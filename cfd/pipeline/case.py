"""OpenFOAM case writer (tested with OpenFOAM v2406, openfoam.com flavour).

Domain: a wedge of artery lumen (default 120 deg, matching Sanhita's CFD
sector) containing one sector of the stent. Geometry is meshed in mm and
then scaled to metres with transformPoints.

Boundary patches after meshing:
    inlet, outlet      planar ends (from blockMesh)
    arteryWall         artery wall (rigid, no-slip)
    symA, symB         sector side planes (symmetry)
    stent              stent struts (+ denticles), no-slip
"""
from __future__ import annotations

import math
import os
from dataclasses import dataclass

import numpy as np

HEADER = """FoamFile
{{
    version     2.0;
    format      ascii;
    class       {cls};
    object      {obj};
}}
"""

# Mesh presets: background cell size (mm) and refinement levels.
# Cell size at level L = bg / 2**L.
PRESETS = {
    #           bg(mm)  wall  stent  (stent cell size)
    "coarse": dict(bg=0.080, wall=1, stent=3),   # 10.0 um at the stent
    "medium": dict(bg=0.060, wall=1, stent=3),   #  7.5 um
    "fine":   dict(bg=0.050, wall=2, stent=4),   #  3.1 um
}


@dataclass
class Physics:
    rho: float = 1060.0           # kg/m3
    # Carreau model (Cho & Kensey 1991 values, the usual choice for blood)
    mu0: float = 0.056            # Pa s
    mu_inf: float = 0.00345       # Pa s
    lam: float = 3.313            # s
    n: float = 0.3568
    newtonian: bool = False
    mu_newtonian: float = 0.0035  # Pa s, used if newtonian=True


@dataclass
class Domain:
    x0: float                    # stent window start (mm)
    x1: float                    # stent window end (mm)
    theta0_deg: float
    theta1_deg: float
    radius: float = 1.59         # lumen radius (mm); struts (r 1.50-1.60) embed 10 um in the wall
    l_in: float = 1.0            # inlet extension upstream of the stent (mm)
    l_out: float = 6.0           # outlet extension downstream (mm)

    @property
    def xa(self):
        return self.x0 - self.l_in

    @property
    def xb(self):
        return self.x1 + self.l_out


def _w(path, text):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w") as fh:
        fh.write(text)


# --------------------------------------------------------------------------
# Meshing dictionaries
# --------------------------------------------------------------------------
def write_mesh_dicts(case: str, dom: Domain, preset: str, nproc: int, has_stent: bool = True):
    pr = PRESETS[preset]
    bg = pr["bg"]
    th = np.radians(np.linspace(dom.theta0_deg, dom.theta1_deg, 721))
    ys, zs = dom.radius * np.cos(th), dom.radius * np.sin(th)
    ymin, ymax = min(ys.min(), 0.0) - 2 * bg, ys.max() + 2 * bg
    zmin, zmax = min(zs.min(), 0.0) - 2 * bg, zs.max() + 2 * bg
    nx = math.ceil((dom.xb - dom.xa) / bg)
    ny = math.ceil((ymax - ymin) / bg)
    nz = math.ceil((zmax - zmin) / bg)

    _w(f"{case}/system/blockMeshDict", HEADER.format(cls="dictionary", obj="blockMeshDict") + f"""
convertToMeters 1;   // mm here; scaled to m after snappyHexMesh
vertices
(
    ({dom.xa} {ymin} {zmin}) ({dom.xb} {ymin} {zmin}) ({dom.xb} {ymax} {zmin}) ({dom.xa} {ymax} {zmin})
    ({dom.xa} {ymin} {zmax}) ({dom.xb} {ymin} {zmax}) ({dom.xb} {ymax} {zmax}) ({dom.xa} {ymax} {zmax})
);
blocks ( hex (0 1 2 3 4 5 6 7) ({nx} {ny} {nz}) simpleGrading (1 1 1) );
boundary
(
    inlet  {{ type patch; faces ((0 4 7 3)); }}
    outlet {{ type patch; faces ((1 2 6 5)); }}
    box    {{ type wall;  faces ((0 1 5 4) (3 7 6 2) (0 3 2 1) (4 5 6 7)); }}
);
""")

    _w(f"{case}/system/surfaceFeatureExtractDict", HEADER.format(cls="dictionary", obj="surfaceFeatureExtractDict") + """
lumen.stl { extractionMethod extractFromSurface; extractFromSurfaceCoeffs { includedAngle 150; } writeObj no; }
""" + ("stent.stl { extractionMethod extractFromSurface; extractFromSurfaceCoeffs { includedAngle 150; } writeObj no; }\n"
       if has_stent else ""))

    L = pr["stent"]
    lw = pr["wall"]
    # point inside the fluid, away from struts: mid-radius, mid-angle, in the inlet extension
    tm = math.radians(0.5 * (dom.theta0_deg + dom.theta1_deg) + 0.37)
    loc = (dom.xa + 0.5 * dom.l_in + 0.0013, 0.61 * dom.radius * math.cos(tm), 0.61 * dom.radius * math.sin(tm))
    if not has_stent:          # validation case: refine the wall one level below the stent level
        lw = max(lw, L - 1)
    st_geo = "stent.stl { type triSurfaceMesh; name stent; }" if has_stent else ""
    st_feat = f'{{ file "stent.eMesh"; level {L}; }}' if has_stent else ""
    st_surf = f"stent {{ level ({L} {L}); patchInfo {{ type wall; }} }}" if has_stent else ""
    st_reg = (f"stent {{ mode distance; levels ((0.015 {L}) (0.05 {max(L - 1, 0)}) (0.12 {max(L - 2, 0)})); }}"
              if has_stent else "")
    _w(f"{case}/system/snappyHexMeshDict", HEADER.format(cls="dictionary", obj="snappyHexMeshDict") + f"""
castellatedMesh true;
snap            true;
addLayers       false;

geometry
{{
    lumen.stl
    {{
        type triSurfaceMesh;
        name lumen;
        regions {{ wall {{ name arteryWall; }} symA {{ name symA; }} symB {{ name symB; }} caps {{ name caps; }} }}
    }}
    {st_geo}
}}

castellatedMeshControls
{{
    maxLocalCells 4000000;
    maxGlobalCells 30000000;
    minRefinementCells 0;
    maxLoadUnbalance 0.10;
    nCellsBetweenLevels 3;
    features
    (
        {{ file "lumen.eMesh"; level {lw}; }}
        {st_feat}
    );
    refinementSurfaces
    {{
        lumen
        {{
            level ({lw} {lw});
            regions
            {{
                wall {{ level ({lw} {lw}); patchInfo {{ type wall; }} }}
                symA {{ level (0 0); patchInfo {{ type symmetry; }} }}
                symB {{ level (0 0); patchInfo {{ type symmetry; }} }}
                caps {{ level (0 0); patchInfo {{ type patch; }} }}
            }}
        }}
        {st_surf}
    }}
    resolveFeatureAngle 30;
    refinementRegions
    {{
        {st_reg}
    }}
    locationInMesh ({loc[0]:.6f} {loc[1]:.6f} {loc[2]:.6f});
    allowFreeStandingZoneFaces true;
}}

snapControls
{{
    nSmoothPatch 3;
    tolerance 2.0;
    nSolveIter 100;
    nRelaxIter 5;
    nFeatureSnapIter 15;
    implicitFeatureSnap false;
    explicitFeatureSnap true;
    multiRegionFeatureSnap false;
}}

addLayersControls
{{
    relativeSizes true;
    layers {{ }}
    expansionRatio 1.2;
    finalLayerThickness 0.5;
    minThickness 0.1;
    nGrow 0;
    featureAngle 60;
    nRelaxIter 3;
    nSmoothSurfaceNormals 1;
    nSmoothNormals 3;
    nSmoothThickness 10;
    maxFaceThicknessRatio 0.5;
    maxThicknessToMedialRatio 0.3;
    minMedialAxisAngle 90;
    nBufferCellsNoExtrude 0;
    nLayerIter 50;
}}

meshQualityControls
{{
    maxNonOrtho 65;
    maxBoundarySkewness 20;
    maxInternalSkewness 4;
    maxConcave 80;
    minVol 1e-30;
    minTetQuality 1e-15;
    minArea -1;
    minTwist 0.02;
    minDeterminant 0.001;
    minFaceWeight 0.05;
    minVolRatio 0.01;
    minTriangleTwist -1;
    nSmoothScale 4;
    errorReduction 0.75;
    relaxed {{ maxNonOrtho 75; }}
}}

mergeTolerance 1e-6;
""")
    _w(f"{case}/system/decomposeParDict", HEADER.format(cls="dictionary", obj="decomposeParDict") + f"""
numberOfSubdomains {nproc};
method hierarchical;   // scotch is not available in every OpenFOAM build
hierarchicalCoeffs {{ n ({nproc} 1 1); delta 0.001; order xyz; }}
""")
    return {"background_cells": nx * ny * nz, "bg_mm": bg,
            "stent_cell_um": 1e3 * bg / 2 ** L, "wall_cell_um": 1e3 * bg / 2 ** lw}


# --------------------------------------------------------------------------
# Inlet: time-varying parabolic profile via timeVaryingMappedFixedValue
# --------------------------------------------------------------------------
def write_inlet_data(case: str, dom: Domain, times: np.ndarray, q_mlmin: np.ndarray):
    """Parabolic (Poiseuille) profile u(r,t) = 2*Umean(t)*(1 - r^2/R^2) on the
    inlet plane. Umean is computed for the FULL vessel cross-section, so the
    wedge automatically carries its share of the flow."""
    R = dom.radius * 1e-3
    xa = dom.xa * 1e-3
    pts = [(xa, 0.0, 0.0)]
    rr = np.linspace(R / 30, R * 1.05, 32)
    tt = np.radians(np.linspace(dom.theta0_deg - 2, dom.theta1_deg + 2, 63))
    for r in rr:
        for t in tt:
            pts.append((xa, r * math.cos(t), r * math.sin(t)))
    pts = np.array(pts)
    r_pts = np.hypot(pts[:, 1], pts[:, 2])
    shape = np.clip(2.0 * (1.0 - (r_pts / R) ** 2), 0.0, None)
    base = f"{case}/constant/boundaryData/inlet"
    os.makedirs(base, exist_ok=True)
    with open(f"{base}/points", "w") as fh:
        fh.write(f"{len(pts)}\n(\n" + "\n".join(f"({p[0]:.9e} {p[1]:.9e} {p[2]:.9e})" for p in pts) + "\n)\n")
    area = math.pi * R ** 2
    for t, q in zip(times, q_mlmin):
        um = q * 1e-6 / 60.0 / area
        d = f"{base}/{_tname(t)}"
        os.makedirs(d, exist_ok=True)
        with open(f"{d}/U", "w") as fh:
            fh.write(f"{len(pts)}\n(\n" + "\n".join(f"({um * s:.9e} 0 0)" for s in shape) + "\n)\n")


def _tname(t: float) -> str:
    s = f"{t:.6f}".rstrip("0").rstrip(".")
    return s if s else "0"


# --------------------------------------------------------------------------
# Fields, physics, solver controls
# --------------------------------------------------------------------------
WALLS = ("arteryWall", "stent", "box")
SYMS = ("symA", "symB")


def write_fields(case: str, has_stent: bool = True):
    fos_patches = "(arteryWall stent)" if has_stent else "(arteryWall)"
    cd = os.path.join(case, "system/controlDict")
    if os.path.exists(cd):     # function objects must not list a patch that does not exist
        txt = open(cd).read().replace("patches (arteryWall stent);", f"patches {fos_patches};")
        open(cd, "w").write(txt)
    U_bc = "\n".join(f"    {w} {{ type noSlip; }}" for w in WALLS)
    U_bc += "\n" + "\n".join(f"    {s} {{ type symmetry; }}" for s in SYMS)
    U_bc += "\n    caps { type noSlip; }"
    _w(f"{case}/0/U", HEADER.format(cls="volVectorField", obj="U") + f"""
dimensions [0 1 -1 0 0 0 0];
internalField uniform (0 0 0);
boundaryField
{{
    inlet
    {{
        type            timeVaryingMappedFixedValue;
        offset          (0 0 0);
        setAverage      off;
        mapMethod       planarInterpolation;
        value           uniform (0 0 0);
    }}
    outlet {{ type inletOutlet; inletValue uniform (0 0 0); value uniform (0 0 0); }}
{U_bc}
}}
""")
    p_bc = "\n".join(f"    {w} {{ type zeroGradient; }}" for w in WALLS)
    p_bc += "\n" + "\n".join(f"    {s} {{ type symmetry; }}" for s in SYMS)
    p_bc += "\n    caps { type zeroGradient; }"
    _w(f"{case}/0/p", HEADER.format(cls="volScalarField", obj="p") + f"""
dimensions [0 2 -2 0 0 0 0];   // kinematic pressure p/rho
internalField uniform 0;
boundaryField
{{
    inlet  {{ type zeroGradient; }}
    outlet {{ type fixedValue; value uniform 0; }}
{p_bc}
}}
""")


def write_physics(case: str, ph: Physics):
    if ph.newtonian:
        tp = f"transportModel Newtonian;\nnu {ph.mu_newtonian / ph.rho:.6e};\n"
    else:
        tp = f"""transportModel BirdCarreau;
BirdCarreauCoeffs
{{
    nu0   [0 2 -1 0 0 0 0] {ph.mu0 / ph.rho:.6e};
    nuInf [0 2 -1 0 0 0 0] {ph.mu_inf / ph.rho:.6e};
    k     [0 0 1 0 0 0 0]  {ph.lam};
    n     [0 0 0 0 0 0 0]  {ph.n};
}}
"""
    _w(f"{case}/constant/transportProperties", HEADER.format(cls="dictionary", obj="transportProperties") + tp)
    _w(f"{case}/constant/turbulenceProperties", HEADER.format(cls="dictionary", obj="turbulenceProperties")
       + "simulationType laminar;\n")


def write_schemes(case: str, steady: bool):
    _w(f"{case}/system/fvSchemes", HEADER.format(cls="dictionary", obj="fvSchemes") + f"""
ddtSchemes      {{ default {'steadyState' if steady else 'backward'}; }}
gradSchemes     {{ default Gauss linear; grad(U) cellLimited Gauss linear 1; }}
divSchemes
{{
    default         none;
    div(phi,U)      {'bounded ' if steady else ''}Gauss linearUpwind grad(U);
    div((nuEff*dev2(T(grad(U))))) Gauss linear;
}}
laplacianSchemes {{ default Gauss linear limited corrected 0.5; }}
interpolationSchemes {{ default linear; }}
snGradSchemes   {{ default limited corrected 0.5; }}
wallDist        {{ method meshWave; }}
""")
    _w(f"{case}/system/fvSolution", HEADER.format(cls="dictionary", obj="fvSolution") + """
solvers
{
    p      { solver GAMG; smoother GaussSeidel; tolerance 1e-7; relTol 0.01; }
    pFinal { $p; relTol 0; }
    "(U|UFinal)" { solver smoothSolver; smoother symGaussSeidel; tolerance 1e-8; relTol 0.1; }
}
PIMPLE
{
    momentumPredictor yes;
    nOuterCorrectors 1;
    nCorrectors 2;
    nNonOrthogonalCorrectors 1;
}
SIMPLE
{
    nNonOrthogonalCorrectors 1;
    consistent no;
    residualControl { p 1e-5; U 1e-6; }
}
relaxationFactors
{
    fields { p 0.3; }
    equations { U 0.7; }
}
""")


def _function_objects(avg_start: float | None) -> str:
    fos = """
    wallShearStress
    {
        type wallShearStress; libs ("libfieldFunctionObjects.so");
        patches (arteryWall stent);
        executeControl timeStep; writeControl writeTime;
    }
    magWSS
    {
        type mag; libs ("libfieldFunctionObjects.so");
        field wallShearStress; result magWSS;
        executeControl timeStep; writeControl writeTime;
    }
    inletFlow
    {
        type surfaceFieldValue; libs ("libfieldFunctionObjects.so");
        regionType patch; name inlet; operation sum; fields (phi);
        writeFields false; writeControl timeStep; writeInterval 5;
    }
    outletFlow { $inletFlow; name outlet; }
    inletPressure
    {
        type surfaceFieldValue; libs ("libfieldFunctionObjects.so");
        regionType patch; name inlet; operation areaAverage; fields (p);
        writeFields false; writeControl timeStep; writeInterval 5;
    }
"""
    if avg_start is not None:
        fos += f"""
    cycleAverage
    {{
        type fieldAverage; libs ("libfieldFunctionObjects.so");
        timeStart {avg_start};
        restartOnRestart false; restartOnOutput false; periodicRestart false;
        executeControl timeStep; writeControl writeTime;
        fields
        (
            wallShearStress {{ mean on; prime2Mean off; base time; }}
            magWSS          {{ mean on; prime2Mean off; base time; }}
            U               {{ mean on; prime2Mean off; base time; }}
            p               {{ mean on; prime2Mean off; base time; }}
        );
    }}
"""
    return "functions\n{" + fos + "}\n"


def write_control_pulsatile(case: str, period: float, n_cycles: int, max_co: float, snapshots: int):
    end = period * n_cycles
    wi = period / snapshots
    _w(f"{case}/system/controlDict", HEADER.format(cls="dictionary", obj="controlDict") + f"""
application     pimpleFoam;
startFrom       latestTime;
startTime       0;
stopAt          endTime;
endTime         {end:.6f};
deltaT          1e-5;
adjustTimeStep  yes;
maxCo           {max_co};
maxDeltaT       {period / 100:.6f};
writeControl    adjustableRunTime;
writeInterval   {wi:.6f};
purgeWrite      {snapshots};
writeFormat     binary;
writePrecision  8;
writeCompression off;
timeFormat      general;
timePrecision   8;
runTimeModifiable true;
""" + _function_objects(avg_start=period * (n_cycles - 1)))


def write_control_steady(case: str, n_iter: int):
    _w(f"{case}/system/controlDict", HEADER.format(cls="dictionary", obj="controlDict") + f"""
application     simpleFoam;
startFrom       latestTime;
startTime       0;
stopAt          endTime;
endTime         {n_iter};
deltaT          1;
writeControl    timeStep;
writeInterval   {n_iter};
purgeWrite      1;
writeFormat     binary;
writePrecision  8;
timeFormat      general;
runTimeModifiable true;
""" + _function_objects(avg_start=None))
