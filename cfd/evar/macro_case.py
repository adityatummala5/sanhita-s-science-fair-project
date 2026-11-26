"""OpenFOAM case writer for vessel-scale EVAR runs (straight-tube validation
and the bifurcated endograft). OpenFOAM v2406 (openfoam.com).

Quarter model: the background blockMesh box spans x in [x_in, x_out],
y in [0, y_max], z in [0, z_max]. snappyHexMesh keeps the part of the box
inside the lumen surface, so the box faces become the patches
    inlet  (x = x_in)   Womersley velocity profile (timeVaryingMappedFixedValue)
    outlet (x = x_out)  p = 0, zero-gradient velocity (backflow allowed)
    symY   (y = 0)      symmetry plane between the two limbs
    symZ   (z = 0)      symmetry plane (bifurcation plane)
and the lumen surface becomes the patch `wall` (rigid graft fabric, no-slip).
The mesh is built in mm and scaled to m.
"""
from __future__ import annotations

import math
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "pipeline"))
import case as fc  # noqa: E402  (coronary pipeline: physics/schemes writers, HEADER)

HEADER = fc.HEADER

PRESETS = {
    #           bg(mm)  wall level  layers
    "coarse": dict(bg=1.0, wall=2, layers=3),     # 0.25 mm at the wall
    "medium": dict(bg=0.8, wall=2, layers=4),     # 0.20 mm
    "fine":   dict(bg=0.6, wall=2, layers=5),     # 0.15 mm
}


def _w(path, text):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w") as fh:
        fh.write(text)


def write_mesh_dicts(case, box, preset, nproc, location):
    pr = PRESETS[preset]
    bg = pr["bg"]
    (x0, x1), (y0, y1), (z0, z1) = box
    nx, ny, nz = (math.ceil((b - a) / bg) for a, b in box)
    _w(f"{case}/system/blockMeshDict", HEADER.format(cls="dictionary", obj="blockMeshDict") + f"""
convertToMeters 1;
vertices
(
    ({x0} {y0} {z0}) ({x1} {y0} {z0}) ({x1} {y1} {z0}) ({x0} {y1} {z0})
    ({x0} {y0} {z1}) ({x1} {y0} {z1}) ({x1} {y1} {z1}) ({x0} {y1} {z1})
);
blocks ( hex (0 1 2 3 4 5 6 7) ({nx} {ny} {nz}) simpleGrading (1 1 1) );
boundary
(
    inlet  {{ type patch;    faces ((0 4 7 3)); }}
    outlet {{ type patch;    faces ((1 2 6 5)); }}
    symY   {{ type symmetry; faces ((0 1 5 4)); }}
    symZ   {{ type symmetry; faces ((0 3 2 1)); }}
    box    {{ type wall;     faces ((3 7 6 2) (4 5 6 7)); }}
);
""")
    L = pr["wall"]
    _w(f"{case}/system/snappyHexMeshDict", HEADER.format(cls="dictionary", obj="snappyHexMeshDict") + f"""
castellatedMesh true;
snap            true;
addLayers       {'true' if pr['layers'] else 'false'};
geometry {{ lumen.stl {{ type triSurfaceMesh; name wall; }} }}
castellatedMeshControls
{{
    maxLocalCells 6000000; maxGlobalCells 40000000; minRefinementCells 0;
    maxLoadUnbalance 0.10; nCellsBetweenLevels 3;
    features ();
    refinementSurfaces {{ wall {{ level ({L} {L}); patchInfo {{ type wall; }} }} }}
    resolveFeatureAngle 30;
    refinementRegions {{ wall {{ mode distance; levels ((0.8 {L}) (2.5 {L - 1})); }} }}
    locationInMesh ({location[0]} {location[1]} {location[2]});
    allowFreeStandingZoneFaces true;
}}
snapControls
{{
    nSmoothPatch 3; tolerance 2.0; nSolveIter 100; nRelaxIter 5;
    nFeatureSnapIter 10; implicitFeatureSnap true; explicitFeatureSnap false; multiRegionFeatureSnap false;
}}
addLayersControls
{{
    relativeSizes true;
    layers {{ wall {{ nSurfaceLayers {pr['layers']}; }} }}
    expansionRatio 1.25; finalLayerThickness 0.45; minThickness 0.08; nGrow 0;
    featureAngle 120; slipFeatureAngle 30; nRelaxIter 5; nSmoothSurfaceNormals 1; nSmoothNormals 3;
    nSmoothThickness 10; maxFaceThicknessRatio 0.5; maxThicknessToMedialRatio 0.3;
    minMedialAxisAngle 90; nBufferCellsNoExtrude 0; nLayerIter 50;
}}
meshQualityControls
{{
    maxNonOrtho 65; maxBoundarySkewness 20; maxInternalSkewness 4; maxConcave 80;
    minVol 1e-30; minTetQuality 1e-15; minArea -1; minTwist 0.02; minDeterminant 0.001;
    minFaceWeight 0.05; minVolRatio 0.01; minTriangleTwist -1; nSmoothScale 4; errorReduction 0.75;
    relaxed {{ maxNonOrtho 75; }}
}}
mergeTolerance 1e-6;
""")
    _w(f"{case}/system/decomposeParDict", HEADER.format(cls="dictionary", obj="decomposeParDict") + f"""
numberOfSubdomains {nproc};
method hierarchical;
hierarchicalCoeffs {{ n ({nproc} 1 1); delta 0.001; order xyz; }}
""")
    return {"background_cells": nx * ny * nz, "bg_mm": bg, "wall_cell_mm": bg / 2 ** L}


def _tname(t: float) -> str:
    s = f"{t:.6f}".rstrip("0").rstrip(".")
    return s if s else "0"


def write_inlet_womersley(case, W, times, x_in_mm, R_mm, n_r=40, n_t=24):
    """Womersley profile at the inlet plane for each time in `times` (s).
    Points cover the +y,+z quadrant of the inlet (m)."""
    R = R_mm * 1e-3
    pts = [(0.0, 0.0)]
    for r in np.linspace(R / n_r, R * 1.04, n_r):
        for th in np.linspace(-0.05, math.pi / 2 + 0.05, n_t):
            pts.append((r * math.cos(th), r * math.sin(th)))
    pts = np.array(pts)
    rr = np.hypot(pts[:, 0], pts[:, 1])
    inside = rr <= R
    base = f"{case}/constant/boundaryData/inlet"
    os.makedirs(base, exist_ok=True)
    with open(f"{base}/points", "w") as fh:
        fh.write(f"{len(pts)}\n(\n" + "\n".join(f"({x_in_mm * 1e-3:.9e} {p[0]:.9e} {p[1]:.9e})" for p in pts) + "\n)\n")
    for t in times:
        u = np.zeros(len(pts))
        u[inside] = W.velocity(np.minimum(rr[inside], R), t)
        d = f"{base}/{_tname(t)}"
        os.makedirs(d, exist_ok=True)
        with open(f"{d}/U", "w") as fh:
            fh.write(f"{len(pts)}\n(\n" + "\n".join(f"({v:.9e} 0 0)" for v in u) + "\n)\n")


def write_inlet_poiseuille(case, q_m3_s, x_in_mm, R_mm):
    """Steady Poiseuille inlet (for the mean-flow initialisation run)."""
    class _P:
        def velocity(self, r, t):
            R = R_mm * 1e-3
            return 2 * q_m3_s / (math.pi * R ** 2) * (1 - (r / R) ** 2)
    write_inlet_womersley(case, _P(), [0.0, 1e6], x_in_mm, R_mm)


def write_fields(case):
    _w(f"{case}/0/U", HEADER.format(cls="volVectorField", obj="U") + """
dimensions [0 1 -1 0 0 0 0];
internalField uniform (0 0 0);
boundaryField
{
    inlet
    {
        type timeVaryingMappedFixedValue; offset (0 0 0); setAverage off;
        mapMethod planarInterpolation; value uniform (0 0 0);
    }
    outlet { type zeroGradient; }
    wall   { type noSlip; }
    box    { type noSlip; }
    symY   { type symmetry; }
    symZ   { type symmetry; }
}
""")
    _w(f"{case}/0/p", HEADER.format(cls="volScalarField", obj="p") + """
dimensions [0 2 -2 0 0 0 0];
internalField uniform 0;
boundaryField
{
    inlet  { type zeroGradient; }
    outlet { type fixedValue; value uniform 0; }
    wall   { type zeroGradient; }
    box    { type zeroGradient; }
    symY   { type symmetry; }
    symZ   { type symmetry; }
}
""")


def _fos(avg_start, probes_m=None):
    fos = """
    wallShearStress { type wallShearStress; libs ("libfieldFunctionObjects.so"); patches (wall);
                      executeControl timeStep; writeControl writeTime; }
    magWSS { type mag; libs ("libfieldFunctionObjects.so"); field wallShearStress; result magWSS;
             executeControl timeStep; writeControl writeTime; }
    inletFlow { type surfaceFieldValue; libs ("libfieldFunctionObjects.so"); regionType patch; name inlet;
                operation sum; fields (phi); writeFields false; writeControl timeStep; writeInterval 5; }
    outletFlow { $inletFlow; name outlet; }
    inletPressure { type surfaceFieldValue; libs ("libfieldFunctionObjects.so"); regionType patch; name inlet;
                    operation areaAverage; fields (p); writeFields false; writeControl timeStep; writeInterval 5; }
"""
    if probes_m:
        pl = " ".join(f"({x:.6e} {y:.6e} {z:.6e})" for x, y, z in probes_m)
        fos += f"""
    wssProbes {{ type patchProbes; libs ("libsampling.so"); patches (wall); fields (wallShearStress);
                probeLocations ({pl}); executeControl timeStep; writeControl timeStep; writeInterval 2; }}
"""
    if avg_start is not None:
        fos += f"""
    cycleAverage
    {{
        type fieldAverage; libs ("libfieldFunctionObjects.so"); timeStart {avg_start};
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


def write_control_pulsatile(case, period, n_cycles, max_co, snapshots, t_start=0.0, probes_m=None):
    end = t_start + period * n_cycles
    _w(f"{case}/system/controlDict", HEADER.format(cls="dictionary", obj="controlDict") + f"""
application     pimpleFoam;
startFrom       latestTime;
startTime       {t_start};
stopAt          endTime;
endTime         {end:.6f};
deltaT          1e-4;
adjustTimeStep  yes;
maxCo           {max_co};
maxDeltaT       {period / 400:.6f};
writeControl    adjustableRunTime;
writeInterval   {period / snapshots:.6f};
purgeWrite      {snapshots};
writeFormat     binary;
writePrecision  8;
writeCompression off;
timeFormat      general;
timePrecision   8;
runTimeModifiable true;
""" + _fos(avg_start=t_start + period * (n_cycles - 1), probes_m=probes_m))


def write_control_steady(case, n_iter):
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
""" + _fos(avg_start=None))


def write_physics(case, ph):
    fc.write_physics(case, ph)


def write_schemes(case, steady):
    fc.write_schemes(case, steady)
    if not steady:      # extra outer correctors: pulsatile flow with strong reversal
        p = f"{case}/system/fvSolution"
        s = open(p).read().replace("nOuterCorrectors 1;", "nOuterCorrectors 2;")
        open(p, "w").write(s)
