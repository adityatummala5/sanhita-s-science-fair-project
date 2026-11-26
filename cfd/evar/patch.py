"""Wall-patch model: Stokes-regime blood flow over a periodic textured graft surface.

Why a patch: the graft is 12-25 mm across but textures are 10-100 um, so the
texture cannot be resolved in the vessel-scale mesh. The patch is a small
periodic box on the graft wall, driven by the shear-stress waveform computed
at a chosen site in the vessel-scale model.

Domain  x (streamwise, along the site's principal shear direction) in [0, Lx],
        z (spanwise) in [0, Lz], both periodic (cyclic) and integer multiples
        of the texture period; y from the textured wall h(x, z) to a lid at
        y = Hlid above the nominal surface.
Drive   the lid moves at U(t) (Couette). Hlid = 0.3 mm is much thinner than the
        oscillatory boundary layer of the graft limb (~1 mm), so the near-wall
        velocity profile is close to linear; and much thicker than the texture
        disturbance (which decays over ~ period/2pi ~ 15-50 um).
Regime  texture Reynolds number gamma * h^2 / nu ~ 0.01-0.1 and patch Womersley
        number Hlid * sqrt(omega/nu) ~ 0.2-0.4, so the flow is quasi-steady
        Stokes flow. Each design is therefore solved as steady flow at a set of
        lid speeds (shear levels); the cardiac cycle is rebuilt by interpolating
        between levels with the site waveform. `--verify-pulsatile` checks this
        against a full time-dependent run.
Mesh    blockMesh box in (x, eta, z), eta in [0, 1] graded toward the wall,
        then every point is moved to y = h + eta (Hlid - h): a surface-following
        structured hex mesh. Periodicity is exact because h is periodic.
Calibration  for the smooth wall the lid speed giving wall stress tau is
        U = Hlid * gamma(tau) with gamma from the Carreau law tau = mu(gamma) gamma.
        The same U(t) drives every textured design, so all designs see the same
        outer flow and differ only by their surface.
"""
from __future__ import annotations

import glob
import json
import math
import os
import re
import shutil
import sys
from dataclasses import dataclass

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "pipeline"))
import case as fc  # noqa: E402

HEADER = fc.HEADER


@dataclass
class PatchMesh:
    h_lid_mm: float = 0.30
    dx_um: float = 3.0          # horizontal cell size
    ny: int = 40                # cells across the gap
    grading: float = 30.0       # last/first cell height (first cell ~1 um for a 0.3 mm gap)


def carreau_mu(g, ph):
    return ph.mu_inf + (ph.mu0 - ph.mu_inf) * (1 + (ph.lam * g) ** 2) ** ((ph.n - 1) / 2)


def gamma_from_tau(tau, ph):
    """Invert tau = mu(gamma) * gamma (monotonic for n > 0) by bisection, vectorised."""
    tau = np.abs(np.asarray(tau, float))
    lo, hi = np.zeros_like(tau), np.full_like(tau, 1e5)
    for _ in range(80):
        mid = 0.5 * (lo + hi)
        f = carreau_mu(mid, ph) * mid
        lo = np.where(f < tau, mid, lo)
        hi = np.where(f >= tau, mid, hi)
    return 0.5 * (lo + hi)


def _w(path, text):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w") as fh:
        fh.write(text)


# ---------------------------------------------------------------------------
def write_block(cdir, Lx, Lz, pm: PatchMesh):
    nx = max(4, round(Lx * 1e3 / pm.dx_um))
    nz = max(4, round(Lz * 1e3 / pm.dx_um))
    _w(f"{cdir}/system/blockMeshDict", HEADER.format(cls="dictionary", obj="blockMeshDict") + f"""
convertToMeters 1;
vertices
(
    (0 0 0) ({Lx} 0 0) ({Lx} 1 0) (0 1 0)
    (0 0 {Lz}) ({Lx} 0 {Lz}) ({Lx} 1 {Lz}) (0 1 {Lz})
);
blocks ( hex (0 1 2 3 4 5 6 7) ({nx} {pm.ny} {nz}) simpleGrading (1 {pm.grading} 1) );
boundary
(
    graft {{ type wall; faces ((0 1 5 4)); }}
    lid   {{ type wall; faces ((3 7 6 2)); }}
    xLo   {{ type cyclic; neighbourPatch xHi; faces ((0 4 7 3)); }}
    xHi   {{ type cyclic; neighbourPatch xLo; faces ((1 2 6 5)); }}
    zLo   {{ type cyclic; neighbourPatch zHi; faces ((0 3 2 1)); }}
    zHi   {{ type cyclic; neighbourPatch zLo; faces ((4 5 6 7)); }}
);
""")
    return nx, nz


def map_points(cdir, tex, pm: PatchMesh):
    """Move blockMesh points (x, eta, z) [mm] onto the textured wall and scale to metres."""
    p = os.path.join(cdir, "constant/polyMesh/points")
    txt = open(p).read()
    m = re.search(r"\n(\d+)\s*\n\(\n", txt)
    n = int(m.group(1))
    body = txt[m.end():]
    arr = np.array(re.findall(r"\(([-\d.eE+]+) ([-\d.eE+]+) ([-\d.eE+]+)\)", body[: body.rfind(")")]), float)
    assert len(arr) == n, (len(arr), n)
    X, E, Z = arr[:, 0], arr[:, 1], arr[:, 2]
    Lx, Lz = tex.domain()
    h = tex.height(X, Z)          # periodic: x = Lx gives the same value as x = 0
    Y = h + E * (pm.h_lid_mm - h)
    out = np.c_[X, Y, Z] * 1e-3
    head = txt[: m.start()]
    with open(p, "w") as fh:
        fh.write(head + f"\n{n}\n(\n" + "\n".join(f"({a:.10e} {b:.10e} {c:.10e})" for a, b, c in out) + "\n)\n")
    # initial guess for the solves: linear Couette profile between wall and lid (cell order = blockMesh order)
    nx, nz = round(Lx * 1e3 / pm.dx_um), round(Lz * 1e3 / pm.dx_um)
    nx, nz = max(4, nx), max(4, nz)
    E3 = E.reshape(nz + 1, pm.ny + 1, nx + 1)
    Ec = 0.125 * sum(E3[a:a + nz, b:b + pm.ny, c:c + nx] for a in (0, 1) for b in (0, 1) for c in (0, 1))
    np.save(os.path.join(cdir, "eta_cells.npy"), Ec.ravel())


def write_common(cdir, ph, nproc=1):
    fc.write_physics(cdir, ph)
    _w(f"{cdir}/system/fvSchemes", HEADER.format(cls="dictionary", obj="fvSchemes") + """
ddtSchemes      { default steadyState; }
gradSchemes     { default Gauss linear; }
divSchemes      { default none; div(phi,U) bounded Gauss linear; div((nuEff*dev2(T(grad(U))))) Gauss linear; }
laplacianSchemes { default Gauss linear corrected; }
interpolationSchemes { default linear; }
snGradSchemes   { default corrected; }
""")
    _w(f"{cdir}/system/fvSolution", HEADER.format(cls="dictionary", obj="fvSolution") + """
solvers
{
    p { solver GAMG; smoother GaussSeidel; tolerance 1e-9; relTol 0.02; }
    U { solver smoothSolver; smoother symGaussSeidel; tolerance 1e-10; relTol 0.05; }
}
SIMPLE
{
    nNonOrthogonalCorrectors 2;
    consistent yes;
    pRefCell 0;
    pRefValue 0;
    residualControl { U 1e-7; }
}
relaxationFactors { equations { U 0.98; ".*" 0.98; } fields { p 1.0; } }
""")
    _w(f"{cdir}/system/decomposeParDict", HEADER.format(cls="dictionary", obj="decomposeParDict") + f"""
numberOfSubdomains {nproc};
method hierarchical;
hierarchicalCoeffs {{ n ({nproc} 1 1); delta 0.001; order xyz; }}
""")


def read_internal_U(path):
    txt = open(path).read()
    i = txt.index("internalField")
    m = re.search(r"nonuniform List<vector>\s*(\d+)\s*\(", txt[i:])
    n = int(m.group(1))
    body = txt[i + m.end():]
    vals = np.array(re.findall(r"\(([-\d.eE+]+) ([-\d.eE+]+) ([-\d.eE+]+)\)", body[: body.index("\n)\n")]), float)
    assert len(vals) == n
    return vals


def write_fields(cdir, u_lid, eta_cells=None, internal=None):
    if internal is not None:
        vals = "\n".join(f"({a:.8e} {b:.8e} {c:.8e})" for a, b, c in internal)
        internal = f"nonuniform List<vector>\n{len(internal)}\n(\n{vals}\n)"
    elif eta_cells is not None:
        vals = "\n".join(f"({u_lid * e:.8e} 0 0)" for e in eta_cells)
        internal = f"nonuniform List<vector>\n{len(eta_cells)}\n(\n{vals}\n)"
    else:
        internal = f"uniform ({0.5 * u_lid:.6e} 0 0)"
    cyc = "\n".join(f"    {c} {{ type cyclic; }}" for c in ("xLo", "xHi", "zLo", "zHi"))
    _w(f"{cdir}/0/U", HEADER.format(cls="volVectorField", obj="U") + f"""
dimensions [0 1 -1 0 0 0 0];
internalField {internal};
boundaryField
{{
    graft {{ type noSlip; }}
    lid   {{ type fixedValue; value uniform ({u_lid:.8e} 0 0); }}
{cyc}
}}
""")
    _w(f"{cdir}/0/p", HEADER.format(cls="volScalarField", obj="p") + f"""
dimensions [0 2 -2 0 0 0 0];
internalField uniform 0;
boundaryField
{{
    graft {{ type zeroGradient; }}
    lid   {{ type zeroGradient; }}
{cyc}
}}
""")


def write_control(cdir, n_iter, ascii_points=True):
    _w(f"{cdir}/system/controlDict", HEADER.format(cls="dictionary", obj="controlDict") + f"""
application     simpleFoam;
startFrom       latestTime;
startTime       0;
stopAt          endTime;
endTime         {n_iter};
deltaT          1;
writeControl    timeStep;
writeInterval   {n_iter};
purgeWrite      1;
writeFormat     {'ascii' if ascii_points else 'binary'};
writePrecision  10;
timeFormat      general;
runTimeModifiable true;
functions
{{
    wallShearStress {{ type wallShearStress; libs ("libfieldFunctionObjects.so"); patches (graft lid);
                      executeControl timeStep; executeInterval 25; writeControl writeTime; }}
    wallStress {{ type surfaceFieldValue; libs ("libfieldFunctionObjects.so"); regionType patch; name graft;
                 operation areaAverage; fields (wallShearStress); writeFields false;
                 executeControl timeStep; executeInterval 25; writeControl timeStep; writeInterval 25; }}
    writeCellVolumes {{ type writeCellVolumes; libs ("libfieldFunctionObjects.so"); writeControl writeTime; }}
    writeCellCentres {{ type writeCellCentres; libs ("libfieldFunctionObjects.so"); writeControl writeTime; }}
}}
""")


def write_pulsatile(cdir, t, U_t, period, n_cycles, eta_cells):
    """Time-dependent check of the quasi-steady assumption: same patch, lid
    driven by the full U(t) table, pimpleFoam with a fixed step of T/400."""
    rows = "\n".join(f"        ({k * period + a:.6f} ({u:.8e} 0 0))"
                     for k in range(n_cycles + 1) for a, u in zip(t, U_t))
    cyc = "\n".join(f"    {c} {{ type cyclic; }}" for c in ("xLo", "xHi", "zLo", "zHi"))
    vals = "\n".join(f"({U_t[0] * e:.8e} 0 0)" for e in eta_cells)
    _w(f"{cdir}/0/U", HEADER.format(cls="volVectorField", obj="U") + f"""
dimensions [0 1 -1 0 0 0 0];
internalField nonuniform List<vector>
{len(eta_cells)}
(
{vals}
);
boundaryField
{{
    graft {{ type noSlip; }}
    lid
    {{
        type uniformFixedValue;
        uniformValue table
        (
{rows}
        );
    }}
{cyc}
}}
""")
    _w(f"{cdir}/0/p", HEADER.format(cls="volScalarField", obj="p") + f"""
dimensions [0 2 -2 0 0 0 0];
internalField uniform 0;
boundaryField
{{
    graft {{ type zeroGradient; }}
    lid   {{ type zeroGradient; }}
{cyc}
}}
""")
    s = open(f"{cdir}/system/fvSchemes").read().replace("default steadyState;", "default Euler;")
    open(f"{cdir}/system/fvSchemes", "w").write(s)
    s = open(f"{cdir}/system/fvSolution").read()
    s = s.replace("SIMPLE\n{", "PIMPLE\n{\n    nOuterCorrectors 3;\n    nCorrectors 2;\n    momentumPredictor yes;")
    s = s.replace('U { solver', 'pFinal { solver GAMG; smoother GaussSeidel; tolerance 1e-9; relTol 0; }\n    "(U|UFinal)" { solver')
    open(f"{cdir}/system/fvSolution", "w").write(s)
    end = period * n_cycles
    _w(f"{cdir}/system/controlDict", HEADER.format(cls="dictionary", obj="controlDict") + f"""
application     pimpleFoam;
startFrom       startTime;
startTime       0;
stopAt          endTime;
endTime         {end:.6f};
deltaT          {period / 400:.6f};
adjustTimeStep  no;
writeControl    runTime;
writeInterval   {period:.6f};
purgeWrite      1;
writeFormat     binary;
writePrecision  10;
timeFormat      general;
runTimeModifiable true;
functions
{{
    wallShearStress {{ type wallShearStress; libs ("libfieldFunctionObjects.so"); patches (graft);
                      executeControl timeStep; writeControl writeTime; }}
    magWSS {{ type mag; libs ("libfieldFunctionObjects.so"); field wallShearStress; result magWSS;
             executeControl timeStep; writeControl writeTime; }}
    cycleAverage
    {{
        type fieldAverage; libs ("libfieldFunctionObjects.so"); timeStart {period * (n_cycles - 1):.6f};
        restartOnRestart false; restartOnOutput false; periodicRestart false;
        executeControl timeStep; writeControl writeTime;
        fields ( wallShearStress {{ mean on; prime2Mean off; base time; }} magWSS {{ mean on; prime2Mean off; base time; }} );
    }}
}}
""")
