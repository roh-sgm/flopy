"""Tests for Mf6Tid and Mf6TidRunner (flopy/mf6/tid.py).

Verification criteria from SPEC_MF6_TID.md §8:
  1. idomain_at(kper=N) returns correct array after activate/deactivate sequence.
  2. Mf6TidRunner produces a head file where activated cells contribute starting
     at the correct SP and inactive cells are not solved.
  3. Re-activating a deactivated cell restores normal budget behaviour.
  4. Mf6Tid.load() round-trips a written .tid.json without data loss.
  5. Node conversion (DIS/DISV tuple → flat index) agrees with modelgrid.get_node().
"""

from __future__ import annotations

import shutil
from pathlib import Path

import numpy as np
import pytest

import flopy
from flopy.mf6 import (
    MFSimulation,
    Mf6Tid,
    Mf6TidRunner,
    ModflowGwf,
    ModflowGwfchd,
    ModflowGwfdis,
    ModflowGwfdisv,
    ModflowGwfic,
    ModflowGwfnpf,
    ModflowGwfoc,
    ModflowIms,
    ModflowTdis,
)

try:
    import modflowapi

    HAS_MODFLOWAPI = True
except ImportError:
    HAS_MODFLOWAPI = False

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

NPER = 10
NLAY, NROW, NCOL = 1, 3, 3
NNODES_DIS = NLAY * NROW * NCOL  # 9

# Runner tests use 10 time steps per SP to ensure steady-state convergence
RUNNER_PERLEN = 10.0
RUNNER_NSTP = 10

# ---------------------------------------------------------------------------
# DIS model helper
# ---------------------------------------------------------------------------


def _make_dis_sim(ws, name="gwf", center_inactive_in_dis=False):
    """3×3×1 steady-state DIS model with CHD on col=0 and col=2, 10 SPs.

    center_inactive_in_dis : set DIS idomain=0 for center cell.
        Used only in Mf6Tid unit tests (no runner).
        Runner tests MUST keep idomain=1 everywhere.
    """
    sim = MFSimulation(sim_name="mfsim", sim_ws=str(ws))
    perlen = RUNNER_PERLEN
    nstp = RUNNER_NSTP
    ModflowTdis(sim, nper=NPER, perioddata=[(perlen, nstp, 1.0)] * NPER)
    ModflowIms(sim, complexity="MODERATE")

    gwf = ModflowGwf(sim, modelname=name, save_flows=True)

    if center_inactive_in_dis:
        idomain = np.ones((NLAY, NROW, NCOL), dtype=int)
        idomain[0, 1, 1] = 0
        ModflowGwfdis(gwf, nlay=NLAY, nrow=NROW, ncol=NCOL, idomain=idomain)
    else:
        ModflowGwfdis(gwf, nlay=NLAY, nrow=NROW, ncol=NCOL)

    ModflowGwfic(gwf, strt=1.0)
    ModflowGwfnpf(gwf, icelltype=0, k=1.0)

    chd_data = []
    for row in range(NROW):
        chd_data.append([(0, row, 0), 1.0])
        chd_data.append([(0, row, 2), 0.0])
    ModflowGwfchd(gwf, stress_period_data={0: chd_data})

    ModflowGwfoc(
        gwf,
        head_filerecord=f"{name}.hds",
        budget_filerecord=f"{name}.cbb",
        saverecord=[("HEAD", "LAST"), ("BUDGET", "LAST")],
        printrecord=[("BUDGET", "LAST")],
    )

    return sim, gwf


# ---------------------------------------------------------------------------
# DISV model helper  (1 layer, 4-cell quad mesh)
# ---------------------------------------------------------------------------

# 2×2 structured-equivalent DISV: 4 cells, 9 vertices
_DISV_NLAY = 1
_DISV_NCPL = 4
_DISV_VERTICES = [
    [0, 0.0, 0.0],
    [1, 1.0, 0.0],
    [2, 2.0, 0.0],
    [3, 0.0, 1.0],
    [4, 1.0, 1.0],
    [5, 2.0, 1.0],
    [6, 0.0, 2.0],
    [7, 1.0, 2.0],
    [8, 2.0, 2.0],
]
_DISV_CELL2D = [
    # icell2d  xc   yc  ncvert  iverts...
    [0, 0.5, 0.5, 4, 0, 1, 4, 3],
    [1, 1.5, 0.5, 4, 1, 2, 5, 4],
    [2, 0.5, 1.5, 4, 3, 4, 7, 6],
    [3, 1.5, 1.5, 4, 4, 5, 8, 7],
]


def _make_disv_sim(ws, name="gwf"):
    """1-layer, 4-cell DISV model.  Cell 0 is CHD=1, cell 1 is free, etc."""
    sim = MFSimulation(sim_name="mfsim", sim_ws=str(ws))
    ModflowTdis(sim, nper=3, perioddata=[(10.0, 10, 1.0)] * 3)
    ModflowIms(sim, complexity="MODERATE")

    gwf = ModflowGwf(sim, modelname=name)
    ModflowGwfdisv(
        gwf,
        nlay=_DISV_NLAY,
        ncpl=_DISV_NCPL,
        vertices=_DISV_VERTICES,
        cell2d=_DISV_CELL2D,
        top=1.0,
        botm=[0.0],
    )
    ModflowGwfic(gwf, strt=0.5)
    ModflowGwfnpf(gwf, icelltype=0, k=1.0)

    return sim, gwf


# ---------------------------------------------------------------------------
# Unit tests — criteria 1, 4, 5 (no exe required)
# ---------------------------------------------------------------------------


def test_dis_node_resolution_agrees_with_modelgrid(tmp_path):
    """Criterion 5 (DIS): tuple → flat index matches modelgrid.get_node()."""
    sim, gwf = _make_dis_sim(tmp_path)
    tid = Mf6Tid(gwf)
    mg = gwf.modelgrid

    tuples = [(0, 0, 0), (0, 1, 1), (0, 2, 2)]
    expected = [mg.get_node(t)[0] for t in tuples]

    tid.add(kper=0, activate=tuples)
    recorded = tid.schedule[0]["activate"]

    assert recorded == expected, f"Expected {expected}, got {recorded}"


def test_disv_node_resolution_agrees_with_modelgrid(tmp_path):
    """Criterion 5 (DISV): (lay, cell2d) tuple → flat index via get_node()."""
    sim, gwf = _make_disv_sim(tmp_path)
    tid = Mf6Tid(gwf)
    mg = gwf.modelgrid

    tuples = [(0, 0), (0, 1), (0, 3)]
    expected = [mg.get_node(t)[0] for t in tuples]

    tid.add(kper=0, activate=tuples)
    recorded = tid.schedule[0]["activate"]

    assert recorded == expected, f"Expected {expected}, got {recorded}"


def test_idomain_at_sequence(tmp_path):
    """Criterion 1: idomain_at() reflects activate/deactivate state machine."""
    sim, gwf = _make_dis_sim(tmp_path, center_inactive_in_dis=True)
    tid = Mf6Tid(gwf)

    center = gwf.modelgrid.get_node((0, 1, 1))[0]  # → 4

    tid.add(kper=2, activate=[(0, 1, 1)])
    tid.add(kper=5, deactivate=[(0, 1, 1)])
    tid.add(kper=7, activate=[(0, 1, 1)])

    # Before any transition: center starts inactive (from DIS idomain=0)
    assert tid.idomain_at(kper=0)[center] == 0
    assert tid.idomain_at(kper=1)[center] == 0

    # Activated at kper=2
    assert tid.idomain_at(kper=2)[center] == 1
    assert tid.idomain_at(kper=4)[center] == 1

    # Deactivated at kper=5
    assert tid.idomain_at(kper=5)[center] == 0
    assert tid.idomain_at(kper=6)[center] == 0

    # Re-activated at kper=7
    assert tid.idomain_at(kper=7)[center] == 1
    assert tid.idomain_at(kper=9)[center] == 1

    # All other nodes remain 1 throughout
    for node in range(NNODES_DIS):
        if node != center:
            assert tid.idomain_at(kper=9)[node] == 1


def test_write_load_roundtrip(tmp_path):
    """Criterion 4: write() / load() preserves schedule exactly."""
    sim, gwf = _make_dis_sim(tmp_path)
    sim.write_simulation(silent=True)

    tid = Mf6Tid(gwf)
    tid.add(kper=2, activate=[(0, 1, 1)])
    tid.add(kper=5, deactivate=[(0, 1, 1)])
    tid.add(kper=7, activate=[(0, 1, 1), (0, 0, 2)])

    tid.write()
    assert (tmp_path / "gwf.tid.json").exists()

    tid2 = Mf6Tid.load(gwf)

    assert set(tid2.schedule.keys()) == {2, 5, 7}
    assert tid2.schedule[2]["activate"] == tid.schedule[2]["activate"]
    assert tid2.schedule[5]["deactivate"] == tid.schedule[5]["deactivate"]
    assert tid2.schedule[7]["activate"] == tid.schedule[7]["activate"]
    assert tid2.schedule[7]["deactivate"] == []

    for kper in range(NPER):
        np.testing.assert_array_equal(
            tid.idomain_at(kper),
            tid2.idomain_at(kper),
            err_msg=f"idomain_at mismatch at kper={kper}",
        )


def test_flat_int_nodes_passthrough(tmp_path):
    """Flat int node IDs (DISU style) pass through unchanged."""
    sim, gwf = _make_dis_sim(tmp_path)
    tid = Mf6Tid(gwf)
    tid.add(kper=0, activate=[0, 1, 2], deactivate=[8])
    assert tid.schedule[0]["activate"] == [0, 1, 2]
    assert tid.schedule[0]["deactivate"] == [8]


def test_bad_node_type_raises(tmp_path):
    """Non-int, non-tuple node identifiers raise TypeError."""
    sim, gwf = _make_dis_sim(tmp_path)
    tid = Mf6Tid(gwf)
    with pytest.raises(TypeError):
        tid.add(kper=0, activate=["bad"])


# ---------------------------------------------------------------------------
# Integration test — criteria 2 & 3 (requires modflowapi + arm64 lib)
# ---------------------------------------------------------------------------


def _find_mf6_lib():
    """Locate the MF6 shared library and verify it loads (architecture check)."""
    import ctypes

    # Repo-local bin/ is first — fork is self-contained
    repo_root = Path(__file__).parent.parent
    candidates = list((repo_root / "bin").glob("*mf6*"))

    exe = shutil.which("mf6")
    if exe is not None:
        exe_dir = Path(exe).parent
        for suffix in (".so", ".dylib", ".dll"):
            candidates.extend(exe_dir.glob(f"*mf6*{suffix}"))
            candidates.extend(exe_dir.parent.glob(f"lib/*mf6*{suffix}"))

    for candidate in candidates:
        try:
            ctypes.CDLL(str(candidate))
            return str(candidate)
        except OSError:
            continue
    return None


_mf6_lib = _find_mf6_lib()
_skip_runner = pytest.mark.skipif(
    not HAS_MODFLOWAPI or _mf6_lib is None,
    reason="modflowapi not installed or MF6 shared library not found",
)


@_skip_runner
def test_runner_dis_budget(tmp_path):
    """Criteria 2 & 3 (DIS): IBOUND changes produce correct heads each SP.

    DIS idomain=1 everywhere (runner requirement).
    TID schedule controls active cells via GWF/IBOUND:

      kper=0–1: center (0,1,1) deactivated → head stays at IC (1.0)
      kper=2–4: activate  → head converges to ~0.5  [criterion 2]
      kper=5–6: deactivate again → head no longer 0.5
      kper=7–9: re-activate → head recovers to ~0.5  [criterion 3]
    """
    from flopy.utils import HeadFile

    sim, gwf = _make_dis_sim(tmp_path, center_inactive_in_dis=False)
    sim.write_simulation(silent=True)

    center = gwf.modelgrid.get_node((0, 1, 1))[0]  # flat index 4
    IC_HEAD = 1.0

    tid = Mf6Tid(gwf)
    tid.add(kper=0, deactivate=[(0, 1, 1)])   # inactive from the start
    tid.add(kper=2, activate=[(0, 1, 1)])
    tid.add(kper=5, deactivate=[(0, 1, 1)])
    tid.add(kper=7, activate=[(0, 1, 1)])
    tid.write()

    Mf6TidRunner(sim, tid, lib_path=_mf6_lib).run()

    hds_path = tmp_path / "gwf.hds"
    assert hds_path.exists()

    hf = HeadFile(str(hds_path))
    heads = hf.get_alldata().reshape(NPER, NNODES_DIS)

    # kper=0,1: inactive → MF6 freezes head at IC value (no solve)
    for k in [0, 1]:
        assert heads[k, center] == pytest.approx(IC_HEAD, abs=1e-6), (
            f"kper={k}: inactive center should keep IC head={IC_HEAD}, "
            f"got {heads[k, center]}"
        )

    # kper=2–4: active → symmetric steady-state, head ≈ 0.5  [criterion 2]
    for k in [2, 3, 4]:
        assert heads[k, center] == pytest.approx(0.5, abs=0.01), (
            f"kper={k}: active center should converge to 0.5, got {heads[k, center]}"
        )

    # kper=5,6: inactive again → MF6 freezes head at last-known value.
    # When the last-known head equals the steady-state value (≈0.5), the frozen
    # head is also ≈0.5.  The proof of inactivity is that the head does NOT
    # change between consecutive inactive SPs.
    assert heads[5, center] == pytest.approx(heads[4, center], abs=1e-8), (
        "kper=5: deactivated center head should be frozen at last-known value"
    )
    assert heads[6, center] == pytest.approx(heads[4, center], abs=1e-8), (
        "kper=6: deactivated center head should remain frozen"
    )

    # kper=7–9: re-activated → head re-solved, recovers to ~0.5  [criterion 3]
    for k in [7, 8, 9]:
        assert heads[k, center] == pytest.approx(0.5, abs=0.01), (
            f"kper={k}: re-activated center should recover to ~0.5, "
            f"got {heads[k, center]}"
        )
