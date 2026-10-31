"""
Optional end-to-end USG-Transport 2.7 executable validation (Stage 3 Card 9).

These tests prove that FloPy-authored-from-scratch USG-T input not only writes
but actually *runs* under the USG-Transport executable and produces stable
outputs. They are intentionally kept out of the default focused suite
(`autotest/test_usg_transport.py`, which stays light) and only run when the
executable is available.

Executable selection (single contract for the whole executable tier):

- Set ``USGT_EXE`` to a **USG-Transport 2.7** executable (a name on ``PATH`` or
  an absolute path). If unset, it defaults to ``mfusg_gsi``.
- ``@requires_exe(USGT_EXE)`` skips every test here cleanly when the executable
  cannot be resolved, so the default suite is unaffected.

The real-model ``Ex1..Ex9`` run tests in `test_usg_transport.py` use the **same**
``USGT_EXE`` contract (they ``load + write + run`` real USG-T models, including
the Ex7 multi-species transport model). This file adds the complementary
*from-scratch authoring → execution* direction.
"""

import os
import re
import shutil
from pathlib import Path

import numpy as np
import pytest
from modflow_devtools.markers import requires_exe, requires_pkg

from flopy.mfusg import (
    MfUsg,
    MfUsgBas,
    MfUsgBcf,
    MfUsgBct,
    MfUsgDpf,
    MfUsgEvt,
    MfUsgLpf,
    MfUsgOc,
    MfUsgPcb,
    MfUsgRch,
    MfUsgSms,
    MfUsgTib,
)
from flopy.modflow import ModflowChd, ModflowDis, ModflowRch
from flopy.utils import (
    CellBudgetFile,
    HeadFile,
    HeadUFile,
    MfusgListBudget,
    MfusgTransportListBudget,
)

USGT_EXE = os.environ.get("USGT_EXE", "mfusg_gsi")


@requires_exe(USGT_EXE)
def test_usgt_exe_minimal_flow_from_scratch(function_tmpdir):
    """A 1-D CHD-driven steady flow model authored from scratch runs and
    reproduces the analytical linear head gradient with a closed budget."""
    ml = MfUsg(
        modelname="flow",
        model_ws=str(function_tmpdir),
        exe_name=USGT_EXE,
        structured=True,
    )
    ModflowDis(
        ml, nlay=1, nrow=1, ncol=5, nper=1, perlen=1.0, nstp=1, steady=True,
        delr=10.0, delc=10.0, top=10.0, botm=0.0,
    )
    MfUsgBas(ml, ibound=1, strt=5.0)
    MfUsgLpf(ml, laytyp=0, hk=1.0, ipakcb=0)
    MfUsgSms(ml, linmeth=1)  # PCGU; linmeth=2/XMD does not converge this trivial system
    MfUsgOc(ml, stress_period_data={(0, 0): ["save head"]})
    ModflowChd(ml, stress_period_data={0: [[0, 0, 0, 8.0, 8.0],
                                           [0, 0, 4, 2.0, 2.0]]})
    ml.write_input()

    success, _ = ml.run_model(silent=True)
    assert success, "USG-T run did not terminate normally"

    heads = HeadFile(os.path.join(ml.model_ws, "flow.hds")).get_data().ravel()
    # CHD at cells 0 and 4 (8.0 and 2.0) -> exact linear interior [6.5, 5.0, 3.5]
    assert np.allclose(heads, [8.0, 6.5, 5.0, 3.5, 2.0], atol=1e-3)

    inc, _cum = MfusgListBudget(
        os.path.join(ml.model_ws, "flow.list")
    ).get_budget()
    assert abs(inc["PERCENT_DISCREPANCY"][-1]) < 0.1


@requires_exe(USGT_EXE)
def test_usgt_exe_minimal_transport_from_scratch(function_tmpdir):
    """A from-scratch BCT transport model with a PCB concentration source runs,
    produces a concentration (`.con`) output file, and closes the
    species-isolated transport mass budget.

    (Concentration *values* are not asserted here — the `.con` is USG node-based
    binary; per-value comparison is covered by the real Ex transport models.)"""
    ml = MfUsg(
        modelname="tran",
        model_ws=str(function_tmpdir),
        exe_name=USGT_EXE,
        structured=True,
    )
    ModflowDis(
        ml, nlay=1, nrow=1, ncol=5, nper=1, perlen=100.0, nstp=10, steady=False,
        delr=10.0, delc=10.0, top=10.0, botm=0.0,
    )
    MfUsgBas(ml, ibound=1, strt=5.0)
    MfUsgLpf(ml, laytyp=0, hk=1.0, ss=1.0e-5, ipakcb=0)
    MfUsgSms(ml, linmeth=1)
    MfUsgBct(ml, itrnsp=1, mcomp=1, prsity=0.2, conc=0.0, idisp=0)
    # Prescribed concentration = 1.0 at the inflow cell (k,i,j,iSpec,conc)
    MfUsgPcb(ml, stress_period_data={0: [[0, 0, 0, 1, 1.0]]})
    MfUsgOc(ml, stress_period_data={(0, 0): ["save head", "save concentration"]})
    ModflowChd(ml, stress_period_data={0: [[0, 0, 0, 8.0, 8.0],
                                           [0, 0, 4, 2.0, 2.0]]})
    ml.write_input()

    success, _ = ml.run_model(silent=True)
    assert success, "USG-T transport run did not terminate normally"

    # Concentration output is produced (USG node-based binary; value comparison
    # is covered by the real Ex models).
    assert os.path.isfile(os.path.join(ml.model_ws, "tran.con"))

    # The transport mass budget closes and is species-isolated (text listing,
    # the reliable cross-check). The PCB source shows up as PRESCRIBED_CONCS_IN.
    inc, _cum = MfusgTransportListBudget(
        os.path.join(ml.model_ws, "tran.list"), species=1
    ).get_budget()
    assert inc is not None and len(inc) > 0
    assert "PRESCRIBED_CONCS_IN" in inc.dtype.names
    assert abs(inc["PERCENT_DISCREPANCY"][-1]) < 0.1


@requires_exe(USGT_EXE)
def test_usgt_exe_tib_prescribed_head_from_scratch(function_tmpdir):
    """A from-scratch semantic TIB package runs under USG-T 2.7.

    A 2-period steady model: period 1 is an unconstrained linear gradient
    between two CHD cells; period 2 uses TIB ``NIBM1`` + ``HEAD`` to prescribe an
    interior head, bending the gradient to a different exact solution. This
    proves the executable reads and applies the FloPy-authored semantic TIB
    record (not just that the file is syntactically accepted)."""
    ml = MfUsg(
        modelname="tibflow",
        model_ws=str(function_tmpdir),
        exe_name=USGT_EXE,
        structured=True,
    )
    ModflowDis(
        ml, nlay=1, nrow=1, ncol=5, nper=2, perlen=1.0, nstp=1, steady=True,
        delr=10.0, delc=10.0, top=10.0, botm=0.0,
    )
    MfUsgBas(ml, ibound=1, strt=5.0)
    MfUsgLpf(ml, laytyp=0, hk=1.0, ipakcb=0)
    MfUsgSms(ml, linmeth=1)
    MfUsgOc(
        ml, stress_period_data={(0, 0): ["save head"], (1, 0): ["save head"]}
    )
    chd = [[0, 0, 0, 8.0, 8.0], [0, 0, 4, 2.0, 2.0]]
    ModflowChd(ml, stress_period_data={0: chd, 1: chd})
    # Period 2: prescribe interior node 2 (0-based) at 6.0 via NIBM1 + HEAD.
    MfUsgTib(ml, stress_period_data={1: {"ibm1": [(2, 6.0)]}})
    ml.write_input()

    success, _ = ml.run_model(silent=True)
    assert success, "USG-T TIB run did not terminate normally"

    heads = HeadFile(os.path.join(ml.model_ws, "tibflow.hds")).get_alldata()
    # Period 1: unconstrained linear gradient between the CHD cells.
    assert np.allclose(heads[0].ravel(), [8.0, 6.5, 5.0, 3.5, 2.0], atol=1e-3)
    # Period 2: TIB pins node 2 at 6.0 -> exact piecewise-linear [8,7,6,4,2].
    assert np.allclose(heads[-1].ravel(), [8.0, 7.0, 6.0, 4.0, 2.0], atol=1e-3)

    inc, _cum = MfusgListBudget(
        os.path.join(ml.model_ws, "tibflow.list")
    ).get_budget()
    assert abs(inc["PERCENT_DISCREPANCY"][-1]) < 0.1


@requires_exe(USGT_EXE)
def test_usgt_exe_evt_from_scratch(function_tmpdir):
    """A from-scratch EVT package runs under USG-T 2.7.

    The same steady CHD flow model plus a NEVTOP=1 EVT (surface above the head
    range, finite extinction depth) terminates normally, reports an ET term in
    the volumetric budget, and closes — proving the executable reads and applies
    the FloPy-authored EVT input."""
    ml = MfUsg(
        modelname="evtflow",
        model_ws=str(function_tmpdir),
        exe_name=USGT_EXE,
        structured=True,
    )
    ModflowDis(
        ml, nlay=1, nrow=1, ncol=5, nper=1, perlen=1.0, nstp=1, steady=True,
        delr=10.0, delc=10.0, top=10.0, botm=0.0,
    )
    MfUsgBas(ml, ibound=1, strt=5.0)
    MfUsgLpf(ml, laytyp=0, hk=1.0, ipakcb=0)
    MfUsgSms(ml, linmeth=1)
    MfUsgOc(ml, stress_period_data={(0, 0): ["save head", "print budget"]})
    ModflowChd(ml, stress_period_data={0: [[0, 0, 0, 8.0, 8.0],
                                           [0, 0, 4, 2.0, 2.0]]})
    MfUsgEvt(ml, nevtop=1, ipakcb=0, surf=10.0, exdp=8.0, evtr=1.0e-3)
    ml.write_input()

    success, _ = ml.run_model(silent=True)
    assert success, "USG-T EVT run did not terminate normally"

    inc, _cum = MfusgListBudget(
        os.path.join(ml.model_ws, "evtflow.list")
    ).get_budget()
    # ET is reported in the budget (ET removes water -> ET_OUT) and it closes.
    assert any("ET" in name for name in inc.dtype.names)
    assert abs(inc["PERCENT_DISCREPANCY"][-1]) < 0.1


@requires_exe(USGT_EXE)
def test_usgt_exe_oc_atsa_from_scratch(function_tmpdir):
    """OC adaptive time-stepping (ATSA) runs under USG-T 2.7. A transient model
    with ``MfUsgOc(atsa=1)`` exercises the ATS path: the Fortran overrides NSTP
    and uses its built-in DELTAT/TMINAT/TMAXAT/TADJAT/TCUTAT (glo2basu1.f), so the
    OC file only needs the ``ATSA`` keyword. With the CHD boundaries held fixed
    the run reaches the analytical steady gradient.

    (FASTFORWARD/FASTFORWARDC/BOOTSTRAPPING are not exe-smoked here — they read
    GWF/CLN/DDF heads (or concentrations) from an external file produced by a
    prior run, i.e. a two-stage fixture; see USGT_STAGE4_18_EXE_VERIFICATION.md
    for the manual tier.)"""
    ml = MfUsg(
        modelname="atsa",
        model_ws=str(function_tmpdir),
        exe_name=USGT_EXE,
        structured=True,
    )
    ModflowDis(
        ml, nlay=1, nrow=1, ncol=5, nper=1, perlen=100.0, nstp=1, steady=False,
        delr=10.0, delc=10.0, top=10.0, botm=0.0,
    )
    MfUsgBas(ml, ibound=1, strt=5.0)
    MfUsgLpf(ml, laytyp=0, hk=1.0, ss=1.0e-5, ipakcb=0)
    MfUsgSms(ml, linmeth=1)
    MfUsgOc(ml, atsa=1, stress_period_data={(0, 0): ["save head"]})
    ModflowChd(ml, stress_period_data={0: [[0, 0, 0, 8.0, 8.0],
                                           [0, 0, 4, 2.0, 2.0]]})
    ml.write_input()
    assert "ATSA" in (function_tmpdir / "atsa.oc").read_text()

    success, _ = ml.run_model(silent=True)
    assert success, "USG-T ATSA run did not terminate normally"

    # transient model with fixed CHD reaches the analytical linear gradient
    heads = HeadFile(os.path.join(ml.model_ws, "atsa.hds")).get_data().ravel()
    assert np.allclose(heads, [8.0, 6.5, 5.0, 3.5, 2.0], atol=1e-3)


@requires_exe(USGT_EXE)
def test_usgt_exe_dpf_layer_order_matches_fortran_contract(function_tmpdir):
    """A two-layer DPF written by FloPy (default DPF unit) runs exactly like a
    hand-written DPF in the GWF2DPFU1AR order (PHIF, DDFTR, SC1IM per layer,
    gwf2dpf1u1.f DO 200). Confined LPF layers, closed system: the mobile heads
    respond only to the dual-domain exchange, so a layer mix-up changes them."""
    phif, ddftr, sc1im = (0.20, 0.05), (1.0e-1, 1.0e-4), (1.0e-4, 1.0e-3)
    hnewim = (9.0, 7.0)

    def build(ws):
        ml = MfUsg(
            modelname="dpf", model_ws=str(ws), exe_name=USGT_EXE, structured=True
        )
        ModflowDis(
            ml,
            nlay=2,
            nrow=1,
            ncol=3,
            nper=1,
            perlen=10.0,
            nstp=10,
            steady=False,
            delr=10.0,
            delc=10.0,
            top=10.0,
            botm=[5.0, 0.0],
        )
        MfUsgBas(ml, ibound=1, strt=8.0)
        MfUsgLpf(ml, laytyp=0, hk=1.0, vka=1.0, ss=1.0e-4, ipakcb=0)
        MfUsgSms(ml, linmeth=1)
        MfUsgOc(ml, stress_period_data={(0, k): ["save head"] for k in range(10)})
        return ml

    written = build(function_tmpdir / "written")
    MfUsgDpf(
        written,
        iboundim=1,
        hnewim=[[[h] * 3] for h in hnewim],
        phif=[[[v] * 3] for v in phif],
        ddftr=[[[v] * 3] for v in ddftr],
        sc1im=[[[v] * 3] for v in sc1im],
    )
    written.write_input()

    manual = build(function_tmpdir / "manual")
    (function_tmpdir / "manual").mkdir(exist_ok=True)
    (function_tmpdir / "manual" / "manual.dpf").write_text(
        "         0         0         0\n"
        "CONSTANT          1 #iboundim layer 1\n"
        "CONSTANT          1 #iboundim layer 2\n"
        "CONSTANT        9.0 #hnewim layer 1\n"
        "CONSTANT        7.0 #hnewim layer 2\n"
        "CONSTANT        0.2 #phif layer 1\n"
        "CONSTANT        0.1 #ddftr layer 1\n"
        "CONSTANT     1.0E-4 #sc1im layer 1\n"
        "CONSTANT       0.05 #phif layer 2\n"
        "CONSTANT     1.0E-4 #ddftr layer 2\n"
        "CONSTANT     1.0E-3 #sc1im layer 2\n"
    )
    # the package only registers the hand-written file in the name file
    MfUsgDpf(manual, filenames="manual.dpf")
    manual.write_input(SelPackList=["DIS", "BAS6", "LPF", "SMS", "OC"])
    assert "manual.dpf" in (function_tmpdir / "manual" / "dpf.nam").read_text()

    results = []
    for ml in (written, manual):
        success, buff = ml.run_model(silent=True, report=True)
        assert success, f"USG-T DPF run in {ml.model_ws} did not terminate normally"
        assert any("normal termination" in line.lower() for line in buff)
        hds = HeadFile(os.path.join(ml.model_ws, "dpf.hds"))
        results.append((hds.get_times(), hds.get_alldata()))

    (times_w, heads_w), (times_m, heads_m) = results
    assert heads_w.shape == heads_m.shape == (10, 2, 1, 3)
    assert times_w == times_m
    assert np.all(np.isfinite(heads_w)) and np.all(np.isfinite(heads_m))
    # the exchange is active: mobile heads move away from the initial 8.0
    assert np.abs(heads_m[-1] - 8.0).max() > 0.5
    np.testing.assert_allclose(heads_w, heads_m, rtol=0, atol=1e-6)


# LPF read-order checks (READLPF/SGWF2LPFU1S/G): the LPF file written by FloPy
# must run exactly like a hand-written file in the Fortran order. Forcing is
# K-sensitive (CHD top-left 26 and bottom-right 24, recharge on top) and each
# case has a perturbed hand-written variant as a sensitivity control. SMS is
# Newton/XMD without option keywords (the SMS options line is out of scope).
_LPF_EXE_CASES = {
    "f06_mixed_richards": (
        {"laytyp": [0, 5], "richards": True, "transient": False},
        ("hk", 2, 0.5),
    ),
    "f08_steady_vkcb_wetdry": (
        {"laytyp": [1, 0], "laywet": [1, 0], "laycbd": [1, 0], "transient": False},
        ("vkcb", 1, 10.0),
    ),
    "sy_thickstrt": (
        {"laytyp": [-1, 1], "thickstrt": True, "transient": True},
        ("hk", 2, 0.5),
    ),
    "f09_disu_ikcflag_m1": (
        {"laytyp": [0], "disu": True, "ikcflag": -1, "transient": False},
        ("ksat", 0, 0.5),
    ),
}


def _lpf_exe_model(ws, case, flow_file=None):
    from autotest.test_usg_transport import _LPF_KSAT, _LPF_VALUES, _lpf_disu_line

    nlay, disu = len(case["laytyp"]), case.get("disu", False)
    steady = not case.get("transient", True)
    ml = MfUsg(
        modelname="lpf", model_ws=str(ws), exe_name=USGT_EXE, structured=not disu
    )
    if disu:
        _lpf_disu_line(ml, nlay, steady)
        ncol = 3
    else:
        laycbd = case.get("laycbd", [0] * nlay)
        ModflowDis(
            ml,
            nlay=nlay,
            nrow=1,
            ncol=3,
            nper=1,
            nstp=2 if not steady else 1,
            steady=steady,
            delr=10.0,
            delc=10.0,
            top=30.0,
            botm=[30.0 - 10.0 * (i + 1) for i in range(nlay + sum(laycbd))],
            laycbd=laycbd,
        )
        ncol = 3
    MfUsgBas(ml, ibound=1, strt=25.0, richards=case.get("richards", False))
    values = {
        n: [np.full(ncol if disu else (1, ncol), v[k]) for k in range(nlay)]
        for n, v in _LPF_VALUES.items()
    }
    lpf = MfUsgLpf(
        ml,
        laytyp=case["laytyp"],
        laywet=case.get("laywet", 0),
        ikcflag=case.get("ikcflag", 0),
        thickstrt=case.get("thickstrt", False),
        ksat=_LPF_KSAT if disu else 1.0,
        **values,
    )
    MfUsgSms(
        ml,
        hclose=1e-5,
        hiclose=1e-7,
        mxiter=250,
        iter1=600,
        iprsms=1,
        nonlinmeth=1,
        linmeth=1,
        theta=0.7,
        akappa=0.07,
        gamma=0.1,
        amomentum=0.0,
        numtrack=200,
        btol=1.1,
        breduc=0.2,
        reslim=1.0,
        iacl=1,
        norder=0,
        level=7,
        north=14,
        iredsys=0,
        rrctol=0.0,
        idroptol=1,
        epsrn=1e-3,
    )
    if disu:
        chd = [[0, 26.0, 26.0], [(nlay - 1) * ncol + ncol - 1, 24.0, 24.0]]
        ModflowRch(ml, rech={0: np.full(ncol, 0.05, dtype=np.float32)}, nrchop=1)
    else:
        chd = [[0, 0, 0, 26.0, 26.0], [nlay - 1, 0, ncol - 1, 24.0, 24.0]]
        ModflowRch(ml, rech=0.05, nrchop=1)
    ModflowChd(ml, stress_period_data={0: chd})
    nstp = 1 if steady else 2
    MfUsgOc(ml, stress_period_data={(0, k): ["save head"] for k in range(nstp)})
    if flow_file is None:
        ml.write_input()
    else:
        (ws / "manual.lpf").write_text(flow_file)
        lpf.file_name[0] = "manual.lpf"
        ml.write_input(SelPackList=[p.name[0] for p in ml.packagelist if p is not lpf])
    return ml


def _lpf_exe_heads(ml, disu):
    success, buff = ml.run_model(silent=True, report=True)
    assert success, f"USG-T run in {ml.model_ws} did not terminate normally"
    assert any("normal termination" in line.lower() for line in buff)
    path = os.path.join(ml.model_ws, "lpf.hds")
    if disu:
        hds = HeadUFile(path)
        times = hds.get_times()
        return times, np.array([np.concatenate(hds.get_data(totim=t)) for t in times])
    hds = HeadFile(path)
    return hds.get_times(), hds.get_alldata()


@requires_exe(USGT_EXE)
@pytest.mark.parametrize("name", list(_LPF_EXE_CASES))
def test_usgt_exe_lpf_written_matches_fortran_order(function_tmpdir, name):
    """FloPy-written LPF = hand-written LPF in the READLPF order (all cells and
    time steps; identical records, so atol=1e-6 only absorbs formatting), and
    a perturbed hand-written variant shows the heads depend on the records the
    defect used to shift or omit."""
    from autotest.test_usg_transport import _LPF_KSAT, _LPF_VALUES, _lpf_manual_text

    case, (prop, layer, factor) = _LPF_EXE_CASES[name]
    disu = case.get("disu", False)
    manual_text = _lpf_manual_text(case)
    if prop == "ksat":
        old = " ".join(f"{v:.2f}" for v in _LPF_KSAT)
        new = " ".join(f"{v * factor:.4f}" for v in _LPF_KSAT)
    else:
        v = _LPF_VALUES[prop][layer - 1]
        old, new = (
            f"CONSTANT {v:.6E} #{prop} layer {layer}",
            f"CONSTANT {v * factor:.6E} #{prop} layer {layer}",
        )
    assert manual_text.count(old) == 1
    perturbed_text = manual_text.replace(old, new)

    t_w, h_w = _lpf_exe_heads(_lpf_exe_model(function_tmpdir / "written", case), disu)
    t_m, h_m = _lpf_exe_heads(
        _lpf_exe_model(function_tmpdir / "manual", case, manual_text), disu
    )
    t_p, h_p = _lpf_exe_heads(
        _lpf_exe_model(function_tmpdir / "perturbed", case, perturbed_text), disu
    )

    assert h_w.shape == h_m.shape == h_p.shape
    assert t_w == t_m == t_p
    assert np.all(np.isfinite(h_w)) and np.all(np.isfinite(h_m))
    np.testing.assert_allclose(h_w, h_m, rtol=0, atol=1e-6)
    # sensitivity control: the perturbed record moves the heads well beyond atol
    assert np.abs(h_p - h_m).max() > 1e-3


# BCF read-order checks (READBCF/SGWF2BCFU1S/G): the BCF file written by FloPy
# must run exactly like a hand-written file in the Fortran order; a perturbed
# hand-written variant is the sensitivity control. Same K-sensitive forcing
# and SMS settings as the LPF checks.
_BCF_EXE_CASES = {
    "b01_intercellt": (
        {"laycon": [3, 0], "intercellt": [1, 0], "transient": False},
        ("hy", 1, 0.5),
    ),
    "b05_richards_steady": (
        {"laycon": [5, 5], "richards": True, "transient": False},
        ("hy", 2, 0.5),
    ),
    "b04_b05_richards_transient": (
        {"laycon": [3, 5], "richards": True, "transient": True},
        ("hy", 2, 0.5),
    ),
    "b08_hy_after_promotion": (
        {"laycon": [0, 2], "richards": True, "transient": False},
        ("hy", 2, 0.5),
    ),
    "b02_f35_disu_ikcflag_m1": (
        {"laycon": [0], "disu": True, "ikcflag": -1, "transient": False},
        ("ksat", 0, 0.5),
    ),
}


def _bcf_exe_model(ws, case, flow_file=None):
    from autotest.test_usg_transport import _BCF_KSAT, _BCF_VALUES, _lpf_disu_line

    nlay, disu = len(case["laycon"]), case.get("disu", False)
    steady = not case.get("transient", True)
    ml = MfUsg(
        modelname="bcf", model_ws=str(ws), exe_name=USGT_EXE, structured=not disu
    )
    ncol = 3
    if disu:
        _lpf_disu_line(ml, nlay, steady)
    else:
        ModflowDis(
            ml,
            nlay=nlay,
            nrow=1,
            ncol=ncol,
            nper=1,
            nstp=1 if steady else 2,
            steady=steady,
            delr=10.0,
            delc=10.0,
            top=30.0,
            botm=[30.0 - 10.0 * (k + 1) for k in range(nlay)],
        )
    MfUsgBas(ml, ibound=1, strt=25.0, richards=case.get("richards", False))
    values = {
        n: [np.full(ncol if disu else (1, ncol), v[k]) for k in range(nlay)]
        for n, v in _BCF_VALUES.items()
    }
    values["vcont"] = values["vcont"][: max(nlay - 1, 1)]
    if not case.get("richards", False) and 5 not in case["laycon"]:
        for name in ("alpha", "beta", "sr", "brook", "bp"):
            del values[name]  # optional arguments, only for Richards layers
    bcf = MfUsgBcf(
        ml,
        laycon=case["laycon"],
        intercellt=case.get("intercellt", 0),
        ikcflag=case.get("ikcflag", 0),
        ksat=_BCF_KSAT if disu else 1.0,
        **values,
    )
    MfUsgSms(
        ml,
        hclose=1e-5,
        hiclose=1e-7,
        mxiter=250,
        iter1=600,
        iprsms=1,
        nonlinmeth=1,
        linmeth=1,
        theta=0.7,
        akappa=0.07,
        gamma=0.1,
        amomentum=0.0,
        numtrack=200,
        btol=1.1,
        breduc=0.2,
        reslim=1.0,
        iacl=1,
        norder=0,
        level=7,
        north=14,
        iredsys=0,
        rrctol=0.0,
        idroptol=1,
        epsrn=1e-3,
    )
    if disu:
        chd = [[0, 26.0, 26.0], [(nlay - 1) * ncol + ncol - 1, 24.0, 24.0]]
        ModflowRch(ml, rech={0: np.full(ncol, 0.05, dtype=np.float32)}, nrchop=1)
    else:
        chd = [[0, 0, 0, 26.0, 26.0], [nlay - 1, 0, ncol - 1, 24.0, 24.0]]
        ModflowRch(ml, rech=0.05, nrchop=1)
    ModflowChd(ml, stress_period_data={0: chd})
    nstp = 1 if steady else 2
    MfUsgOc(ml, stress_period_data={(0, k): ["save head"] for k in range(nstp)})
    if flow_file is None:
        ml.write_input()
    else:
        (ws / "manual.bcf").write_text(flow_file)
        bcf.file_name[0] = "manual.bcf"
        ml.write_input(SelPackList=[p.name[0] for p in ml.packagelist if p is not bcf])
    return ml


def _run_heads(ml, name, disu):
    success, buff = ml.run_model(silent=True, report=True)
    assert success, f"USG-T run in {ml.model_ws} did not terminate normally"
    assert any("normal termination" in line.lower() for line in buff)
    path = os.path.join(ml.model_ws, name)
    if disu:
        hds = HeadUFile(path)
        times = hds.get_times()
        return times, np.array([np.concatenate(hds.get_data(totim=t)) for t in times])
    hds = HeadFile(path)
    return hds.get_times(), hds.get_alldata()


def _perturb_record(text, prop, layer, factor):
    """Scale the data line that follows a hand-written array control record."""
    lines = text.splitlines()
    tag = "#ksat" if prop == "ksat" else f"#{prop} layer {layer}"
    (idx,) = [i for i, ln in enumerate(lines) if ln.endswith(tag)]
    values = [float(v) * factor for v in lines[idx + 1].split()]
    lines[idx + 1] = " ".join(f"{v:.6E}" for v in values)
    return "\n".join(lines) + "\n"


@requires_exe(USGT_EXE)
@pytest.mark.parametrize("name", list(_BCF_EXE_CASES))
def test_usgt_exe_bcf_written_matches_fortran_order(function_tmpdir, name):
    """FloPy-written BCF = hand-written BCF in the READBCF order (all cells and
    time steps; identical records, so atol=1e-6 only absorbs formatting), and
    a perturbed hand-written variant shows the heads depend on the records the
    defect used to shift, omit or mislabel."""
    from autotest.test_usg_transport import _bcf_manual_text

    case, (prop, layer, factor) = _BCF_EXE_CASES[name]
    case = dict(case, ncol=3)
    disu = case.get("disu", False)
    manual_text = _bcf_manual_text(case)
    perturbed_text = _perturb_record(manual_text, prop, layer, factor)

    t_w, h_w = _run_heads(
        _bcf_exe_model(function_tmpdir / "written", case), "bcf.hds", disu
    )
    t_m, h_m = _run_heads(
        _bcf_exe_model(function_tmpdir / "manual", case, manual_text), "bcf.hds", disu
    )
    t_p, h_p = _run_heads(
        _bcf_exe_model(function_tmpdir / "perturbed", case, perturbed_text),
        "bcf.hds",
        disu,
    )

    assert h_w.shape == h_m.shape == h_p.shape
    assert t_w == t_m == t_p
    assert np.all(np.isfinite(h_w)) and np.all(np.isfinite(h_m))
    np.testing.assert_allclose(h_w, h_m, rtol=0, atol=1e-6)
    assert np.abs(h_p - h_m).max() > 1e-3


@requires_exe(USGT_EXE)
@pytest.mark.parametrize(
    "subdir, namefile, heads, disu",
    [
        ("Ex3_CLN_Conduit/Conduit", "Conduit.nam", "Conduit.hds", False),
        ("Ex4_Dual_Domain", "DualDomain.nam", "DualDomain.hds", False),
        ("Ex8_Lake", "lak_usg_01.nam", "lak_usg_01.hds", True),
    ],
    ids=["ex3-fixed-2layers", "ex4-fixed-1layer", "ex8-free-disu"],
)
def test_usgt_exe_bcf_example_rewrite_keeps_heads(
    function_tmpdir, subdir, namefile, heads, disu
):
    """A real BCF example loaded and rewritten by FloPy gives the same heads as
    the original files (copied to a new folder; the examples stay untouched).
    atol=1e-5 m is two orders below the smallest HCLOSE of these models
    (1e-3 for Ex3); identical records are expected to give identical heads."""
    import shutil

    from autotest.conftest import get_example_data_path

    original = function_tmpdir / "original"
    shutil.copytree(get_example_data_path() / "mfusg_transport" / subdir, original)
    ml = MfUsg.load(namefile, model_ws=str(original), exe_name=USGT_EXE, check=False)
    ml.change_model_ws(str(function_tmpdir / "rewritten"))
    ml.write_input()

    import flopy

    results = []
    for ws in (original, function_tmpdir / "rewritten"):
        success, buff = flopy.mbase.run_model(
            USGT_EXE, namefile, model_ws=str(ws), silent=True, report=True
        )
        assert success and any("normal termination" in b.lower() for b in buff)
        path = os.path.join(ws, heads)
        if disu:
            hds = HeadUFile(path)
            times = hds.get_times()
            data = np.array([np.concatenate(hds.get_data(totim=t)) for t in times])
        else:
            hds = HeadFile(path)
            times, data = hds.get_times(), hds.get_alldata()
        results.append((times, data))
    (t_o, h_o), (t_r, h_r) = results
    assert h_o.shape == h_r.shape and t_o == t_r
    assert np.all(np.isfinite(h_o)) and np.all(np.isfinite(h_r))
    np.testing.assert_allclose(h_r, h_o, rtol=0, atol=1e-5)


@requires_exe(USGT_EXE)
def test_usgt_exe_bcf_fixed_format_reals_echoed(function_tmpdir):
    """READBCF echoes the fixed-format item 1 reals (HDRY with G13.5, WETFCT
    with F10.5): HDRY=-1e30 and WETFCT=1.2345678 written by FloPy must be read
    unchanged, not scaled by the F10.3 implied decimals."""
    ml = MfUsg(modelname="fx", model_ws=str(function_tmpdir), exe_name=USGT_EXE)
    ModflowDis(
        ml,
        nlay=1,
        nrow=1,
        ncol=3,
        nper=1,
        steady=True,
        delr=10.0,
        delc=10.0,
        top=30.0,
        botm=0.0,
    )
    MfUsgBas(ml, ibound=1, strt=25.0, ifrefm=False)
    MfUsgBcf(ml, laycon=3, hdry=-1.0e30, iwdflg=1, wetfct=1.2345678, iwetit=1, hy=1.0)
    MfUsgSms(
        ml,
        hclose=1e-5,
        hiclose=1e-7,
        mxiter=250,
        iter1=600,
        iprsms=1,
        nonlinmeth=1,
        linmeth=1,
        theta=0.7,
        akappa=0.07,
        gamma=0.1,
        amomentum=0.0,
        numtrack=200,
        btol=1.1,
        breduc=0.2,
        reslim=1.0,
        iacl=1,
        norder=0,
        level=7,
        north=14,
        iredsys=0,
        rrctol=0.0,
        idroptol=1,
        epsrn=1e-3,
    )
    ModflowChd(
        ml, stress_period_data={0: [[0, 0, 0, 26.0, 26.0], [0, 0, 2, 24.0, 24.0]]}
    )
    MfUsgOc(ml, stress_period_data={(0, 0): ["save head"]})
    ml.write_input()
    success, buff = ml.run_model(silent=True, report=True)
    assert success and any("normal termination" in b.lower() for b in buff)

    listing = (function_tmpdir / "fx.list").read_text().splitlines()
    hdry = [ln for ln in listing if "HEAD AT CELLS THAT CONVERT TO DRY=" in ln]
    wetfct = [ln for ln in listing if "WETTING FACTOR=" in ln]
    assert len(hdry) == 1 and len(wetfct) == 1
    # G13.5 keeps five significant digits, F10.5 five decimals
    assert np.isclose(float(hdry[0].split("=")[1]), -1.0e30, rtol=1e-4)
    assert np.isclose(
        float(wetfct[0].split("=")[1].split()[0]), 1.2345678, rtol=0, atol=5e-6
    )


# DPF budget characterization with constant heads (USG-T 2.8.1). In
# GWF2DPFU1FM the mobile-immobile exchange also applies at CHD cells
# (IBOUNDIM!=0 and IBOUND!=0) and GWF2DPFU1BDS reports their IMMOBILE STORAGE,
# while the CONSTANT HEAD term (GWF2BCFU1BDCH/BDADJ) only sums the flows
# between cells. In this small LPF/Richards/DPF case the printed residual
# D = IN - OUT therefore coincides, within tolerance, with the signed
# IMMOBILE STORAGE q summed under the CHD cells. These tests characterize the
# engine behaviour (they pass before and after the documentation change);
# they do not claim that the printed budget is closed or that the result
# holds for every DPF model.
#
# Tolerances: the listing prints budget rates with four decimals (rounding
# up to 5e-5 per printed value) and the IDPFCB record is single precision
# (about 3e-6 on values near 20), so D and q are compared at 2e-4 L3/T.
# Heads of runs with identical inputs are compared at 1e-6.
_DPFCHD_TOL = 2.0e-4
_DPFCHD_NSTP = 5
_DPFCHD_CHD = [(0, 0, 0), (1, 0, 2)]  # (layer, row, column), zero-based
_DPFCHD_DPF = {  # per layer (1, 2); every layer has LAYCON 5 in this case
    "hnewim": (23.5, 23.0),
    "phif": (0.2, 0.1),
    "ddftr": (0.01, 0.02),
    "sc1im": (1.0e-4, 2.0e-4),
    "sc2im": (0.3, 0.25),
    "alphaim": (1.0, 0.8),
    "betaim": (2.0, 2.5),
    "srim": (0.1, 0.15),
    "brookim": (3.5, 4.0),
}


def _dpfchd_iboundim(no_immobile_under_chd):
    iboundim = np.ones((2, 1, 3), dtype=int)
    if no_immobile_under_chd:
        for lay, row, col in _DPFCHD_CHD:
            iboundim[lay, row, col] = 0
    return iboundim


def _dpfchd_manual_dpf(iboundim):
    """DPF written by hand in the GWF2DPFU1AR read order (gwf2dpf1u1.f):
    item 1, IBOUNDIM and HNEWIM for every layer, then for each layer PHIF,
    DDFTR, SC1IM, SC2IM (LAYCON!=0) and the Richards arrays (LAYCON=5)."""
    lines = ["        61         0         0"]
    for k in range(2):
        lines += [
            f"INTERNAL 1 (FREE) -1 #iboundim layer {k + 1}",
            " ".join(str(v) for v in iboundim[k, 0]),
        ]
    lines += [
        f"CONSTANT {_DPFCHD_DPF['hnewim'][k]} #hnewim layer {k + 1}" for k in range(2)
    ]
    for k in range(2):
        for name in (
            "phif",
            "ddftr",
            "sc1im",
            "sc2im",
            "alphaim",
            "betaim",
            "srim",
            "brookim",
        ):
            lines.append(f"CONSTANT {_DPFCHD_DPF[name][k]} #{name} layer {k + 1}")
    return "\n".join(lines) + "\n"


def _dpfchd_run(ws, no_immobile_under_chd=False, manual=False):
    """LPF LAYTYP [5,5] with BAS RICHARDS, 1x3x2 grid, CHD 26 (layer 1,
    column 1) and 24 (layer 2, column 3), recharge 0.05, five 1-day steps,
    and DPF with immobile Richards arrays. Returns per-step listing rates,
    the IMMOBILE STORAGE records and the heads."""
    ml = MfUsg(modelname="dc", model_ws=str(ws), exe_name=USGT_EXE)
    ModflowDis(
        ml,
        nlay=2,
        nrow=1,
        ncol=3,
        nper=1,
        nstp=_DPFCHD_NSTP,
        perlen=5.0,
        steady=False,
        delr=10.0,
        delc=10.0,
        top=30.0,
        botm=[20.0, 10.0],
    )
    MfUsgBas(ml, ibound=1, strt=25.0, richards=True)
    lpf_values = {
        "hk": (21.0, 22.0),
        "vka": (3.1, 3.2),
        "ss": (1.1e-5, 1.2e-5),
        "sy": (0.051, 0.052),
        "alpha": (0.71, 0.72),
        "beta": (8.1, 8.2),
        "sr": (0.091, 0.092),
        "brook": (10.1, 10.2),
    }
    MfUsgLpf(
        ml,
        laytyp=[5, 5],
        **{n: [[[v[k]] * 3] for k in range(2)] for n, v in lpf_values.items()},
    )
    MfUsgSms(
        ml,
        hclose=1e-5,
        hiclose=1e-7,
        mxiter=250,
        iter1=600,
        iprsms=1,
        nonlinmeth=1,
        linmeth=1,
        theta=0.7,
        akappa=0.07,
        gamma=0.1,
        amomentum=0.0,
        numtrack=200,
        btol=1.1,
        breduc=0.2,
        reslim=1.0,
        iacl=1,
        norder=0,
        level=7,
        north=14,
        iredsys=0,
        rrctol=0.0,
        idroptol=1,
        epsrn=1e-3,
    )
    ModflowChd(
        ml, stress_period_data={0: [[0, 0, 0, 26.0, 26.0], [1, 0, 2, 24.0, 24.0]]}
    )
    ModflowRch(ml, rech=0.05, nrchop=1)
    MfUsgOc(
        ml,
        stress_period_data={
            (0, k): ["save head", "save budget", "print budget"]
            for k in range(_DPFCHD_NSTP)
        },
    )
    iboundim = _dpfchd_iboundim(no_immobile_under_chd)
    dpf = MfUsgDpf(
        ml,
        ipakcb=61,
        iboundim=iboundim,
        **{n: [[[v[k]] * 3] for k in range(2)] for n, v in _DPFCHD_DPF.items()},
    )
    ml.add_output_file(61, fname="dc.dpf.cbc", binflag=True, package="DPF")
    if manual:
        (ws / "manual.dpf").write_text(_dpfchd_manual_dpf(iboundim))
        dpf.file_name[0] = "manual.dpf"
        ml.write_input(SelPackList=[p.name[0] for p in ml.packagelist if p is not dpf])
    else:
        ml.write_input()

    success, buff = ml.run_model(silent=True, report=True)
    assert success and any("normal termination" in b.lower() for b in buff)
    listing = (ws / "dc.list").read_text()
    assert "FAILED TO CONVERGE" not in listing.upper()

    rates, _ = MfusgListBudget(str(ws / "dc.list")).get_budget()
    assert len(rates) == _DPFCHD_NSTP
    assert np.allclose(rates["totim"], np.arange(1.0, _DPFCHD_NSTP + 1.0))
    cbc = CellBudgetFile(str(ws / "dc.dpf.cbc"), precision="single")
    assert any("IMMOBILE STORAGE" in t.decode().strip() for t in cbc.textlist)
    immobile = np.array(
        [np.asarray(a).reshape(2, 1, 3) for a in cbc.get_data(text="IMMOBILE STORAGE")]
    )
    hds = HeadFile(str(ws / "dc.hds"))
    heads = hds.get_alldata()
    assert immobile.shape == (_DPFCHD_NSTP, 2, 1, 3)
    assert heads.shape == (_DPFCHD_NSTP, 2, 1, 3)
    assert np.allclose(hds.get_times(), np.arange(1.0, _DPFCHD_NSTP + 1.0))
    assert np.all(np.isfinite(heads)) and np.all(np.isfinite(immobile))
    return rates, immobile, heads


def _dpfchd_under_chd(immobile):
    return np.array(
        [sum(step[lay, row, col] for lay, row, col in _DPFCHD_CHD) for step in immobile]
    )


@requires_exe(USGT_EXE)
def test_usgt_exe_dpf_chd_residual_equals_immobile_storage_under_chd(
    function_tmpdir,
):
    """D = IN - OUT is non-trivial and coincides, within tolerance, with the
    signed IMMOBILE STORAGE q summed under the CHD cells in every step."""
    rates, immobile, _ = _dpfchd_run(function_tmpdir)
    d = np.asarray(rates["IN-OUT"], dtype=float)
    q = _dpfchd_under_chd(immobile)
    assert len(d) == len(q) == _DPFCHD_NSTP
    assert np.all(np.abs(d) > 1.0)  # far above the tolerance: not a trivial case
    np.testing.assert_allclose(d, q, rtol=0, atol=_DPFCHD_TOL)
    # the printed percentage uses the mean of total in and total out
    total_in = np.asarray(rates["TOTAL_IN"], dtype=float)
    total_out = np.asarray(rates["TOTAL_OUT"], dtype=float)
    percent = 100.0 * d / ((total_in + total_out) / 2.0)
    np.testing.assert_allclose(
        rates["PERCENT_DISCREPANCY"], percent, rtol=0, atol=0.01
    )  # printed with two decimals


@requires_exe(USGT_EXE)
def test_usgt_exe_dpf_no_immobile_under_chd_closes_budget(function_tmpdir):
    """Control for this case only: without immobile cells under the CHD cells
    (a different physical model there) the printed budget closes and the
    mobile heads coincide with the original run."""
    rates0, _, heads0 = _dpfchd_run(function_tmpdir / "original")
    rates1, immobile1, heads1 = _dpfchd_run(
        function_tmpdir / "no_immobile_under_chd", no_immobile_under_chd=True
    )
    assert len(rates1) == len(rates0) == _DPFCHD_NSTP
    assert np.all(np.abs(np.asarray(rates0["IN-OUT"], dtype=float)) > 1.0)
    assert np.all(np.abs(np.asarray(rates1["IN-OUT"], dtype=float)) <= _DPFCHD_TOL)
    assert np.all(np.abs(_dpfchd_under_chd(immobile1)) == 0.0)
    np.testing.assert_allclose(heads1, heads0, rtol=0, atol=1e-6)


@requires_exe(USGT_EXE)
def test_usgt_exe_dpf_chd_manual_dpf_matches_writer(function_tmpdir):
    """A DPF file written by hand in the GWF2DPFU1AR order gives the same heads
    and per-cell IMMOBILE STORAGE as the file written by MfUsgDpf."""
    rates_w, immobile_w, heads_w = _dpfchd_run(function_tmpdir / "written")
    rates_m, immobile_m, heads_m = _dpfchd_run(function_tmpdir / "manual", manual=True)
    assert len(rates_w) == len(rates_m) == _DPFCHD_NSTP
    np.testing.assert_allclose(heads_m, heads_w, rtol=0, atol=1e-6)
    np.testing.assert_allclose(immobile_m, immobile_w, rtol=0, atol=_DPFCHD_TOL)
    np.testing.assert_allclose(
        rates_m["IN-OUT"], rates_w["IN-OUT"], rtol=0, atol=_DPFCHD_TOL
    )


# D7 stage 1: HNEWIM precision through the executable. With DDFTR=0 the
# immobile heads keep HNEWIM (HNEWIM = HOLDIM), and with DPOUT/DPIO the saved
# IMMOBILE HEAD record is REAL*8 (ULASAV8), so the saved values are the
# values the executable read. Saved files are read with struct, not FloPy.

_D7X_HNEWIM = np.array(
    [
        [23.1234567890123, 23.1234567890124, 22.75],
        [21.987654321098765, 22.000000123456789, 21.5],
    ]
)  # the first two values of each layer are equal in float32
_D7X_NOT_REAL = 23.1234567890123


def _d7x_model(ws, options=("dpio",)):
    ml = MfUsg(modelname="d7", model_ws=str(ws), exe_name=USGT_EXE)
    ModflowDis(
        ml,
        nlay=2,
        nrow=1,
        ncol=3,
        nper=1,
        nstp=1,
        perlen=1.0,
        steady=False,
        delr=10.0,
        delc=10.0,
        top=30.0,
        botm=[20.0, 10.0],
    )
    flags = {"dpin": "double_prec", "dpout": "double_out", "dpio": "double_io"}
    MfUsgBas(ml, ibound=1, strt=25.0, **{flags[o]: True for o in options})
    MfUsgLpf(ml, laytyp=0, hk=10.0, vka=1.0, ss=1.0e-5)
    MfUsgSms(
        ml,
        hclose=1e-9,
        hiclose=1e-11,
        mxiter=250,
        iter1=600,
        iprsms=1,
        nonlinmeth=1,
        linmeth=1,
        theta=0.7,
        akappa=0.07,
        gamma=0.1,
        amomentum=0.0,
        numtrack=200,
        btol=1.1,
        breduc=0.2,
        reslim=1.0,
        iacl=1,
        norder=0,
        level=7,
        north=14,
        iredsys=0,
        rrctol=0.0,
        idroptol=1,
        epsrn=1e-3,
    )
    ModflowChd(
        ml, stress_period_data={0: [[0, 0, 0, 26.0, 26.0], [0, 0, 2, 24.0, 24.0]]}
    )
    MfUsgOc(ml, stress_period_data={(0, 0): ["save head"]})
    ml.add_output_file(72, fname="d7.imh", binflag=True, package="DPF")
    return ml


def _d7x_dpf_kwargs():
    return {
        "idpfhd": 72,
        "iboundim": 1,
        "phif": 0.2,
        "ddftr": 0.0,
        "sc1im": 1.0e-4,
    }


def _d7x_manual(ws, hnewim_records, options=("dpio",), files=()):
    """Hand-written DPF (GWF2DPFU1AR order, LAYCON 0); other packages by
    FloPy. files: (name, bytes, unit) registered as DATA(BINARY)."""
    ml = _d7x_model(ws, options)
    MfUsgDpf(ml, filenames="manual.dpf", **_d7x_dpf_kwargs())
    lines = [
        "         0        72         0",
        "CONSTANT 1 #iboundim layer 1",
        "CONSTANT 1 #iboundim layer 2",
    ]
    lines += hnewim_records
    for k in (1, 2):
        lines += [
            f"CONSTANT 0.2 #phif layer {k}",
            f"CONSTANT 0.0 #ddftr layer {k}",
            f"CONSTANT 1.0E-4 #sc1im layer {k}",
        ]
    (ws / "manual.dpf").write_text("\n".join(lines) + "\n")
    for name, data, unit in files:
        (ws / name).write_bytes(data)
        ml.add_external(name, unit, binflag=True)
    ml.write_input(SelPackList=["DIS", "BAS6", "LPF", "SMS", "CHD", "OC"])
    return ml


def _d7x_binary(values=_D7X_HNEWIM):
    """U1DREL8 stream layout: per layer KSTP,KPER (int32), PERTIM,TOTIM
    (REAL*8), TEXT (16), NSTRT,NNDLAY,ILAY (int32), then REAL*8 values."""
    import struct

    out, nstrt = b"", 1
    for k, row in enumerate(values, start=1):
        out += struct.pack(
            "<iidd16siii",
            1,
            1,
            0.0,
            0.0,
            b"  IMMOBILE HEAD ",
            nstrt,
            nstrt + len(row) - 1,
            k,
        )
        out += struct.pack(f"<{len(row)}d", *row)
        nstrt += len(row)
    return out


def _d7x_run_saved(ml, fname="d7.imh", unstructured=False):
    """Run and return the saved REAL*8 IMMOBILE HEAD values per layer
    (NCOL,NROW in the header, or NSTRT,NNDLAY for unstructured grids)."""
    import struct

    success, buff = ml.run_model(silent=True, report=True)
    assert success and any("normal termination" in ln.lower() for ln in buff)
    data = (Path(ml.model_ws) / fname).read_bytes()
    layers, off = [], 0
    while off < len(data):
        _, _, _, _, text, n1, n2, _ = struct.unpack_from("<iidd16siii", data, off)
        assert text.strip().startswith(b"IMMOBILE HEAD")
        n = n2 - n1 + 1 if unstructured else n1 * n2
        layers.append(np.array(struct.unpack_from(f"<{n}d", data, off + 52)))
        off += 52 + 8 * n
    return layers if unstructured else np.array(layers)


@requires_exe(USGT_EXE)
@pytest.mark.parametrize("options", [("dpio",), ("dpin", "dpout")])
def test_usgt_exe_hnewim_double_written_by_flopy(function_tmpdir, options):
    ml = _d7x_model(function_tmpdir, options)
    dpf = MfUsgDpf(ml, hnewim=_D7X_HNEWIM.reshape(2, 1, 3), **_d7x_dpf_kwargs())
    assert dpf.hnewim.dtype == np.float64
    ml.write_input()
    saved = _d7x_run_saved(ml)
    np.testing.assert_array_equal(saved, _D7X_HNEWIM)
    assert saved[0, 0] != saved[0, 1]


@requires_exe(USGT_EXE)
def test_usgt_exe_hnewim_text_matches_manual_binary(function_tmpdir):
    """FloPy text, a hand-written REAL*8 binary record, and the binary file
    loaded and rewritten by FloPy all reach the executable identically."""
    written = _d7x_model(function_tmpdir / "text")
    MfUsgDpf(written, hnewim=_D7X_HNEWIM.reshape(2, 1, 3), **_d7x_dpf_kwargs())
    written.write_input()
    (function_tmpdir / "bin").mkdir()
    manual = _d7x_manual(
        function_tmpdir / "bin",
        [f"EXTERNAL -71 1.0 #hnewim layer {k}" for k in (1, 2)],
        files=[("h.bin", _d7x_binary(), 71)],
    )
    saved_text = _d7x_run_saved(written)
    saved_bin = _d7x_run_saved(manual)
    np.testing.assert_array_equal(saved_text, _D7X_HNEWIM)
    np.testing.assert_array_equal(saved_bin, saved_text)

    reloaded = MfUsg.load("d7.nam", model_ws=str(function_tmpdir / "bin"), check=False)
    loaded = reloaded.get_package("DPF").hnewim
    assert loaded.dtype == np.float64
    np.testing.assert_array_equal(loaded.array.reshape(2, 3), _D7X_HNEWIM)
    rewrite = _d7x_model(function_tmpdir / "rewrite")
    MfUsgDpf(rewrite, hnewim=loaded, **_d7x_dpf_kwargs())  # loaded Util3d
    rewrite.write_input()
    rewritten = (function_tmpdir / "rewrite" / "d7.dpf").read_text()
    assert "EXTERNAL" not in rewritten and "(FREE)" in rewritten
    np.testing.assert_array_equal(_d7x_run_saved(rewrite), _D7X_HNEWIM)


@requires_exe(USGT_EXE)
def test_usgt_exe_hnewim_uniform_not_real_written_internal(function_tmpdir):
    ml = _d7x_model(function_tmpdir)
    MfUsgDpf(ml, hnewim=_D7X_NOT_REAL, **_d7x_dpf_kwargs())
    ml.write_input()
    assert "CONSTANT" not in [
        ln.split()[0]
        for ln in (function_tmpdir / "d7.dpf").read_text().splitlines()
        if "hnewim" in ln
    ]
    np.testing.assert_array_equal(_d7x_run_saved(ml), np.full((2, 3), _D7X_NOT_REAL))


@requires_exe(USGT_EXE)
def test_usgt_exe_hnewim_assigned_multiplier_is_real(function_tmpdir):
    ml = _d7x_model(function_tmpdir)
    dpf = MfUsgDpf(ml, hnewim=_D7X_HNEWIM.reshape(2, 1, 3), **_d7x_dpf_kwargs())
    dpf.hnewim.cnstnt = 0.1
    ml.write_input()
    saved = _d7x_run_saved(ml)
    np.testing.assert_array_equal(saved, dpf.hnewim.array.reshape(2, 3))
    np.testing.assert_array_equal(saved, _D7X_HNEWIM * np.float64(np.float32(0.1)))
    assert np.any(saved != _D7X_HNEWIM * 0.1)


@requires_exe(USGT_EXE)
def test_usgt_exe_hnewim_binary_multiplier_not_applied(function_tmpdir):
    ml = _d7x_manual(
        function_tmpdir,
        [f"EXTERNAL -71 2.0 #hnewim layer {k}" for k in (1, 2)],
        files=[("h.bin", _d7x_binary(), 71)],
    )
    np.testing.assert_array_equal(_d7x_run_saved(ml), _D7X_HNEWIM)
    with pytest.warns(UserWarning, match="not applied"):
        loaded = MfUsg.load("d7.nam", model_ws=str(function_tmpdir), check=False)
    np.testing.assert_array_equal(
        loaded.get_package("DPF").hnewim.array.reshape(2, 3), _D7X_HNEWIM
    )


@requires_exe(USGT_EXE)
def test_usgt_exe_hnewim_fixed_format_records(function_tmpdir):
    """HNEWIM written with fixed control records (I10,F10.0,A20,I10) in a
    free-format model; U1DREL8 detects them by the missing keyword."""
    ml = _d7x_model(function_tmpdir)
    ml.array_free_format = False
    dpf = MfUsgDpf(ml, hnewim=_D7X_HNEWIM.reshape(2, 1, 3), **_d7x_dpf_kwargs())
    ml.array_free_format = True
    ml.write_input()
    records = [
        ln
        for ln in (function_tmpdir / "d7.dpf").read_text().splitlines()
        if "hnewim" in ln
    ]
    assert all(ln[20:40].strip() == "(3E24.16)" and int(ln[0:10]) > 0 for ln in records)
    assert dpf.hnewim.dtype == np.float64
    np.testing.assert_array_equal(_d7x_run_saved(ml), _D7X_HNEWIM)


@requires_exe(USGT_EXE)
def test_usgt_exe_hnewim_without_dpin_is_single(function_tmpdir):
    """Control: without DPIN/DPIO FloPy keeps float32 and the executable reads
    REAL (saved in double with DPOUT)."""
    ml = _d7x_model(function_tmpdir, ("dpout",))
    dpf = MfUsgDpf(ml, hnewim=_D7X_HNEWIM.reshape(2, 1, 3), **_d7x_dpf_kwargs())
    assert dpf.hnewim.dtype == np.float32
    ml.write_input()
    assert "(3E15.6)" in (function_tmpdir / "d7.dpf").read_text()
    saved = _d7x_run_saved(ml)
    np.testing.assert_allclose(saved, _D7X_HNEWIM, rtol=1e-6)
    assert saved[0, 0] == saved[0, 1]


@requires_exe(USGT_EXE)
def test_usgt_exe_hnewim_disu_irregular_double(function_tmpdir):
    """DISU with 3 and 2 nodes per layer: U1DREL8 reads NDSLAY values per
    layer; the saved records carry NSTRT,NNDLAY."""
    from flopy.mfusg import MfUsgDisU

    values = np.array(
        [
            23.1234567890123,
            23.1234567890124,
            22.75,
            21.987654321098765,
            22.000000123456789,
        ]
    )
    ml = MfUsg(
        modelname="u",
        model_ws=str(function_tmpdir),
        structured=False,
        exe_name=USGT_EXE,
    )
    MfUsgDisU(
        ml,
        nodes=5,
        nlay=2,
        njag=15,
        nodelay=[3, 2],
        iac=np.array([3, 4, 2, 3, 3]),
        ja=np.array([0, 1, 3, 1, 0, 2, 4, 2, 1, 3, 0, 4, 4, 1, 3]),
        ivsd=0,
        cl12=5.0,
        fahl=10.0,
        top=[30.0] * 3 + [20.0] * 2,
        bot=[20.0] * 3 + [10.0] * 2,
        area=100.0,
        nper=1,
        perlen=1.0,
        nstp=1,
        steady=False,
    )
    MfUsgBas(ml, ibound=1, strt=25.0, double_io=True)
    MfUsgLpf(ml, laytyp=0, hk=10.0, vka=1.0, ss=1.0e-5)
    MfUsgSms(
        ml, hclose=1e-9, hiclose=1e-11, mxiter=250, iter1=600, nonlinmeth=1, linmeth=1
    )
    MfUsgOc(ml, stress_period_data={(0, 0): ["save head"]})
    ml.add_output_file(72, fname="u.imh", binflag=True, package="DPF")
    dpf = MfUsgDpf(ml, hnewim=[values[:3], values[3:]], **_d7x_dpf_kwargs())
    assert [u.shape for u in dpf.hnewim.util_2ds] == [(3,), (2,)]
    ml.write_input()
    saved = _d7x_run_saved(ml, "u.imh", unstructured=True)
    np.testing.assert_array_equal(saved[0], values[:3])
    np.testing.assert_array_equal(saved[1], values[3:])


@requires_exe(USGT_EXE)
def test_usgt_exe_hnewim_layer_replaced_with_double_util2d(function_tmpdir):
    """A layer replaced by a float64 Util2d after construction reaches the
    executable with values that are distinct only in float64."""
    from flopy.utils import Util2d

    ml = _d7x_model(function_tmpdir, ("dpin", "dpout"))
    dpf = MfUsgDpf(ml, hnewim=0.0, **_d7x_dpf_kwargs())
    for k in range(2):
        dpf.hnewim[k] = Util2d(
            ml, (1, 3), np.float64, _D7X_HNEWIM[k : k + 1], name="replacement"
        )
    ml.write_input()
    saved = _d7x_run_saved(ml)
    np.testing.assert_array_equal(saved, _D7X_HNEWIM)
    assert saved[0, 0] != saved[0, 1]  # equal in float32


# D7 stage 2: binary HNEWIM output (EXTERNAL with a negative unit and a
# DATA(BINARY) file) read by the executable; DDFTR=0 and REAL*8 IMMOBILE HEAD
# output show the values read.


def _d7x_binary_dpf(ml, hnewim=_D7X_HNEWIM.reshape(2, 1, 3), layers=(0, 1)):
    dpf = MfUsgDpf(ml, hnewim=hnewim, **_d7x_dpf_kwargs())
    for k in layers:
        dpf.hnewim[k].format.binary = True
    return dpf


def _d7x_bin_size(ws, fname="d7.hnewim.bin"):
    return (Path(ws) / fname).stat().st_size


@requires_exe(USGT_EXE)
@pytest.mark.parametrize("options", [("dpio",), ("dpin", "dpout")])
def test_usgt_exe_hnewim_binary_double(function_tmpdir, options):
    ml = _d7x_model(function_tmpdir, options)
    _d7x_binary_dpf(ml)
    ml.write_input()
    assert _d7x_bin_size(function_tmpdir) == 2 * (52 + 24)
    saved = _d7x_run_saved(ml)
    np.testing.assert_array_equal(saved, _D7X_HNEWIM)
    assert saved[0, 0] != saved[0, 1]


@requires_exe(USGT_EXE)
def test_usgt_exe_hnewim_binary_single(function_tmpdir):
    ml = _d7x_model(function_tmpdir, ("dpout",))
    _d7x_binary_dpf(ml)
    ml.write_input()
    assert _d7x_bin_size(function_tmpdir) == 2 * (44 + 12)
    saved = _d7x_run_saved(ml)
    np.testing.assert_array_equal(saved, _D7X_HNEWIM.astype(np.float32))
    assert saved[0, 0] == saved[0, 1]


@requires_exe(USGT_EXE)
@pytest.mark.parametrize("options", [("dpio",), ("dpout",)])
def test_usgt_exe_hnewim_binary_multiplier_applied_once(function_tmpdir, options):
    ml = _d7x_model(function_tmpdir, options)
    dpf = _d7x_binary_dpf(ml)
    dpf.hnewim[1].cnstnt = 0.1
    ml.write_input()
    np.testing.assert_array_equal(_d7x_run_saved(ml), dpf.hnewim.array.reshape(2, 3))


@requires_exe(USGT_EXE)
@pytest.mark.parametrize(
    "value, options", [(_D7X_NOT_REAL, ("dpio",)), (22.75, ("dpout",))]
)
def test_usgt_exe_hnewim_binary_scalar(function_tmpdir, value, options):
    ml = _d7x_model(function_tmpdir, options)
    dpf = _d7x_binary_dpf(ml, hnewim=value)
    ml.write_input()
    np.testing.assert_array_equal(_d7x_run_saved(ml), dpf.hnewim.array.reshape(2, 3))


@requires_exe(USGT_EXE)
def test_usgt_exe_hnewim_binary_mixed_layers(function_tmpdir):
    ml = _d7x_model(function_tmpdir)
    _d7x_binary_dpf(ml, layers=(1,))
    ml.write_input()
    assert _d7x_bin_size(function_tmpdir) == 52 + 24
    np.testing.assert_array_equal(_d7x_run_saved(ml), _D7X_HNEWIM)


@requires_exe(USGT_EXE)
def test_usgt_exe_hnewim_binary_fixed_format(function_tmpdir):
    ml = _d7x_model(function_tmpdir)
    ml.array_free_format = False
    _d7x_binary_dpf(ml)
    ml.array_free_format = True
    ml.write_input()
    lines = [
        ln
        for ln in (function_tmpdir / "d7.dpf").read_text().splitlines()
        if "hnewim" in ln
    ]
    assert all(int(ln[0:10]) < 0 for ln in lines)
    np.testing.assert_array_equal(_d7x_run_saved(ml), _D7X_HNEWIM)


@requires_exe(USGT_EXE)
def test_usgt_exe_hnewim_binary_matches_manual_text(function_tmpdir):
    binary = _d7x_model(function_tmpdir / "bin")
    _d7x_binary_dpf(binary)
    binary.write_input()
    (function_tmpdir / "text").mkdir()
    records = []
    for k in (0, 1):
        records += [
            f"INTERNAL 1.0 (FREE) -1 #hnewim layer {k + 1}",
            " ".join(repr(float(v)) for v in _D7X_HNEWIM[k]),
        ]
    manual = _d7x_manual(function_tmpdir / "text", records)
    np.testing.assert_array_equal(_d7x_run_saved(binary), _d7x_run_saved(manual))


@requires_exe(USGT_EXE)
def test_usgt_exe_hnewim_binary_disu_irregular(function_tmpdir):
    from flopy.mfusg import MfUsgDisU

    values = np.array(
        [
            23.1234567890123,
            23.1234567890124,
            22.75,
            21.987654321098765,
            22.000000123456789,
        ]
    )
    ml = MfUsg(
        modelname="u",
        model_ws=str(function_tmpdir),
        structured=False,
        exe_name=USGT_EXE,
    )
    MfUsgDisU(
        ml,
        nodes=5,
        nlay=2,
        njag=15,
        nodelay=[3, 2],
        iac=np.array([3, 4, 2, 3, 3]),
        ja=np.array([0, 1, 3, 1, 0, 2, 4, 2, 1, 3, 0, 4, 4, 1, 3]),
        ivsd=0,
        cl12=5.0,
        fahl=10.0,
        top=[30.0] * 3 + [20.0] * 2,
        bot=[20.0] * 3 + [10.0] * 2,
        area=100.0,
        nper=1,
        perlen=1.0,
        nstp=1,
        steady=False,
    )
    MfUsgBas(ml, ibound=1, strt=25.0, double_io=True)
    MfUsgLpf(ml, laytyp=0, hk=10.0, vka=1.0, ss=1.0e-5)
    MfUsgSms(
        ml, hclose=1e-9, hiclose=1e-11, mxiter=250, iter1=600, nonlinmeth=1, linmeth=1
    )
    MfUsgOc(ml, stress_period_data={(0, 0): ["save head"]})
    ml.add_output_file(72, fname="u.imh", binflag=True, package="DPF")
    _d7x_binary_dpf(ml, hnewim=[values[:3], values[3:]])
    ml.write_input()
    assert _d7x_bin_size(function_tmpdir, "u.hnewim.bin") == (52 + 24) + (52 + 16)
    saved = _d7x_run_saved(ml, "u.imh", unstructured=True)
    np.testing.assert_array_equal(np.concatenate(saved), values)


@requires_exe(USGT_EXE)
def test_usgt_exe_hnewim_binary_foreign_file_and_new_folder(function_tmpdir):
    """With a foreign d7.hnewim.bin registered on unit 777, the run uses the
    free name, leaves the foreign file and entry, and a later workspace
    change writes a new file without touching the first one."""
    ml = _d7x_model(function_tmpdir / "a")
    (function_tmpdir / "a" / "d7.hnewim.bin").write_bytes(b"FOREIGN")
    ml.add_external("d7.hnewim.bin", 777, binflag=True)
    _d7x_binary_dpf(ml)
    ml.remove_external(unit=777)  # the run must not open the foreign file
    ml.add_external("foreign.lst", 777)
    (function_tmpdir / "a" / "foreign.lst").write_text("x\n")
    ml.write_input()
    assert (function_tmpdir / "a" / "d7.hnewim.bin").read_bytes() == b"FOREIGN"
    np.testing.assert_array_equal(_d7x_run_saved(ml), _D7X_HNEWIM)
    first = (function_tmpdir / "a" / "d7.hnewim.bin").read_bytes()
    ml.change_model_ws(str(function_tmpdir / "b"))
    (function_tmpdir / "b" / "foreign.lst").write_text("x\n")
    ml.write_input()
    np.testing.assert_array_equal(_d7x_run_saved(ml), _D7X_HNEWIM)
    assert (function_tmpdir / "a" / "d7.hnewim.bin").read_bytes() == first


@requires_exe(USGT_EXE)
def test_usgt_exe_hnewim_binary_own_path_replaced_by_link(function_tmpdir):
    """The own file is kept as preserved.bin and its path is replaced by a
    link to a foreign file: the next write uses a new file, the foreign and
    preserved files are unchanged and the executable reads the new file."""
    ws = function_tmpdir
    ml = _d7x_model(ws)
    _d7x_binary_dpf(ml)
    ml.write_input()
    os.rename(ws / "d7.hnewim.bin", ws / "preserved.bin")
    (ws / "foreign.bin").write_bytes(b"FOREIGN")
    os.symlink("foreign.bin", ws / "d7.hnewim.bin")
    preserved = (ws / "preserved.bin").read_bytes()
    ml.write_input()
    assert (ws / "foreign.bin").read_bytes() == b"FOREIGN"
    assert (ws / "preserved.bin").read_bytes() == preserved
    assert "d7.hnewim_1.bin" in (ws / "d7.nam").read_text()
    np.testing.assert_array_equal(_d7x_run_saved(ml), _D7X_HNEWIM)


# G02: RIV on DISU. The second RIV header field is IRIVCB (GWF2RIV7U1AR):
# > 0 saves RIVER LEAKAGE on that unit, < 0 prints the rates when ICBCFL != 0,
# 0 neither. A 3-node DISU line with CHD at node 1 and RIV at node 3.


def _rivx_model(ws, ipakcb, manual=None, riv=None):
    """manual: None (FloPy writes RIV) or the text of a hand-written RIV;
    riv: optional MfUsgRiv keyword arguments for the FloPy-written RIV."""
    from flopy.mfusg import MfUsgChd, MfUsgDisU, MfUsgRiv

    ml = MfUsg(modelname="r", model_ws=str(ws), structured=False, exe_name=USGT_EXE)
    MfUsgDisU(
        ml,
        nodes=3,
        nlay=1,
        njag=7,
        iac=np.array([2, 3, 2]),
        ja=np.array([0, 1, 1, 0, 2, 2, 1]),
        ivsd=0,
        cl12=5.0,
        fahl=10.0,
        area=100.0,
        top=[10.0] * 3,
        bot=[0.0] * 3,
        nper=1,
        perlen=1.0,
        nstp=1,
        steady=True,
    )
    MfUsgBas(ml, ibound=1, strt=10.0)
    MfUsgLpf(ml, laytyp=0, hk=1.0, ipakcb=0)
    MfUsgSms(
        ml, hclose=1e-9, hiclose=1e-11, mxiter=100, iter1=100, nonlinmeth=1, linmeth=1
    )
    MfUsgChd(ml, stress_period_data={0: [[0, 10.0, 10.0]]})
    MfUsgOc(ml, stress_period_data={(0, 0): ["save head", "save budget"]})
    if manual is None:
        kwargs = riv or {"stress_period_data": {0: [[2, 12.0, 5.0, 8.0]]}}
        MfUsgRiv(ml, ipakcb=ipakcb, **kwargs)
        ml.write_input()
    else:
        MfUsgRiv(ml, ipakcb=ipakcb, filenames="manual.riv")  # nam + CBC unit only
        Path(ws).mkdir(parents=True, exist_ok=True)
        (Path(ws) / "manual.riv").write_text(manual)
        ml.write_input(SelPackList=["DISU", "BAS6", "LPF", "SMS", "CHD", "OC"])
    return ml


def _rivx_run(ml):
    success, buff = ml.run_model(silent=True, report=True)
    assert success and any("normal termination" in b.lower() for b in buff)
    lst = (Path(ml.model_ws) / "r.list").read_text()
    assert "FAILED TO CONVERGE" not in lst.upper()
    heads = HeadUFile(str(Path(ml.model_ws) / "r.hds")).get_data()[0]
    heads = np.asarray(heads, dtype=float).ravel()
    assert heads.shape == (3,) and np.all(np.isfinite(heads))
    return heads, lst


def _rivx_budget_in(lst):
    line = next(ln for ln in lst.splitlines() if "RIVER LEAKAGE =" in ln)
    return float(line.split("=")[1].split()[0])


@requires_exe(USGT_EXE)
@pytest.mark.parametrize("unit", [55, 61])
def test_usgt_exe_riv_disu_cbc_on_requested_unit(function_tmpdir, unit):
    written = _rivx_model(function_tmpdir / "w", unit)
    header = (function_tmpdir / "w" / "r.riv").read_text().splitlines()[1].split()
    assert header[:2] == ["1", str(unit)]
    manual = _rivx_model(
        function_tmpdir / "m",
        unit,
        manual=f"# manual RIV\n 1 {unit}\n 1\n 3 12.0 5.0 8.0\n",
    )
    results = []
    for ml in (written, manual):
        heads, lst = _rivx_run(ml)
        assert f'SAVING "   RIVER LEAKAGE" ON UNIT  {unit}' in lst
        nam = (Path(ml.model_ws) / "r.nam").read_text().split("\n")
        cbc = [
            ln.split()[2]
            for ln in nam
            if ln.startswith("DATA(BINARY)") and int(ln.split()[1]) == unit
        ]
        assert cbc == ["r.cbc"]
        assert not list(Path(ml.model_ws).glob("fort.*"))  # no unregistered unit
        cbf = CellBudgetFile(str(Path(ml.model_ws) / cbc[0]), precision="single")
        assert [n.strip() for n in cbf.get_unique_record_names()] == [b"RIVER LEAKAGE"]
        assert cbf.get_times() == [1.0] and cbf.get_kstpkper() == [(0, 0)]
        q = cbf.get_data(text="RIVER LEAKAGE")[0]
        assert q["node"].tolist() == [3] and np.isfinite(q["q"][0]) and q["q"][0] > 0
        # Q = C (stage - h); float32 storage in the budget file
        np.testing.assert_allclose(q["q"][0], 5.0 * (12.0 - heads[2]), rtol=1e-6)
        assert abs(_rivx_budget_in(lst) - q["q"][0]) < 1e-4  # listing: 4 decimals
        results.append((heads, q["q"][0]))
    # identical input values -> identical results
    np.testing.assert_array_equal(results[0][0], results[1][0])
    assert results[0][1] == results[1][1]


@requires_exe(USGT_EXE)
@pytest.mark.parametrize("unit", [0, -1])
def test_usgt_exe_riv_disu_irivcb_zero_and_negative(function_tmpdir, unit):
    reference, _ = _rivx_run(_rivx_model(function_tmpdir / "ref", 55))
    ml = _rivx_model(function_tmpdir / "w", unit)
    header = (function_tmpdir / "w" / "r.riv").read_text().splitlines()[1].split()
    assert header[:2] == ["1", str(unit)]
    heads, lst = _rivx_run(ml)
    np.testing.assert_array_equal(heads, reference)
    nam = (function_tmpdir / "w" / "r.nam").read_text()
    assert "r.cbc" not in nam and not list((function_tmpdir / "w").glob("*.cbc"))
    assert "CELL-BY-CELL FLOWS WILL BE SAVED" not in lst
    assert 'SAVING "   RIVER LEAKAGE"' not in lst
    assert ("CELL-BY-CELL FLOWS WILL BE PRINTED" in lst) == (unit < 0)
    assert ("RIVER LEAKAGE   PERIOD    1   STEP   1" in lst) == (unit < 0)
    assert abs(_rivx_budget_in(lst) - 5.0 * (12.0 - heads[2])) < 1e-4


@requires_exe(USGT_EXE)
def test_usgt_exe_riv_disu_declared_aux_irch(function_tmpdir):
    """A hand-written RIV with 'AUX IRCH' and the same file loaded and
    rewritten by FloPy give the same heads, RIV flow and CBC (with the AUX
    value saved by COMPACT BUDGET AUX)."""
    from flopy.mfusg import MfUsgRiv

    manual = _rivx_model(
        function_tmpdir / "m",
        55,
        manual="# manual RIV\n 1 55 AUX IRCH\n 1 0\n 3 12.0 5.0 8.0 7.0\n",
    )
    loaded = MfUsgRiv.load(
        str(function_tmpdir / "m" / "manual.riv"),
        _rivx_model(function_tmpdir / "l", 0, manual="x\n"),
        ext_unit_dict={},
    )
    rewritten = _rivx_model(
        function_tmpdir / "w",
        loaded.ipakcb,
        riv={
            "stress_period_data": loaded.stress_period_data.data,
            "dtype": loaded.dtype,
            "options": loaded.options,
        },
    )
    lines = (function_tmpdir / "w" / "r.riv").read_text().splitlines()
    assert lines[1].split() == ["1", "55", "AUX", "IRCH"]
    assert lines[3].split() == [
        "3",
        "12.000000",
        "5.000000e+00",
        "8.000000",
        "7.000000e+00",
    ]
    results = []
    for ml in (manual, rewritten):
        heads, lst = _rivx_run(ml)
        assert "AUXILIARY RIVER VARIABLE: IRCH" in lst
        cbf = CellBudgetFile(str(Path(ml.model_ws) / "r.cbc"), precision="single")
        assert cbf.get_times() == [1.0] and cbf.get_kstpkper() == [(0, 0)]
        q = cbf.get_data(text="RIVER LEAKAGE")[0]
        assert q.dtype.names == ("node", "q", "IRCH")
        assert q["node"].tolist() == [3] and q["IRCH"][0] == 7.0
        assert np.isfinite(q["q"][0]) and q["q"][0] > 0
        # q and h are single precision; 12 - h3 (~0.18) amplifies the head
        # rounding (~6e-8 relative) about 66 times
        np.testing.assert_allclose(q["q"][0], 5.0 * (12.0 - heads[2]), rtol=1e-5)
        results.append((heads, q["q"][0]))
    np.testing.assert_array_equal(results[0][0], results[1][0])
    assert results[0][1] == results[1][1]


# G03: EVS/RTS name-file entries kept as opaque time-series files. The fixture
# (hand-written, see test_usg_transport._g03_fixture) has analytical rates:
# recharge 0.3 on [0,4] and 0.4 on [4,10]; ET 0.1 on [0,6] and 0.02 on [6,10];
# ATS ends time steps at the series times (TEND).


def _g03x_run(ws, expect_success=True):
    from flopy.mbase import run_model

    success, buff = run_model(
        USGT_EXE, "t.nam", model_ws=str(ws), silent=True, report=True
    )
    ok = success and any("normal termination" in b.lower() for b in buff)
    assert ok == expect_success
    if not expect_success:
        return None
    lst = (Path(ws) / "t.list").read_text()
    assert "FAILED TO CONVERGE" not in lst.upper()
    steps = []
    marker = "VOLUMETRIC BUDGET FOR ENTIRE MODEL AT END OF TIME STEP"
    for block in lst.split(marker)[1:]:
        budget_in, rout = block.split("OUT:", 1)
        rec = {}
        for side, part in (("in", budget_in), ("out", rout.split("IN - OUT")[0])):
            for item in ("STORAGE", "CONSTANT HEAD", "ET", "RECHARGE"):
                m = re.search(rf"^\s+{item} =\s+(\S+)\s+{item} =\s+(\S+)", part, re.M)
                if m:
                    rec[f"cum_{side}_{item}"] = float(m.group(1))
                    rec[f"rate_{side}_{item}"] = float(m.group(2))
        rec["pct"] = float(re.search(r"PERCENT DISCREPANCY =\s+(\S+)", block).group(1))
        m = re.search(r"TOTAL TIME\s+\S+\s+\S+\s+\S+\s+(\S+)", block)
        rec["totim"] = float(m.group(1))
        steps.append(rec)
    hf = HeadUFile(str(Path(ws) / "t.hds"))
    times = hf.get_times()
    heads = np.asarray(hf.get_data(totim=times[-1])[0], dtype=float).ravel()
    assert times == [10.0] and heads.shape == (3,) and np.all(np.isfinite(heads))
    assert all(abs(r["pct"]) <= 0.01 for r in steps)
    return steps, heads


@requires_exe(USGT_EXE)
def test_usgt_exe_evsrts_control_analytic(function_tmpdir):
    """Hand-written control: ATS step ends, recharge/ET rates and cumulative
    volumes match the analytical values of the RTS/EVS series."""
    from autotest.test_usg_transport import _g03_fixture

    steps, heads = _g03x_run(_g03_fixture(function_tmpdir))
    assert [r["totim"] for r in steps] == [1.0, 3.0, 4.0, 6.0, 10.0]
    assert [r["rate_in_RECHARGE"] for r in steps] == [0.3, 0.3, 0.3, 0.4, 0.4]
    assert [r["rate_out_ET"] for r in steps] == [0.1, 0.1, 0.1, 0.1, 0.02]
    # list file values have 5 significant digits
    np.testing.assert_allclose(
        [steps[-1]["cum_in_RECHARGE"], steps[-1]["cum_out_ET"]], [3.6, 0.68], rtol=1e-5
    )
    assert heads[0] == 10.0 and np.all(heads[1:] > 10.0)


@requires_exe(USGT_EXE)
def test_usgt_exe_evsrts_sensitivity(function_tmpdir):
    """Without the series (plain RCH/EVT) the rates change; with the RCH/EVT
    options but without the EVS/RTS entries USG-T stops (reads unit 0)."""
    from autotest.test_usg_transport import (
        _G03_EVT_PLAIN,
        _G03_RCH_PLAIN,
        _g03_fixture,
    )

    plain = _g03_fixture(function_tmpdir / "plain")
    (plain / "t.rch").write_text(_G03_RCH_PLAIN)
    (plain / "t.evt").write_text(_G03_EVT_PLAIN)
    steps, _ = _g03x_run(plain)
    assert {r["rate_in_RECHARGE"] for r in steps} == {0.2}
    assert {r["rate_out_ET"] for r in steps} == {0.04}
    assert [r["totim"] for r in steps] != [1.0, 3.0, 4.0, 6.0, 10.0]
    dropped = _g03_fixture(function_tmpdir / "dropped")
    nam = [
        ln
        for ln in (dropped / "t.nam").read_text().splitlines()
        if not ln.startswith(("RTS", "EVS"))
    ]
    (dropped / "t.nam").write_text("\n".join(nam) + "\n")
    _g03x_run(dropped, expect_success=False)


@requires_exe(USGT_EXE)
@pytest.mark.parametrize("series", ["both", "rts", "evs"])
@pytest.mark.parametrize("target", ["same_ws", "new_ws"])
def test_usgt_exe_evsrts_rewrite_matches_control(function_tmpdir, series, target):
    """FloPy load (load_only without RCH/EVT) + rewrite runs like the
    hand-written control. In a new workspace the test copies t.rch/t.evt
    itself: entries skipped by load_only are not copied by FloPy."""
    from autotest.test_usg_transport import _G03_BASE, _g03_fixture

    kw = {"rts": None} if series == "evs" else {}
    if series == "rts":
        kw["evs"] = None
    control = _g03x_run(_g03_fixture(function_tmpdir / "control", **kw))
    src = _g03_fixture(function_tmpdir / "src", **kw)
    ml = MfUsg.load("t.nam", model_ws=str(src), load_only=list(_G03_BASE), check=False)
    if target == "new_ws":
        new = function_tmpdir / "new"
        ml.change_model_ws(str(new))
        ml.write_input()
        for name in ("t.rch", "t.evt"):
            if (src / name).exists():
                shutil.copy2(src / name, new / name)
    else:
        ml.write_input()
    steps, heads = _g03x_run(ml.model_ws)
    assert steps == control[0]
    np.testing.assert_allclose(heads, control[1], rtol=0, atol=1e-6)


@requires_exe(USGT_EXE)
def test_usgt_exe_evsrts_quoted_path(function_tmpdir):
    """A quoted RTS name with spaces is read by USG-T and kept by FloPy."""
    from autotest.test_usg_transport import _G03_BASE, _g03_fixture

    control = _g03x_run(_g03_fixture(function_tmpdir / "control"))
    line = "RTS               40  'ts dir/r s.rts' FLUSH"
    src = _g03_fixture(function_tmpdir / "src", rts=line)
    manual = _g03x_run(_g03_fixture(function_tmpdir / "manual", rts=line))
    assert manual[0] == control[0]
    ml = MfUsg.load("t.nam", model_ws=str(src), load_only=list(_G03_BASE), check=False)
    new = function_tmpdir / "new"
    ml.change_model_ws(str(new))
    ml.write_input()
    assert line in (new / "t.nam").read_text()
    for name in ("t.rch", "t.evt"):
        shutil.copy2(src / name, new / name)
    steps, _ = _g03x_run(new)
    assert steps == control[0]


# G06: TRTOSS stress periods. Hand-written fixtures from
# test_usg_transport._g06_model (5 cells, CHD h=10, WEL -0.5, strt 20,
# confined, ATS): heads decay to the steady line; a TRTOSS period ends one
# step after the storage rate falls below the fraction of the total rate.


def _g06x_run(ws, grid):
    from flopy.mbase import run_model

    success, buff = run_model(
        USGT_EXE, "t.nam", model_ws=str(ws), silent=True, report=True
    )
    assert success and any("normal termination" in b.lower() for b in buff)
    lst = (Path(ws) / "t.list").read_text()
    assert "FAILED TO CONVERGE" not in lst.upper()
    steps = []
    marker = "VOLUMETRIC BUDGET FOR ENTIRE MODEL AT END OF TIME STEP"
    for block in lst.split(marker)[1:]:
        head = re.match(r"\s*(\d+) IN STRESS PERIOD\s+(\d+)", block)
        kstp, kper = map(int, head.groups())
        budget_in, rout = block.split("OUT:", 1)
        rout = rout.split("IN - OUT")[0]

        def rate(part, item):
            m = re.search(rf"^\s+{item} =\s+\S+\s+{item} =\s+(\S+)", part, re.M)
            return float(m.group(1)) if m else 0.0

        sto = rate(budget_in, "STORAGE") + rate(rout, "STORAGE")
        tot = rate(budget_in, "TOTAL IN") + rate(rout, "TOTAL OUT")
        m = re.search(r"TOTAL TIME\s+\S+\s+\S+\s+\S+\s+(\S+)", block)
        pct = float(re.search(r"PERCENT DISCREPANCY =\s+(\S+)", block).group(1))
        steps.append(
            {
                "kstp": kstp,
                "kper": kper,
                "totim": float(m.group(1)),
                "sto": sto,
                "tot": tot,
                "stfrc": sto / tot,
                "pct": pct,
                "chd_in": rate(budget_in, "CONSTANT HEAD"),
                "chd_out": rate(rout, "CONSTANT HEAD"),
            }
        )
    assert all(abs(s["pct"]) <= 0.01 for s in steps)
    reached = [
        tuple(map(int, m))
        for m in re.findall(
            r"TRANSIENT REACHED STEADY-STATE CONDITION AT TIME STEP\s+(\d+) FOR "
            r"STRESS PERIOD\s+(\d+)",
            lst,
        )
    ]
    hfile = Path(ws) / "t.hds"
    if grid == "dis":
        hf = HeadFile(str(hfile))
        heads = [hf.get_data(totim=t).ravel() for t in hf.get_times()]
    else:
        hf = HeadUFile(str(hfile))
        heads = [np.asarray(hf.get_data(totim=t)[0], float) for t in hf.get_times()]
    heads = np.array(heads, dtype=float)
    assert heads.shape[1] == 5 and np.all(np.isfinite(heads))
    return {
        "steps": steps,
        "reached": reached,
        "heads": heads,
        "times": hf.get_times(),
        "kstpkper": hf.get_kstpkper(),
    }


def _g06x_sp(run, kper):
    return [s for s in run["steps"] if s["kper"] == kper]


@requires_exe(USGT_EXE)
@pytest.mark.parametrize("grid", ["dis", "disu"])
def test_usgt_exe_trtoss_criterion_and_sensitivity(function_tmpdir, grid):
    """The period ends at the step after the first one with storage
    fraction below STOFRAC (message at that step); 1e-6 ends later."""
    from autotest.test_usg_transport import _G06_E3, _g06_model

    runs = {}
    for sf in ("1.0E-3", "1.0E-6"):
        periods = [f"100000.0   1   1.0   TRTOSS {sf}"] + _G06_E3[1:]
        ws = _g06_model(function_tmpdir / sf, grid, periods=periods)
        runs[sf] = _g06x_run(ws, grid)
        sp1 = _g06x_sp(runs[sf], 1)
        n = len(sp1)
        assert sp1[-1]["kstp"] == n and runs[sf]["reached"] == [(n, 1)]
        assert sp1[-2]["stfrc"] < float(sf) <= sp1[-3]["stfrc"]
        assert sp1[-1]["totim"] < 1000.0
        assert [len(_g06x_sp(runs[sf], k)) > 0 for k in (2, 3)] == [True, True]
    e3, e6 = (_g06x_sp(runs[s], 1) for s in ("1.0E-3", "1.0E-6"))
    assert len(e6) > len(e3) and e6[-1]["totim"] > e3[-1]["totim"]


@requires_exe(USGT_EXE)
def test_usgt_exe_trtoss_without_ats_runs_like_tr(function_tmpdir):
    """Without ATS the TRTOSS period is not ended early."""
    from autotest.test_usg_transport import _G06_E3, _g06_model

    tail = _G06_E3[1:]
    first = ["100000.0  20   1.5   TRTOSS 1.0E-3", "100000.0  20   1.5   TR"]
    a, b = (
        _g06x_run(
            _g06_model(function_tmpdir / f"m{i}", ats=False, periods=[rec] + tail),
            "dis",
        )
        for i, rec in enumerate(first)
    )
    assert a["reached"] == [] and a["steps"] == b["steps"]
    np.testing.assert_array_equal(a["heads"], b["heads"])
    assert _g06x_sp(a, 1)[-1]["totim"] == 100000.0


_G06X_CASES = {
    "e3": [
        "100000.0   1   1.0   TRTOSS 1.0E-3",
        "    10.0   1   1.0   TR",
        "    10.0   1   1.0   TR",
    ],
    "all": ["100000.0   1   1.0   TRTOSS 1.0E-3"] * 3,
    "short": [
        "     5.0   1   1.0   TRTOSS 1.0E-3",
        "    10.0   1   1.0   TR",
        "    10.0   1   1.0   TR",
    ],
    "digits": [
        "100000.0   1   1.0   TRTOSS 1.23456789D-3",
        "    10.0   1   1.0   TR",
        "    10.0   1   1.0   TR",
    ],
}


@requires_exe(USGT_EXE)
@pytest.mark.parametrize("grid", ["dis", "disu"])
@pytest.mark.parametrize("flow", ["lpf", "bcf"])
@pytest.mark.parametrize("case", list(_G06X_CASES))
def test_usgt_exe_trtoss_rewrite_matches_control(function_tmpdir, grid, flow, case):
    """FloPy load + rewrite (OC skipped by load_only and kept as written)
    runs exactly like the hand-written control: steps, times, termination
    messages, budgets and heads."""
    from autotest.test_usg_transport import _g06_model, _g06_records

    periods = _G06X_CASES[case]
    control = _g06x_run(
        _g06_model(function_tmpdir / "control", grid, flow, periods), grid
    )
    ws = _g06_model(function_tmpdir / "rewrite", grid, flow, periods)
    disname = "DIS" if grid == "dis" else "DISU"
    ml = MfUsg.load(
        "t.nam",
        model_ws=str(ws),
        check=False,
        load_only=[
            disname,
            "BAS6",
            "LPF" if flow == "lpf" else "BCF6",
            "SMS",
            "CHD",
            "WEL",
        ],
    )
    ml.write_input()
    records = _g06_records(ws / ("t.dis" if grid == "dis" else "t.disu"))
    assert [r[0] for r in records] == [p.split()[3] for p in periods]
    for r, p in zip(records, periods):
        if r[0] == "TRTOSS":
            assert np.float32(r[1]) == np.float32(p.split()[4].replace("D", "E"))
    run = _g06x_run(ws, grid)
    assert run["steps"] == control["steps"] and run["reached"] == control["reached"]
    assert run["times"] == control["times"]
    np.testing.assert_allclose(run["heads"], control["heads"], rtol=0, atol=1e-6)


@requires_exe(USGT_EXE)
@pytest.mark.parametrize("grid", ["dis", "disu"])
def test_usgt_exe_trtoss_short_is_not_ss(function_tmpdir, grid):
    """A TRTOSS period that does not reach its criterion runs like TR and
    differs from SS: writing it as SS would change the heads."""
    from autotest.test_usg_transport import _G06_E3, _g06_model

    runs = {}
    for i, f in enumerate(("TRTOSS 1.0E-3", "TR", "SS")):
        periods = [f"     5.0   1   1.0   {f}"] + _G06_E3[1:]
        ws = _g06_model(function_tmpdir / f"m{i}", grid, periods=periods)
        runs[f] = _g06x_run(ws, grid)
    trtoss, tr, ss = runs.values()
    assert trtoss["reached"] == [] and trtoss["steps"] == tr["steps"]
    end1 = [k for k, kk in enumerate(trtoss["kstpkper"]) if kk[1] == 0][-1]
    end1_ss = [k for k, kk in enumerate(ss["kstpkper"]) if kk[1] == 0][-1]
    assert np.max(np.abs(trtoss["heads"][end1] - ss["heads"][end1_ss])) > 0.1


# G04: SMS options, NUMTRACK records and REAL values (hand-written models from
# test_usg_transport._g04_model). The engine echo (E15.6), option messages,
# outer iterations per step, heads and budgets of a FloPy rewrite are
# compared with the hand-written control.

_G04X_ECHO = [
    "HCLOSE",
    "HICLOSE",
    "MXITER",
    "ITER1",
    "NONLINMETH",
    "LINMETH",
    "THETA",
    "KAPPA",
    "GAMMA",
    "AMOMENTUM",
    "NUMTRACK",
    "BTOL",
    "BREDUC",
    "RES_LIM",
    "ITRUNCNEWTON",
]
_G04X_MSG = {
    "SHIFT": "SOLUTION VECTOR WILL BE SHIFTED",
    "DAMPBOT": "BOTTOM DAMPING",
    "SOLVEACTIVE": "ONLY ACTIVE NODES",
    "CUTOFF": "WITH MAXCUTOFF",
    "TRUNCNEWTON": "TRUNCATED NEWTON OPTION",
    "SIMPLE": "SIMPLE OPTION",
}


def _g04x_run(ws):
    from flopy.mbase import run_model

    ws = Path(ws)
    success, buff = run_model(
        USGT_EXE, "t.nam", model_ws=str(ws), silent=True, report=True
    )
    lst = (ws / "t.list").read_text()
    out = {
        "normal": success and any("normal termination" in b.lower() for b in buff),
        "echo": {
            k: m.group(1)
            for k in _G04X_ECHO
            if (m := re.search(rf"\({k}\)\s*=\s*(\S+)", lst))
        },
        "messages": sorted(k for k, v in _G04X_MSG.items() if v in lst),
        "outer": [int(x) for x in re.findall(r"TOTAL OF\s+(\d+)OUTER ITERATIONS", lst)],
        "failed": "FAILED TO CONVERGE" in lst.upper(),
        "pct": [float(x) for x in re.findall(r"PERCENT DISCREPANCY =\s+(\S+)", lst)],
        "times": [],
        "heads": None,
    }
    hds = ws / "t.hds"
    if out["normal"] and hds.exists():
        hf = HeadFile(str(hds)) if (ws / "t.dis").exists() else HeadUFile(str(hds))
        out["times"] = hf.get_times()
        out["heads"] = np.array(
            [np.asarray(hf.get_data(totim=t), float).ravel() for t in hf.get_times()]
        )
    return out


def _g04x_rewrite(ws):
    disname = "DIS" if (Path(ws) / "t.dis").exists() else "DISU"
    ml = MfUsg.load("t.nam", model_ws=str(ws), check=False, load_only=[disname, "SMS"])
    ml.write_input()


_G04X_XMD = " 1 0 7 2 0 0.0 0 1.0E-3\n"
# (grid, flow, laytyp, SMS text, converged): converged=False marks cases whose
# expected engine behaviour is an abort or dry cells (equivalence of
# behaviour, not a usable model)
_G04X_CASES = {
    "dis_lpf_picard_xmd_shift": (
        "dis",
        "lpf",
        1,
        " 1.0E-6 1.0E-6 100 200 1 0 1 SHIFT\n" + _G04X_XMD,
        True,
    ),
    "dis_lpf_newton_nt0": (
        "dis",
        "lpf",
        4,
        " 1.0E-6 1.0E-6 100 200 1 1 1\n 0.7 0.1 0.2 0.001 0\n" + _G04X_XMD,
        True,
    ),
    "dis_lpf_newton_nt20": (
        "dis",
        "lpf",
        4,
        " 1.0E-6 1.0E-6 100 200 1 1 1\n"
        " 0.7 0.1 0.2 0.001 20 1.0E4 0.2 100.0\n" + _G04X_XMD,
        True,
    ),
    "dis_lpf_newton_digits": (
        "dis",
        "lpf",
        4,
        " 1.23456789E-6 9.87654321D-7 100 200 1 1 1"
        " truncatednewtoncutoff 0.7654321\n"
        " 0.712345678 0.1234567 0.2345678 0.0012345678 20"
        " 1.2345678D4 0.2345678 123.45678\n"
        " 1 0 7 2 0 0.0 0 1.23456789E-3\n",
        True,
    ),
    "dis_lpf_simple_mixed": (
        "dis",
        "lpf",
        1,
        "Simple\n 1.0E-6 1.0E-6 100 200 1 0 1 Shift solveactive\n",
        True,
    ),
    "dis_lpf_simple_pcgu": (
        "dis",
        "lpf",
        1,
        "SIMPLE\n 1.0E-6 1.0E-6 100 200 1 0 2\n\n",
        True,
    ),
    "disu_lpf_picard_shift_lower": (
        "disu",
        "lpf",
        1,
        " 1.0E-6 1.0E-6 100 200 1 0 1 shift dampbot\n" + _G04X_XMD,
        True,
    ),
    "disu_lpf_newton_nt0": (
        "disu",
        "lpf",
        4,
        " 1.0E-6 1.0E-6 100 200 1 1 1\n 0.7 0.1 0.2 0.001 0\n" + _G04X_XMD,
        True,
    ),
    "dis_bcf_picard_digits": (
        "dis",
        "bcf",
        0,
        " 1.23456789E-6 1.0E-6 100 200 1 0 1"
        " SOLVEACTIVE\n 1 0 7 2 0 0.0 0 1.23456789E-3\n",
        True,
    ),
    "dis_lpf_pcgu_shift_dry": (
        "dis",
        "lpf",
        1,
        " 1.0E-3 1.0E-2 100 3 1 0 2 SHIFT\n BCGS 1 0 0 1.0E-2 1.0\n",
        False,
    ),
    "dis_lpf_pcgu_shift_abort": (
        "dis",
        "lpf",
        1,
        " 1.0E-4 1.0E-3 100 5 1 0 2 SHIFT\n CG 1 0 0 1.0E-3 1.0\n",
        False,
    ),
}


@requires_exe(USGT_EXE)
@pytest.mark.parametrize("case", list(_G04X_CASES))
def test_usgt_exe_sms_rewrite_matches_control(function_tmpdir, case):
    from autotest.test_usg_transport import _g04_model

    grid, flow, laytyp, text, converged = _G04X_CASES[case]
    ws0 = _g04_model(function_tmpdir / "control", text, grid, flow, laytyp)
    control = _g04x_run(ws0)
    ws = _g04_model(function_tmpdir / "rewrite", text, grid, flow, laytyp)
    _g04x_rewrite(ws)
    run = _g04x_run(ws)
    for key in ("normal", "echo", "messages", "outer", "failed", "pct", "times"):
        assert run[key] == control[key], key
    if converged:
        assert control["normal"] and not control["failed"]
        assert all(abs(p) <= 0.01 for p in control["pct"])
        assert np.all(np.abs(control["heads"]) < 1.0e29)
        np.testing.assert_array_equal(run["heads"], control["heads"])
    elif control["heads"] is not None:  # same dry cells (1e30)
        np.testing.assert_array_equal(run["heads"], control["heads"])


@requires_exe(USGT_EXE)
def test_usgt_exe_sms_shift_is_read(function_tmpdir):
    """SHIFT is read by the engine and changes this PCGU run (engine
    behaviour observed with CHD: dry cells), unlike DAMPBOT or no option;
    a rewrite must therefore keep SHIFT."""
    from autotest.test_usg_transport import _g04_model

    runs = {}
    for opt in ("", " SHIFT", " DAMPBOT"):
        text = f" 1.0E-3 1.0E-2 100 3 1 0 2{opt}\n BCGS 1 0 0 1.0E-2 1.0\n"
        runs[opt] = _g04x_run(_g04_model(function_tmpdir / f"m{len(runs)}", text))
    assert runs[" SHIFT"]["messages"] == ["SHIFT"]
    assert runs[" SHIFT"]["outer"] != runs[""]["outer"]
    assert runs[" DAMPBOT"]["outer"] == runs[""]["outer"]
    assert np.max(np.abs(runs[" SHIFT"]["heads"])) >= 1.0e29  # dry cells
    assert np.max(np.abs(runs[""]["heads"])) < 1.0e29


@requires_exe(USGT_EXE)
def test_usgt_exe_sms_option_separators(function_tmpdir):
    """Record 1b with commas, tabs, mixed case and two cutoffs (both read by
    USG-T, in order): the FloPy rewrite reads the same."""
    from autotest.test_usg_transport import _g04_model

    text = (
        " 1.0E-6 1.0E-6 100 200 1 1 1 SHIFT,dampbot\tsolveactive "
        "TRUNCATEDNEWTONCUTOFF 0.9,truncatednewtoncutoff 0.7654321\n"
        " 0.7 0.1 0.2 0.001 0\n" + _G04X_XMD
    )
    runs = []
    for kind in ("control", "rewrite"):
        ws = _g04_model(function_tmpdir / kind, text, "dis", "lpf", 4)
        if kind == "rewrite":
            _g04x_rewrite(ws)
        run = _g04x_run(ws)
        run["cutoffs"] = re.findall(
            r"MAXCUTOFF\s*=\s*(\S+)", (ws / "t.list").read_text()
        )
        runs.append(run)
    control, rewrite = runs
    assert control["normal"] and not control["failed"]
    assert control["cutoffs"] == ["0.90000E+00", "0.76543E+00"]
    assert set(control["messages"]) >= {"SHIFT", "DAMPBOT", "SOLVEACTIVE", "CUTOFF"}
    for key in ("echo", "messages", "outer", "pct", "times", "cutoffs"):
        assert rewrite[key] == control[key], key
    np.testing.assert_array_equal(rewrite["heads"], control["heads"])


# G05: OC with ATS. Hand-written controls against a FloPy load -> write of the
# OC file only (the other files are the same); numeric ATS records are
# written as words, so the list-file echo of the ATS input differs.
_G05X_HDR = "HEAD SAVE UNIT 51\nCOMPACT BUDGET AUX\n"
_G05X_ACT = "  save head\n  save budget\n  print budget\n"
_G05X_NUM = " 0 1 1 1\n 0 0 1 0\n"
_G05X_P2 = "     1.0   1   1.0   SS\n    10.0   5   1.0   TR\n"
_G05X_P3 = _G05X_P2 + "    10.0   5   1.0   TR\n"
_G05X_LONG = "     1.0   1   1.0   SS\n  1000.0   5   1.0   TR\n"


def _g05x_words(first, timot="", p2=""):
    return (
        f"# OC\n{first}\n{timot}"
        + _G05X_HDR
        + "period 1\n"
        + _G05X_ACT
        + "period 2\n"
        + p2
        + _G05X_ACT
    )


def _g05x_num(header, timot, rec1, rec2, flags=_G05X_NUM):
    return f"# OC\n 0 0 51 0 ATSA{header}\n{timot}{rec1}\n{flags}{rec2}\n{flags}"


_G05X_T3 = "3.5 6.25 11.0\n"
_G05X_R1 = " 1.0 1.0E-5 10.0 2.0 5.0"
_G05X_R2 = " 0.5 1.0E-5 2.0 2.0 5.0"
_G05X_STEPS = "".join(
    f"period {k} step {s}\n" + _G05X_ACT
    for k, n in ((1, 1), (2, 5))
    for s in range(1, n + 1)
)
# name: (grid, periods, control OC[, hand-written words control for the rewrite])
_G05X_CASES = {
    "dis_noats_words": ("dis", _G05X_P2, "# OC\n" + _G05X_HDR + _G05X_STEPS),
    "dis_noats_numeric": ("dis", _G05X_P2, "# OC\n 0 0 51 0\n" + _G05X_NUM * 6),
    "disu_noats_words": ("disu", _G05X_P2, "# OC\n" + _G05X_HDR + _G05X_STEPS),
    "dis_ats_dec_timot": (
        "dis",
        _G05X_P2,
        _g05x_words("ATSA NPTIMES 3", _G05X_T3, "  DELTAT 0.5\n  TMAXAT 2.0\n"),
    ),
    "dis_ats_int_timot": (
        "dis",
        _G05X_P2,
        _g05x_words("ATSA NPTIMES 3", "3 6 11\n", "  DELTAT 0.5\n  TMAXAT 2.0\n"),
    ),
    "dis_ats_timot_forms": (
        "dis",
        _G05X_P2,
        _g05x_words(
            "atsa nptimes 4",
            "3.5E0,\n\n6.25d0 , 8.5\t1.1D1\n",
            "  DELTAT 0.5\n  TMAXAT 2.0\n",
        ),
    ),
    "disu_ats_dec_timot": (
        "disu",
        _G05X_P2,
        _g05x_words("ATSA NPTIMES 3", _G05X_T3, "  DELTAT 0.5\n  TMAXAT 2.0\n"),
    ),
    "dis_ats_param_digits": (
        "dis",
        _G05X_P2,
        _g05x_words(
            "ATSA", "", "  DELTAT 0.123456789\n  TMAXAT 1.7654321\n  TADJAT 1.33333\n"
        ),
    ),
    # ATS actions separated by a comma, a tab and blanks around a comma
    "dis_ats_comma_params": (
        "dis",
        _G05X_P2,
        _g05x_words(
            "ATSA",
            "",
            "  DELTAT,0.123456789\n  TMAXAT\t1.7654321\n  tadjat , 1.33333\n",
        ),
    ),
    "dis_ats_npstps": (
        "dis",
        _G05X_P2,
        _g05x_words("ATSA NPSTPS 2", "", "  DELTAT 0.5\n  TMAXAT 2.0\n"),
    ),
    "dis_long_timot_npstps": (
        "dis",
        _G05X_LONG,
        _g05x_words(
            "ATSA NPTIMES 60 NPSTPS 4",
            "".join(
                " ".join(str(10 * (6 * r + c + 1)) for c in range(6)) + "\n"
                for r in range(10)
            ),
        ),
    ),
    "dis_old_flopy_inline_timot": (
        "dis",
        _G05X_LONG,
        "# OC\nATSA NPTIMES 60 "
        + " ".join(str(10 * (i + 1)) for i in range(60))
        + " NPSTPS 4 \n"
        + " ".join(str(10 * (i + 1)) for i in range(60))
        + " \n"
        + _G05X_HDR
        + "period 1\n"
        + _G05X_ACT
        + "period 2\n"
        + _G05X_ACT,
    ),
    "dis_ats_ddref": (
        "dis",
        _G05X_P2,
        _g05x_words("ATSA", "", "  DELTAT 0.5\n").replace(
            "period 2\n", "period 2 ddreference\n"
        ),
    ),
    "dis_numeric_ats_dec": (
        "dis",
        _G05X_P2,
        _g05x_num(
            " NPTIMES 3", _G05X_T3, _G05X_R1 + " 3 1000 0 0", _G05X_R2 + " 3 1000 0 0"
        ),
    ),
    "disu_numeric_ats_dec": (
        "disu",
        _G05X_P2,
        _g05x_num(
            " NPTIMES 3", _G05X_T3, _G05X_R1 + " 3 1000 0 0", _G05X_R2 + " 3 1000 0 0"
        ),
    ),
    "dis_numeric_npstps_fields": (
        "dis",
        _G05X_P2,
        _g05x_num("", "", _G05X_R1 + " 0 2 0 0", _G05X_R2 + " 0 2 0 0"),
    ),
    "dis_numeric_nptimes_no_atsa": (
        "dis",
        _G05X_P2,
        "# OC\n 0 0 51 0 NPTIMES 3\n3 6 11\n"
        + _G05X_R1
        + " 3 1000 0 0\n"
        + _G05X_NUM
        + _G05X_R2
        + " 3 1000 0 0\n"
        + _G05X_NUM,
    ),
    "dis_numeric_digits": (
        "dis",
        _G05X_P2,
        _g05x_num(
            "",
            "",
            _G05X_R1 + " 0 1000 0 0",
            " 0.123456789 1.0E-5 1.7654321 1.33333 5.0 0 1000 0 0",
        ),
    ),
    "dis_numeric_limits": (
        "dis",
        _G05X_P2,
        _g05x_num("", "", _G05X_R1 + " 0 1000 0 0", " 0.0 0.0 0.0 0.5 0.9 0 1000 0 0"),
    ),
    "dis_numeric_ss_perlen5": (
        "dis",
        "     5.0   1   1.0   SS\n    10.0   5   1.0   TR\n",
        _g05x_num(
            " NPTIMES 2",
            "8.5 15.0\n",
            _G05X_R1 + " 2 1000 0 0",
            _G05X_R2 + " 2 1000 0 0",
        ),
    ),
    "dis_numeric_ss_perlen01": (
        "dis",
        "     0.1   1   1.0   SS\n    10.0   5   1.0   TR\n",
        _g05x_num("", "", _G05X_R1 + " 0 1000 0 0", _G05X_R2 + " 0 1000 0 0"),
    ),
    "dis_numeric_three_periods": (
        "dis",
        _G05X_P3,
        _g05x_num(
            " NPTIMES 4",
            "3.5 6.25 11.0 17.5\n",
            _G05X_R1 + " 4 1000 0 0",
            _G05X_R2 + " 4 1000 0 0",
        )
        + " 1.0 1.0E-5 3.0 2.0 5.0 4 1000 0 0\n"
        + _G05X_NUM,
    ),
    "dis_numeric_incode_layers": (
        "dis",
        _G05X_P2,
        _g05x_num(
            "",
            "",
            _G05X_R1 + " 0 1000 0 0",
            _G05X_R2 + " 0 1000 0 0",
            flags=" 1 1 1 1\n 0 0 1 0\n",
        ),
    ),
    # NPSTPS 0 in every record is written as no NPSTPS; the rewrite is checked
    # against a hand-written words control (a run of the numeric file would
    # evaluate MOD(KSTP,0), whose result depends on the platform)
    "dis_numeric_npstps0": (
        "dis",
        _G05X_P2,
        _g05x_num(" NPTIMES 3", _G05X_T3, _G05X_R1 + " 3 0 0 0", _G05X_R2 + " 3 0 0 0"),
        "# OC\nATSA NPTIMES 3\n" + _G05X_T3 + "HEAD SAVE UNIT 51\nperiod 1\n"
        "  DELTAT 1.0\n  TMINAT 1.0E-5\n  TMAXAT 10.0\n  TADJAT 2.0\n  TCUTAT 5.0\n"
        "  print budget\n  save budget\n  save head\nperiod 2\n  DELTAT 0.5\n"
        "  TMINAT 1.0E-5\n  TMAXAT 2.0\n  TADJAT 2.0\n  TCUTAT 5.0\n"
        "  print budget\n  save budget\n  save head\n",
    ),
}


def _g05x_run(ws):
    """Run the model in ws; outputs the OC asked for, read from new files."""
    from flopy.mbase import run_model

    ws = Path(ws)
    success, buff = run_model(
        USGT_EXE, "t.nam", model_ws=str(ws), silent=True, report=True
    )
    lst = (ws / "t.list").read_text()
    out = {
        "normal": success and any("normal termination" in b.lower() for b in buff),
        "failed": "FAILED TO CONVERGE" in lst.upper(),
        "pct": [float(x) for x in re.findall(r"PERCENT DISCREPANCY =\s+(\S+)", lst)],
        "budget_at": re.findall(
            r"VOLUMETRIC BUDGET FOR ENTIRE MODEL AT END OF "
            r"TIME STEP\s+(\d+) IN STRESS PERIOD\s+(\d+)",
            lst,
        ),
        "delt": re.findall(r"TIME STEP LENGTH\s+\S+\s+\S+\s+\S+\s+(\S+)", lst),
        "hds": [],
        "heads": None,
        "cbc": [],
    }
    hds = ws / "t.hds"
    if out["normal"] and hds.exists():
        hf = HeadFile(str(hds)) if (ws / "t.dis").exists() else HeadUFile(str(hds))
        rec = hf.recordarray
        out["hds"] = [(int(r["kstp"]), int(r["kper"]), float(r["totim"])) for r in rec]
        out["heads"] = np.array(
            [np.asarray(hf.get_data(idx=i), float).ravel() for i in range(len(rec))]
        )
    if out["normal"] and (ws / "t.cbc").exists():
        cf = CellBudgetFile(str(ws / "t.cbc"))
        out["cbc"] = sorted(set(zip(cf.get_kstpkper(), cf.get_times())))
    return out


@pytest.mark.parametrize("case", list(_G05X_CASES))
def test_usgt_exe_oc_ats_rewrite_matches_control(function_tmpdir, case):
    """A FloPy load -> write of the OC file gives the same output times
    (heads, budgets, list file), time steps, budgets, convergence and heads
    as the hand-written control."""
    from autotest.test_usg_transport import _g05_dis, _g05_model

    grid, periods, text, *words = _G05X_CASES[case]
    ws0 = _g05_model(
        function_tmpdir / "control", words[0] if words else text, grid, periods
    )
    control = _g05x_run(ws0)
    assert control["normal"] and not control["failed"] and control["hds"]
    assert np.all(np.isfinite(control["heads"]))
    ws = _g05_model(function_tmpdir / "rewrite", text, grid, periods)
    oc = MfUsgOc.load(str(ws / "t.oc"), _g05_dis(ws))
    oc.fn_path = str(ws / "t.oc")
    oc.write_file()
    assert "{" not in (ws / "t.oc").read_text()
    rewrite = _g05x_run(ws)
    for key in ("normal", "failed", "pct", "budget_at", "delt", "hds", "cbc"):
        assert rewrite[key] == control[key], key
    np.testing.assert_array_equal(rewrite["heads"], control["heads"])


def test_usgt_exe_oc_ats_contract(function_tmpdir):
    """The engine behaviour the G05 rules rely on: decimal output times are
    reached exactly with variable steps; NPTIMES as first word stops USG-T;
    the word form applies neither the ATS1RP limits nor the TIMOT limit of
    DELTAT."""
    from autotest.test_usg_transport import _g05_model

    def run(name, text):
        return _g05x_run(_g05_model(function_tmpdir / name, text, "dis", _G05X_P2))

    dec = run("dec", _G05X_CASES["dis_ats_dec_timot"][2])
    assert [t for _, _, t in dec["hds"]] == [1.0, 3.5, 6.25, 11.0]
    assert len(set(dec["delt"])) > 3
    assert not run("first", _g05x_words("NPTIMES 3", "3 6 11\n"))["normal"]
    num = run(
        "num",
        _g05x_num(
            " NPTIMES 3",
            _G05X_T3,
            _G05X_R1 + " 3 1000 0 0",
            " 5.0 1.0E-5 10.0 2.0 5.0 3 1000 0 0",
        ),
    )
    words = run(
        "words",
        _g05x_words("ATSA NPTIMES 3", _G05X_T3, "  DELTAT 5.0\n  TMAXAT 10.0\n"),
    )
    assert num["hds"] != words["hds"]
    assert not run(
        "limits", _g05x_words("ATSA", "", "  DELTAT 0.5\n  TADJAT 0.5\n  TCUTAT 0.9\n")
    )["normal"]


# G01 (partial): SEEPELEV and RTS zones with USG-T. Hand-written controls
# (test_usg_transport._g01_seep_text on the 5-cell _g05_model; a one-zone RTS
# variant of _g03_fixture) against the RCH file written by FloPy, after
# MfUsgRch.load or built with the MfUsgRch API.


def _g01x_run(ws):
    """_g05x_run plus the RECHARGE list-file lines and RECHARGE cbc records."""
    out = _g05x_run(ws)
    lst = (Path(ws) / "t.list").read_text()
    out["rch"] = re.findall(r"^\s+RECHARGE =\s+(\S+)\s+RECHARGE =\s+(\S+)", lst, re.M)
    cbc = CellBudgetFile(str(Path(ws) / "t.cbc"))
    out["rch_cbc"] = np.array(
        [np.asarray(r, float) for r in cbc.get_data(text="RECHARGE")]
    )
    return out


def _g01x_same(a, b):
    for key in ("normal", "failed", "pct", "budget_at", "delt", "hds", "cbc", "rch"):
        assert a[key] == b[key], key
    np.testing.assert_array_equal(a["heads"], b["heads"])
    np.testing.assert_array_equal(a["rch_cbc"], b["rch_cbc"])


@requires_exe(USGT_EXE)
@pytest.mark.parametrize("grid", ["dis", "disu"])
@pytest.mark.parametrize("source", ["load_write", "api"])
def test_usgt_exe_rch_seepelev_matches_control(function_tmpdir, grid, source):
    """The RCH file FloPy writes with SELEV runs like the hand-written
    control: heads, RECHARGE budget (list and cbc), steps and convergence.
    Period 2 of the control repeats SELEV (INSELEV 1) or reuses it
    (INSELEV -1); both controls give the same results."""
    from autotest.test_usg_transport import (
        _G01_RECH,
        _G01_SELEV,
        _g01_seep_model,
        _g01_seep_text,
        _g05_dis,
    )

    control = _g01x_run(
        _g01_seep_model(function_tmpdir / "control", _g01_seep_text(), grid)
    )
    assert control["normal"] and not control["failed"] and control["hds"]
    reuse = _g01x_run(
        _g01_seep_model(function_tmpdir / "reuse", _g01_seep_text(p2="-1"), grid)
    )
    _g01x_same(reuse, control)
    ws = _g01_seep_model(function_tmpdir / source, _g01_seep_text(), grid)
    if source == "load_write":
        rch = MfUsgRch.load(str(ws / "t.rch"), _g05_dis(ws), check=False)
    else:

        def arr(values):
            return np.array([float(v) for v in values]).reshape(1, 5)

        s1, s3 = (arr(s) for s in _G01_SELEV)
        rch = MfUsgRch(
            _g05_dis(ws),
            nrchop=3,
            ipakcb=50,
            rech={0: arr(_G01_RECH[0]), 1: arr(_G01_RECH[1])},
            seepelev=1,
            selev={0: s1, 1: s1, 2: s3},
        )
    rch.fn_path = str(ws / "t.rch")
    rch.write_file(check=False)
    assert (ws / "t.rch").read_text() != _g01_seep_text()
    _g01x_same(_g01x_run(ws), control)


@requires_exe(USGT_EXE)
@pytest.mark.parametrize("grid", ["dis", "disu"])
def test_usgt_exe_rch_seepelev_sensitivity(function_tmpdir, grid):
    """SELEV changes the results: one period 1 value changed (also used in
    period 2) or the option removed give other heads and recharge."""
    from autotest.test_usg_transport import _G01_SELEV, _g01_seep_model, _g01_seep_text

    def run(name, text):
        out = _g01x_run(_g01_seep_model(function_tmpdir / name, text, grid))
        assert out["normal"] and not out["failed"]
        return out

    control = run("control", _g01_seep_text())
    changed = list(_G01_SELEV[0])
    changed[2] = "10.7"
    pert = run("pert", _g01_seep_text(selev1=changed))
    plain = run("plain", _g01_seep_text(seepelev=False))
    for other in (pert, plain):
        assert other["hds"] == control["hds"]
        for i in (0, 1, 5):  # period 1, first and last step of period 2
            assert not np.array_equal(other["heads"][i], control["heads"][i])
        assert other["rch"][0] != control["rch"][0]
    # without SELEV all recharge enters (cell 1 is CHD): 8 in period 1
    assert float(plain["rch"][0][1]) == 8.0 > float(control["rch"][0][1])


_G01X_RTS = "0.0 4.0 1.0 1.0e-3\n4.0 10.0 1.0 3.0e-3\n10.0 1.0e6 1.0 0.0\n"


def _g01x_rts(ws, zones="1 0 1"):
    """_g03_fixture without EVS, RCH with RTS 1 and one-zone series."""
    from autotest.test_usg_transport import _g03_fixture

    ws = _g03_fixture(ws, evs=None)
    (ws / "t.rch").write_text(
        "# RCH\n 1 0 RTS 1\n 0 INRCHZONES 1\nCONSTANT 0.0\n"
        f"INTERNAL 1 (FREE) -1\n {zones}\n"
    )
    (ws / "ts" / "rech.rts").write_text(_G01X_RTS)
    return ws


@requires_exe(USGT_EXE)
@pytest.mark.parametrize("source", ["load_write", "api"])
def test_usgt_exe_rch_rts_zones_match_control(function_tmpdir, source):
    """One RTS zone on node 3 (node 1 is CHD): recharge 0.1 then 0.3. The
    RCH file FloPy writes (RTS 1, IZNRCH 1 0 1) runs like the control;
    other zones give other rates."""
    control = _g03x_run(_g01x_rts(function_tmpdir / "control"))
    assert {(r["totim"] <= 4.0, r["rate_in_RECHARGE"]) for r in control[0]} == {
        (True, 0.1),
        (False, 0.3),
    }
    other = _g03x_run(_g01x_rts(function_tmpdir / "zones", "0 1 1"))
    assert other[0] != control[0]
    ws = _g01x_rts(function_tmpdir / source)
    ml = MfUsg.load("t.nam", model_ws=str(ws), check=False, load_only=["DISU"])
    if source == "load_write":
        rch = MfUsgRch.load(str(ws / "t.rch"), ml, check=False)
    else:
        rch = MfUsgRch(
            ml, nrchop=1, rech=0.0, mxrtzones=1, iznrch=np.array([[1, 0, 1]])
        )
    rch.fn_path = str(ws / "t.rch")
    rch.write_file(check=False)
    lines = (ws / "t.rch").read_text().splitlines()
    assert lines[1].split() == ["1", "0", "RTS", "1"]
    steps, heads = _g03x_run(ws)
    assert steps == control[0]
    np.testing.assert_array_equal(heads, control[1])


# VORONOI-USG: DISU arrays from flopy.utils.voronoi.get_disu5_gridprops run
# with USG-T on irregular voronoi grids of a 1000 x 500 rectangle
# (test_usg_transport._vor_rect). Tolerances (HCLOSE 1e-9, heads and budgets
# saved in single precision, DISU and CHD values written as float32):
# heads 1e-4, flows rel 1e-4. MfUsgDisU converts ja to one-based in place,
# so each model gets a copy.


def _vorx_run(ws, g, chd, hk, vka, strt=95.0):
    """Steady LPF (LAYCON 0) model on the DISU arrays g; heads per node and
    the cell budget file."""
    from flopy.mfusg import MfUsgDisU

    m = MfUsg(modelname="v", model_ws=str(ws), exe_name=USGT_EXE, structured=False)
    ja = g["ja"].copy()
    MfUsgDisU(m, **{**g, "ja": ja}, nper=1, perlen=1.0, nstp=1, steady=True)
    assert ja[0] == 1 and g["ja"][0] == 0
    MfUsgBas(m, ibound=1, strt=strt, structured=False)
    MfUsgLpf(m, laytyp=0, hk=hk, vka=vka, ipakcb=50)
    MfUsgSms(m, hclose=1e-9, hiclose=1e-10, mxiter=50, iter1=1000, linmeth=1)
    MfUsgOc(m, stress_period_data={(0, 0): ["save head", "save budget"]})
    ModflowChd(m, stress_period_data={0: [[nd, h, h] for nd, h in chd]})
    m.write_input()
    success, buff = m.run_model(silent=True, report=True)
    assert success and any("normal termination" in b.lower() for b in buff)
    assert "FAILED TO CONVERGE" not in (Path(ws) / "v.list").read_text().upper()
    heads = np.concatenate(HeadUFile(str(Path(ws) / "v.hds")).get_data())
    return heads, CellBudgetFile(str(Path(ws) / "v.cbc"))


def _vorx_chd(cbc, nodes):
    """Net CONSTANT HEAD flow per node (positive into the aquifer)."""
    rec = cbc.get_data(text="CONSTANT HEAD")[0]
    q = np.zeros(nodes)
    q[rec["node"] - 1] = rec["q"]
    return q


def _vorx_x(vor, x):
    return np.flatnonzero(np.abs(vor.points[:, 0] - x) <= 1e-6)


@requires_pkg("scipy")
@requires_exe("triangle", USGT_EXE)
def test_usgt_exe_voronoi_disu5_linear_flow(tmp_path):
    """CHD on x=0 and x=1000, no-flow sides: the linear head and the inflow
    K b W dh / L are exact on a voronoi grid. Horizontal FAHL given as a
    width instead of an area changes the inflow."""
    from autotest.test_usg_transport import _vor_rect

    _, vor = _vor_rect(tmp_path / "mesh")
    ncpl = vor.ncpl
    g = vor.get_disu5_gridprops(20.0, [0.0])
    left, right = _vorx_x(vor, 0.0), _vorx_x(vor, 1000.0)
    assert len(left) >= 2 and len(right) >= 2
    chd = [(n, 100.0) for n in left] + [(n, 90.0) for n in right]
    heads, cbc = _vorx_run(tmp_path / "e1", g, chd, 5.0, 5.0)
    exact = 100.0 - 10.0 * vor.points[:, 0] / 1000.0
    assert np.max(np.abs(heads - exact)) <= 1e-4
    q = _vorx_chd(cbc, ncpl)
    assert q[left].sum() == pytest.approx(5.0 * 20.0 * 500.0 * 10.0 / 1000.0, rel=1e-4)
    assert q[right].sum() == pytest.approx(-q[left].sum(), rel=1e-4)
    width = {**g, "fahl": np.where(g["ivc"] == 0, g["fahl"] / 20.0, g["fahl"])}
    _, cbcw = _vorx_run(tmp_path / "e1_width", width, chd, 5.0, 5.0)
    assert _vorx_chd(cbcw, ncpl)[left].sum() == pytest.approx(25.0, rel=1e-4)


@requires_pkg("scipy", "shapely")
@requires_exe("triangle", USGT_EXE)
def test_usgt_exe_voronoi_disu5_vertical_series(tmp_path):
    """Three layers of variable thickness; layers 1 and 3 all CHD (50, 30),
    layer 2 with HK 0: every column is a series of vertical conductances
    A / (b1/(2 Kv1) + b2/(2 Kv2)) (USG-T computes no flow between two CHD
    cells). The cell areas come from shapely."""
    from shapely.geometry import Polygon

    from autotest.test_usg_transport import _vor_rect

    _, vor = _vor_rect(tmp_path / "mesh", seed=11)
    ncpl = vor.ncpl
    x, y = vor.points[:, 0], vor.points[:, 1]
    top = 40.0 + 5.0 * np.sin(x / 150.0)
    bot1 = 25.0 + 3.0 * np.cos(y / 80.0)
    bot2 = 10.0 + 2.0 * np.sin(x / 70.0)
    g = vor.get_disu5_gridprops(top, [bot1, bot2, 0.0])
    chd = [(n, 50.0) for n in range(ncpl)] + [(2 * ncpl + n, 30.0) for n in range(ncpl)]
    kv = (0.1, 0.3, 0.2)
    heads, cbc = _vorx_run(
        tmp_path / "e2",
        g,
        chd,
        [np.full(ncpl, h) for h in (1.0, 0.0, 2.0)],
        [np.full(ncpl, k) for k in kv],
        strt=40.0,
    )
    area = np.array([Polygon(np.asarray(vor.verts)[iv]).area for iv in vor.iverts])
    b1, b2, b3 = top - bot1, bot1 - bot2, bot2
    r12 = 0.5 * b1 / kv[0] + 0.5 * b2 / kv[1]
    r23 = 0.5 * b2 / kv[1] + 0.5 * b3 / kv[2]
    expected = area * 20.0 / (r12 + r23)
    q = _vorx_chd(cbc, 3 * ncpl)
    np.testing.assert_allclose(q[:ncpl], expected, rtol=1e-4)
    np.testing.assert_allclose(-q[2 * ncpl :], expected, rtol=1e-4)
    np.testing.assert_allclose(
        heads[ncpl : 2 * ncpl], 50.0 - 20.0 * r12 / (r12 + r23), atol=1e-4
    )


def _vorx_reference(vor, top, botm, hk, vka, chd):
    """Independent numpy solution of the USG-T LPF equations (LAYCON 0,
    LAYAVG 0: CBCK12/CBCKISO harmonic mean weighted by CL1 and CL2;
    SGWF2LPFU1VCOND for vertical flow) with the geometry from shapely: shared
    boundaries and distances from the cell centers to the face lines."""
    from shapely.geometry import LineString, Point, Polygon
    from shapely.ops import linemerge

    polys = [Polygon(np.asarray(vor.verts)[iv]) for iv in vor.iverts]
    ncpl, nlay = len(polys), len(botm)
    tops = np.vstack([top, np.array(botm[:-1])])
    thk = tops - np.array(botm)
    nodes = nlay * ncpl
    cond = np.zeros((nodes, nodes))
    pts = vor.points
    for i in range(ncpl):
        for j in range(i + 1, ncpl):
            if not polys[i].intersects(polys[j]):
                continue
            shared = polys[i].boundary.intersection(polys[j].boundary)
            lines = [
                s
                for s in getattr(shared, "geoms", [shared])
                if s.geom_type == "LineString"
            ]
            if not lines:
                continue
            face = linemerge(lines)
            if face.length <= 1e-9:
                continue
            (x1, y1), (x2, y2) = face.coords[0], face.coords[-1]
            line = LineString(
                [
                    (x1 - 1e4 * (x2 - x1), y1 - 1e4 * (y2 - y1)),
                    (x1 + 1e4 * (x2 - x1), y1 + 1e4 * (y2 - y1)),
                ]
            )
            cl1, cl2 = line.distance(Point(pts[i])), line.distance(Point(pts[j]))
            for k in range(nlay):
                t1, t2 = hk[k] * thk[k, i], hk[k] * thk[k, j]
                n, m = k * ncpl + i, k * ncpl + j
                cond[n, m] = cond[m, n] = face.length * t1 * t2 / (t1 * cl2 + t2 * cl1)
    for i in range(ncpl):
        for k in range(nlay - 1):
            c = polys[i].area / (
                0.5 * thk[k, i] / vka[k] + 0.5 * thk[k + 1, i] / vka[k + 1]
            )
            n, m = k * ncpl + i, (k + 1) * ncpl + i
            cond[n, m] = cond[m, n] = c
    mat = np.diag(cond.sum(axis=1)) - cond
    rhs = np.zeros(nodes)
    for n, h in chd:
        mat[n, :] = 0.0
        mat[n, n] = 1.0
        rhs[n] = h
    return np.linalg.solve(mat, rhs)


@requires_pkg("scipy", "shapely")
@requires_exe("triangle", USGT_EXE)
def test_usgt_exe_voronoi_disu5_two_layers_reference(tmp_path):
    """Two layers of variable thickness with horizontal and vertical flow:
    USG-T heads equal the independent reference within 1e-4."""
    from autotest.test_usg_transport import _vor_rect

    _, vor = _vor_rect(tmp_path / "mesh", seed=3, n=25, max_area=8000.0)
    ncpl = vor.ncpl
    x, y = vor.points[:, 0], vor.points[:, 1]
    top = 40.0 + 0.01 * x + 4.0 * np.sin(y / 60.0)
    bot1 = 18.0 - 0.005 * x + 2.0 * np.cos(x / 90.0)
    botm = [bot1, np.zeros(ncpl)]
    hk, vka = [3.0, 0.5], [0.05, 0.2]
    chd = [(n, 100.0) for n in _vorx_x(vor, 0.0)] + [
        (ncpl + n, 90.0) for n in _vorx_x(vor, 1000.0)
    ]
    g = vor.get_disu5_gridprops(top, botm)
    heads, _ = _vorx_run(
        tmp_path / "e3",
        g,
        chd,
        [np.full(ncpl, h) for h in hk],
        [np.full(ncpl, k) for k in vka],
    )
    ref = _vorx_reference(vor, top, botm, hk, vka, chd)
    assert np.max(np.abs(heads - ref)) <= 1e-4
    assert (
        np.ptp(heads[:ncpl]) > 1.0 and np.max(np.abs(heads[:ncpl] - heads[ncpl:])) > 0.1
    )
