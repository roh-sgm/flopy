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

import numpy as np
import pytest
from modflow_devtools.markers import requires_exe

from flopy.mfusg import (
    MfUsg,
    MfUsgBas,
    MfUsgBct,
    MfUsgEvt,
    MfUsgLpf,
    MfUsgOc,
    MfUsgPcb,
    MfUsgSms,
    MfUsgTib,
)
from flopy.modflow import ModflowChd, ModflowDis
from flopy.utils import HeadFile, MfusgListBudget, MfusgTransportListBudget

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
