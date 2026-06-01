import io
import os
import shutil
from pathlib import Path

import numpy as np
import pytest
from flaky import flaky
from modflow_devtools.markers import requires_exe

from autotest.conftest import get_example_data_path
from flopy.mfusg import (
    MfUsg,
    MfUsgBas,
    MfUsgBcf,
    MfUsgBct,
    MfUsgCln,
    MfUsgDdf,
    MfUsgDis,
    MfUsgDisU,
    MfUsgDpf,
    MfUsgDpt,
    MfUsgDrn,
    MfUsgDrt,
    MfUsgEvt,
    MfUsgGhb,
    MfUsgGnc,
    MfUsgGsf,
    MfUsgLak,
    MfUsgLpf,
    MfUsgMdt,
    MfUsgOc,
    MfUsgPcb,
    MfUsgQrt,
    MfUsgRch,
    MfUsgSgb,
    MfUsgSms,
    MfUsgTvm,
    MfUsgWel,
)
from flopy.modflow import (
    #    ModflowBas,
    #    ModflowDis,
    ModflowChd,
    ModflowFhb,
    ModflowGhb,
)
from flopy.utils import Util2d, Util3d

# USG-Transport 2.7 executable for the optional real-model run tests. Resolved
# via the USGT_EXE environment variable (executable name on PATH or an absolute
# path), defaulting to "mfusg_gsi". @requires_exe skips these tests cleanly when
# the executable cannot be resolved, so the default suite stays light.
USGT_EXE = os.environ.get("USGT_EXE", "mfusg_gsi")


def test_mfusgcln():
    ml = MfUsg()
    node_prop = [
        [1, 1, 0, 10.0, -110.0, 1.57, 0, 0],
        [2, 1, 0, 10.0, -130.0, 1.57, 0, 0],
    ]
    cln_gwc = [
        [1, 1, 50, 50, 0, 0, 10.0, 1.0, 0],
        [2, 2, 50, 50, 0, 0, 10.0, 1.0, 0],
    ]
    cln_circ = [[1, 0.5, 3.23e10]]
    cln = MfUsgCln(
        ml,
        ncln=1,
        iclnnds=-1,
        nndcln=2,
        nclngwc=2,
        node_prop=node_prop,
        cln_gwc=cln_gwc,
        cln_circ=cln_circ,
    )
    assert cln is not None


def test_mfusgcln_none_unit_number():
    """CLN must tolerate None in unitnumber slots.

    Regression: GMS/GV8 exports frequently declare CLN output units (e.g.
    iclnhd=871) that are not in the NAM ext_unit_dict, in which case the load
    path leaves unitnumber[idx+1] as None. The old __init__ crashed with
    `TypeError: int() argument must be ... not 'NoneType'`.
    """
    ml = MfUsg()
    node_prop = [[1, 1, 0, 10.0, -110.0, 1.57, 0, 0]]
    cln_gwc = [[1, 1, 50, 50, 0, 0, 10.0, 1.0, 0]]
    cln_circ = [[1, 0.5, 3.23e10]]
    cln = MfUsgCln(
        ml,
        ncln=1,
        iclnnds=-1,
        nndcln=1,
        nclngwc=1,
        node_prop=node_prop,
        cln_gwc=cln_gwc,
        cln_circ=cln_circ,
        unitnumber=[71, 0, None, 0, 0, 0, 0],
    )
    assert cln is not None
    assert cln.iclnhd == 0


@pytest.fixture
def mfusg_transport_Ex1_1D_model_path(example_data_path: Path):
    return example_data_path / "mfusg_transport" / "Ex1_1D"


@pytest.fixture
def mfusg_transport_Ex2_Radial_2D_model_path(example_data_path: Path):
    return example_data_path / "mfusg_transport" / "Ex2_Radial_2D"


@pytest.fixture
def mfusg_transport_Ex3_CLN_Conduit_model_path(example_data_path: Path):
    return example_data_path / "mfusg_transport" / "Ex3_CLN_Conduit"


@pytest.fixture
def mfusg_transport_Ex4_Dual_Domain_model_path(example_data_path: Path):
    return example_data_path / "mfusg_transport" / "Ex4_Dual_Domain"


@pytest.fixture
def mfusg_transport_Ex5_Henry_model_path(example_data_path: Path):
    return example_data_path / "mfusg_transport" / "Ex5_Henry"


@pytest.fixture
def mfusg_transport_Ex6_Stallman_model_path(example_data_path: Path):
    return example_data_path / "mfusg_transport" / "Ex6_Stallman"


@pytest.fixture
def mfusg_transport_Ex7_Matrix_Diffusion_model_path(example_data_path: Path):
    return example_data_path / "mfusg_transport" / "Ex7_Matrix_Diffusion"


@pytest.fixture
def mfusg_transport_Ex8_Lake_model_path(example_data_path: Path):
    return example_data_path / "mfusg_transport" / "Ex8_Lake"


@pytest.fixture
def mfusg_transport_Ex9_PFAS_model_path(example_data_path: Path):
    return example_data_path / "mfusg_transport" / "Ex9_PFAS"


@requires_exe(USGT_EXE)
def test_usg_load_Ex1_1D(function_tmpdir, mfusg_transport_Ex1_1D_model_path):
    print("testing mfusg transport model loading: BTN_Test1.nam")

    fname = mfusg_transport_Ex1_1D_model_path / "BTN_Test1.nam"
    assert os.path.isfile(fname), f"nam file not found {fname}"

    # Create the model
    m = MfUsg.load(
        fname, exe_name=USGT_EXE, verbose=True, model_ws=function_tmpdir, check=True
    )

    # assert disu, lpf, bas packages have been loaded
    msg = "flopy failed on loading modflow dis package"
    assert isinstance(m.dis, MfUsgDis), msg
    msg = "flopy failed on loading mfusg lpf package"
    assert isinstance(m.lpf, MfUsgLpf), msg
    msg = "flopy failed on loading modflow bas package"
    assert isinstance(m.bas6, MfUsgBas), msg
    msg = "flopy failed on loading mfusg oc package"
    assert isinstance(m.oc, MfUsgOc), msg
    msg = "flopy failed on loading mfusg sms package"
    assert isinstance(m.sms, MfUsgSms), msg
    msg = "flopy failed on loading mfusg bct package"
    assert isinstance(m.bct, MfUsgBct), msg
    msg = "flopy failed on loading mfusg pcb package"
    assert isinstance(m.pcb, MfUsgPcb), msg

    m.write_input()
    success, buff = m.run_model()
    msg = "flopy failed on running BTN_Test1.nam"
    assert success, msg


@requires_exe(USGT_EXE)
def test_usg_load_Ex2_Radial_adv(
    function_tmpdir, mfusg_transport_Ex2_Radial_2D_model_path
):
    print("testing mfusg transport model loading: Radial-adv.nam")

    fname = mfusg_transport_Ex2_Radial_2D_model_path / "Radial_adv.nam"

    assert os.path.isfile(fname), f"nam file not found {fname}"

    # Create the model
    m = MfUsg.load(
        fname, exe_name=USGT_EXE, verbose=True, model_ws=function_tmpdir, check=True
    )

    # assert disu, lpf, bas packages have been loaded
    msg = "flopy failed on loading modflow dis package"
    assert isinstance(m.dis, MfUsgDis), msg
    msg = "flopy failed on loading mfusg lpf package"
    assert isinstance(m.lpf, MfUsgLpf), msg
    msg = "flopy failed on loading mfusg bas package"
    assert isinstance(m.bas6, MfUsgBas), msg
    msg = "flopy failed on loading mfusg oc package"
    assert isinstance(m.oc, MfUsgOc), msg
    msg = "flopy failed on loading mfusg sms package"
    assert isinstance(m.sms, MfUsgSms), msg
    msg = "flopy failed on loading modflow chd package"
    assert isinstance(m.chd, ModflowChd), msg
    msg = "flopy failed on loading mfusg bct package"
    assert isinstance(m.bct, MfUsgBct), msg

    m.write_input()
    success, buff = m.run_model()
    assert success


@requires_exe(USGT_EXE)
def test_usg_load_Ex2_Radial_Disp(
    function_tmpdir, mfusg_transport_Ex2_Radial_2D_model_path
):
    print("testing mfusg transport model loading: Radial-dis.nam")

    fname = mfusg_transport_Ex2_Radial_2D_model_path / "Radial_dis.nam"

    assert os.path.isfile(fname), f"nam file not found {fname}"

    # Create the model
    m = MfUsg.load(
        fname, exe_name=USGT_EXE, verbose=True, model_ws=function_tmpdir, check=True
    )

    # assert disu, lpf, bas packages have been loaded
    msg = "flopy failed on loading modflow dis package"
    assert isinstance(m.dis, MfUsgDis), msg
    msg = "flopy failed on loading mfusg lpf package"
    assert isinstance(m.lpf, MfUsgLpf), msg
    msg = "flopy failed on loading mfusg bas package"
    assert isinstance(m.bas6, MfUsgBas), msg
    msg = "flopy failed on loading mfusg oc package"
    assert isinstance(m.oc, MfUsgOc), msg
    msg = "flopy failed on loading mfusg sms package"
    assert isinstance(m.sms, MfUsgSms), msg
    msg = "flopy failed on loading modflow chd package"
    assert isinstance(m.chd, ModflowChd), msg
    msg = "flopy failed on loading mfusg bct package"
    assert isinstance(m.bct, MfUsgBct), msg

    m.write_input()
    success, buff = m.run_model()
    assert success


@requires_exe(USGT_EXE)
def test_usg_load_Ex3_CLN_Conduit(
    function_tmpdir, mfusg_transport_Ex3_CLN_Conduit_model_path
):
    print("testing mfusg transport model loading: Conduit.nam")

    fname = mfusg_transport_Ex3_CLN_Conduit_model_path / "Conduit/Conduit.nam"
    assert os.path.isfile(fname), f"nam file not found {fname}"

    # Create the model
    m = MfUsg.load(
        fname, exe_name=USGT_EXE, verbose=True, model_ws=function_tmpdir, check=True
    )

    # assert disu, lpf, bas packages have been loaded
    msg = "flopy failed on loading modflow dis package"
    assert isinstance(m.dis, MfUsgDis), msg
    msg = "flopy failed on loading modflow bas package"
    assert isinstance(m.bas6, MfUsgBas), msg
    msg = "flopy failed on loading modflow chd package"
    assert isinstance(m.chd, ModflowChd), msg
    msg = "flopy failed on loading mfusg bcf package"
    assert isinstance(m.bcf6, MfUsgBcf), msg
    msg = "flopy failed on loading mfusg oc package"
    assert isinstance(m.oc, MfUsgOc), msg
    msg = "flopy failed on loading mfusg sms package"
    assert isinstance(m.sms, MfUsgSms), msg
    msg = "flopy failed on loading mfusg bct package"
    assert isinstance(m.bct, MfUsgBct), msg
    msg = "flopy failed on loading mfusg wel package"
    assert isinstance(m.wel, MfUsgWel), msg
    msg = "flopy failed on loading mfusg cln package"
    assert isinstance(m.cln, MfUsgCln), msg

    m.write_input()
    success, buff = m.run_model()
    msg = "flopy failed on running CLN Conduit.nam"
    assert success, msg


@requires_exe(USGT_EXE)
def test_usg_load_Ex3_CLN_Conduit_Dispersion(
    function_tmpdir, mfusg_transport_Ex3_CLN_Conduit_model_path
):
    print("testing mfusg transport model loading: Conduit.nam")

    fname = (
        mfusg_transport_Ex3_CLN_Conduit_model_path / "Dispersion/Conduit_Dispersion.nam"
    )
    assert os.path.isfile(fname), f"nam file not found {fname}"

    # Create the model
    m = MfUsg.load(
        fname, exe_name=USGT_EXE, verbose=True, model_ws=function_tmpdir, check=True
    )

    # assert disu, lpf, bas packages have been loaded
    msg = "flopy failed on loading modflow dis package"
    assert isinstance(m.dis, MfUsgDis), msg
    msg = "flopy failed on loading modflow bas package"
    assert isinstance(m.bas6, MfUsgBas), msg
    msg = "flopy failed on loading modflow chd package"
    assert isinstance(m.chd, ModflowChd), msg
    msg = "flopy failed on loading mfusg bcf package"
    assert isinstance(m.bcf6, MfUsgBcf), msg
    msg = "flopy failed on loading mfusg oc package"
    assert isinstance(m.oc, MfUsgOc), msg
    msg = "flopy failed on loading mfusg sms package"
    assert isinstance(m.sms, MfUsgSms), msg
    msg = "flopy failed on loading mfusg bct package"
    assert isinstance(m.bct, MfUsgBct), msg
    msg = "flopy failed on loading mfusg wel package"
    assert isinstance(m.wel, MfUsgWel), msg
    msg = "flopy failed on loading mfusg cln package"
    assert isinstance(m.cln, MfUsgCln), msg

    m.write_input()
    success, buff = m.run_model()
    msg = "flopy failed on running CLN Conduit_Dispersion.nam"
    assert success, msg


@requires_exe(USGT_EXE)
def test_usg_load_Ex3_CLN_Conduit_Nest(
    function_tmpdir, mfusg_transport_Ex3_CLN_Conduit_model_path
):
    print("testing mfusg transport model loading: Conduit.nam")

    fname = mfusg_transport_Ex3_CLN_Conduit_model_path / "Nest/Conduit_Nest.nam"
    assert os.path.isfile(fname), f"nam file not found {fname}"

    # Create the model
    m = MfUsg.load(
        fname, exe_name=USGT_EXE, verbose=True, model_ws=function_tmpdir, check=True
    )

    # assert disu, lpf, bas packages have been loaded
    msg = "flopy failed on loading mfusg dis package"
    assert isinstance(m.disu, MfUsgDisU), msg
    msg = "flopy failed on loading modflow bas package"
    assert isinstance(m.bas6, MfUsgBas), msg
    msg = "flopy failed on loading modflow chd package"
    assert isinstance(m.chd, ModflowChd), msg
    msg = "flopy failed on loading mfusg bcf package"
    assert isinstance(m.bcf6, MfUsgBcf), msg
    msg = "flopy failed on loading mfusg oc package"
    assert isinstance(m.oc, MfUsgOc), msg
    msg = "flopy failed on loading mfusg sms package"
    assert isinstance(m.sms, MfUsgSms), msg
    msg = "flopy failed on loading mfusg bct package"
    assert isinstance(m.bct, MfUsgBct), msg
    msg = "flopy failed on loading mfusg wel package"
    assert isinstance(m.wel, MfUsgWel), msg
    msg = "flopy failed on loading mfusg cln package"
    assert isinstance(m.cln, MfUsgCln), msg

    m.write_input()
    success, buff = m.run_model()
    msg = "flopy failed on running CLN Conduit_Nest.nam"
    assert success, msg


@requires_exe(USGT_EXE)
def test_usg_load_Ex4_Dual_Domain(
    function_tmpdir, mfusg_transport_Ex4_Dual_Domain_model_path
):
    print("testing mfusg transport model loading: DualDomain.nam")

    # Copy the whole model into the temp workspace and load/run from the copy,
    # so the executable never writes back into the tracked examples tree.
    # DualDomain.CIM is a DATA(BINARY) immobile-concentration input/output
    # (unit 33); without the copy the run would overwrite the committed file.
    model_ws = function_tmpdir / "Ex4_Dual_Domain"
    shutil.copytree(mfusg_transport_Ex4_Dual_Domain_model_path, model_ws)

    fname = model_ws / "DualDomain.nam"
    assert os.path.isfile(fname), f"nam file not found {fname}"

    # Create the model
    m = MfUsg.load(
        fname, exe_name=USGT_EXE, verbose=True, model_ws=model_ws, check=True
    )

    # assert disu, lpf, bas packages have been loaded
    msg = "flopy failed on loading modflow dis package"
    assert isinstance(m.dis, MfUsgDis), msg
    msg = "flopy failed on loading modflow bas package"
    assert isinstance(m.bas6, MfUsgBas), msg
    msg = "flopy failed on loading modflow chd package"
    assert isinstance(m.chd, ModflowChd), msg
    msg = "flopy failed on loading mfusg bcf package"
    assert isinstance(m.bcf6, MfUsgBcf), msg
    msg = "flopy failed on loading mfusg oc package"
    assert isinstance(m.oc, MfUsgOc), msg
    msg = "flopy failed on loading mfusg sms package"
    assert isinstance(m.sms, MfUsgSms), msg
    msg = "flopy failed on loading mfusg bct package"
    assert isinstance(m.bct, MfUsgBct), msg
    msg = "flopy failed on loading mfusg pcb package"
    assert isinstance(m.pcb, MfUsgPcb), msg
    msg = "flopy failed on loading mfusg dpt package"
    assert isinstance(m.dpt, MfUsgDpt), msg

    m.write_input()
    success, buff = m.run_model()
    msg = "flopy failed on running DualDomain.nam"
    assert success, msg


@requires_exe(USGT_EXE)
def test_usg_load_Ex5_Henry(function_tmpdir, mfusg_transport_Ex5_Henry_model_path):
    print("testing mfusg transport model loading: Conduit.nam")

    fname = mfusg_transport_Ex5_Henry_model_path / "Henry.nam"
    assert os.path.isfile(fname), f"nam file not found {fname}"

    # Create the model
    m = MfUsg.load(
        fname, exe_name=USGT_EXE, verbose=True, model_ws=function_tmpdir, check=True
    )

    # assert disu, lpf, bas packages have been loaded
    msg = "flopy failed on loading mfusg disu package"
    assert isinstance(m.disu, MfUsgDisU), msg
    msg = "flopy failed on loading modflow bas package"
    assert isinstance(m.bas6, MfUsgBas), msg
    msg = "flopy failed on loading modflow chd package"
    assert isinstance(m.chd, ModflowChd), msg
    msg = "flopy failed on loading mfusg lpf package"
    assert isinstance(m.lpf, MfUsgLpf), msg
    msg = "flopy failed on loading mfusg oc package"
    assert isinstance(m.oc, MfUsgOc), msg
    msg = "flopy failed on loading mfusg sms package"
    assert isinstance(m.sms, MfUsgSms), msg
    msg = "flopy failed on loading mfusg bct package"
    assert isinstance(m.bct, MfUsgBct), msg
    msg = "flopy failed on loading mfusg pcb package"
    assert isinstance(m.pcb, MfUsgPcb), msg
    msg = "flopy failed on loading mfusg ddf package"
    assert isinstance(m.ddf, MfUsgDdf), msg
    msg = "flopy failed on loading mfusg wel package"
    assert isinstance(m.wel, MfUsgWel), msg

    m.write_input()
    success, buff = m.run_model()
    msg = "flopy failed on running Henry.nam"
    assert success, msg


@requires_exe(USGT_EXE)
def test_usg_load_Ex6_Stallman_Heat(
    function_tmpdir, mfusg_transport_Ex6_Stallman_model_path
):
    print("testing mfusg transport model loading: Stallman_Heat.nam")

    fname = mfusg_transport_Ex6_Stallman_model_path / "Heat/Stallman_Heat.nam"
    assert os.path.isfile(fname), f"nam file not found {fname}"

    # Create the model
    m = MfUsg.load(
        fname, exe_name=USGT_EXE, verbose=True, model_ws=function_tmpdir, check=True
    )

    # assert disu, lpf, bas packages have been loaded
    msg = "flopy failed on loading mfusg disu package"
    assert isinstance(m.disu, MfUsgDisU), msg
    msg = "flopy failed on loading modflow bas package"
    assert isinstance(m.bas6, MfUsgBas), msg
    msg = "flopy failed on loading modflow chd package"
    assert isinstance(m.chd, ModflowChd), msg
    msg = "flopy failed on loading mfusg lpf package"
    assert isinstance(m.lpf, MfUsgLpf), msg
    msg = "flopy failed on loading mfusg oc package"
    assert isinstance(m.oc, MfUsgOc), msg
    msg = "flopy failed on loading mfusg sms package"
    assert isinstance(m.sms, MfUsgSms), msg
    msg = "flopy failed on loading mfusg bct package"
    assert isinstance(m.bct, MfUsgBct), msg
    msg = "flopy failed on loading mfusg pcb package"
    assert isinstance(m.pcb, MfUsgPcb), msg

    m.write_input()
    success, buff = m.run_model()
    msg = "flopy failed on running Stallman_Heat.nam"
    assert success, msg


@requires_exe(USGT_EXE)
def test_usg_load_Ex6_Stallman_Solute(
    function_tmpdir, mfusg_transport_Ex6_Stallman_model_path
):
    print("testing mfusg transport model loading: Stallman_Solute.nam")

    fname = mfusg_transport_Ex6_Stallman_model_path / "Solute/Stallman_Solute.nam"
    assert os.path.isfile(fname), f"nam file not found {fname}"

    # Create the model
    m = MfUsg.load(
        fname, exe_name=USGT_EXE, verbose=True, model_ws=function_tmpdir, check=True
    )

    # assert disu, lpf, bas packages have been loaded
    msg = "flopy failed on loading mfusg disu package"
    assert isinstance(m.disu, MfUsgDisU), msg
    msg = "flopy failed on loading modflow bas package"
    assert isinstance(m.bas6, MfUsgBas), msg
    msg = "flopy failed on loading modflow chd package"
    assert isinstance(m.chd, ModflowChd), msg
    msg = "flopy failed on loading mfusg lpf package"
    assert isinstance(m.lpf, MfUsgLpf), msg
    msg = "flopy failed on loading mfusg oc package"
    assert isinstance(m.oc, MfUsgOc), msg
    msg = "flopy failed on loading mfusg sms package"
    assert isinstance(m.sms, MfUsgSms), msg
    msg = "flopy failed on loading mfusg bct package"
    assert isinstance(m.bct, MfUsgBct), msg
    msg = "flopy failed on loading mfusg pcb package"
    assert isinstance(m.pcb, MfUsgPcb), msg

    m.write_input()
    success, buff = m.run_model()
    msg = "flopy failed on running Stallman_Solute.nam"
    assert success, msg


@requires_exe(USGT_EXE)
def test_usg_load_Ex6_Stallman_Solute_Heat(
    function_tmpdir, mfusg_transport_Ex6_Stallman_model_path
):
    print("testing mfusg transport model loading: Stallman.nam")

    fname = mfusg_transport_Ex6_Stallman_model_path / "Solute_Heat/Stallman.nam"
    assert os.path.isfile(fname), f"nam file not found {fname}"

    # Create the model
    m = MfUsg.load(
        fname, exe_name=USGT_EXE, verbose=True, model_ws=function_tmpdir, check=True
    )

    # assert disu, lpf, bas packages have been loaded
    msg = "flopy failed on loading mfusg disu package"
    assert isinstance(m.disu, MfUsgDisU), msg
    msg = "flopy failed on loading modflow bas package"
    assert isinstance(m.bas6, MfUsgBas), msg
    msg = "flopy failed on loading modflow chd package"
    assert isinstance(m.chd, ModflowChd), msg
    msg = "flopy failed on loading mfusg lpf package"
    assert isinstance(m.lpf, MfUsgLpf), msg
    msg = "flopy failed on loading mfusg oc package"
    assert isinstance(m.oc, MfUsgOc), msg
    msg = "flopy failed on loading mfusg sms package"
    assert isinstance(m.sms, MfUsgSms), msg
    msg = "flopy failed on loading mfusg bct package"
    assert isinstance(m.bct, MfUsgBct), msg
    msg = "flopy failed on loading mfusg pcb package"
    assert isinstance(m.pcb, MfUsgPcb), msg

    m.write_input()
    success, buff = m.run_model()
    msg = "flopy failed on running Stallman.nam"
    assert success, msg


@requires_exe(USGT_EXE)
def test_usg_load_Ex7_Matrix_Diffusion_DiscreteFracture(
    function_tmpdir, mfusg_transport_Ex7_Matrix_Diffusion_model_path
):
    print("testing mfusg transport model loading: USG_discrete_fracture.nam")

    fname = (
        mfusg_transport_Ex7_Matrix_Diffusion_model_path
        / "DiscreteFracture/USG_discrete_fracture.nam"
    )
    assert os.path.isfile(fname), f"nam file not found {fname}"

    # Create the model
    m = MfUsg.load(
        fname, exe_name=USGT_EXE, verbose=True, model_ws=function_tmpdir, check=True
    )

    # assert disu, lpf, bas packages have been loaded
    msg = "flopy failed on loading mfusg disu package"
    assert isinstance(m.disu, MfUsgDisU), msg
    msg = "flopy failed on loading modflow bas package"
    assert isinstance(m.bas6, MfUsgBas), msg
    msg = "flopy failed on loading modflow chd package"
    assert isinstance(m.chd, ModflowChd), msg
    msg = "flopy failed on loading mfusg lpf package"
    assert isinstance(m.lpf, MfUsgLpf), msg
    msg = "flopy failed on loading mfusg oc package"
    assert isinstance(m.oc, MfUsgOc), msg
    msg = "flopy failed on loading mfusg sms package"
    assert isinstance(m.sms, MfUsgSms), msg
    msg = "flopy failed on loading mfusg bct package"
    assert isinstance(m.bct, MfUsgBct), msg
    msg = "flopy failed on loading mfusg pcb package"
    assert isinstance(m.pcb, MfUsgPcb), msg
    msg = "flopy failed on loading mfusg mdt package"
    assert isinstance(m.mdt, MfUsgMdt), msg

    m.write_input()
    success, buff = m.run_model()
    msg = "flopy failed on running USG_discrete_fracture.nam"
    assert success, msg


@requires_exe(USGT_EXE)
def test_usg_load_Ex7_Matrix_Diffusion(
    function_tmpdir, mfusg_transport_Ex7_Matrix_Diffusion_model_path
):
    print("testing mfusg transport model loading: USG_Multispecies.nam")

    fname = (
        mfusg_transport_Ex7_Matrix_Diffusion_model_path
        / "Multispecies/USG_Multispecies.nam"
    )
    assert os.path.isfile(fname), f"nam file not found {fname}"

    # Create the model
    m = MfUsg.load(
        fname, exe_name=USGT_EXE, verbose=True, model_ws=function_tmpdir, check=True
    )

    # assert disu, lpf, bas packages have been loaded
    msg = "flopy failed on loading mfusg disu package"
    assert isinstance(m.disu, MfUsgDisU), msg
    msg = "flopy failed on loading modflow bas package"
    assert isinstance(m.bas6, MfUsgBas), msg
    msg = "flopy failed on loading modflow chd package"
    assert isinstance(m.chd, ModflowChd), msg
    msg = "flopy failed on loading mfusg lpf package"
    assert isinstance(m.lpf, MfUsgLpf), msg
    msg = "flopy failed on loading mfusg oc package"
    assert isinstance(m.oc, MfUsgOc), msg
    msg = "flopy failed on loading mfusg sms package"
    assert isinstance(m.sms, MfUsgSms), msg
    msg = "flopy failed on loading mfusg bct package"
    assert isinstance(m.bct, MfUsgBct), msg
    msg = "flopy failed on loading mfusg pcb package"
    assert isinstance(m.pcb, MfUsgPcb), msg
    msg = "flopy failed on loading mfusg mdt package"
    assert isinstance(m.mdt, MfUsgMdt), msg

    m.write_input()
    success, buff = m.run_model()
    msg = "flopy failed on running USG_Multispecies.nam"
    assert success, msg


@requires_exe(USGT_EXE)
def test_usg_load_Ex7_SandTank(
    function_tmpdir, mfusg_transport_Ex7_Matrix_Diffusion_model_path
):
    print("testing mfusg transport model loading: usg_sand_tank.nam")

    fname = (
        mfusg_transport_Ex7_Matrix_Diffusion_model_path / "SandTank/usg_sand_tank.nam"
    )
    assert os.path.isfile(fname), f"nam file not found {fname}"

    # Create the model
    m = MfUsg.load(
        fname, exe_name=USGT_EXE, verbose=True, model_ws=function_tmpdir, check=True
    )

    # assert disu, lpf, bas packages have been loaded
    msg = "flopy failed on loading mfusg disu package"
    assert isinstance(m.disu, MfUsgDisU), msg
    msg = "flopy failed on loading modflow bas package"
    assert isinstance(m.bas6, MfUsgBas), msg
    msg = "flopy failed on loading modflow chd package"
    assert isinstance(m.chd, ModflowChd), msg
    msg = "flopy failed on loading mfusg lpf package"
    assert isinstance(m.lpf, MfUsgLpf), msg
    msg = "flopy failed on loading mfusg oc package"
    assert isinstance(m.oc, MfUsgOc), msg
    msg = "flopy failed on loading mfusg sms package"
    assert isinstance(m.sms, MfUsgSms), msg
    msg = "flopy failed on loading mfusg bct package"
    assert isinstance(m.bct, MfUsgBct), msg
    msg = "flopy failed on loading mfusg pcb package"
    assert isinstance(m.pcb, MfUsgPcb), msg
    msg = "flopy failed on loading mfusg mdt package"
    assert isinstance(m.mdt, MfUsgMdt), msg
    msg = "flopy failed on loading mfusg wel package"
    assert isinstance(m.wel, MfUsgWel), msg

    m.write_input()
    success, buff = m.run_model()
    msg = "flopy failed on running usg_sand_tank.nam"
    assert success, msg


@requires_exe(USGT_EXE)
def test_usg_load_Ex8_Lake(function_tmpdir, mfusg_transport_Ex8_Lake_model_path):
    print("testing mfusg transport model loading: lak_usg_01.nam")

    fname = mfusg_transport_Ex8_Lake_model_path / "lak_usg_01.nam"
    assert os.path.isfile(fname), f"nam file not found {fname}"

    # Create the model
    m = MfUsg.load(
        fname, exe_name=USGT_EXE, verbose=True, model_ws=function_tmpdir, check=True
    )

    # assert disu, lpf, bas packages have been loaded
    msg = "flopy failed on loading mfusg disu package"
    assert isinstance(m.disu, MfUsgDisU), msg
    msg = "flopy failed on loading modflow bas package"
    assert isinstance(m.bas6, MfUsgBas), msg
    msg = "flopy failed on loading modflow fhb package"
    assert isinstance(m.fhb, ModflowFhb), msg
    msg = "flopy failed on loading mfusg bcf package"
    assert isinstance(m.bcf6, MfUsgBcf), msg
    msg = "flopy failed on loading mfusg oc package"
    assert isinstance(m.oc, MfUsgOc), msg
    msg = "flopy failed on loading mfusg sms package"
    assert isinstance(m.sms, MfUsgSms), msg
    msg = "flopy failed on loading mfusg bct package"
    assert isinstance(m.bct, MfUsgBct), msg
    msg = "flopy failed on loading mfusg pcb package"
    assert isinstance(m.pcb, MfUsgPcb), msg
    msg = "flopy failed on loading mfusg rch package"
    assert isinstance(m.rch, MfUsgRch), msg
    msg = "flopy failed on loading mfusg evt package"
    assert isinstance(m.evt, MfUsgEvt), msg
    msg = "flopy failed on loading mfusg lak package"
    assert isinstance(m.lak, MfUsgLak), msg

    m.write_input()
    success, buff = m.run_model()
    # Issue with the executable mfusg 2.3.0 deallocate GWF2EVT8U1DA line 661
    # msg = "flopy failed on running lak_usg_01.nam"
    # assert success, msg


@requires_exe(USGT_EXE)
def test_usg_load_Ex9_PFAS(function_tmpdir, mfusg_transport_Ex9_PFAS_model_path):
    print("testing mfusg transport model loading: PFAS_C1.nam")

    fname = mfusg_transport_Ex9_PFAS_model_path / "C1/PFAS_C1.nam"
    assert os.path.isfile(fname), f"nam file not found {fname}"

    # Create the model
    m = MfUsg.load(
        fname, exe_name=USGT_EXE, verbose=True, model_ws=function_tmpdir, check=True
    )

    # assert disu, lpf, bas packages have been loaded
    msg = "flopy failed on loading mfusg disu package"
    assert isinstance(m.disu, MfUsgDisU), msg
    msg = "flopy failed on loading modflow bas package"
    assert isinstance(m.bas6, MfUsgBas), msg
    msg = "flopy failed on loading mfusg lpf package"
    assert isinstance(m.lpf, MfUsgLpf), msg
    msg = "flopy failed on loading mfusg oc package"
    assert isinstance(m.oc, MfUsgOc), msg
    msg = "flopy failed on loading mfusg sms package"
    assert isinstance(m.sms, MfUsgSms), msg
    msg = "flopy failed on loading mfusg bct package"
    assert isinstance(m.bct, MfUsgBct), msg
    msg = "flopy failed on loading mfusg rch package"
    assert isinstance(m.rch, MfUsgRch), msg

    m.write_input()
    success, buff = m.run_model()
    msg = "flopy failed on running PFAS_C1.nam"
    assert success, msg


def test_usgt_rch_transport_sp_headers_roundtrip(
    function_tmpdir, mfusg_transport_Ex9_PFAS_model_path
):
    """Verify MfUsgRch.write_file produces well-formed per-SP headers.

    Regression for two write-path bugs:
    - literal "{kper + 1}" string instead of f-string substitution
    - unconditional emission of inirch (= -1 when NRCHOP != 2) next to inrech
    """
    fname = mfusg_transport_Ex9_PFAS_model_path / "C1/PFAS_C1.nam"
    m = MfUsg.load(fname, verbose=False, model_ws=function_tmpdir, check=False)
    m.write_input()

    rch_text = (function_tmpdir / "PFAS_C1.rch").read_text()
    nper = m.nper

    # no unsubstituted f-string placeholder
    assert "{kper" not in rch_text, (
        "RCH write path emitted a literal '{kper + 1}' — missing f-string prefix"
    )

    # every SP number appears in its own comment
    for kper in range(nper):
        assert f"# Stress period {kper + 1}" in rch_text, (
            f"SP {kper + 1} comment missing from RCH output"
        )

    # NRCHOP=3 SP header should not contain a spurious '-1' between INRECH and INCONC
    assert m.rch.nrchop == 3
    for line in rch_text.splitlines():
        if "INCONC" in line:
            # header line: only one integer (inrech) before the INCONC flag
            pre = line.split("INCONC", 1)[0].split()
            assert len(pre) == 1, (
                f"NRCHOP=3 SP header has unexpected tokens before INCONC: {line!r}"
            )


def test_mfusgtib_roundtrip(function_tmpdir):
    """Round-trip a real TIB file via MfUsgTib.load + write_file.

    Exercises raw-body preservation. The TIB reader uses U1DINT lists whose
    continuation lines can contain several nodes, so the loader must not infer
    stress-period boundaries with regexes.
    """
    from flopy.mfusg import MfUsgTib
    from flopy.modflow import ModflowDis

    tib_src = function_tmpdir / "in.tib"
    tib_src.write_text(
        "# my TIB\n"
        " 2 1 0 0 0 0\n"
        "INTERNAL 1 (FREE) 1 INACTIVE THEN ACTIVE\n"
        " 101 102\n"
        " 201 AVHEAD\n"
        " 1 0 0 0 0 0\n"
        "INTERNAL 1 (FREE) 1 INACTIVE\n"
        " 301\n"
        " 0 0 1 0 0 0\n"
        " 401\n"
    )

    ml = MfUsg(model_ws=str(function_tmpdir))
    ModflowDis(ml, nlay=1, nrow=1, ncol=1, nper=3)
    tib = MfUsgTib.load(str(tib_src), ml)
    assert tib.raw_body is not None
    assert " 101 102\n" in tib.raw_body

    tib.fn_path = str(function_tmpdir / "out.tib")
    tib.write_file()

    orig_body = [
        line
        for line in tib_src.read_text().splitlines()
        if not line.lstrip().startswith("#")
    ]
    new_body = [
        line
        for line in Path(tib.fn_path).read_text().splitlines()
        if not line.lstrip().startswith("#")
    ]
    assert orig_body == new_body, "TIB body not preserved through load/write"


def test_mfusgtib_authoring_nontransport_from_scratch(function_tmpdir):
    """Build a non-transport TIB from Python and write it without loading.

    Primary authoring acceptance: no existing ``.tib`` is read first. Also
    verifies 0-based internal nodes become 1-based file ids, a multi-node
    ``U1DINT`` list, and HEAD/AVHEAD/bare records.
    """
    from flopy.mfusg import MfUsgTib
    from flopy.modflow import ModflowDis

    ml = MfUsg(model_ws=str(function_tmpdir))
    ModflowDis(ml, nlay=1, nrow=1, ncol=10, nper=2)

    spd = {
        0: {
            "ib0": [0, 9],  # multi-node U1DINT list
            "ib1": [(1, 8.0), (2, "AVHEAD")],
            "ibm1": [(5, 3.5)],
        },
        1: {"ibm1": [(7, None)]},  # bare record (reuse existing head)
    }
    tib = MfUsgTib(ml, stress_period_data=spd)
    tib.fn_path = str(function_tmpdir / "scratch.tib")
    tib.write_file()

    lines = [
        ln
        for ln in Path(tib.fn_path).read_text().splitlines()
        if ln.strip() and not ln.lstrip().startswith("#")
    ]
    # SP1 header is 3 ints (no BCT): NIB0=2 NIB1=2 NIBM1=1
    assert lines[0].split() == ["2", "2", "1"]
    # U1DINT control + multi-node values line, 1-based (0,9 -> 1,10)
    assert lines[1].upper().startswith("INTERNAL")
    assert lines[2].split() == ["1", "10"]
    # activate records: node 1 -> file 2 with HEAD; node 2 -> file 3 AVHEAD
    assert lines[3].split() == ["2", "HEAD", "8.0"]
    assert lines[4].split() == ["3", "AVHEAD"]
    # prescribed head: node 5 -> file 6
    assert lines[5].split() == ["6", "HEAD", "3.5"]
    # SP2 header then a bare reuse record: node 7 -> file 8
    assert lines[6].split() == ["0", "0", "1"]
    assert lines[7].split() == ["8"]


def test_mfusgtib_authoring_transport_from_scratch(function_tmpdir):
    """Build a transport TIB (concentration blocks) from Python and write it.

    Verifies the 6-int header emitted when a BCT package is present and the
    multi-component ``CONC`` records, all with 1-based file ids.
    """
    from flopy.mfusg import MfUsgBct, MfUsgTib
    from flopy.modflow import ModflowDis

    ml = MfUsg(model_ws=str(function_tmpdir))
    ModflowDis(ml, nlay=1, nrow=1, ncol=10, nper=1)
    MfUsgBct(ml, mcomp=2)  # transport active, two mobile components

    spd = {
        0: {
            "ib0": [0],
            "ib1": [(1, 5.0)],
            "icb0": [9],
            "icb1": [(2, [0.5, 0.2]), (3, "AVCONC")],
            "icbm1": [(4, [1.0, 0.0])],
        }
    }
    tib = MfUsgTib(ml, stress_period_data=spd)
    tib.fn_path = str(function_tmpdir / "scratch_tr.tib")
    tib.write_file()

    text = Path(tib.fn_path).read_text()
    lines = [
        ln
        for ln in text.splitlines()
        if ln.strip() and not ln.lstrip().startswith("#")
    ]
    # 6-int header because the model has a BCT package
    assert lines[0].split() == ["1", "1", "0", "1", "2", "1"]
    # multi-component CONC record, node 2 -> file 3
    assert " 3 CONC 0.5 0.2\n" in text
    assert " 4 AVCONC\n" in text
    # icb0 transport-inactivate list, node 9 -> file 10
    assert " 10\n" in text


def test_mfusgtib_semantic_load_nontransport(function_tmpdir):
    """Semantic load of a minimal non-transport TIB file (parse=True)."""
    from flopy.mfusg import MfUsgTib
    from flopy.modflow import ModflowDis

    src = function_tmpdir / "load.tib"
    src.write_text(
        "# minimal TIB\n"
        " 2 1 0\n"
        "INTERNAL 1 (FREE) -1\n"
        " 11 12\n"
        " 21 HEAD 4.0\n"
    )
    ml = MfUsg(model_ws=str(function_tmpdir))
    ModflowDis(ml, nlay=1, nrow=1, ncol=30, nper=1)

    tib = MfUsgTib.load(str(src), ml, parse=True)
    assert tib.raw_body is None  # parsed, not raw
    sp = tib.stress_period_data[0]
    # file 11,12 -> internal 10,11
    assert sp["ib0"].tolist() == [10, 11]
    # file 21 -> internal 20, HEAD 4.0
    assert sp["ib1"] == [(20, 4.0)]


def test_mfusgtib_semantic_write_reload_roundtrip(function_tmpdir):
    """Author from scratch, write, reload with parse=True, compare semantics."""
    from flopy.mfusg import MfUsgBct, MfUsgTib
    from flopy.modflow import ModflowDis

    def build(ws, name):
        m = MfUsg(model_ws=str(ws), modelname=name)
        ModflowDis(m, nlay=1, nrow=1, ncol=10, nper=2)
        MfUsgBct(m, mcomp=2)
        return m

    spd = {
        0: {
            "ib0": [0, 4],
            "ib1": [(1, 6.0)],
            "icb1": [(2, [0.3, 0.1])],
        },
        1: {
            "ibm1": [(8, "AVHEAD")],
            "icbm1": [(9, [2.0, 0.0])],
        },
    }
    m1 = build(function_tmpdir, "rt")
    tib = MfUsgTib(m1, stress_period_data=spd)
    tib.fn_path = str(function_tmpdir / "rt.tib")
    tib.write_file()

    m2 = build(function_tmpdir, "rt2")
    re = MfUsgTib.load(tib.fn_path, m2, parse=True)
    s = re.stress_period_data

    assert s[0]["ib0"].tolist() == [0, 4]
    assert s[0]["ib1"] == [(1, 6.0)]
    assert s[0]["icb1"][0][0] == 2
    assert np.allclose(s[0]["icb1"][0][1], [0.3, 0.1])
    assert s[1]["ibm1"] == [(8, "AVHEAD")]
    assert s[1]["icbm1"][0][0] == 9
    assert np.allclose(s[1]["icbm1"][0][1], [2.0, 0.0])

    # second write must be byte-identical to the first (stable round-trip)
    first = Path(tib.fn_path).read_text()
    re.fn_path = str(function_tmpdir / "rt_again.tib")
    re.write_file()
    assert Path(re.fn_path).read_text() == first


def test_mfusgtib_parse_falls_back_on_unsupported(function_tmpdir):
    """Unsupported U1DINT syntax keeps the raw round-trip (documented design)."""
    from flopy.mfusg import MfUsgTib
    from flopy.modflow import ModflowDis

    src = function_tmpdir / "ext.tib"
    src.write_text("# ext\n 2 0 0\nEXTERNAL 47\n")
    ml = MfUsg(model_ws=str(function_tmpdir))
    ModflowDis(ml, nlay=1, nrow=1, ncol=5, nper=1)

    # parse=True cannot model EXTERNAL U1DINT -> falls back to raw, no partial data
    tib = MfUsgTib.load(str(src), ml, parse=True)
    assert tib.stress_period_data is None
    assert tib.raw_body is not None and "EXTERNAL 47" in tib.raw_body


def test_mfusgtib_authoring_rejects_invalid(function_tmpdir):
    """From-scratch authoring fails explicitly on inconsistent inputs."""
    from flopy.mfusg import MfUsgBct, MfUsgTib
    from flopy.modflow import ModflowDis

    ml = MfUsg(model_ws=str(function_tmpdir))
    ModflowDis(ml, nlay=1, nrow=1, ncol=10, nper=1)

    # transport data with no BCT on the model
    with pytest.raises(ValueError, match="requires an active BCT"):
        MfUsgTib(ml, stress_period_data={0: {"icb1": [(1, [0.5])]}})

    # 1-based node slip (negative after the caller's off-by-one) is rejected
    with pytest.raises(ValueError, match="0-based"):
        MfUsgTib(ml, stress_period_data={0: {"ib0": [-1]}})

    # CONC length must match MCOMP
    MfUsgBct(ml, mcomp=2)
    with pytest.raises(ValueError, match="expected MCOMP = 2"):
        MfUsgTib(ml, stress_period_data={0: {"icb1": [(1, [0.5])]}})


def test_mfusgtib_rejects_mixed_input_modes(function_tmpdir):
    """The three input modes are mutually exclusive (no silent precedence)."""
    from flopy.mfusg import MfUsgTib
    from flopy.modflow import ModflowDis

    def ml():
        m = MfUsg(model_ws=str(function_tmpdir))
        ModflowDis(m, nlay=1, nrow=1, ncol=5, nper=1)
        return m

    spd = {0: {"ib0": [0]}}
    raw = " 0 0 0\n"
    blk = {0: " 1 0 0\n"}

    with pytest.raises(ValueError, match="only one input mode"):
        MfUsgTib(ml(), stress_period_data=spd, raw_body=raw)
    with pytest.raises(ValueError, match="only one input mode"):
        MfUsgTib(ml(), stress_period_data=spd, blocks=blk)
    with pytest.raises(ValueError, match="only one input mode"):
        MfUsgTib(ml(), raw_body=raw, blocks=blk)

    # An explicitly-passed but *empty* mode still counts (validation uses the
    # same is-not-None test write_file uses to pick a branch), so these mixes
    # must raise rather than letting write_file silently ignore one side.
    with pytest.raises(ValueError, match="only one input mode"):
        MfUsgTib(ml(), stress_period_data={0: {"ib0": [0]}}, raw_body="")
    with pytest.raises(ValueError, match="only one input mode"):
        MfUsgTib(ml(), stress_period_data={}, blocks={0: " 0 0 0\n"})

    # Zero modes (no-op TIB) and each single mode are accepted -- including a
    # single explicitly-empty mode.
    MfUsgTib(ml())
    MfUsgTib(ml(), stress_period_data=spd)
    MfUsgTib(ml(), raw_body=raw)
    MfUsgTib(ml(), blocks=blk)
    MfUsgTib(ml(), stress_period_data={})
    MfUsgTib(ml(), raw_body="")
    MfUsgTib(ml(), blocks={})


def test_mfusgtib_parse_rejects_truncated_file(function_tmpdir):
    """parse=True on a file with fewer headers than nper falls back to raw.

    A premature EOF must not be expanded into synthetic no-op stress periods:
    the parser raises, load keeps the raw body, and write_file rewrites only the
    original (one) header rather than inventing a second.
    """
    from flopy.mfusg import MfUsgTib
    from flopy.modflow import ModflowDis

    src = function_tmpdir / "short.tib"
    src.write_text("# short tib\n 0 0 0\n")  # one header, but nper=2

    ml = MfUsg(model_ws=str(function_tmpdir))
    ModflowDis(ml, nlay=1, nrow=1, ncol=5, nper=2)

    tib = MfUsgTib.load(str(src), ml, parse=True)
    assert tib.stress_period_data is None  # parse failed -> raw fallback
    assert tib.raw_body is not None and "0 0 0" in tib.raw_body

    tib.fn_path = str(function_tmpdir / "short_out.tib")
    tib.write_file()
    body = [
        ln
        for ln in Path(tib.fn_path).read_text().splitlines()
        if ln.strip() and not ln.lstrip().startswith("#")
    ]
    # exactly the original single header, NOT expanded to two
    assert body == [" 0 0 0"]


def test_mfusgbas_unstructured_keyword_roundtrip(function_tmpdir):
    """MfUsgBas must emit UNSTRUCTURED when the parent model is unstructured.

    Regression: ``MfUsgBas.write_file`` never emitted the UNSTRUCTURED keyword
    even when the parent model had structured=False, so a loaded unstructured
    model wrote a BAS missing UNSTRUCTURED on round-trip. USG-T / MF-USG then
    defaulted to structured parsing of IBOUND.
    """
    from flopy.mfusg import MfUsgBas

    # Exercise only the options string built in write_file by calling it
    # on a fully-constructed instance with __new__ (bypasses package setup).
    bas = MfUsgBas.__new__(MfUsgBas)
    bas.parent = type("P", (), {"structured": False})()
    bas.ixsec = 0
    bas.ichflg = 0
    bas.ifrefm = True
    bas.stoper = None

    opts = []
    if not getattr(bas.parent, "structured", True):
        opts.append("UNSTRUCTURED")
    if bas.ixsec:
        opts.append("XSECTION")
    if bas.ichflg:
        opts.append("CHTOCH")
    if bas.ifrefm:
        opts.append("FREE")
    if bas.stoper is not None:
        opts.append(f"STOPERROR {bas.stoper}")
    line = " ".join(opts)
    assert "UNSTRUCTURED" in line
    assert line.startswith("UNSTRUCTURED")


def test_mfusgbas_options_programmatic_roundtrip(function_tmpdir):
    """BAS options can be authored from scratch and survive load/write."""
    from flopy.mfusg import MfUsgBas
    from flopy.modflow import ModflowDis

    ml = MfUsg(structured=False, model_ws=str(function_tmpdir))
    ModflowDis(ml, nlay=1, nrow=1, ncol=1, nper=1)
    bas = MfUsgBas(
        ml,
        ibound=1,
        strt=10.0,
        structured=False,
        ifrefm=True,
        iprintfv=True,
        iprinttime=True,
        converge=True,
        richards=True,
        double_prec=True,
        double_out=True,
        double_io=True,
        sy_all=True,
        ishowp=True,
        stoper=0.01,
    )
    bas.fn_path = str(function_tmpdir / "created.bas")
    bas.write_file(check=False)

    option_line = Path(bas.fn_path).read_text().splitlines()[1]
    for token in (
        "PRINTFV",
        "CONVERGE",
        "UNSTRUCTURED",
        "FREE",
        "PRINTTIME",
        "SHOWPROGRESS",
        "RICHARDS",
        "DPIN",
        "DPOUT",
        "DPIO",
        "SY-ALL",
        "STOPERROR",
    ):
        assert token in option_line

    ml2 = MfUsg(structured=False)
    ModflowDis(ml2, nlay=1, nrow=1, ncol=1, nper=1)
    bas2 = MfUsgBas.load(str(bas.fn_path), ml2, check=False)
    assert bas2.iprintfv
    assert bas2.iprinttime
    assert bas2.converge
    assert bas2.richards
    assert bas2.double_prec
    assert bas2.double_out
    assert bas2.double_io
    assert bas2.sy_all
    assert bas2.ishowp
    assert np.isclose(bas2.stoper, 0.01)
    assert not bas2.structured


# ============================================================
# Tests for fork additions: MfUsgChd, MfUsgRiv, MfUsgEts,
#                          MfusgTransportListBudget
# ============================================================


def test_mfusgchd_roundtrip(function_tmpdir):
    """MfUsgChd load + write round-trip on a minimal unstructured CHD file.

    Checks:
    - SP with data loads correct node / shead / ehead values
    - SP with -1 (reuse) copies previous SP's data
    - Write emits the correct number of data rows per SP
    """
    from flopy.mfusg import MfUsgChd
    from flopy.modflow import ModflowDis

    chd_in = function_tmpdir / "test.chd"
    chd_in.write_text(
        "# MfUsgChd test\n"
        "         3\n"
        " 3    Stress Period 1\n"
        " 101   15.000000  14.500000\n"
        " 102   20.000000  19.500000\n"
        " 103   25.000000  24.500000\n"
        " -1    Stress Period 2\n"
    )

    ml = MfUsg(structured=False, model_ws=str(function_tmpdir))
    ModflowDis(ml, nlay=1, nrow=1, ncol=1, nper=2)
    chd = MfUsgChd.load(str(chd_in), ml, nper=2, ext_unit_dict={})

    sp0 = chd.stress_period_data[0]
    assert len(sp0) == 3
    assert list(sp0["node"]) == [100, 101, 102]
    assert np.isclose(sp0["shead"][0], 15.0, atol=0.01)
    assert np.isclose(sp0["ehead"][2], 24.5, atol=0.01)

    # Reuse SP copies previous data
    sp1 = chd.stress_period_data[1]
    assert len(sp1) == 3
    assert list(sp1["node"]) == [100, 101, 102]

    # Write and verify both SPs appear in output
    chd_out = function_tmpdir / "out.chd"
    chd.fn_path = str(chd_out)
    chd.write_file()
    text = chd_out.read_text()
    assert "Stress Period 1" in text
    assert "Stress Period 2" in text
    assert "\n 101" in text

    # Re-load and verify data survives round-trip
    ml2 = MfUsg(structured=False)
    ModflowDis(ml2, nlay=1, nrow=1, ncol=1, nper=2)
    chd2 = MfUsgChd.load(str(chd_out), ml2, nper=2, ext_unit_dict={})
    sp0b = chd2.stress_period_data[0]
    assert np.isclose(sp0b["shead"][0], 15.0, atol=0.01)
    assert list(sp0b["node"]) == [100, 101, 102]


def test_mfusgchd_aux_roundtrip(function_tmpdir):
    """MfUsgChd with AUX concentration variable: load + write round-trip."""
    from flopy.mfusg import MfUsgChd
    from flopy.modflow import ModflowDis

    chd_in = function_tmpdir / "aux.chd"
    chd_in.write_text(
        "# MfUsgChd aux test\n"
        "         2 AUX C01\n"
        " 2    Stress Period 1\n"
        " 101   15.000000  14.500000  1.000000e-01\n"
        " 102   20.000000  19.500000  2.000000e-01\n"
        " -1    Stress Period 2\n"
    )

    ml = MfUsg(structured=False, model_ws=str(function_tmpdir))
    ModflowDis(ml, nlay=1, nrow=1, ncol=1, nper=2)
    chd = MfUsgChd.load(str(chd_in), ml, nper=2, ext_unit_dict={})

    # Field name is preserved as-is from file ('C01', not lowercased)
    assert "C01" in chd.dtype.names
    sp0 = chd.stress_period_data[0]
    assert np.isclose(sp0["C01"][0], 0.1, atol=1e-5)
    assert np.isclose(sp0["C01"][1], 0.2, atol=1e-5)

    # Write and verify AUX keyword appears in header
    chd_out = function_tmpdir / "out_aux.chd"
    chd.fn_path = str(chd_out)
    chd.write_file()
    lines = chd_out.read_text().splitlines()
    header = next(l for l in lines if not l.startswith("#"))
    assert "AUX" in header.upper()

    # Values survive write → re-load
    ml2 = MfUsg(structured=False)
    ModflowDis(ml2, nlay=1, nrow=1, ncol=1, nper=2)
    chd2 = MfUsgChd.load(str(chd_out), ml2, nper=2, ext_unit_dict={})
    assert np.isclose(chd2.stress_period_data[0]["C01"][0], 0.1, atol=1e-4)


def test_mfusgriv_roundtrip(function_tmpdir):
    """MfUsgRiv load + write: IRDFLAG, node/stage/cond/rbot, reuse SP."""
    from flopy.mfusg import MfUsgRiv
    from flopy.modflow import ModflowDis

    riv_in = function_tmpdir / "test.riv"
    riv_in.write_text(
        "# MfUsgRiv test\n"
        " 2 50\n"
        " 2 0    Stress Period 1\n"
        " 101  15.000000  1.000000e-04  12.000000\n"
        " 102  20.000000  2.000000e-04  17.000000\n"
        " -1 0    Stress Period 2\n"
    )

    ml = MfUsg(structured=False, model_ws=str(function_tmpdir))
    ModflowDis(ml, nlay=1, nrow=1, ncol=1, nper=2)
    riv = MfUsgRiv.load(str(riv_in), ml, nper=2, ext_unit_dict={})

    assert riv.irdflag == 50
    sp0 = riv.stress_period_data[0]
    assert len(sp0) == 2
    assert sp0["node"][0] == 100
    assert np.isclose(sp0["stage"][0], 15.0, atol=0.01)
    assert np.isclose(sp0["rbot"][1], 17.0, atol=0.01)

    # Reuse SP copies previous data
    sp1 = riv.stress_period_data[1]
    assert len(sp1) == 2
    assert sp1["node"][0] == 100

    # Write and verify both SPs in output, then re-load
    riv_out = function_tmpdir / "out.riv"
    riv.fn_path = str(riv_out)
    riv.write_file()
    text = riv_out.read_text()
    assert "Stress Period 1" in text
    assert "Stress Period 2" in text
    assert "\n 101" in text

    ml2 = MfUsg(structured=False)
    ModflowDis(ml2, nlay=1, nrow=1, ncol=1, nper=2)
    riv2 = MfUsgRiv.load(str(riv_out), ml2, nper=2, ext_unit_dict={})
    assert np.isclose(riv2.stress_period_data[0]["stage"][0], 15.0, atol=0.01)


def test_mfusgriv_irch_detection(function_tmpdir):
    """MfUsgRiv auto-detects trailing irch column on load and preserves it on write."""
    from flopy.mfusg import MfUsgRiv
    from flopy.modflow import ModflowDis

    riv_in = function_tmpdir / "irch.riv"
    riv_in.write_text(
        "# MfUsgRiv irch test\n"
        " 2 50\n"
        " 2 0    Stress Period 1\n"
        " 101  15.000000  1.000000e-04  12.000000  5\n"
        " 102  20.000000  2.000000e-04  17.000000  7\n"
        " -1 0    Stress Period 2\n"
    )

    ml = MfUsg(structured=False, model_ws=str(function_tmpdir))
    ModflowDis(ml, nlay=1, nrow=1, ncol=1, nper=2)
    riv = MfUsgRiv.load(str(riv_in), ml, nper=2, ext_unit_dict={})

    assert "irch" in riv.dtype.names
    sp0 = riv.stress_period_data[0]
    assert sp0["irch"][0] == 5
    assert sp0["irch"][1] == 7

    # irch must NOT appear as AUX in options
    assert not any("irch" in o.lower() for o in riv.options)

    # Round-trip: irch values survive write → load
    riv_out = function_tmpdir / "out_irch.riv"
    riv.fn_path = str(riv_out)
    riv.write_file()

    ml2 = MfUsg(structured=False)
    ModflowDis(ml2, nlay=1, nrow=1, ncol=1, nper=2)
    riv2 = MfUsgRiv.load(str(riv_out), ml2, nper=2, ext_unit_dict={})
    assert "irch" in riv2.dtype.names
    assert riv2.stress_period_data[0]["irch"][0] == 5
    assert riv2.stress_period_data[0]["irch"][1] == 7


def test_mfusgriv_aux_and_irch(function_tmpdir):
    """MfUsgRiv with both AUX concentration and trailing irch."""
    from flopy.mfusg import MfUsgRiv

    riv_in = function_tmpdir / "aux_irch.riv"
    riv_in.write_text(
        "# MfUsgRiv aux+irch test\n"
        " 1 50 AUX C01\n"
        " 1 0    Stress Period 1\n"
        " 101  15.000000  1.000000e-04  12.000000  5.000000e-02  3\n"
    )

    ml = MfUsg(structured=False)
    riv = MfUsgRiv.load(str(riv_in), ml, nper=1, ext_unit_dict={})

    # Field names are preserved as-is from file header (uppercase 'C01')
    assert "C01" in riv.dtype.names
    assert "irch" in riv.dtype.names
    sp0 = riv.stress_period_data[0]
    assert np.isclose(sp0["C01"][0], 0.05, atol=1e-5)
    assert sp0["irch"][0] == 3


def test_mfusgets_construction(function_tmpdir):
    """MfUsgEts can be constructed and attached to a structured model."""
    from flopy.mfusg import MfUsgEts
    from flopy.modflow import ModflowDis

    ml = MfUsg(structured=True, model_ws=str(function_tmpdir))
    ModflowDis(ml, nlay=1, nrow=2, ncol=2, nper=1)
    ets = MfUsgEts(ml, netsop=1, evtr=1.2e-4, netseg=1)

    assert ml.ets is not None
    assert ml.ets.netseg == 1
    assert ml.ets.netsop == 1


def test_mfusgets_write(function_tmpdir):
    """MfUsgEts with netseg=2 writes a file with PXDP/PETM segment arrays."""
    from flopy.mfusg import MfUsgEts
    from flopy.modflow import ModflowDis

    ml = MfUsg(structured=True, model_ws=str(function_tmpdir))
    ModflowDis(ml, nlay=1, nrow=2, ncol=2, nper=1)
    ets = MfUsgEts(ml, netsop=1, evtr=1.2e-4, netseg=2, pxdp=[0.5], petm=[0.5])
    ml.write_input()

    ets_file = Path(ets.fn_path)
    assert ets_file.exists(), f"ETS file not written at {ets_file}"
    content = ets_file.read_text()
    # netseg=2 → PXDP and PETM segment arrays must be written
    assert "pxdp" in content.lower() and "petm" in content.lower()


def test_mfusgets_parameterized_write_fails_explicitly(function_tmpdir):
    """Programmatic ETS parameters are not silently written as incomplete files."""
    from flopy.mfusg import MfUsgEts
    from flopy.modflow import ModflowDis

    ml = MfUsg(structured=True, model_ws=str(function_tmpdir))
    ModflowDis(ml, nlay=1, nrow=2, ncol=2, nper=1)
    ets = MfUsgEts(ml, netsop=1, evtr=1.2e-4, npets=1)

    with pytest.raises(NotImplementedError, match="parameter"):
        ets.write_file()


def test_modflow_name_file_preserves_input_external_paths(function_tmpdir):
    """NAM writing preserves external input subdirs and rebases output files."""
    from flopy.modflow import Modflow

    ml = Modflow(modelname="nam_paths", model_ws=str(function_tmpdir))
    ml.add_external("arrays/recharge.ref", unit=101, binflag=False, output=False)
    ml.add_external("arrays/binary.ref", unit=102, binflag=True, output=False)
    ml.add_external(
        "old_external_outputs/heads.hds", unit=103, binflag=True, output=True
    )
    ml.add_output("old_outputs/nam_paths.cbc", unit=201, binflag=True)
    ml.write_name_file()

    entries = {}
    for line in (function_tmpdir / "nam_paths.nam").read_text().splitlines():
        parts = line.split()
        if parts and parts[0].startswith("DATA"):
            entries[int(parts[1])] = parts[2]

    assert entries[101] == "arrays/recharge.ref"
    assert entries[102] == "arrays/binary.ref"
    assert entries[103] == "heads.hds"
    assert entries[201] == "nam_paths.cbc"


def test_mfusg_boundary_programmatic_creation_uses_zero_based_nodes(function_tmpdir):
    """USG-T boundary packages can be authored from scratch with 0-based nodes."""
    from flopy.mfusg import MfUsgChd, MfUsgDrn, MfUsgGhb, MfUsgRiv
    from flopy.modflow import ModflowDis

    ml = MfUsg(structured=False, model_ws=str(function_tmpdir))
    ModflowDis(ml, nlay=1, nrow=1, ncol=1, nper=2)

    chd_dtype = MfUsgChd.get_default_dtype(structured=False)
    chd_data = {
        0: np.array(
            [(0, 15.0, 14.5), (4, 20.0, 19.5)],
            dtype=chd_dtype,
        ).view(np.recarray)
    }
    chd = MfUsgChd(ml, stress_period_data=chd_data)
    chd.fn_path = str(function_tmpdir / "created.chd")
    chd.write_file()
    assert "\n 1   15.000000" in Path(chd.fn_path).read_text()
    assert "\n 5   20.000000" in Path(chd.fn_path).read_text()

    riv_dtype = np.dtype(
        [
            ("node", int),
            ("stage", np.float64),
            ("cond", np.float32),
            ("rbot", np.float32),
            ("C01", np.float32),
            ("irch", int),
        ]
    )
    riv = MfUsgRiv(
        ml,
        stress_period_data={
            0: np.array(
                [(2, 10.0, 1.0e-4, 9.0, 0.25, 7)],
                dtype=riv_dtype,
            ).view(np.recarray)
        },
        dtype=riv_dtype,
    )
    riv.fn_path = str(function_tmpdir / "created.riv")
    riv.write_file()
    riv_text = Path(riv.fn_path).read_text()
    assert "AUX C01" in riv_text
    assert "\n 3  10.000000" in riv_text
    assert riv_text.splitlines()[3].split()[-1] == "7"

    ghb_dtype = np.dtype(
        [
            ("node", int),
            ("bhead", np.float32),
            ("cond", np.float32),
            ("C01", np.float32),
        ]
    )
    ghb = MfUsgGhb(
        ml,
        stress_period_data={
            0: np.array([(8, 5.0, 1.0e2, 0.15)], dtype=ghb_dtype).view(np.recarray)
        },
        dtype=ghb_dtype,
    )
    ghb.fn_path = str(function_tmpdir / "created.ghb")
    ghb.write_file()
    ghb_text = Path(ghb.fn_path).read_text()
    assert "AUX C01" in ghb_text
    assert "\n 9  5.000000" in ghb_text

    drn_dtype = MfUsgDrn.get_default_dtype(structured=False)
    drn = MfUsgDrn(
        ml,
        stress_period_data={
            0: np.array([(10, 2.5, 1.0e1)], dtype=drn_dtype).view(np.recarray)
        },
    )
    drn.fn_path = str(function_tmpdir / "created.drn")
    drn.write_file()
    assert "\n 11  2.500000" in Path(drn.fn_path).read_text()

    ml2 = MfUsg(structured=False)
    ModflowDis(ml2, nlay=1, nrow=1, ncol=1, nper=2)
    created_chd = MfUsgChd.load(chd.fn_path, ml2, nper=2)
    assert list(created_chd.stress_period_data[0]["node"]) == [0, 4]


def test_mfusgwel_programmatic_cln_itmpcln_roundtrip(function_tmpdir):
    """WEL authoring writes the USG-T ITMP/NP/ITMPCLN stress-period header."""
    from flopy.mfusg import MfUsgWel
    from flopy.modflow import ModflowDis

    ml = MfUsg(structured=False, model_ws=str(function_tmpdir))
    ModflowDis(ml, nlay=1, nrow=1, ncol=1, nper=3)

    dtype = np.dtype([("node", int), ("flux", np.float32), ("C01", np.float32)])
    stress_period_data = {
        0: np.array([(0, -100.0, 0.10)], dtype=dtype).view(np.recarray)
    }
    cln_stress_period_data = {
        0: np.array([(2, -10.0, 0.20)], dtype=dtype).view(np.recarray),
        1: np.array([(3, -20.0, 0.30)], dtype=dtype).view(np.recarray),
    }
    wel = MfUsgWel(
        ml,
        stress_period_data=stress_period_data,
        cln_stress_period_data=cln_stress_period_data,
        dtype=dtype,
        cln_dtype=dtype,
    )
    wel.fn_path = str(function_tmpdir / "created.wel")
    wel.write_file()

    lines = [
        line
        for line in Path(wel.fn_path).read_text().splitlines()
        if line.strip() and not line.lstrip().startswith("#")
    ]
    assert "aux C01" in lines[0]
    assert lines[1].split()[:3] == ["1", "0", "1"]
    t = lines[2].split()
    assert t[0] == "1"
    assert np.isclose(float(t[1]), -100.0)
    assert np.isclose(float(t[2]), 0.10)
    t = lines[3].split()
    assert t[0] == "3"
    assert np.isclose(float(t[1]), -10.0)
    assert np.isclose(float(t[2]), 0.20)
    assert lines[4].split()[:3] == ["0", "0", "1"]
    t = lines[5].split()
    assert t[0] == "4"
    assert np.isclose(float(t[1]), -20.0)
    assert np.isclose(float(t[2]), 0.30)
    assert lines[6].split()[:3] == ["0", "0", "0"]

    ml2 = MfUsg(structured=False)
    ModflowDis(ml2, nlay=1, nrow=1, ncol=1, nper=3)
    wel2 = MfUsgWel.load(str(wel.fn_path), ml2, nper=3, check=False)
    assert list(wel2.stress_period_data[0]["node"]) == [0]
    assert list(wel2.cln_stress_period_data[0]["node"]) == [2]
    assert list(wel2.cln_stress_period_data[1]["node"]) == [3]
    assert np.isclose(wel2.cln_stress_period_data[1]["c01"][0], 0.30)


def test_mfusgcln_general_section_processccf_roundtrip(function_tmpdir):
    """CLN supports USG-T PROCESSCCF, GENERAL_SEC, and ISHAPE authoring."""
    from flopy.modflow import ModflowDis

    ml = MfUsg(model_ws=str(function_tmpdir))
    ModflowDis(ml, nlay=1, nrow=1, ncol=1, nper=1)

    cln = MfUsgCln(
        ml,
        ncln=1,
        iclnnds=-1,
        nndcln=1,
        nclngwc=1,
        node_prop=[[1, 3, 1, 0, 10.0, -1.0, 0.0, 0, 0]],
        cln_gwc=[[1, 1, 1, 1, 0, 0.0, 10.0, 1.0, 0]],
        nconduityp=1,
        cln_circ=[[1, 0.5, 100.0]],
        processccf=True,
        iclngwcb=902,
        ngenshptyp=1,
        ngentabrows=2,
        cln_gen=[(1, 3.0, [[0.0, 0.0, 0.0, 0.0], [1.0, 2.0, 3.0, 4.0]])],
    )
    cln.fn_path = str(function_tmpdir / "created.cln")
    cln.write_file()

    lines = [
        line
        for line in Path(cln.fn_path).read_text().splitlines()
        if line.strip() and not line.lstrip().startswith("#")
    ]
    assert "PROCESSCCF 902" in lines[0]
    assert "GENERAL_SEC 1 2" in lines[1]
    assert any(line.split()[:3] == ["1", "3", "1"] for line in lines)

    ml2 = MfUsg()
    ModflowDis(ml2, nlay=1, nrow=1, ncol=1, nper=1)
    cln2 = MfUsgCln.load(str(cln.fn_path), ml2)
    assert cln2.processccf
    assert cln2.iclngwcb == 902
    assert cln2.ngenshptyp == 1
    assert cln2.ngentabrows == 2
    assert "ishape" in cln2.node_prop.dtype.names
    assert cln2.node_prop["ishape"][0] == 3
    assert cln2.cln_gen[0][0] == 1
    assert np.allclose(cln2.cln_gen[0][2], [[0, 0, 0, 0], [1, 2, 3, 4]])


def test_mfusgdpf_programmatic_tabrich_and_sc2im_roundtrip(function_tmpdir):
    """DPF authoring honors FRAHK, IUZONTABIM, IDPF, and conditional SC2IM."""
    from flopy.modflow import ModflowDis

    ml = MfUsg(model_ws=str(function_tmpdir))
    ModflowDis(ml, nlay=2, nrow=1, ncol=1, nper=1, steady=False)
    MfUsgBas(ml, ibound=1, strt=1.0, richards=True)
    MfUsgBcf(
        ml,
        laycon=[0, 5],
        ipakcb=0,
        tabrich=True,
        nuzones=1,
        nutabrows=2,
        tran=1.0,
        hy=1.0,
        kv=1.0,
        sf1=1.0e-5,
        sf2=0.15,
    )

    dpf = MfUsgDpf(
        ml,
        frahk=True,
        iuzontabim=[7, 8],
        iboundim=1,
        hnewim=2.0,
        phif=0.1,
        ddftr=0.01,
        sc1im=1.0e-5,
        sc2im=[[[0.0]], [[0.25]]],
    )
    assert ml.idpf == 1
    dpf.fn_path = str(function_tmpdir / "created.dpf")
    dpf.write_file()

    text = Path(dpf.fn_path).read_text()
    assert "FRAHK" in text.splitlines()[0]
    assert "#iuzontabim" in text
    assert "#sc2im layer 1" not in text
    assert "#sc2im layer 2" in text

    ml2 = MfUsg()
    ModflowDis(ml2, nlay=2, nrow=1, ncol=1, nper=1, steady=False)
    MfUsgBas(ml2, ibound=1, strt=1.0, richards=True)
    MfUsgBcf(
        ml2,
        laycon=[0, 5],
        ipakcb=0,
        tabrich=True,
        nuzones=1,
        nutabrows=2,
        tran=1.0,
        hy=1.0,
        kv=1.0,
        sf1=1.0e-5,
        sf2=0.15,
    )
    dpf2 = MfUsgDpf.load(str(dpf.fn_path), ml2)
    assert ml2.idpf == 1
    assert dpf2.frahk
    assert list(dpf2.iuzontabim.array) == [7, 8]
    assert np.allclose(dpf2.sc2im[0].array, 0.0)
    assert np.allclose(dpf2.sc2im[1].array, 0.25)


def test_mfusgdpf_programmatic_richards_immobile_roundtrip(function_tmpdir):
    """DPF writes immobile Richards arrays for LAYCON=5 when TABRICH is off."""
    from flopy.modflow import ModflowDis

    ml = MfUsg(model_ws=str(function_tmpdir))
    ModflowDis(ml, nlay=1, nrow=1, ncol=1, nper=1, steady=False)
    MfUsgBas(ml, ibound=1, strt=1.0, richards=True)
    MfUsgBcf(
        ml,
        laycon=[5],
        ipakcb=0,
        bubblept=True,
        tran=1.0,
        hy=1.0,
        kv=1.0,
        sf1=1.0e-5,
        sf2=0.15,
    )

    dpf = MfUsgDpf(
        ml,
        iboundim=1,
        hnewim=2.0,
        phif=0.1,
        ddftr=0.01,
        sc1im=1.0e-5,
        sc2im=0.25,
        alphaim=0.2,
        betaim=4.0,
        srim=0.1,
        brookim=5.0,
        bpim=-0.25,
    )
    dpf.fn_path = str(function_tmpdir / "richards.dpf")
    dpf.write_file()

    text = Path(dpf.fn_path).read_text()
    assert "#alphaim layer 1" in text
    assert "#betaim layer 1" in text
    assert "#srim layer 1" in text
    assert "#brookim layer 1" in text
    assert "#bpim layer 1" in text

    ml2 = MfUsg()
    ModflowDis(ml2, nlay=1, nrow=1, ncol=1, nper=1, steady=False)
    MfUsgBas(ml2, ibound=1, strt=1.0, richards=True)
    MfUsgBcf(
        ml2,
        laycon=[5],
        ipakcb=0,
        bubblept=True,
        tran=1.0,
        hy=1.0,
        kv=1.0,
        sf1=1.0e-5,
        sf2=0.15,
    )
    dpf2 = MfUsgDpf.load(str(dpf.fn_path), ml2)
    assert np.allclose(dpf2.alphaim[0].array, 0.2)
    assert np.allclose(dpf2.betaim[0].array, 4.0)
    assert np.allclose(dpf2.srim[0].array, 0.1)
    assert np.allclose(dpf2.brookim[0].array, 5.0)
    assert np.allclose(dpf2.bpim[0].array, -0.25)


def test_mfusghfb_static_fortran_layout(function_tmpdir):
    """Static non-parametric HFB writes only NHFBNP rows after the header."""
    from flopy.mfusg import MfUsgHfb
    from flopy.modflow import ModflowDis

    dtype = MfUsgHfb.get_default_dtype(structured=False)
    hfb_data = np.array(
        [(0, 1, 1.0e-4), (2, 3, 2.0e-4)],
        dtype=dtype,
    ).view(np.recarray)

    ml = MfUsg(structured=False, model_ws=str(function_tmpdir))
    ModflowDis(ml, nlay=1, nrow=1, ncol=1, nper=3)
    hfb = MfUsgHfb(ml, hfb_data=hfb_data, transient=False)

    hfb_out = function_tmpdir / "static.hfb"
    hfb.fn_path = str(hfb_out)
    hfb.write_file()

    lines = [
        line
        for line in hfb_out.read_text().splitlines()
        if line.strip() and not line.lstrip().startswith("#")
    ]
    assert len(lines) == 3
    assert lines[0].split()[:3] == ["0", "0", "2"]
    assert lines[1].split()[:2] == ["1", "2"]
    assert lines[2].split()[:2] == ["3", "4"]

    ml2 = MfUsg(structured=False)
    ModflowDis(ml2, nlay=1, nrow=1, ncol=1, nper=3)
    hfb2 = MfUsgHfb.load(str(hfb_out), ml2, nper=3, ext_unit_dict={})
    assert list(hfb2.hfb_data["node1"]) == [0, 2]
    assert list(hfb2.hfb_data["node2"]) == [1, 3]


def test_mfusghfb_structured_static_fortran_layout(function_tmpdir):
    """Structured HFB keeps k/i/j internal indices 0-based and writes 1-based."""
    from flopy.mfusg import MfUsgHfb
    from flopy.modflow import ModflowDis

    dtype = MfUsgHfb.get_default_dtype(structured=True)
    hfb_data = np.array(
        [(0, 0, 0, 0, 1, 1.0e-4)],
        dtype=dtype,
    ).view(np.recarray)

    ml = MfUsg(structured=True, model_ws=str(function_tmpdir))
    ModflowDis(ml, nlay=1, nrow=1, ncol=2, nper=1)
    hfb = MfUsgHfb(ml, hfb_data=hfb_data, transient=False)

    hfb_out = function_tmpdir / "structured.hfb"
    hfb.fn_path = str(hfb_out)
    hfb.write_file()

    lines = [
        line
        for line in hfb_out.read_text().splitlines()
        if line.strip() and not line.lstrip().startswith("#")
    ]
    assert lines[1].split()[:5] == ["1", "1", "1", "1", "2"]

    ml2 = MfUsg(structured=True)
    ModflowDis(ml2, nlay=1, nrow=1, ncol=2, nper=1)
    hfb2 = MfUsgHfb.load(str(hfb_out), ml2, nper=1, ext_unit_dict={})
    assert hfb2.hfb_data[0]["k"] == 0
    assert hfb2.hfb_data[0]["icol2"] == 1


def test_mfusghfb_transient_fortran_layout(function_tmpdir):
    """Transient HFB uses IHFBRD as a flag and -1 for stress-period reuse."""
    from flopy.mfusg import MfUsgHfb
    from flopy.modflow import ModflowDis

    dtype = MfUsgHfb.get_default_dtype(structured=False)
    sp0 = np.array(
        [(0, 1, 1.0e-4), (2, 3, 2.0e-4)],
        dtype=dtype,
    ).view(np.recarray)
    sp1 = np.array(
        [(4, 5, 3.0e-4), (6, 7, 4.0e-4)],
        dtype=dtype,
    ).view(np.recarray)

    ml = MfUsg(structured=False, model_ws=str(function_tmpdir))
    ModflowDis(ml, nlay=1, nrow=1, ncol=1, nper=3)
    hfb = MfUsgHfb(
        ml,
        hfb_data=sp0,
        transient=True,
        stress_period_data={1: sp1},
    )

    hfb_out = function_tmpdir / "transient.hfb"
    hfb.fn_path = str(hfb_out)
    hfb.write_file()

    lines = [
        line
        for line in hfb_out.read_text().splitlines()
        if line.strip() and not line.lstrip().startswith("#")
    ]
    assert "TRANSIENT_HFB" in lines[0]
    assert lines[1].split() == ["1"]
    assert lines[2].split()[:2] == ["1", "2"]
    assert lines[4].split() == ["1"]
    assert lines[5].split()[:2] == ["5", "6"]
    assert lines[7].split() == ["-1"]

    ml2 = MfUsg(structured=False)
    ModflowDis(ml2, nlay=1, nrow=1, ncol=1, nper=3)
    hfb2 = MfUsgHfb.load(str(hfb_out), ml2, nper=3, ext_unit_dict={})
    assert hfb2.transient
    assert set(hfb2.stress_period_data) == {0, 1}
    assert list(hfb2.stress_period_data[0]["node1"]) == [0, 2]
    assert list(hfb2.stress_period_data[1]["node1"]) == [4, 6]

    hfb_zero = function_tmpdir / "transient_zero.hfb"
    hfb_zero.write_text(
        "# HFB with IHFBRD=0 reuse\n"
        "         0         0         2  TRANSIENT_HFB\n"
        "1\n"
        "1 2 1.0e-4\n"
        "3 4 2.0e-4\n"
        "0\n"
        "-1\n"
    )
    ml3 = MfUsg(structured=False)
    ModflowDis(ml3, nlay=1, nrow=1, ncol=1, nper=3)
    hfb3 = MfUsgHfb.load(str(hfb_zero), ml3, nper=3, ext_unit_dict={})
    assert set(hfb3.stress_period_data) == {0}
    assert list(hfb3.stress_period_data[0]["node1"]) == [0, 2]


def _make_usgt_lst_old_format():
    """Return a minimal valid USG-T listing file (old format: VOLUMETRIC BUDGET
    for both flow and transport) with 1 SP and 2 transport species."""
    budget_block = (
        " {bkey} AT END OF TIME STEP    1 STRESS PERIOD   1\n"
        " IN:\n"
        " WELLS =    {val:.6E}     WELLS =    {val:.6E}\n"
        " -------------------------------------------------------\n"
        " OUT:\n"
        " WELLS =    {val:.6E}     WELLS =    {val:.6E}\n"
        " PERCENT DISCREPANCY =    0.000000E+00"
        "     PERCENT DISCREPANCY =    0.000000E+00\n"
    )
    time_block = (
        " TIME SUMMARY AT END OF TIME STEP    1 IN STRESS PERIOD    1\n"
        "                SECONDS     MINUTES      HOURS       DAYS        YEARS\n"
        " -----------------------------------------------------------------------\n"
        " TIME STEP          86400.     1440.0     24.000     1.0000    2.7397E-03\n"
        " STRESS PERIOD      86400.     1440.0     24.000     1.0000    2.7397E-03\n"
        " TOTAL              86400.     1440.0     24.000     1.0000    2.7397E-03\n"
    )

    lines = []
    # Flow section
    lines.append(" IN FLOW TIME STEP    1   STRESS PERIOD    1\n")
    lines.append("\n")
    lines.append(budget_block.format(bkey="VOLUMETRIC BUDGET FOR ENTIRE MODEL", val=100.0))
    lines.append("\n")
    lines.append(time_block)
    lines.append("\n")

    # Transport section
    lines.append(" TRANSPORT SOLUTION COMPLETE FOR ALL SPECIES\n")
    for sp_num, val in [(1, 0.1), (2, 0.2)]:
        lines.append(f" TRANSPORT OUTPUT FOR COMPONENT SPECIES NUMBER  {sp_num:4d}\n")
        lines.append(budget_block.format(bkey="VOLUMETRIC BUDGET FOR ENTIRE MODEL", val=val))
        lines.append("\n")
        lines.append(time_block)
        lines.append("\n")

    return "".join(lines)


def _make_usgt_lst_new_format():
    """Return a minimal valid USG-T listing file (new format: MASS BUDGET for
    transport) with 1 SP and 2 transport species."""
    flow_budget_block = (
        " VOLUMETRIC BUDGET FOR ENTIRE MODEL AT END OF TIME STEP    1 STRESS PERIOD   1\n"
        " IN:\n"
        " WELLS =    1.000000E+02     WELLS =    1.000000E+02\n"
        " -------------------------------------------------------\n"
        " OUT:\n"
        " WELLS =    1.000000E+02     WELLS =    1.000000E+02\n"
        " PERCENT DISCREPANCY =    0.000000E+00"
        "     PERCENT DISCREPANCY =    0.000000E+00\n"
    )
    mass_budget_block = (
        " MASS BUDGET FOR ENTIRE MODEL AT END OF TIME STEP    1 STRESS PERIOD   1\n"
        " IN:\n"
        " WELLS =    {val:.6E}     WELLS =    {val:.6E}\n"
        " -------------------------------------------------------\n"
        " OUT:\n"
        " WELLS =    {val:.6E}     WELLS =    {val:.6E}\n"
        " PERCENT DISCREPANCY =    0.000000E+00"
        "     PERCENT DISCREPANCY =    0.000000E+00\n"
    )
    time_block = (
        " TIME SUMMARY AT END OF TIME STEP    1 IN STRESS PERIOD    1\n"
        "                SECONDS     MINUTES      HOURS       DAYS        YEARS\n"
        " -----------------------------------------------------------------------\n"
        " TIME STEP          86400.     1440.0     24.000     1.0000    2.7397E-03\n"
        " STRESS PERIOD      86400.     1440.0     24.000     1.0000    2.7397E-03\n"
        " TOTAL              86400.     1440.0     24.000     1.0000    2.7397E-03\n"
    )

    lines = []
    lines.append(" IN FLOW TIME STEP    1   STRESS PERIOD    1\n")
    lines.append("\n")
    lines.append(flow_budget_block)
    lines.append("\n")
    lines.append(time_block)
    lines.append("\n")
    lines.append(" TRANSPORT SOLUTION COMPLETE FOR ALL SPECIES\n")
    for sp_num, val in [(1, 0.1), (2, 0.2)]:
        lines.append(f" TRANSPORT OUTPUT FOR COMPONENT SPECIES NUMBER  {sp_num:4d}\n")
        lines.append(mass_budget_block.format(val=val))
        lines.append("\n")
        lines.append(time_block)
        lines.append("\n")

    return "".join(lines)


def test_mfusg_transport_list_budget_old_format(function_tmpdir):
    """MfusgTransportListBudget reads species budgets from old-format LST.

    Old format: transport blocks use 'VOLUMETRIC BUDGET', same keyword as flow.
    State machine must use in_transport + current_species to disambiguate.
    """
    from flopy.utils import MfusgTransportListBudget

    lst_path = function_tmpdir / "old.lst"
    lst_path.write_text(_make_usgt_lst_old_format())

    # Species 1 — 1 SP → idx_map length must be 1
    reader1 = MfusgTransportListBudget(str(lst_path), species=1)
    assert len(reader1.idx_map) == 1, (
        f"Expected 1 budget block for species 1, got {len(reader1.idx_map)}"
    )
    inc1, cum1 = reader1.get_budget()
    assert len(inc1) == 1
    assert np.isclose(inc1["WELLS_IN"][0], 0.1, atol=1e-4)

    # Species 2 — independent reader
    reader2 = MfusgTransportListBudget(str(lst_path), species=2)
    assert len(reader2.idx_map) == 1
    inc2, cum2 = reader2.get_budget()
    assert np.isclose(inc2["WELLS_IN"][0], 0.2, atol=1e-4)


def test_mfusg_transport_list_budget_new_format(function_tmpdir):
    """MfusgTransportListBudget reads species budgets from new-format LST.

    New format: transport blocks use 'MASS BUDGET' keyword.
    The flow block still uses 'VOLUMETRIC BUDGET' and must not be indexed.
    """
    from flopy.utils import MfusgTransportListBudget

    lst_path = function_tmpdir / "new.lst"
    lst_path.write_text(_make_usgt_lst_new_format())

    reader1 = MfusgTransportListBudget(str(lst_path), species=1)
    assert len(reader1.idx_map) == 1
    inc1, _ = reader1.get_budget()
    assert np.isclose(inc1["WELLS_IN"][0], 0.1, atol=1e-4)

    reader2 = MfusgTransportListBudget(str(lst_path), species=2)
    assert len(reader2.idx_map) == 1
    inc2, _ = reader2.get_budget()
    assert np.isclose(inc2["WELLS_IN"][0], 0.2, atol=1e-4)


def test_mfusg_transport_list_budget_species_isolation(function_tmpdir):
    """Each species reader returns only its own blocks; neither reads the other's."""
    from flopy.utils import MfusgTransportListBudget

    # Use old format (more complex state machine)
    lst_path = function_tmpdir / "isolation.lst"
    lst_path.write_text(_make_usgt_lst_old_format())

    reader1 = MfusgTransportListBudget(str(lst_path), species=1)
    reader2 = MfusgTransportListBudget(str(lst_path), species=2)

    inc1, _ = reader1.get_budget()
    inc2, _ = reader2.get_budget()

    # Values must differ — each reader got its own block
    assert not np.isclose(inc1["WELLS_IN"][0], inc2["WELLS_IN"][0], atol=1e-4), (
        "Species 1 and 2 readers returned identical values — isolation failed"
    )
    assert np.isclose(inc1["WELLS_IN"][0], 0.1, atol=1e-4)
    assert np.isclose(inc2["WELLS_IN"][0], 0.2, atol=1e-4)


# ---------------------------------------------------------------------------
# MfUsgGhb tests
# ---------------------------------------------------------------------------

def test_mfusgghb_roundtrip(function_tmpdir):
    """Basic GHB: load → inspect → write → reload preserves data."""
    from flopy.modflow import ModflowDis

    ghb_in = function_tmpdir / "test.ghb"
    ghb_in.write_text(
        "# MfUsgGhb test\n"
        "         2 0\n"
        " 2 0    Stress Period 1\n"
        " 101   5.000000  1.000000e+02\n"
        " 102   4.500000  5.000000e+01\n"
        " -1 0    Stress Period 2\n"
    )
    ml = MfUsg(structured=False, model_ws=str(function_tmpdir))
    ModflowDis(ml, nlay=1, nrow=1, ncol=1, nper=2)
    ghb = MfUsgGhb.load(str(ghb_in), ml, nper=2, ext_unit_dict={})

    sp0 = ghb.stress_period_data[0]
    assert len(sp0) == 2
    assert list(sp0["node"]) == [100, 101]
    assert np.isclose(sp0["bhead"][0], 5.0, atol=0.01)
    assert np.isclose(sp0["cond"][1], 50.0, atol=0.1)

    # Reuse SP should copy SP0 data
    sp1 = ghb.stress_period_data[1]
    assert len(sp1) == 2

    ghb_out = function_tmpdir / "out.ghb"
    ghb.fn_path = str(ghb_out)
    ghb.write_file()
    text = ghb_out.read_text()
    assert "Stress Period 1" in text
    assert "Stress Period 2" in text
    assert "\n 101" in text

    ml2 = MfUsg(structured=False, model_ws=str(function_tmpdir))
    ModflowDis(ml2, nlay=1, nrow=1, ncol=1, nper=2)
    ghb2 = MfUsgGhb.load(str(ghb_out), ml2, nper=2, ext_unit_dict={})
    assert np.isclose(ghb2.stress_period_data[0]["bhead"][0], 5.0, atol=0.01)
    assert list(ghb2.stress_period_data[0]["node"]) == [100, 101]


def test_mfusgghb_aux_roundtrip(function_tmpdir):
    """GHB with one AUX concentration field round-trips correctly."""
    from flopy.modflow import ModflowDis

    ghb_in = function_tmpdir / "aux.ghb"
    ghb_in.write_text(
        "# MfUsgGhb aux test\n"
        " 1 0 AUX C01\n"
        " 1 0    Stress Period 1\n"
        " 101   6.000000  2.000000e+02  1.500000e-01\n"
    )
    ml = MfUsg(structured=False, model_ws=str(function_tmpdir))
    ModflowDis(ml, nlay=1, nrow=1, ncol=1, nper=1)
    ghb = MfUsgGhb.load(str(ghb_in), ml, nper=1, ext_unit_dict={})

    sp0 = ghb.stress_period_data[0]
    assert "C01" in ghb.dtype.names
    assert np.isclose(sp0["C01"][0], 0.15, atol=1e-5)

    ghb_out = function_tmpdir / "aux_out.ghb"
    ghb.fn_path = str(ghb_out)
    ghb.write_file()
    content = ghb_out.read_text()
    assert "AUX C01" in content

    ml2 = MfUsg(structured=False, model_ws=str(function_tmpdir))
    ModflowDis(ml2, nlay=1, nrow=1, ncol=1, nper=1)
    ghb2 = MfUsgGhb.load(str(ghb_out), ml2, nper=1, ext_unit_dict={})
    assert np.isclose(ghb2.stress_period_data[0]["C01"][0], 0.15, atol=1e-5)


# ---------------------------------------------------------------------------
# MfUsgDrn tests
# ---------------------------------------------------------------------------

def test_mfusgdrn_roundtrip(function_tmpdir):
    """Basic DRN: load → inspect → write → reload preserves data."""
    from flopy.modflow import ModflowDis

    drn_in = function_tmpdir / "test.drn"
    drn_in.write_text(
        "# MfUsgDrn test\n"
        "         3 0\n"
        " 3 0    Stress Period 1\n"
        " 101   2.500000  1.000000e+01\n"
        " 102   3.000000  2.000000e+01\n"
        " 103   1.800000  5.000000e+00\n"
        " -1 0    Stress Period 2\n"
    )
    ml = MfUsg(structured=False, model_ws=str(function_tmpdir))
    ModflowDis(ml, nlay=1, nrow=1, ncol=1, nper=2)
    drn = MfUsgDrn.load(str(drn_in), ml, nper=2, ext_unit_dict={})

    sp0 = drn.stress_period_data[0]
    assert len(sp0) == 3
    assert list(sp0["node"]) == [100, 101, 102]
    assert np.isclose(sp0["elev"][0], 2.5, atol=0.01)
    assert np.isclose(sp0["cond"][2], 5.0, atol=0.1)

    drn_out = function_tmpdir / "out.drn"
    drn.fn_path = str(drn_out)
    drn.write_file()
    text = drn_out.read_text()
    assert "Stress Period 1" in text
    assert "Stress Period 2" in text
    assert "\n 101" in text

    ml2 = MfUsg(structured=False, model_ws=str(function_tmpdir))
    ModflowDis(ml2, nlay=1, nrow=1, ncol=1, nper=2)
    drn2 = MfUsgDrn.load(str(drn_out), ml2, nper=2, ext_unit_dict={})
    assert np.isclose(drn2.stress_period_data[0]["elev"][0], 2.5, atol=0.01)
    assert list(drn2.stress_period_data[0]["node"]) == [100, 101, 102]


def test_mfusgdrn_aux_roundtrip(function_tmpdir):
    """DRN with AUX concentration field round-trips correctly."""
    from flopy.modflow import ModflowDis

    drn_in = function_tmpdir / "aux.drn"
    drn_in.write_text(
        "# MfUsgDrn aux test\n"
        " 1 0 AUX C01\n"
        " 1 0    Stress Period 1\n"
        " 201   1.200000  3.000000e+01  2.500000e-02\n"
    )
    ml = MfUsg(structured=False, model_ws=str(function_tmpdir))
    ModflowDis(ml, nlay=1, nrow=1, ncol=1, nper=1)
    drn = MfUsgDrn.load(str(drn_in), ml, nper=1, ext_unit_dict={})

    assert "C01" in drn.dtype.names
    assert np.isclose(drn.stress_period_data[0]["C01"][0], 0.025, atol=1e-6)

    drn_out = function_tmpdir / "aux_out.drn"
    drn.fn_path = str(drn_out)
    drn.write_file()
    assert "AUX C01" in drn_out.read_text()


# ---------------------------------------------------------------------------
# MfUsgTvm tests
# ---------------------------------------------------------------------------

def test_mfusgtvm_roundtrip(function_tmpdir):
    """TVM semantic load/write preserves global controls and SP boundaries."""
    from flopy.modflow import ModflowDis

    tvm_content = (
        "# MODFLOW-USG Time-Variant Materials (TVM) Package\n"
        "  1  -1  -1  0  0  -1  0\n"
        "  0  0  0  0  0  0    Stress Period 1\n"
        "  0  0  0  0  0  0    Stress Period 2\n"
        "  0  0  0  0  0  0    Stress Period 3\n"
    )
    tvm_in = function_tmpdir / "test.tvm"
    tvm_in.write_text(tvm_content)

    ml = MfUsg(structured=False, model_ws=str(function_tmpdir))
    ModflowDis(ml, nlay=1, nrow=1, ncol=1, nper=3)
    tvm = MfUsgTvm.load(str(tvm_in), ml, nper=3, ext_unit_dict={})

    assert tvm.itvmprint == 1
    assert tvm.tvmlogbasehk == -1.0
    assert tvm.tvmlogbasevka == -1.0
    assert sorted(tvm.stress_period_data) == [0, 1, 2]

    tvm_out = function_tmpdir / "out.tvm"
    tvm.fn_path = str(tvm_out)
    tvm.write_file()
    written = tvm_out.read_text()
    assert "Stress period number 1 start" in written
    assert "Stress period number 1 end" in written
    assert "Stress period number 2" in written
    assert "Stress period number 3" in written


def test_mfusgtvm_missing_sp_gets_zeros(function_tmpdir):
    """TVM SP blocks missing from the file are written as all-zero headers."""
    from flopy.modflow import ModflowDis

    # File only has 2 SPs but model has 3
    tvm_in = function_tmpdir / "short.tvm"
    tvm_in.write_text(
        "# TVM\n"
        "  0  0  0  0  0  0  0\n"
        "  0  0  0  0  0  0    Stress Period 1\n"
        "  0  0  0  0  0  0    Stress Period 2\n"
    )
    ml = MfUsg(structured=False, model_ws=str(function_tmpdir))
    ModflowDis(ml, nlay=1, nrow=1, ncol=1, nper=3)
    tvm = MfUsgTvm.load(str(tvm_in), ml, nper=3, ext_unit_dict={})
    tvm_out = function_tmpdir / "out.tvm"
    tvm.fn_path = str(tvm_out)
    tvm.write_file()
    written = tvm_out.read_text()
    # SP 3 must be emitted with zeros even though it wasn't in the source
    assert "Stress period number 3" in written


# ---------------------------------------------------------------------------
# MfUsgGsf tests
# ---------------------------------------------------------------------------

# Minimal valid GSF: 1-layer, 1 triangular cell, 3 vertices.
# Format follows UnstructuredGrid.from_gridspec():
#   UNSTRUCTURED
#   NNODES [...]
#   NVERTS
#   X Y Z   (for each vertex)
#   NODENO XC YC ZC LAY NVERTS V1 V2 ...  (for each node, 1-based vertex indices)
_MINIMAL_GSF_LINES = [
    "UNSTRUCTURED\n",
    "1 1 1 1\n",       # nnodes=1 (only first token used)
    "3\n",             # 3 vertices
    "0.0 0.0 10.0\n",  # vertex 1
    "1.0 0.0 10.0\n",  # vertex 2
    "0.5 1.0 10.0\n",  # vertex 3
    "1 0.5 0.333 5.0 1 3 1 2 3\n",  # node 1: xc=0.5, yc=0.333, lay=1, 3 verts (1-based)
]


def test_mfusggsf_load_stores_lines(function_tmpdir):
    """MfUsgGsf.load() stores all file lines verbatim."""
    from flopy.modflow import ModflowDis

    gsf_in = function_tmpdir / "test.gsf"
    gsf_in.write_text("".join(_MINIMAL_GSF_LINES))

    ml = MfUsg(structured=False, model_ws=str(function_tmpdir))
    ModflowDis(ml, nlay=1, nrow=1, ncol=1, nper=1)
    gsf = MfUsgGsf.load(str(gsf_in), ml, ext_unit_dict={})

    assert len(gsf.lines) == len(_MINIMAL_GSF_LINES)
    assert "UNSTRUCTURED" in gsf.lines[0].upper()


def test_mfusggsf_text_roundtrip(function_tmpdir):
    """write_file() reproduces the exact content that was loaded."""
    from flopy.modflow import ModflowDis

    gsf_in = function_tmpdir / "in.gsf"
    original_text = "".join(_MINIMAL_GSF_LINES)
    gsf_in.write_text(original_text)

    ml = MfUsg(structured=False, model_ws=str(function_tmpdir))
    ModflowDis(ml, nlay=1, nrow=1, ncol=1, nper=1)
    gsf = MfUsgGsf.load(str(gsf_in), ml, ext_unit_dict={})

    gsf_out = function_tmpdir / "out.gsf"
    gsf.fn_path = str(gsf_out)
    gsf.write_file()

    # Every line of the original must appear in the written file
    written_lines = gsf_out.read_text().splitlines()
    for ln in original_text.splitlines():
        assert ln in written_lines, f"Line missing from written GSF: {ln!r}"


def test_mfusggsf_to_grid(function_tmpdir):
    """to_grid() returns an UnstructuredGrid parsed from the synthetic GSF."""
    from flopy.discretization import UnstructuredGrid
    from flopy.modflow import ModflowDis

    gsf_in = function_tmpdir / "grid.gsf"
    gsf_in.write_text("".join(_MINIMAL_GSF_LINES))

    ml = MfUsg(structured=False, model_ws=str(function_tmpdir))
    ModflowDis(ml, nlay=1, nrow=1, ncol=1, nper=1)
    gsf = MfUsgGsf.load(str(gsf_in), ml, ext_unit_dict={})

    # Point fn_path to the input file so to_grid() reads it without needing write_file first
    gsf.fn_path = str(gsf_in)

    grid = gsf.to_grid()
    assert isinstance(grid, UnstructuredGrid)
    assert len(grid.xcellcenters) == 1          # 1 node
    assert np.isclose(grid.xcellcenters[0], 0.5, atol=1e-6)
    assert np.isclose(grid.ycellcenters[0], 0.333, atol=1e-3)


# Semantic equivalent of _MINIMAL_GSF_LINES (0-based vertex/node references).
_MINIMAL_GSF_VERTICES = [(0.0, 0.0, 10.0), (1.0, 0.0, 10.0), (0.5, 1.0, 10.0)]
_MINIMAL_GSF_NODES = [
    {"node": 0, "xc": 0.5, "yc": 0.333, "zc": 5.0, "layer": 0, "vertices": [0, 1, 2]}
]


def test_mfusggsf_semantic_load(function_tmpdir):
    """load(parse=True) parses a minimal GSF into 0-based semantic data."""
    from flopy.modflow import ModflowDis

    gsf_in = function_tmpdir / "in.gsf"
    gsf_in.write_text("".join(_MINIMAL_GSF_LINES))

    ml = MfUsg(structured=False, model_ws=str(function_tmpdir))
    ModflowDis(ml, nlay=1, nrow=1, ncol=1, nper=1)
    gsf = MfUsgGsf.load(str(gsf_in), ml, ext_unit_dict={}, parse=True)

    assert gsf.lines is None  # parsed, not raw
    assert gsf.nnodes == 1 and gsf.nlay == 1
    # file vertex ids 1,2,3 -> internal 0,1,2
    assert gsf.node_data[0]["vertices"] == [0, 1, 2]
    # file node 1 -> internal 0; file layer 1 -> internal 0
    assert gsf.node_data[0]["node"] == 0
    assert gsf.node_data[0]["layer"] == 0
    assert np.allclose(gsf.vertices, _MINIMAL_GSF_VERTICES)


def test_mfusggsf_authoring_from_scratch(function_tmpdir):
    """Build a GSF from Python data (no file loaded); file ids are 1-based."""
    from flopy.discretization import UnstructuredGrid
    from flopy.modflow import ModflowDis

    ml = MfUsg(structured=False, model_ws=str(function_tmpdir))
    ModflowDis(ml, nlay=1, nrow=1, ncol=1, nper=1)
    gsf = MfUsgGsf(
        ml, vertices=_MINIMAL_GSF_VERTICES, node_data=_MINIMAL_GSF_NODES
    )
    gsf.fn_path = str(function_tmpdir / "scratch.gsf")
    gsf.write_file()

    lines = [
        ln for ln in Path(gsf.fn_path).read_text().splitlines() if ln.strip()
    ]
    assert lines[0].upper().startswith("UNSTRUCTURED")
    assert lines[1].split()[:2] == ["1", "1"]  # nnodes nlay
    assert lines[2].split()[0] == "3"  # nverts
    node_line = lines[-1].split()
    assert node_line[0] == "1"  # internal node 0 -> file 1
    assert node_line[4] == "1"  # internal layer 0 -> file 1
    assert node_line[5] == "3"  # nvert
    assert node_line[6:9] == ["1", "2", "3"]  # internal verts 0,1,2 -> file 1,2,3

    # The authored file parses back to a real UnstructuredGrid.
    grid = gsf.to_grid()
    assert isinstance(grid, UnstructuredGrid)
    assert len(grid.xcellcenters) == 1
    assert np.isclose(grid.xcellcenters[0], 0.5, atol=1e-6)


def test_mfusggsf_write_reload_roundtrip(function_tmpdir):
    """Author -> write -> load(parse=True) preserves semantics and 0-based ids."""
    from flopy.modflow import ModflowDis

    def ml(name):
        m = MfUsg(structured=False, model_ws=str(function_tmpdir), modelname=name)
        ModflowDis(m, nlay=1, nrow=1, ncol=1, nper=1)
        return m

    gsf = MfUsgGsf(
        ml("w"), vertices=_MINIMAL_GSF_VERTICES, node_data=_MINIMAL_GSF_NODES
    )
    gsf.fn_path = str(function_tmpdir / "rt.gsf")
    gsf.write_file()

    re = MfUsgGsf.load(gsf.fn_path, ml("w2"), parse=True)
    assert np.allclose(re.vertices, _MINIMAL_GSF_VERTICES)
    assert re.node_data[0]["vertices"] == [0, 1, 2]
    assert re.node_data[0]["layer"] == 0

    # Second write is byte-identical (stable round-trip).
    first = Path(gsf.fn_path).read_text()
    re.fn_path = str(function_tmpdir / "rt_again.gsf")
    re.write_file()
    assert Path(re.fn_path).read_text() == first


def test_mfusggsf_rejects_invalid_and_mixed_modes(function_tmpdir):
    """Invalid vertex refs and ambiguous input modes fail explicitly."""
    from flopy.modflow import ModflowDis

    def ml():
        m = MfUsg(structured=False, model_ws=str(function_tmpdir))
        ModflowDis(m, nlay=1, nrow=1, ncol=1, nper=1)
        return m

    # vertex reference out of range
    with pytest.raises(ValueError, match="out of range"):
        MfUsgGsf(
            ml(),
            vertices=_MINIMAL_GSF_VERTICES,
            node_data=[{"xc": 0.5, "yc": 0.3, "layer": 0, "vertices": [0, 1, 9]}],
        )

    # raw lines + semantic data together (no silent precedence)
    with pytest.raises(ValueError, match="not both"):
        MfUsgGsf(
            ml(),
            lines=["UNSTRUCTURED\n"],
            vertices=_MINIMAL_GSF_VERTICES,
            node_data=_MINIMAL_GSF_NODES,
        )

    # semantic mode needs both vertices and node_data
    with pytest.raises(ValueError, match="requires both"):
        MfUsgGsf(ml(), vertices=_MINIMAL_GSF_VERTICES)


def test_mfusggsf_from_grid(function_tmpdir):
    """from_grid builds a GSF from an UnstructuredGrid (requires per-vertex z)."""
    from flopy.discretization import UnstructuredGrid
    from flopy.modflow import ModflowDis

    def ml():
        m = MfUsg(structured=False, model_ws=str(function_tmpdir))
        ModflowDis(m, nlay=1, nrow=1, ncol=1, nper=1)
        return m

    gsf_in = function_tmpdir / "src.gsf"
    gsf_in.write_text("".join(_MINIMAL_GSF_LINES))
    grid = UnstructuredGrid.from_gridspec(str(gsf_in), split_vertices=True)

    # UnstructuredGrid drops per-vertex z, so from_grid must demand it.
    with pytest.raises(ValueError, match="zverts"):
        MfUsgGsf.from_grid(ml(), grid)

    gsf = MfUsgGsf.from_grid(ml(), grid, zverts=[10.0, 10.0, 10.0])
    assert gsf.nnodes == 1 and gsf.nlay == 1
    assert gsf.node_data[0]["vertices"] == [0, 1, 2]  # grid iverts are 0-based
    gsf.fn_path = str(function_tmpdir / "from_grid.gsf")
    gsf.write_file()
    assert isinstance(gsf.to_grid(), UnstructuredGrid)


def test_mfusggsf_unstructured_gwf_header(function_tmpdir):
    """An 'UNSTRUCTURED GWF' header loads, writes, and round-trips to_grid().

    Regression: from_gridspec previously rejected the valid two-token header
    'UNSTRUCTURED GWF' due to an operator-precedence bug.
    """
    from flopy.discretization import UnstructuredGrid
    from flopy.modflow import ModflowDis

    text = "UNSTRUCTURED GWF\n" + "".join(_MINIMAL_GSF_LINES[1:])
    src = function_tmpdir / "gwf.gsf"
    src.write_text(text)

    # The grid parser itself must accept the two-token header.
    assert isinstance(
        UnstructuredGrid.from_gridspec(str(src)), UnstructuredGrid
    )

    ml = MfUsg(structured=False, model_ws=str(function_tmpdir))
    ModflowDis(ml, nlay=1, nrow=1, ncol=1, nper=1)
    gsf = MfUsgGsf.load(str(src), ml, parse=True)
    assert gsf.header == "UNSTRUCTURED GWF"
    gsf.fn_path = str(function_tmpdir / "gwf_out.gsf")
    gsf.write_file()
    assert "UNSTRUCTURED GWF" in Path(gsf.fn_path).read_text()
    assert isinstance(gsf.to_grid(), UnstructuredGrid)


def test_mfusggsf_parse_rejects_trailing_content(function_tmpdir):
    """parse=True with extra non-comment content falls back to the raw lines."""
    from flopy.modflow import ModflowDis

    src = function_tmpdir / "trailing.gsf"
    src.write_text("".join(_MINIMAL_GSF_LINES) + "EXTRA UNEXPECTED LINE\n")

    ml = MfUsg(structured=False, model_ws=str(function_tmpdir))
    ModflowDis(ml, nlay=1, nrow=1, ncol=1, nper=1)
    gsf = MfUsgGsf.load(str(src), ml, parse=True)

    assert gsf.node_data is None  # not parsed semantically
    assert gsf.lines is not None  # raw fallback
    gsf.fn_path = str(function_tmpdir / "trailing_out.gsf")
    gsf.write_file()
    assert "EXTRA UNEXPECTED LINE" in Path(gsf.fn_path).read_text()


def test_mfusggsf_from_grid_top_bottom(function_tmpdir):
    """from_grid with top/bottom elevations writes the USG-T doubled-vertex GSF."""
    from flopy.discretization import UnstructuredGrid
    from flopy.modflow import ModflowDis

    src = function_tmpdir / "tri.gsf"
    src.write_text("".join(_MINIMAL_GSF_LINES))
    grid = UnstructuredGrid.from_gridspec(str(src))
    nverts = grid.verts.shape[0]  # 3

    ml = MfUsg(structured=False, model_ws=str(function_tmpdir))
    ModflowDis(ml, nlay=1, nrow=1, ncol=1, nper=1)
    gsf = MfUsgGsf.from_grid(
        ml, grid, top_zverts=[10.0, 10.0, 10.0], bot_zverts=[0.0, 0.0, 0.0]
    )
    gsf.fn_path = str(function_tmpdir / "tb.gsf")
    gsf.write_file()

    lines = [ln for ln in Path(gsf.fn_path).read_text().splitlines() if ln.strip()]
    assert lines[2].split()[0] == str(2 * nverts)  # NVERTS doubled
    node = lines[-1].split()
    assert node[5] == str(2 * 3)  # 2 * ncell_vertices
    assert node[6:9] == ["1", "2", "3"]  # top ids, 1-based
    # bottom ids = top ids + totalverts (nverts)
    assert node[9:12] == [str(1 + nverts), str(2 + nverts), str(3 + nverts)]

    # split_vertices=True reconstructs the correct top/botm
    grid2 = gsf.to_grid()
    assert isinstance(grid2, UnstructuredGrid)
    assert np.allclose(grid2.top, 10.0)
    assert np.allclose(grid2.botm, 0.0)

    # top_zverts without bot_zverts is rejected
    with pytest.raises(ValueError, match="together"):
        MfUsgGsf.from_grid(ml, grid, top_zverts=[1.0, 1.0, 1.0])


def test_mfusggsf_from_disv_gridprops(function_tmpdir):
    """from_disv_gridprops builds a single-layer GSF (notebook workflow)."""
    from flopy.discretization import UnstructuredGrid
    from flopy.modflow import ModflowDis

    # 6 vertices; cell 0 is a square with a closing duplicate vertex.
    disv = {
        "vertices": [
            (0, 0.0, 0.0),
            (1, 1.0, 0.0),
            (2, 1.0, 1.0),
            (3, 0.0, 1.0),
            (4, 2.0, 0.0),
            (5, 2.0, 1.0),
        ],
        "cell2d": [
            [0, 0.5, 0.5, 5, 0, 1, 2, 3, 0],  # closing dup -> 4 unique verts
            [1, 1.5, 0.5, 3, 1, 4, 5],  # triangle, 3 unique verts
        ],
    }
    def ml():
        m = MfUsg(structured=False, model_ws=str(function_tmpdir))
        ModflowDis(m, nlay=1, nrow=1, ncol=1, nper=1)
        return m

    gsf = MfUsgGsf.from_disv_gridprops(ml(), disv, top=10.0, botm=0.0)
    assert gsf.nlay == 1
    assert len(gsf.vertices) == 12  # 2 * 6 (top + bottom)
    # square cell: closing dup dropped, bottom ids shifted by nvert (6)
    assert gsf.node_data[0]["vertices"] == [0, 1, 2, 3, 6, 7, 8, 9]
    gsf.fn_path = str(function_tmpdir / "disv.gsf")
    gsf.write_file()
    assert isinstance(gsf.to_grid(), UnstructuredGrid)

    # a degenerate cell (< 3 unique vertices) is rejected by default,
    disv_deg = {
        "vertices": disv["vertices"],
        "cell2d": [[0, 0.5, 0.5, 5, 0, 1, 2, 3, 0], [1, 1.5, 0.5, 2, 4, 5]],
    }
    with pytest.raises(ValueError, match="unique vertices"):
        MfUsgGsf.from_disv_gridprops(ml(), disv_deg, top=10.0, botm=0.0)
    # ... or skipped (with consecutive node renumbering) on request
    gsf2 = MfUsgGsf.from_disv_gridprops(
        ml(), disv_deg, top=10.0, botm=0.0, skip_degenerate=True
    )
    assert gsf2.nnodes == 1
    assert gsf2.node_data[0]["node"] == 0


def test_mfusggsf_vertex_modes(function_tmpdir):
    """shared vs cell vertex modes (the two GRIDGEN2GSF layouts).

    Two adjacent quad cells share an edge (vertices 1 and 2). The parsimonious
    'shared' mode reuses those vertex ids across the neighbouring cells; the
    non-parsimonious 'cell' mode gives every cell its own unique vertices
    (8 per quad), with the top half then the bottom half so that
    from_gridspec(split_vertices=True) recovers top/botm.
    """
    from flopy.discretization import UnstructuredGrid
    from flopy.modflow import ModflowDis

    def ml():
        m = MfUsg(structured=False, model_ws=str(function_tmpdir))
        ModflowDis(m, nlay=1, nrow=1, ncol=1, nper=1)
        return m

    disv = {
        "vertices": [
            (0, 0.0, 0.0),
            (1, 1.0, 0.0),
            (2, 1.0, 1.0),
            (3, 0.0, 1.0),
            (4, 2.0, 0.0),
            (5, 2.0, 1.0),
        ],
        "cell2d": [
            [0, 0.5, 0.5, 4, 0, 1, 2, 3],
            [1, 1.5, 0.5, 4, 1, 4, 5, 2],  # shares verts 1, 2 with cell 0
        ],
    }

    # parsimonious / shared: neighbouring cells reuse vertex ids
    shared = MfUsgGsf.from_disv_gridprops(
        ml(), disv, top=10.0, botm=0.0, vertex_mode="parsimonious"
    )
    assert len(shared.vertices) == 12  # 2 * 6 shared
    common = set(shared.node_data[0]["vertices"]) & set(
        shared.node_data[1]["vertices"]
    )
    assert common  # ids are shared between neighbours

    # non-parsimonious / cell: no shared ids, 8 unique vertices per quad
    cell = MfUsgGsf.from_disv_gridprops(
        ml(), disv, top=10.0, botm=0.0, vertex_mode="nonparsimonious"
    )
    assert len(cell.vertices) == 2 * 8  # 8 unique per quad cell
    assert not (
        set(cell.node_data[0]["vertices"]) & set(cell.node_data[1]["vertices"])
    )
    v0 = cell.node_data[0]["vertices"]
    assert len(v0) == 8  # quad -> 8 vertices
    # top half then bottom half (split_vertices convention)
    assert v0[:4] == [0, 1, 2, 3] and v0[4:] == [4, 5, 6, 7]

    # both modes reconstruct correct top/botm through from_gridspec
    for gsf, name in ((shared, "shared"), (cell, "cell")):
        gsf.fn_path = str(function_tmpdir / f"{name}.gsf")
        gsf.write_file()
        grid = gsf.to_grid()
        assert isinstance(grid, UnstructuredGrid)
        assert np.allclose(grid.top, 10.0)
        assert np.allclose(grid.botm, 0.0)

    # cell mode requires real top/bottom elevations (not a single surface)
    src = function_tmpdir / "tri.gsf"
    src.write_text("".join(_MINIMAL_GSF_LINES))
    tgrid = UnstructuredGrid.from_gridspec(str(src))
    nv = tgrid.verts.shape[0]
    with pytest.raises(ValueError, match="vertex_mode='cell'"):
        MfUsgGsf.from_grid(ml(), tgrid, zverts=[1.0] * nv, vertex_mode="cell")


def test_mfusggsf_hardening_rejects(function_tmpdir):
    """Semantic constructor rejects bad header, duplicate nodes, and small nlay."""
    from flopy.modflow import ModflowDis

    def ml():
        m = MfUsg(structured=False, model_ws=str(function_tmpdir))
        ModflowDis(m, nlay=1, nrow=1, ncol=1, nper=1)
        return m

    v = _MINIMAL_GSF_VERTICES
    nd = [{"node": 0, "xc": 0.5, "yc": 0.3, "layer": 0, "vertices": [0, 1, 2]}]

    with pytest.raises(ValueError, match="header"):
        MfUsgGsf(ml(), vertices=v, node_data=nd, header="STRUCTURED")
    with pytest.raises(ValueError, match="contiguous and ordered"):
        MfUsgGsf(ml(), vertices=v, node_data=[nd[0], dict(nd[0])])
    with pytest.raises(ValueError, match="nlay"):
        MfUsgGsf(
            ml(),
            vertices=v,
            node_data=[{"xc": 0.5, "yc": 0.3, "layer": 2, "vertices": [0, 1, 2]}],
            nlay=1,
        )


# Minimal GSF body shared by the strict-parse tests (vertices + one node).
_GSF_BODY = "3\n0 0 10\n1 0 10\n0.5 1 10\n1 0.5 0.333 5 1 3 1 2 3\n"


def _gsf_parse_mode(function_tmpdir, text, name):
    """Load `text` with parse=True; return 'semantic' or 'raw'."""
    from flopy.modflow import ModflowDis

    p = function_tmpdir / f"{name}.gsf"
    p.write_text(text)
    m = MfUsg(structured=False, model_ws=str(function_tmpdir), modelname=name)
    ModflowDis(m, nlay=1, nrow=1, ncol=3, nper=1)
    gsf = MfUsgGsf.load(str(p), m, parse=True)
    return "raw" if gsf.node_data is None else "semantic"


def test_mfusggsf_parse_header_strict(function_tmpdir):
    """parse=True accepts only UNSTRUCTURED / UNSTRUCTURED GWF; else raw fallback."""
    assert (
        _gsf_parse_mode(function_tmpdir, "UNSTRUCTURED\n1 1 1 1\n" + _GSF_BODY, "a")
        == "semantic"
    )
    assert (
        _gsf_parse_mode(
            function_tmpdir, "UNSTRUCTURED GWF\n1 1 1 1\n" + _GSF_BODY, "b"
        )
        == "semantic"
    )
    # extra token between UNSTRUCTURED and GWF is rejected -> raw round-trip
    assert (
        _gsf_parse_mode(
            function_tmpdir, "UNSTRUCTURED EXTRA GWF\n1 1 1 1\n" + _GSF_BODY, "c"
        )
        == "raw"
    )


def test_mfusggsf_iz_ic_flags(function_tmpdir):
    """IZ/IC line-2 flags: authoring requires (1,1); parse allows omitted or 1 1."""
    from flopy.modflow import ModflowDis

    def ml():
        m = MfUsg(structured=False, model_ws=str(function_tmpdir))
        ModflowDis(m, nlay=1, nrow=1, ncol=3, nper=1)
        return m

    v = _MINIMAL_GSF_VERTICES
    nd = [{"node": 0, "xc": 0.5, "yc": 0.3, "layer": 0, "vertices": [0, 1, 2]}]

    # authoring: omitted (assumed) and explicit (1,1) are accepted
    assert MfUsgGsf(ml(), vertices=v, node_data=nd).extra_header == (1, 1)
    assert (
        MfUsgGsf(ml(), vertices=v, node_data=nd, extra_header=(1, 1)).extra_header
        == (1, 1)
    )
    # authoring: anything other than (1, 1), or a bad length, is rejected
    for bad in ((0, 1), (1, 0), (1, 1, 9)):
        with pytest.raises(ValueError, match="IZ IC"):
            MfUsgGsf(ml(), vertices=v, node_data=nd, extra_header=bad)

    # parse: 'nnode nlay' (omitted -> assumed) and 'nnode nlay 1 1' are semantic
    assert (
        _gsf_parse_mode(function_tmpdir, "UNSTRUCTURED\n1 1\n" + _GSF_BODY, "p2")
        == "semantic"
    )
    assert (
        _gsf_parse_mode(function_tmpdir, "UNSTRUCTURED\n1 1 1 1\n" + _GSF_BODY, "p4")
        == "semantic"
    )
    # parse: 0 1, 1 0, or an odd length -> raw fallback
    for flags, name in (("0 1", "f01"), ("1 0", "f10"), ("1", "f3"), ("1 1 9", "f5")):
        text = f"UNSTRUCTURED\n1 1 {flags}\n" + _GSF_BODY
        assert _gsf_parse_mode(function_tmpdir, text, name) == "raw"


def test_mfusggsf_inode_validation(function_tmpdir):
    """Node ids must be 0..nnodes-1 ordered (authoring error / parse raw fallback)."""
    from flopy.modflow import ModflowDis

    def ml():
        m = MfUsg(structured=False, model_ws=str(function_tmpdir))
        ModflowDis(m, nlay=1, nrow=1, ncol=3, nper=1)
        return m

    v = _MINIMAL_GSF_VERTICES

    def rec(node):
        return {"node": node, "xc": 0.5, "yc": 0.3, "layer": 0, "vertices": [0, 1, 2]}

    # authoring: gap, duplicate, and reorder all fail explicitly
    for ids in ([0, 2], [0, 0], [1, 0]):
        with pytest.raises(ValueError, match="contiguous and ordered"):
            MfUsgGsf(ml(), vertices=v, node_data=[rec(i) for i in ids])

    # parse: node numbers out of order (2 then 1) -> raw fallback
    body = (
        "3\n0 0 10\n1 0 10\n0.5 1 10\n"
        "2 0.5 0.3 5 1 3 1 2 3\n"
        "1 0.6 0.3 5 1 3 1 2 3\n"
    )
    assert (
        _gsf_parse_mode(function_tmpdir, "UNSTRUCTURED\n2 1 1 1\n" + body, "ord")
        == "raw"
    )


# ---------------------------------------------------------------------------
# gridgen_to_gsf utility tests (flopy/mfusg/gridgen2gsf.py)
# ---------------------------------------------------------------------------

# Two adjacent quad cells sharing the edge vertices 1 and 2.
_GG_DISV = {
    "vertices": [
        (0, 0.0, 0.0),
        (1, 1.0, 0.0),
        (2, 1.0, 1.0),
        (3, 0.0, 1.0),
        (4, 2.0, 0.0),
        (5, 2.0, 1.0),
    ],
    "cell2d": [
        [0, 0.5, 0.5, 4, 0, 1, 2, 3],
        [1, 1.5, 0.5, 4, 1, 4, 5, 2],
    ],
}


def _gg_model(function_tmpdir, name="gg"):
    from flopy.modflow import ModflowDis

    m = MfUsg(structured=False, model_ws=str(function_tmpdir), modelname=name)
    ModflowDis(m, nlay=1, nrow=1, ncol=1, nper=1)
    return m


def test_gridgen_to_gsf_disv_modes(function_tmpdir):
    """gridgen_to_gsf from disv_gridprops: shared shares ids, cell does not."""
    from flopy.discretization import UnstructuredGrid
    from flopy.mfusg import MfUsgGsf, gridgen_to_gsf

    shared = gridgen_to_gsf(
        _gg_model(function_tmpdir, "s"),
        _GG_DISV,
        top=10.0,
        botm=0.0,
        vertex_mode="parsimonious",
    )
    assert isinstance(shared, MfUsgGsf)
    assert len(shared.vertices) == 12  # 2 * 6 shared
    assert set(shared.node_data[0]["vertices"]) & set(
        shared.node_data[1]["vertices"]
    )  # neighbours reuse ids

    cell = gridgen_to_gsf(
        _gg_model(function_tmpdir, "c"),
        _GG_DISV,
        top=10.0,
        botm=0.0,
        vertex_mode="nonparsimonious",
    )
    assert len(cell.vertices) == 2 * 8  # unique per quad
    assert not (
        set(cell.node_data[0]["vertices"]) & set(cell.node_data[1]["vertices"])
    )
    assert len(cell.node_data[0]["vertices"]) == 8  # quad -> 8 vertices

    # both reconstruct top/botm through to_grid (split_vertices)
    for gsf, name in ((shared, "s"), (cell, "c")):
        gsf.fn_path = str(function_tmpdir / f"gg_{name}.gsf")
        gsf.write_file()
        grid = gsf.to_grid()
        assert isinstance(grid, UnstructuredGrid)
        assert np.allclose(grid.top, 10.0) and np.allclose(grid.botm, 0.0)


def test_gridgen_to_gsf_source_types(function_tmpdir):
    """gridgen_to_gsf accepts a Gridgen-like object and an UnstructuredGrid."""
    from flopy.discretization import UnstructuredGrid
    from flopy.mfusg import gridgen_to_gsf

    # flopy Gridgen-like object: only get_gridprops_disv() is required
    class _FakeGridgen:
        def get_gridprops_disv(self):
            return _GG_DISV

    g = gridgen_to_gsf(_gg_model(function_tmpdir, "fg"), _FakeGridgen(), top=5.0)
    assert g.nnodes == 2

    # UnstructuredGrid source (top/botm broadcast to per-vertex surfaces)
    src = function_tmpdir / "tri.gsf"
    src.write_text("".join(_MINIMAL_GSF_LINES))
    ug = UnstructuredGrid.from_gridspec(str(src))
    gu = gridgen_to_gsf(_gg_model(function_tmpdir, "ug"), ug, top=20.0, botm=2.0)
    gu.fn_path = str(function_tmpdir / "gg_ug.gsf")
    gu.write_file()
    grid = gu.to_grid()
    assert np.allclose(grid.top, 20.0) and np.allclose(grid.botm, 2.0)


def test_gridgen_to_gsf_validation(function_tmpdir):
    """gridgen_to_gsf fails explicitly on bad mode/source/geometry."""
    from flopy.mfusg import gridgen_to_gsf

    with pytest.raises(ValueError, match="vertex_mode"):
        gridgen_to_gsf(_gg_model(function_tmpdir), _GG_DISV, vertex_mode="bogus")

    with pytest.raises(TypeError, match="source"):
        gridgen_to_gsf(_gg_model(function_tmpdir), 12345)

    with pytest.raises(ValueError, match="vertices.*cell2d|cell2d"):
        gridgen_to_gsf(_gg_model(function_tmpdir), {"vertices": []})

    # degenerate cell: raises by default, dropped + renumbered with the flag
    disv_deg = {
        "vertices": _GG_DISV["vertices"],
        "cell2d": [[0, 0.5, 0.5, 4, 0, 1, 2, 3], [1, 1.5, 0.5, 2, 4, 5]],
    }
    with pytest.raises(ValueError, match="unique vertices"):
        gridgen_to_gsf(_gg_model(function_tmpdir), disv_deg, top=10.0, botm=0.0)
    g = gridgen_to_gsf(
        _gg_model(function_tmpdir), disv_deg, top=10.0, botm=0.0, skip_degenerate=True
    )
    assert g.nnodes == 1


# ---------------------------------------------------------------------------
# MfUsgSgb tests (Specified Gradient Boundary, glo2sgbu1.f)
# ---------------------------------------------------------------------------

def test_mfusgsgb_roundtrip(function_tmpdir):
    """SGB: load → inspect → write → reload preserves data and -1 reuse."""
    from flopy.modflow import ModflowDis

    sgb_in = function_tmpdir / "test.sgb"
    sgb_in.write_text(
        "# MfUsgSgb test\n"
        "         2 0\n"
        " 2 0    Stress Period 1\n"
        " 101   1.000000e-02\n"
        " 102   2.500000e-02\n"
        " -1 0    Stress Period 2\n"
    )
    ml = MfUsg(structured=False, model_ws=str(function_tmpdir))
    ModflowDis(ml, nlay=1, nrow=1, ncol=1, nper=2)
    sgb = MfUsgSgb.load(str(sgb_in), ml, nper=2, ext_unit_dict={})

    sp0 = sgb.stress_period_data[0]
    assert len(sp0) == 2
    assert list(sp0["node"]) == [100, 101]
    assert np.isclose(sp0["gradient"][0], 0.01, atol=1e-6)
    assert np.isclose(sp0["gradient"][1], 0.025, atol=1e-6)

    # -1 reuse copies SP0 data into SP1
    sp1 = sgb.stress_period_data[1]
    assert len(sp1) == 2
    assert list(sp1["node"]) == [100, 101]

    sgb_out = function_tmpdir / "out.sgb"
    sgb.fn_path = str(sgb_out)
    sgb.write_file()
    text = sgb_out.read_text()
    assert "Stress Period 1" in text
    assert "Stress Period 2" in text
    assert "\n 101" in text

    ml2 = MfUsg(structured=False, model_ws=str(function_tmpdir))
    ModflowDis(ml2, nlay=1, nrow=1, ncol=1, nper=2)
    sgb2 = MfUsgSgb.load(str(sgb_out), ml2, nper=2, ext_unit_dict={})
    assert list(sgb2.stress_period_data[0]["node"]) == [100, 101]
    assert np.isclose(sgb2.stress_period_data[0]["gradient"][1], 0.025, atol=1e-6)


def test_mfusgsgb_aux_roundtrip(function_tmpdir):
    """SGB with one AUX concentration field round-trips correctly."""
    from flopy.modflow import ModflowDis

    sgb_in = function_tmpdir / "aux.sgb"
    sgb_in.write_text(
        "# MfUsgSgb aux test\n"
        " 1 0 AUX C01\n"
        " 1 0    Stress Period 1\n"
        " 101   3.000000e-02  1.500000e-01\n"
    )
    ml = MfUsg(structured=False, model_ws=str(function_tmpdir))
    ModflowDis(ml, nlay=1, nrow=1, ncol=1, nper=1)
    sgb = MfUsgSgb.load(str(sgb_in), ml, nper=1, ext_unit_dict={})

    assert "C01" in sgb.dtype.names
    assert np.isclose(sgb.stress_period_data[0]["C01"][0], 0.15, atol=1e-5)

    sgb_out = function_tmpdir / "aux_out.sgb"
    sgb.fn_path = str(sgb_out)
    sgb.write_file()
    assert "AUX C01" in sgb_out.read_text()

    ml2 = MfUsg(structured=False, model_ws=str(function_tmpdir))
    ModflowDis(ml2, nlay=1, nrow=1, ncol=1, nper=1)
    sgb2 = MfUsgSgb.load(str(sgb_out), ml2, nper=1, ext_unit_dict={})
    assert np.isclose(sgb2.stress_period_data[0]["C01"][0], 0.15, atol=1e-5)


def test_mfusgsgb_programmatic_authoring(function_tmpdir):
    """SGB authored from scratch uses 0-based nodes; file is written 1-based."""
    from flopy.modflow import ModflowDis

    ml = MfUsg(structured=False, model_ws=str(function_tmpdir))
    ModflowDis(ml, nlay=1, nrow=1, ncol=1, nper=1)

    dtype = MfUsgSgb.get_default_dtype()
    spd = {0: np.array([(0, 0.01), (4, 0.05)], dtype=dtype).view(np.recarray)}
    sgb = MfUsgSgb(ml, stress_period_data=spd)
    sgb.fn_path = str(function_tmpdir / "created.sgb")
    sgb.write_file()
    text = Path(sgb.fn_path).read_text()
    assert "\n 1  1.000000e-02" in text
    assert "\n 5  5.000000e-02" in text

    ml2 = MfUsg(structured=False, model_ws=str(function_tmpdir))
    ModflowDis(ml2, nlay=1, nrow=1, ncol=1, nper=1)
    sgb2 = MfUsgSgb.load(sgb.fn_path, ml2, nper=1)
    assert list(sgb2.stress_period_data[0]["node"]) == [0, 4]
    assert np.isclose(sgb2.stress_period_data[0]["gradient"][1], 0.05, atol=1e-6)


def test_mfusgsgb_nam_registry():
    """MfUsg.load() registry maps the SGB package key to MfUsgSgb."""
    ml = MfUsg(structured=False)
    assert ml.mfnam_packages["sgb"] is MfUsgSgb


def test_mfusgsgb_parameters_fail_explicitly(function_tmpdir):
    """Named SGB parameters (NPSGB>0 / per-SP NP>0) raise NotImplementedError."""
    from flopy.modflow import ModflowDis

    # PARAMETER prefix in the header (NPSGB > 0)
    sgb_param = function_tmpdir / "param.sgb"
    sgb_param.write_text(
        "# param header\n"
        " PARAMETER 1 5 2 0\n"
        " 0 0    Stress Period 1\n"
    )
    ml = MfUsg(structured=False, model_ws=str(function_tmpdir))
    ModflowDis(ml, nlay=1, nrow=1, ncol=1, nper=1)
    with pytest.raises(NotImplementedError):
        MfUsgSgb.load(str(sgb_param), ml, nper=1, ext_unit_dict={})

    # active parameter declared in a stress period (NP > 0)
    sgb_sp = function_tmpdir / "sp_param.sgb"
    sgb_sp.write_text(
        "# sp param\n"
        "         2 0\n"
        " 1 1    Stress Period 1\n"
        " 101   1.000000e-02\n"
    )
    ml2 = MfUsg(structured=False, model_ws=str(function_tmpdir))
    ModflowDis(ml2, nlay=1, nrow=1, ncol=1, nper=1)
    with pytest.raises(NotImplementedError):
        MfUsgSgb.load(str(sgb_sp), ml2, nper=1, ext_unit_dict={})


# ---------------------------------------------------------------------------
# MfUsgQrt tests (Sink with Return Flow, gwf2QRT8u.f)
# ---------------------------------------------------------------------------

def test_mfusgqrt_minimal_authoring(function_tmpdir):
    """Minimal QRT authored from scratch: 0-based API, 1-based file, reload."""
    from flopy.modflow import ModflowDis

    ml = MfUsg(structured=False, model_ws=str(function_tmpdir))
    ModflowDis(ml, nlay=1, nrow=1, ncol=1, nper=1)

    dtype = MfUsgQrt.get_default_dtype(returnflow=True, changec=False)
    spd = {0: np.array([(0, -100.0, 0.75)], dtype=dtype).view(np.recarray)}
    recips = {0: [[9, 10]]}  # 0-based recipient nodes
    qrt = MfUsgQrt(
        ml, stress_period_data=spd, recipient_nodes=recips, options=["RETURNFLOW"]
    )
    qrt.fn_path = str(function_tmpdir / "created.qrt")
    qrt.write_file()
    text = Path(qrt.fn_path).read_text()

    # Sink node 0 -> 1 in file; recipients 9,10 -> 10,11 in file
    assert "\n 1  -1.000000e+02  2  7.500000e-01" in text
    assert "INTERNAL" in text
    assert "\n 10 11\n" in text

    ml2 = MfUsg(structured=False, model_ws=str(function_tmpdir))
    ModflowDis(ml2, nlay=1, nrow=1, ncol=1, nper=1)
    qrt2 = MfUsgQrt.load(qrt.fn_path, ml2, nper=1)
    rec = qrt2.stress_period_data[0]
    assert list(rec["node"]) == [0]
    assert np.isclose(rec["q"][0], -100.0)
    assert np.isclose(rec["rfprop"][0], 0.75)
    assert qrt2.recipient_nodes[0][0] == [9, 10]


def test_mfusgqrt_returnflow_concentration(function_tmpdir):
    """QRT CHANGEC + AUX concentration round-trips (line order: rfprop, ichng, aux)."""
    from flopy.modflow import ModflowDis

    ml = MfUsg(structured=False, model_ws=str(function_tmpdir))
    ModflowDis(ml, nlay=1, nrow=1, ncol=1, nper=1)

    opts = ["RETURNFLOW", "CHANGEC", "AUX C01"]
    dtype = MfUsgQrt.get_default_dtype(
        returnflow=True, changec=True, aux_names=["C01"]
    )
    spd = {0: np.array([(3, -25.0, 0.5, 2, 0.15)], dtype=dtype).view(np.recarray)}
    recips = {0: [[7]]}
    qrt = MfUsgQrt(
        ml, stress_period_data=spd, recipient_nodes=recips, options=opts
    )
    qrt.fn_path = str(function_tmpdir / "conc.qrt")
    qrt.write_file()
    text = qrt.fn_path
    content = Path(text).read_text()
    assert "RETURNFLOW CHANGEC AUX C01" in content

    ml2 = MfUsg(structured=False, model_ws=str(function_tmpdir))
    ModflowDis(ml2, nlay=1, nrow=1, ncol=1, nper=1)
    qrt2 = MfUsgQrt.load(text, ml2, nper=1)
    rec = qrt2.stress_period_data[0]
    assert "iqchngtyp" in qrt2.dtype.names
    assert "C01" in qrt2.dtype.names
    assert rec["iqchngtyp"][0] == 2
    assert np.isclose(rec["rfprop"][0], 0.5)
    assert np.isclose(rec["C01"][0], 0.15)
    assert qrt2.recipient_nodes[0][0] == [7]


def test_mfusgqrt_multi_recipient_and_pure_sink(function_tmpdir):
    """QRT mixing a multi-recipient sink with a no-return (NumRT=0) sink."""
    from flopy.modflow import ModflowDis

    ml = MfUsg(structured=False, model_ws=str(function_tmpdir))
    ModflowDis(ml, nlay=1, nrow=1, ncol=1, nper=1)

    dtype = MfUsgQrt.get_default_dtype(returnflow=True)
    spd = {
        0: np.array(
            [(0, -200.0, 0.9), (4, -30.0, 0.0)], dtype=dtype
        ).view(np.recarray)
    }
    recips = {0: [[10, 11, 12], []]}
    qrt = MfUsgQrt(
        ml, stress_period_data=spd, recipient_nodes=recips, options=["RETURNFLOW"]
    )
    qrt.fn_path = str(function_tmpdir / "multi.qrt")
    qrt.write_file()

    ml2 = MfUsg(structured=False, model_ws=str(function_tmpdir))
    ModflowDis(ml2, nlay=1, nrow=1, ncol=1, nper=1)
    qrt2 = MfUsgQrt.load(qrt.fn_path, ml2, nper=1)
    assert qrt2.recipient_nodes[0][0] == [10, 11, 12]
    assert qrt2.recipient_nodes[0][1] == []
    # MXRTCELLS in the header equals the total recipient nodes in the SP (3)
    header = Path(qrt.fn_path).read_text().splitlines()[1]
    assert header.split()[1] == "3"


def test_mfusgqrt_nam_registry():
    """MfUsg.load() registry maps the QRT package key to MfUsgQrt."""
    ml = MfUsg(structured=False)
    assert ml.mfnam_packages["qrt"] is MfUsgQrt


def test_mfusgqrt_unsupported_modes_fail_explicitly(function_tmpdir):
    """NPQRT>0 and TRANSIENTQ raise NotImplementedError rather than partial-write."""
    from flopy.modflow import ModflowDis

    ml = MfUsg(structured=False, model_ws=str(function_tmpdir))
    ModflowDis(ml, nlay=1, nrow=1, ncol=1, nper=1)

    # NPQRT > 0 (named parameters)
    qrt_param = function_tmpdir / "param.qrt"
    qrt_param.write_text(
        "# qrt params\n"
        "        10        10 0 1 5 RETURNFLOW\n"
        " 0    Stress Period 1\n"
    )
    with pytest.raises(NotImplementedError):
        MfUsgQrt.load(str(qrt_param), ml, nper=1, ext_unit_dict={})

    # TRANSIENTQ option
    qrt_tq = function_tmpdir / "tq.qrt"
    qrt_tq.write_text(
        "# qrt transientq\n"
        "        10        10 0 0 0 RETURNFLOW TRANSIENTQ 5\n"
        " 0    Stress Period 1\n"
    )
    with pytest.raises(NotImplementedError):
        MfUsgQrt.load(str(qrt_tq), ml, nper=1, ext_unit_dict={})


# ---------------------------------------------------------------------------
# MfUsgDrt tests (Drain Return DRT8, gwf2drt8u.f)
# ---------------------------------------------------------------------------

def test_mfusgdrt_inline_single_recipient(function_tmpdir):
    """DRT with a single inline recipient (NR>0): 0-based API, 1-based file."""
    from flopy.modflow import ModflowDis

    ml = MfUsg(structured=False, model_ws=str(function_tmpdir))
    ModflowDis(ml, nlay=1, nrow=1, ncol=1, nper=1)

    dtype = MfUsgDrt.get_usg_dtype(returnflow=True)
    spd = {0: np.array([(0, 5.0, 100.0, 0.7)], dtype=dtype).view(np.recarray)}
    recips = {0: [[8]]}
    drt = MfUsgDrt(
        ml, stress_period_data=spd, recipient_nodes=recips, options=["RETURNFLOW"]
    )
    drt.fn_path = str(function_tmpdir / "created.drt")
    drt.write_file()
    text = Path(drt.fn_path).read_text()
    # drain node 0->1, recipient node 8->9, inline (no INTERNAL block)
    assert "\n 1  5.000000e+00  1.000000e+02  9  7.000000e-01" in text
    assert "INTERNAL" not in text

    ml2 = MfUsg(structured=False, model_ws=str(function_tmpdir))
    ModflowDis(ml2, nlay=1, nrow=1, ncol=1, nper=1)
    drt2 = MfUsgDrt.load(drt.fn_path, ml2, nper=1)
    rec = drt2.stress_period_data[0]
    assert list(rec["node"]) == [0]
    assert np.isclose(rec["elev"][0], 5.0)
    assert np.isclose(rec["cond"][0], 100.0)
    assert np.isclose(rec["rfprop"][0], 0.7)
    assert drt2.recipient_nodes[0][0] == [8]


def test_mfusgdrt_changec_concentration(function_tmpdir):
    """DRT with CHANGEC (IDCHNGTYP) + AUX concentration round-trips."""
    from flopy.modflow import ModflowDis

    ml = MfUsg(structured=False, model_ws=str(function_tmpdir))
    ModflowDis(ml, nlay=1, nrow=1, ncol=1, nper=1)

    opts = ["RETURNFLOW", "CHANGEC", "AUX C01"]
    dtype = MfUsgDrt.get_usg_dtype(
        returnflow=True, changec=True, aux_names=["C01"]
    )
    spd = {0: np.array([(2, 6.0, 75.0, 0.4, 3, 0.2)], dtype=dtype).view(np.recarray)}
    recips = {0: [[5]]}
    drt = MfUsgDrt(
        ml, stress_period_data=spd, recipient_nodes=recips, options=opts
    )
    drt.fn_path = str(function_tmpdir / "conc.drt")
    drt.write_file()
    assert "RETURNFLOW CHANGEC AUX C01" in Path(drt.fn_path).read_text()

    ml2 = MfUsg(structured=False, model_ws=str(function_tmpdir))
    ModflowDis(ml2, nlay=1, nrow=1, ncol=1, nper=1)
    drt2 = MfUsgDrt.load(drt.fn_path, ml2, nper=1)
    rec = drt2.stress_period_data[0]
    assert "idchngtyp" in drt2.dtype.names
    assert rec["idchngtyp"][0] == 3
    assert np.isclose(rec["C01"][0], 0.2)
    assert drt2.recipient_nodes[0][0] == [5]


def test_mfusgdrt_spread_multi_node(function_tmpdir):
    """DRT spreading ground (NR<0): multiple recipients via a U1DINT block."""
    from flopy.modflow import ModflowDis

    ml = MfUsg(structured=False, model_ws=str(function_tmpdir))
    ModflowDis(ml, nlay=1, nrow=1, ncol=1, nper=1)

    dtype = MfUsgDrt.get_usg_dtype(returnflow=True)
    spd = {0: np.array([(0, 5.0, 100.0, 0.9)], dtype=dtype).view(np.recarray)}
    recips = {0: [[10, 11, 12]]}
    drt = MfUsgDrt(
        ml, stress_period_data=spd, recipient_nodes=recips, options=["RETURNFLOW"]
    )
    drt.fn_path = str(function_tmpdir / "spread.drt")
    drt.write_file()
    content = Path(drt.fn_path).read_text()
    lines = content.splitlines()
    # NR = -3 (spreading), SPREAD 3 in header, U1DINT block with 1-based nodes
    assert lines[1].split()[-2:] == ["SPREAD", "3"]
    assert "1.000000e+02  -3" in content  # NR = -3 on the drain line
    assert "INTERNAL" in content
    assert "\n 11 12 13\n" in content

    ml2 = MfUsg(structured=False, model_ws=str(function_tmpdir))
    ModflowDis(ml2, nlay=1, nrow=1, ncol=1, nper=1)
    drt2 = MfUsgDrt.load(drt.fn_path, ml2, nper=1)
    assert drt2.recipient_nodes[0][0] == [10, 11, 12]


def test_mfusgdrt_stress_period_reuse(function_tmpdir):
    """DRT -1 reuse copies the previous stress period's drains."""
    from flopy.modflow import ModflowDis

    drt_in = function_tmpdir / "reuse.drt"
    drt_in.write_text(
        "# drt reuse\n"
        "         1 0 0 0 RETURNFLOW\n"
        " 1 0    Stress Period 1\n"
        " 1  5.000000e+00  1.000000e+02  9  7.000000e-01\n"
        " -1    Stress Period 2\n"
    )
    ml = MfUsg(structured=False, model_ws=str(function_tmpdir))
    ModflowDis(ml, nlay=1, nrow=1, ncol=1, nper=2)
    drt = MfUsgDrt.load(str(drt_in), ml, nper=2, ext_unit_dict={})
    assert list(drt.stress_period_data[0]["node"]) == [0]
    assert list(drt.stress_period_data[1]["node"]) == [0]
    assert drt.recipient_nodes[1][0] == [8]


def test_mfusgdrt_nam_registry():
    """MfUsg.load() registry maps the DRT package key to MfUsgDrt."""
    ml = MfUsg(structured=False)
    assert ml.mfnam_packages["drt"] is MfUsgDrt


def test_mfusgdrt_parameters_fail_explicitly(function_tmpdir):
    """Named DRT parameters (NPDRT>0) raise NotImplementedError on load."""
    from flopy.modflow import ModflowDis

    drt_param = function_tmpdir / "param.drt"
    drt_param.write_text(
        "# drt params\n"
        "        10 0 1 5 RETURNFLOW\n"
        " 0    Stress Period 1\n"
    )
    ml = MfUsg(structured=False, model_ws=str(function_tmpdir))
    ModflowDis(ml, nlay=1, nrow=1, ncol=1, nper=1)
    with pytest.raises(NotImplementedError):
        MfUsgDrt.load(str(drt_param), ml, nper=1, ext_unit_dict={})


# ---------------------------------------------------------------------------
# BCF / LPF TABRICH tests (items 1c IUZONTAB + 1d RETCRVS, gwf2bcf-lpf-u1.f)
# ---------------------------------------------------------------------------

_RETCRVS_2x2 = np.array(
    [
        [[0.0, 1.0, 1.0], [1.0, 0.5, 0.10]],  # zone 1: (caphead, sat, relperm)
        [[0.0, 1.0, 1.0], [2.0, 0.3, 0.05]],  # zone 2
    ]
)


def test_mfusgbcf_tabrich_authoring_roundtrip(function_tmpdir):
    """BCF TABRICH: author IUZONTAB + RETCRVS from scratch, write, reload."""
    from flopy.modflow import ModflowDis

    ml = MfUsg(model_ws=str(function_tmpdir))
    ModflowDis(ml, nlay=2, nrow=1, ncol=1, nper=1, steady=True)
    MfUsgBas(ml, ibound=1, strt=1.0, richards=True)
    bcf = MfUsgBcf(
        ml,
        laycon=[0, 5],
        ipakcb=0,
        tabrich=True,
        nuzones=2,
        nutabrows=2,
        iuzontab=[1, 2],
        retcrvs=_RETCRVS_2x2,
        tran=1.0,
        hy=1.0,
        kv=1.0,
        sf1=1.0e-5,
        sf2=0.15,
    )
    bcf.fn_path = str(function_tmpdir / "tabrich.bcf")
    bcf.write_file()
    text = Path(bcf.fn_path).read_text()
    assert "TABRICH" in text.splitlines()[0]  # BCF writes item 1 first (no heading)
    assert "#iuzontab" in text

    ml2 = MfUsg(model_ws=str(function_tmpdir))
    ModflowDis(ml2, nlay=2, nrow=1, ncol=1, nper=1, steady=True)
    MfUsgBas(ml2, ibound=1, strt=1.0, richards=True)
    bcf2 = MfUsgBcf.load(bcf.fn_path, ml2)
    assert bcf2.tabrich and bcf2.nuzones == 2 and bcf2.nutabrows == 2
    assert list(bcf2.iuzontab.array) == [1, 2]
    assert bcf2.retcrvs.shape == (2, 2, 3)
    assert np.allclose(bcf2.retcrvs, _RETCRVS_2x2)


def test_mfusglpf_tabrich_authoring_roundtrip(function_tmpdir):
    """LPF TABRICH: author IUZONTAB + RETCRVS, skip Richards arrays, reload."""
    from flopy.modflow import ModflowDis

    ml = MfUsg(model_ws=str(function_tmpdir))
    ModflowDis(ml, nlay=2, nrow=1, ncol=1, nper=1, steady=True)
    MfUsgBas(ml, ibound=1, strt=1.0, richards=True)
    lpf = MfUsgLpf(
        ml,
        laytyp=[0, 5],
        ipakcb=0,
        tabrich=True,
        nuzones=2,
        nutabrows=2,
        iuzontab=[1, 2],
        retcrvs=_RETCRVS_2x2,
        hk=1.0,
        vka=1.0,
        ss=1.0e-5,
        sy=0.15,
    )
    lpf.fn_path = str(function_tmpdir / "tabrich.lpf")
    lpf.write_file(check=False)
    text = Path(lpf.fn_path).read_text()
    assert "TABRICH 2 2" in text
    assert "#iuzontab" in text
    # Richards per-layer arrays must NOT be written under TABRICH
    assert "richards alpha" not in text

    ml2 = MfUsg(model_ws=str(function_tmpdir))
    ModflowDis(ml2, nlay=2, nrow=1, ncol=1, nper=1, steady=True)
    MfUsgBas(ml2, ibound=1, strt=1.0, richards=True)
    lpf2 = MfUsgLpf.load(lpf.fn_path, ml2, check=False)
    assert lpf2.tabrich and lpf2.nuzones == 2 and lpf2.nutabrows == 2
    assert list(lpf2.iuzontab.array) == [1, 2]
    assert np.allclose(lpf2.retcrvs, _RETCRVS_2x2)


def test_mfusgbcf_tabrich_incomplete_write_fails(function_tmpdir):
    """tabrich=True without iuzontab/retcrvs: construction OK, write fails."""
    from flopy.modflow import ModflowDis

    ml = MfUsg(model_ws=str(function_tmpdir))
    ModflowDis(ml, nlay=2, nrow=1, ncol=1, nper=1, steady=True)
    MfUsgBas(ml, ibound=1, strt=1.0, richards=True)
    # Construction must NOT raise (BCF can be used purely as model setup).
    bcf = MfUsgBcf(
        ml,
        laycon=[0, 5],
        ipakcb=0,
        tabrich=True,
        nuzones=1,
        nutabrows=2,
        tran=1.0,
        hy=1.0,
        kv=1.0,
        sf1=1.0e-5,
        sf2=0.15,
    )
    bcf.fn_path = str(function_tmpdir / "bad.bcf")
    with pytest.raises(ValueError):
        bcf.write_file()


def test_tabrich_retcrvs_shape_validation(function_tmpdir):
    """A mis-shaped RETCRVS raises ValueError at construction."""
    from flopy.modflow import ModflowDis

    ml = MfUsg(model_ws=str(function_tmpdir))
    ModflowDis(ml, nlay=2, nrow=1, ncol=1, nper=1, steady=True)
    MfUsgBas(ml, ibound=1, strt=1.0, richards=True)
    with pytest.raises(ValueError):
        MfUsgBcf(
            ml,
            laycon=[0, 5],
            ipakcb=0,
            tabrich=True,
            nuzones=2,
            nutabrows=2,
            iuzontab=[1, 2],
            retcrvs=np.zeros((2, 3, 3)),  # nutabrows mismatch (expects (2,2,3))
            tran=1.0,
            hy=1.0,
            kv=1.0,
            sf1=1.0e-5,
            sf2=0.15,
        )


# ---------------------------------------------------------------------------
# Priority-2 closeout: BAS niche options, ETS NETSOP=2/IESFACTOR, HFB params
# ---------------------------------------------------------------------------

def test_mfusgbas_richards_hp_and_ihm_roundtrip(function_tmpdir):
    """BAS RICHARDS_HP and IHM [IUIHM] options author and round-trip."""
    from flopy.modflow import ModflowDis

    ml = MfUsg(model_ws=str(function_tmpdir))
    ModflowDis(ml, nlay=1, nrow=1, ncol=1, nper=1)
    bas = MfUsgBas(ml, ibound=1, strt=1.0, richards_hp=True, ihm=True, iuihm=44)
    # RICHARDS_HP implies Richards mode for dependent packages (BCF/LPF LAYTYP=5)
    assert bas.richards is True
    bas.fn_path = str(function_tmpdir / "rhp.bas")
    bas.write_file(check=False)
    opt_line = next(ln for ln in Path(bas.fn_path).read_text().splitlines()
                    if "RICHARDS" in ln)
    assert "RICHARDS_HP" in opt_line
    assert "IHM 44" in opt_line

    ml2 = MfUsg(model_ws=str(function_tmpdir))
    ModflowDis(ml2, nlay=1, nrow=1, ncol=1, nper=1)
    bas2 = MfUsgBas.load(bas.fn_path, ml2, check=False)
    assert bas2.richards_hp is True
    assert bas2.richards is True  # effective Richards mode preserved
    assert bas2.ihm is True
    assert bas2.iuihm == 44


def test_mfusgets_netsop2_authoring_roundtrip(function_tmpdir):
    """ETS NETSOP=2 (layer-indicator IEVT) authors and round-trips."""
    from flopy.mfusg import MfUsgEts
    from flopy.modflow import ModflowDis

    ml = MfUsg(structured=True, model_ws=str(function_tmpdir))
    ModflowDis(ml, nlay=2, nrow=2, ncol=2, nper=1)
    ets = MfUsgEts(ml, netsop=2, evtr=1.2e-4, netseg=1, ievt=1)
    ml.write_input()

    content = Path(ets.fn_path).read_text()
    # Item 2a: NETSOP IETSCB NPETS NETSEG IESFACTOR -> NETSOP field is 2
    assert content.splitlines()[1].split()[0] == "2"

    ml2 = MfUsg(structured=True, model_ws=str(function_tmpdir))
    ModflowDis(ml2, nlay=2, nrow=2, ncol=2, nper=1)
    ets2 = MfUsgEts.load(ets.fn_path, ml2, nper=1, ext_unit_dict={})
    assert ets2.netsop == 2


def test_mfusgets_iesfactor_authoring(function_tmpdir):
    """ETS IESFACTOR=1 writes the transport ESFACTOR record when transport is on."""
    from flopy.mfusg import MfUsgEts
    from flopy.modflow import ModflowDis

    ml = MfUsg(structured=True, model_ws=str(function_tmpdir))
    ModflowDis(ml, nlay=1, nrow=2, ncol=2, nper=1)
    # Emulate an active transport simulation (BCT would set these).
    ml.itrnsp = 1
    ml.mcomp = 1
    ets = MfUsgEts(
        ml, netsop=1, evtr=1.2e-4, netseg=1, iesfactor=1, esfactor=[2.5]
    )
    ets.fn_path = str(function_tmpdir / "ies.ets")
    ets.write_file()
    content = Path(ets.fn_path).read_text()
    # Item 2a IESFACTOR field is 1; ESFACTOR(MCOMP) record written
    assert content.splitlines()[1].split()[4] == "1"
    assert "2.500000" in content


def test_mfusghfb_parameterized_fails_explicitly(function_tmpdir):
    """Parameterized HFB (NPHFB>0) fails explicitly on load rather than partial."""
    from flopy.mfusg import MfUsgHfb
    from flopy.modflow import ModflowDis

    # Header: NPHFB MXFB NHFBNP -> NPHFB=1 (named parameters), unsupported
    hfb_param = function_tmpdir / "param.hfb"
    hfb_param.write_text(
        "# parameterized HFB\n"
        "         1         1         0\n"
    )
    ml = MfUsg(structured=False, model_ws=str(function_tmpdir))
    ModflowDis(ml, nlay=1, nrow=1, ncol=1, nper=1)
    with pytest.raises(NotImplementedError):
        MfUsgHfb.load(str(hfb_param), ml, ext_unit_dict={})


def test_mfusgdpt_aw_adsorbim_fails_explicitly(function_tmpdir):
    """DPT immobile-domain air-water adsorption (A-W_ADSORBIM) fails explicitly.

    The sub-mode reads extra function indices and arrays; without support it
    would silently shift all later reads, so load raises instead.
    """
    from flopy.mfusg import MfUsgDpt
    from flopy.modflow import ModflowDis

    # Two option-line forms: the bare keyword and the form with function
    # indices (IAREA_FNIM=5, IKAWI_FNIM=4 = tabular, which would otherwise read
    # a zone map + tabular area arrays). Both must raise at the option line,
    # before any extra array read, so nothing downstream shifts.
    for header in (
        " 0 0 0 0 0 0 0 A-W_ADSORBIM\n",
        " 0 0 0 0 0 0 0 A-W_ADSORBIM 5 4\n",
    ):
        dpt_file = function_tmpdir / "aw.dpt"
        dpt_file.write_text("# DPT immobile air-water adsorption\n" + header)
        ml = MfUsg(structured=False, model_ws=str(function_tmpdir))
        ModflowDis(ml, nlay=1, nrow=1, ncol=1, nper=1)
        with pytest.raises(NotImplementedError, match="A-W_ADSORBIM"):
            MfUsgDpt.load(str(dpt_file), ml, ext_unit_dict={})


# ---------------------------------------------------------------------------
# Priority-3 review: from-scratch authoring tests for BCT / DDF
# (previously only exercised via real-model round-trips)
# ---------------------------------------------------------------------------

def _minimal_transport_model(ws, nrow=2, ncol=2):
    """A minimal structured MfUsg + DIS + BAS for BCT/DDF authoring tests."""
    from flopy.modflow import ModflowDis

    ml = MfUsg(model_ws=str(ws))
    ModflowDis(ml, nlay=1, nrow=nrow, ncol=ncol, nper=1, steady=True)
    MfUsgBas(ml, ibound=1, strt=1.0)
    return ml


def test_mfusgbct_minimal_authoring_roundtrip(function_tmpdir):
    """BCT single-species (IDISP=1) authored from scratch, written, reloaded."""
    ml = _minimal_transport_model(function_tmpdir)
    bct = MfUsgBct(ml, mcomp=1, idisp=1, prsity=0.2, conc=0.0)
    bct.fn_path = str(function_tmpdir / "min.bct")
    bct.write_file()

    ml2 = _minimal_transport_model(function_tmpdir)
    bct2 = MfUsgBct.load(bct.fn_path, ml2)
    assert bct2.mcomp == 1
    assert bct2.idisp == 1


def test_mfusgbct_idisp2_authoring_roundtrip(function_tmpdir):
    """BCT IDISP=2 (full dispersion tensor DLX/DLY/DLZ/DTXY/DTYZ/DTXZ) authoring."""
    ml = _minimal_transport_model(function_tmpdir)
    bct = MfUsgBct(
        ml, mcomp=1, idisp=2, prsity=0.2,
        dlx=1.0, dly=1.0, dlz=0.1, dtxy=0.1, dtyz=0.1, dtxz=0.1, conc=0.0,
    )
    bct.fn_path = str(function_tmpdir / "idisp2.bct")
    bct.write_file()

    ml2 = _minimal_transport_model(function_tmpdir)
    bct2 = MfUsgBct.load(bct.fn_path, ml2)
    assert bct2.idisp == 2


def test_mfusgbct_multispecies_authoring_roundtrip(function_tmpdir):
    """BCT multi-species (MCOMP=2) authoring round-trips with species isolated."""
    ml = _minimal_transport_model(function_tmpdir)
    bct = MfUsgBct(ml, mcomp=2, idisp=1, prsity=0.2, conc=0.0)
    bct.fn_path = str(function_tmpdir / "multi.bct")
    bct.write_file()

    ml2 = _minimal_transport_model(function_tmpdir)
    bct2 = MfUsgBct.load(bct.fn_path, ml2)
    assert bct2.mcomp == 2


def test_mfusgddf_nonlinear_table_authoring_roundtrip(function_tmpdir):
    """DDF NONLINEAR density table authored from scratch, written, reloaded."""
    ml = _minimal_transport_model(function_tmpdir)
    table = [(0.0, 1000.0), (17.5, 1012.5), (35.0, 1025.0)]
    ddf = MfUsgDdf(ml, rhofresh=1000.0, rhostd=1025.0, cstd=35.0,
                   nonlinear=True, density_table=table)
    ddf.fn_path = str(function_tmpdir / "nl.ddf")
    ddf.write_file()
    content = Path(ddf.fn_path).read_text()
    assert "NONLINEAR" in content

    ml2 = _minimal_transport_model(function_tmpdir)
    ddf2 = MfUsgDdf.load(ddf.fn_path, ml2)
    assert ddf2.nonlinear
    assert len(ddf2.density_table) == 3
    assert np.isclose(ddf2.density_table[1][1], 1012.5, atol=1e-3)


# ===========================================================================
# Phase 2 hardening (USGT_PHASE2_REVIEW.md)
# ===========================================================================

# --- P0: MfUsgDrt structured authoring contract ---------------------------

def test_mfusgdrt_structured_construction_raises(function_tmpdir):
    """MfUsgDrt is unstructured-only; structured construction fails clearly."""
    from flopy.modflow import ModflowDis

    ml = MfUsg(structured=True, model_ws=str(function_tmpdir))
    ModflowDis(ml, nlay=1, nrow=2, ncol=2, nper=1)
    with pytest.raises(NotImplementedError, match="unstructured"):
        MfUsgDrt(ml, stress_period_data={0: [(0, 5.0, 100.0)]})


def test_modflowdrt_structured_authoring_from_scratch(function_tmpdir):
    """Classic structured DRT (base ModflowDrt) authors a valid file from scratch."""
    from flopy.modflow import ModflowDis, ModflowDrt

    ml = MfUsg(structured=True, model_ws=str(function_tmpdir))
    ModflowDis(ml, nlay=1, nrow=2, ncol=2, nper=1)
    dtype = ModflowDrt.get_default_dtype(structured=True)
    spd = {0: np.array([(0, 0, 0, 5.0, 100.0, 0, 0, 0, 0.0)],
                       dtype=dtype).view(np.recarray)}
    drt = ModflowDrt(ml, stress_period_data=spd)
    drt.fn_path = str(function_tmpdir / "classic.drt")
    drt.write_file()
    assert Path(drt.fn_path).read_text().strip() != ""


def test_mfusgdrt_load_structured_delegates_to_base(function_tmpdir):
    """Loading a structured DRT via the registry class returns a base object."""
    from flopy.modflow import ModflowDis, ModflowDrt

    drt_in = function_tmpdir / "classic.drt"
    drt_in.write_text(
        "# classic structured DRT\n"
        "         1         0\n"
        " 1 0    Stress Period 1\n"
        " 1 1 1  5.000000  1.000000e+02\n"
    )
    ml = MfUsg(structured=True, model_ws=str(function_tmpdir))
    ModflowDis(ml, nlay=1, nrow=2, ncol=2, nper=1)
    drt = MfUsgDrt.load(str(drt_in), ml, nper=1, ext_unit_dict={})
    # Delegation yields a base ModflowDrt, not the unstructured MfUsgDrt
    assert isinstance(drt, ModflowDrt)
    assert not isinstance(drt, MfUsgDrt)


# --- P0: SGB/QRT/DRT main-list controls (SFAC / OPEN-CLOSE / EXTERNAL) ------

def _usgt_unstructured_model(ws, nper=1):
    from flopy.modflow import ModflowDis

    ml = MfUsg(structured=False, model_ws=str(ws))
    ModflowDis(ml, nlay=1, nrow=1, ncol=1, nper=nper)
    return ml


def test_mfusgsgb_sfac_and_open_close(function_tmpdir):
    """SGB honors SFAC (inert on gradient per Fortran) and OPEN/CLOSE rows."""
    # SFAC: USG-T scales an internal dummy column for SGB, not the gradient.
    sfac_file = function_tmpdir / "sfac.sgb"
    sfac_file.write_text(
        "# sgb sfac\n         1 0\n 1 0    SP1\n SFAC 5.0\n 101   1.000000e-02\n"
    )
    sgb = MfUsgSgb.load(str(sfac_file), _usgt_unstructured_model(function_tmpdir),
                        nper=1, ext_unit_dict={})
    assert np.isclose(sgb.stress_period_data[0]["gradient"][0], 0.01)

    # OPEN/CLOSE: rows live in a separate file resolved against model_ws.
    (function_tmpdir / "sgb_rows.dat").write_text(
        " 201   2.000000e-02\n 202   3.000000e-02\n"
    )
    oc_file = function_tmpdir / "oc.sgb"
    oc_file.write_text(
        "# sgb open/close\n         2 0\n 2 0    SP1\n OPEN/CLOSE sgb_rows.dat\n"
    )
    sgb2 = MfUsgSgb.load(str(oc_file), _usgt_unstructured_model(function_tmpdir),
                         nper=1, ext_unit_dict={})
    assert list(sgb2.stress_period_data[0]["node"]) == [200, 201]


def test_mfusgqrt_sfac_with_recipients_expanded_write(function_tmpdir):
    """QRT SFAC scales Q, keeps recipients, and writes expanded (no SFAC)."""
    sfac_file = function_tmpdir / "sfac.qrt"
    sfac_file.write_text(
        "# qrt sfac\n        1         1 0 0 0 RETURNFLOW\n 1    SP1\n"
        " SFAC 2.0\n 5  -5.000000e+01  1  8.000000e-01\nINTERNAL  1  (FREE)  -1\n 10\n"
    )
    qrt = MfUsgQrt.load(str(sfac_file), _usgt_unstructured_model(function_tmpdir),
                        nper=1, ext_unit_dict={})
    rec = qrt.stress_period_data[0]
    assert np.isclose(rec["q"][0], -100.0)  # -50 * SFAC 2.0
    assert qrt.recipient_nodes[0][0] == [9]

    out = function_tmpdir / "out.qrt"
    qrt.fn_path = str(out)
    qrt.write_file()
    text = out.read_text()
    assert "SFAC" not in text  # expanded valid write
    qrt2 = MfUsgQrt.load(str(out), _usgt_unstructured_model(function_tmpdir),
                         nper=1, ext_unit_dict={})
    assert np.isclose(qrt2.stress_period_data[0]["q"][0], -100.0)
    assert qrt2.recipient_nodes[0][0] == [9]


def test_mfusgdrt_sfac_single_and_spread(function_tmpdir):
    """DRT SFAC scales COND for both inline-single and spreading recipients."""
    sfac_file = function_tmpdir / "sfac.drt"
    sfac_file.write_text(
        "# drt sfac\n         2 0 0 0 RETURNFLOW\n 2    SP1\n SFAC 3.0\n"
        " 1  5.000000e+00  1.000000e+01  9  7.000000e-01\n"
        " 2  4.000000e+00  2.000000e+01  -2  5.000000e-01\n"
        "INTERNAL  1  (FREE)  -1\n 11 12\n"
    )
    drt = MfUsgDrt.load(str(sfac_file), _usgt_unstructured_model(function_tmpdir),
                        nper=1, ext_unit_dict={})
    rec = drt.stress_period_data[0]
    assert np.isclose(rec["cond"][0], 30.0)  # 10 * 3.0
    assert np.isclose(rec["cond"][1], 60.0)  # 20 * 3.0
    assert drt.recipient_nodes[0][0] == [8]        # inline single (NR=9 -> 8)
    assert drt.recipient_nodes[0][1] == [10, 11]   # spreading (NR=-2 -> nodes 11,12)


def test_usgt_list_external_without_dict_fails(function_tmpdir):
    """EXTERNAL list input fails explicitly (NotImplementedError) when the
    unit cannot be resolved, for each list-package family."""
    cases = [
        (MfUsgSgb, "ext.sgb", "# sgb\n         1 0\n 1 0    SP1\n EXTERNAL 77\n"),
        (MfUsgQrt, "ext.qrt",
         "# qrt\n        1         0 0 0 0\n 1    SP1\n EXTERNAL 77\n"),
        (MfUsgDrt, "ext.drt", "# drt\n         1 0 0 0\n 1    SP1\n EXTERNAL 77\n"),
    ]
    for cls, name, text in cases:
        p = function_tmpdir / name
        p.write_text(text)
        with pytest.raises(NotImplementedError, match="EXTERNAL"):
            cls.load(str(p), _usgt_unstructured_model(function_tmpdir),
                     nper=1, ext_unit_dict={})


# --- P1: BAS IHM optional integer parsing ---------------------------------

def test_mfusgbas_ihm_optional_unit_parsing(function_tmpdir):
    """BAS IHM loads with or without the optional IUIHM unit (no IndexError)."""
    from flopy.modflow import ModflowDis

    def _write_bas(optline, name):
        p = function_tmpdir / name
        p.write_text(
            f"# test\n{optline}\nCONSTANT          1\n"
            f"        -999.99\nCONSTANT       1.0\n"
        )
        return p

    def _model():
        ml = MfUsg(structured=False, model_ws=str(function_tmpdir))
        ModflowDis(ml, nlay=1, nrow=1, ncol=1, nper=1)
        return ml

    # bare IHM -> iuihm=0
    b = MfUsgBas.load(str(_write_bas("FREE IHM", "bare.bas")), _model(), check=False)
    assert b.ihm is True and b.iuihm == 0

    # IHM <int> -> iuihm=int
    b = MfUsgBas.load(str(_write_bas("FREE IHM 123", "unit.bas")), _model(),
                      check=False)
    assert b.ihm is True and b.iuihm == 123

    # IHM followed by another option -> iuihm=0, the other option still parsed
    b = MfUsgBas.load(str(_write_bas("FREE IHM SY-ALL", "opt.bas")), _model(),
                      check=False)
    assert b.ihm is True and b.iuihm == 0 and b.sy_all is True

    # programmatic write normalizes to "IHM <iuihm>" and round-trips
    for iuihm in (0, 99):
        ml = _model()
        bas = MfUsgBas(ml, ibound=1, strt=1.0, ihm=True, iuihm=iuihm)
        bas.fn_path = str(function_tmpdir / f"prog_{iuihm}.bas")
        bas.write_file(check=False)
        assert f"IHM {iuihm}" in Path(bas.fn_path).read_text()
        bas2 = MfUsgBas.load(bas.fn_path, _model(), check=False)
        assert bas2.ihm is True and bas2.iuihm == iuihm


# --- P1: QRT/DRT recipient_nodes validation -------------------------------

def test_mfusgqrt_recipient_count_mismatch_fails(function_tmpdir):
    """QRT recipient_nodes shorter/longer than the sink list fails on write."""
    ml = _usgt_unstructured_model(function_tmpdir)
    dtype = MfUsgQrt.get_default_dtype(returnflow=True)
    spd = {0: np.array([(0, -100.0, 0.5), (4, -50.0, 0.5)],
                       dtype=dtype).view(np.recarray)}

    too_short = MfUsgQrt(ml, stress_period_data=spd,
                         recipient_nodes={0: [[9]]}, options=["RETURNFLOW"])
    too_short.fn_path = str(function_tmpdir / "short.qrt")
    with pytest.raises(ValueError, match="recipient_nodes"):
        too_short.write_file()

    ml2 = _usgt_unstructured_model(function_tmpdir)
    too_long = MfUsgQrt(ml2, stress_period_data=spd,
                        recipient_nodes={0: [[9], [10], [11]]},
                        options=["RETURNFLOW"])
    too_long.fn_path = str(function_tmpdir / "long.qrt")
    with pytest.raises(ValueError, match="recipient_nodes"):
        too_long.write_file()


def test_mfusgdrt_recipient_count_mismatch_fails(function_tmpdir):
    """DRT recipient_nodes shorter/longer than the drain list fails on write."""
    ml = _usgt_unstructured_model(function_tmpdir)
    dtype = MfUsgDrt.get_usg_dtype(returnflow=True)
    spd = {0: np.array([(0, 5.0, 100.0, 0.5), (4, 4.0, 50.0, 0.5)],
                       dtype=dtype).view(np.recarray)}

    too_short = MfUsgDrt(ml, stress_period_data=spd,
                         recipient_nodes={0: [[9]]}, options=["RETURNFLOW"])
    too_short.fn_path = str(function_tmpdir / "short.drt")
    with pytest.raises(ValueError, match="recipient_nodes"):
        too_short.write_file()

    ml2 = _usgt_unstructured_model(function_tmpdir)
    too_long = MfUsgDrt(ml2, stress_period_data=spd,
                        recipient_nodes={0: [[9], [10], [11]]},
                        options=["RETURNFLOW"])
    too_long.fn_path = str(function_tmpdir / "long.drt")
    with pytest.raises(ValueError, match="recipient_nodes"):
        too_long.write_file()


def test_mfusgqrt_zero_recipients_when_omitted(function_tmpdir):
    """QRT with RETURNFLOW but omitted recipient_nodes => all-zero recipients."""
    ml = _usgt_unstructured_model(function_tmpdir)
    dtype = MfUsgQrt.get_default_dtype(returnflow=True)
    spd = {0: np.array([(0, -100.0, 0.0), (4, -50.0, 0.0)],
                       dtype=dtype).view(np.recarray)}
    qrt = MfUsgQrt(ml, stress_period_data=spd, options=["RETURNFLOW"])
    qrt.fn_path = str(function_tmpdir / "zero.qrt")
    qrt.write_file()  # must not raise

    qrt2 = MfUsgQrt.load(qrt.fn_path, _usgt_unstructured_model(function_tmpdir),
                         nper=1, ext_unit_dict={})
    assert qrt2.recipient_nodes[0] == [[], []]


# --- P2: TABRICH node-count contract --------------------------------------

def test_tabrich_node_count_contract(function_tmpdir):
    """_tabrich.node_count: structured product, DISU nodes, DIS fallback, error."""
    from flopy.mfusg._tabrich import node_count
    from flopy.modflow import ModflowDis

    # structured DIS -> nlay*nrow*ncol
    ml = MfUsg(structured=True, model_ws=str(function_tmpdir))
    ModflowDis(ml, nlay=2, nrow=1, ncol=3, nper=1)
    assert node_count(ml) == 6

    # structured=False with a classic DIS -> fallback to grid product
    # (previously raised AttributeError assuming DISU)
    ml2 = MfUsg(structured=False, model_ws=str(function_tmpdir))
    ModflowDis(ml2, nlay=1, nrow=1, ncol=2, nper=1)
    assert node_count(ml2) == 2

    # neither DISU nor DIS -> explicit, actionable error
    ml3 = MfUsg(structured=False, model_ws=str(function_tmpdir))
    with pytest.raises(ValueError, match="node count"):
        node_count(ml3)


# ===========================================================================
# Phase 2 polish pass (USGT_PHASE2_REREVIEW.md)
# ===========================================================================

def test_usgt_list_open_close_quoted_filenames(function_tmpdir):
    """OPEN/CLOSE accepts single/double quotes and single-quoted names w/ spaces."""
    # SGB: single-quoted name
    (function_tmpdir / "rows.dat").write_text(" 201   2.000000e-02\n")
    p = function_tmpdir / "q.sgb"
    p.write_text("# sgb\n         1 0\n 1 0 SP1\n OPEN/CLOSE 'rows.dat'\n")
    sgb = MfUsgSgb.load(str(p), _usgt_unstructured_model(function_tmpdir),
                        nper=1, ext_unit_dict={})
    assert list(sgb.stress_period_data[0]["node"]) == [200]
    # expanded write emits inline rows; it must not preserve OPEN/CLOSE
    out = function_tmpdir / "out.sgb"
    sgb.fn_path = str(out)
    sgb.write_file()
    text = out.read_text()
    assert "OPEN/CLOSE" not in text
    assert "\n 201  2.000000e-02" in text

    # DRT: single-quoted name WITH SPACES, with a recipient
    (function_tmpdir / "rows with spaces.dat").write_text(
        " 1  5.000000e+00  1.000000e+01  9  7.000000e-01\n"
    )
    p2 = function_tmpdir / "sp.drt"
    p2.write_text(
        "# drt\n         1 0 0 0 RETURNFLOW\n 1 SP1\n"
        " OPEN/CLOSE 'rows with spaces.dat'\n"
    )
    drt = MfUsgDrt.load(str(p2), _usgt_unstructured_model(function_tmpdir),
                        nper=1, ext_unit_dict={})
    assert list(drt.stress_period_data[0]["node"]) == [0]
    assert drt.recipient_nodes[0][0] == [8]

    # SGB: double-quoted name
    p3 = function_tmpdir / "dq.sgb"
    p3.write_text('# sgb\n         1 0\n 1 0 SP1\n OPEN/CLOSE "rows.dat"\n')
    sgb3 = MfUsgSgb.load(str(p3), _usgt_unstructured_model(function_tmpdir),
                         nper=1, ext_unit_dict={})
    assert list(sgb3.stress_period_data[0]["node"]) == [200]


def test_usgt_list_external_positive_with_ext_unit_dict(function_tmpdir):
    """EXTERNAL via ext_unit_dict loads for SGB/QRT/DRT (one file begins with SFAC)."""
    from flopy.utils.mfreadnam import NamData

    def eud(unit, fname):
        # Minimal NAM entry: only filename/filetype are consulted on this path.
        return {unit: NamData("DATA", fname, None, {})}

    # SGB: external file begins with SFAC (inert on the SGB gradient).
    (function_tmpdir / "s_ext.dat").write_text(" SFAC 5.0\n 101   1.000000e-02\n")
    p = function_tmpdir / "e.sgb"
    p.write_text("# sgb\n         1 0\n 1 0 SP1\n EXTERNAL 77\n")
    sgb = MfUsgSgb.load(str(p), _usgt_unstructured_model(function_tmpdir),
                        nper=1, ext_unit_dict=eud(77, "s_ext.dat"))
    assert list(sgb.stress_period_data[0]["node"]) == [100]
    assert np.isclose(sgb.stress_period_data[0]["gradient"][0], 0.01)

    # QRT: external file begins with SFAC (scales Q: -50 * 2 = -100) + recipients.
    (function_tmpdir / "q_ext.dat").write_text(
        " SFAC 2.0\n 5  -5.000000e+01  1  8.000000e-01\n"
        "INTERNAL  1  (FREE)  -1\n 10\n"
    )
    p = function_tmpdir / "e.qrt"
    p.write_text("# qrt\n        1         1 0 0 0 RETURNFLOW\n 1 SP1\n EXTERNAL 78\n")
    qrt = MfUsgQrt.load(str(p), _usgt_unstructured_model(function_tmpdir),
                        nper=1, ext_unit_dict=eud(78, "q_ext.dat"))
    assert np.isclose(qrt.stress_period_data[0]["q"][0], -100.0)
    assert qrt.recipient_nodes[0][0] == [9]

    # DRT: external file with a spreading (NR<0) U1DINT block.
    (function_tmpdir / "d_ext.dat").write_text(
        " 1  5.000000e+00  1.000000e+01  -2  7.000000e-01\n"
        "INTERNAL  1  (FREE)  -1\n 11 12\n"
    )
    p = function_tmpdir / "e.drt"
    p.write_text("# drt\n         1 0 0 0 RETURNFLOW\n 1 SP1\n EXTERNAL 79\n")
    drt = MfUsgDrt.load(str(p), _usgt_unstructured_model(function_tmpdir),
                        nper=1, ext_unit_dict=eud(79, "d_ext.dat"))
    assert drt.recipient_nodes[0][0] == [10, 11]


def test_mfusgdrt_zero_recipients_when_omitted(function_tmpdir):
    """DRT with RETURNFLOW but omitted recipient_nodes => all-zero recipients."""
    ml = _usgt_unstructured_model(function_tmpdir)
    dtype = MfUsgDrt.get_usg_dtype(returnflow=True)
    spd = {0: np.array([(0, 5.0, 100.0, 0.0), (4, 4.0, 50.0, 0.0)],
                       dtype=dtype).view(np.recarray)}
    drt = MfUsgDrt(ml, stress_period_data=spd, options=["RETURNFLOW"])
    drt.fn_path = str(function_tmpdir / "zero.drt")
    drt.write_file()  # must not raise

    drt2 = MfUsgDrt.load(drt.fn_path, _usgt_unstructured_model(function_tmpdir),
                         nper=1, ext_unit_dict={})
    assert drt2.recipient_nodes[0] == [[], []]


# --- Stage 3 Card 2: parameter strategy (Expanded valid write for ETS) ----

def test_mfusgets_parameterized_load_expands_to_npets0(function_tmpdir):
    """A parameterized ETS file loads, expands the parameter to a concrete
    ETSR array, and writes valid non-parametric input (NPETS=0, no PARAMETER).

    This is the documented `Expanded valid write` policy: ETS reads MODFLOW
    array-parameter syntax but does not preserve it on output.
    """
    from flopy.mfusg import MfUsgEts
    from flopy.modflow import ModflowDis

    ml = MfUsg(structured=True, model_ws=str(function_tmpdir))
    ModflowDis(ml, nlay=1, nrow=2, ncol=2, nper=1)
    param_ets = function_tmpdir / "param.ets"
    param_ets.write_text(
        "# parameterized ETS (NPETS=1)\n"
        "PARAMETER 1\n"
        "1 0 1 1 0\n"            # NETSOP IETSCB NPETS NETSEG IESFACTOR
        "etsrate ets 5.0E-4 1\n"  # param: name type value nclu
        "NONE ALL\n"             # cluster: no multiplier, all cells
        "0 1 0\n"                 # SP1: INSURF INETSR INEXDP (INETSR=1 param)
        "CONSTANT 10.0\n"        # ETSS surface
        "etsrate\n"              # ETSR via parameter "etsrate"
        "CONSTANT 5.0\n"         # ETSX extinction depth
    )
    ets = MfUsgEts.load(str(param_ets), ml, nper=1)
    # Parameter expanded on load -> NPETS reset to 0, ETSR filled with parval.
    assert ets.npets == 0
    assert np.allclose(ets.evtr[0].array, 5.0e-4)

    out = function_tmpdir / "expanded.ets"
    ets.fn_path = str(out)
    ets.write_file()
    content = out.read_text()
    assert "PARAMETER" not in content
    item2a = next(ln for ln in content.splitlines() if not ln.startswith("#"))
    assert item2a.split()[2] == "0"  # NPETS field expanded to 0


# --- Stage 3 Card 7: recipient-node U1DINT controls (INTERNAL/CONSTANT) -----

def test_mfusgdrt_recipient_constant_u1dint(function_tmpdir):
    """A spreading recipient list written as CONSTANT expands to that node."""
    drt_in = function_tmpdir / "const.drt"
    drt_in.write_text(
        "# drt CONSTANT recipients\n"
        "         1 0 0 0 RETURNFLOW\n 1 SP1\n"
        " 1  5.000000e+00  1.000000e+01  -3  7.000000e-01\n"
        "CONSTANT 5\n"  # 3 recipients, all node 5 (1-based) -> node 4 (0-based)
    )
    drt = MfUsgDrt.load(str(drt_in), _usgt_unstructured_model(function_tmpdir),
                        nper=1, ext_unit_dict={})
    assert drt.recipient_nodes[0][0] == [4, 4, 4]


def test_mfusg_recipient_external_u1dint_unsupported(function_tmpdir):
    """EXTERNAL/OPEN-CLOSE recipient U1DINT lists fail explicitly (rare; deferred).

    The Fortran U1DINT technically accepts these controls, but recipient lists
    are short inline blocks in practice, so they are documented unsupported.
    """
    drt_ext = function_tmpdir / "ext.drt"
    drt_ext.write_text(
        "# drt EXTERNAL recipients\n"
        "         1 0 0 0 RETURNFLOW\n 1 SP1\n"
        " 1  5.000000e+00  1.000000e+01  -2  7.000000e-01\n"
        "EXTERNAL 88\n"
    )
    with pytest.raises(NotImplementedError, match="recipient"):
        MfUsgDrt.load(str(drt_ext), _usgt_unstructured_model(function_tmpdir),
                      nper=1, ext_unit_dict={})


# --- Stage 3 Card 8: plain-checkmark hardening (PCB, EVT transport) --------

def test_mfusgpcb_authoring_roundtrip(function_tmpdir):
    """PCB (Prescribed Concentration Boundary) authored from scratch: 0-based
    nodes internal, 1-based in the file, species + concentration preserved."""
    from flopy.mfusg import MfUsgPcb

    ml = _usgt_unstructured_model(function_tmpdir)
    dtype = MfUsgPcb.get_default_dtype(structured=False)  # (node, iSpec, conc)
    spd = {0: np.array([(4, 1, 0.25), (9, 1, 0.5)], dtype=dtype).view(np.recarray)}
    pcb = MfUsgPcb(ml, stress_period_data=spd)
    pcb.fn_path = str(function_tmpdir / "created.pcb")
    pcb.write_file()
    text = Path(pcb.fn_path).read_text()
    assert "\n         5         1" in text  # node 4 -> 5 (1-based)
    assert "\n        10         1" in text  # node 9 -> 10

    ml2 = _usgt_unstructured_model(function_tmpdir)
    pcb2 = MfUsgPcb.load(pcb.fn_path, ml2, nper=1)
    assert list(pcb2.stress_period_data[0]["node"]) == [4, 9]
    assert np.isclose(pcb2.stress_period_data[0]["conc"][1], 0.5)


def test_mfusgevt_transport_etfactor_authoring(function_tmpdir):
    """EVT transport ET factor (IETFACTOR/ETFACTOR) authors when transport is on."""
    from flopy.mfusg import MfUsgEvt
    from flopy.modflow import ModflowDis

    ml = MfUsg(structured=True, model_ws=str(function_tmpdir))
    ModflowDis(ml, nlay=1, nrow=2, ncol=2, nper=1)
    ml.itrnsp = 1   # emulate active transport (BCT would set this)
    ml.mcomp = 1
    evt = MfUsgEvt(ml, nevtop=1, evtr=1.0e-4, ietfactor=1, etfactor=[2.5])
    evt.fn_path = str(function_tmpdir / "transport.evt")
    evt.write_file()
    content = Path(evt.fn_path).read_text()
    assert "2.5" in content  # ETFACTOR(MCOMP) record written
    assert evt.ietfactor == 1


def test_mfusgoc_atsa_authoring_roundtrip(function_tmpdir):
    """OC ATS adaptive time-stepping (ATSA) authors from scratch and round-trips."""
    from flopy.mfusg import MfUsgOc
    from flopy.modflow import ModflowDis

    ml = MfUsg(structured=True, model_ws=str(function_tmpdir))
    ModflowDis(ml, nlay=1, nrow=2, ncol=2, nper=1, nstp=1)
    MfUsgBas(ml, ibound=1, strt=1.0)
    oc = MfUsgOc(ml, atsa=1, stress_period_data={(0, 0): ["save head"]})
    oc.fn_path = str(function_tmpdir / "ats.oc")
    oc.write_file()
    assert "ATSA" in Path(oc.fn_path).read_text().splitlines()[1]

    ml2 = MfUsg(structured=True, model_ws=str(function_tmpdir))
    ModflowDis(ml2, nlay=1, nrow=2, ncol=2, nper=1, nstp=1)
    MfUsgBas(ml2, ibound=1, strt=1.0)
    oc2 = MfUsgOc.load(oc.fn_path, ml2, nper=1)
    assert oc2.atsa == 1
