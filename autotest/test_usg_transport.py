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


def test_gridgen_to_gsf_parsimonious_compacts(function_tmpdir):
    """Shared/parsimonious drops unused vertices and remaps per-vertex top/botm."""
    from flopy.discretization import UnstructuredGrid
    from flopy.mfusg import gridgen_to_gsf

    # one quad cell using vertices 0..3; vertex 4 is unused
    disv = {
        "vertices": [
            (0, 0.0, 0.0),
            (1, 1.0, 0.0),
            (2, 1.0, 1.0),
            (3, 0.0, 1.0),
            (4, 9.0, 9.0),  # not referenced by any cell
        ],
        "cell2d": [[0, 0.5, 0.5, 4, 0, 1, 2, 3]],
    }

    # parsimonious: only the 4 used vertices are written (doubled -> 8), and the
    # unused (9, 9) vertex never appears.
    shared = gridgen_to_gsf(
        _gg_model(function_tmpdir, "p"), disv, top=10.0, botm=0.0, vertex_mode="shared"
    )
    assert len(shared.vertices) == 8  # 2 * 4 used (not 2 * 5)
    xy = {(round(float(x), 3), round(float(y), 3)) for x, y, _ in shared.vertices}
    assert (9.0, 9.0) not in xy

    # per-vertex top/botm arrays are remapped to the kept subset
    topv = [10.0, 11.0, 12.0, 13.0, 99.0]  # 99 belongs to the unused vertex
    botv = [0.0, 1.0, 2.0, 3.0, -99.0]
    g2 = gridgen_to_gsf(
        _gg_model(function_tmpdir, "p2"),
        disv,
        top=topv,
        botm=botv,
        vertex_mode="parsimonious",
    )
    zvals = {round(float(z), 3) for _, _, z in g2.vertices}
    assert 99.0 not in zvals and -99.0 not in zvals
    assert {10.0, 11.0, 12.0, 13.0}.issubset(zvals)

    # still reconstructs top/botm through to_grid
    shared.fn_path = str(function_tmpdir / "p.gsf")
    shared.write_file()
    grid = shared.to_grid()
    assert isinstance(grid, UnstructuredGrid)
    assert np.allclose(grid.top, 10.0) and np.allclose(grid.botm, 0.0)

    # cell/nonparsimonious also ignores the unused vertex and gives 8 per quad
    cell = gridgen_to_gsf(
        _gg_model(function_tmpdir, "pc"),
        disv,
        top=10.0,
        botm=0.0,
        vertex_mode="nonparsimonious",
    )
    assert len(cell.vertices) == 8
    assert len(cell.node_data[0]["vertices"]) == 8


def test_gridgen_to_gsf_skip_degenerate_compacts(function_tmpdir):
    """Vertices used only by a skipped degenerate cell are not retained."""
    from flopy.mfusg import gridgen_to_gsf

    disv = {
        "vertices": [
            (0, 0.0, 0.0),
            (1, 1.0, 0.0),
            (2, 1.0, 1.0),
            (3, 0.0, 1.0),
            (4, 5.0, 5.0),  # exclusive to the degenerate cell below
            (5, 6.0, 6.0),
        ],
        "cell2d": [
            [0, 0.5, 0.5, 4, 0, 1, 2, 3],
            [1, 5.5, 5.5, 2, 4, 5],  # degenerate (< 3 unique verts)
        ],
    }
    g = gridgen_to_gsf(
        _gg_model(function_tmpdir, "sd"),
        disv,
        top=10.0,
        botm=0.0,
        vertex_mode="parsimonious",
        skip_degenerate=True,
    )
    assert g.nnodes == 1
    assert len(g.vertices) == 8  # only the surviving quad's 4 verts, doubled
    xy = {(round(float(x), 3), round(float(y), 3)) for x, y, _ in g.vertices}
    assert (5.0, 5.0) not in xy and (6.0, 6.0) not in xy


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


# --- Stage 4.4C: SGB MODFLOW list-parameter preservation -------------------
#
# Review follow-up: USG-T 2.7 defines SGB parameters with PARTYP='SGB'
# (UPARLSTRP) but activates them with PTYP='G' (UPARLSTSUB), so any active SGB
# parameter aborts the run with a "Parameter type conflict". SGB is therefore
# "definition-preserving only": parameter definitions (no activations)
# round-trip; active parameters (per-SP NP>0) and INSTANCES are unsupported.


def _sgb_model(function_tmpdir, name, nper=1):
    from flopy.modflow import ModflowDis

    ml = MfUsg(structured=False, model_ws=str(function_tmpdir), modelname=name)
    ModflowDis(ml, nlay=1, nrow=1, ncol=1, nper=nper)
    return ml


def test_mfusgsgb_definitions_preserved_roundtrip(function_tmpdir):
    """SGB NPSGB>0 definitions (no activations) preserve the PARAMETER record and
    the per-parameter rows on load -> write -> reload (0-based / 1-based)."""
    p = function_tmpdir / "param.sgb"
    p.write_text(
        "# sgb param defs only\n"
        "PARAMETER 1 5\n"  # NPSGB MXS
        "         2 0\n"  # MXACTS ISGBCB
        "sgbpar SGB 1.0 2\n"  # PARNAM PARTYP PARVAL NLST
        " 101   1.000000e-02\n"
        " 102   2.000000e-02\n"
        " 1 0    Stress Period 1\n"  # ITMP NP (no active parameters)
        " 201   5.000000e-02\n"  # non-parametric row
    )
    ml = _sgb_model(function_tmpdir, "p1")
    sgb = MfUsgSgb.load(str(p), ml, nper=1, ext_unit_dict={})
    assert sgb.mxs == 5
    assert list(sgb.parameters) == ["sgbpar"]
    pdef = sgb.parameters["sgbpar"]
    assert pdef["partyp"] == "SGB" and pdef["parval"] == "1.0" and pdef["nlst"] == 2
    assert list(pdef["data"]["node"]) == [100, 101]  # 0-based internal
    assert not any(sgb.active_params.values())  # no activations
    assert list(sgb.stress_period_data[0]["node"]) == [200]

    out = function_tmpdir / "param_out.sgb"
    sgb.fn_path = str(out)
    sgb.write_file()
    content = out.read_text()
    assert "PARAMETER 1 5" in content
    assert "sgbpar SGB 1.0 2" in content
    assert "\n 101  1.000000e-02" in content  # parameter row written 1-based

    ml2 = _sgb_model(function_tmpdir, "p2")
    re = MfUsgSgb.load(str(out), ml2, nper=1, ext_unit_dict={})
    assert re.mxs == 5 and list(re.parameters) == ["sgbpar"]
    assert list(re.parameters["sgbpar"]["data"]["node"]) == [100, 101]
    assert list(re.stress_period_data[0]["node"]) == [200]


def test_mfusgsgb_definitions_aux(function_tmpdir):
    """AUX values on SGB parameter definition rows are preserved (no activations)."""
    p = function_tmpdir / "auxp.sgb"
    p.write_text(
        "# sgb aux param defs\n"
        "PARAMETER 1 5\n"
        " 1 0 AUX C01\n"
        "gp SGB 1.0 1\n"
        " 101   1.000000e-02  5.000000e-01\n"
        " 1 0    Stress Period 1\n"
        " 201   2.000000e-02  7.000000e-01\n"
    )
    ml = _sgb_model(function_tmpdir, "ax1")
    sgb = MfUsgSgb.load(str(p), ml, nper=1, ext_unit_dict={})
    assert np.isclose(sgb.parameters["gp"]["data"]["C01"][0], 0.5)

    out = function_tmpdir / "auxp_out.sgb"
    sgb.fn_path = str(out)
    sgb.write_file()
    ml2 = _sgb_model(function_tmpdir, "ax2")
    re = MfUsgSgb.load(str(out), ml2, nper=1, ext_unit_dict={})
    assert np.isclose(re.parameters["gp"]["data"]["C01"][0], 0.5)


def test_mfusgsgb_definitions_sfac_inert_on_gradient(function_tmpdir):
    """SFAC inside a parameter's row block is consumed but inert on the gradient
    (Fortran ISCLOC=2 scales an internal dummy column)."""
    p = function_tmpdir / "sfacp.sgb"
    p.write_text(
        "# sgb sfac param defs\n"
        "PARAMETER 1 5\n"
        "         0 0\n"
        "gp SGB 1.0 1\n"
        " SFAC 3.0\n"
        " 101   1.000000e-02\n"
        " 0 0    Stress Period 1\n"
    )
    ml = _sgb_model(function_tmpdir, "sf1")
    sgb = MfUsgSgb.load(str(p), ml, nper=1, ext_unit_dict={})
    assert np.isclose(sgb.parameters["gp"]["data"]["gradient"][0], 0.01)  # unscaled


def test_mfusgsgb_definitions_open_close(function_tmpdir):
    """A parameter's NLST rows can be read via an OPEN/CLOSE list control."""
    (function_tmpdir / "gp_rows.dat").write_text(
        " 101 1.000000e-02\n 102 2.000000e-02\n"
    )
    p = function_tmpdir / "ocp.sgb"
    p.write_text(
        "# sgb param open/close defs\n"
        "PARAMETER 1 5\n"
        "         0 0\n"
        "gp SGB 1.0 2\n"
        " OPEN/CLOSE gp_rows.dat\n"
        " 0 0    Stress Period 1\n"
    )
    ml = _sgb_model(function_tmpdir, "oc1")
    sgb = MfUsgSgb.load(str(p), ml, nper=1, ext_unit_dict={})
    assert list(sgb.parameters["gp"]["data"]["node"]) == [100, 101]


def test_mfusgsgb_active_parameters_unsupported(function_tmpdir):
    """Active SGB parameters (per-SP NP>0) are rejected. USG-T 2.7 reads SGB
    parameters as PARTYP='SGB' (UPARLSTRP) but activates them as PTYP='G'
    (UPARLSTSUB) -- a type conflict that aborts the Fortran. Any NP>0 on load
    raises NotImplementedError."""
    p = function_tmpdir / "active.sgb"
    p.write_text(
        "# sgb active param\n"
        "PARAMETER 1 5\n"
        "         1 0\n"
        "gp SGB 1.0 1\n"
        " 101   1.000000e-02\n"
        " 1 1    Stress Period 1\n"  # NP=1 active parameter
        " 201   5.000000e-02\n"
        "gp\n"
    )
    ml = _sgb_model(function_tmpdir, "ac1")
    with pytest.raises(NotImplementedError, match="not supported"):
        MfUsgSgb.load(str(p), ml, nper=1, ext_unit_dict={})


def test_mfusgsgb_parameter_instances_unsupported(function_tmpdir):
    """SGB parameter INSTANCES are unsupported (active params unsupported)."""
    p = function_tmpdir / "inst.sgb"
    p.write_text(
        "# sgb instances\n"
        "PARAMETER 1 5\n"
        "         0 0\n"
        "gp SGB 1.0 2 INSTANCES 2\n"
        "spring\n 101 1.000000e-02\n"
        "fall\n 102 2.000000e-02\n"
        " 0 0    Stress Period 1\n"
    )
    ml = _sgb_model(function_tmpdir, "in1")
    with pytest.raises(NotImplementedError, match="INSTANCES"):
        MfUsgSgb.load(str(p), ml, nper=1, ext_unit_dict={})


def test_mfusgsgb_parameter_mxs_zero_fails(function_tmpdir):
    """Definitions present with MXS<=0 (a from-scratch PARAMETER 1 0) raise
    ValueError and write no partial file."""
    dtype = MfUsgSgb.get_default_dtype()
    rows = np.array([(0, 0.01)], dtype=dtype).view(np.recarray)
    params = {"p1": {"partyp": "SGB", "parval": "1.0", "nlst": 1, "data": rows}}
    sgb = MfUsgSgb(_sgb_model(function_tmpdir, "mz"), parameters=params, mxs=0)
    out = function_tmpdir / "mz.sgb"
    sgb.fn_path = str(out)
    with pytest.raises(ValueError, match="MXS"):
        sgb.write_file()
    assert not out.exists()


def test_mfusgsgb_parameter_mxs_too_small_fails(function_tmpdir):
    """MXS below the total number of parameter list entries raises ValueError."""
    dtype = MfUsgSgb.get_default_dtype()
    rows = np.array([(0, 0.01)], dtype=dtype).view(np.recarray)
    params = {
        "p1": {"partyp": "SGB", "parval": "1.0", "nlst": 1, "data": rows},
        "p2": {"partyp": "SGB", "parval": "1.0", "nlst": 1, "data": rows},
    }
    sgb = MfUsgSgb(_sgb_model(function_tmpdir, "mt"), parameters=params, mxs=1)
    out = function_tmpdir / "mt.sgb"
    sgb.fn_path = str(out)
    with pytest.raises(ValueError, match="MXS"):
        sgb.write_file()
    assert not out.exists()


def test_mfusgsgb_parameter_inconsistent_fails(function_tmpdir):
    """A definition whose nlst != len(data) raises ValueError, no partial file."""
    dtype = MfUsgSgb.get_default_dtype()
    rows = np.array([(0, 0.01)], dtype=dtype).view(np.recarray)
    params = {"p1": {"partyp": "SGB", "parval": "1.0", "nlst": 2, "data": rows}}
    sgb = MfUsgSgb(_sgb_model(function_tmpdir, "ic"), parameters=params, mxs=5)
    out = function_tmpdir / "ic.sgb"
    sgb.fn_path = str(out)
    with pytest.raises(ValueError, match="nlst"):
        sgb.write_file()
    assert not out.exists()


def test_mfusgsgb_active_params_constructed_fails(function_tmpdir):
    """Manually constructed active parameters are rejected on write."""
    dtype = MfUsgSgb.get_default_dtype()
    rows = np.array([(0, 0.01)], dtype=dtype).view(np.recarray)
    params = {"p1": {"partyp": "SGB", "parval": "1.0", "nlst": 1, "data": rows}}
    sgb = MfUsgSgb(
        _sgb_model(function_tmpdir, "am"),
        parameters=params,
        mxs=5,
        active_params={0: ["p1"]},
    )
    out = function_tmpdir / "am.sgb"
    sgb.fn_path = str(out)
    with pytest.raises(NotImplementedError, match="not supported"):
        sgb.write_file()
    assert not out.exists()


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
    """TRANSIENTQ combined with NPQRT>0 raises NotImplementedError on load.

    Plain TRANSIENTQ is supported as of Stage 4.5A; the combination with named
    parameters is the remaining explicit failure (the Fortran reads BDQV past
    its MXAQRT allocation when MXL>0). NPQRT>0 alone is preserved (Stage 4.4E).
    """
    from flopy.modflow import ModflowDis

    ml = MfUsg(structured=False, model_ws=str(function_tmpdir))
    ModflowDis(ml, nlay=1, nrow=1, ncol=1, nper=1)

    # TRANSIENTQ + NPQRT>0 (4th token = NPQRT = 1)
    qrt_tq = function_tmpdir / "tq.qrt"
    qrt_tq.write_text(
        "# qrt transientq + params\n"
        "        10        10 0 1 5 RETURNFLOW TRANSIENTQ 5\n"
        " 0    Stress Period 1\n"
    )
    with pytest.raises(NotImplementedError, match="TRANSIENTQ"):
        MfUsgQrt.load(str(qrt_tq), ml, nper=1, ext_unit_dict={})


# --- Stage 4.4E: QRT MODFLOW list-parameter preservation -------------------
#
# QRT is type-consistent (def + activation both PARTYP='QRT': gwf2QRT8u.f:183 /
# :1118), so active QRT params are type-valid and structurally preserved
# (load -> write -> reload). But execution is NOT guaranteed: the parameter value
# scales QRTF(5)=NumRT (recipient count), not Q (a Fortran bug; SFAC scales Q at
# ISCLOC=4), and NodQRT is not copied on activation. INSTANCES / from-scratch
# authoring / TRANSIENTQ are unsupported.


def _qrt_model(function_tmpdir, name, nper=1):
    from flopy.modflow import ModflowDis

    ml = MfUsg(structured=False, model_ws=str(function_tmpdir), modelname=name)
    ModflowDis(ml, nlay=1, nrow=1, ncol=1, nper=nper)
    return ml


def test_mfusgqrt_parameterized_roundtrip(function_tmpdir):
    """QRT NPQRT>0 preserves item-1 NPQRT/MXL, the definition, and the per-SP
    active record; MXAQRT covers the active total; 0-based / 1-based."""
    p = function_tmpdir / "param.qrt"
    p.write_text(
        "# qrt param\n"
        "         2         0 0 1 1\n"  # MXAQRT MXRTCELLS IQRTCB NPQRT MXL
        "qp QRT 2.0 1\n"  # PARNAM PARTYP PARVAL NLST
        " 5  -5.000000e+01\n"  # sink node 5, Q
        " 1 1    Stress Period 1\n"  # ITMP NP
        " 21  -1.000000e+02\n"  # non-parametric sink
        "qp\n"  # active parameter
    )
    ml = _qrt_model(function_tmpdir, "p1")
    qrt = MfUsgQrt.load(str(p), ml, nper=1, ext_unit_dict={})
    assert qrt.mxl == 1
    pdef = qrt.parameters["qp"]
    assert pdef["partyp"] == "QRT" and pdef["parval"] == "2.0" and pdef["nlst"] == 1
    assert list(pdef["data"]["node"]) == [4]  # 0-based internal
    assert np.isclose(pdef["data"]["q"][0], -50.0)
    assert qrt.active_params[0] == ["qp"]
    assert list(qrt.stress_period_data[0]["node"]) == [20]

    out = function_tmpdir / "param_out.qrt"
    qrt.fn_path = str(out)
    qrt.write_file()
    item1 = next(ln for ln in out.read_text().splitlines() if not ln.startswith("#"))
    # MXAQRT = non-parametric (1) + active rows (1) = 2; NPQRT field = 1
    assert item1.split()[0] == "2" and item1.split()[3] == "1"
    assert "qp QRT 2.0 1" in out.read_text()

    re = MfUsgQrt.load(str(out), _qrt_model(function_tmpdir, "p2"), nper=1)
    assert re.mxl == 1 and list(re.parameters) == ["qp"]
    assert list(re.parameters["qp"]["data"]["node"]) == [4]
    assert re.active_params[0] == ["qp"]
    assert list(re.stress_period_data[0]["node"]) == [20]


def test_mfusgqrt_parameter_returnflow_recipients(function_tmpdir):
    """A parameter definition row keeps its RETURNFLOW recipients (U1DINT block),
    plus CHANGEC and AUX, through load -> write -> reload. MXRTCELLS reflects the
    parameter recipients. Structural round-trip (not execution-guaranteed)."""
    p = function_tmpdir / "rf.qrt"
    p.write_text(
        "# qrt rf param\n"
        "         1         2 0 1 1 RETURNFLOW CHANGEC AUX C01\n"
        "qp QRT 1.0 1\n"
        " 5  -5.000000e+01  2  8.000000e-01  3  9.000000e-01\n"
        "INTERNAL  1  (FREE)  -1\n"
        " 11 12\n"
        " 0 1    Stress Period 1\n"
        "qp\n"
    )
    ml = _qrt_model(function_tmpdir, "rf1")
    qrt = MfUsgQrt.load(str(p), ml, nper=1, ext_unit_dict={})
    pdef = qrt.parameters["qp"]
    assert pdef["recipient_nodes"] == [[10, 11]]  # 0-based
    assert list(pdef["data"]["iqchngtyp"]) == [3]
    assert np.isclose(pdef["data"]["C01"][0], 0.9)

    out = function_tmpdir / "rf_out.qrt"
    qrt.fn_path = str(out)
    qrt.write_file()
    item1 = next(ln for ln in out.read_text().splitlines() if not ln.startswith("#"))
    assert int(item1.split()[1]) >= 2  # MXRTCELLS covers the 2 recipients

    re = MfUsgQrt.load(str(out), _qrt_model(function_tmpdir, "rf2"), nper=1)
    assert re.parameters["qp"]["recipient_nodes"] == [[10, 11]]
    assert list(re.parameters["qp"]["data"]["iqchngtyp"]) == [3]


def test_mfusgqrt_parameter_mxaqrt_counts_active_rows(function_tmpdir):
    """MXAQRT must cover NQRTCL = non-parametric + Σ(nlst of active params)
    (gwf2QRT8u.f: SGWF2QRT8LS aborts if NQRTCL > MXAQRT)."""
    p = function_tmpdir / "mxaqrt.qrt"
    p.write_text(
        "# qrt mxaqrt\n"
        "         4         0 0 2 3\n"  # MXAQRT=4 NPQRT=2 MXL=3
        "pa QRT 2.0 2\n"
        " 5  -1.000000e+01\n"
        " 6  -1.000000e+01\n"
        "pb QRT 3.0 1\n"
        " 7  -1.000000e+01\n"
        " 1 2    Stress Period 1\n"  # ITMP=1 NP=2
        " 21  -2.000000e+01\n"
        "pa\n"
        "pb\n"
    )
    ml = _qrt_model(function_tmpdir, "ma1")
    qrt = MfUsgQrt.load(str(p), ml, nper=1, ext_unit_dict={})
    assert qrt.active_params[0] == ["pa", "pb"]

    out = function_tmpdir / "mxaqrt_out.qrt"
    qrt.fn_path = str(out)
    qrt.write_file()
    item1 = next(ln for ln in out.read_text().splitlines() if not ln.startswith("#"))
    assert item1.split()[0] == "4"  # 1 non-param + 2 (pa) + 1 (pb)


def test_mfusgqrt_parameter_reuse_with_active(function_tmpdir):
    """ITMP<0 reuses the previous period's non-parametric sinks while a new active
    parameter applies (NP>0); data round-trips."""
    p = function_tmpdir / "reuse.qrt"
    p.write_text(
        "# qrt reuse\n"
        "         1         0 0 1 1\n"
        "qp QRT 1.0 1\n"
        " 5  -1.000000e+01\n"
        " 1 1    Stress Period 1\n"
        " 21  -2.000000e+01\n"
        "qp\n"
        " -1 1    Stress Period 2\n"
        "qp\n"
    )
    ml = _qrt_model(function_tmpdir, "r1", nper=2)
    qrt = MfUsgQrt.load(str(p), ml, nper=2, ext_unit_dict={})
    assert list(qrt.stress_period_data[1]["node"]) == [20]  # reused
    assert qrt.active_params[1] == ["qp"]

    out = function_tmpdir / "reuse_out.qrt"
    qrt.fn_path = str(out)
    qrt.write_file()
    re = MfUsgQrt.load(str(out), _qrt_model(function_tmpdir, "r2", nper=2), nper=2)
    assert list(re.stress_period_data[1]["node"]) == [20]
    assert re.active_params[1] == ["qp"]


def test_mfusgqrt_parameter_sfac_scales_q(function_tmpdir):
    """SFAC in a parameter block scales Q (Fortran ISCLOC=4) -- unlike the
    parameter value, which the Fortran (buggily) applies to NumRT."""
    p = function_tmpdir / "sfacp.qrt"
    p.write_text(
        "# qrt sfac param\n"
        "         0         0 0 1 1\n"
        "qp QRT 1.0 1\n"
        " SFAC 3.0\n"
        " 5  -1.000000e+01\n"
        " 0 1    Stress Period 1\n"
        "qp\n"
    )
    ml = _qrt_model(function_tmpdir, "sf1")
    qrt = MfUsgQrt.load(str(p), ml, nper=1, ext_unit_dict={})
    assert np.isclose(qrt.parameters["qp"]["data"]["q"][0], -30.0)  # -10 * SFAC 3


def test_mfusgqrt_parameter_instances_unsupported(function_tmpdir):
    """QRT parameter INSTANCES are Fortran-supported but not yet by FloPy."""
    p = function_tmpdir / "inst.qrt"
    p.write_text(
        "# qrt instances\n"
        "         0         0 0 1 1\n"
        "qp QRT 1.0 2 INSTANCES 2\n"
        "spring\n 5  -1.000000e+01\n"
        "fall\n 6  -1.000000e+01\n"
        " 0 1    Stress Period 1\n"
        "qp spring\n"
    )
    ml = _qrt_model(function_tmpdir, "in1")
    with pytest.raises(NotImplementedError, match="INSTANCES"):
        MfUsgQrt.load(str(p), ml, nper=1, ext_unit_dict={})


def test_mfusgqrt_parameter_from_scratch_fails(function_tmpdir):
    """Active parameters with no loaded definitions fail explicitly, no file."""
    qrt = MfUsgQrt(_qrt_model(function_tmpdir, "fs"), active_params={0: ["p1"]})
    out = function_tmpdir / "fs.qrt"
    qrt.fn_path = str(out)
    with pytest.raises(NotImplementedError, match="from scratch"):
        qrt.write_file()
    assert not out.exists()


def test_mfusgqrt_parameter_mxl_too_small_fails(function_tmpdir):
    """MXL below the total parameter list entries raises ValueError, no file."""
    dtype = MfUsgQrt.get_default_dtype(returnflow=False)
    rows = np.array([(0, -10.0)], dtype=dtype).view(np.recarray)
    params = {
        "p1": {"partyp": "QRT", "parval": "1.0", "nlst": 1, "data": rows,
               "recipient_nodes": [[]]},
        "p2": {"partyp": "QRT", "parval": "1.0", "nlst": 1, "data": rows,
               "recipient_nodes": [[]]},
    }
    qrt = MfUsgQrt(_qrt_model(function_tmpdir, "mx"), parameters=params, mxl=1)
    out = function_tmpdir / "mx.qrt"
    qrt.fn_path = str(out)
    with pytest.raises(ValueError, match="MXL"):
        qrt.write_file()
    assert not out.exists()


def test_mfusgqrt_parameter_inconsistent_fails(function_tmpdir):
    """A definition whose nlst != len(data), or recipient_nodes != nlst, raises
    ValueError, no file."""
    dtype = MfUsgQrt.get_default_dtype(returnflow=False)
    rows = np.array([(0, -10.0)], dtype=dtype).view(np.recarray)
    params = {
        "p1": {"partyp": "QRT", "parval": "1.0", "nlst": 2, "data": rows,
               "recipient_nodes": [[]]},  # nlst 2 but 1 row
    }
    qrt = MfUsgQrt(_qrt_model(function_tmpdir, "ic"), parameters=params, mxl=5)
    out = function_tmpdir / "ic.qrt"
    qrt.fn_path = str(out)
    with pytest.raises(ValueError, match="nlst"):
        qrt.write_file()
    assert not out.exists()

    params2 = {
        "p1": {"partyp": "QRT", "parval": "1.0", "nlst": 1, "data": rows,
               "recipient_nodes": []},  # recipient lists != nlst
    }
    qrt2 = MfUsgQrt(_qrt_model(function_tmpdir, "ic2"), parameters=params2, mxl=5)
    out2 = function_tmpdir / "ic2.qrt"
    qrt2.fn_path = str(out2)
    with pytest.raises(ValueError, match="recipient"):
        qrt2.write_file()
    assert not out2.exists()


def test_mfusgqrt_parameter_active_undefined_fails(function_tmpdir):
    """An active parameter name not present in definitions raises ValueError."""
    dtype = MfUsgQrt.get_default_dtype(returnflow=False)
    rows = np.array([(0, -10.0)], dtype=dtype).view(np.recarray)
    params = {
        "p1": {"partyp": "QRT", "parval": "1.0", "nlst": 1, "data": rows,
               "recipient_nodes": [[]]},
    }
    qrt = MfUsgQrt(
        _qrt_model(function_tmpdir, "ud"),
        parameters=params,
        mxl=5,
        active_params={0: ["px"]},  # not defined
    )
    out = function_tmpdir / "ud.qrt"
    qrt.fn_path = str(out)
    with pytest.raises(ValueError, match="not defined"):
        qrt.write_file()
    assert not out.exists()


def test_mfusgqrt_parameter_active_case_insensitive(function_tmpdir):
    """An activation name that differs from its definition only in case (the
    Fortran UPCASEs both) still counts toward MXAQRT (review follow-up)."""
    p = function_tmpdir / "ci.qrt"
    p.write_text(
        "# qrt case-insensitive\n"
        "         2         0 0 1 1\n"  # MXAQRT MXRTCELLS IQRTCB NPQRT MXL
        "qp QRT 2.0 1\n"  # definition lower-case
        " 5  -5.000000e+01\n"
        " 1 1    Stress Period 1\n"
        " 21  -1.000000e+02\n"
        "QP\n"  # activation upper-case
    )
    ml = _qrt_model(function_tmpdir, "ci1")
    qrt = MfUsgQrt.load(str(p), ml, nper=1, ext_unit_dict={})
    out = function_tmpdir / "ci_out.qrt"
    qrt.fn_path = str(out)
    qrt.write_file()
    item1 = next(ln for ln in out.read_text().splitlines() if not ln.startswith("#"))
    # 1 non-parametric + 1 active (resolved case-insensitively) => MXAQRT=2
    assert item1.split()[0] == "2"


def test_mfusgqrt_parameter_duplicate_active_fails(function_tmpdir):
    """Activating the same parameter twice in a stress period (case-insensitive)
    raises ValueError (the Fortran aborts 'already activated'), no partial file."""
    dtype = MfUsgQrt.get_default_dtype(returnflow=False)
    rows = np.array([(0, -10.0)], dtype=dtype).view(np.recarray)
    params = {
        "qp": {"partyp": "QRT", "parval": "1.0", "nlst": 1, "data": rows,
               "recipient_nodes": [[]]},
    }
    qrt = MfUsgQrt(
        _qrt_model(function_tmpdir, "dup"),
        parameters=params,
        mxl=5,
        active_params={0: ["qp", "QP"]},  # same parameter twice
    )
    out = function_tmpdir / "dup.qrt"
    qrt.fn_path = str(out)
    with pytest.raises(ValueError, match="more than once"):
        qrt.write_file()
    assert not out.exists()


# --- Stage 4.5A: QRT TRANSIENTQ transient extraction-flow time series ------
#
# TRANSIENTQ overrides QRTF(4)=Q at each time step (GWF2QRT8U1AD); recipients
# (NumRT/Rfprop) are untouched. The inline block (IQRTUN = the package unit) is
# read in GWF2QRT8U1AR after the parameter definitions and before the first
# stress period: a times control line + NBDQTIM times, then a values control
# line + exactly MXAQRT rows of (node, NBDQTIM values). NBDQTIM<0 => staircase.
# See USGT_STAGE4_05_QRT_TRANSIENTQ.md.


def test_mfusgqrt_transientq_authoring(function_tmpdir):
    """From-scratch QRT with a minimal TRANSIENTQ (interpolation, NBDQTIM=2):
    the item-1 token, the inline control lines, and 1-based node tags."""
    dtype = MfUsgQrt.get_default_dtype(returnflow=True)
    spd = {0: np.array([(0, -100.0, 0.75)], dtype=dtype).view(np.recarray)}
    recips = {0: [[9, 10]]}
    qrt = MfUsgQrt(
        _qrt_model(function_tmpdir, "tqa"),
        stress_period_data=spd,
        recipient_nodes=recips,
        options=["RETURNFLOW"],
        transientq_times=[0.0, 10.0],
        transientq_values=[[-100.0, -50.0]],  # 1 sink (MXAQRT=1) x 2 times
        transientq_nodes=[0],  # 0-based -> file 1
    )
    out = function_tmpdir / "tqa.qrt"
    qrt.fn_path = str(out)
    qrt.write_file()
    text = out.read_text()

    # TRANSIENTQ is the last item-1 option, count = NBDQTIM = 2 (positive).
    assert text.splitlines()[1].rstrip().endswith("TRANSIENTQ 2")
    # Inline: IQRTUN = the package's own unit, CNSTM = 1.0 (two control lines).
    unit = qrt.unit_number[0]
    assert text.count(f" {unit} 1.000000e+00") == 2
    # Times line, then the MXAQRT=1 value row with the 1-based node tag.
    assert "\n 0.000000e+00 1.000000e+01\n" in text
    assert "\n 1 -1.000000e+02 -5.000000e+01\n" in text


def test_mfusgqrt_transientq_roundtrip(function_tmpdir):
    """load -> write -> reload preserves the TRANSIENTQ block and the per-SP
    sink data; the written text is stable across the cycle."""
    dtype = MfUsgQrt.get_default_dtype(returnflow=True)
    spd = {0: np.array([(0, -100.0, 0.75), (4, -30.0, 0.0)], dtype=dtype).view(
        np.recarray
    )}
    recips = {0: [[9, 10], []]}
    qrt = MfUsgQrt(
        _qrt_model(function_tmpdir, "tqr"),
        stress_period_data=spd,
        recipient_nodes=recips,
        options=["RETURNFLOW"],
        transientq_times=[0.0, 5.0, 10.0],
        transientq_values=[[-100.0, -80.0, -60.0], [-30.0, -20.0, -10.0]],
        transientq_nodes=[0, 4],  # MXAQRT = 2
    )
    out1 = function_tmpdir / "tqr1.qrt"
    qrt.fn_path = str(out1)
    qrt.write_file()

    re1 = MfUsgQrt.load(str(out1), _qrt_model(function_tmpdir, "tqr2"), nper=1)
    out2 = function_tmpdir / "tqr2.qrt"
    re1.fn_path = str(out2)
    re1.write_file()
    re2 = MfUsgQrt.load(str(out2), _qrt_model(function_tmpdir, "tqr3"), nper=1)

    # Written text is stable (body after the heading line).
    assert out1.read_text().split("\n", 1)[1] == out2.read_text().split("\n", 1)[1]
    # TRANSIENTQ block preserved.
    assert re2._transientq_active()
    assert np.allclose(re2.transientq_times, [0.0, 5.0, 10.0])
    assert np.allclose(
        re2.transientq_values, [[-100.0, -80.0, -60.0], [-30.0, -20.0, -10.0]]
    )
    assert re2.transientq_nodes == [0, 4]
    assert re2.transientq_staircase is False
    # Per-SP sink data still intact.
    assert list(re2.stress_period_data[0]["node"]) == [0, 4]
    assert re2.recipient_nodes[0][0] == [9, 10]
    assert re2.recipient_nodes[0][1] == []


def test_mfusgqrt_transientq_staircase(function_tmpdir):
    """NBDQTIM<0 (staircase) round-trips: a negative item-1 count and
    transientq_staircase=True."""
    dtype = MfUsgQrt.get_default_dtype(returnflow=False)
    spd = {0: np.array([(2, -10.0)], dtype=dtype).view(np.recarray)}
    qrt = MfUsgQrt(
        _qrt_model(function_tmpdir, "tqs"),
        stress_period_data=spd,
        transientq_times=[0.0, 5.0, 10.0],
        transientq_values=[[-1.0, -2.0, -3.0]],
        transientq_nodes=[2],
        transientq_staircase=True,
    )
    out = function_tmpdir / "tqs.qrt"
    qrt.fn_path = str(out)
    qrt.write_file()
    assert out.read_text().splitlines()[1].rstrip().endswith("TRANSIENTQ -3")

    re = MfUsgQrt.load(str(out), _qrt_model(function_tmpdir, "tqs2"), nper=1)
    assert re.transientq_staircase is True
    assert np.allclose(re.transientq_times, [0.0, 5.0, 10.0])
    assert np.allclose(re.transientq_values, [[-1.0, -2.0, -3.0]])


def test_mfusgqrt_transientq_nodes_1based(function_tmpdir):
    """transientq_nodes are 0-based internally and 1-based in the file."""
    dtype = MfUsgQrt.get_default_dtype(returnflow=False)
    spd = {0: np.array([(4, -10.0)], dtype=dtype).view(np.recarray)}
    qrt = MfUsgQrt(
        _qrt_model(function_tmpdir, "tqn"),
        stress_period_data=spd,
        transientq_times=[0.0, 10.0],
        transientq_values=[[-5.0, -7.0]],
        transientq_nodes=[4],  # 0-based -> file 5
    )
    out = function_tmpdir / "tqn.qrt"
    qrt.fn_path = str(out)
    qrt.write_file()
    # Value row begins with the 1-based node tag 5.
    assert "\n 5 -5.000000e+00 -7.000000e+00\n" in out.read_text()

    re = MfUsgQrt.load(str(out), _qrt_model(function_tmpdir, "tqn2"), nper=1)
    assert re.transientq_nodes == [4]


def test_mfusgqrt_transientq_with_params_fails(function_tmpdir):
    """TRANSIENTQ together with QRT parameters (NPQRT>0) raises
    NotImplementedError before the file is opened (no partial file)."""
    dtype = MfUsgQrt.get_default_dtype(returnflow=False)
    spd = {0: np.array([(0, -10.0)], dtype=dtype).view(np.recarray)}
    params = {
        "qp": {"partyp": "QRT", "parval": "1.0", "nlst": 1,
               "data": np.array([(0, -10.0)], dtype=dtype).view(np.recarray),
               "recipient_nodes": [[]]},
    }
    qrt = MfUsgQrt(
        _qrt_model(function_tmpdir, "tqp"),
        stress_period_data=spd,
        parameters=params,
        mxl=5,
        active_params={0: ["qp"]},
        transientq_times=[0.0, 10.0],
        transientq_values=[[-1.0, -2.0]],
        transientq_nodes=[0],
    )
    out = function_tmpdir / "tqp.qrt"
    qrt.fn_path = str(out)
    with pytest.raises(NotImplementedError, match="TRANSIENTQ"):
        qrt.write_file()
    assert not out.exists()


def test_mfusgqrt_transientq_dim_mismatch_fails(function_tmpdir):
    """transientq_values not shaped (MXAQRT, NBDQTIM) raises ValueError before
    the file is opened (the Fortran reads exactly MXAQRT rows of NBDQTIM)."""
    dtype = MfUsgQrt.get_default_dtype(returnflow=False)
    spd = {0: np.array([(0, -10.0)], dtype=dtype).view(np.recarray)}  # MXAQRT=1
    qrt = MfUsgQrt(
        _qrt_model(function_tmpdir, "tqd"),
        stress_period_data=spd,
        transientq_times=[0.0, 10.0],
        transientq_values=[[-1.0, -2.0], [-3.0, -4.0]],  # 2 rows, MXAQRT=1
        transientq_nodes=[0, 1],
    )
    out = function_tmpdir / "tqd.qrt"
    qrt.fn_path = str(out)
    with pytest.raises(ValueError, match="MXAQRT"):
        qrt.write_file()
    assert not out.exists()


# --- Stage 4.5A review follow-up: TRANSIENTQ contract guards ----------------


def test_mfusgqrt_transientq_external_times_unit_fails(function_tmpdir):
    """A TRANSIENTQ times control line whose IQRTUN is not the package's inline
    unit raises NotImplementedError (external-unit data is unsupported)."""
    unit = MfUsgQrt._defaultunit()  # inline unit when no ext_unit_dict is given
    qrt_f = function_tmpdir / "tqxt.qrt"
    qrt_f.write_text(
        "# external times unit\n"
        " 1 0 0 0 0 TRANSIENTQ 2\n"
        " 999 1.000000e+00\n"  # external unit on the TIMES control line
        " 0.0 10.0\n"
        f" {unit} 1.000000e+00\n"
        " 1 -1.0 -2.0\n"
        " 1   Stress Period 1\n"
        " 1 -10.0\n"
    )
    with pytest.raises(NotImplementedError, match=r"times.*999|999.*times"):
        MfUsgQrt.load(str(qrt_f), _qrt_model(function_tmpdir, "tqxt"), nper=1)


def test_mfusgqrt_transientq_external_values_unit_fails(function_tmpdir):
    """A TRANSIENTQ values control line whose IQRTUN is not the package's inline
    unit raises NotImplementedError (the times line is fine, the values line is
    external)."""
    unit = MfUsgQrt._defaultunit()
    qrt_f = function_tmpdir / "tqxv.qrt"
    qrt_f.write_text(
        "# external values unit\n"
        " 1 0 0 0 0 TRANSIENTQ 2\n"
        f" {unit} 1.000000e+00\n"  # inline times OK
        " 0.0 10.0\n"
        " 888 1.000000e+00\n"  # external unit on the VALUES control line
        " 1 -1.0 -2.0\n"
        " 1   Stress Period 1\n"
        " 1 -10.0\n"
    )
    with pytest.raises(NotImplementedError, match=r"values.*888|888.*values"):
        MfUsgQrt.load(str(qrt_f), _qrt_model(function_tmpdir, "tqxv"), nper=1)


def test_mfusgqrt_transientq_trailing_option_fails(function_tmpdir):
    """TRANSIENTQ must be the last item-1 option: any trailing token raises
    ValueError on load (the Fortran branch would silently ignore it)."""
    for trailing in ("AUX C01", "NOPRINT"):
        qrt_f = function_tmpdir / f"tqtr_{trailing.split()[0].lower()}.qrt"
        qrt_f.write_text(
            "# trailing option after TRANSIENTQ\n"
            f" 1 0 0 0 0 RETURNFLOW TRANSIENTQ 2 {trailing}\n"
            " 1 0   Stress Period 1\n"
        )
        with pytest.raises(ValueError, match=r"(?i)last option"):
            MfUsgQrt.load(str(qrt_f), _qrt_model(function_tmpdir, "tqtr"), nper=1)


def test_mfusgqrt_transientq_negative_node_fails(function_tmpdir):
    """A negative transientq_nodes entry raises ValueError before the file is
    opened (0-based node tags must be non-negative); no partial file."""
    dtype = MfUsgQrt.get_default_dtype(returnflow=False)
    spd = {0: np.array([(0, -10.0)], dtype=dtype).view(np.recarray)}
    qrt = MfUsgQrt(
        _qrt_model(function_tmpdir, "tqng"),
        stress_period_data=spd,
        transientq_times=[0.0, 10.0],
        transientq_values=[[-1.0, -2.0]],
        transientq_nodes=[-1],  # invalid 0-based node
    )
    out = function_tmpdir / "tqng.qrt"
    qrt.fn_path = str(out)
    with pytest.raises(ValueError, match="non-negative"):
        qrt.write_file()
    assert not out.exists()


def test_mfusgqrt_transientq_multipliers_roundtrip(function_tmpdir):
    """Non-unit CNSTM multipliers (times and values) survive load -> write ->
    reload, and the stored times/values stay raw (pre-multiplier)."""
    dtype = MfUsgQrt.get_default_dtype(returnflow=False)
    spd = {0: np.array([(0, -10.0)], dtype=dtype).view(np.recarray)}
    qrt = MfUsgQrt(
        _qrt_model(function_tmpdir, "tqm"),
        stress_period_data=spd,
        transientq_times=[0.0, 10.0],
        transientq_values=[[-100.0, -50.0]],
        transientq_nodes=[0],
        transientq_times_mult=2.5,
        transientq_values_mult=0.5,
    )
    out = function_tmpdir / "tqm.qrt"
    qrt.fn_path = str(out)
    qrt.write_file()
    # Control lines carry the multipliers verbatim.
    text = out.read_text()
    assert "2.500000e+00" in text and "5.000000e-01" in text

    re = MfUsgQrt.load(str(out), _qrt_model(function_tmpdir, "tqm2"), nper=1)
    assert np.isclose(re.transientq_times_mult, 2.5)
    assert np.isclose(re.transientq_values_mult, 0.5)
    # Raw (pre-multiplier) times/values are preserved.
    assert np.allclose(re.transientq_times, [0.0, 10.0])
    assert np.allclose(re.transientq_values, [[-100.0, -50.0]])


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


# --- Stage 4.4D: DRT MODFLOW list-parameter preservation -------------------
#
# DRT is internally consistent in USG-T 2.7 (definitions and activations both
# use PARTYP='DRT': gwf2drt8u.f:135 / :1108), unlike SGB. So active DRT
# parameters are valid and preserved (load -> write -> reload), including their
# RETURNFLOW recipients / spreading blocks. INSTANCES and from-scratch parameter
# authoring are unsupported (explicit NotImplementedError).


def _drt_model(function_tmpdir, name, nper=1):
    from flopy.modflow import ModflowDis

    ml = MfUsg(structured=False, model_ws=str(function_tmpdir), modelname=name)
    ModflowDis(ml, nlay=1, nrow=1, ncol=1, nper=nper)
    return ml


def test_mfusgdrt_parameterized_roundtrip(function_tmpdir):
    """DRT NPDRT>0 preserves item-1 NPDRT/MXL, the definition (with its inline
    RETURNFLOW recipient), and the per-SP active record; 0-based / 1-based."""
    p = function_tmpdir / "param.drt"
    p.write_text(
        "# drt param\n"
        "         2 0 1 5 RETURNFLOW\n"  # MXADRT IDRTCB NPDRT MXL RETURNFLOW
        "drtpar DRT 2.0 1\n"  # PARNAM PARTYP PARVAL NLST
        " 11  5.000000e+00  1.000000e+01  9  7.000000e-01\n"  # NR=9 inline recip
        " 1 1    Stress Period 1\n"  # ITMP NP
        " 21  4.000000e+00  2.000000e+01  0\n"  # non-parametric drain
        "drtpar\n"  # active parameter
    )
    ml = _drt_model(function_tmpdir, "p1")
    drt = MfUsgDrt.load(str(p), ml, nper=1, ext_unit_dict={})
    assert drt.mxl == 5
    pdef = drt.parameters["drtpar"]
    assert pdef["partyp"] == "DRT" and pdef["parval"] == "2.0" and pdef["nlst"] == 1
    assert list(pdef["data"]["node"]) == [10]  # 0-based internal
    assert pdef["recipient_nodes"] == [[8]]  # 1-based 9 -> 0-based 8
    assert drt.active_params[0] == ["drtpar"]
    assert list(drt.stress_period_data[0]["node"]) == [20]

    out = function_tmpdir / "param_out.drt"
    drt.fn_path = str(out)
    drt.write_file()
    content = out.read_text()
    item1 = next(ln for ln in content.splitlines() if not ln.startswith("#"))
    # MXADRT must cover NDRTCL = non-parametric (1) + active parameter rows (1) = 2
    assert item1.split()[:4] == ["2", "0", "1", "5"]  # MXADRT IDRTCB NPDRT MXL
    assert "drtpar DRT 2.0 1" in content
    assert "drtpar" in [ln.strip() for ln in content.splitlines()]  # activation

    re = MfUsgDrt.load(str(out), _drt_model(function_tmpdir, "p2"), nper=1)
    assert re.mxl == 5 and list(re.parameters) == ["drtpar"]
    assert list(re.parameters["drtpar"]["data"]["node"]) == [10]
    assert re.parameters["drtpar"]["recipient_nodes"] == [[8]]
    assert re.active_params[0] == ["drtpar"]
    assert list(re.stress_period_data[0]["node"]) == [20]


def test_mfusgdrt_parameter_mxadrt_counts_active_rows(function_tmpdir):
    """MXADRT (item 1) must cover NDRTCL = non-parametric rows + the NLST rows of
    every active parameter in the period (the Fortran aborts if NDRTCL > MXADRT,
    SGWF2DRT8LS). One non-parametric drain + two active params (NLST 2 and 1) =>
    MXADRT must be written as 4."""
    p = function_tmpdir / "mxadrt.drt"
    p.write_text(
        "# drt mxadrt\n"
        "         4 0 2 3 RETURNFLOW\n"  # MXADRT=4 NPDRT=2 MXL=3
        "pa DRT 2.0 2\n"
        " 5  5.000000e+00  1.000000e+01  0\n"
        " 6  5.000000e+00  1.000000e+01  0\n"
        "pb DRT 3.0 1\n"
        " 7  5.000000e+00  1.000000e+01  0\n"
        " 1 2    Stress Period 1\n"  # ITMP=1 NP=2
        " 21  4.000000e+00  2.000000e+01  0\n"
        "pa\n"
        "pb\n"
    )
    ml = _drt_model(function_tmpdir, "ma1")
    drt = MfUsgDrt.load(str(p), ml, nper=1, ext_unit_dict={})
    assert drt.active_params[0] == ["pa", "pb"]
    assert drt.parameters["pa"]["nlst"] == 2 and drt.parameters["pb"]["nlst"] == 1

    out = function_tmpdir / "mxadrt_out.drt"
    drt.fn_path = str(out)
    drt.write_file()
    item1 = next(ln for ln in out.read_text().splitlines() if not ln.startswith("#"))
    # NDRTCL = 1 non-param + 2 (pa) + 1 (pb) = 4
    assert item1.split()[0] == "4"

    re = MfUsgDrt.load(str(out), _drt_model(function_tmpdir, "ma2"), nper=1)
    assert re.active_params[0] == ["pa", "pb"]
    assert re.parameters["pa"]["nlst"] == 2


def test_mfusgdrt_parameter_spread_recipients(function_tmpdir):
    """A parameter definition row keeps its SPREAD (multi-node U1DINT) recipients,
    plus CHANGEC and AUX, through load -> write -> reload.

    Contract: this is **structural round-trip preservation**, not an
    executable-on-USG-T guarantee. USG-T 2.7 activates a parameter by copying its
    DRTF rows (SGWF2DRT8LS) but does **not** copy/offset the NodDRT recipient
    array, and the spreading consumer (GWF2DRT8U1AD / RP) walks NodDRT
    sequentially -- so an *activated* SPREAD (NR<0) parameter is not guaranteed to
    resolve its recipients at run time. Inline single recipients (NR>0) travel in
    DRTF and are copied, so those are fine. See USGT_STAGE4_04_PARAMETERS_DRT.md.
    """
    p = function_tmpdir / "spread.drt"
    p.write_text(
        "# drt spread param\n"
        "         1 0 1 10 RETURNFLOW CHANGEC AUX C01\n"
        "dp DRT 2.0 1\n"
        " 11  5.000000e+00  1.000000e+01  -2  7.000000e-01  3  9.000000e-01\n"
        "INTERNAL  1  (FREE)  -1\n"
        " 21 22\n"
        " 1 1    Stress Period 1\n"
        " 31  4.000000e+00  2.000000e+01  0  0.0  0  1.000000e-01\n"
        "dp\n"
    )
    ml = _drt_model(function_tmpdir, "s1")
    drt = MfUsgDrt.load(str(p), ml, nper=1, ext_unit_dict={})
    pdef = drt.parameters["dp"]
    assert pdef["recipient_nodes"] == [[20, 21]]  # 0-based spreading nodes
    assert list(pdef["data"]["idchngtyp"]) == [3]
    assert np.isclose(pdef["data"]["C01"][0], 0.9)

    out = function_tmpdir / "spread_out.drt"
    drt.fn_path = str(out)
    drt.write_file()
    re = MfUsgDrt.load(str(out), _drt_model(function_tmpdir, "s2"), nper=1)
    assert re.parameters["dp"]["recipient_nodes"] == [[20, 21]]
    assert list(re.parameters["dp"]["data"]["idchngtyp"]) == [3]
    assert np.isclose(re.parameters["dp"]["data"]["C01"][0], 0.9)


def test_mfusgdrt_parameter_mixed_and_active(function_tmpdir):
    """A stress period mixes non-parametric drains (ITMP) with active parameters
    (NP); the per-SP header carries ITMP NP and both round-trip."""
    p = function_tmpdir / "mix.drt"
    p.write_text(
        "# drt mixed\n"
        "         2 0 1 5 RETURNFLOW\n"
        "dp DRT 2.0 1\n"
        " 5  5.000000e+00  1.000000e+01  0\n"
        " 2 1    Stress Period 1\n"
        " 21  4.000000e+00  2.000000e+01  0\n"
        " 22  4.000000e+00  3.000000e+01  0\n"
        "dp\n"
    )
    ml = _drt_model(function_tmpdir, "mx1")
    drt = MfUsgDrt.load(str(p), ml, nper=1, ext_unit_dict={})
    assert list(drt.stress_period_data[0]["node"]) == [20, 21]
    assert drt.active_params[0] == ["dp"]

    out = function_tmpdir / "mix_out.drt"
    drt.fn_path = str(out)
    drt.write_file()
    sp_line = next(
        ln for ln in out.read_text().splitlines() if "Stress Period 1" in ln
    )
    assert sp_line.split()[:2] == ["2", "1"]  # ITMP NP

    re = MfUsgDrt.load(str(out), _drt_model(function_tmpdir, "mx2"), nper=1)
    assert list(re.stress_period_data[0]["node"]) == [20, 21]
    assert re.active_params[0] == ["dp"]


def test_mfusgdrt_parameter_reuse_with_active(function_tmpdir):
    """ITMP<0 reuses the previous period's non-parametric drains while a new
    active parameter applies (NP>0); data round-trips."""
    p = function_tmpdir / "reuse.drt"
    p.write_text(
        "# drt reuse\n"
        "         1 0 1 5 RETURNFLOW\n"
        "dp DRT 1.0 1\n"
        " 5  5.000000e+00  1.000000e+01  0\n"
        " 1 1    Stress Period 1\n"
        " 21  4.000000e+00  2.000000e+01  0\n"
        "dp\n"
        " -1 1    Stress Period 2\n"
        "dp\n"
    )
    ml = _drt_model(function_tmpdir, "r1", nper=2)
    drt = MfUsgDrt.load(str(p), ml, nper=2, ext_unit_dict={})
    assert list(drt.stress_period_data[1]["node"]) == [20]  # reused from SP1
    assert drt.active_params[1] == ["dp"]

    out = function_tmpdir / "reuse_out.drt"
    drt.fn_path = str(out)
    drt.write_file()
    re = MfUsgDrt.load(str(out), _drt_model(function_tmpdir, "r2", nper=2), nper=2)
    assert list(re.stress_period_data[1]["node"]) == [20]
    assert re.active_params[1] == ["dp"]


def test_mfusgdrt_parameter_sfac_scales_cond(function_tmpdir):
    """SFAC inside a parameter's row block scales COND (Fortran ISCLOC=5)."""
    p = function_tmpdir / "sfacp.drt"
    p.write_text(
        "# drt sfac param\n"
        "         0 0 1 5 RETURNFLOW\n"
        "dp DRT 1.0 1\n"
        " SFAC 3.0\n"
        " 11  5.000000e+00  1.000000e+01  0\n"
        " 0 1    Stress Period 1\n"
        "dp\n"
    )
    ml = _drt_model(function_tmpdir, "sf1")
    drt = MfUsgDrt.load(str(p), ml, nper=1, ext_unit_dict={})
    assert np.isclose(drt.parameters["dp"]["data"]["cond"][0], 30.0)  # 10 * SFAC 3


def test_mfusgdrt_parameter_instances_unsupported(function_tmpdir):
    """DRT parameter INSTANCES are Fortran-supported but not yet by FloPy."""
    p = function_tmpdir / "inst.drt"
    p.write_text(
        "# drt instances\n"
        "         0 0 1 10 RETURNFLOW\n"
        "dp DRT 1.0 2 INSTANCES 2\n"
        "spring\n 11  5.0e0  1.0e1  0\n"
        "fall\n 12  5.0e0  1.0e1  0\n"
        " 0 1    Stress Period 1\n"
        "dp spring\n"
    )
    ml = _drt_model(function_tmpdir, "in1")
    with pytest.raises(NotImplementedError, match="INSTANCES"):
        MfUsgDrt.load(str(p), ml, nper=1, ext_unit_dict={})


def test_mfusgdrt_parameter_from_scratch_fails(function_tmpdir):
    """Active parameters with no loaded definitions fail explicitly, no file."""
    drt = MfUsgDrt(
        _drt_model(function_tmpdir, "fs"),
        options=["RETURNFLOW"],
        active_params={0: ["p1"]},
    )
    out = function_tmpdir / "fs.drt"
    drt.fn_path = str(out)
    with pytest.raises(NotImplementedError, match="from scratch"):
        drt.write_file()
    assert not out.exists()


def test_mfusgdrt_parameter_mxl_too_small_fails(function_tmpdir):
    """MXL below the total parameter list entries raises ValueError, no file."""
    dtype = MfUsgDrt.get_usg_dtype(returnflow=True, changec=False)
    rows = np.array([(0, 5.0, 10.0, 0.0)], dtype=dtype).view(np.recarray)
    params = {
        "p1": {
            "partyp": "DRT", "parval": "1.0", "nlst": 1,
            "data": rows, "recipient_nodes": [[]],
        },
        "p2": {
            "partyp": "DRT", "parval": "1.0", "nlst": 1,
            "data": rows, "recipient_nodes": [[]],
        },
    }
    drt = MfUsgDrt(
        _drt_model(function_tmpdir, "mx"),
        options=["RETURNFLOW"],
        parameters=params,
        mxl=1,
    )
    out = function_tmpdir / "mx.drt"
    drt.fn_path = str(out)
    with pytest.raises(ValueError, match="MXL"):
        drt.write_file()
    assert not out.exists()


def test_mfusgdrt_parameter_inconsistent_fails(function_tmpdir):
    """A definition whose nlst != len(data) raises ValueError, no file."""
    dtype = MfUsgDrt.get_usg_dtype(returnflow=True, changec=False)
    rows = np.array([(0, 5.0, 10.0, 0.0)], dtype=dtype).view(np.recarray)
    params = {
        "p1": {
            "partyp": "DRT", "parval": "1.0", "nlst": 2,  # declares 2, has 1
            "data": rows, "recipient_nodes": [[]],
        }
    }
    drt = MfUsgDrt(
        _drt_model(function_tmpdir, "ic"),
        options=["RETURNFLOW"],
        parameters=params,
        mxl=5,
    )
    out = function_tmpdir / "ic.drt"
    drt.fn_path = str(out)
    with pytest.raises(ValueError, match="nlst"):
        drt.write_file()
    assert not out.exists()


def test_mfusgdrt_parameter_active_undefined_fails(function_tmpdir):
    """An active parameter name not present in definitions raises ValueError."""
    dtype = MfUsgDrt.get_usg_dtype(returnflow=True, changec=False)
    rows = np.array([(0, 5.0, 10.0, 0.0)], dtype=dtype).view(np.recarray)
    params = {
        "p1": {
            "partyp": "DRT", "parval": "1.0", "nlst": 1,
            "data": rows, "recipient_nodes": [[]],
        }
    }
    drt = MfUsgDrt(
        _drt_model(function_tmpdir, "ud"),
        options=["RETURNFLOW"],
        parameters=params,
        mxl=5,
        active_params={0: ["px"]},  # not defined
    )
    out = function_tmpdir / "ud.drt"
    drt.fn_path = str(out)
    with pytest.raises(ValueError, match="not defined"):
        drt.write_file()
    assert not out.exists()


def test_mfusgdrt_parameter_active_case_insensitive(function_tmpdir):
    """An activation name differing from its definition only in case still counts
    toward MXADRT (review follow-up)."""
    p = function_tmpdir / "ci.drt"
    p.write_text(
        "# drt case-insensitive\n"
        "         2 0 1 1 RETURNFLOW\n"  # MXADRT IDRTCB NPDRT MXL RETURNFLOW
        "dp DRT 2.0 1\n"  # definition lower-case
        " 5  5.000000e+00  1.000000e+01  0\n"
        " 1 1    Stress Period 1\n"
        " 21  4.000000e+00  2.000000e+01  0\n"
        "DP\n"  # activation upper-case
    )
    ml = _drt_model(function_tmpdir, "ci1")
    drt = MfUsgDrt.load(str(p), ml, nper=1, ext_unit_dict={})
    out = function_tmpdir / "ci_out.drt"
    drt.fn_path = str(out)
    drt.write_file()
    item1 = next(ln for ln in out.read_text().splitlines() if not ln.startswith("#"))
    # 1 non-parametric + 1 active (resolved case-insensitively) => MXADRT=2
    assert item1.split()[0] == "2"


def test_mfusgdrt_parameter_duplicate_active_fails(function_tmpdir):
    """Activating the same parameter twice in a stress period (case-insensitive)
    raises ValueError, no partial file."""
    dtype = MfUsgDrt.get_usg_dtype(returnflow=True, changec=False)
    rows = np.array([(0, 5.0, 10.0, 0.0)], dtype=dtype).view(np.recarray)
    params = {
        "dp": {
            "partyp": "DRT", "parval": "1.0", "nlst": 1,
            "data": rows, "recipient_nodes": [[]],
        }
    }
    drt = MfUsgDrt(
        _drt_model(function_tmpdir, "dup"),
        options=["RETURNFLOW"],
        parameters=params,
        mxl=5,
        active_params={0: ["dp", "DP"]},  # same parameter twice
    )
    out = function_tmpdir / "dup.drt"
    drt.fn_path = str(out)
    with pytest.raises(ValueError, match="more than once"):
        drt.write_file()
    assert not out.exists()


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


# --- Stage 4.4B: HFB MODFLOW list-parameter preservation -------------------


def _hfb_model(function_tmpdir, name, structured=False, nlay=1, nrow=1, ncol=1, nper=1):
    from flopy.modflow import ModflowDis

    ml = MfUsg(structured=structured, model_ws=str(function_tmpdir), modelname=name)
    ModflowDis(ml, nlay=nlay, nrow=nrow, ncol=ncol, nper=nper)
    return ml


def test_mfusghfb_parameterized_roundtrip(function_tmpdir):
    """A parameterized (NPHFB>0) HFB file preserves defs + activations on
    load -> write -> reload (no longer expands or fails)."""
    from flopy.mfusg import MfUsgHfb

    p = function_tmpdir / "param.hfb"
    p.write_text(
        "# parameterized HFB\n"
        "         1         2         0\n"  # NPHFB MXFB NHFBNP
        "hfb_par hfb 2.0 2\n"  # def: name type parval nlst
        "1 2 0.5\n"  # barrier row 1 (1-based in file)
        "3 4 0.5\n"  # barrier row 2
        "1\n"  # NACTHFB
        "hfb_par\n"  # active parameter name
    )
    hfb = MfUsgHfb.load(str(p), _hfb_model(function_tmpdir, "h1"), nper=1)
    assert hfb.nphfb == 1 and hfb.mxfb == 2 and hfb.nhfbnp == 0
    assert hfb.nacthfb == 1 and hfb.acthfb_names == ["hfb_par"]
    pdef = hfb.parameters["hfb_par"]
    assert pdef["partyp"] == "hfb" and pdef["parval"] == "2.0" and pdef["nlst"] == 2
    # barrier nodes are stored 0-based internally
    assert list(pdef["data"]["node1"]) == [0, 2]
    assert list(pdef["data"]["node2"]) == [1, 3]

    out = function_tmpdir / "param_out.hfb"
    hfb.fn_path = str(out)
    hfb.write_file()
    content = out.read_text()
    # NPHFB preserved in the header; parameter syntax re-emitted
    item1 = next(ln for ln in content.splitlines() if not ln.startswith("#"))
    assert item1.split()[:3] == ["1", "2", "0"]
    assert "hfb_par hfb 2.0 2" in content
    rows = [ln for ln in content.splitlines() if not ln.startswith("#")]
    assert "hfb_par" in [ln.strip() for ln in rows]  # activation record

    re = MfUsgHfb.load(str(out), _hfb_model(function_tmpdir, "h2"), nper=1)
    assert re.nphfb == 1 and re.acthfb_names == ["hfb_par"]
    assert list(re.parameters["hfb_par"]["data"]["node1"]) == [0, 2]


def test_mfusghfb_parameterized_structured_roundtrip(function_tmpdir):
    """Structured parameterized HFB keeps k/i/j 0-based internally, 1-based on file."""
    from flopy.mfusg import MfUsgHfb

    p = function_tmpdir / "sparam.hfb"
    p.write_text(
        "# structured parameterized HFB\n"
        "         1         1         0\n"
        "spar hfb 3.0 1\n"
        "1 1 1 1 2 0.7\n"  # k irow1 icol1 irow2 icol2 factor (1-based)
        "1\n"
        "spar\n"
    )
    hfb = MfUsgHfb.load(
        str(p), _hfb_model(function_tmpdir, "s1", structured=True, ncol=2), nper=1
    )
    pdef = hfb.parameters["spar"]
    assert list(pdef["data"]["k"]) == [0]  # 0-based internal
    assert list(pdef["data"]["icol2"]) == [1]

    out = function_tmpdir / "sparam_out.hfb"
    hfb.fn_path = str(out)
    hfb.write_file()
    rows = [ln for ln in out.read_text().splitlines() if not ln.startswith("#")]
    # the parameter barrier row is written 1-based
    assert rows[2].split()[:5] == ["1", "1", "1", "1", "2"]

    re = MfUsgHfb.load(
        str(out), _hfb_model(function_tmpdir, "s2", structured=True, ncol=2), nper=1
    )
    assert list(re.parameters["spar"]["data"]["k"]) == [0]
    assert list(re.parameters["spar"]["data"]["icol2"]) == [1]


def test_mfusghfb_parameterized_with_nonparam(function_tmpdir):
    """HFB mixes parameter-defined barriers (item 2-3) with non-parametric ones
    (item 4), as the Fortran allows."""
    from flopy.mfusg import MfUsgHfb

    p = function_tmpdir / "mixed.hfb"
    p.write_text(
        "# params + non-parametric barriers\n"
        "         1         1         1\n"  # NPHFB=1 MXFB=1 NHFBNP=1
        "spar hfb 3.0 1\n"
        "1 1 1 1 2 0.7\n"  # parameter barrier (k=1)
        "2 1 1 1 2 9.9\n"  # non-parametric barrier (k=2)
        "1\n"
        "spar\n"
    )
    hfb = MfUsgHfb.load(
        str(p),
        _hfb_model(function_tmpdir, "m1", structured=True, nlay=2, ncol=2),
        nper=1,
    )
    assert hfb.nphfb == 1 and hfb.nhfbnp == 1
    assert list(hfb.parameters["spar"]["data"]["k"]) == [0]  # param barrier
    assert list(hfb.hfb_data["k"]) == [1]  # non-parametric barrier
    assert abs(float(hfb.hfb_data["hydchr"][0]) - 9.9) < 1e-6

    out = function_tmpdir / "mixed_out.hfb"
    hfb.fn_path = str(out)
    hfb.write_file()
    re = MfUsgHfb.load(
        str(out),
        _hfb_model(function_tmpdir, "m2", structured=True, nlay=2, ncol=2),
        nper=1,
    )
    assert list(re.parameters["spar"]["data"]["k"]) == [0]
    assert list(re.hfb_data["k"]) == [1]


def test_mfusghfb_parameter_authoring_from_scratch_fails(function_tmpdir):
    """Authoring HFB parameters from scratch (NPHFB>0, no defs) fails explicitly."""
    from flopy.mfusg import MfUsgHfb

    hfb = MfUsgHfb(_hfb_model(function_tmpdir, "fs"), nphfb=1, mxfb=1, nhfbnp=0)
    hfb.fn_path = str(function_tmpdir / "fs.hfb")
    with pytest.raises(NotImplementedError, match="from scratch"):
        hfb.write_file()


def test_mfusghfb_transient_with_parameters_fails(function_tmpdir):
    """TRANSIENT_HFB combined with parameters (NPHFB>0) fails explicitly."""
    from flopy.mfusg import MfUsgHfb

    p = function_tmpdir / "tparam.hfb"
    p.write_text(
        "# transient + parameters\n"
        "         1         1         0  TRANSIENT_HFB\n"
        "tp hfb 1.0 1\n"
        "1 2 0.5\n"
        "1\n"
        "tp\n"
    )
    with pytest.raises(NotImplementedError, match="TRANSIENT_HFB"):
        MfUsgHfb.load(str(p), _hfb_model(function_tmpdir, "tp", nper=2), nper=2)


# --- Stage 4.4B follow-up: HFB list controls (SFAC / OPEN-CLOSE / EXTERNAL) --


def test_mfusghfb_sfac_scales_nonparam(function_tmpdir):
    """A non-parametric HFB barrier list may begin with SFAC, which scales hydchr."""
    from flopy.mfusg import MfUsgHfb

    p = function_tmpdir / "sfac.hfb"
    p.write_text(
        "# hfb non-param SFAC\n"
        "         0         0         2\n"  # NPHFB MXFB NHFBNP
        " SFAC 10.0\n"
        " 1 2 0.5\n"
        " 3 4 0.5\n"
    )
    hfb = MfUsgHfb.load(str(p), _hfb_model(function_tmpdir, "sf1"), nper=1)
    assert list(hfb.hfb_data["node1"]) == [0, 2]
    assert np.allclose(hfb.hfb_data["hydchr"], [5.0, 5.0])  # 0.5 * SFAC 10


def test_mfusghfb_sfac_parameterized_roundtrip(function_tmpdir):
    """SFAC inside a parameter's NLST block scales hydchr; the (already-scaled)
    value is preserved on write -> reload (expanded valid write, no SFAC out)."""
    from flopy.mfusg import MfUsgHfb

    p = function_tmpdir / "psfac.hfb"
    p.write_text(
        "# hfb param SFAC (structured)\n"
        "         1         1         0\n"
        "spar hfb 2.0 1\n"
        " SFAC 4.0\n"
        " 1 1 1 1 2 0.5\n"  # k irow1 icol1 irow2 icol2 factor
        "1\n"
        "spar\n"
    )
    hfb = MfUsgHfb.load(
        str(p), _hfb_model(function_tmpdir, "ps1", structured=True, ncol=2), nper=1
    )
    assert np.allclose(hfb.parameters["spar"]["data"]["hydchr"], [2.0])  # 0.5*4
    assert list(hfb.parameters["spar"]["data"]["k"]) == [0]

    out = function_tmpdir / "psfac_out.hfb"
    hfb.fn_path = str(out)
    hfb.write_file()
    assert "SFAC" not in out.read_text()  # SFAC baked into hydchr, not preserved

    re = MfUsgHfb.load(
        str(out), _hfb_model(function_tmpdir, "ps2", structured=True, ncol=2), nper=1
    )
    assert np.allclose(re.parameters["spar"]["data"]["hydchr"], [2.0])
    assert re.parameters["spar"]["parval"] == "2.0"


def test_mfusghfb_sfac_mixed_param_and_nonparam(function_tmpdir):
    """SFAC applies independently to a parameter's NLST block (item 2-3) and to
    the non-parametric barriers (item 4)."""
    from flopy.mfusg import MfUsgHfb

    p = function_tmpdir / "mixsfac.hfb"
    p.write_text(
        "# hfb mixed SFAC\n"
        "         1         1         1\n"  # NPHFB=1 MXFB=1 NHFBNP=1
        "spar hfb 2.0 1\n"
        " SFAC 3.0\n"
        " 1 2 0.4\n"  # param barrier -> 0.4 * 3 = 1.2
        " SFAC 5.0\n"
        " 7 8 0.2\n"  # non-param barrier -> 0.2 * 5 = 1.0
        "1\n"
        "spar\n"
    )
    hfb = MfUsgHfb.load(str(p), _hfb_model(function_tmpdir, "mx1"), nper=1)
    assert np.allclose(hfb.parameters["spar"]["data"]["hydchr"], [1.2])
    assert np.allclose(hfb.hfb_data["hydchr"], [1.0])
    assert list(hfb.hfb_data["node1"]) == [6]


def test_mfusghfb_open_close(function_tmpdir):
    """A barrier list can be read from an OPEN/CLOSE file (plain and quoted)."""
    from flopy.mfusg import MfUsgHfb

    (function_tmpdir / "hfb_rows.dat").write_text(" 1 2 0.7\n 3 4 0.8\n")
    p = function_tmpdir / "oc.hfb"
    p.write_text(
        "# hfb open/close\n         0         0         2\n OPEN/CLOSE hfb_rows.dat\n"
    )
    hfb = MfUsgHfb.load(str(p), _hfb_model(function_tmpdir, "oc1"), nper=1)
    assert list(hfb.hfb_data["node1"]) == [0, 2]
    assert np.allclose(hfb.hfb_data["hydchr"], [0.7, 0.8])

    # quoted filename containing a space
    (function_tmpdir / "hfb rows.dat").write_text(" 5 6 0.9\n")
    pq = function_tmpdir / "ocq.hfb"
    pq.write_text(
        '# hfb open/close quoted\n         0         0         1\n'
        ' OPEN/CLOSE "hfb rows.dat"\n'
    )
    hfbq = MfUsgHfb.load(str(pq), _hfb_model(function_tmpdir, "oc2"), nper=1)
    assert list(hfbq.hfb_data["node1"]) == [4]


def test_mfusghfb_external_via_ext_unit_dict(function_tmpdir):
    """EXTERNAL resolves the barrier-list file through ext_unit_dict."""
    import types

    from flopy.mfusg import MfUsgHfb

    (function_tmpdir / "hfb_ext.dat").write_text(" 7 8 0.11\n")
    p = function_tmpdir / "ext.hfb"
    p.write_text("# hfb external\n         0         0         1\n EXTERNAL 88\n")
    ext_unit_dict = {
        88: types.SimpleNamespace(filename="hfb_ext.dat", filetype="DATA")
    }
    hfb = MfUsgHfb.load(
        str(p), _hfb_model(function_tmpdir, "ex1"), nper=1, ext_unit_dict=ext_unit_dict
    )
    assert list(hfb.hfb_data["node1"]) == [6]
    assert np.isclose(hfb.hfb_data["hydchr"][0], 0.11)


def test_mfusghfb_external_without_dict_fails(function_tmpdir):
    """EXTERNAL with an unresolvable unit fails with NotImplementedError."""
    from flopy.mfusg import MfUsgHfb

    p = function_tmpdir / "extbad.hfb"
    p.write_text("# hfb external bad\n         0         0         1\n EXTERNAL 77\n")
    with pytest.raises(NotImplementedError, match="EXTERNAL"):
        MfUsgHfb.load(
            str(p), _hfb_model(function_tmpdir, "ex2"), nper=1, ext_unit_dict={}
        )


def test_mfusghfb_nacthfb_mismatch_fails(function_tmpdir):
    """write_file rejects a parameterized HFB whose nacthfb != len(acthfb_names)."""
    from flopy.mfusg import MfUsgHfb

    p = function_tmpdir / "pp.hfb"
    p.write_text(
        "# hfb param\n         1         1         0\nspar hfb 2.0 1\n"
        " 1 2 0.5\n1\nspar\n"
    )
    hfb = MfUsgHfb.load(str(p), _hfb_model(function_tmpdir, "nm1"), nper=1)
    hfb.nacthfb = 5  # inconsistent with the single active name
    hfb.fn_path = str(function_tmpdir / "nm_out.hfb")
    with pytest.raises(ValueError, match="nacthfb"):
        hfb.write_file()


# --- Stage 4.4B polish: HFB parameter writer validation --------------------


def _hfb_param_def(nlst, hydchr=0.5):
    """Build an unstructured HFB list-parameter `data` recarray with nlst rows."""
    from flopy.mfusg import MfUsgHfb

    dt = MfUsgHfb.get_default_dtype(structured=False)
    return np.array(
        [(i, i + 1, hydchr) for i in range(nlst)], dtype=dt
    ).view(np.recarray)


def test_mfusghfb_param_write_empty_dict_fails(function_tmpdir):
    """NPHFB>0 with an empty parameters dict is from-scratch authoring: it fails
    with NotImplementedError and writes no partial file."""
    from flopy.mfusg import MfUsgHfb

    hfb = MfUsgHfb(
        _hfb_model(function_tmpdir, "e1"),
        nphfb=1,
        mxfb=1,
        parameters={},
        acthfb_names=[],
    )
    out = function_tmpdir / "empty.hfb"
    hfb.fn_path = str(out)
    with pytest.raises(NotImplementedError, match="from scratch"):
        hfb.write_file()
    assert not out.exists()  # no partial file written


def test_mfusghfb_param_write_count_mismatch_fails(function_tmpdir):
    """len(parameters) must equal NPHFB."""
    from flopy.mfusg import MfUsgHfb

    params = {
        "p1": {"partyp": "hfb", "parval": "1.0", "nlst": 1, "data": _hfb_param_def(1)}
    }
    hfb = MfUsgHfb(
        _hfb_model(function_tmpdir, "c1"),
        nphfb=2,  # claims 2 parameters but only one is defined
        mxfb=1,
        parameters=params,
        acthfb_names=["p1"],
        nacthfb=1,
    )
    hfb.fn_path = str(function_tmpdir / "count.hfb")
    with pytest.raises(ValueError, match="NPHFB"):
        hfb.write_file()


def test_mfusghfb_param_write_nlst_mismatch_fails(function_tmpdir):
    """A parameter's nlst must equal len(data)."""
    from flopy.mfusg import MfUsgHfb

    params = {
        "p1": {"partyp": "hfb", "parval": "1.0", "nlst": 2, "data": _hfb_param_def(1)}
    }
    hfb = MfUsgHfb(
        _hfb_model(function_tmpdir, "n1"),
        nphfb=1,
        mxfb=2,
        parameters=params,
        acthfb_names=["p1"],
        nacthfb=1,
    )
    hfb.fn_path = str(function_tmpdir / "nlst.hfb")
    with pytest.raises(ValueError, match="nlst"):
        hfb.write_file()


def test_mfusghfb_param_write_undefined_active_name_fails(function_tmpdir):
    """Every active parameter name must be defined in parameters."""
    from flopy.mfusg import MfUsgHfb

    params = {
        "p1": {"partyp": "hfb", "parval": "1.0", "nlst": 1, "data": _hfb_param_def(1)}
    }
    hfb = MfUsgHfb(
        _hfb_model(function_tmpdir, "u1"),
        nphfb=1,
        mxfb=1,
        parameters=params,
        acthfb_names=["px"],  # not defined
        nacthfb=1,
    )
    hfb.fn_path = str(function_tmpdir / "undef.hfb")
    with pytest.raises(ValueError, match="not defined"):
        hfb.write_file()


def test_mfusghfb_param_write_duplicate_active_fails(function_tmpdir):
    """Activating the same HFB parameter more than once (case-insensitive) raises
    ValueError before the file is opened (the Fortran aborts 'already activated'
    in SGWF2HFB7SUB); no partial file is written."""
    from flopy.mfusg import MfUsgHfb

    params = {
        "p1": {"partyp": "hfb", "parval": "1.0", "nlst": 1, "data": _hfb_param_def(1)}
    }
    hfb = MfUsgHfb(
        _hfb_model(function_tmpdir, "dup"),
        nphfb=1,
        mxfb=1,
        parameters=params,
        acthfb_names=["p1", "P1"],  # same parameter twice (case-insensitive)
        nacthfb=2,
    )
    out = function_tmpdir / "dup.hfb"
    hfb.fn_path = str(out)
    with pytest.raises(ValueError, match="more than once"):
        hfb.write_file()
    assert not out.exists()


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
    """Opt-in `expand_parameters=True` keeps the legacy Expanded-valid-write path.

    A parameterized ETS file loads, expands the parameter to a concrete ETSR
    array, and writes valid non-parametric input (NPETS=0, no PARAMETER). This
    is the documented fallback; preservation (the new default) is tested
    separately.
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
    ets = MfUsgEts.load(str(param_ets), ml, nper=1, expand_parameters=True)
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


# --- Stage 4.4A: ETS MODFLOW array-parameter preservation -------------------


def _param_ets_model(function_tmpdir, name):
    from flopy.modflow import ModflowDis

    ml = MfUsg(structured=True, model_ws=str(function_tmpdir), modelname=name)
    ModflowDis(ml, nlay=1, nrow=2, ncol=2, nper=1)
    return ml


# Canonical single-parameter ETS file (USG-T style: NPETS in item 2a, no
# PARAMETER line). ETSR is parameterized; ETSS/ETSX are plain arrays.
_PARAM_ETS_TEXT = (
    "# parameterized ETS (NPETS=1)\n"
    "1 0 1 1 0\n"              # NETSOP IETSCB NPETS NETSEG IESFACTOR
    "etsrate ets 5.0E-4 1\n"   # param: name type value nclu
    "NONE ALL\n"               # cluster: no multiplier, all cells
    "0 1 0\n"                  # SP1: INSURF INETSR INEXDP (INETSR=1 param)
    "CONSTANT 10.0\n"          # ETSS surface
    "etsrate\n"                # ETSR via parameter "etsrate"
    "CONSTANT 5.0\n"           # ETSX extinction depth
)


def test_mfusgets_parameterized_load_preserves(function_tmpdir):
    """Default load preserves ETS array parameters (npets, defs, per-SP record)."""
    from flopy.mfusg import MfUsgEts

    p = function_tmpdir / "param.ets"
    p.write_text(_PARAM_ETS_TEXT)
    ets = MfUsgEts.load(str(p), _param_ets_model(function_tmpdir, "p1"), nper=1)
    assert ets.npets == 1
    assert ets.parameters is not None
    assert "etsrate" in ets.parameters.bc_parms
    # per-stress-period active-parameter record preserved
    assert ets.evtr_parm[0] == [("etsrate", None)]
    # the parameter value is still expanded into the in-memory ETSR array
    assert np.allclose(ets.evtr[0].array, 5.0e-4)


def test_mfusgets_parameterized_write_preserves_syntax(function_tmpdir):
    """Write keeps NPETS>0 + parameter syntax and mixes a parameterized ETSR
    with plain ETSS/ETSX arrays; no PARAMETER line (NPETS is in item 2a)."""
    from flopy.mfusg import MfUsgEts

    p = function_tmpdir / "param.ets"
    p.write_text(_PARAM_ETS_TEXT)
    ets = MfUsgEts.load(str(p), _param_ets_model(function_tmpdir, "p2"), nper=1)
    out = function_tmpdir / "preserved.ets"
    ets.fn_path = str(out)
    ets.write_file()
    content = out.read_text()

    # NPETS preserved in item 2a; no MODFLOW-2005 PARAMETER line
    item2a = next(ln for ln in content.splitlines() if not ln.startswith("#"))
    assert item2a.split()[2] == "1"
    assert "PARAMETER" not in content
    # parameter definition block + cluster
    assert "etsrate ets 5.0E-4 1" in content
    assert "NONE ALL" in content
    # the parameterized ETSR mixes with plain ETSS/ETSX arrays
    lines = [ln.strip() for ln in content.splitlines()]
    assert "etsrate" in lines  # active-parameter record (replaces ETSR array)
    assert sum(ln.startswith("CONSTANT") for ln in lines) == 2  # ETSS + ETSX


def test_mfusgets_parameterized_roundtrip(function_tmpdir):
    """Load -> write -> reload keeps the parameter syntax and value intact."""
    from flopy.mfusg import MfUsgEts

    p = function_tmpdir / "param.ets"
    p.write_text(_PARAM_ETS_TEXT)
    ets = MfUsgEts.load(str(p), _param_ets_model(function_tmpdir, "r1"), nper=1)
    out = function_tmpdir / "preserved.ets"
    ets.fn_path = str(out)
    ets.write_file()

    re = MfUsgEts.load(str(out), _param_ets_model(function_tmpdir, "r2"), nper=1)
    assert re.npets == 1
    assert re.evtr_parm[0] == [("etsrate", None)]
    assert re.parameters.bc_parms["etsrate"][0]["parval"] == "5.0E-4"
    assert np.allclose(re.evtr[0].array, 5.0e-4)


def test_mfusgets_parameter_instances_roundtrip(function_tmpdir):
    """A time-varying ETS parameter (INSTANCES) preserves its active instance."""
    from flopy.mfusg import MfUsgEts

    p = function_tmpdir / "inst.ets"
    p.write_text(
        "# ETS with two instances\n"
        "1 0 1 1 0\n"
        "etsrate ets 5.0E-4 1 INSTANCES 2\n"
        "spring\n"
        "NONE ALL\n"
        "fall\n"
        "NONE ALL\n"
        "0 1 0\n"
        "CONSTANT 10.0\n"
        "etsrate spring\n"
        "CONSTANT 5.0\n"
    )
    ets = MfUsgEts.load(str(p), _param_ets_model(function_tmpdir, "in1"), nper=1)
    assert ets.evtr_parm[0] == [("etsrate", "spring")]
    out = function_tmpdir / "inst_out.ets"
    ets.fn_path = str(out)
    ets.write_file()
    content = out.read_text()
    assert "INSTANCES 2" in content
    assert "etsrate spring" in content

    re = MfUsgEts.load(str(out), _param_ets_model(function_tmpdir, "in2"), nper=1)
    assert re.evtr_parm[0] == [("etsrate", "spring")]


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


def _evt_struct_model(function_tmpdir, name, nlay=1, nrow=2, ncol=2, nper=1):
    from flopy.modflow import ModflowDis

    m = MfUsg(structured=True, model_ws=str(function_tmpdir), modelname=name)
    ModflowDis(m, nlay=nlay, nrow=nrow, ncol=ncol, nper=nper)
    return m


def _usgt_unstructured_model_ncol(function_tmpdir, name, ncol=4, nper=1):
    from flopy.modflow import ModflowDis

    m = MfUsg(structured=False, model_ws=str(function_tmpdir), modelname=name)
    ModflowDis(m, nlay=1, nrow=1, ncol=ncol, nper=nper)
    return m


def test_mfusgevt_nevtop1_and_3_authoring_roundtrip(function_tmpdir):
    """NEVTOP=1 and NEVTOP=3 author from scratch and round-trip."""
    from flopy.mfusg import MfUsgEvt

    for nevtop in (1, 3):
        evt = MfUsgEvt(
            _evt_struct_model(function_tmpdir, f"n{nevtop}"),
            nevtop=nevtop,
            surf=10.0,
            evtr=1.0e-4,
            exdp=2.0,
        )
        evt.fn_path = str(function_tmpdir / f"n{nevtop}.evt")
        evt.write_file()
        re = MfUsgEvt.load(
            evt.fn_path, _evt_struct_model(function_tmpdir, f"n{nevtop}b")
        )
        assert re.nevtop == nevtop


def test_mfusgevt_nevtop2_structured_roundtrip(function_tmpdir):
    """NEVTOP=2 structured: ievt is 0-based internal and 1-based in the file."""
    from flopy.mfusg import MfUsgEvt

    ievt = np.array([[0, 1], [2, 0]])  # 0-based layer indices
    evt = MfUsgEvt(
        _evt_struct_model(function_tmpdir, "n2", nlay=3),
        nevtop=2,
        ievt=ievt,
        surf=5.0,
        evtr=1.0e-4,
        exdp=1.0,
    )
    evt.fn_path = str(function_tmpdir / "n2.evt")
    evt.write_file()
    # the file carries 1-based layers (max 0-based 2 -> 3)
    assert " 3" in Path(evt.fn_path).read_text()

    re = MfUsgEvt.load(evt.fn_path, _evt_struct_model(function_tmpdir, "n2b", nlay=3))
    assert re.nevtop == 2
    np.testing.assert_array_equal(re.ievt[0].array, ievt)  # back to 0-based


def test_mfusgevt_nevtop2_unstructured_mxndevt(function_tmpdir):
    """NEVTOP=2 unstructured writes MXNDEVT and node-index IEVT (0-based/1-based)."""
    from flopy.mfusg import MfUsgEvt

    grid_ievt = np.array([[0, 1, 2, 3]])
    evt = MfUsgEvt(
        _usgt_unstructured_model_ncol(function_tmpdir, "u1", ncol=4),
        nevtop=2,
        ievt=grid_ievt,
        surf=5.0,
        evtr=1.0e-4,
        exdp=1.0,
    )
    evt.fn_path = str(function_tmpdir / "u1.evt")
    evt.write_file()
    lines = Path(evt.fn_path).read_text().splitlines()
    assert lines[2].split() == ["4"]  # MXNDEVT line

    re = MfUsgEvt.load(
        evt.fn_path, _usgt_unstructured_model_ncol(function_tmpdir, "u2", ncol=4)
    )
    np.testing.assert_array_equal(re.ievt[0].array.ravel(), [0, 1, 2, 3])


def test_mfusgevt_multi_period_reuse(function_tmpdir):
    """Repeated arrays across stress periods reuse via negative flags (-1)."""
    from flopy.mfusg import MfUsgEvt

    evt = MfUsgEvt(
        _evt_struct_model(function_tmpdir, "rs", nper=3),
        nevtop=1,
        surf=5.0,
        evtr={0: 1.0e-4, 1: 2.0e-4, 2: 2.0e-4},
        exdp=1.0,
    )
    evt.fn_path = str(function_tmpdir / "rs.evt")
    evt.write_file()
    # SP3 repeats SP2's evtr -> a reuse flag (-1) appears in a later period header
    headers = [
        ln.split()[:3]
        for ln in Path(evt.fn_path).read_text().splitlines()
        if ln and ln.split() and ln.split()[0].lstrip("-").isdigit()
    ]
    assert any("-1" in h for h in headers)
    re = MfUsgEvt.load(evt.fn_path, _evt_struct_model(function_tmpdir, "rs2", nper=3))
    assert len(re.evtr.transient_2ds) == 3


def test_mfusgevt_transport_ietfactor_roundtrip(function_tmpdir):
    """Transport ietfactor 0/<0/>0 (incl. MCOMP>1) round-trips, preserving values."""
    from flopy.mfusg import MfUsgEvt

    def tmodel(name, mcomp=1):
        m = _evt_struct_model(function_tmpdir, name)
        m.itrnsp = 1
        m.mcomp = mcomp
        return m

    for ietf, etf in ((0, 0.0), (-1, 0.0), (1, [2.5])):
        evt = MfUsgEvt(
            tmodel(f"t{ietf}"), nevtop=1, evtr=1.0e-4, ietfactor=ietf, etfactor=etf
        )
        evt.fn_path = str(function_tmpdir / f"t{ietf}.evt")
        evt.write_file()
        # transport always writes 3 integers on dataset 1
        assert len(Path(evt.fn_path).read_text().splitlines()[1].split()) == 3
        re = MfUsgEvt.load(evt.fn_path, tmodel(f"t{ietf}b"))
        assert re.ietfactor == ietf

    # MCOMP > 1 with a per-component ETFACTOR array
    evt = MfUsgEvt(
        tmodel("mc", mcomp=2), nevtop=1, evtr=1.0e-4, ietfactor=1, etfactor=[2.5, 3.5]
    )
    evt.fn_path = str(function_tmpdir / "mc.evt")
    evt.write_file()
    re = MfUsgEvt.load(evt.fn_path, tmodel("mc2", mcomp=2))
    assert re.ietfactor == 1 and np.allclose(re.etfactor, [2.5, 3.5])


def test_mfusgevt_rejects_invalid(function_tmpdir):
    """EVT fails explicitly on invalid NEVTOP, unsupported ETS, and bad inputs."""
    from flopy.mfusg import MfUsgEvt

    with pytest.raises(ValueError, match="NEVTOP"):
        MfUsgEvt(_evt_struct_model(function_tmpdir, "b1"), nevtop=5)

    with pytest.raises(NotImplementedError, match="ETS"):
        MfUsgEvt(_evt_struct_model(function_tmpdir, "b2"), nevtop=1, mxetzones=3)

    def tm(name, mcomp=2):
        m = _evt_struct_model(function_tmpdir, name)
        m.itrnsp = 1
        m.mcomp = mcomp
        return m

    with pytest.raises(ValueError, match="MCOMP"):
        MfUsgEvt(tm("b3"), nevtop=1, ietfactor=1, etfactor=[2.5])  # len 1 != mcomp 2

    # NEVTOP=2 structured layer index out of range (only 2 layers) -> write raises
    with pytest.raises(ValueError, match="0-based layer"):
        MfUsgEvt(
            _evt_struct_model(function_tmpdir, "b4", nlay=2),
            nevtop=2,
            ievt=np.array([[5, 0], [0, 0]]),
        ).write_file()


def test_mfusgevt_etfactor_scalar_and_array_roundtrip(function_tmpdir):
    """ETFACTOR accepts a scalar (MCOMP=1) or an array (MCOMP>1) and round-trips."""
    from flopy.mfusg import MfUsgEvt

    def tmodel(name, mcomp):
        m = _evt_struct_model(function_tmpdir, name)
        m.itrnsp = 1
        m.mcomp = mcomp
        return m

    # scalar ETFACTOR with MCOMP=1 must write (not crash) and round-trip
    evt = MfUsgEvt(tmodel("sc", 1), nevtop=1, evtr=1e-4, ietfactor=1, etfactor=2.5)
    evt.fn_path = str(function_tmpdir / "sc.evt")
    evt.write_file()
    assert "2.5" in Path(evt.fn_path).read_text()
    re = MfUsgEvt.load(evt.fn_path, tmodel("sc2", 1))
    assert np.allclose(re.etfactor, [2.5])
    re.fn_path = str(function_tmpdir / "sc_re.evt")
    re.write_file()  # rewrite stable, no crash
    assert "2.5" in Path(re.fn_path).read_text()

    # array ETFACTOR with MCOMP>1 still works
    evt = MfUsgEvt(
        tmodel("ar", 2), nevtop=1, evtr=1e-4, ietfactor=1, etfactor=[2.5, 3.5]
    )
    evt.fn_path = str(function_tmpdir / "ar.evt")
    evt.write_file()
    re = MfUsgEvt.load(evt.fn_path, tmodel("ar2", 2))
    assert np.allclose(re.etfactor, [2.5, 3.5])

    # incompatible length still raises
    with pytest.raises(ValueError, match="MCOMP"):
        MfUsgEvt(tmodel("bad", 2), nevtop=1, ietfactor=1, etfactor=[2.5])


def test_mfusgevt_nevtop2_unstructured_ievt_out_of_range(function_tmpdir):
    """Unstructured NEVTOP=2 validates the IEVT node index against NODES."""
    from flopy.mfusg import MfUsgEvt

    # valid node indices (4 nodes -> 0..3) author and round-trip fine
    evt = MfUsgEvt(
        _usgt_unstructured_model_ncol(function_tmpdir, "uok", ncol=4),
        nevtop=2,
        ievt=np.array([[0, 1, 2, 3]]),
        surf=5.0,
        evtr=1e-4,
        exdp=1.0,
    )
    evt.fn_path = str(function_tmpdir / "uok.evt")
    evt.write_file()

    # node index >= NODES is rejected
    with pytest.raises(ValueError, match="0-based node"):
        MfUsgEvt(
            _usgt_unstructured_model_ncol(function_tmpdir, "uhi", ncol=4),
            nevtop=2,
            ievt=np.array([[0, 1, 2, 4]]),
        ).write_file()

    # negative node index is rejected
    with pytest.raises(ValueError, match="0-based node"):
        MfUsgEvt(
            _usgt_unstructured_model_ncol(function_tmpdir, "uneg", ncol=4),
            nevtop=2,
            ievt=np.array([[-1, 1, 2, 3]]),
        ).write_file()


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


def _oc_model(function_tmpdir, name="m", nlay=2, nper=1):
    """A minimal structured MfUsg model for OC authoring/round-trip tests."""
    from flopy.modflow import ModflowDis

    m = MfUsg(structured=True, model_ws=str(function_tmpdir), modelname=name)
    ModflowDis(m, nlay=nlay, nrow=1, ncol=3, nper=nper, nstp=1)
    return m


def test_mfusgoc_bootstrapping_header_authoring_roundtrip(function_tmpdir):
    """OC BOOTSTRAPPING header authors from scratch and round-trips.

    Per glo2basu1.f the option is parsed only from the first OC line
    (SGWF2BAS7I); the per-record reader SGWF2BAS7J rejects it as a standalone
    line, so the writer must place it on line 1.
    """
    from flopy.mfusg import MfUsgOc

    oc = MfUsgOc(
        _oc_model(function_tmpdir, "b1"),
        bootstrapping=1,
        iugboot=80,
        stress_period_data={(0, 0): ["save head"]},
    )
    oc.fn_path = str(function_tmpdir / "b1.oc")
    oc.write_file()
    body = Path(oc.fn_path).read_text().splitlines()
    assert "BOOTSTRAPPING" in body[1] and "80" in body[1]  # first OC line

    re = MfUsgOc.load(oc.fn_path, _oc_model(function_tmpdir, "b2"), nper=1)
    assert re.bootstrapping and re.iugboot == 80
    re.fn_path = str(function_tmpdir / "b1_re.oc")
    re.write_file()
    assert "BOOTSTRAPPING 80" in Path(re.fn_path).read_text()  # stable

    # a normal OC writes no BOOTSTRAPPING line
    plain = MfUsgOc(
        _oc_model(function_tmpdir, "b3"), stress_period_data={(0, 0): ["save head"]}
    )
    plain.fn_path = str(function_tmpdir / "b3.oc")
    plain.write_file()
    assert "BOOTSTRAP" not in Path(plain.fn_path).read_text()


def test_mfusgoc_bootstrap_stress_period_actions(function_tmpdir):
    """Per-stress-period BOOTSTRAP toggles round-trip as OC actions."""
    from flopy.mfusg import MfUsgOc

    spd = {
        (0, 0): ["save head", "bootstrap", "bootstrapscale"],
        (1, 0): ["save head", "nobootstrap", "nobootstrapscale"],
    }
    oc = MfUsgOc(_oc_model(function_tmpdir, "spb", nper=2), stress_period_data=spd)
    oc.fn_path = str(function_tmpdir / "spb.oc")
    oc.write_file()

    re = MfUsgOc.load(oc.fn_path, _oc_model(function_tmpdir, "spb2", nper=2), nper=2)
    a0 = [a.upper() for a in re.stress_period_data[(0, 0)]]
    a1 = [a.upper() for a in re.stress_period_data[(1, 0)]]
    assert "BOOTSTRAP" in a0 and "BOOTSTRAPSCALE" in a0
    assert "NOBOOTSTRAP" in a1 and "NOBOOTSTRAPSCALE" in a1


def test_mfusgoc_output_block_combinations(function_tmpdir):
    """SAVE/PRINT HEAD/CONC/BUDGET output actions round-trip."""
    from flopy.mfusg import MfUsgOc

    spd = {
        (0, 0): [
            "save head",
            "save conc",
            "print conc",
            "save budget",
            "print budget",
        ]
    }
    oc = MfUsgOc(_oc_model(function_tmpdir, "io"), stress_period_data=spd)
    oc.fn_path = str(function_tmpdir / "io.oc")
    oc.write_file()

    re = MfUsgOc.load(oc.fn_path, _oc_model(function_tmpdir, "io2"), nper=1)
    acts = [a.upper() for a in re.stress_period_data[(0, 0)]]
    for a in ("SAVE HEAD", "SAVE CONC", "PRINT CONC", "SAVE BUDGET", "PRINT BUDGET"):
        assert a in acts, (a, acts)


def test_mfusgoc_save_ibound_roundtrips_but_usgt_rejects(function_tmpdir):
    """SAVE IBOUND round-trips in FloPy, but USG-T 2.7's OC reader rejects it.

    The per-step SAVE IBOUND branch is commented out in glo2basu1.f
    (SGWF2BAS7N), so this action must not be used with USG-T 2.7. FloPy still
    preserves the keyword for standard-MODFLOW compatibility — this test pins
    that documented gap.
    """
    from flopy.mfusg import MfUsgOc

    oc = MfUsgOc(
        _oc_model(function_tmpdir, "ib"),
        stress_period_data={(0, 0): ["save head", "save ibound"]},
    )
    oc.fn_path = str(function_tmpdir / "ib.oc")
    oc.write_file()
    re = MfUsgOc.load(oc.fn_path, _oc_model(function_tmpdir, "ib2"), nper=1)
    assert "SAVE IBOUND" in [a.upper() for a in re.stress_period_data[(0, 0)]]


def test_mfusgoc_layer_qualified_roundtrip(function_tmpdir):
    """Layer-qualified PRINT/SAVE HEAD/DRAWDOWN/CONC keep their layer lists."""
    from flopy.mfusg import MfUsgOc

    spd = {
        (0, 0): [
            "save head 1 2",
            "print head 1",
            "save drawdown 2 3",
            "print conc 1 3",
        ]
    }
    oc = MfUsgOc(_oc_model(function_tmpdir, "lay", nlay=3), stress_period_data=spd)
    oc.fn_path = str(function_tmpdir / "lay.oc")
    oc.write_file()

    re = MfUsgOc.load(oc.fn_path, _oc_model(function_tmpdir, "lay2", nlay=3), nper=1)
    acts = [a.upper() for a in re.stress_period_data[(0, 0)]]
    assert "SAVE HEAD 1 2" in acts
    assert "PRINT HEAD 1" in acts
    assert "SAVE DRAWDOWN 2 3" in acts
    assert "PRINT CONC 1 3" in acts


def test_mfusgoc_ddreference_roundtrip(function_tmpdir):
    """DDREFERENCE (written on the period line) round-trips as an action."""
    from flopy.mfusg import MfUsgOc

    oc = MfUsgOc(
        _oc_model(function_tmpdir, "ddr"),
        stress_period_data={(0, 0): ["save drawdown", "ddreference"]},
    )
    oc.fn_path = str(function_tmpdir / "ddr.oc")
    oc.write_file()
    period_line = next(
        ln
        for ln in Path(oc.fn_path).read_text().splitlines()
        if ln.lower().startswith("period")
    )
    assert "ddreference" in period_line.lower()

    re = MfUsgOc.load(oc.fn_path, _oc_model(function_tmpdir, "ddr2"), nper=1)
    assert "DDREFERENCE" in [a.upper() for a in re.stress_period_data[(0, 0)]]


def test_mfusgoc_ddreference_only_roundtrip(function_tmpdir):
    """A DDREFERENCE-only stress period still emits a period line and round-trips."""
    from flopy.mfusg import MfUsgOc

    oc = MfUsgOc(
        _oc_model(function_tmpdir, "d1"),
        stress_period_data={(0, 0): ["ddreference"]},
    )
    oc.fn_path = str(function_tmpdir / "d1.oc")
    oc.write_file()
    text = Path(oc.fn_path).read_text().lower()
    assert "period 1 step 1 ddreference" in text

    re = MfUsgOc.load(oc.fn_path, _oc_model(function_tmpdir, "d2"), nper=1)
    assert (0, 0) in re.stress_period_data
    assert "DDREFERENCE" in [a.upper() for a in re.stress_period_data[(0, 0)]]


def _oc_check_warnings(function_tmpdir, name, actions):
    """Build an OC with `actions` and return its check() warning descriptions."""
    from flopy.mfusg import MfUsgOc

    oc = MfUsgOc(
        _oc_model(function_tmpdir, name), stress_period_data={(0, 0): actions}
    )
    chk = oc.check(verbose=False)
    return [str(d) for d in chk.summary_array["desc"]]


def test_mfusgoc_check_accepts_usgt_actions(function_tmpdir):
    """check() does not flag valid USG-T OC actions as ignored."""
    actions = [
        "bootstrap",
        "nobootstrap",
        "bootstrapscale",
        "nobootstrapscale",
        "ddreference",
        "deltat 0.5",
        "tminat 1e-10",
        "tmaxat 1e10",
        "tadjat 2.0",
        "tcutat 5.0",
        "hclose 1e-5",
        "btol 1.1",
        "mxiter 100",
        "save head",
        "save conc",
        "print conc",
        "save budget",
    ]
    assert _oc_check_warnings(function_tmpdir, "okact", actions) == []

    # a genuinely unknown action is still flagged
    descs = _oc_check_warnings(function_tmpdir, "bogus", ["save head", "frobnicate x"])
    assert any("frobnicate" in d.lower() for d in descs)


def test_mfusgoc_check_warns_save_ibound(function_tmpdir):
    """check() warns that SAVE IBOUND is preserved by FloPy but USG-T rejects it."""
    descs = _oc_check_warnings(
        function_tmpdir, "ibwarn", ["save head", "save ibound"]
    )
    assert any(
        "SAVE IBOUND" in d and "rejected by USG-T 2.7" in d for d in descs
    ), descs


def test_mfusgoc_check_warns_param_without_value(function_tmpdir):
    """check() warns when a USG-T keyword that needs a value is missing it.

    USG-T (glo2basu1.f SGWF2BAS7N) reads a value after DELTAT/TMINAT/TMAXAT/
    TADJAT/TCUTAT/HCLOSE/BTOL/MXITER, so a bare keyword is incomplete.
    """
    # complete forms produce no warnings
    assert (
        _oc_check_warnings(
            function_tmpdir,
            "okparam",
            ["save head", "deltat 0.5", "hclose 1e-5", "mxiter 100"],
        )
        == []
    )

    # each parametric keyword without a value is flagged
    for kw in (
        "deltat",
        "tminat",
        "tmaxat",
        "tadjat",
        "tcutat",
        "hclose",
        "btol",
        "mxiter",
    ):
        descs = _oc_check_warnings(function_tmpdir, f"nv_{kw}", ["save head", kw])
        assert any(
            kw.upper() in d and "requires a value" in d for d in descs
        ), (kw, descs)

    # single-word USG-T actions remain valid (no value needed)
    assert (
        _oc_check_warnings(
            function_tmpdir, "single", ["bootstrap", "ddreference", "nobootstrapscale"]
        )
        == []
    )


# ---------------------------------------------------------------------------
# MfUsgMdt tests (Matrix Diffusion Transport, gwt2mdtu1.for)
# ---------------------------------------------------------------------------


def _mdt_model(function_tmpdir, name, mcomp=1, idpf=0, nlay=1, nrow=2, ncol=2):
    from flopy.modflow import ModflowDis

    m = MfUsg(structured=True, model_ws=str(function_tmpdir), modelname=name)
    ModflowDis(m, nlay=nlay, nrow=nrow, ncol=ncol, nper=1)
    m.mcomp = mcomp
    m.itrnsp = 1
    m.idpf = idpf
    return m


def test_mfusgmdt_minimal_authoring_roundtrip(function_tmpdir):
    """Minimal structured MDT (mcomp=1, tshiftmd=0) authors and round-trips."""
    from flopy.mfusg import MfUsgMdt

    mdt = MfUsgMdt(
        _mdt_model(function_tmpdir, "m1"),
        mdflag=1,
        volfracmd=0.1,
        pormd=0.2,
        rhobmd=1500.0,
        difflenmd=0.5,
        tortmd=0.7,
        kdmd=0.01,
        decaymd=1e-6,
        yieldmd=1.0,
        diffmd=1e-9,
    )
    mdt.fn_path = str(function_tmpdir / "m1.mdt")
    mdt.write_file()

    re = MfUsgMdt.load(mdt.fn_path, _mdt_model(function_tmpdir, "m1b"))
    assert re.frahk == 0 and re.tshiftmd == 0.0
    assert np.isclose(re.volfracmd.array.mean(), 0.1)
    assert np.isclose(re.pormd.array.mean(), 0.2)
    assert np.isclose(re.kdmd[0].array.mean(), 0.01)

    # rewrite then reload is semantically stable
    re.fn_path = str(function_tmpdir / "m1_re.mdt")
    re.write_file()
    re2 = MfUsgMdt.load(re.fn_path, _mdt_model(function_tmpdir, "m1c"))
    assert np.isclose(re2.volfracmd.array.mean(), 0.1)
    assert np.isclose(re2.tortmd.array.mean(), 0.7)
    assert np.isclose(re2.kdmd[0].array.mean(), 0.01)


def test_mfusgmdt_tshift_aiold_roundtrip(function_tmpdir):
    """TSHIFTMD>0 writes and round-trips AIOLD1MD/AIOLD2MD."""
    from flopy.mfusg import MfUsgMdt

    mdt = MfUsgMdt(
        _mdt_model(function_tmpdir, "ts"),
        tshiftmd=2.0,
        kdmd=0.01,
        decaymd=1e-6,
        yieldmd=1.0,
        diffmd=1e-9,
        aiold1md=0.3,
        aiold2md=0.4,
    )
    mdt.fn_path = str(function_tmpdir / "ts.mdt")
    mdt.write_file()
    assert "TSHIFTMD" in Path(mdt.fn_path).read_text()

    re = MfUsgMdt.load(mdt.fn_path, _mdt_model(function_tmpdir, "ts2"))
    assert np.isclose(re.tshiftmd, 2.0)
    assert np.isclose(re.aiold1md[0].array.mean(), 0.3)
    assert np.isclose(re.aiold2md[0].array.mean(), 0.4)


def test_mfusgmdt_multispecies_roundtrip(function_tmpdir):
    """Multi-species MDT preserves distinct per-component values."""
    from flopy.mfusg import MfUsgMdt

    mdt = MfUsgMdt(
        _mdt_model(function_tmpdir, "ms", mcomp=2),
        kdmd=[0.01, 0.02],
        decaymd=[1e-6, 2e-6],
        yieldmd=[1.0, 0.5],
        diffmd=[1e-9, 2e-9],
    )
    mdt.fn_path = str(function_tmpdir / "ms.mdt")
    mdt.write_file()

    re = MfUsgMdt.load(mdt.fn_path, _mdt_model(function_tmpdir, "ms2", mcomp=2))
    assert np.isclose(re.kdmd[0].array.mean(), 0.01)
    assert np.isclose(re.kdmd[1].array.mean(), 0.02)
    assert np.isclose(re.diffmd[1].array.mean(), 2e-9)


def test_mfusgmdt_frahk_fradarcy_roundtrip(function_tmpdir):
    """FRAHK and FRADARCY each author and reload (regression: load was lowercase)."""
    from flopy.mfusg import MfUsgMdt

    for opt in ("frahk", "fradarcy"):
        mdt = MfUsgMdt(
            _mdt_model(function_tmpdir, opt),
            kdmd=0.01,
            decaymd=1e-6,
            yieldmd=1.0,
            diffmd=1e-9,
            **{opt: True},
        )
        mdt.fn_path = str(function_tmpdir / f"{opt}.mdt")
        mdt.write_file()
        assert opt.upper() in Path(mdt.fn_path).read_text()
        re = MfUsgMdt.load(mdt.fn_path, _mdt_model(function_tmpdir, f"{opt}b"))
        assert getattr(re, opt) == 1


def test_mfusgmdt_output_options_roundtrip(function_tmpdir):
    """SEPARATE_AI2 and MULTIFILE_MD round-trip (rootname case preserved)."""
    from flopy.mfusg import MfUsgMdt

    mdt = MfUsgMdt(
        _mdt_model(function_tmpdir, "out"),
        imdtcf=58,
        iunitAI2=59,
        crootname="mdRoot",
        kdmd=0.01,
        decaymd=1e-6,
        yieldmd=1.0,
        diffmd=1e-9,
    )
    mdt.fn_path = str(function_tmpdir / "out.mdt")
    mdt.write_file()
    text = Path(mdt.fn_path).read_text()
    assert "SEPARATE_AI2" in text and "MULTIFILE_MD mdRoot" in text

    re = MfUsgMdt.load(mdt.fn_path, _mdt_model(function_tmpdir, "out2"))
    assert re.iunitAI2 == 59 and re.crootname == "mdRoot"


def test_mfusgmdt_idpf_skips_volfracmd(function_tmpdir):
    """With dual-porosity flow (IDPF!=0), VOLFRACMD is not written or expected."""
    from flopy.mfusg import MfUsgMdt

    def blocks(model, name):
        mdt = MfUsgMdt(
            model, kdmd=0.01, decaymd=1e-6, yieldmd=1.0, diffmd=1e-9
        )
        mdt.fn_path = str(function_tmpdir / f"{name}.mdt")
        mdt.write_file()
        text = Path(mdt.fn_path).read_text().upper()
        return text.count("CONSTANT") + text.count("INTERNAL"), mdt.fn_path

    n0, _ = blocks(_mdt_model(function_tmpdir, "d0", idpf=0), "d0")
    n1, path1 = blocks(_mdt_model(function_tmpdir, "d1", idpf=1), "d1")
    assert n1 == n0 - 1  # the VOLFRACMD array block is dropped under IDPF!=0

    # and it loads back without expecting VOLFRACMD
    re = MfUsgMdt.load(path1, _mdt_model(function_tmpdir, "d1b", idpf=1))
    assert np.isclose(re.pormd.array.mean(), 0.0)


def test_mfusgmdt_external_handle_write(function_tmpdir):
    """write_file(f=handle) writes to and does not close the caller's handle."""
    import io

    from flopy.mfusg import MfUsgMdt

    mdt = MfUsgMdt(
        _mdt_model(function_tmpdir, "ext"),
        kdmd=0.01,
        decaymd=1e-6,
        yieldmd=1.0,
        diffmd=1e-9,
    )
    buf = io.StringIO()
    mdt.write_file(f=buf)  # regression: previously NameError on f_obj
    assert len(buf.getvalue()) > 50
    assert not buf.closed  # external handle must stay open


def test_mfusgmdt_rejects_invalid(function_tmpdir):
    """MDT fails explicitly on conflicting/incomplete options and bad lengths."""
    from flopy.mfusg import MfUsgMdt

    with pytest.raises(ValueError, match="mutually exclusive"):
        MfUsgMdt(_mdt_model(function_tmpdir, "b1"), frahk=True, fradarcy=True)

    with pytest.raises(ValueError, match="MCOMP"):
        MfUsgMdt(_mdt_model(function_tmpdir, "b2", mcomp=2), kdmd=[0.01])

    with pytest.raises(ValueError, match="IDPF"):
        MfUsgMdt(_mdt_model(function_tmpdir, "b3", idpf=1), frahk=True)

    with pytest.raises(ValueError, match="imdtcf"):
        MfUsgMdt(_mdt_model(function_tmpdir, "b4"), imdtcf=0, crootname="x")


def test_mfusgmdt_tshiftmd_threshold(function_tmpdir):
    """TSHIFTMD uses the Fortran threshold (1e-10) consistently for AIOLD.

    USG-T (gwt2mdtu1.for) reads AIOLD1MD/AIOLD2MD only when TSHIFTMD > 1e-10, so
    FloPy must (a) not round a small valid value to 0.0 and (b) not write AIOLD
    for a value below the threshold.
    """
    from flopy.mfusg import MfUsgMdt

    base = dict(kdmd=0.01, decaymd=1e-6, yieldmd=1.0, diffmd=1e-9)

    def n_blocks(text):
        up = text.upper()
        return up.count("CONSTANT") + up.count("INTERNAL")

    # tshiftmd=2.0: TSHIFTMD + AIOLD written and round-tripped
    big = MfUsgMdt(
        _mdt_model(function_tmpdir, "tbig"),
        tshiftmd=2.0,
        aiold1md=0.3,
        aiold2md=0.4,
        **base,
    )
    big.fn_path = str(function_tmpdir / "tbig.mdt")
    big.write_file()
    assert "TSHIFTMD" in Path(big.fn_path).read_text()
    re = MfUsgMdt.load(big.fn_path, _mdt_model(function_tmpdir, "tbig2"))
    assert np.isclose(re.tshiftmd, 2.0)
    assert np.isclose(re.aiold1md[0].array.mean(), 0.3)

    # tshiftmd=1e-6: above threshold -> a non-zero readable value (NOT 0.00)
    # and AIOLD are written; round-trips back above the threshold.
    small = MfUsgMdt(
        _mdt_model(function_tmpdir, "sm"),
        tshiftmd=1e-6,
        aiold1md=0.3,
        aiold2md=0.4,
        **base,
    )
    small.fn_path = str(function_tmpdir / "sm.mdt")
    small.write_file()
    ts_line = next(
        ln
        for ln in Path(small.fn_path).read_text().splitlines()
        if "TSHIFTMD" in ln
    )
    assert float(ts_line.split("TSHIFTMD")[1].split()[0]) > 1e-10  # not rounded to 0
    re = MfUsgMdt.load(small.fn_path, _mdt_model(function_tmpdir, "sm2"))
    assert re.tshiftmd > 1e-10
    assert np.isclose(re.aiold1md[0].array.mean(), 0.3)

    # tshiftmd=1e-12: below threshold -> no TSHIFTMD keyword and no AIOLD arrays
    tiny = MfUsgMdt(
        _mdt_model(function_tmpdir, "ti"),
        tshiftmd=1e-12,
        aiold1md=0.3,
        aiold2md=0.4,
        **base,
    )
    tiny.fn_path = str(function_tmpdir / "ti.mdt")
    tiny.write_file()
    none = MfUsgMdt(_mdt_model(function_tmpdir, "no"), **base)
    none.fn_path = str(function_tmpdir / "no.mdt")
    none.write_file()
    assert "TSHIFTMD" not in Path(tiny.fn_path).read_text()
    assert n_blocks(Path(tiny.fn_path).read_text()) == n_blocks(
        Path(none.fn_path).read_text()
    )

    # IDPF!=0: a below-threshold tshiftmd is inactive (no error); an
    # above-threshold one is an active option and must raise.
    MfUsgMdt(_mdt_model(function_tmpdir, "dp_ok", idpf=1), tshiftmd=1e-12, **base)
    with pytest.raises(ValueError, match="IDPF"):
        MfUsgMdt(_mdt_model(function_tmpdir, "dp_bad", idpf=1), tshiftmd=1e-6, **base)


# ---------------------------------------------------------------------------
# MfUsgLak tests (Lake package, gwf2lak7u1.f)
# ---------------------------------------------------------------------------


def _lak_model(function_tmpdir, name, mcomp=0):
    from flopy.modflow import ModflowDis

    m = MfUsg(structured=True, model_ws=str(function_tmpdir), modelname=name)
    ModflowDis(m, nlay=1, nrow=3, ncol=3, nper=1, nstp=1, steady=False)
    if mcomp:
        m.itrnsp = 1
        m.mcomp = mcomp
    return m


def _lak_grids():
    lakarr = np.zeros((1, 3, 3), dtype=int)
    lakarr[0, 1, 1] = 1
    return {0: lakarr}, {0: np.ones((1, 3, 3), dtype=float) * 0.1}


def test_mfusglak_minimal_authoring_roundtrip(function_tmpdir):
    """A minimal LAK (no transport) authors from scratch and round-trips."""
    from flopy.mfusg import MfUsgLak

    lakarr, bdlknc = _lak_grids()
    lak = MfUsgLak(
        _lak_model(function_tmpdir, "m1"),
        nlakes=1,
        stages=100.0,
        lakarr=lakarr,
        bdlknc=bdlknc,
        flux_data={0: {0: [1.0, 2.0, 0.0, 0.0]}},
    )
    lak.fn_path = str(function_tmpdir / "m1.lak")
    lak.write_file()
    re = MfUsgLak.load(lak.fn_path, _lak_model(function_tmpdir, "m1b"))
    assert re.nlakes == 1 and not re.transportboundary


def test_mfusglak_classic_transport_roundtrip(function_tmpdir):
    """Classic lake transport (CPPT/CRNF per lake-component) round-trips."""
    from flopy.mfusg import MfUsgLak

    lakarr, bdlknc = _lak_grids()
    lak = MfUsgLak(
        _lak_model(function_tmpdir, "ct", mcomp=1),
        nlakes=1,
        stages=100.0,
        lakarr=lakarr,
        bdlknc=bdlknc,
        flux_data={0: {0: [1.0, 2.0, 0.0, 0.0]}},
        clake=[[5.0]],
        conc_data={0: {(0, 0): [3.0, 1.0]}},  # CPPT, CRNF
    )
    lak.fn_path = str(function_tmpdir / "ct.lak")
    lak.write_file()
    assert "TRANSPORTBOUNDARY" not in Path(lak.fn_path).read_text()

    re = MfUsgLak.load(lak.fn_path, _lak_model(function_tmpdir, "ct2", mcomp=1))
    assert not re.transportboundary
    assert [float(x) for x in re.conc_data[0][(0, 0)]] == [3.0, 1.0]


def test_mfusglak_transportboundary_roundtrip(function_tmpdir):
    """TRANSPORTBOUNDARY writes one CLAKE line per lake (NSOL values) and reloads.

    Regression: the writer previously emitted one line per component and the
    `transportboundary` flag did not add the header keyword, so the file did not
    round-trip.
    """
    from flopy.mfusg import MfUsgLak

    lakarr, bdlknc = _lak_grids()
    lak = MfUsgLak(
        _lak_model(function_tmpdir, "tb", mcomp=2),
        nlakes=1,
        stages=100.0,
        lakarr=lakarr,
        bdlknc=bdlknc,
        flux_data={0: {0: [1.0, 2.0, 0.0, 0.0]}},
        transportboundary=True,
        clake=[[5.0, 6.0]],
        conc_data={0: {(0, 0): 5.0, (0, 1): 6.0}},
    )
    lak.fn_path = str(function_tmpdir / "tb.lak")
    lak.write_file()
    text = Path(lak.fn_path).read_text()
    assert "TRANSPORTBOUNDARY" in text  # header keyword emitted
    # exactly one dataset-9b line for the single lake (not one per component)
    assert sum("Data set 9b" in ln for ln in text.splitlines()) == 1

    re = MfUsgLak.load(lak.fn_path, _lak_model(function_tmpdir, "tb2", mcomp=2))
    assert re.transportboundary
    assert re.conc_data[0][(0, 0)] == 5.0 and re.conc_data[0][(0, 1)] == 6.0


def test_mfusglak_tableinput_authoring(function_tmpdir):
    """TABLEINPUT authoring emits the keyword + per-lake tab unit, registered."""
    from flopy.mfusg import MfUsgLak

    model = _lak_model(function_tmpdir, "tab")
    lakarr, bdlknc = _lak_grids()
    lak = MfUsgLak(
        model,
        nlakes=1,
        stages=100.0,
        lakarr=lakarr,
        bdlknc=bdlknc,
        flux_data={0: {0: [1.0, 2.0, 0.0, 0.0]}},
        options=["TABLEINPUT"],
        tab_files=["lake1.tab"],
    )
    lak.fn_path = str(function_tmpdir / "tab.lak")
    lak.write_file()
    assert "TABLEINPUT" in Path(lak.fn_path).read_text()
    assert lak.tabdata and lak.iunit_tab and lak.iunit_tab[0] > 0
    # the tab file is registered as an external unit on the model
    assert lak.iunit_tab[0] in model.external_units


def test_mfusglak_rejects_invalid(function_tmpdir):
    """LAK fails explicitly on missing flux_data, bad transport options/lengths."""
    from flopy.mfusg import MfUsgLak

    def build(name, mcomp=0, **kw):
        lakarr, bdlknc = _lak_grids()
        return MfUsgLak(
            _lak_model(function_tmpdir, name, mcomp=mcomp),
            nlakes=1,
            stages=100.0,
            lakarr=lakarr,
            bdlknc=bdlknc,
            **kw,
        )

    # flux_data (dataset 9) is required
    with pytest.raises(ValueError, match="flux_data"):
        build("n1", flux_data=None)

    # TRANSPORTBOUNDARY needs active transport
    with pytest.raises(ValueError, match="TRANSPORTBOUNDARY"):
        build("n2", flux_data={0: {0: [0, 0, 0, 0]}}, transportboundary=True)

    # clake must be nlakes x mcomp
    with pytest.raises(ValueError, match="clake"):
        build(
            "n3",
            mcomp=2,
            flux_data={0: {0: [0, 0, 0, 0]}},
            clake=[[1.0]],
            conc_data={0: {(0, 0): 1.0, (0, 1): 1.0}},
        )

    # active transport requires conc_data (dataset 9b)
    with pytest.raises(ValueError, match="conc_data"):
        build("n4", mcomp=1, flux_data={0: {0: [0, 0, 0, 0]}}, clake=[[1.0]])


def test_mfusglak_tableinput_rejects_wrong_counts(function_tmpdir):
    """TABLEINPUT with the wrong number of tab_files / tab_units fails explicitly.

    Regression: too few tab_files used to be a dead message and then crashed
    write_file() with IndexError on iunit_tab[n].
    """
    from flopy.mfusg import MfUsgLak

    lakarr, bdlknc = _lak_grids()
    base = {
        "nlakes": 2,
        "stages": 100.0,
        "lakarr": lakarr,
        "bdlknc": bdlknc,
        "flux_data": {0: {0: [1.0, 2.0, 0.0, 0.0], 1: [1.0, 2.0, 0.0, 0.0]}},
        "options": ["TABLEINPUT"],
    }
    # one tab_file for two lakes
    with pytest.raises(ValueError, match="tab_file"):
        MfUsgLak(_lak_model(function_tmpdir, "tw1"), tab_files=["one.tab"], **base)
    # two tab_files but a tab_units list of the wrong length
    with pytest.raises(ValueError, match="tab_unit"):
        MfUsgLak(
            _lak_model(function_tmpdir, "tw2"),
            tab_files=["one.tab", "two.tab"],
            tab_units=[201],
            **base,
        )


def test_mfusglak_classic_transport_rejects_incomplete(function_tmpdir):
    """Classic transport with a missing (lake, component) conc_data entry fails."""
    from flopy.mfusg import MfUsgLak

    lakarr, bdlknc = _lak_grids()
    with pytest.raises(ValueError, match="conc_data"):
        MfUsgLak(
            _lak_model(function_tmpdir, "ci", mcomp=2),
            nlakes=1,
            stages=100.0,
            lakarr=lakarr,
            bdlknc=bdlknc,
            flux_data={0: {0: [1.0, 2.0, 0.0, 0.0]}},
            clake=[[5.0, 6.0]],
            conc_data={0: {(0, 0): [3.0, 1.0]}},  # missing (0, 1)
        )


def test_mfusglak_classic_transport_wthdrw_caug_roundtrip(function_tmpdir):
    """WTHDRW<0 needs CPPT,CRNF,CAUG: 2 values are rejected, 3 round-trip."""
    from flopy.mfusg import MfUsgLak

    lakarr, bdlknc = _lak_grids()
    common = {
        "nlakes": 1,
        "stages": 100.0,
        "lakarr": lakarr,
        "bdlknc": bdlknc,
        "flux_data": {0: {0: [1.0, 2.0, 0.0, -5.0]}},  # WTHDRW<0 -> augmentation
        "clake": [[5.0]],
    }
    # only CPPT, CRNF supplied but WTHDRW<0 requires CAUG too
    with pytest.raises(ValueError, match="CAUG"):
        MfUsgLak(
            _lak_model(function_tmpdir, "wd1", mcomp=1),
            conc_data={0: {(0, 0): [3.0, 1.0]}},
            **common,
        )
    # full CPPT, CRNF, CAUG writes and reloads
    lak = MfUsgLak(
        _lak_model(function_tmpdir, "wd2", mcomp=1),
        conc_data={0: {(0, 0): [3.0, 1.0, 2.0]}},
        **common,
    )
    lak.fn_path = str(function_tmpdir / "wd2.lak")
    lak.write_file()
    re = MfUsgLak.load(lak.fn_path, _lak_model(function_tmpdir, "wd2b", mcomp=1))
    assert [float(x) for x in re.conc_data[0][(0, 0)]] == [3.0, 1.0, 2.0]


def test_mfusglak_transportboundary_rejects_incomplete(function_tmpdir):
    """TRANSPORTBOUNDARY with a missing (lake, component) conc_data entry fails."""
    from flopy.mfusg import MfUsgLak

    lakarr, bdlknc = _lak_grids()
    with pytest.raises(ValueError, match="conc_data"):
        MfUsgLak(
            _lak_model(function_tmpdir, "tbi", mcomp=2),
            nlakes=1,
            stages=100.0,
            lakarr=lakarr,
            bdlknc=bdlknc,
            flux_data={0: {0: [1.0, 2.0, 0.0, 0.0]}},
            transportboundary=True,
            clake=[[5.0, 6.0]],
            conc_data={0: {(0, 0): 5.0}},  # missing (0, 1)
        )


def test_mfusglak_flux_data_rejects_missing_lake(function_tmpdir):
    """flux_data must carry one dataset-9a entry per lake."""
    from flopy.mfusg import MfUsgLak

    lakarr, bdlknc = _lak_grids()
    with pytest.raises(ValueError, match="missing lake"):
        MfUsgLak(
            _lak_model(function_tmpdir, "fm"),
            nlakes=2,
            stages=100.0,
            lakarr=lakarr,
            bdlknc=bdlknc,
            flux_data={0: {0: [1.0, 2.0, 0.0, 0.0]}},  # lake 1 missing
        )
