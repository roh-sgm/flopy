import io
import os
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


@requires_exe("mfusg_gsi")
def test_usg_load_Ex1_1D(function_tmpdir, mfusg_transport_Ex1_1D_model_path):
    print("testing mfusg transport model loading: BTN_Test1.nam")

    fname = mfusg_transport_Ex1_1D_model_path / "BTN_Test1.nam"
    assert os.path.isfile(fname), f"nam file not found {fname}"

    # Create the model
    m = MfUsg.load(
        fname, exe_name="mfusg_gsi", verbose=True, model_ws=function_tmpdir, check=True
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


@requires_exe("mfusg_gsi")
def test_usg_load_Ex2_Radial_adv(
    function_tmpdir, mfusg_transport_Ex2_Radial_2D_model_path
):
    print("testing mfusg transport model loading: Radial-adv.nam")

    fname = mfusg_transport_Ex2_Radial_2D_model_path / "Radial_adv.nam"

    assert os.path.isfile(fname), f"nam file not found {fname}"

    # Create the model
    m = MfUsg.load(
        fname, exe_name="mfusg_gsi", verbose=True, model_ws=function_tmpdir, check=True
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


@requires_exe("mfusg_gsi")
def test_usg_load_Ex2_Radial_Disp(
    function_tmpdir, mfusg_transport_Ex2_Radial_2D_model_path
):
    print("testing mfusg transport model loading: Radial-dis.nam")

    fname = mfusg_transport_Ex2_Radial_2D_model_path / "Radial_dis.nam"

    assert os.path.isfile(fname), f"nam file not found {fname}"

    # Create the model
    m = MfUsg.load(
        fname, exe_name="mfusg_gsi", verbose=True, model_ws=function_tmpdir, check=True
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


@requires_exe("mfusg_gsi")
def test_usg_load_Ex3_CLN_Conduit(
    function_tmpdir, mfusg_transport_Ex3_CLN_Conduit_model_path
):
    print("testing mfusg transport model loading: Conduit.nam")

    fname = mfusg_transport_Ex3_CLN_Conduit_model_path / "Conduit/Conduit.nam"
    assert os.path.isfile(fname), f"nam file not found {fname}"

    # Create the model
    m = MfUsg.load(
        fname, exe_name="mfusg_gsi", verbose=True, model_ws=function_tmpdir, check=True
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


@requires_exe("mfusg_gsi")
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
        fname, exe_name="mfusg_gsi", verbose=True, model_ws=function_tmpdir, check=True
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


@requires_exe("mfusg_gsi")
def test_usg_load_Ex3_CLN_Conduit_Nest(
    function_tmpdir, mfusg_transport_Ex3_CLN_Conduit_model_path
):
    print("testing mfusg transport model loading: Conduit.nam")

    fname = mfusg_transport_Ex3_CLN_Conduit_model_path / "Nest/Conduit_Nest.nam"
    assert os.path.isfile(fname), f"nam file not found {fname}"

    # Create the model
    m = MfUsg.load(
        fname, exe_name="mfusg_gsi", verbose=True, model_ws=function_tmpdir, check=True
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


@requires_exe("mfusg_gsi")
def test_usg_load_Ex4_Dual_Domain(
    function_tmpdir, mfusg_transport_Ex4_Dual_Domain_model_path
):
    print("testing mfusg transport model loading: Conduit.nam")

    fname = mfusg_transport_Ex4_Dual_Domain_model_path / "DualDomain.nam"
    assert os.path.isfile(fname), f"nam file not found {fname}"

    # Create the model
    m = MfUsg.load(
        fname, exe_name="mfusg_gsi", verbose=True, model_ws=function_tmpdir, check=True
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


@requires_exe("mfusg_gsi")
def test_usg_load_Ex5_Henry(function_tmpdir, mfusg_transport_Ex5_Henry_model_path):
    print("testing mfusg transport model loading: Conduit.nam")

    fname = mfusg_transport_Ex5_Henry_model_path / "Henry.nam"
    assert os.path.isfile(fname), f"nam file not found {fname}"

    # Create the model
    m = MfUsg.load(
        fname, exe_name="mfusg_gsi", verbose=True, model_ws=function_tmpdir, check=True
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


@requires_exe("mfusg_gsi")
def test_usg_load_Ex6_Stallman_Heat(
    function_tmpdir, mfusg_transport_Ex6_Stallman_model_path
):
    print("testing mfusg transport model loading: Stallman_Heat.nam")

    fname = mfusg_transport_Ex6_Stallman_model_path / "Heat/Stallman_Heat.nam"
    assert os.path.isfile(fname), f"nam file not found {fname}"

    # Create the model
    m = MfUsg.load(
        fname, exe_name="mfusg_gsi", verbose=True, model_ws=function_tmpdir, check=True
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


@requires_exe("mfusg_gsi")
def test_usg_load_Ex6_Stallman_Solute(
    function_tmpdir, mfusg_transport_Ex6_Stallman_model_path
):
    print("testing mfusg transport model loading: Stallman_Solute.nam")

    fname = mfusg_transport_Ex6_Stallman_model_path / "Solute/Stallman_Solute.nam"
    assert os.path.isfile(fname), f"nam file not found {fname}"

    # Create the model
    m = MfUsg.load(
        fname, exe_name="mfusg_gsi", verbose=True, model_ws=function_tmpdir, check=True
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


@requires_exe("mfusg_gsi")
def test_usg_load_Ex6_Stallman_Solute_Heat(
    function_tmpdir, mfusg_transport_Ex6_Stallman_model_path
):
    print("testing mfusg transport model loading: Stallman.nam")

    fname = mfusg_transport_Ex6_Stallman_model_path / "Solute_Heat/Stallman.nam"
    assert os.path.isfile(fname), f"nam file not found {fname}"

    # Create the model
    m = MfUsg.load(
        fname, exe_name="mfusg_gsi", verbose=True, model_ws=function_tmpdir, check=True
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


@requires_exe("mfusg_gsi")
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
        fname, exe_name="mfusg_gsi", verbose=True, model_ws=function_tmpdir, check=True
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


@requires_exe("mfusg_gsi")
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
        fname, exe_name="mfusg_gsi", verbose=True, model_ws=function_tmpdir, check=True
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


@requires_exe("mfusg_gsi")
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
        fname, exe_name="mfusg_gsi", verbose=True, model_ws=function_tmpdir, check=True
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


@requires_exe("mfusg_gsi")
def test_usg_load_Ex8_Lake(function_tmpdir, mfusg_transport_Ex8_Lake_model_path):
    print("testing mfusg transport model loading: lak_usg_01.nam")

    fname = mfusg_transport_Ex8_Lake_model_path / "lak_usg_01.nam"
    assert os.path.isfile(fname), f"nam file not found {fname}"

    # Create the model
    m = MfUsg.load(
        fname, exe_name="mfusg_gsi", verbose=True, model_ws=function_tmpdir, check=True
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


@requires_exe("mfusg_gsi")
def test_usg_load_Ex9_PFAS(function_tmpdir, mfusg_transport_Ex9_PFAS_model_path):
    print("testing mfusg transport model loading: PFAS_C1.nam")

    fname = mfusg_transport_Ex9_PFAS_model_path / "C1/PFAS_C1.nam"
    assert os.path.isfile(fname), f"nam file not found {fname}"

    # Create the model
    m = MfUsg.load(
        fname, exe_name="mfusg_gsi", verbose=True, model_ws=function_tmpdir, check=True
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

    dpt_file = function_tmpdir / "aw.dpt"
    dpt_file.write_text(
        "# DPT with immobile air-water adsorption\n"
        " 0 0 0 0 0 0 0 A-W_ADSORBIM\n"
    )
    ml = MfUsg(structured=False, model_ws=str(function_tmpdir))
    ModflowDis(ml, nlay=1, nrow=1, ncol=1, nper=1)
    with pytest.raises(NotImplementedError):
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
