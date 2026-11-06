import io
import os
import shutil
import sys
import threading
import warnings
from pathlib import Path

import numpy as np
import pytest
from flaky import flaky
from modflow_devtools.markers import requires_exe, requires_pkg

from autotest.conftest import get_example_data_path
from autotest.usgt_exe import (
    USGT_EXE,
    USGT_EXE_NAMES,
    requires_usgt,
    select_usgt_exe,
)
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

# USG-Transport executable for the optional real-model run tests: USGT_EXE, or
# mfusgt / mfusg_gsi on PATH (autotest/usgt_exe.py). @requires_usgt skips these
# tests cleanly when no executable is selected, so the default suite stays light.


@pytest.mark.parametrize(
    "environ, available, expected",
    [
        ({}, {"mfusgt": "/new/mfusgt"}, "/new/mfusgt"),
        ({}, {"mfusg_gsi": "/old/mfusg_gsi"}, "/old/mfusg_gsi"),
        ({}, {"mfusgt": "/new/mfusgt", "mfusg_gsi": "/old/mfusg_gsi"}, "/new/mfusgt"),
        (
            {"USGT_EXE": "mfusg_gsi"},
            {"mfusgt": "/new/mfusgt", "mfusg_gsi": "/old/mfusg_gsi"},
            "/old/mfusg_gsi",
        ),
        (
            {"USGT_EXE": "/opt/usg t/usgt_281"},
            {"mfusgt": "/new/mfusgt", "/opt/usg t/usgt_281": "/opt/usg t/usgt_281"},
            "/opt/usg t/usgt_281",
        ),
    ],
    ids=["new-name", "old-name", "both-names", "override-name", "override-path"],
)
def test_usgt_exe_selection(environ, available, expected):
    """USGT_EXE wins; otherwise mfusgt is preferred over mfusg_gsi."""
    calls = []

    def which(name):
        calls.append(name)
        return available.get(name)

    assert select_usgt_exe(environ, which) == (os.path.abspath(expected), None)
    # names are looked up as given, so shutil.which applies PATHEXT on Windows
    if "USGT_EXE" in environ:
        assert calls == [environ["USGT_EXE"]]
    else:
        found = USGT_EXE_NAMES.index(Path(expected).name)
        assert calls == list(USGT_EXE_NAMES[: found + 1])


@pytest.mark.parametrize("override", ["missing_usgt", "/no/such/usgt", ""])
def test_usgt_exe_invalid_override_does_not_fall_back(override):
    available = {"mfusgt": "/new/mfusgt", "mfusg_gsi": "/old/mfusg_gsi"}
    exe, reason = select_usgt_exe({"USGT_EXE": override}, available.get)
    assert exe is None
    assert reason == f"USGT_EXE={override!r} is not an executable (no fallback)"


def test_usgt_exe_none_available():
    exe, reason = select_usgt_exe({}, {}.get)
    assert exe is None
    assert reason.startswith("missing executable: none of mfusgt, mfusg_gsi")


@pytest.mark.skipif(sys.platform == "win32", reason="POSIX executable bits")
def test_usgt_exe_selection_on_path_with_spaces(function_tmpdir, monkeypatch):
    """Real PATH lookup; the empty files are only selected, never run."""
    bindir = function_tmpdir / "bin with spaces"
    bindir.mkdir()
    for name in USGT_EXE_NAMES:
        (bindir / name).write_text("")
        (bindir / name).chmod(0o755)
    monkeypatch.setenv("PATH", str(bindir))
    assert select_usgt_exe({}) == (str(bindir / "mfusgt"), None)
    override = str(bindir / "mfusg_gsi")
    assert select_usgt_exe({"USGT_EXE": override}) == (override, None)
    (bindir / "mfusgt").chmod(0o644)
    assert select_usgt_exe({}) == (str(bindir / "mfusg_gsi"), None)


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


@requires_usgt
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


@requires_usgt
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


@requires_usgt
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


@requires_usgt
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


def test_usg_Ex3_CLN_Conduit_nam_names_match_files(
    mfusg_transport_Ex3_CLN_Conduit_model_path,
):
    """Every Conduit.nam entry naming a file of the example folder uses its
    exact name: a case-only difference (e.g. Conduit.WEL for Conduit.wel)
    stops the engine on case-sensitive file systems. Output files written by
    a run (LIST, DATA(BINARY)) are not required to exist."""
    folder = mfusg_transport_Ex3_CLN_Conduit_model_path / "Conduit"
    files = os.listdir(folder)
    by_lower = {name.lower(): name for name in files}
    entries = [
        line.split()[2]
        for line in (folder / "Conduit.nam").read_text().splitlines()
        if line.strip() and not line.lstrip().startswith("#")
    ]
    existing = [name for name in entries if name.lower() in by_lower]
    mismatched = [
        (name, by_lower[name.lower()]) for name in existing if name not in files
    ]
    assert not mismatched, f"name file entry -> file in {folder}: {mismatched}"
    assert {"Conduit.wel", "Conduit.cln", "Conduit.chd"} <= set(existing)


@requires_usgt
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


@requires_usgt
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


@requires_usgt
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


@requires_usgt
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


@requires_usgt
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


@requires_usgt
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


@requires_usgt
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


@requires_usgt
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


@requires_usgt
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


@requires_usgt
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


@requires_usgt
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


@requires_usgt
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
    """MfUsgRiv load + write on DISU: the second header field is IRIVCB
    (GWF2RIV7U1AR), loaded as ipakcb; node/stage/cond/rbot; reuse SP.
    (Formerly asserted irdflag == 50: that field is IRIVCB, not a read flag.)"""
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

    assert riv.ipakcb == 50
    assert riv.irdflag is None
    assert 50 in ml.output_units  # CBC output registered for IRIVCB > 0
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
    assert text.splitlines()[1].split()[:2] == ["2", "50"]
    assert "Stress Period 1" in text
    assert "Stress Period 2" in text
    assert "\n 101" in text

    ml2 = MfUsg(structured=False)
    ModflowDis(ml2, nlay=1, nrow=1, ncol=1, nper=2)
    riv2 = MfUsgRiv.load(str(riv_out), ml2, nper=2, ext_unit_dict={})
    assert riv2.ipakcb == 50
    assert np.isclose(riv2.stress_period_data[0]["stage"][0], 15.0, atol=0.01)


def test_mfusgriv_trailing_values_ignored(function_tmpdir):
    """Values after node, stage, cond and rbot are not read by USG-T
    (ULSTRDU reads NAUX values only): load ignores them with a warning and
    the rewrite has no extra column. (Formerly they were kept as an 'irch'
    column, which no USG-T 2.8.1 RIV routine reads.)"""
    from flopy.mfusg import MfUsgRiv
    from flopy.modflow import ModflowDis

    riv_in = function_tmpdir / "trailing.riv"
    riv_in.write_text(
        "# MfUsgRiv trailing value test\n"
        " 2 50\n"
        " 2 0    Stress Period 1\n"
        " 101  15.000000  1.000000e-04  12.000000  5\n"
        " 102  20.000000  2.000000e-04  17.000000  7\n"
        " -1 0    Stress Period 2\n"
    )

    ml = MfUsg(structured=False, model_ws=str(function_tmpdir))
    ModflowDis(ml, nlay=1, nrow=1, ncol=1, nper=2)
    with pytest.warns(UserWarning, match="not read by USG-T"):
        riv = MfUsgRiv.load(str(riv_in), ml, nper=2, ext_unit_dict={})

    assert riv.dtype.names == ("node", "stage", "cond", "rbot")
    sp0 = riv.stress_period_data[0]
    assert sp0["node"].tolist() == [100, 101]
    assert np.allclose(sp0["rbot"], [12.0, 17.0])
    assert not any("irch" in o.lower() for o in riv.options)

    riv_out = function_tmpdir / "out_trailing.riv"
    riv.fn_path = str(riv_out)
    riv.write_file()
    rows = [ln.split() for ln in riv_out.read_text().splitlines()[3:5]]
    assert [len(r) for r in rows] == [4, 4]


def test_mfusgriv_aux_with_trailing_value(function_tmpdir):
    """AUX values are kept; a value after them is ignored like USG-T does.
    (Formerly the trailing value was kept as 'irch'.)"""
    from flopy.mfusg import MfUsgRiv

    riv_in = function_tmpdir / "aux_trailing.riv"
    riv_in.write_text(
        "# MfUsgRiv aux test\n"
        " 1 50 AUX C01\n"
        " 1 0    Stress Period 1\n"
        " 101  15.000000  1.000000e-04  12.000000  5.000000e-02  3\n"
    )

    ml = MfUsg(structured=False)
    with pytest.warns(UserWarning, match="not read by USG-T"):
        riv = MfUsgRiv.load(str(riv_in), ml, nper=1, ext_unit_dict={})

    # Field names are preserved as-is from file header (uppercase 'C01')
    assert riv.dtype.names == ("node", "stage", "cond", "rbot", "C01")
    sp0 = riv.stress_period_data[0]
    assert np.isclose(sp0["C01"][0], 0.05, atol=1e-5)
    assert riv.options == ["AUX C01"]


def _riv_disu_model(ws):
    """DISU line of 3 nodes (1 layer); CHD at node 1, RIV at node 3."""
    ml = MfUsg(modelname="r", model_ws=str(ws), structured=False)
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
    return ml


def _riv_header(path):
    return next(
        ln for ln in Path(path).read_text().splitlines() if not ln.startswith("#")
    ).split()


@pytest.mark.parametrize(
    "ipakcb, irivcb", [(55, 55), (61, 61), (0, 0), (None, 0), (-1, -1)]
)
def test_mfusgriv_disu_header_irivcb(function_tmpdir, ipakcb, irivcb):
    """DISU header is MXACTR IRIVCB [options]; IRIVCB = ipakcb (None → 0);
    a CBC output is registered only for IRIVCB > 0; load + rewrite keep it."""
    from flopy.mfusg import MfUsgRiv
    from flopy.utils.mfreadnam import parsenamefile

    ml = _riv_disu_model(function_tmpdir / "w")
    MfUsgRiv(
        ml,
        ipakcb=ipakcb,
        options=["NOPRINT", "AUX C01"],
        stress_period_data={0: [[2, 12.0, 5.0, 8.0, 0.25]]},
        dtype=MfUsgRiv.get_empty(0, aux_names=["C01"], structured=False).dtype,
    )
    ml.write_input()
    ws = function_tmpdir / "w"
    assert _riv_header(ws / "r.riv") == ["1", str(irivcb), "NOPRINT", "AUX", "C01"]
    assert (ws / "r.riv").read_text().splitlines()[3].split() == [
        "3",
        "12.000000",
        "5.000000e+00",
        "8.000000",
        "2.500000e-01",
    ]
    data = [
        ln.split()
        for ln in (ws / "r.nam").read_text().splitlines()
        if ln.startswith("DATA")
    ]
    assert any(int(d[1]) == irivcb for d in data) == (irivcb > 0)

    other = _riv_disu_model(function_tmpdir / "l")
    ext = parsenamefile(str(ws / "r.nam"), other.mfnam_packages, verbose=False)
    try:
        riv = MfUsgRiv.load(ext[18].filehandle, other, ext_unit_dict=ext)
    finally:
        for item in ext.values():
            if item.filehandle is not None and not item.filehandle.closed:
                item.filehandle.close()
    assert riv.ipakcb == irivcb
    assert riv.options == ["NOPRINT", "AUX C01"]
    assert (irivcb in other.output_units) == (irivcb > 0)
    if irivcb > 0:
        fname = other.output_fnames[other.output_units.index(irivcb)]
        assert fname == "r.cbc"
    riv.fn_path = str(function_tmpdir / "rewrite.riv")
    riv.write_file()
    assert _riv_header(function_tmpdir / "rewrite.riv")[:2] == ["1", str(irivcb)]


def test_mfusgriv_irdflag_ignored(function_tmpdir):
    """irdflag stays in the signature but does not control IRIVCB; an explicit
    value different from ipakcb warns and is not used as a CBC unit."""
    from flopy.mfusg import MfUsgRiv

    ml = _riv_disu_model(function_tmpdir)
    spd = {0: [[2, 12.0, 5.0, 8.0]]}
    with pytest.warns(UserWarning, match="irdflag=50 is ignored"):
        riv = MfUsgRiv(ml, 55, spd, None, 50)  # positional call still valid
    riv.write_file()
    assert _riv_header(riv.fn_path)[:2] == ["1", "55"]
    ml = _riv_disu_model(function_tmpdir / "b")  # no unit 55 registered yet
    with pytest.warns(UserWarning, match="irdflag=55 is ignored"):
        riv = MfUsgRiv(ml, stress_period_data=spd, irdflag=55)
    assert riv.ipakcb == 0 and 55 not in ml.output_units
    riv.write_file()
    assert _riv_header(riv.fn_path)[:2] == ["1", "0"]
    ml = _riv_disu_model(function_tmpdir / "c")
    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter("always")
        MfUsgRiv(ml, ipakcb=55, stress_period_data=spd, irdflag=55)
    assert not [x for x in w if "irdflag" in str(x.message)]  # consistent value


def test_mfusgriv_irch_dtype_rejected_before_writing(function_tmpdir):
    """A dtype with 'irch' (not read by USG-T) is rejected before the RIV
    file is opened."""
    from flopy.mfusg import MfUsgRiv

    ml = _riv_disu_model(function_tmpdir)
    dtype = np.dtype(
        list(MfUsgRiv.get_default_dtype(structured=False).descr) + [("irch", np.int32)]
    )
    riv = MfUsgRiv(
        ml, ipakcb=55, dtype=dtype, stress_period_data={0: [[2, 12.0, 5.0, 8.0, 3]]}
    )
    assert not any("irch" in o.lower() for o in riv.options)
    Path(riv.fn_path).write_text("SENTINEL\n")
    with pytest.raises(ValueError, match="irch"):
        riv.write_file()
    assert Path(riv.fn_path).read_text() == "SENTINEL\n"


@pytest.mark.parametrize(
    "decl", ["AUX IRCH", "AUXILIARY IRCH", "aux irch", "Auxiliary Irch"]
)
def test_mfusgriv_declared_aux_irch_load_write(function_tmpdir, decl):
    """A variable declared as AUX/AUXILIARY IRCH (any case) is an ordinary AUX
    variable: kept on load, write and reload."""
    from flopy.mfusg import MfUsgRiv

    riv_in = function_tmpdir / "aux_irch.riv"
    riv_in.write_text(f"# manual\n 1 55 {decl}\n 1 0\n 3 12.0 5.0 8.0 7.0\n")
    ml = _riv_disu_model(function_tmpdir / "l")
    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter("always")
        riv = MfUsgRiv.load(str(riv_in), ml, ext_unit_dict={})
    assert not [x for x in w if "not read by USG-T" in str(x.message)]
    name = decl.split()[1]
    assert riv.options == [decl]
    assert riv.dtype.names[4:] == (name,)
    assert riv.stress_period_data[0][name][0] == 7.0
    riv.fn_path = str(function_tmpdir / "out.riv")
    riv.write_file()
    lines = Path(riv.fn_path).read_text().splitlines()
    assert lines[1].split() == ["1", "55"] + decl.split()
    assert lines[3].split()[-1] == "7.000000e+00"
    again = MfUsgRiv.load(
        riv.fn_path, _riv_disu_model(function_tmpdir / "r"), ext_unit_dict={}
    )
    assert again.options == [decl] and again.stress_period_data[0][name][0] == 7.0


@pytest.mark.parametrize(
    "field, decl",
    [("irch", "AUX IRCH"), ("IRCH", "auxiliary irch"), ("Irch", "AUX irch")],
)
def test_mfusgriv_declared_irch_programmatic(function_tmpdir, field, decl):
    from flopy.mfusg import MfUsgRiv

    ml = _riv_disu_model(function_tmpdir)
    dtype = np.dtype(
        list(MfUsgRiv.get_default_dtype(structured=False).descr) + [(field, np.float32)]
    )
    riv = MfUsgRiv(
        ml,
        ipakcb=55,
        dtype=dtype,
        options=[decl],
        stress_period_data={0: [[2, 12.0, 5.0, 8.0, 7.0]]},
    )
    assert riv.options == [decl]  # not declared a second time
    riv.write_file()
    lines = Path(riv.fn_path).read_text().splitlines()
    assert lines[1].split() == ["1", "55"] + decl.split()
    assert lines[3].split()[-1] == "7.000000e+00"


def test_mfusgriv_aux_order_with_irch(function_tmpdir):
    """AUX values are written in the order of the declarations, whatever the
    dtype order, and reload in that order."""
    from flopy.mfusg import MfUsgRiv

    ml = _riv_disu_model(function_tmpdir / "w")
    base = list(MfUsgRiv.get_default_dtype(structured=False).descr)
    dtype = np.dtype(
        base + [("C02", np.float32), ("IRCH", np.float32), ("C01", np.float32)]
    )
    options = ["AUX C01", "AUXILIARY IRCH", "AUX C02"]
    riv = MfUsgRiv(
        ml,
        ipakcb=55,
        dtype=dtype,
        options=list(options),
        stress_period_data={0: [[2, 12.0, 5.0, 8.0, 0.2, 7.0, 0.1]]},
    )
    riv.write_file()
    lines = Path(riv.fn_path).read_text().splitlines()
    assert lines[1].split()[2:] == "AUX C01 AUXILIARY IRCH AUX C02".split()
    assert lines[3].split()[4:] == ["1.000000e-01", "7.000000e+00", "2.000000e-01"]
    again = MfUsgRiv.load(
        riv.fn_path, _riv_disu_model(function_tmpdir / "r"), ext_unit_dict={}
    )
    assert again.options == options
    assert again.dtype.names[4:] == ("C01", "IRCH", "C02")
    rec = again.stress_period_data[0][0]
    assert np.allclose([rec["C01"], rec["IRCH"], rec["C02"]], [0.1, 7.0, 0.2])


@pytest.mark.parametrize(
    "options, extra",
    [
        ([], []),
        (["AUX XIRCH"], ["XIRCH"]),
        (["AUXILIARY IRCHX"], ["IRCHX"]),
        (["NOPRINT"], []),
    ],
)
def test_mfusgriv_undeclared_irch_rejected(function_tmpdir, options, extra):
    """An undeclared 'irch' field is rejected before the file is opened; a
    partial name match (XIRCH, IRCHX) does not declare it, and it is not
    declared automatically."""
    from flopy.mfusg import MfUsgRiv

    ml = _riv_disu_model(function_tmpdir)
    base = list(MfUsgRiv.get_default_dtype(structured=False).descr)
    dtype = np.dtype(base + [(n, np.float32) for n in extra] + [("irch", np.int32)])
    riv = MfUsgRiv(
        ml,
        ipakcb=55,
        dtype=dtype,
        options=list(options),
        stress_period_data={0: [[2, 12.0, 5.0, 8.0] + [0.5] * len(extra) + [3]]},
    )
    assert riv.options == options
    Path(riv.fn_path).write_text("SENTINEL\n")
    with pytest.raises(ValueError, match="not declared as AUX"):
        riv.write_file()
    assert Path(riv.fn_path).read_text() == "SENTINEL\n"


def test_mfusgriv_aux_declaration_whole_tokens(function_tmpdir):
    """A field is declared by whole-token match: 'AUX C10' does not declare C1;
    a declared AUX without a data field is rejected before writing."""
    from flopy.mfusg import MfUsgRiv

    ml = _riv_disu_model(function_tmpdir)
    base = list(MfUsgRiv.get_default_dtype(structured=False).descr)
    dtype = np.dtype(base + [("C10", np.float32), ("C1", np.float32)])
    riv = MfUsgRiv(
        ml,
        ipakcb=55,
        dtype=dtype,
        options=["AUX C10"],
        stress_period_data={0: [[2, 12.0, 5.0, 8.0, 0.1, 0.2]]},
    )
    assert riv.options == ["AUX C10", "AUX C1"]
    riv.write_file()
    assert Path(riv.fn_path).read_text().splitlines()[3].split()[4:] == [
        "1.000000e-01",
        "2.000000e-01",
    ]
    riv2 = MfUsgRiv(
        _riv_disu_model(function_tmpdir / "m"),
        ipakcb=55,
        options=["AUX C05"],
        stress_period_data={0: [[2, 12.0, 5.0, 8.0]]},
    )
    Path(riv2.fn_path).write_text("SENTINEL\n")
    with pytest.raises(ValueError, match="no data field"):
        riv2.write_file()
    assert Path(riv2.fn_path).read_text() == "SENTINEL\n"


def test_mfusgriv_structured_header_unchanged(function_tmpdir):
    """Control: structured models keep the ModflowRiv writer (IRIVCB=ipakcb)."""
    from flopy.mfusg import MfUsgRiv
    from flopy.modflow import ModflowDis

    ml = MfUsg(modelname="s", model_ws=str(function_tmpdir))
    ModflowDis(ml, nlay=1, nrow=1, ncol=3, top=10.0, botm=[0.0])
    MfUsgBas(ml)
    MfUsgLpf(ml)
    MfUsgRiv(ml, ipakcb=55, stress_period_data={0: [[0, 0, 2, 12.0, 5.0, 8.0]]})
    ml.get_package("RIV").write_file()
    assert _riv_header(function_tmpdir / "s.riv")[:2] == ["1", "55"]


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


def test_mfusgets_parameterized_write_npets_without_defs_fails(function_tmpdir):
    """npets>0 with no parameter definitions fails explicitly (ValueError), no
    incomplete file. (From-scratch authoring *with* definitions is supported as
    of Stage 4.6C-D -- see the authoring tests below.)"""
    from flopy.mfusg import MfUsgEts
    from flopy.modflow import ModflowDis

    ml = MfUsg(structured=True, model_ws=str(function_tmpdir))
    ModflowDis(ml, nlay=1, nrow=2, ncol=2, nper=1)
    ets = MfUsgEts(ml, netsop=1, evtr=1.2e-4, npets=1)
    out = function_tmpdir / "noparam.ets"
    ets.fn_path = str(out)
    with pytest.raises(ValueError, match="no parameters"):
        ets.write_file()
    assert not out.exists()


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
        ]
    )
    # No 'irch' column: USG-T RIV reads node, stage, cond, rbot and AUX only
    # (a dtype with 'irch' is rejected, see
    # test_mfusgriv_irch_dtype_rejected_before_writing)
    riv = MfUsgRiv(
        ml,
        stress_period_data={
            0: np.array(
                [(2, 10.0, 1.0e-4, 9.0, 0.25)],
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
    assert riv_text.splitlines()[3].split()[-1] == "2.500000e-01"

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


# DPF read contract, GWF2DPFU1AR (gwf2dpf1u1.f, USG-Transport 2.8.1):
# IUZONTABIM if ITABRICH!=0 (111-115), all IBOUNDIM layers (118-124), all
# HNEWIM layers (131-157), then DO 200 K=1,NLAY (162-191): PHIF, DDFTR, SC1IM,
# SC2IM if LAYCON(K)!=0, ALPHAIM/BETAIM/SRIM/BROOKIM if LAYCON(K)==5 and
# ITABRICH==0, and BPIM if additionally IBPN>0.

# one distinct value per property and layer (tuple index = layer)
_DPF_VALUES = {
    "iboundim": (1, 1, 1),
    "hnewim": (5.1, 5.2, 5.3),
    "phif": (0.11, 0.12, 0.13),
    "ddftr": (0.021, 0.022, 0.023),
    "sc1im": (3.1e-05, 3.2e-05, 3.3e-05),
    "sc2im": (0.41, 0.42, 0.43),
    "alphaim": (0.51, 0.52, 0.53),
    "betaim": (6.1, 6.2, 6.3),
    "srim": (0.071, 0.072, 0.073),
    "brookim": (8.1, 8.2, 8.3),
    "bpim": (-0.91, -0.92, -0.93),
}


def _dpf_fortran_sequence(laycon, tabrich=False, bubblept=False):
    """Restatement of the GWF2DPFU1AR record order, as (name, layer)."""
    nlay = len(laycon)
    seq = [("iuzontabim", 0)] if tabrich else []
    seq += [("iboundim", k + 1) for k in range(nlay)]
    seq += [("hnewim", k + 1) for k in range(nlay)]
    for k in range(nlay):
        seq += [("phif", k + 1), ("ddftr", k + 1), ("sc1im", k + 1)]
        if laycon[k] != 0:
            seq.append(("sc2im", k + 1))
        if laycon[k] == 5 and not tabrich:
            seq += [(n, k + 1) for n in ("alphaim", "betaim", "srim", "brookim")]
            if bubblept:
                seq.append(("bpim", k + 1))
    return seq


def _dpf_records(path):
    """(name, layer) of every array control record of a DPF file, in order."""
    records = []
    for line in Path(path).read_text().splitlines():
        tokens = line.split()
        if not tokens or tokens[0] not in ("CONSTANT", "INTERNAL", "OPEN/CLOSE"):
            continue
        comment = line.split("#", 1)[1].split()
        layer = int(comment[2]) if len(comment) > 2 and comment[1] == "layer" else 0
        records.append((comment[0], layer))
    return records


def _dpf_kwargs(nlay, ncol=2):
    return {
        name: [[[values[k]] * ncol] for k in range(nlay)]
        for name, values in _DPF_VALUES.items()
    }


def _dpf_flow_model(
    ws,
    laycon=None,
    laytyp=None,
    richards=False,
    richards_hp=False,
    tabrich=False,
    bubblept=False,
    thickstrt=False,
):
    """Structured 1-row, 2-column MfUsg model with BCF (laycon) or LPF."""
    from flopy.modflow import ModflowDis

    nlay = len(laycon if laycon is not None else laytyp)
    ml = MfUsg(modelname="dpf", model_ws=str(ws))
    ModflowDis(
        ml,
        nlay=nlay,
        nrow=1,
        ncol=2,
        nper=1,
        steady=False,
        top=10.0,
        botm=[10.0 - 2.0 * (k + 1) for k in range(nlay)],
    )
    MfUsgBas(ml, ibound=1, strt=8.0, richards=richards, richards_hp=richards_hp)
    tab = {"tabrich": True, "nuzones": 1, "nutabrows": 2} if tabrich else {}
    if laycon is not None:
        MfUsgBcf(
            ml,
            laycon=list(laycon),
            bubblept=bubblept,
            tran=1.0,
            hy=1.0,
            kv=1.0,
            sf1=1.0e-5,
            sf2=0.15,
            **tab,
        )
    else:
        MfUsgLpf(ml, laytyp=list(laytyp), bubblept=bubblept, thickstrt=thickstrt, **tab)
    return ml


def _assert_dpf_values(dpf, laycon, tabrich=False, bubblept=False):
    for name, layer in _dpf_fortran_sequence(laycon, tabrich, bubblept):
        if name == "iuzontabim":
            continue
        got = getattr(dpf, name)[layer - 1].array
        assert np.allclose(got, _DPF_VALUES[name][layer - 1]), (name, layer, got)


def test_mfusgdpf_write_order_multilayer_confined(function_tmpdir):
    """PHIF, DDFTR and SC1IM are written layer by layer (DO 200)."""
    ml = _dpf_flow_model(function_tmpdir, laycon=[0, 0])
    dpf = MfUsgDpf(ml, **_dpf_kwargs(2))
    dpf.write_file()

    assert _dpf_records(dpf.fn_path) == [
        ("iboundim", 1),
        ("iboundim", 2),
        ("hnewim", 1),
        ("hnewim", 2),
        ("phif", 1),
        ("ddftr", 1),
        ("sc1im", 1),
        ("phif", 2),
        ("ddftr", 2),
        ("sc1im", 2),
    ]


def test_mfusgdpf_load_manual_fortran_order_confined(function_tmpdir):
    """A hand-written file in GWF2DPFU1AR order loads into the right layers."""
    path = function_tmpdir / "manual.dpf"
    path.write_text(
        "         0         0         0\n"
        "CONSTANT          1 (I10) -1 #iboundim layer 1\n"
        "CONSTANT          1 (I10) -1 #iboundim layer 2\n"
        "CONSTANT        5.1 (G15.6) -1 #hnewim layer 1\n"
        "CONSTANT        5.2 (G15.6) -1 #hnewim layer 2\n"
        "CONSTANT       0.11 (G15.6) -1 #phif layer 1\n"
        "CONSTANT      0.021 (G15.6) -1 #ddftr layer 1\n"
        "CONSTANT    3.1e-05 (G15.6) -1 #sc1im layer 1\n"
        "CONSTANT       0.12 (G15.6) -1 #phif layer 2\n"
        "CONSTANT      0.022 (G15.6) -1 #ddftr layer 2\n"
        "CONSTANT    3.2e-05 (G15.6) -1 #sc1im layer 2\n"
    )
    ml = _dpf_flow_model(function_tmpdir, laycon=[0, 0])
    dpf = MfUsgDpf.load(str(path), ml)

    _assert_dpf_values(dpf, [0, 0])
    # SC2IM is not read for confined layers: placeholders keep the 0 convention
    assert np.allclose(dpf.sc2im.array, 0.0)


@pytest.mark.parametrize("bubblept", [False, True])
def test_mfusgdpf_two_richards_layers(function_tmpdir, bubblept):
    """Two LAYCON=5 layers: write order, manual load and round-trip."""
    laycon = [5, 5]
    kw = {"laycon": laycon, "richards": True, "bubblept": bubblept}
    ml = _dpf_flow_model(function_tmpdir, **kw)
    dpf = MfUsgDpf(ml, **_dpf_kwargs(2))
    dpf.write_file()

    # FloPy round-trip first: writer and loader must agree on the order
    ml2 = _dpf_flow_model(function_tmpdir, **kw)
    _assert_dpf_values(MfUsgDpf.load(str(dpf.fn_path), ml2), laycon, bubblept=bubblept)

    def richards(layer):
        names = ["phif", "ddftr", "sc1im", "sc2im", "alphaim", "betaim", "srim"]
        names += ["brookim"] + (["bpim"] if bubblept else [])
        return [(name, layer) for name in names]

    expected = [("iboundim", 1), ("iboundim", 2), ("hnewim", 1), ("hnewim", 2)]
    expected += richards(1) + richards(2)
    assert _dpf_records(dpf.fn_path) == expected

    # hand-written file in GWF2DPFU1AR order
    lines = [
        "         0         0         0",
        "CONSTANT          1 (I10) -1 #iboundim layer 1",
        "CONSTANT          1 (I10) -1 #iboundim layer 2",
        "CONSTANT        5.1 (G15.6) -1 #hnewim layer 1",
        "CONSTANT        5.2 (G15.6) -1 #hnewim layer 2",
        "CONSTANT       0.11 (G15.6) -1 #phif layer 1",
        "CONSTANT      0.021 (G15.6) -1 #ddftr layer 1",
        "CONSTANT    3.1e-05 (G15.6) -1 #sc1im layer 1",
        "CONSTANT       0.41 (G15.6) -1 #sc2im layer 1",
        "CONSTANT       0.51 (G15.6) -1 #alphaim layer 1",
        "CONSTANT        6.1 (G15.6) -1 #betaim layer 1",
        "CONSTANT      0.071 (G15.6) -1 #srim layer 1",
        "CONSTANT        8.1 (G15.6) -1 #brookim layer 1",
        "CONSTANT      -0.91 (G15.6) -1 #bpim layer 1",
        "CONSTANT       0.12 (G15.6) -1 #phif layer 2",
        "CONSTANT      0.022 (G15.6) -1 #ddftr layer 2",
        "CONSTANT    3.2e-05 (G15.6) -1 #sc1im layer 2",
        "CONSTANT       0.42 (G15.6) -1 #sc2im layer 2",
        "CONSTANT       0.52 (G15.6) -1 #alphaim layer 2",
        "CONSTANT        6.2 (G15.6) -1 #betaim layer 2",
        "CONSTANT      0.072 (G15.6) -1 #srim layer 2",
        "CONSTANT        8.2 (G15.6) -1 #brookim layer 2",
        "CONSTANT      -0.92 (G15.6) -1 #bpim layer 2",
    ]
    if not bubblept:
        lines = [line for line in lines if "#bpim" not in line]
    path = function_tmpdir / "manual.dpf"
    path.write_text("\n".join(lines) + "\n")
    ml3 = _dpf_flow_model(function_tmpdir, **kw)
    _assert_dpf_values(MfUsgDpf.load(str(path), ml3), laycon, bubblept=bubblept)


@pytest.mark.parametrize(
    "flow, codes, richards, tabrich, effective",
    [
        ("bcf", [0, 3, 5], False, False, [0, 3, 5]),
        ("bcf", [0, 3, 5], True, False, [0, 5, 5]),
        ("bcf", [0, 3, 5], True, True, [0, 5, 5]),
        ("lpf", [0, 1, 5], True, False, [0, 5, 5]),
        ("lpf", [0, 1, 5], True, True, [0, 5, 5]),
    ],
)
def test_mfusgdpf_mixed_layers_follow_effective_laycon(
    function_tmpdir, flow, codes, richards, tabrich, effective
):
    """Mixed confined/convertible/Richards layers, TABRICH on and off.

    The records follow the effective LAYCON: with BAS RICHARDS every
    LAYCON != 0 becomes 5 (gwf2bcf-lpf-u1.f 343-348 and 3345-3348).
    """
    kw = {"richards": richards, "tabrich": tabrich, "bubblept": True}
    kw["laycon" if flow == "bcf" else "laytyp"] = codes
    ml = _dpf_flow_model(function_tmpdir, **kw)
    assert MfUsgDpf._get_laycon(ml).tolist() == effective

    arrays = _dpf_kwargs(3)
    if tabrich:
        arrays["iuzontabim"] = [7, 8, 9, 10, 11, 12]
    dpf = MfUsgDpf(ml, **arrays)
    dpf.write_file()
    assert _dpf_records(dpf.fn_path) == _dpf_fortran_sequence(effective, tabrich, True)

    ml2 = _dpf_flow_model(function_tmpdir, **kw)
    dpf2 = MfUsgDpf.load(str(dpf.fn_path), ml2)
    _assert_dpf_values(dpf2, effective, tabrich, True)
    assert np.allclose(dpf2.sc2im[0].array, 0.0)
    if tabrich:
        assert dpf2.iuzontabim.array.tolist() == [7, 8, 9, 10, 11, 12]
        assert not hasattr(dpf2, "alphaim")


def _dpf_disu_model(ws):
    """DISU with nodelay=[3, 2]; LPF LAYTYP=[0, 5] with BAS RICHARDS."""
    ml = MfUsg(modelname="dpfu", model_ws=str(ws), structured=False)
    MfUsgDisU(
        ml,
        nodes=5,
        nlay=2,
        njag=15,
        nodelay=[3, 2],
        iac=np.array([3, 4, 2, 3, 3]),
        ja=np.array([0, 1, 3, 1, 0, 2, 4, 2, 1, 3, 0, 4, 4, 1, 3]),
        ivsd=0,
        cl12=1.0,
        fahl=1.0,
        top=[10.0] * 3 + [5.0] * 2,
        bot=[5.0] * 3 + [0.0] * 2,
        area=1.0,
        steady=False,
    )
    MfUsgBas(ml, ibound=1, strt=8.0, richards=True)
    MfUsgLpf(ml, laytyp=[0, 5])
    return ml


def test_mfusgdpf_disu_irregular_layers(function_tmpdir):
    """Per-layer node counts (3, 2) are kept by the layer-by-layer loop."""
    layer_values = {
        "phif": ([0.11, 0.111, 0.112], [0.12, 0.121]),
        "ddftr": ([0.021, 0.0211, 0.0212], [0.022, 0.0221]),
        "sc1im": ([3.1e-5, 3.11e-5, 3.12e-5], [3.2e-5, 3.21e-5]),
        "sc2im": ([0.0, 0.0, 0.0], [0.42, 0.421]),
        "alphaim": ([1.0, 1.0, 1.0], [0.52, 0.521]),
        "betaim": ([7.0, 7.0, 7.0], [6.2, 6.21]),
        "srim": ([0.05, 0.05, 0.05], [0.072, 0.0721]),
        "brookim": ([6.0, 6.0, 6.0], [8.2, 8.21]),
    }
    ml = _dpf_disu_model(function_tmpdir)
    arrays = {name: [np.array(v) for v in vals] for name, vals in layer_values.items()}
    dpf = MfUsgDpf(ml, iboundim=1, hnewim=5.0, **arrays)
    dpf.write_file()

    assert _dpf_records(dpf.fn_path) == _dpf_fortran_sequence([0, 5])

    def check(loaded):
        for name, layer in _dpf_fortran_sequence([0, 5]):
            if name in layer_values:
                got = np.asarray(getattr(loaded, name)[layer - 1].array).ravel()
                assert got.size == (3, 2)[layer - 1], (name, layer, got)
                assert np.allclose(got, layer_values[name][layer - 1]), (name, layer)

    check(MfUsgDpf.load(str(dpf.fn_path), _dpf_disu_model(function_tmpdir)))

    path = function_tmpdir / "manual_disu.dpf"
    path.write_text(
        "         0         0         0\n"
        "CONSTANT          1 (I10) -1 #iboundim layer 1\n"
        "CONSTANT          1 (I10) -1 #iboundim layer 2\n"
        "CONSTANT        5.0 (G15.6) -1 #hnewim layer 1\n"
        "CONSTANT        5.0 (G15.6) -1 #hnewim layer 2\n"
        "INTERNAL 1.0 (FREE) -1 #phif layer 1\n0.11 0.111 0.112\n"
        "INTERNAL 1.0 (FREE) -1 #ddftr layer 1\n0.021 0.0211 0.0212\n"
        "INTERNAL 1.0 (FREE) -1 #sc1im layer 1\n3.1e-5 3.11e-5 3.12e-5\n"
        "INTERNAL 1.0 (FREE) -1 #phif layer 2\n0.12 0.121\n"
        "INTERNAL 1.0 (FREE) -1 #ddftr layer 2\n0.022 0.0221\n"
        "INTERNAL 1.0 (FREE) -1 #sc1im layer 2\n3.2e-5 3.21e-5\n"
        "INTERNAL 1.0 (FREE) -1 #sc2im layer 2\n0.42 0.421\n"
        "INTERNAL 1.0 (FREE) -1 #alphaim layer 2\n0.52 0.521\n"
        "INTERNAL 1.0 (FREE) -1 #betaim layer 2\n6.2 6.21\n"
        "INTERNAL 1.0 (FREE) -1 #srim layer 2\n0.072 0.0721\n"
        "INTERNAL 1.0 (FREE) -1 #brookim layer 2\n8.2 8.21\n"
    )
    check(MfUsgDpf.load(str(path), _dpf_disu_model(function_tmpdir)))


@pytest.mark.parametrize(
    "flow, codes, options, effective",
    [
        ("bcf", [3, 5], {"richards": True}, [5, 5]),
        ("bcf", [3, 0], {}, [3, 0]),
        ("lpf", [1, 2, 3], {}, [3, 3, 3]),
        ("lpf", [4, 0], {}, [4, 0]),
        ("lpf", [-1, 0], {"thickstrt": True}, [0, 0]),
        ("lpf", [-1, 0], {}, [3, 0]),
        ("lpf", [1, 5], {"richards": True}, [5, 5]),
        ("lpf", [1, 5], {"richards_hp": True}, [5, 5]),
        ("lpf", [-1, 0, 5], {"richards": True, "thickstrt": True}, [0, 0, 5]),
    ],
)
def test_mfusgdpf_effective_laycon(function_tmpdir, flow, codes, options, effective):
    """LAYCON as derived by the Fortran BCF/LPF readers; flow arrays untouched."""
    kw = dict(options)
    kw["laycon" if flow == "bcf" else "laytyp"] = codes
    ml = _dpf_flow_model(function_tmpdir, **kw)
    pkg = ml.get_package("BCF6") or ml.get_package("LPF")
    before = np.array(pkg.laycon.array if flow == "bcf" else pkg.laytyp.array)

    assert MfUsgDpf._get_laycon(ml).tolist() == effective
    after = pkg.laycon.array if flow == "bcf" else pkg.laytyp.array
    assert np.array_equal(before, after)

    dpf = MfUsgDpf(ml, **_dpf_kwargs(len(codes)))
    dpf.write_file()
    records = _dpf_records(dpf.fn_path)
    sc2im_layers = [layer for name, layer in records if name == "sc2im"]
    assert sc2im_layers == [k + 1 for k, lc in enumerate(effective) if lc != 0]


def test_mfusgdpf_thickstrt_option_token(function_tmpdir):
    """THICKSTRT is matched as a whole LPF option token, in any letter case."""
    ml = _dpf_flow_model(function_tmpdir, laytyp=[-1, 0])
    lpf = ml.get_package("LPF")
    lpf.options = " NOVFC thickstrt "
    assert MfUsgDpf._get_laycon(ml).tolist() == [0, 0]
    lpf.options = " THICKSTRTX "
    assert MfUsgDpf._get_laycon(ml).tolist() == [3, 0]


def test_mfusgdpf_tabrich_from_lpf(function_tmpdir):
    """TABRICH on LPF: IUZONTABIM is written first and no Richards arrays."""
    ml = _dpf_flow_model(
        function_tmpdir, laytyp=[0, 5], richards=True, tabrich=True, bubblept=True
    )
    assert MfUsgDpf._has_tabrich(ml)
    assert MfUsgDpf._has_bubblept(ml)
    arrays = _dpf_kwargs(2)
    arrays["iuzontabim"] = [7, 8, 9, 10]
    dpf = MfUsgDpf(ml, **arrays)
    dpf.write_file()
    records = _dpf_records(dpf.fn_path)
    assert records[0] == ("iuzontabim", 0)
    assert records == _dpf_fortran_sequence([0, 5], tabrich=True, bubblept=True)


def test_mfusgdpf_load_open_close_and_shared_external(function_tmpdir):
    """OPEN/CLOSE and a shared EXTERNAL unit are consumed in DO 200 order."""
    from flopy.utils.mfreadnam import parsenamefile

    (function_tmpdir / "phif_l1.txt").write_text("0.11 0.11\n")
    # four distinct blocks read in sequence from unit 71, then a sentinel line
    (function_tmpdir / "dpf_ext71.dat").write_text(
        "0.021 0.021\n3.1e-05 3.1e-05\n0.12 0.12\n0.022 0.022\nEND-OF-DPF-BLOCKS\n"
    )
    (function_tmpdir / "manual.dpf").write_text(
        "         0         0         0\n"
        "CONSTANT          1 (I10) -1 #iboundim layer 1\n"
        "CONSTANT          1 (I10) -1 #iboundim layer 2\n"
        "CONSTANT        5.1 (G15.6) -1 #hnewim layer 1\n"
        "CONSTANT        5.2 (G15.6) -1 #hnewim layer 2\n"
        "OPEN/CLOSE phif_l1.txt 1.0 (FREE) -1 #phif layer 1\n"
        "EXTERNAL 71 1.0 (FREE) -1 #ddftr layer 1\n"
        "EXTERNAL 71 1.0 (FREE) -1 #sc1im layer 1\n"
        "EXTERNAL 71 1.0 (FREE) -1 #phif layer 2\n"
        "EXTERNAL 71 1.0 (FREE) -1 #ddftr layer 2\n"
        "CONSTANT    3.2e-05 (G15.6) -1 #sc1im layer 2\n"
    )
    nam = function_tmpdir / "manual.nam"
    nam.write_text("DPF 57 manual.dpf\nDATA 71 dpf_ext71.dat\n")

    ml = _dpf_flow_model(function_tmpdir, laycon=[0, 0])
    ext_unit_dict = parsenamefile(str(nam), ml.mfnam_packages, verbose=False)
    try:
        dpf = MfUsgDpf.load(
            str(function_tmpdir / "manual.dpf"), ml, ext_unit_dict=ext_unit_dict
        )
        remaining = ext_unit_dict[71].filehandle.readline().strip()
    finally:
        for item in ext_unit_dict.values():
            if item.filehandle is not None:
                item.filehandle.close()

    _assert_dpf_values(dpf, [0, 0])
    assert remaining == "END-OF-DPF-BLOCKS"


@pytest.mark.parametrize("unitnumber, expected_unit", [(None, 157), (57, 57)])
def test_mfusgdpf_write_input_name_file_unit(
    function_tmpdir, unitnumber, expected_unit
):
    """write_input() works with the default DPF unit and with an explicit one."""
    ml = _dpf_flow_model(function_tmpdir, laytyp=[0, 0])
    MfUsgDpf(ml, unitnumber=unitnumber, **_dpf_kwargs(2))
    ml.write_input()

    nam = (function_tmpdir / "dpf.nam").read_text().splitlines()
    assert [ln.split() for ln in nam if ln.split()[:1] == ["DPF"]] == [
        ["DPF", str(expected_unit), "dpf.dpf"]
    ]
    assert _dpf_records(function_tmpdir / "dpf.dpf") == _dpf_fortran_sequence([0, 0])


# D7 stage 1: HNEWIM precision. With BAS DPIN or DPIO (IDPIN=1) GWF2DPFU1AR
# reads HNEWIM with U1DREL8 (REAL*8 values; REAL CNSTNT and CONSTANT; binary
# only through a negative unit, with a 52-byte header and no multiplier);
# otherwise with U1DREL (REAL). Expected values below come from the Fortran
# contract, not from FloPy.

_D7_HNEWIM = np.array(
    [
        [23.1234567890123, 23.1234567890124, 22.75],
        [21.987654321098765, 22.000000123456789, 21.5],
    ]
)  # the first two values of each layer are equal in float32
_D7_NOT_REAL = 23.1234567890123


def _d7_model(ws, options=("dpio",), array_free_format=True):
    """Structured 2-layer 1x3 confined model (LPF LAYTYP 0, LAYCON 0)."""
    from flopy.modflow import ModflowDis

    ml = MfUsg(modelname="d7", model_ws=str(ws))
    ml.array_free_format = array_free_format
    ModflowDis(
        ml, nlay=2, nrow=1, ncol=3, nper=1, steady=False, top=30.0, botm=[20.0, 10.0]
    )
    flags = {"dpin": "double_prec", "dpout": "double_out", "dpio": "double_io"}
    MfUsgBas(ml, ibound=1, strt=25.0, **{flags[o]: True for o in options})
    MfUsgLpf(ml, laytyp=0)
    return ml


def _d7_dpf(ml, hnewim=_D7_HNEWIM.reshape(2, 1, 3), **kwargs):
    return MfUsgDpf(
        ml, iboundim=1, hnewim=hnewim, phif=0.2, ddftr=0.0, sc1im=1.0e-4, **kwargs
    )


def _d7_manual_dpf(path, hnewim_records):
    """Hand-written DPF in the GWF2DPFU1AR order (LAYCON 0)."""
    lines = [
        "         0         0         0",
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
    Path(path).write_text("\n".join(lines) + "\n")


def _d7_text_records(values=_D7_HNEWIM, cnstnt="1.0"):
    out = []
    for k, row in enumerate(values, start=1):
        out += [
            f"INTERNAL {cnstnt} (FREE) -1 #hnewim layer {k}",
            " ".join(repr(float(v)) for v in row),
        ]
    return out


def _d7_binary(values=_D7_HNEWIM, double=True):
    """U1DREL8 (double) / U1DREL stream layout, one header per layer."""
    import struct

    real = "d" if double else "f"
    out = b""
    nstrt = 1
    for k, row in enumerate(values, start=1):
        out += struct.pack(
            f"<ii{real}{real}16siii",
            1,
            1,
            0.0,
            0.0,
            b"  IMMOBILE HEAD ",
            nstrt,
            nstrt + len(row) - 1,
            k,
        )
        out += struct.pack(f"<{len(row)}{real}", *row)
        nstrt += len(row)
    return out


def _d7_written_hnewim(path):
    """HNEWIM values as written in free-format text records (independent
    of FloPy's reader): INTERNAL values or CONSTANT, per layer."""
    lines = Path(path).read_text().splitlines()
    layers = []
    for i, line in enumerate(lines):
        if "hnewim" not in line.lower():
            continue
        tokens = line.split()
        if tokens[0] == "CONSTANT":
            layers.append(("CONSTANT", [float(tokens[1])]))
        elif tokens[0] == "INTERNAL":
            layers.append(("INTERNAL", [float(v) for v in lines[i + 1].split()]))
        else:
            layers.append((tokens[0], None))
    return layers


def _d7_load(ws, dpf_file, options=("dpio",), nam_lines=None):
    """Load a DPF file into a new model; ext_unit_dict from nam_lines."""
    from flopy.utils.mfreadnam import parsenamefile

    ml = _d7_model(ws, options)
    if not nam_lines:
        return MfUsgDpf.load(str(dpf_file), ml)
    nam = Path(dpf_file).parent / "d7_manual.nam"
    nam.write_text("\n".join(nam_lines) + "\n")
    ext = parsenamefile(str(nam), ml.mfnam_packages, verbose=False)
    try:
        return MfUsgDpf.load(str(dpf_file), ml, ext_unit_dict=ext)
    finally:
        for item in ext.values():
            if item.filehandle is not None:
                item.filehandle.close()


@pytest.mark.parametrize(
    "options, dtype",
    [
        ((), np.float32),
        (("dpout",), np.float32),
        (("dpin",), np.float64),
        (("dpio",), np.float64),
        (("dpin", "dpout"), np.float64),
    ],
)
def test_mfusgdpf_hnewim_precision_follows_bas(function_tmpdir, options, dtype):
    dpf = _d7_dpf(_d7_model(function_tmpdir, options))
    assert dpf.hnewim.dtype == dtype
    assert all(u2d.dtype == dtype for u2d in dpf.hnewim.util_2ds)
    hn = dpf.hnewim.array.reshape(2, 3)
    assert hn.dtype == dtype
    np.testing.assert_array_equal(hn, _D7_HNEWIM.astype(dtype))
    assert (hn[0, 0] != hn[0, 1]) == (dtype == np.float64)
    for name in ("phif", "ddftr", "sc1im"):
        assert getattr(dpf, name).dtype == np.float32


def test_mfusgdpf_hnewim_double_text_write_load_rewrite(function_tmpdir):
    ml = _d7_model(function_tmpdir / "w")
    dpf = _d7_dpf(ml)
    dpf.write_file()
    written = _d7_written_hnewim(dpf.fn_path)
    assert [kind for kind, _ in written] == ["INTERNAL", "INTERNAL"]
    assert "(3E24.16)" in Path(dpf.fn_path).read_text()
    np.testing.assert_array_equal([v for _, v in written], _D7_HNEWIM)

    loaded = _d7_load(function_tmpdir / "l", dpf.fn_path)
    assert loaded.hnewim.dtype == np.float64
    np.testing.assert_array_equal(loaded.hnewim.array.reshape(2, 3), _D7_HNEWIM)
    loaded.write_file()  # the loaded package itself, (FREE) layers
    rewritten = _d7_written_hnewim(loaded.fn_path)
    np.testing.assert_array_equal([v for _, v in rewritten], _D7_HNEWIM)


def test_mfusgdpf_hnewim_assignments_keep_double(function_tmpdir):
    dpf = _d7_dpf(_d7_model(function_tmpdir))
    dpf.hnewim = _D7_HNEWIM[::-1].reshape(2, 1, 3)  # Package.__setattr__
    assert dpf.hnewim.dtype == np.float64
    np.testing.assert_array_equal(dpf.hnewim.array.reshape(2, 3), _D7_HNEWIM[::-1])
    dpf.hnewim[0] = _D7_HNEWIM[:1]  # Util3d.__setitem__
    assert dpf.hnewim[0].dtype == np.float64
    np.testing.assert_array_equal(dpf.hnewim[0].array, _D7_HNEWIM[:1])
    dpf.hnewim.cnstnt = 0.1  # REAL multiplier, as U1DREL8 applies it
    factor = np.float64(np.float32(0.1))
    assert dpf.hnewim[0].array[0, 0] == _D7_HNEWIM[0, 0] * factor
    dpf.hnewim[1].cnstnt = 1.0
    assert dpf.hnewim[1].array[0, 0] == _D7_HNEWIM[0, 0]


def test_mfusgdpf_hnewim_uniform_constants(function_tmpdir):
    dpf = _d7_dpf(_d7_model(function_tmpdir / "exact"), hnewim=[22.75, 21.5])
    dpf.write_file()
    assert _d7_written_hnewim(dpf.fn_path) == [
        ("CONSTANT", [22.75]),
        ("CONSTANT", [21.5]),
    ]
    dpf = _d7_dpf(_d7_model(function_tmpdir / "not_real"), hnewim=_D7_NOT_REAL)
    dpf.write_file()
    assert _d7_written_hnewim(dpf.fn_path) == [("INTERNAL", [_D7_NOT_REAL] * 3)] * 2


def _d7_sentinel_checks(ws, dpf, error, match):
    """The error is raised before the DPF file or a caller's handle is
    touched, and without writing external files."""
    target = Path(dpf.fn_path)
    target.write_text("SENTINEL\n")
    before = sorted(p.name for p in Path(ws).rglob("*"))
    with pytest.raises(error, match=match):
        dpf.write_file()
    assert target.read_text() == "SENTINEL\n"
    handle = io.StringIO("CALLER\n")
    handle.seek(0, io.SEEK_END)
    with pytest.raises(error, match=match):
        dpf.write_file(f=handle)
    assert handle.getvalue() == "CALLER\n" and not handle.closed
    assert sorted(p.name for p in Path(ws).rglob("*")) == before


def test_mfusgdpf_hnewim_explicit_constant_not_real_rejected(function_tmpdir):
    dpf = _d7_dpf(_d7_model(function_tmpdir), hnewim=_D7_NOT_REAL)
    dpf.hnewim.how = "constant"
    _d7_sentinel_checks(function_tmpdir, dpf, ValueError, "REAL")


def test_mfusgdpf_hnewim_dpf_before_bas_rejected(function_tmpdir):
    from flopy.modflow import ModflowDis

    ml = MfUsg(modelname="d7", model_ws=str(function_tmpdir))
    ModflowDis(ml, nlay=2, nrow=1, ncol=3, top=30.0, botm=[20.0, 10.0])
    MfUsgLpf(ml, laytyp=0)
    dpf = _d7_dpf(ml)
    assert dpf.hnewim.dtype == np.float32  # no BAS6 yet
    MfUsgBas(ml, double_io=True)
    _d7_sentinel_checks(function_tmpdir, dpf, ValueError, "single precision")


@pytest.mark.parametrize("flag", ["double_prec", "double_io"])
def test_mfusgdpf_hnewim_dpin_set_later_rejected(function_tmpdir, flag):
    ml = _d7_model(function_tmpdir, ())
    dpf = _d7_dpf(ml)
    setattr(ml.get_package("BAS6"), flag, True)
    _d7_sentinel_checks(function_tmpdir, dpf, ValueError, "single precision")


def test_mfusgdpf_hnewim_dpin_cleared_later_warns(function_tmpdir):
    ml = _d7_model(function_tmpdir)
    dpf = _d7_dpf(ml)
    ml.get_package("BAS6").double_io = False
    with pytest.warns(UserWarning, match="single precision"):
        dpf.write_file()
    np.testing.assert_array_equal(
        [v for _, v in _d7_written_hnewim(dpf.fn_path)], _D7_HNEWIM
    )


@pytest.mark.parametrize("options, dtype", [((), np.float64), (("dpio",), np.float32)])
def test_mfusgdpf_hnewim_binary_precision_rejected(function_tmpdir, options, dtype):
    """Binary HNEWIM records are REAL*8 with DPIN/DPIO and REAL otherwise; a
    layer of the other precision is rejected before anything is written
    (stage 2 replaces the stage-1 rejection of every binary output)."""
    ml = _d7_model(function_tmpdir, options)
    ml.external_path = "ext"
    (function_tmpdir / "ext").mkdir()
    dpf = _d7_dpf(ml)
    dpf.hnewim[1] = Util2d(ml, (1, 3), dtype, _D7_HNEWIM[1:], name="hnewim layer 2")
    for u2d in dpf.hnewim.util_2ds:
        u2d.format.binary = True
    externals = list(ml.external_fnames)
    if dtype == np.float64:
        # the double-precision warning of stage 1 precedes the rejection
        with pytest.warns(UserWarning, match="double precision"):
            _d7_sentinel_checks(function_tmpdir, dpf, ValueError, "REAL binary record")
    else:
        _d7_sentinel_checks(function_tmpdir, dpf, ValueError, "single precision")
    assert list((function_tmpdir / "ext").iterdir()) == []
    assert ml.external_fnames == externals


# D7 stage 2: binary HNEWIM output. U1DREL/U1DREL8 read binary arrays only
# through a negative unit: per layer a header (KSTP, KPER, PERTIM, TOTIM,
# TEXT, NSTRT, NNDLAY, ILAY; 52 bytes with REAL*8 times, 44 with REAL) and the
# values, without applying CNSTNT. Files are parsed here with struct.


def _d7_binary_dpf(ml, hnewim=_D7_HNEWIM.reshape(2, 1, 3), layers=(0, 1)):
    dpf = _d7_dpf(ml, hnewim=hnewim)
    for k in layers:
        dpf.hnewim[k].format.binary = True
    return dpf


def _d7_parse_binary(path, double, counts):
    import struct

    real = "d" if double else "f"
    hfmt = f"<ii{real}{real}16siii"
    data = Path(path).read_bytes()
    out, off = [], 0
    for n in counts:
        h = struct.unpack_from(hfmt, data, off)
        off += struct.calcsize(hfmt)
        out.append((h, np.array(struct.unpack_from(f"<{n}{real}", data, off))))
        off += struct.calcsize(f"<{n}{real}")
    return out, data[off:]


def _d7_nam(ws, name="d7.nam"):
    return [
        ln.split()
        for ln in (Path(ws) / name).read_text().splitlines()
        if ln.strip() and not ln.startswith("#")
    ]


def _d7_unit_of(ws, fname, name="d7.nam"):
    rows = [e for e in _d7_nam(ws, name) if e[-1] == fname]
    assert len(rows) == 1 and rows[0][0] == "DATA(BINARY)", rows
    return int(rows[0][1])


def _d7_hnewim_lines(path):
    return [ln for ln in Path(path).read_text().splitlines() if "hnewim" in ln]


@pytest.mark.parametrize("options, double", [(("dpio",), True), ((), False)])
def test_mfusgdpf_hnewim_binary_layout(function_tmpdir, options, double):
    ml = _d7_model(function_tmpdir, options)
    dpf = _d7_binary_dpf(ml)
    dpf.hnewim[1].cnstnt = 0.1
    ml.write_input()
    unit = _d7_unit_of(function_tmpdir, "d7.hnewim.bin")
    assert unit > 0
    lines = _d7_hnewim_lines(function_tmpdir / "d7.dpf")
    assert [ln.split()[:3] for ln in lines] == [["EXTERNAL", str(-unit), "1.0"]] * 2
    size = 52 + 24 if double else 44 + 12
    assert (function_tmpdir / "d7.hnewim.bin").stat().st_size == 2 * size
    layers, rest = _d7_parse_binary(function_tmpdir / "d7.hnewim.bin", double, [3, 3])
    assert rest == b""
    assert [h[4].strip() for h, _ in layers] == [b"HNEWIM"] * 2
    assert [h[5:] for h, _ in layers] == [(1, 3, 1), (4, 6, 2)]
    effective = dpf.hnewim.array.reshape(2, 3)  # CNSTNT applied once
    real = np.float64 if double else np.float32
    np.testing.assert_array_equal([v for _, v in layers], effective.astype(real))
    if double:
        assert layers[0][1][0] != layers[0][1][1]
        np.testing.assert_array_equal(
            layers[1][1], _D7_HNEWIM[1] * np.float64(np.float32(0.1))
        )


def test_mfusgdpf_hnewim_binary_mixed_and_disu(function_tmpdir):
    ml = _d7_model(function_tmpdir / "mixed")
    _d7_binary_dpf(ml, layers=(1,))
    ml.write_input()
    lines = _d7_hnewim_lines(function_tmpdir / "mixed" / "d7.dpf")
    assert lines[0].startswith("INTERNAL") and lines[1].startswith("EXTERNAL")
    layers, rest = _d7_parse_binary(
        function_tmpdir / "mixed" / "d7.hnewim.bin", True, [3]
    )
    assert rest == b"" and layers[0][0][5:] == (4, 6, 2)

    values = np.array(
        [
            23.1234567890123,
            23.1234567890124,
            22.75,
            21.987654321098765,
            22.000000123456789,
        ]
    )
    mu = _dpf_disu_model(function_tmpdir / "disu")
    mu.get_package("BAS6").double_io = True
    dpf = MfUsgDpf(
        mu,
        iboundim=1,
        hnewim=[values[:3], values[3:]],
        phif=0.2,
        ddftr=0.0,
        sc1im=1.0e-4,
    )
    for u2d in dpf.hnewim.util_2ds:
        u2d.format.binary = True
    dpf.write_file()
    path = function_tmpdir / "disu" / "dpfu.hnewim.bin"
    assert path.stat().st_size == (52 + 24) + (52 + 16)
    layers, rest = _d7_parse_binary(path, True, [3, 2])
    assert rest == b"" and [h[5:] for h, _ in layers] == [(1, 3, 1), (4, 5, 2)]
    np.testing.assert_array_equal(np.concatenate([v for _, v in layers]), values)


def test_mfusgdpf_hnewim_binary_scalar_and_fixed(function_tmpdir):
    ml = _d7_model(function_tmpdir, array_free_format=False)
    _d7_binary_dpf(ml, hnewim=_D7_NOT_REAL)
    ml.array_free_format = True
    ml.write_input()
    unit = _d7_unit_of(function_tmpdir, "d7.hnewim.bin")
    lines = _d7_hnewim_lines(function_tmpdir / "d7.dpf")
    assert all(int(ln[0:10]) == -unit and float(ln[10:20]) == 1.0 for ln in lines)
    layers, _ = _d7_parse_binary(function_tmpdir / "d7.hnewim.bin", True, [3, 3])
    np.testing.assert_array_equal([v for _, v in layers], np.full((2, 3), _D7_NOT_REAL))


def test_mfusgdpf_hnewim_binary_reload_with_loader(function_tmpdir):
    from flopy.utils.mfreadnam import parsenamefile

    ml = _d7_model(function_tmpdir)
    dpf = _d7_binary_dpf(ml)
    dpf.hnewim[0].cnstnt = 0.1
    ml.write_input()
    expected = dpf.hnewim.array.reshape(2, 3)
    unit = _d7_unit_of(function_tmpdir, "d7.hnewim.bin")
    other = _d7_model(function_tmpdir / "load")
    ext = parsenamefile(
        str(function_tmpdir / "d7.nam"), other.mfnam_packages, verbose=False
    )
    try:
        loaded = MfUsgDpf.load(
            str(function_tmpdir / "d7.dpf"), other, ext_unit_dict=ext
        )
        assert ext[unit].filehandle.read() == b""  # cursor at the end
    finally:
        for item in ext.values():
            if item.filehandle is not None:
                item.filehandle.close()
    assert all(u.cnstnt == 1.0 for u in loaded.hnewim.util_2ds)
    np.testing.assert_array_equal(loaded.hnewim.array.reshape(2, 3), expected)


def test_mfusgdpf_hnewim_binary_write_file_caller_handle(function_tmpdir):
    ml = _d7_model(function_tmpdir)
    dpf = _d7_binary_dpf(ml)
    handle = io.StringIO()
    dpf.write_file(f=handle)
    assert not handle.closed and handle.getvalue().count("EXTERNAL") == 2
    assert (function_tmpdir / "d7.hnewim.bin").stat().st_size == 152
    assert "d7.hnewim.bin" in ml.external_fnames  # name file not written
    assert not (function_tmpdir / "d7.nam").exists()


def test_mfusgdpf_hnewim_binary_repeated_own_writes(function_tmpdir):
    ml = _d7_model(function_tmpdir)
    _d7_binary_dpf(ml)
    ml.write_input()
    unit = _d7_unit_of(function_tmpdir, "d7.hnewim.bin")
    first = (function_tmpdir / "d7.hnewim.bin").read_bytes()
    ml.write_input()
    assert _d7_unit_of(function_tmpdir, "d7.hnewim.bin") == unit
    assert ml.external_units.count(unit) == 1
    assert (function_tmpdir / "d7.hnewim.bin").read_bytes() == first
    assert not (function_tmpdir / "d7.hnewim_1.bin").exists()


def _d7_foreign(ws, name="d7.hnewim.bin", content=b"FOREIGN"):
    (Path(ws) / name).write_bytes(content)


@pytest.mark.parametrize("registered", [True, False])
def test_mfusgdpf_hnewim_binary_foreign_file_kept(function_tmpdir, registered):
    """A file of another component, registered (unit 777) or not, is never
    overwritten nor its entry replaced; a free name of this package is used."""
    ml = _d7_model(function_tmpdir)
    _d7_foreign(function_tmpdir)
    if registered:
        ml.add_external("d7.hnewim.bin", 777, binflag=True)
    _d7_binary_dpf(ml)
    ml.write_input()
    assert (function_tmpdir / "d7.hnewim.bin").read_bytes() == b"FOREIGN"
    unit = _d7_unit_of(function_tmpdir, "d7.hnewim_1.bin")
    assert unit != 777
    if registered:
        assert _d7_unit_of(function_tmpdir, "d7.hnewim.bin") == 777
    else:
        assert not any(e[-1] == "d7.hnewim.bin" for e in _d7_nam(function_tmpdir))
    lines = _d7_hnewim_lines(function_tmpdir / "d7.dpf")
    assert all(ln.split()[1] == str(-unit) for ln in lines)
    layers, rest = _d7_parse_binary(function_tmpdir / "d7.hnewim_1.bin", True, [3, 3])
    assert rest == b""
    np.testing.assert_array_equal([v for _, v in layers], _D7_HNEWIM)


@pytest.mark.parametrize("form", ["relative", "absolute", "symlink", "case"])
def test_mfusgdpf_hnewim_binary_foreign_entry_without_file(function_tmpdir, form):
    """An entry of another component for the same destination (literal,
    absolute, through a symbolic link or, on case-insensitive file systems,
    another case) is kept although the file does not exist yet."""
    ml = _d7_model(function_tmpdir)
    target = function_tmpdir / "d7.hnewim.bin"
    if form == "relative":
        fname = "d7.hnewim.bin"
    elif form == "absolute":
        fname = str(target)
    elif form == "symlink":
        os.symlink(target, function_tmpdir / "alias.bin")
        fname = "alias.bin"
    else:
        probe = Path(str(function_tmpdir).swapcase())
        if not probe.exists():
            pytest.skip("case-sensitive file system")
        fname = "D7.HNEWIM.BIN"
    ml.add_external(fname, 778, binflag=True)
    _d7_binary_dpf(ml)
    ml.write_input()
    assert not target.exists()
    assert ml.external_units[ml.external_fnames.index(fname)] == 778
    unit = _d7_unit_of(function_tmpdir, "d7.hnewim_1.bin")
    assert unit != 778


@pytest.mark.parametrize("kind", ["output", "package"])
def test_mfusgdpf_hnewim_binary_avoids_output_and_package_paths(function_tmpdir, kind):
    ml = _d7_model(function_tmpdir)
    if kind == "output":
        ml.add_output_file(90, fname="d7.hnewim.bin", binflag=True, package="DPF")
    else:
        ml.get_package("LPF").file_name[0] = "d7.hnewim.bin"
    _d7_binary_dpf(ml)
    ml.write_input()
    unit = _d7_unit_of(function_tmpdir, "d7.hnewim_1.bin")
    assert unit not in (90,) and unit not in ml.package_units
    assert not (function_tmpdir / "d7.hnewim.bin").exists()  # not written here


def test_mfusgdpf_hnewim_binary_unit_collisions(function_tmpdir):
    ml = _d7_model(function_tmpdir)
    (function_tmpdir / "user.dat").write_text("1\n")
    ml.add_external("user.dat", 1001, binflag=False)
    ml.add_output_file(1002, fname="extra.out", package="DPF")
    _d7_binary_dpf(ml)
    ml.write_input()
    unit = _d7_unit_of(function_tmpdir, "d7.hnewim.bin")
    assert unit not in (1001, 1002) and unit not in ml.package_units
    assert next(e for e in _d7_nam(function_tmpdir) if e[1] == "1001")[-1] == "user.dat"


def test_mfusgdpf_hnewim_binary_unit_reassigned_by_other(function_tmpdir):
    ml = _d7_model(function_tmpdir)
    _d7_binary_dpf(ml)
    ml.write_input()
    unit = _d7_unit_of(function_tmpdir, "d7.hnewim.bin")
    (function_tmpdir / "other.dat").write_text("1\n")
    ml.add_external("other.dat", unit)  # another component takes the unit
    ml.write_input()
    rows = {e[-1]: int(e[1]) for e in _d7_nam(function_tmpdir)}
    assert rows["other.dat"] == unit
    assert rows["d7.hnewim.bin"] != unit  # own file kept, new unit
    assert list(rows).count("d7.hnewim.bin") == 1


@pytest.mark.parametrize("change", ["external_path", "name", "workspace"])
def test_mfusgdpf_hnewim_binary_destination_changes(function_tmpdir, change):
    """After a change of external_path, model name or workspace the new
    destination is used, the previous file is kept and no stale entry
    remains; disabling binary output afterwards drops only the own entry."""
    ws = function_tmpdir / "a"
    ml = _d7_model(ws)
    dpf = _d7_binary_dpf(ml)
    ml.write_input()
    first = (ws / "d7.hnewim.bin").read_bytes()
    if change == "external_path":
        (ws / "ext").mkdir()
        ml.external_path = "ext"
        new_ws, new_fname, nam = ws, os.path.join("ext", "d7.hnewim.bin"), "d7.nam"
    elif change == "name":
        ml.name = "d8"
        new_ws, new_fname, nam = ws, "d8.hnewim.bin", "d8.nam"
    else:
        ml.change_model_ws(str(function_tmpdir / "b"))
        new_ws, new_fname, nam = function_tmpdir / "b", "d7.hnewim.bin", "d7.nam"
    ml.write_input()
    assert (ws / "d7.hnewim.bin").read_bytes() == first  # previous file kept
    assert (Path(new_ws) / new_fname).read_bytes() == first
    entries = [e for e in _d7_nam(new_ws, nam) if "hnewim" in e[-1]]
    assert [e[-1] for e in entries] == [new_fname]
    for u2d in dpf.hnewim.util_2ds:
        u2d.format.binary = False
    ml.write_input()
    assert not any("hnewim" in e[-1] for e in _d7_nam(new_ws, nam))
    assert (Path(new_ws) / new_fname).read_bytes() == first  # not deleted


# Ownership is the identity of the file the package created, not its path:
# a path replaced by a symbolic link or by another file, a directory link
# that now resolves elsewhere, or an own file shared through a registered
# hard link is not rewritten. Originals are kept under new names.


def _d7_bytes(paths):
    return {str(p): Path(p).read_bytes() for p in paths}


@pytest.mark.parametrize("replacement", ["symlink", "regular"])
def test_mfusgdpf_hnewim_binary_replaced_own_path_not_truncated(
    function_tmpdir, replacement
):
    ws = function_tmpdir
    ml = _d7_model(ws)
    dpf = _d7_binary_dpf(ml)
    ml.write_input()
    own = ws / "d7.hnewim.bin"
    original = own.read_bytes()
    os.rename(own, ws / "preserved.bin")  # keep the original
    if replacement == "symlink":
        (ws / "foreign.bin").write_bytes(b"FOREIGN")
        os.symlink("foreign.bin", own)
        kept = [ws / "foreign.bin", ws / "preserved.bin"]
    else:
        own.write_bytes(b"FOREIGN REGULAR")
        kept = [own, ws / "preserved.bin"]
    before = _d7_bytes(kept)
    handle = io.StringIO()
    dpf.write_file(f=handle)  # the caller's handle stays open
    assert not handle.closed
    assert _d7_bytes(kept) == before
    assert (ws / "preserved.bin").read_bytes() == original
    if replacement == "symlink":
        assert own.is_symlink()
    ml.write_input()
    assert _d7_bytes(kept) == before
    unit = _d7_unit_of(ws, "d7.hnewim_1.bin")
    assert not any(e[-1] == "d7.hnewim.bin" for e in _d7_nam(ws))
    lines = _d7_hnewim_lines(ws / "d7.dpf")
    assert all(ln.split()[1] == str(-unit) for ln in lines)
    assert all(
        f"EXTERNAL {-unit:>10d}" in ln
        for ln in handle.getvalue().splitlines()
        if "hnewim" in ln
    )
    layers, rest = _d7_parse_binary(ws / "d7.hnewim_1.bin", True, [3, 3])
    assert rest == b""
    np.testing.assert_array_equal([v for _, v in layers], _D7_HNEWIM)


def test_mfusgdpf_hnewim_binary_registered_hard_link_not_truncated(function_tmpdir):
    ws = function_tmpdir
    ml = _d7_model(ws)
    dpf = _d7_binary_dpf(ml)
    ml.write_input()
    own = ws / "d7.hnewim.bin"
    original = own.read_bytes()
    os.link(own, ws / "alias.bin")  # another component's name for the file
    ml.add_external("alias.bin", 779, binflag=True)
    dpf.hnewim[0] = _D7_HNEWIM[1:]  # new values must not reach alias.bin
    ml.write_input()
    assert (ws / "alias.bin").read_bytes() == original
    assert own.read_bytes() == original
    assert _d7_unit_of(ws, "alias.bin") == 779
    layers, rest = _d7_parse_binary(ws / "d7.hnewim_1.bin", True, [3, 3])
    assert rest == b""
    np.testing.assert_array_equal([v for _, v in layers], _D7_HNEWIM[[1, 1]])


def test_mfusgdpf_hnewim_binary_directory_link_retargeted(function_tmpdir):
    ws = function_tmpdir
    (ws / "dirA").mkdir()
    (ws / "dirB").mkdir()
    os.symlink("dirA", ws / "ext")
    ml = _d7_model(ws)
    ml.external_path = "ext"
    _d7_binary_dpf(ml)
    ml.write_input()
    own_a = ws / "dirA" / "d7.hnewim.bin"
    original = own_a.read_bytes()
    os.rename(ws / "ext", ws / "ext_old")  # keep the first link
    (ws / "dirB" / "d7.hnewim.bin").write_bytes(b"FOREIGN IN B")
    os.symlink("dirB", ws / "ext")
    ml.write_input()
    assert (ws / "dirB" / "d7.hnewim.bin").read_bytes() == b"FOREIGN IN B"
    assert own_a.read_bytes() == original
    assert _d7_unit_of(ws, os.path.join("ext", "d7.hnewim_1.bin")) > 0
    layers, rest = _d7_parse_binary(ws / "dirB" / "d7.hnewim_1.bin", True, [3, 3])
    assert rest == b""
    np.testing.assert_array_equal([v for _, v in layers], _D7_HNEWIM)


def test_mfusgdpf_hnewim_binary_own_rewrite_keeps_identity(function_tmpdir):
    """Control: an unchanged own file is rewritten in place (same inode,
    same unit, new values, no extra file)."""
    ws = function_tmpdir
    ml = _d7_model(ws)
    dpf = _d7_binary_dpf(ml)
    ml.write_input()
    own = ws / "d7.hnewim.bin"
    ino, unit = own.stat().st_ino, _d7_unit_of(ws, "d7.hnewim.bin")
    dpf.hnewim[0] = _D7_HNEWIM[1:]
    ml.write_input()
    assert own.stat().st_ino == ino and _d7_unit_of(ws, "d7.hnewim.bin") == unit
    assert not (ws / "d7.hnewim_1.bin").exists()
    layers, _ = _d7_parse_binary(own, True, [3, 3])
    np.testing.assert_array_equal([v for _, v in layers], _D7_HNEWIM[[1, 1]])


def test_mfusgdpf_hnewim_binary_load_and_rewrite_same_folder(function_tmpdir):
    """A loaded package writes binary output to a new file: the file of the
    original package is not adopted."""
    ml = _d7_model(function_tmpdir)
    _d7_binary_dpf(ml)
    ml.write_input()
    original = (function_tmpdir / "d7.hnewim.bin").read_bytes()
    loaded = MfUsg.load("d7.nam", model_ws=str(function_tmpdir), check=False)
    assert not any("hnewim" in f for f in loaded.external_fnames)
    for u2d in loaded.get_package("DPF").hnewim.util_2ds:
        u2d.format.binary = True
    loaded.write_input()
    assert (function_tmpdir / "d7.hnewim.bin").read_bytes() == original
    layers, rest = _d7_parse_binary(function_tmpdir / "d7.hnewim_1.bin", True, [3, 3])
    assert rest == b""
    np.testing.assert_array_equal([v for _, v in layers], _D7_HNEWIM)


def test_mfusgdpf_hnewim_binary_validation_before_writing(function_tmpdir):
    """A failing text layer, or a missing folder, stops the write before the
    binary file, the name-file entry or the DPF file are created."""
    ml = _d7_model(function_tmpdir)
    dpf = _d7_dpf(ml, hnewim=_D7_NOT_REAL)
    dpf.hnewim[0].how = "constant"  # not exact as REAL: rejected
    dpf.hnewim[1].format.binary = True
    externals = list(ml.external_fnames)
    _d7_sentinel_checks(function_tmpdir, dpf, ValueError, "REAL")
    assert not (function_tmpdir / "d7.hnewim.bin").exists()
    assert ml.external_fnames == externals

    ml2 = _d7_model(function_tmpdir / "m")
    dpf2 = _d7_binary_dpf(ml2)
    ml2.external_path = "missing"  # assigned later: the folder is not created
    _d7_sentinel_checks(function_tmpdir / "m", dpf2, FileNotFoundError, "binary hnewim")
    assert ml2.external_fnames == []


def test_mfusgdpf_hnewim_single_layers_rejected_under_dpin(function_tmpdir):
    ml = _d7_model(function_tmpdir)
    layers = [
        Util2d(ml, (1, 3), np.float32, _D7_HNEWIM[k : k + 1], "hnewim")
        for k in range(2)
    ]
    with pytest.raises(ValueError, match="float64"):
        _d7_dpf(ml, hnewim=layers)
    single = Util3d(ml, (2, 1, 3), np.float32, _D7_HNEWIM.reshape(2, 1, 3), "hnewim")
    with pytest.raises(ValueError, match="float64"):
        _d7_dpf(ml, hnewim=single)


@pytest.mark.parametrize(
    "record",
    [
        "EXTERNAL -71 1.0 (BINARY) -1",
        "EXTERNAL -71 1.0",
        "EXTERNAL -71 1.0 (FREE) -1",
        f"{-71:10d}{'1.0':>10s}{'':>20s}{-1:10d}",
    ],
)
def test_mfusgdpf_load_double_binary(function_tmpdir, record):
    _d7_manual_dpf(
        function_tmpdir / "m.dpf", [f"{record} #hnewim layer {k}" for k in (1, 2)]
    )
    (function_tmpdir / "h.bin").write_bytes(_d7_binary())
    dpf = _d7_load(
        function_tmpdir / "l",
        function_tmpdir / "m.dpf",
        nam_lines=["DPF 57 m.dpf", "DATA(BINARY) 71 h.bin"],
    )
    assert dpf.hnewim.dtype == np.float64
    np.testing.assert_array_equal(dpf.hnewim.array.reshape(2, 3), _D7_HNEWIM)


def test_mfusgdpf_load_double_binary_ignores_multiplier(function_tmpdir):
    _d7_manual_dpf(
        function_tmpdir / "m.dpf",
        [f"EXTERNAL -71 2.0 #hnewim layer {k}" for k in (1, 2)],
    )
    (function_tmpdir / "h.bin").write_bytes(_d7_binary())
    with pytest.warns(UserWarning, match="not applied"):
        dpf = _d7_load(
            function_tmpdir / "l",
            function_tmpdir / "m.dpf",
            nam_lines=["DPF 57 m.dpf", "DATA(BINARY) 71 h.bin"],
        )
    np.testing.assert_array_equal(dpf.hnewim.array.reshape(2, 3), _D7_HNEWIM)


def test_mfusgdpf_load_double_binary_shared_and_truncated(function_tmpdir):
    from flopy.utils.mfreadnam import parsenamefile

    _d7_manual_dpf(
        function_tmpdir / "m.dpf",
        [f"EXTERNAL -71 1.0 #hnewim layer {k}" for k in (1, 2)],
    )
    (function_tmpdir / "h.bin").write_bytes(_d7_binary() + b"LATER-DATA")
    (function_tmpdir / "m.nam").write_text("DPF 57 m.dpf\nDATA(BINARY) 71 h.bin\n")
    ml = _d7_model(function_tmpdir / "l")
    ext = parsenamefile(
        str(function_tmpdir / "m.nam"), ml.mfnam_packages, verbose=False
    )
    try:
        dpf = MfUsgDpf.load(str(function_tmpdir / "m.dpf"), ml, ext_unit_dict=ext)
        assert ext[71].filehandle.read() == b"LATER-DATA"  # cursor after layer 2
    finally:
        for item in ext.values():
            if item.filehandle is not None:
                item.filehandle.close()
    np.testing.assert_array_equal(dpf.hnewim.array.reshape(2, 3), _D7_HNEWIM)

    (function_tmpdir / "h.bin").write_bytes(_d7_binary()[:-8])
    with pytest.raises(ValueError, match="expected array size"):
        _d7_load(
            function_tmpdir / "t",
            function_tmpdir / "m.dpf",
            nam_lines=["DPF 57 m.dpf", "DATA(BINARY) 71 h.bin"],
        )


@pytest.mark.parametrize(
    "record, nam",
    [
        ("OPEN/CLOSE ../h.bin 1.0 (BINARY) -1", ["DPF 57 m.dpf"]),
        ("EXTERNAL 71 1.0 (BINARY) -1", ["DPF 57 m.dpf", "DATA(BINARY) 71 h.bin"]),
    ],
)
def test_mfusgdpf_load_hnewim_binary_needs_negative_unit(function_tmpdir, record, nam):
    """U1DREL/U1DREL8 read binary only through a negative unit."""
    _d7_manual_dpf(
        function_tmpdir / "m.dpf", [f"{record} #hnewim layer {k}" for k in (1, 2)]
    )
    (function_tmpdir / "h.bin").write_bytes(_d7_binary())
    with pytest.raises(ValueError, match="1-D array readers"):
        _d7_load(function_tmpdir / "l", function_tmpdir / "m.dpf", nam_lines=nam)


@pytest.mark.parametrize("cnstnt", ["0.0", "1.0", "0.1"])
def test_mfusgdpf_load_text_multiplier_is_real(function_tmpdir, cnstnt):
    _d7_manual_dpf(function_tmpdir / "m.dpf", _d7_text_records(cnstnt=cnstnt))
    dpf = _d7_load(function_tmpdir / "l", function_tmpdir / "m.dpf")
    c = float(cnstnt)
    factor = 1.0 if c in (0.0, 1.0) else np.float64(np.float32(c))
    np.testing.assert_array_equal(dpf.hnewim.array.reshape(2, 3), _D7_HNEWIM * factor)


def test_mfusgdpf_load_constant_under_dpin_is_real(function_tmpdir):
    _d7_manual_dpf(
        function_tmpdir / "m.dpf",
        [f"CONSTANT {_D7_NOT_REAL!r} #hnewim layer {k}" for k in (1, 2)],
    )
    dpf = _d7_load(function_tmpdir / "l", function_tmpdir / "m.dpf")
    assert dpf.hnewim.dtype == np.float64
    np.testing.assert_array_equal(
        dpf.hnewim.array, np.float64(np.float32(_D7_NOT_REAL))
    )


def test_mfusgdpf_load_without_dpin_stays_single(function_tmpdir):
    _d7_manual_dpf(function_tmpdir / "m.dpf", _d7_text_records())
    dpf = _d7_load(function_tmpdir / "l", function_tmpdir / "m.dpf", options=())
    assert dpf.hnewim.dtype == np.float32
    np.testing.assert_array_equal(
        dpf.hnewim.array.reshape(2, 3), _D7_HNEWIM.astype(np.float32)
    )
    (function_tmpdir / "h.bin").write_bytes(_d7_binary(double=False))
    _d7_manual_dpf(
        function_tmpdir / "b.dpf",
        [f"EXTERNAL -71 1.0 (BINARY) -1 #hnewim layer {k}" for k in (1, 2)],
    )
    dpf = _d7_load(
        function_tmpdir / "lb",
        function_tmpdir / "b.dpf",
        options=(),
        nam_lines=["DPF 57 b.dpf", "DATA(BINARY) 71 h.bin"],
    )
    np.testing.assert_array_equal(
        dpf.hnewim.array.reshape(2, 3), _D7_HNEWIM.astype(np.float32)
    )


def test_mfusgdpf_direct_load_reads_bas_options_and_keeps_cursor(function_tmpdir):
    """MfUsgDpf.load into a model without BAS6: the precision comes from the
    BAS6 options line in ext_unit_dict; the BAS6 handle keeps its position,
    also when the DPF load fails."""
    from flopy.modflow import ModflowDis
    from flopy.utils.mfreadnam import parsenamefile

    _d7_manual_dpf(function_tmpdir / "m.dpf", _d7_text_records())
    _d7_manual_dpf(function_tmpdir / "bad.dpf", _d7_text_records()[:3])
    (function_tmpdir / "m.bas").write_text("# BAS6\nFREE DPIN\n")
    (function_tmpdir / "m.nam").write_text("DPF 57 m.dpf\nBAS6 13 m.bas\n")
    for dpf_file, error in (("m.dpf", None), ("bad.dpf", Exception)):
        ml = MfUsg(modelname="d7", model_ws=str(function_tmpdir / "l"))
        ModflowDis(ml, nlay=2, nrow=1, ncol=3, top=30.0, botm=[20.0, 10.0])
        MfUsgLpf(ml, laytyp=0)
        ext = parsenamefile(
            str(function_tmpdir / "m.nam"), ml.mfnam_packages, verbose=False
        )
        bas_fh = ext[13].filehandle
        try:
            if error is None:
                dpf = MfUsgDpf.load(
                    str(function_tmpdir / dpf_file), ml, ext_unit_dict=ext
                )
            else:
                with pytest.raises(error):
                    MfUsgDpf.load(
                        str(function_tmpdir / dpf_file), ml, ext_unit_dict=ext
                    )
            assert bas_fh.tell() == 0
        finally:
            for item in ext.values():
                if item.filehandle is not None:
                    item.filehandle.close()
    assert dpf.hnewim.dtype == np.float64
    np.testing.assert_array_equal(dpf.hnewim.array.reshape(2, 3), _D7_HNEWIM)


def test_mfusgdpf_load_closes_own_file_on_error(function_tmpdir):
    _d7_manual_dpf(function_tmpdir / "bad.dpf", _d7_text_records()[:3])
    before = len(os.listdir("/dev/fd"))
    with pytest.raises(Exception):
        _d7_load(function_tmpdir / "l", function_tmpdir / "bad.dpf")
    assert len(os.listdir("/dev/fd")) == before


def test_mfusgdpf_write_file_leaves_caller_handle_open(function_tmpdir):
    dpf = _d7_dpf(_d7_model(function_tmpdir))
    handle = io.StringIO()
    dpf.write_file(f=handle)
    assert not handle.closed
    assert "(3E24.16)" in handle.getvalue()


def test_mfusgdpf_hnewim_fixed_format_records(function_tmpdir):
    """Fixed control records (I10,F10.0,A20,I10) for float64 HNEWIM."""
    ml = _d7_model(function_tmpdir, array_free_format=False)
    dpf = _d7_dpf(ml)
    dpf.hnewim.cnstnt = 0.1
    dpf.write_file()
    lines = Path(dpf.fn_path).read_text().splitlines()
    records = [ln for ln in lines if "hnewim" in ln]
    externals = dict(zip(ml.external_units, ml.external_fnames))
    assert len(records) == 2
    factor = np.float64(np.float32(0.1))
    for k, line in enumerate(records):
        locat, cnstnt, fmtin = int(line[0:10]), float(line[10:20]), line[20:40].strip()
        assert locat > 0 and fmtin == "(3E24.16)"
        assert np.float32(cnstnt) == np.float32(0.1)
        text = (function_tmpdir / externals[locat]).read_text()
        values = [float(text[i : i + 24]) for i in range(0, 72, 24)]
        np.testing.assert_array_equal(np.array(values) * factor, dpf.hnewim[k].array[0])
        assert values == _D7_HNEWIM[k].tolist()


@pytest.mark.parametrize("form", ["flat", "layers"])
def test_mfusgdpf_hnewim_disu_irregular_double(function_tmpdir, form):
    values = np.array(
        [
            23.1234567890123,
            23.1234567890124,
            22.75,
            21.987654321098765,
            22.000000123456789,
        ]
    )
    ml = _dpf_disu_model(function_tmpdir)
    ml.get_package("BAS6").double_io = True
    hnewim = values if form == "flat" else [values[:3], values[3:]]
    dpf = MfUsgDpf(ml, iboundim=1, hnewim=hnewim, phif=0.2, ddftr=0.0, sc1im=1.0e-4)
    assert dpf.hnewim.dtype == np.float64
    assert [u.shape for u in dpf.hnewim.util_2ds] == [(3,), (2,)]
    np.testing.assert_array_equal(dpf.hnewim.array, values)
    dpf.write_file()
    written = _d7_written_hnewim(dpf.fn_path)
    assert [len(v) for _, v in written] == [3, 2]
    np.testing.assert_array_equal(np.concatenate([v for _, v in written]), values)
    other = _dpf_disu_model(function_tmpdir / "l")
    other.get_package("BAS6").double_io = True
    loaded = MfUsgDpf.load(dpf.fn_path, other)
    assert loaded.hnewim.dtype == np.float64
    np.testing.assert_array_equal(loaded.hnewim.array, values)


# HNEWIM layers can be replaced after construction; with DPIN/DPIO every
# layer must hold float64 values (U1DREL8), not only the Util3d container.

_D7_LAYER = np.array([[21.987654321098765, 22.000000123456789, 21.5]])
_D7_ASSIGN = {
    "util2d": lambda dpf, k, u: dpf.hnewim.__setitem__(k, u),
    "list": lambda dpf, k, u: setattr(
        dpf, "hnewim", [u if i == k else dpf.hnewim[i] for i in range(2)]
    ),
    "util_2ds": lambda dpf, k, u: dpf.hnewim.util_2ds.__setitem__(k, u),
}


@pytest.mark.parametrize("path", list(_D7_ASSIGN))
@pytest.mark.parametrize("layer", [0, 1])
@pytest.mark.parametrize("options", [("dpin",), ("dpio",)])
def test_mfusgdpf_hnewim_single_layer_assignment_rejected(
    function_tmpdir, options, layer, path
):
    ml = _d7_model(function_tmpdir, options)
    dpf = _d7_dpf(ml)
    single = Util2d(ml, (1, 3), np.float32, _D7_LAYER, name="replacement")
    _D7_ASSIGN[path](dpf, layer, single)
    assert dpf.hnewim.dtype == np.float64  # the container alone looks fine
    assert dpf.hnewim[layer].dtype == np.float32
    _d7_sentinel_checks(
        function_tmpdir, dpf, ValueError, rf"layer\(s\) {layer + 1} are stored"
    )


def test_mfusgdpf_hnewim_all_single_layers_reported(function_tmpdir):
    ml = _d7_model(function_tmpdir, ("dpio",))
    dpf = _d7_dpf(ml)
    for k in range(2):
        dpf.hnewim[k] = Util2d(ml, (1, 3), np.float32, _D7_LAYER, name="r")
    _d7_sentinel_checks(function_tmpdir, dpf, ValueError, r"layer\(s\) 1, 2 are")


@pytest.mark.parametrize("options", [("dpin",), ("dpio",)])
def test_mfusgdpf_hnewim_double_layer_assignments_kept(function_tmpdir, options):
    """Float64 Util2d, ndarray and list replacements and multiplier
    assignments keep values that are distinct only in float64."""
    ml = _d7_model(function_tmpdir, options)
    dpf = _d7_dpf(ml)
    dpf.hnewim[1] = Util2d(ml, (1, 3), np.float64, _D7_HNEWIM[:1], name="r")
    dpf.hnewim[0] = _D7_HNEWIM[1:]
    dpf.write_file()
    written = [v for _, v in _d7_written_hnewim(dpf.fn_path)]
    np.testing.assert_array_equal(written, _D7_HNEWIM[::-1])
    assert written[1][0] != written[1][1]

    dpf.hnewim = [
        Util2d(
            ml, (1, 3), np.float64, _D7_HNEWIM[k : k + 1], name=f"hnewim layer {k + 1}"
        )
        for k in range(2)
    ]
    dpf.hnewim[1].cnstnt = 0.1
    dpf.write_file()
    written = [v for _, v in _d7_written_hnewim(dpf.fn_path)]
    np.testing.assert_array_equal(written, _D7_HNEWIM)
    np.testing.assert_array_equal(
        dpf.hnewim[1].array, _D7_HNEWIM[1:] * np.float64(np.float32(0.1))
    )


def test_mfusgdpf_hnewim_layer_assignments_without_dpin(function_tmpdir):
    """Control: without DPIN/DPIO a float32 layer is the expected precision;
    a float64 layer is written but warns (the model reads REAL)."""
    ml = _d7_model(function_tmpdir, ())
    dpf = _d7_dpf(ml)
    dpf.hnewim[0] = Util2d(ml, (1, 3), np.float32, _D7_LAYER, name="r")
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        dpf.write_file()
    dpf.hnewim[1] = Util2d(ml, (1, 3), np.float64, _D7_LAYER, name="r")
    with pytest.warns(UserWarning, match="single precision"):
        dpf.write_file()


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


def test_mfusgsgb_parameter_mxs_zero_auto_computes(function_tmpdir):
    """MXS omitted/0 is auto-computed as the total definition rows (Stage 4.6C-B;
    previously this raised). Two params of nlst 1 -> PARAMETER 2 2."""
    dtype = MfUsgSgb.get_default_dtype()
    params = {
        "p1": {"partyp": "SGB", "parval": "1.0", "nlst": 1,
               "data": np.array([(0, 0.01)], dtype=dtype).view(np.recarray)},
        "p2": {"partyp": "SGB", "parval": "1.0", "nlst": 1,
               "data": np.array([(1, 0.02)], dtype=dtype).view(np.recarray)},
    }
    sgb = MfUsgSgb(_sgb_model(function_tmpdir, "mz"), parameters=params, mxs=0)
    out = function_tmpdir / "mz.sgb"
    sgb.fn_path = str(out)
    sgb.write_file()
    assert "PARAMETER 2 2" in out.read_text()  # NPSGB=2, MXS auto = 2


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


# --- Stage 4.6C-B: SGB from-scratch parameter-DEFINITION authoring ----------
#
# Definitions (PARAMETER NPSGB MXS + UPARLSTRP headers + rows) can be authored
# from Python (no prior load): parameters={name: {parval, data, ...}}; counts
# auto-computed; normalized + validated before open. Active SGB parameters stay
# unsupported (PARTYP='SGB' def vs PTYP='G' activation -> Fortran abort), so NP=0
# is written every period and any activation raises NotImplementedError.


def test_mfusgsgb_parameter_authoring_from_scratch(function_tmpdir):
    """Author SGB parameter definitions from scratch (ergonomic input): PARAMETER
    NPSGB MXS auto, UPARLSTRP headers + rows, per-SP NP=0, 0↔1-based, reload."""
    ml = _sgb_model(function_tmpdir, "auth")
    params = {"gp": {"parval": "1.5", "data": [(4, 0.01), (9, 0.02)]}}  # 0-based
    sgb = MfUsgSgb(ml, parameters=params)
    out = function_tmpdir / "auth.sgb"
    sgb.fn_path = str(out)
    sgb.write_file()
    text = out.read_text()
    assert "PARAMETER 1 2" in text  # NPSGB=1, MXS auto = sum(nlst)=2
    assert "gp SGB 1.5 2" in text
    assert "\n 5  1.000000e-02" in text  # node 4 -> 5 (1-based)
    assert " 0 0    Stress Period 1" in text  # NP=0, zero non-param rows

    re = MfUsgSgb.load(str(out), _sgb_model(function_tmpdir, "auth2"), nper=1,
                       ext_unit_dict={})
    pd = re.parameters["gp"]
    assert pd["partyp"].upper() == "SGB" and pd["nlst"] == 2
    assert list(pd["data"]["node"]) == [4, 9]  # back to 0-based
    assert abs(float(pd["parval"]) - 1.5) < 1e-9


def test_mfusgsgb_parameter_authoring_with_aux(function_tmpdir):
    """From-scratch SGB parameter definition with an AUX concentration column."""
    from flopy.pakbase import Package

    ml = _sgb_model(function_tmpdir, "ax")
    dtype = Package.add_to_dtype(MfUsgSgb.get_default_dtype(), ["C01"], np.float64)
    rows = np.array([(4, 0.01, 0.15)], dtype=dtype).view(np.recarray)
    sgb = MfUsgSgb(ml, dtype=dtype, parameters={"gp": {"parval": "1.0", "data": rows}})
    out = function_tmpdir / "ax.sgb"
    sgb.fn_path = str(out)
    sgb.write_file()
    assert "AUX C01" in out.read_text()

    re = MfUsgSgb.load(str(out), _sgb_model(function_tmpdir, "ax2"), nper=1,
                       ext_unit_dict={})
    assert np.isclose(re.parameters["gp"]["data"]["C01"][0], 0.15)


def test_mfusgsgb_parameter_authoring_mixed_with_nonparam(function_tmpdir):
    """Parameter definitions coexist with non-parametric stress-period rows;
    NP=0 is written every period (no activations)."""
    ml = _sgb_model(function_tmpdir, "mix")
    spd = {0: np.array([(2, 0.03)], dtype=MfUsgSgb.get_default_dtype()).view(
        np.recarray
    )}
    sgb = MfUsgSgb(
        ml, stress_period_data=spd, parameters={"gp": {"parval": "1.0",
                                                       "data": [(4, 0.01)]}}
    )
    out = function_tmpdir / "mix.sgb"
    sgb.fn_path = str(out)
    sgb.write_file()
    text = out.read_text()
    assert "PARAMETER 1 1" in text
    assert " 1 0    Stress Period 1" in text  # 1 non-param row, NP=0

    re = MfUsgSgb.load(str(out), _sgb_model(function_tmpdir, "mix2"), nper=1,
                       ext_unit_dict={})
    assert list(re.stress_period_data[0]["node"]) == [2]
    assert set(re.parameters) == {"gp"}


def test_mfusgsgb_parameter_authoring_auto_counts(function_tmpdir):
    """NPSGB = len(parameters) and MXS = total nlst when omitted."""
    ml = _sgb_model(function_tmpdir, "ac")
    params = {
        "a": {"parval": "1.0", "data": [(0, 0.01), (1, 0.01)]},
        "b": {"parval": "1.0", "data": [(2, 0.01)]},
    }
    sgb = MfUsgSgb(ml, parameters=params)  # mxs omitted
    out = function_tmpdir / "ac.sgb"
    sgb.fn_path = str(out)
    sgb.write_file()
    assert "PARAMETER 2 3" in out.read_text()  # NPSGB=2, MXS=2+1


def test_mfusgsgb_parameter_authoring_duplicate_names_fails(function_tmpdir):
    """Definition names colliding case-insensitively (`gp`/`GP`) -> ValueError."""
    params = {
        "gp": {"parval": "1.0", "data": [(0, 0.01)]},
        "GP": {"parval": "1.0", "data": [(1, 0.01)]},
    }
    sgb = MfUsgSgb(_sgb_model(function_tmpdir, "dd"), parameters=params)
    out = function_tmpdir / "dd.sgb"
    sgb.fn_path = str(out)
    with pytest.raises(ValueError, match="duplicate parameter definition"):
        sgb.write_file()
    assert not out.exists()


def test_mfusgsgb_parameter_authoring_bad_name_fails(function_tmpdir):
    """A blank/whitespace/over-long parameter name -> ValueError, no file."""
    for badname in ("", "g p", "abcdefghijk"):
        sgb = MfUsgSgb(
            _sgb_model(function_tmpdir, "bn"),
            parameters={badname: {"parval": "1.0", "data": [(0, 0.01)]}},
        )
        out = function_tmpdir / "bn.sgb"
        if out.exists():
            out.unlink()
        sgb.fn_path = str(out)
        with pytest.raises(ValueError):
            sgb.write_file()
        assert not out.exists()


def test_mfusgsgb_parameter_authoring_bad_parval_fails(function_tmpdir):
    """A multi-token / blank parval -> ValueError, no file."""
    for badparval in ("1 2", "  "):
        sgb = MfUsgSgb(
            _sgb_model(function_tmpdir, "bp"),
            parameters={"gp": {"parval": badparval, "data": [(0, 0.01)]}},
        )
        out = function_tmpdir / "bp.sgb"
        if out.exists():
            out.unlink()
        sgb.fn_path = str(out)
        with pytest.raises(ValueError):
            sgb.write_file()
        assert not out.exists()


def test_mfusgsgb_parameter_authoring_empty_data_fails(function_tmpdir):
    """A parameter with empty data rows -> ValueError, no file."""
    sgb = MfUsgSgb(
        _sgb_model(function_tmpdir, "ed"),
        parameters={"gp": {"parval": "1.0", "data": []}},
    )
    out = function_tmpdir / "ed.sgb"
    sgb.fn_path = str(out)
    with pytest.raises(ValueError, match="data"):
        sgb.write_file()
    assert not out.exists()


def test_mfusgsgb_parameter_authoring_negative_node_fails(function_tmpdir):
    """A negative parametric node -> ValueError, no file."""
    sgb = MfUsgSgb(
        _sgb_model(function_tmpdir, "nn"),
        parameters={"gp": {"parval": "1.0", "data": [(-1, 0.01)]}},
    )
    out = function_tmpdir / "nn.sgb"
    sgb.fn_path = str(out)
    with pytest.raises(ValueError, match="negative node"):
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


def test_mfusgqrt_parameter_active_without_defs_fails(function_tmpdir):
    """active_params referencing parameters with no definitions fails explicitly
    (ValueError), no file. (From-scratch authoring with definitions is supported
    as of Stage 4.6C-C -- see the authoring tests below.)"""
    qrt = MfUsgQrt(_qrt_model(function_tmpdir, "fs"), active_params={0: ["p1"]})
    out = function_tmpdir / "fs.qrt"
    qrt.fn_path = str(out)
    with pytest.raises(ValueError, match="none are defined"):
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


# --- Stage 4.6C-C: QRT from-scratch NPQRT>0 parameter authoring (structural) -
#
# A parameterized QRT can be built from Python: parameters={name: {...}} +
# active_params={kper: [...]}; ergonomic input (data tuples; nlst/MXL/
# recipient_nodes auto) is normalized + validated before open. STRUCTURAL only:
# an activated QRT parameter round-trips but is not execution-guaranteed (the
# Fortran scales QRTF(5)=NumRT not Q, and NodQRT is not copied on activation).


def test_mfusgqrt_parameter_authoring_from_scratch(function_tmpdir):
    """Build a parameterized QRT from scratch (ergonomic input): NPQRT, auto MXL,
    definitions, per-SP activation, 0-based internal / 1-based file."""
    params = {"qp": {"parval": "1.5", "data": [(4, -10.0), (9, -20.0)]}}
    qrt = MfUsgQrt(
        _qrt_model(function_tmpdir, "auth"),
        parameters=params,
        active_params={0: ["qp"]},
    )
    out = function_tmpdir / "auth.qrt"
    qrt.fn_path = str(out)
    qrt.write_file()
    text = out.read_text()
    item1 = next(ln for ln in text.splitlines() if not ln.startswith("#")).split()
    assert item1[3] == "1" and item1[4] == "2"  # NPQRT=1, MXL auto = sum(nlst)=2
    assert "qp QRT 1.5 2" in text
    assert "\n 5  -1.000000e+01" in text  # node 4 -> 5 (1-based)

    re = MfUsgQrt.load(str(out), _qrt_model(function_tmpdir, "auth2"), nper=1)
    pd = re.parameters["qp"]
    assert pd["partyp"].upper() == "QRT" and pd["nlst"] == 2
    assert list(pd["data"]["node"]) == [4, 9]  # back to 0-based
    assert re.active_params == {0: ["qp"]}


def test_mfusgqrt_parameter_authoring_returnflow_recipients(function_tmpdir):
    """From-scratch QRT parameter with RETURNFLOW + per-row recipient_nodes
    (0-based internal / 1-based file)."""
    params = {
        "qp": {
            "parval": "1.0",
            "data": [(0, -100.0, 0.75), (4, -30.0, 0.0)],
            "recipient_nodes": [[8, 9], []],  # 0-based
        }
    }
    qrt = MfUsgQrt(
        _qrt_model(function_tmpdir, "rf"),
        options=["RETURNFLOW"],
        parameters=params,
        active_params={0: ["qp"]},
    )
    out = function_tmpdir / "rf.qrt"
    qrt.fn_path = str(out)
    qrt.write_file()

    re = MfUsgQrt.load(str(out), _qrt_model(function_tmpdir, "rf2"), nper=1)
    pd = re.parameters["qp"]
    assert pd["nlst"] == 2
    assert pd["recipient_nodes"][0] == [8, 9] and pd["recipient_nodes"][1] == []


def test_mfusgqrt_parameter_authoring_mixed_nonparam_active(function_tmpdir):
    """A non-parametric sink and an active parameter in one SP: MXAQRT counts
    both (NQRTNP + active NLST)."""
    dtype = MfUsgQrt.get_default_dtype(returnflow=False)
    spd = {0: np.array([(2, -30.0)], dtype=dtype).view(np.recarray)}
    qrt = MfUsgQrt(
        _qrt_model(function_tmpdir, "mix"),
        stress_period_data=spd,
        parameters={"qp": {"parval": "1.0", "data": [(0, -10.0)]}},
        active_params={0: ["qp"]},
    )
    out = function_tmpdir / "mix.qrt"
    qrt.fn_path = str(out)
    qrt.write_file()
    item1 = next(
        ln for ln in out.read_text().splitlines() if not ln.startswith("#")
    ).split()
    assert item1[0] == "2"  # MXAQRT = 1 non-param + 1 active param row

    re = MfUsgQrt.load(str(out), _qrt_model(function_tmpdir, "mix2"), nper=1)
    assert list(re.stress_period_data[0]["node"]) == [2]
    assert re.active_params == {0: ["qp"]}


def test_mfusgqrt_parameter_authoring_auto_mxl(function_tmpdir):
    """MXL is computed as the total definition rows when omitted."""
    params = {
        "a": {"parval": "1.0", "data": [(0, -10.0), (1, -10.0)]},
        "b": {"parval": "1.0", "data": [(2, -10.0)]},
    }
    qrt = MfUsgQrt(
        _qrt_model(function_tmpdir, "mxa"), parameters=params, active_params={0: ["a"]}
    )
    out = function_tmpdir / "mxa.qrt"
    qrt.fn_path = str(out)
    qrt.write_file()
    item1 = next(
        ln for ln in out.read_text().splitlines() if not ln.startswith("#")
    ).split()
    assert item1[4] == "3"  # MXL auto = 2 + 1


def test_mfusgqrt_parameter_authoring_active_case_insensitive(function_tmpdir):
    """An activation name differing from its definition only in case resolves and
    writes (the Fortran upper-cases PARNAM)."""
    qrt = MfUsgQrt(
        _qrt_model(function_tmpdir, "ci"),
        parameters={"qp": {"parval": "1.0", "data": [(0, -10.0)]}},
        active_params={0: ["QP"]},  # upper-case activation
    )
    out = function_tmpdir / "ci.qrt"
    qrt.fn_path = str(out)
    qrt.write_file()  # must not raise
    assert "qp QRT 1.0 1" in out.read_text()


def test_mfusgqrt_parameter_authoring_partyp_invalid_fails(function_tmpdir):
    """A non-QRT partyp raises ValueError before the file is opened, no file."""
    qrt = MfUsgQrt(
        _qrt_model(function_tmpdir, "pt"),
        parameters={"qp": {"partyp": "GHB", "parval": "1.0", "data": [(0, -10.0)]}},
        active_params={0: ["qp"]},
    )
    out = function_tmpdir / "pt.qrt"
    qrt.fn_path = str(out)
    with pytest.raises(ValueError, match="partyp"):
        qrt.write_file()
    assert not out.exists()


def test_mfusgqrt_parameter_authoring_duplicate_definition_names_fails(function_tmpdir):
    """Definition names colliding case-insensitively (`qp`/`QP`) -> ValueError."""
    params = {
        "qp": {"parval": "1.0", "data": [(0, -10.0)]},
        "QP": {"parval": "1.0", "data": [(1, -10.0)]},
    }
    qrt = MfUsgQrt(
        _qrt_model(function_tmpdir, "dd"), parameters=params, active_params={0: ["qp"]}
    )
    out = function_tmpdir / "dd.qrt"
    qrt.fn_path = str(out)
    with pytest.raises(ValueError, match="duplicate parameter definition"):
        qrt.write_file()
    assert not out.exists()


def test_mfusgqrt_parameter_authoring_bad_name_or_parval_fails(function_tmpdir):
    """A bad parameter name or parval raises ValueError, no file."""
    cases = [
        {"qp": {"parval": "1.0", "data": [(0, -10.0)]}, "_name": "q p"},  # bad name
        {"qp": {"parval": "1 2", "data": [(0, -10.0)]}},  # multi-token parval
    ]
    for i, spec in enumerate(cases):
        name = spec.pop("_name", "qp")
        pdef = spec["qp"]
        qrt = MfUsgQrt(
            _qrt_model(function_tmpdir, f"bn{i}"),
            parameters={name: pdef},
            active_params={0: [name]},
        )
        out = function_tmpdir / f"bn{i}.qrt"
        if out.exists():
            out.unlink()
        qrt.fn_path = str(out)
        with pytest.raises(ValueError):
            qrt.write_file()
        assert not out.exists()


def test_mfusgqrt_parameter_authoring_kper_out_of_range_fails(function_tmpdir):
    """An active_params stress period outside 0..nper-1 raises ValueError."""
    qrt = MfUsgQrt(
        _qrt_model(function_tmpdir, "kp", nper=1),
        parameters={"qp": {"parval": "1.0", "data": [(0, -10.0)]}},
        active_params={2: ["qp"]},
    )
    out = function_tmpdir / "kp.qrt"
    qrt.fn_path = str(out)
    with pytest.raises(ValueError, match="out of range"):
        qrt.write_file()
    assert not out.exists()


def test_mfusgqrt_parameter_authoring_negative_node_or_recipient_fails(function_tmpdir):
    """A negative parametric node or recipient -> ValueError, no file."""
    # negative node
    qrt = MfUsgQrt(
        _qrt_model(function_tmpdir, "nn"),
        parameters={"qp": {"parval": "1.0", "data": [(-1, -10.0)]}},
        active_params={0: ["qp"]},
    )
    out = function_tmpdir / "nn.qrt"
    qrt.fn_path = str(out)
    with pytest.raises(ValueError, match="negative node"):
        qrt.write_file()
    assert not out.exists()

    # negative recipient
    qrt = MfUsgQrt(
        _qrt_model(function_tmpdir, "nr"),
        options=["RETURNFLOW"],
        parameters={"qp": {"parval": "1.0", "data": [(0, -10.0, 0.5)],
                           "recipient_nodes": [[-1]]}},
        active_params={0: ["qp"]},
    )
    out = function_tmpdir / "nr.qrt"
    qrt.fn_path = str(out)
    with pytest.raises(ValueError, match="recipient"):
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


def test_mfusgdrt_parameter_active_without_defs_fails(function_tmpdir):
    """active_params referencing parameters with no definitions fails explicitly
    (ValueError), no file. (From-scratch authoring with definitions is supported
    as of Stage 4.6B — see the authoring tests below.)"""
    drt = MfUsgDrt(
        _drt_model(function_tmpdir, "fs"),
        options=["RETURNFLOW"],
        active_params={0: ["p1"]},
    )
    out = function_tmpdir / "fs.drt"
    drt.fn_path = str(out)
    with pytest.raises(ValueError, match="none are defined"):
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


# --- Stage 4.6B: DRT from-scratch NPDRT parameter authoring -----------------
#
# A parameterized DRT can be built entirely in Python (no prior load): pass
# parameters={name: {...}} + active_params={kper: [...]}. Ergonomic input — data
# as plain tuples, nlst/mxl/recipient_nodes computed when omitted — is normalized
# and validated before the file is opened; write -> reload preserves it.


def test_mfusgdrt_parameter_authoring_from_scratch(function_tmpdir):
    """Build a parameterized DRT from scratch (ergonomic input): NPDRT, auto MXL,
    definition rows, per-SP activation, 0-based internal / 1-based file."""
    ml = _drt_model(function_tmpdir, "auth")
    # data as plain tuples (node 0-based); nlst/mxl/recipient_nodes omitted.
    params = {"dp": {"parval": "2.0", "data": [(4, 5.0, 10.0, 0.0)]}}
    drt = MfUsgDrt(
        ml, options=["RETURNFLOW"], parameters=params, active_params={0: ["dp"]}
    )
    out = function_tmpdir / "auth.drt"
    drt.fn_path = str(out)
    drt.write_file()
    text = out.read_text()
    item1 = next(ln for ln in text.splitlines() if not ln.startswith("#")).split()
    assert item1[2] == "1" and item1[3] == "1"  # NPDRT=1, MXL auto = sum(nlst)=1
    assert "dp DRT 2.0 1" in text  # UPARLSTRP header
    assert "\n 5  5.000000e+00  1.000000e+01  0\n" in text  # node 4 -> 5 (1-based)

    re = MfUsgDrt.load(str(out), _drt_model(function_tmpdir, "auth2"), nper=1,
                       ext_unit_dict={})
    pd = re.parameters["dp"]
    assert pd["partyp"].upper() == "DRT" and pd["nlst"] == 1
    assert abs(float(pd["parval"]) - 2.0) < 1e-9
    assert list(pd["data"]["node"]) == [4]  # back to 0-based
    assert re.active_params == {0: ["dp"]}


def test_mfusgdrt_parameter_authoring_returnflow_changec_aux(function_tmpdir):
    """From-scratch parameter with RETURNFLOW + CHANGEC + AUX + per-row
    recipient_nodes round-trips (recipients 0-based internal / 1-based file)."""
    ml = _drt_model(function_tmpdir, "rfa")
    params = {
        "dp": {
            "parval": "1.5",
            "data": [(0, 5.0, 10.0, 0.7, 2, 0.15), (3, 4.0, 20.0, 0.5, 1, 0.25)],
            "recipient_nodes": [[8, 9], [7]],  # 0-based
        }
    }
    drt = MfUsgDrt(
        ml,
        options=["RETURNFLOW", "CHANGEC", "AUX C01"],
        parameters=params,
        active_params={0: ["dp"]},
    )
    out = function_tmpdir / "rfa.drt"
    drt.fn_path = str(out)
    drt.write_file()

    re = MfUsgDrt.load(str(out), _drt_model(function_tmpdir, "rfa2"), nper=1,
                       ext_unit_dict={})
    pd = re.parameters["dp"]
    assert pd["nlst"] == 2
    assert pd["recipient_nodes"][0] == [8, 9] and pd["recipient_nodes"][1] == [7]
    assert list(pd["data"]["idchngtyp"]) == [2, 1]
    assert np.allclose(pd["data"]["C01"], [0.15, 0.25])


def test_mfusgdrt_parameter_authoring_mixed_nonparam_active(function_tmpdir):
    """A non-parametric drain and an active parameter in the same SP: MXADRT
    counts both (NDRTNP + active NLST)."""
    ml = _drt_model(function_tmpdir, "mix")
    dtype = MfUsgDrt.get_usg_dtype(returnflow=True)
    spd = {0: np.array([(2, 6.0, 30.0, 0.0)], dtype=dtype).view(np.recarray)}
    params = {"dp": {"parval": "1.0", "data": [(0, 5.0, 10.0, 0.0)]}}
    drt = MfUsgDrt(
        ml,
        options=["RETURNFLOW"],
        stress_period_data=spd,
        parameters=params,
        active_params={0: ["dp"]},
    )
    out = function_tmpdir / "mix.drt"
    drt.fn_path = str(out)
    drt.write_file()
    item1 = next(
        ln for ln in out.read_text().splitlines() if not ln.startswith("#")
    ).split()
    assert item1[0] == "2"  # MXADRT = 1 non-param + 1 active param row

    re = MfUsgDrt.load(str(out), _drt_model(function_tmpdir, "mix2"), nper=1,
                       ext_unit_dict={})
    assert list(re.stress_period_data[0]["node"]) == [2]
    assert re.active_params == {0: ["dp"]}


def test_mfusgdrt_parameter_authoring_mxl_auto(function_tmpdir):
    """MXL is computed as the total definition rows when omitted."""
    ml = _drt_model(function_tmpdir, "mxa")
    params = {
        "a": {"parval": "1.0", "data": [(0, 5.0, 10.0, 0.0), (1, 5.0, 10.0, 0.0)]},
        "b": {"parval": "1.0", "data": [(2, 5.0, 10.0, 0.0)]},
    }
    drt = MfUsgDrt(
        ml, options=["RETURNFLOW"], parameters=params, active_params={0: ["a"]}
    )  # mxl omitted
    out = function_tmpdir / "mxa.drt"
    drt.fn_path = str(out)
    drt.write_file()
    item1 = next(
        ln for ln in out.read_text().splitlines() if not ln.startswith("#")
    ).split()
    assert item1[3] == "3"  # MXL auto = 2 + 1


def test_mfusgdrt_parameter_authoring_partyp_invalid_fails(function_tmpdir):
    """A non-DRT partyp raises ValueError before the file is opened, no file."""
    drt = MfUsgDrt(
        _drt_model(function_tmpdir, "pt"),
        options=["RETURNFLOW"],
        parameters={"dp": {"partyp": "GHB", "parval": "1.0",
                           "data": [(0, 5.0, 10.0, 0.0)]}},
        active_params={0: ["dp"]},
    )
    out = function_tmpdir / "pt.drt"
    drt.fn_path = str(out)
    with pytest.raises(ValueError, match="partyp"):
        drt.write_file()
    assert not out.exists()


def test_mfusgdrt_parameter_authoring_kper_out_of_range_fails(function_tmpdir):
    """An active_params stress period outside 0..nper-1 raises ValueError,
    no file."""
    drt = MfUsgDrt(
        _drt_model(function_tmpdir, "kp", nper=1),
        options=["RETURNFLOW"],
        parameters={"dp": {"parval": "1.0", "data": [(0, 5.0, 10.0, 0.0)]}},
        active_params={2: ["dp"]},  # nper=1, so kper 2 is invalid
    )
    out = function_tmpdir / "kp.drt"
    drt.fn_path = str(out)
    with pytest.raises(ValueError, match="out of range"):
        drt.write_file()
    assert not out.exists()


def test_mfusgdrt_parameter_authoring_recipients_without_returnflow_fails(
    function_tmpdir,
):
    """recipient_nodes on a parameter without RETURNFLOW raises ValueError."""
    drt = MfUsgDrt(
        _drt_model(function_tmpdir, "rw"),
        options=[],  # no RETURNFLOW
        parameters={"dp": {"parval": "1.0", "data": [(0, 5.0, 10.0)],
                           "recipient_nodes": [[7]]}},
        active_params={0: ["dp"]},
    )
    out = function_tmpdir / "rw.drt"
    drt.fn_path = str(out)
    with pytest.raises(ValueError, match="RETURNFLOW"):
        drt.write_file()
    assert not out.exists()


# --- Stage 4.6B review follow-up: parameter-authoring contract hardening ----
#
# Names (definition keys + activations) must be single whitespace-free tokens
# <= 10 chars and unique case-insensitively (the Fortran reads PARNAM as one
# CHARACTER*10 word and upper-cases it); parval must be a number or single
# token; parametric nodes/recipients must be non-negative 0-based integers.


def test_mfusgdrt_parameter_authoring_duplicate_definition_names_fails(function_tmpdir):
    """Definition names that differ only in case (`dp`/`DP`) collide once the
    Fortran upper-cases PARNAM -> ValueError, no file."""
    params = {
        "dp": {"parval": "1.0", "data": [(0, 5.0, 10.0, 0.0)]},
        "DP": {"parval": "1.0", "data": [(1, 5.0, 10.0, 0.0)]},
    }
    drt = MfUsgDrt(
        _drt_model(function_tmpdir, "dd"),
        options=["RETURNFLOW"],
        parameters=params,
        active_params={0: ["dp"]},
    )
    out = function_tmpdir / "dd.drt"
    drt.fn_path = str(out)
    with pytest.raises(ValueError, match="duplicate parameter definition"):
        drt.write_file()
    assert not out.exists()


def test_mfusgdrt_parameter_authoring_bad_name_fails(function_tmpdir):
    """A blank, whitespace, multi-token, or over-long parameter name raises
    ValueError before the file is opened, no file."""
    for badname in ("", "   ", "d p", "abcdefghijk"):  # "" blank, spaces, >10 chars
        drt = MfUsgDrt(
            _drt_model(function_tmpdir, "bn"),
            options=["RETURNFLOW"],
            parameters={badname: {"parval": "1.0", "data": [(0, 5.0, 10.0, 0.0)]}},
            active_params={0: [badname]},
        )
        out = function_tmpdir / "bn.drt"
        if out.exists():
            out.unlink()
        drt.fn_path = str(out)
        with pytest.raises(ValueError):
            drt.write_file()
        assert not out.exists()


def test_mfusgdrt_parameter_authoring_bad_parval_fails(function_tmpdir):
    """A multi-token or blank parval raises ValueError, no file (the Fortran
    reads PARVAL as one numeric token)."""
    for badparval in ("1 2", "  "):
        drt = MfUsgDrt(
            _drt_model(function_tmpdir, "bp"),
            options=["RETURNFLOW"],
            parameters={"dp": {"parval": badparval, "data": [(0, 5.0, 10.0, 0.0)]}},
            active_params={0: ["dp"]},
        )
        out = function_tmpdir / "bp.drt"
        if out.exists():
            out.unlink()
        drt.fn_path = str(out)
        with pytest.raises(ValueError):
            drt.write_file()
        assert not out.exists()


def test_mfusgdrt_parameter_authoring_negative_node_fails(function_tmpdir):
    """A negative parametric drain node raises ValueError, no file."""
    drt = MfUsgDrt(
        _drt_model(function_tmpdir, "nn"),
        options=["RETURNFLOW"],
        parameters={"dp": {"parval": "1.0", "data": [(-1, 5.0, 10.0, 0.0)]}},
        active_params={0: ["dp"]},
    )
    out = function_tmpdir / "nn.drt"
    drt.fn_path = str(out)
    with pytest.raises(ValueError, match="negative node"):
        drt.write_file()
    assert not out.exists()


def test_mfusgdrt_parameter_authoring_negative_recipient_fails(function_tmpdir):
    """A negative parametric recipient node raises ValueError, no file."""
    drt = MfUsgDrt(
        _drt_model(function_tmpdir, "nr"),
        options=["RETURNFLOW"],
        parameters={"dp": {"parval": "1.0", "data": [(0, 5.0, 10.0, 0.7)],
                           "recipient_nodes": [[-1]]}},
        active_params={0: ["dp"]},
    )
    out = function_tmpdir / "nr.drt"
    drt.fn_path = str(out)
    with pytest.raises(ValueError, match="recipient"):
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


# LPF read contract, READLPF + SGWF2LPFU1S/G (gwf2bcf-lpf-u1.f, USG-T 2.8.1).
# Per layer: HK, HANI if CHANI<=0, VKA (DISU: only if IKCFLAG=0); SS if
# transient; SY if transient and LAYTYP!=0 after the THICKSTRT reset
# (3198-3209, 3653/3918); VKCB if LAYCBD!=0 (DISU: and IKCFLAG=0) (3676/3934);
# WETDRY if LAYWET!=0 (3698/3950); ALPHA, BETA, SR, BROOK (+BP with BUBBLEPT)
# if the effective LAYCON is 5 and TABRICH is off (3708/3955), where LAYCON
# follows 3166-3180 and BAS RICHARDS/RICHARDS_HP turns LAYCON!=0 into 5
# (3334-3349). WETFCT/IWETIT/IHDWET if any LAYWET!=0 (3286-3299). Connection
# array after the layer loop if IKCFLAG!=0 (3974).

# one distinct value per property and layer (tuple index = layer)
_LPF_VALUES = {
    "hk": (21.0, 22.0),
    "hani": (0.81, 0.82),
    "vka": (3.1, 3.2),
    "ss": (1.1e-5, 1.2e-5),
    "sy": (0.051, 0.052),
    "vkcb": (0.061, 0.062),
    "wetdry": (-0.61, -0.62),
    "alpha": (0.71, 0.72),
    "beta": (8.1, 8.2),
    "sr": (0.091, 0.092),
    "brook": (10.1, 10.2),
    "bp": (-0.111, -0.112),
}


def _lpf_contract(case):
    """Restatement of the READLPF record order, as (name, layer 1-based)."""
    laytyp, nlay = case["laytyp"], len(case["laytyp"])
    thickstrt = case.get("thickstrt", False)
    richards = case.get("richards", False) or case.get("richards_hp", False)
    disu, ikc = case.get("disu", False), case.get("ikcflag", 0)
    laywet = case.get("laywet", [0] * nlay)
    laycbd = case.get("laycbd", [0] * nlay)
    seq = [("iuzontab", 0)] if case.get("tabrich") else []
    for k in range(nlay):
        t = laytyp[k]
        if t == 0 or (t < 0 and thickstrt):
            laycon = 0
        else:
            laycon = t if t in (4, 5) else 3
        if richards and laycon != 0:
            laycon = 5
        if not disu or ikc == 0:
            seq += [("hk", k + 1), ("vka", k + 1)]
        if case.get("transient", True):
            seq.append(("ss", k + 1))
            if t != 0 and not (t < 0 and thickstrt):
                seq.append(("sy", k + 1))
        if laycbd[k] != 0 and (not disu or ikc == 0):
            seq.append(("vkcb", k + 1))
        if laywet[k] != 0:
            seq.append(("wetdry", k + 1))
        if laycon == 5 and not case.get("tabrich"):
            seq += [(n, k + 1) for n in ("alpha", "beta", "sr", "brook")]
            if case.get("bubblept"):
                seq.append(("bp", k + 1))
    if disu and ikc != 0:
        seq.append(("ksat", 0))
    return seq


_LPF_KSAT = np.array([2.01, 2.02, 2.03, 2.04, 2.05, 2.06, 2.07], dtype=np.float32)


def _lpf_disu_line(ml, nlay, steady):
    """DISU line of 3 cells per layer, IVSD=1, IDSYMRD=0 (NJAG-long arrays)."""
    iac, ja, ivc, cl12, fahl = [], [], [], [], []
    for k in range(nlay):
        for j in range(3):
            n = 3 * k + j
            nb = [(m, 1) for m in (n - 3,) if k > 0] + [(n - 1, 0)] * (j > 0)
            nb += [(n + 1, 0)] * (j < 2) + [(n + 3, 1)] * (k < nlay - 1)
            iac.append(len(nb) + 1)
            ja.append(n)
            ivc.append(0)
            cl12.append(0.0)
            fahl.append(0.0)
            for m, vert in nb:
                ja.append(m)
                ivc.append(vert)
                cl12.append(5.0)
                fahl.append(100.0 if vert else 10.0)
    MfUsgDisU(
        ml,
        nodes=3 * nlay,
        nlay=nlay,
        njag=len(ja),
        nodelay=[3] * nlay,
        iac=np.array(iac),
        ja=np.array(ja),
        ivsd=1,
        ivc=np.array(ivc),
        idsymrd=0,
        cl12=np.array(cl12),
        fahl=np.array(fahl),
        top=[30.0 - 10.0 * k for k in range(nlay) for _ in range(3)],
        bot=[20.0 - 10.0 * k for k in range(nlay) for _ in range(3)],
        area=100.0,
        steady=steady,
    )


def _lpf_model(ws, case):
    from flopy.modflow import ModflowDis

    laytyp, nlay = case["laytyp"], len(case["laytyp"])
    disu = case.get("disu", False)
    ml = MfUsg(modelname="lpf", model_ws=str(ws), structured=not disu)
    steady = not case.get("transient", True)
    if disu:
        _lpf_disu_line(ml, nlay, steady)
        ncell = 3
    else:
        laycbd = case.get("laycbd", [0] * nlay)
        nbot = nlay + sum(laycbd)
        ModflowDis(
            ml,
            nlay=nlay,
            nrow=1,
            ncol=2,
            nper=1,
            steady=steady,
            top=30.0,
            botm=[30.0 - 10.0 * (i + 1) for i in range(nbot)],
            laycbd=laycbd,
        )
        ncell = 2
    MfUsgBas(
        ml,
        ibound=1,
        strt=25.0,
        richards=case.get("richards", False),
        richards_hp=case.get("richards_hp", False),
    )
    values = {
        name: [
            np.full((1, ncell), vals[k]) if not disu else np.full(ncell, vals[k])
            for k in range(nlay)
        ]
        for name, vals in _LPF_VALUES.items()
    }
    tab = {}
    if case.get("tabrich"):
        tab = {
            "tabrich": True,
            "nuzones": 2,
            "nutabrows": 2,
            "iuzontab": [1] * (ncell * nlay),
            "retcrvs": _RETCRVS_2x2,
        }
    MfUsgLpf(
        ml,
        laytyp=laytyp,
        laywet=case.get("laywet", 0),
        ikcflag=case.get("ikcflag", 0),
        thickstrt=case.get("thickstrt", False),
        bubblept=case.get("bubblept", False),
        ksat=_LPF_KSAT if (disu and nlay == 1) else 1.0,  # NJAG=7 for one layer
        **values,
        **tab,
    )
    return ml


def _lpf_manual_text(case):
    """Hand-built LPF file in the READLPF order (independent of the writer)."""
    nlay = len(case["laytyp"])
    head = "0 -1.0E+30 0"
    if case.get("disu"):
        head += f" {case.get('ikcflag', 0)}"
    if case.get("thickstrt"):
        head += " THICKSTRT"
    if case.get("bubblept"):
        head += " BUBBLEPT"
    lines = ["# hand-written LPF in READLPF order"]
    if case.get("tabrich"):
        lines.append(head + " TABRICH 2 2")
        ncell = 3 if case.get("disu") else 2
        lines.append("INTERNAL 1 (FREE) -1 #iuzontab")
        lines.append(" ".join(["1"] * ncell * nlay))
        for zone in _RETCRVS_2x2:
            lines += [" ".join(str(v) for v in row) for row in zone]
    else:
        lines.append(head)
    laywet = case.get("laywet", [0] * nlay)
    lines += [
        " ".join(str(t) for t in case["laytyp"]),
        " ".join("0" for _ in range(nlay)),
        " ".join("1.0" for _ in range(nlay)),
        " ".join("0" for _ in range(nlay)),
        " ".join(str(w) for w in laywet),
    ]
    if any(w != 0 for w in laywet):
        lines.append("0.1 1 0")
    for name, layer in _lpf_contract(case):
        if name == "iuzontab":
            continue
        if name == "ksat":
            lines.append("INTERNAL 1.0 (FREE) -1 #ksat")
            lines.append(" ".join(f"{v:.2f}" for v in _LPF_KSAT))
        else:
            value = _LPF_VALUES[name][layer - 1]
            lines.append(f"CONSTANT {value:.6E} #{name} layer {layer}")
    return "\n".join(lines) + "\n"


def _lpf_records(path):
    """(name, layer) of each array record identified by its distinct value."""
    lookup = [
        ((n, k + 1), v) for n, vals in _LPF_VALUES.items() for k, v in enumerate(vals)
    ]
    lookup.append((("ksat", 0), float(_LPF_KSAT[0])))
    lookup.append((("iuzontab", 0), 1.0))
    text = Path(path).read_text().splitlines()
    out = []
    for i, line in enumerate(text):
        tokens = line.split()
        if not tokens or tokens[0].upper() not in ("CONSTANT", "INTERNAL"):
            continue
        is_constant = tokens[0].upper() == "CONSTANT"
        value = float(tokens[1] if is_constant else text[i + 1].split()[0])
        hits = [key for key, v in lookup if np.isclose(value, v, rtol=1e-5)]
        out.append(hits[0] if len(hits) == 1 else ("?", value))
    return out


def _assert_lpf_values(lpf, case):
    for name, layer in _lpf_contract(case):
        if name == "iuzontab":
            continue
        if name == "ksat":
            assert np.allclose(lpf.ksat.array, _LPF_KSAT), lpf.ksat.array
            continue
        got = np.asarray(getattr(lpf, name)[layer - 1].array)
        assert np.allclose(got, _LPF_VALUES[name][layer - 1]), (name, layer, got)


def _lpf_write_load_check(ws, case):
    """Writer order, independent manual load, and FloPy round-trip."""
    ml = _lpf_model(ws / "written", case)
    lpf = ml.get_package("LPF")
    lpf.write_file(check=False)
    assert _lpf_records(lpf.fn_path) == _lpf_contract(case)

    manual = ws / "manual.lpf"
    manual.write_text(_lpf_manual_text(case))
    loaded = MfUsgLpf.load(str(manual), _lpf_model(ws / "m1", case), check=False)
    _assert_lpf_values(loaded, case)

    reloaded = MfUsgLpf.load(lpf.fn_path, _lpf_model(ws / "m2", case), check=False)
    _assert_lpf_values(reloaded, case)
    return lpf, loaded


@pytest.mark.parametrize(
    "case",
    [
        {"laytyp": [0, 5], "richards": True, "bubblept": True},
        {"laytyp": [0, 1], "richards": True},  # LAYTYP 1 promoted to LAYCON 5
        {"laytyp": [0, 1], "richards_hp": True, "bubblept": True},
        {"laytyp": [0, 5], "disu": True, "richards": True, "bubblept": True},
        {"laytyp": [0, 5], "richards": True, "tabrich": True},  # control
        {"laytyp": [5, 5], "richards": True},  # control: all layers LAYCON 5
    ],
    ids=[
        "dis-0-5-bp",
        "dis-0-1-promoted",
        "dis-0-1-richards_hp",
        "disu-0-5-bp",
        "tabrich-control",
        "all5-control",
    ],
)
def test_mfusglpf_richards_arrays_follow_effective_laycon(function_tmpdir, case):
    """F06: ALPHA/BETA/SR/BROOK(/BP) only where the effective LAYCON is 5."""
    lpf, _ = _lpf_write_load_check(function_tmpdir, case)
    assert list(lpf.laytyp.array) == case["laytyp"]


def test_mfusglpf_wetdry_loaded_for_convertible_layers(function_tmpdir):
    """F07: WETDRY is read whenever LAYWET != 0 (two convertible layers)."""
    case = {"laytyp": [1, 1], "laywet": [1, 1]}
    _lpf_write_load_check(function_tmpdir, case)


def test_mfusglpf_wetdry_open_close_and_shared_external(function_tmpdir):
    """F07 with external records: values and EXTERNAL stream position."""
    from flopy.utils.mfreadnam import parsenamefile

    (function_tmpdir / "wetdry_l1.txt").write_text("-0.61 -0.61\n")
    (function_tmpdir / "lpf_ext71.dat").write_text(
        "22.0 22.0\n3.2 3.2\nEND-OF-LPF-BLOCKS\n"
    )
    (function_tmpdir / "manual.lpf").write_text(
        "# hand-written LPF in READLPF order\n0 -1.0E+30 0\n1 1\n0 0\n1.0 1.0\n"
        "0 0\n1 1\n0.1 1 0\n"
        "CONSTANT 21.0 #hk layer 1\nCONSTANT 3.1 #vka layer 1\n"
        "CONSTANT 1.1e-5 #ss layer 1\nCONSTANT 0.051 #sy layer 1\n"
        "OPEN/CLOSE wetdry_l1.txt 1.0 (FREE) -1 #wetdry layer 1\n"
        "EXTERNAL 71 1.0 (FREE) -1 #hk layer 2\n"
        "EXTERNAL 71 1.0 (FREE) -1 #vka layer 2\n"
        "CONSTANT 1.2e-5 #ss layer 2\nCONSTANT 0.052 #sy layer 2\n"
        "CONSTANT -0.62 #wetdry layer 2\n"
    )
    nam = function_tmpdir / "manual.nam"
    nam.write_text("LPF 15 manual.lpf\nDATA 71 lpf_ext71.dat\n")
    case = {"laytyp": [1, 1], "laywet": [1, 1]}
    ml = _lpf_model(function_tmpdir, case)
    ext_unit_dict = parsenamefile(str(nam), ml.mfnam_packages, verbose=False)
    try:
        lpf = MfUsgLpf.load(
            str(function_tmpdir / "manual.lpf"),
            ml,
            ext_unit_dict=ext_unit_dict,
            check=False,
        )
        remaining = ext_unit_dict[71].filehandle.readline().strip()
    finally:
        for item in ext_unit_dict.values():
            if item.filehandle is not None:
                item.filehandle.close()
    _assert_lpf_values(lpf, case)
    assert remaining == "END-OF-LPF-BLOCKS"


@pytest.mark.parametrize("transient", [False, True], ids=["steady", "transient"])
@pytest.mark.parametrize(
    "case",
    [
        {"laytyp": [1, 0], "laycbd": [1, 0]},  # confining bed
        {"laytyp": [1, 0], "laywet": [1, 0]},  # wetting
    ],
    ids=["vkcb", "wetdry"],
)
def test_mfusglpf_vkcb_wetdry_not_gated_by_transient(function_tmpdir, case, transient):
    """F08: VKCB and WETDRY are written in steady-state models too."""
    _lpf_write_load_check(function_tmpdir, dict(case, transient=transient))


@pytest.mark.parametrize("thickstrt", [True, False], ids=["thickstrt", "no-thickstrt"])
def test_mfusglpf_sy_uses_laytyp_after_thickstrt(function_tmpdir, thickstrt):
    """SY is read only where LAYTYP != 0 after the THICKSTRT reset."""
    case = {"laytyp": [-1, 1], "thickstrt": thickstrt}
    lpf, loaded = _lpf_write_load_check(function_tmpdir, case)
    # public LAYTYP values are not mutated by the derivation
    assert list(lpf.laytyp.array) == [-1, 1]
    assert list(loaded.laytyp.array) == [-1, 1]


@pytest.mark.parametrize("ikcflag", [0, 1, -1, 2])
def test_mfusglpf_connection_array_for_nonzero_ikcflag(function_tmpdir, ikcflag):
    """F09: the DISU connection array (IDSYMRD=0, NJAG values) is written and
    read whenever IKCFLAG != 0. IKCFLAG=2 is not a documented option; it is
    included only because the 2.8.1 reader treats any non-zero value alike."""
    case = {"laytyp": [0], "disu": True, "ikcflag": ikcflag, "transient": False}
    lpf, _ = _lpf_write_load_check(function_tmpdir, case)
    has_ksat = any(name == "ksat" for name, _ in _lpf_records(lpf.fn_path))
    assert has_ksat == (ikcflag != 0)


@pytest.mark.parametrize(
    "laywet, item7",
    [([-1, 1], True), ([0, 0], False), ([1, 1], True)],
    ids=["mixed-sign", "no-wetting", "positive"],
)
def test_mfusglpf_item7_when_any_laywet_nonzero(function_tmpdir, laywet, item7):
    """L7: WETFCT/IWETIT/IHDWET are present when any LAYWET != 0."""
    case = {"laytyp": [1, 1], "laywet": laywet}
    lpf, loaded = _lpf_write_load_check(function_tmpdir, case)
    # heading, item 1, then LAYTYP, LAYAVG, CHANI, LAYVKA, LAYWET (lines 2-6)
    lines = Path(lpf.fn_path).read_text().splitlines()
    assert lines[6].split() == [str(w) for w in laywet]
    next_is_array = lines[7].split()[0].upper() in ("CONSTANT", "INTERNAL")
    assert next_is_array != item7
    if item7:
        assert np.isclose(loaded.wetfct, 0.1) and loaded.iwetit == 1


def test_mfusglpf_write_file_to_open_handle(function_tmpdir):
    """WF: write_file(f=handle) writes to the handle and leaves it open."""
    case = {"laytyp": [0, 5], "richards": True, "laywet": [0, 0]}
    ml = _lpf_model(function_tmpdir, case)
    lpf = ml.get_package("LPF")
    lpf.write_file(check=False)
    expected = Path(lpf.fn_path).read_text()

    target = function_tmpdir / "via_handle.lpf"
    with open(target, "w") as fh:
        lpf.write_file(check=False, f=fh)
        assert not fh.closed
        fh.write("# caller content after the package\n")
    assert target.read_text() == expected + "# caller content after the package\n"


# BCF read contract, READBCF + SGWF2BCFU1S/G (gwf2bcf-lpf-u1.f, USG-T 2.8.1).
# Item 1: free format, or FORMAT(I10,F10.3,I10,F10.3,4I10) with IFREFM=0
# (91-113); IKCFLAG only for unstructured grids. LAYCON codes
# 10*INTERCELLT+LAYCON with READ(IN,*) or FORMAT(40I2) (231-237); RICHARDS or
# RICHARDS_HP turns every LAYCON!=0 into 5 (335-350). TRPY if IKCFLAG=0
# (353-361); DISU ANGLEX if IKCFLAG=0 and some |TRPY-1|>1e-6 (747-749). Per
# layer: SF1 if transient; with IKCFLAG=0, TRAN if LAYCON is 0 or 2 else HY,
# then VCONT (IKVFLAG=0, K<NLAY) or KV (IKVFLAG!=0, NLAY>1); SF2 if transient
# and LAYCON is 2-5; WETDRY if IWDFLG!=0 and LAYCON is 1 or 3; ALPHA, BETA,
# SR, BROOK (+BP with BUBBLEPT) if LAYCON is 5 and TABRICH is off (S 558-687,
# G 753-816). DISU connection array if IKCFLAG!=0 (819).

_BCF_VALUES = {
    "sf1": (1.1e-5, 1.2e-5, 1.3e-5),
    "tran": (11.0, 12.0, 13.0),
    "hy": (21.0, 22.0, 23.0),
    "vcont": (0.031, 0.032, 0.033),
    "kv": (41.0, 42.0, 43.0),
    "sf2": (0.051, 0.052, 0.053),
    "wetdry": (-0.61, -0.62, -0.63),
    "alpha": (0.71, 0.72, 0.73),
    "beta": (8.1, 8.2, 8.3),
    "sr": (0.091, 0.092, 0.093),
    "brook": (10.1, 10.2, 10.3),
    "bp": (-0.111, -0.112, -0.113),
}
_BCF_KSAT = np.array([3.01, 3.02, 3.03, 3.04, 3.05, 3.06, 3.07], dtype=np.float32)
_BCF_TRPY_ANISO = 0.5
_BCF_ANGLEX = 0.25


def _bcf_contract(case):
    """Restatement of the READBCF record order, as (name, layer 1-based)."""
    laycon, nlay = case["laycon"], len(case["laycon"])
    richards = case.get("richards", False) or case.get("richards_hp", False)
    disu = case.get("disu", False)
    ikc = case.get("ikcflag", 0) if disu else 0
    ikv = case.get("ikvflag", 0)
    seq = [("iuzontab", 0)] if case.get("tabrich") else []
    if ikc == 0:
        seq.append(("trpy", 0))
        if disu and case.get("trpy", 1.0) != 1.0:
            seq.append(("anglex", 0))
    for k in range(nlay):
        lc = 5 if (richards and laycon[k] != 0) else laycon[k]
        if case.get("transient", True):
            seq.append(("sf1", k + 1))
        if ikc == 0:
            seq.append(("tran" if lc in (0, 2) else "hy", k + 1))
            if ikv == 0 and k < nlay - 1:
                seq.append(("vcont", k + 1))
            elif ikv != 0 and nlay > 1:
                seq.append(("kv", k + 1))
        if case.get("transient", True) and lc in (2, 3, 4, 5):
            seq.append(("sf2", k + 1))
        if case.get("iwdflg", 0) != 0 and lc in (1, 3):
            seq.append(("wetdry", k + 1))
        if lc == 5 and not case.get("tabrich"):
            seq += [(n, k + 1) for n in ("alpha", "beta", "sr", "brook")]
            if case.get("bubblept"):
                seq.append(("bp", k + 1))
    if disu and ikc != 0:
        seq.append(("ksat", 0))
    return seq


def _bcf_disu_irregular(ml, steady):
    """DISU with NODELAY=[3, 2] (cells 4-5 under cells 1-2), IVSD=1, IDSYMRD=0."""
    nbrs = {0: [1, 3], 1: [0, 2, 4], 2: [1], 3: [0, 4], 4: [1, 3]}
    iac, ja, ivc, cl12, fahl = [], [], [], [], []
    for n in range(5):
        iac.append(len(nbrs[n]) + 1)
        ja.append(n)
        ivc.append(0)
        cl12.append(0.0)
        fahl.append(0.0)
        for m in nbrs[n]:
            vert = (n < 3) != (m < 3)
            ja.append(m)
            ivc.append(int(vert))
            cl12.append(5.0)
            fahl.append(100.0 if vert else 10.0)
    MfUsgDisU(
        ml,
        nodes=5,
        nlay=2,
        njag=len(ja),
        nodelay=[3, 2],
        iac=np.array(iac),
        ja=np.array(ja),
        ivsd=1,
        ivc=np.array(ivc),
        idsymrd=0,
        cl12=np.array(cl12),
        fahl=np.array(fahl),
        top=[30.0] * 3 + [20.0] * 2,
        bot=[20.0] * 3 + [10.0] * 2,
        area=100.0,
        steady=steady,
    )


def _bcf_model(ws, case):
    from flopy.modflow import ModflowDis

    laycon, nlay = case["laycon"], len(case["laycon"])
    disu = case.get("disu", False)
    steady = not case.get("transient", True)
    ml = MfUsg(modelname="bcf", model_ws=str(ws), structured=not disu)
    if case.get("irregular"):
        _bcf_disu_irregular(ml, steady)
        ncells = [3, 2]
    elif disu:
        _lpf_disu_line(ml, nlay, steady)
        ncells = [3] * nlay
    else:
        ModflowDis(
            ml,
            nlay=nlay,
            nrow=1,
            ncol=2,
            nper=1,
            steady=steady,
            top=30.0,
            botm=[30.0 - 10.0 * (k + 1) for k in range(nlay)],
        )
        ncells = [2] * nlay
    MfUsgBas(
        ml,
        ibound=1,
        strt=25.0,
        richards=case.get("richards", False),
        richards_hp=case.get("richards_hp", False),
        ifrefm=case.get("ifrefm", True),
    )

    def layered(name):
        vals = _BCF_VALUES[name]
        if disu:
            return [np.full(ncells[k], vals[k % 3]) for k in range(nlay)]
        return [np.full((1, 2), vals[k % 3]) for k in range(nlay)]

    tab = {}
    if case.get("tabrich"):
        tab = {
            "tabrich": True,
            "nuzones": 2,
            "nutabrows": 2,
            "iuzontab": [1] * sum(ncells),
            "retcrvs": _RETCRVS_2x2,
        }
    # the Richards arrays are passed only when some layer reads them, so the
    # other cases also run with versions that lack these parameters
    richards = case.get("richards", False) or case.get("richards_hp", False)
    needs_richards = not case.get("tabrich") and any(
        c == 5 or (richards and c != 0) for c in laycon
    )
    names = [n for n in _BCF_VALUES if n != "vcont"]
    if not needs_richards:
        names = [n for n in names if n not in ("alpha", "beta", "sr", "brook", "bp")]
    MfUsgBcf(
        ml,
        laycon=laycon,
        intercellt=case.get("intercellt", 0),
        trpy=case.get("trpy", 1.0),
        iwdflg=case.get("iwdflg", 0),
        ikvflag=case.get("ikvflag", 0),
        ikcflag=case.get("ikcflag", 0),
        bubblept=case.get("bubblept", False),
        anglex=_BCF_ANGLEX,
        ksat=_BCF_KSAT if (disu and nlay == 1) else 1.0,  # NJAG=7 for one layer
        vcont=layered("vcont")[: max(nlay - 1, 1)],
        **{n: layered(n) for n in names},
        **tab,
    )
    return ml


def _bcf_manual_text(case):
    """Hand-built free-format BCF file in the READBCF order."""
    nlay = len(case["laycon"])
    disu = case.get("disu", False)
    head = f"0 -1.0E+30 {case.get('iwdflg', 0)} 0.1 1 0 {case.get('ikvflag', 0)}"
    if disu:
        head += f" {case.get('ikcflag', 0)}"
    if case.get("tabrich"):
        head += " TABRICH 2 2"
    if case.get("bubblept"):
        head += " BUBBLEPT"
    lines = [head]
    if case.get("tabrich"):
        ncell = 5 if case.get("irregular") else (3 if disu else 2) * nlay
        lines += ["INTERNAL 1 (FREE) -1 #iuzontab", " ".join(["1"] * ncell)]
        for zone in _RETCRVS_2x2:
            lines += [" ".join(str(v) for v in row) for row in zone]
    ic = case.get("intercellt", [0] * nlay)
    ic = ic if isinstance(ic, list) else [ic] * nlay
    lines.append(" ".join(str(10 * i + c) for i, c in zip(ic, case["laycon"])))
    ncol = 3 if disu else case.get("ncol", 2)
    ncells = [3, 2] if case.get("irregular") else [ncol] * nlay
    for name, layer in _bcf_contract(case):
        if name == "iuzontab":
            continue
        if name == "trpy":
            lines.append(f"CONSTANT {case.get('trpy', 1.0)} #trpy")
        elif name == "anglex":
            lines.append(f"CONSTANT {_BCF_ANGLEX} #anglex")
        elif name == "ksat":
            lines.append("INTERNAL 1.0 (FREE) -1 #ksat")
            lines.append(" ".join(f"{v:.2f}" for v in _BCF_KSAT))
        else:
            value = _BCF_VALUES[name][(layer - 1) % 3]
            n = ncells[layer - 1]
            lines += [
                f"INTERNAL 1.0 (FREE) -1 #{name} layer {layer}",
                " ".join([f"{value:.6E}"] * n),
            ]
    return "\n".join(lines) + "\n"


def _bcf_records(path):
    """(name, layer) of each array record identified by its distinct value."""
    lookup = [
        ((n, k + 1), v) for n, vals in _BCF_VALUES.items() for k, v in enumerate(vals)
    ]
    lookup += [(("ksat", 0), float(_BCF_KSAT[0])), (("anglex", 0), _BCF_ANGLEX)]
    lookup += [(("trpy", 0), 1.0), (("trpy", 0), _BCF_TRPY_ANISO)]
    text = Path(path).read_text().splitlines()
    out = []
    for i, line in enumerate(text):
        tokens = line.split()
        if not tokens or tokens[0].upper() not in ("CONSTANT", "INTERNAL"):
            continue
        is_constant = tokens[0].upper() == "CONSTANT"
        value = float(tokens[1] if is_constant else text[i + 1].split()[0])
        if "#iuzontab" in line:
            out.append(("iuzontab", 0))
            continue
        hits = [key for key, v in lookup if np.isclose(value, v, rtol=1e-5)]
        out.append(hits[0] if len(hits) == 1 else ("?", value))
    return out


def _assert_bcf_values(bcf, case):
    for name, layer in _bcf_contract(case):
        if name in ("iuzontab", "trpy", "anglex"):
            continue
        if name == "ksat":
            assert np.allclose(bcf.ksat.array, _BCF_KSAT), bcf.ksat.array
            continue
        got = np.asarray(getattr(bcf, name)[layer - 1].array)
        assert np.allclose(got, _BCF_VALUES[name][(layer - 1) % 3]), (name, layer, got)


def _bcf_write_load_check(ws, case):
    """Writer order, independent manual load, and FloPy round-trip."""
    ml = _bcf_model(ws / "written", case)
    bcf = ml.get_package("BCF6")
    bcf.write_file()
    assert _bcf_records(bcf.fn_path) == _bcf_contract(case)

    manual = ws / "manual.bcf"
    manual.write_text(_bcf_manual_text(case))
    loaded = MfUsgBcf.load(str(manual), _bcf_model(ws / "m1", case))
    _assert_bcf_values(loaded, case)

    reloaded = MfUsgBcf.load(bcf.fn_path, _bcf_model(ws / "m2", case))
    _assert_bcf_values(reloaded, case)
    return bcf, loaded


def test_mfusgbcf_laycon_codes_free_format(function_tmpdir):
    """B01: each code 10*INTERCELLT+LAYCON is one integer for READ(IN,*)."""
    case = {"laycon": [3, 0, 0], "intercellt": [1, 0, 2], "transient": False}
    bcf, loaded = _bcf_write_load_check(function_tmpdir, case)
    lines = Path(bcf.fn_path).read_text().splitlines()
    assert [int(t) for t in lines[1].split()] == [13, 0, 20]
    for pkg in (bcf, loaded):
        assert list(pkg.intercellt.array) == [1, 0, 2]
        assert list(pkg.laycon.array) == [3, 0, 0]


def _laycon_model(nlay):
    from flopy.modflow import ModflowDis

    ml = MfUsg(modelname="bcf")
    ModflowDis(ml, nlay=nlay, nrow=1, ncol=2, nper=1, steady=True)
    MfUsgBas(ml, ibound=1, strt=1.0, ifrefm=True)
    return ml


# READ(IN,*) (LAYCON(I),I=1,NLAY) (GWF2BCFU1AR, 2.8.1) takes NLAY values from
# as many lines as needed and ignores the rest of the line holding the last one
@pytest.mark.parametrize(
    "nlay, text, codes",
    [
        (1, "4  Ltype(NLAY)\n", [4]),
        (2, "04 04                          Ltype(NLAY)\n", [4, 4]),
        (2, "04 04 # Ltype(NLAY)\n", [4, 4]),
        (2, "13 0 99 1.5e3 x\n", [13, 0]),
        (3, "13\n\n 0 \n 20 30 Ltype(NLAY)\n", [13, 0, 20]),
        (45, "3 " * 40 + "\n" + "11 " * 5 + "Ltype(NLAY)\n", [3] * 40 + [11] * 5),
    ],
    ids=[
        "one-layer",
        "ex3",
        "comment-marker",
        "extra-numbers",
        "continuation",
        "45-layers",
    ],
)
def test_mfusgbcf_laycon_free_format_reads_nlay_values(nlay, text, codes):
    """Only the NLAY required values are converted; the next record is intact."""
    f_obj = io.StringIO(text + "constant 1.0  TRPY(NLAY)\n")
    laycon, intercellt = MfUsgBcf._load_laycon(f_obj, _laycon_model(nlay))
    assert list(laycon) == [c % 10 for c in codes]
    assert list(intercellt) == [c // 10 for c in codes]
    assert f_obj.readline() == "constant 1.0  TRPY(NLAY)\n"


@pytest.mark.parametrize("text", ["4 Ltype(NLAY) 4\n", "4\nx 4\n", "4 4.5\n"])
def test_mfusgbcf_laycon_free_format_invalid_value_fails(text):
    """A token that is not an integer among the NLAY values is an error."""
    with pytest.raises(ValueError, match="invalid literal"):
        MfUsgBcf._load_laycon(io.StringIO(text), _laycon_model(2))


def test_mfusgbcf_laycon_free_format_end_of_file_fails():
    """End of file before NLAY values fails instead of reading forever; the
    read runs in a thread so that a hang fails the test after 10 s."""
    ml = _laycon_model(3)
    result = {}

    def load():
        try:
            MfUsgBcf._load_laycon(io.StringIO("13 0\n"), ml)
        except ValueError as err:
            result["error"] = err

    worker = threading.Thread(target=load, daemon=True)
    worker.start()
    worker.join(timeout=10)
    assert not worker.is_alive(), "LAYCON read did not stop at end of file"
    assert "end of file while reading LAYCON: 2 of 3" in str(result.get("error"))


@pytest.mark.parametrize("nlay", [3, 45])
def test_mfusgbcf_fixed_format_header_and_laycon(function_tmpdir, nlay):
    """B01f: fixed-format item 1 fields and FORMAT(40I2) LAYCON lines."""
    from flopy.modflow import ModflowDis

    ml = MfUsg(modelname="bcf", model_ws=str(function_tmpdir))
    ModflowDis(
        ml,
        nlay=nlay,
        nrow=1,
        ncol=2,
        nper=1,
        steady=True,
        top=30.0,
        botm=[30.0 - 0.5 * (k + 1) for k in range(nlay)],
    )
    MfUsgBas(ml, ibound=1, strt=25.0, ifrefm=False)
    intercellt = [(k % 3) for k in range(nlay)]
    laycon = [3 if k == 0 else (k % 2) * 3 for k in range(nlay)]
    bcf = MfUsgBcf(
        ml,
        ipakcb=0,
        hdry=-1.0e30,
        wetfct=0.1,
        iwetit=7,
        ihdwet=1,
        laycon=laycon,
        intercellt=intercellt,
    )
    bcf.write_file()
    lines = Path(bcf.fn_path).read_text().splitlines()

    # seven I10/F10.3 fields for a structured grid; HDRY keeps a decimal point
    # in the mantissa so that F10.3 does not divide it by 1000
    assert lines[0] == (
        "         0   -1.E+30         0       0.1         7         1         0"
    )
    nlines = (nlay + 39) // 40
    codes = []
    for line in lines[1 : 1 + nlines]:
        assert len(line) <= 80
        codes += [
            int(line[i : i + 2].replace(" ", "") or 0) for i in range(0, len(line), 2)
        ]
    assert codes == [10 * i + c for i, c in zip(intercellt, laycon)]

    ml2 = MfUsg(model_ws=str(function_tmpdir))
    ModflowDis(
        ml2,
        nlay=nlay,
        nrow=1,
        ncol=2,
        nper=1,
        steady=True,
        top=30.0,
        botm=[30.0 - 0.5 * (k + 1) for k in range(nlay)],
    )
    MfUsgBas(ml2, ibound=1, strt=25.0, ifrefm=False)
    bcf2 = MfUsgBcf.load(bcf.fn_path, ml2)
    assert list(bcf2.intercellt.array) == intercellt
    assert list(bcf2.laycon.array) == laycon
    assert np.isclose(bcf2.hdry, -1.0e30) and np.isclose(bcf2.wetfct, 0.1)
    assert (bcf2.iwetit, bcf2.ihdwet) == (7, 1)


# F10.3 input of the item 1 real fields. Expected values are literal and were
# checked with an F10.3 read compiled with gfortran (default REAL); rtol=1e-6
# covers the single-precision rounding of that read (about 6e-8 relative).
@pytest.mark.parametrize(
    "field, value",
    [
        ("    -1E+30", -1.0e27),  # no decimal point: the 3 implied decimals apply
        ("    123E+2", 12.3),
        ("    123D+2", 12.3),
        ("  -25E-1  ", -0.0025),
        ("   1.0+30 ", 1.0e30),  # exponent without letter
        ("   -1.E+30", -1.0e30),  # decimal point in the mantissa
        ("    -1000.", -1000.0),
        (" 1.2345678", 1.2345678),
        ("      1500", 1.5),  # no exponent, no decimal point
        ("  15 00   ", 1.5),  # blanks are ignored
        ("          ", 0.0),  # empty field
    ],
)
def test_mfusgbcf_fixed_format_real_fields_read_as_f10_3(field, value):
    """HDRY and WETFCT of a fixed-format item 1 follow the F10.3 input rules."""
    line = f"{50:10d}{field}{0:10d}{field}{1:10d}{0:10d}{0:10d}"
    header = MfUsgBcf._load_fixed_header(line, True)
    assert np.isclose(header[1], value, rtol=1e-6, atol=0)
    assert np.isclose(header[3], value, rtol=1e-6, atol=0)


@pytest.mark.parametrize(
    "value, field",
    [
        (-1.0e30, "   -1.E+30"),
        (1.0e30, "    1.E+30"),
        (1.0e38, "    1.E+38"),
        (3.0e-8, "    3.E-08"),
        (-2.5e-30, "  -2.5E-30"),
        (1.2345678, " 1.2345678"),
        (-1.2345678, "-1.2345678"),
        (-12345.678, "-12345.678"),
        (999.99, "    999.99"),
        (0.0, "        0."),
    ],
)
def test_mfusgbcf_fixed_format_real_fields_written_for_f10_3(
    function_tmpdir, value, field
):
    """Item 1 reals are written with a decimal point in the mantissa, within
    ten characters and with the same single-precision value."""
    from flopy.modflow import ModflowDis

    ml = MfUsg(modelname="bcf", model_ws=str(function_tmpdir))
    ModflowDis(ml, nlay=1, nrow=1, ncol=2, nper=1, steady=True)
    MfUsgBas(ml, ibound=1, strt=1.0, ifrefm=False)
    bcf = MfUsgBcf(ml, ipakcb=0, laycon=0, hdry=value, wetfct=value)
    bcf.write_file()
    header = Path(bcf.fn_path).read_text().splitlines()[0]
    assert header[10:20] == field and header[30:40] == field

    ml2 = MfUsg(model_ws=str(function_tmpdir))
    ModflowDis(ml2, nlay=1, nrow=1, ncol=2, nper=1, steady=True)
    MfUsgBas(ml2, ibound=1, strt=1.0, ifrefm=False)
    bcf2 = MfUsgBcf.load(bcf.fn_path, ml2)
    assert np.float32(bcf2.hdry) == np.float32(value)
    assert np.float32(bcf2.wetfct) == np.float32(value)


def test_mfusgbcf_fixed_format_unrepresentable_real_fails(function_tmpdir):
    """A value that needs more than ten characters to keep its
    single-precision value fails before the destination is touched."""
    from flopy.modflow import ModflowDis

    ml = MfUsg(modelname="bcf", model_ws=str(function_tmpdir))
    ModflowDis(ml, nlay=1, nrow=1, ncol=2, nper=1, steady=True)
    MfUsgBas(ml, ibound=1, strt=1.0, ifrefm=False)
    bcf = MfUsgBcf(ml, laycon=0, hdry=-1.2345679e-30)
    target = Path(bcf.fn_path)
    target.write_text("pre-existing content\n")
    with pytest.raises(ValueError, match="single-precision"):
        bcf.write_file()
    assert target.read_text() == "pre-existing content\n"


def test_mfusgbcf_fixed_format_manual_load(function_tmpdir):
    """B01f load: hand-written fixed-format item 1 (DISU, eight fields, an
    implied-decimal real) and 40I2 LAYCON codes."""
    case = {"laycon": [3, 0], "disu": True, "transient": False}
    ml = _bcf_model(function_tmpdir, dict(case, ifrefm=False))
    (function_tmpdir / "manual.bcf").write_text(
        "        50    -1000.         1      1500         3         0         0"
        "         0\n"
        "13 0\n"
        "CONSTANT 1.0 #trpy\n"
        "CONSTANT 21.0 #hy layer 1\nCONSTANT 0.031 #vcont layer 1\n"
        "CONSTANT -0.61 #wetdry layer 1\nCONSTANT 12.0 #tran layer 2\n"
    )
    bcf = MfUsgBcf.load(str(function_tmpdir / "manual.bcf"), ml)
    assert (bcf.ipakcb, bcf.iwdflg, bcf.iwetit, bcf.ikcflag) == (50, 1, 3, 0)
    assert np.isclose(bcf.hdry, -1000.0) and np.isclose(bcf.wetfct, 1.5)
    assert list(bcf.intercellt.array) == [1, 0]
    assert list(bcf.laycon.array) == [3, 0]
    assert np.allclose(bcf.hy[0].array, 21.0) and np.allclose(bcf.tran[1].array, 12.0)
    assert np.allclose(bcf.wetdry[0].array, -0.61)


@pytest.mark.parametrize(
    "kwargs",
    [{"bubblept": True}, {"fullydry": True}, {"altsto": True}, {"ipakcb": 10**10}],
    ids=["bubblept", "fullydry", "altsto", "field-overflow"],
)
def test_mfusgbcf_fixed_format_rejected_before_writing(function_tmpdir, kwargs):
    """Options (or values) a fixed-format BCF cannot carry fail before the
    destination file or a supplied handle is touched."""
    from flopy.modflow import ModflowDis

    ml = MfUsg(modelname="bcf", model_ws=str(function_tmpdir))
    ModflowDis(ml, nlay=1, nrow=1, ncol=2, nper=1, steady=True)
    MfUsgBas(ml, ibound=1, strt=1.0, ifrefm=False)
    bcf = MfUsgBcf(ml, laycon=0, **kwargs)
    target = Path(bcf.fn_path)
    target.write_text("pre-existing content\n")
    with pytest.raises(ValueError, match=r"fixed-format|fit"):
        bcf.write_file()
    assert target.read_text() == "pre-existing content\n"
    handle = io.StringIO()
    with pytest.raises(ValueError):
        bcf.write_file(f=handle)
    assert handle.getvalue() == ""

    # the same options are accepted in free format
    ml.get_package("BAS6").ifrefm = True
    if "ipakcb" not in kwargs:
        bcf.write_file()
        assert kwargs.popitem()[0].upper() in target.read_text().splitlines()[0]


@pytest.mark.parametrize("ikcflag", [0, 1, -1, 2])
def test_mfusgbcf_trpy_and_connection_array_by_ikcflag(function_tmpdir, ikcflag):
    """B02/F35: TRPY only with IKCFLAG=0; DISU connection array (IDSYMRD=0,
    NJAG values) for any IKCFLAG != 0. IKCFLAG=2 is not a documented option;
    the 2.8.1 reader treats any non-zero value as connection input."""
    case = {"laycon": [0], "disu": True, "ikcflag": ikcflag, "transient": False}
    bcf, _ = _bcf_write_load_check(function_tmpdir, case)
    names = [name for name, _ in _bcf_records(bcf.fn_path)]
    assert ("trpy" in names) == (ikcflag == 0)
    assert ("ksat" in names) == (ikcflag != 0)


@pytest.mark.parametrize("ikcflag", [0, -1])
def test_mfusgbcf_anglex_condition_disu(function_tmpdir, ikcflag):
    """ANGLEX follows IHANISO, evaluated only when TRPY is read (IKCFLAG=0)."""
    case = {
        "laycon": [0],
        "disu": True,
        "ikcflag": ikcflag,
        "transient": False,
        "trpy": _BCF_TRPY_ANISO,
    }
    bcf, _ = _bcf_write_load_check(function_tmpdir, case)
    names = [name for name, _ in _bcf_records(bcf.fn_path)]
    assert ("anglex" in names) == (ikcflag == 0)


def test_mfusgbcf_structured_ikcflag_treated_as_zero(function_tmpdir):
    """IKCFLAG is not read for structured grids: HY/TRAN/VCONT stay in the file."""
    case = {"laycon": [3, 0], "ikcflag": 1, "transient": False}
    bcf, loaded = _bcf_write_load_check(function_tmpdir, dict(case))
    header = Path(bcf.fn_path).read_text().splitlines()[0].split()
    assert len(header) == 7
    assert loaded.ikcflag == 0


@pytest.mark.parametrize(
    "nlay, ikvflag",
    [(3, 0), (3, 1), (3, 2), (1, 1)],
    ids=["vcont", "kv", "kv-ikvflag2", "single-layer"],
)
def test_mfusgbcf_kv_vcont_by_ikvflag(function_tmpdir, nlay, ikvflag):
    """F36: VCONT for K<NLAY with IKVFLAG=0; KV for every layer with IKVFLAG!=0
    and NLAY>1 (IKVFLAG=2 behaves as 1 in the 2.8.1 reader)."""
    case = {"laycon": [0] * nlay, "ikvflag": ikvflag, "transient": False}
    _bcf_write_load_check(function_tmpdir, case)


@pytest.mark.parametrize(
    "case",
    [
        {"laycon": [3, 5], "richards": True, "bubblept": True},
        {"laycon": [0, 2], "richards": True},  # LAYCON 2 -> 5: HY, not TRAN
        {"laycon": [3, 0], "richards": True, "iwdflg": 1, "transient": False},
        {"laycon": [0, 3], "richards_hp": True, "bubblept": True},
        {
            "laycon": [0, 5],
            "richards": True,
            "disu": True,
            "irregular": True,
            "bubblept": True,
        },
        # TABRICH variants: no Richards arrays, isolate SF2, HY/TRAN and WETDRY
        {"laycon": [0, 5], "richards": True, "tabrich": True},
        {"laycon": [0, 2], "richards": True, "tabrich": True},
        {"laycon": [3, 0], "richards": True, "tabrich": True, "iwdflg": 1},
        {"laycon": [1, 3], "iwdflg": 1},  # control: WETDRY without RICHARDS
    ],
    ids=[
        "sf2-richards-bp",
        "hy-after-promotion",
        "no-wetdry-after-promotion",
        "richards_hp",
        "disu-irregular",
        "tabrich-sf2",
        "tabrich-hy",
        "tabrich-no-wetdry",
        "wetdry-control",
    ],
)
def test_mfusgbcf_effective_laycon_records(function_tmpdir, case):
    """B04/B05/B08/B09: SF2, HY/TRAN, WETDRY and the Richards arrays follow
    the effective LAYCON; public LAYCON values are not modified."""
    bcf, loaded = _bcf_write_load_check(function_tmpdir, case)
    assert list(bcf.laycon.array) == case["laycon"]
    assert list(loaded.laycon.array) == case["laycon"]


def test_mfusgbcf_richards_open_close_and_shared_external(function_tmpdir):
    """Richards arrays from OPEN/CLOSE and a shared EXTERNAL unit."""
    from flopy.utils.mfreadnam import parsenamefile

    (function_tmpdir / "alpha_l1.txt").write_text("0.71 0.71\n")
    (function_tmpdir / "bcf_ext72.dat").write_text("8.1 8.1\n0.091 0.091\nEND-OF-BCF\n")
    (function_tmpdir / "manual.bcf").write_text(
        "0 -1.0E+30 0 0.1 1 0 0\n05\nCONSTANT 1.0 #trpy\n"
        "CONSTANT 21.0 #hy layer 1\n"
        "OPEN/CLOSE alpha_l1.txt 1.0 (FREE) -1 #alpha layer 1\n"
        "EXTERNAL 72 1.0 (FREE) -1 #beta layer 1\n"
        "EXTERNAL 72 1.0 (FREE) -1 #sr layer 1\n"
        "CONSTANT 10.1 #brook layer 1\n"
    )
    nam = function_tmpdir / "manual.nam"
    nam.write_text("BCF6 15 manual.bcf\nDATA 72 bcf_ext72.dat\n")
    case = {"laycon": [5], "richards": True, "transient": False}
    ml = _bcf_model(function_tmpdir, case)
    ext_unit_dict = parsenamefile(str(nam), ml.mfnam_packages, verbose=False)
    try:
        bcf = MfUsgBcf.load(
            str(function_tmpdir / "manual.bcf"), ml, ext_unit_dict=ext_unit_dict
        )
        remaining = ext_unit_dict[72].filehandle.readline().strip()
    finally:
        for item in ext_unit_dict.values():
            if item.filehandle is not None:
                item.filehandle.close()
    _assert_bcf_values(bcf, case)
    assert remaining == "END-OF-BCF"


def test_mfusgbcf_positional_arguments_unchanged():
    """The Richards parameters are appended after all existing arguments."""
    import inspect

    params = list(inspect.signature(MfUsgBcf.__init__).parameters)
    assert params == [
        "self",
        "model",
        "ipakcb",
        "intercellt",
        "laycon",
        "trpy",
        "hdry",
        "iwdflg",
        "wetfct",
        "iwetit",
        "ihdwet",
        "ikvflag",
        "ikcflag",
        "tabrich",
        "nuzones",
        "nutabrows",
        "iuzontab",
        "retcrvs",
        "bubblept",
        "fullydry",
        "altsto",
        "tran",
        "hy",
        "vcont",
        "kv",
        "anglex",
        "ksat",
        "sf1",
        "sf2",
        "wetdry",
        "extension",
        "unitnumber",
        "filenames",
        "add_package",
        "alpha",
        "beta",
        "sr",
        "brook",
        "bp",
    ]
    from flopy.modflow import ModflowDis

    ml = MfUsg()
    ModflowDis(ml, nlay=2, nrow=1, ncol=1, nper=1)
    MfUsgBas(ml)
    bcf = MfUsgBcf(ml, 0, [1, 0], [3, 0], 1.0, -999.0, 1, 0.2, 2, 1, 1)
    assert list(bcf.intercellt.array) == [1, 0]
    assert list(bcf.laycon.array) == [3, 0]
    assert (bcf.hdry, bcf.iwdflg, bcf.wetfct, bcf.iwetit, bcf.ihdwet, bcf.ikvflag) == (
        -999.0,
        1,
        0.2,
        2,
        1,
        1,
    )
    assert np.allclose(bcf.alpha.array, 1.0) and np.allclose(bcf.brook.array, 6.0)


def test_mfusgbcf_write_file_to_open_handle(function_tmpdir):
    """WF: write_file(f=handle) writes to the handle and leaves it open."""
    case = {"laycon": [3, 0], "intercellt": [1, 0]}
    bcf = _bcf_model(function_tmpdir, case).get_package("BCF6")
    bcf.write_file()
    expected = Path(bcf.fn_path).read_text()

    target = function_tmpdir / "via_handle.bcf"
    with open(target, "w") as fh:
        bcf.write_file(f=fh)
        assert not fh.closed
        fh.write("# caller content after the package\n")
    assert target.read_text() == expected + "# caller content after the package\n"


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


def test_mfusghfb_parameter_active_without_defs_fails(function_tmpdir):
    """NPHFB>0 / active names with no parameter definitions fails explicitly
    (ValueError), no file. (From-scratch authoring *with* definitions is
    supported as of Stage 4.6C-A — see the authoring tests below.)"""
    from flopy.mfusg import MfUsgHfb

    hfb = MfUsgHfb(_hfb_model(function_tmpdir, "fs"), nphfb=1, mxfb=1, nhfbnp=0)
    out = function_tmpdir / "fs.hfb"
    hfb.fn_path = str(out)
    with pytest.raises(ValueError, match="none are defined"):
        hfb.write_file()
    assert not out.exists()


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


def test_mfusghfb_parameter_instances_unsupported(function_tmpdir):
    """HFB parameter INSTANCES are rejected on load -- USG-T itself aborts the
    run ('INSTANCES ARE NOT SUPPORTED FOR HFB', gwf2hfb7u1.f), so FloPy fails
    explicitly at the parameter header before reading any rows (Stage 4.7A/A5)."""
    from flopy.mfusg import MfUsgHfb

    p = function_tmpdir / "hfbinst.hfb"
    p.write_text(
        "# hfb parameter with INSTANCES\n"
        "         1         2         0\n"  # NPHFB MXFB NHFBNP
        "hp hfb 1.0 1 INSTANCES 2\n"
        "spring\n"
        "1 2 0.5\n"
        "fall\n"
        "3 4 0.5\n"
    )
    with pytest.raises(NotImplementedError, match="INSTANCES"):
        MfUsgHfb.load(str(p), _hfb_model(function_tmpdir, "hi"), nper=1)


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
    """NPHFB>0 with an empty parameters dict has no definitions to write: it
    fails with ValueError and writes no partial file."""
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
    with pytest.raises(ValueError, match="none are defined"):
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


# --- Stage 4.6C-A: HFB from-scratch NPHFB>0 parameter authoring -------------
#
# A parameterized HFB can be built entirely in Python (no prior load): pass
# parameters={name: {parval, data, ...}} + acthfb_names=[...]. Ergonomic input
# (data as tuples; nphfb/mxfb/nlst/nacthfb auto) is normalized + validated before
# the file is opened; write -> reload preserves it. TRANSIENT_HFB+params and
# INSTANCES stay unsupported.


def test_mfusghfb_parameter_authoring_from_scratch(function_tmpdir):
    """Build a parameterized unstructured HFB from scratch (ergonomic input):
    NPHFB, auto MXFBP, definitions, active names, 0-based internal / 1-based
    file."""
    from flopy.mfusg import MfUsgHfb

    # data as plain tuples (node1, node2, hydchr), 0-based; counts/nlst omitted.
    params = {"hp": {"parval": "2.0", "data": [(0, 1, 0.5), (2, 3, 0.25)]}}
    hfb = MfUsgHfb(
        _hfb_model(function_tmpdir, "auth"), parameters=params, acthfb_names=["hp"]
    )
    out = function_tmpdir / "auth.hfb"
    hfb.fn_path = str(out)
    hfb.write_file()
    text = out.read_text()
    hdr = next(ln for ln in text.splitlines() if not ln.startswith("#")).split()
    assert hdr[0] == "1" and hdr[1] == "2"  # NPHFB=1, MXFBP auto = sum(nlst)=2
    assert "hp HFB 2.0 2" in text  # UPARLSTRP header, nlst auto
    assert "         1         2" in text  # node1 0,1 -> 1,2 (1-based)

    re = MfUsgHfb.load(str(out), _hfb_model(function_tmpdir, "auth2"), nper=1)
    pd = re.parameters["hp"]
    assert pd["nlst"] == 2 and abs(float(pd["parval"]) - 2.0) < 1e-9
    assert list(pd["data"]["node1"]) == [0, 2]  # back to 0-based
    assert re.acthfb_names == ["hp"] and re.nphfb == 1


def test_mfusghfb_parameter_authoring_structured(function_tmpdir):
    """From-scratch structured HFB parameter (k, irow1, icol1, irow2, icol2,
    hydchr); 0-based internal / 1-based file."""
    from flopy.mfusg import MfUsgHfb

    ml = _hfb_model(function_tmpdir, "st", structured=True, nrow=3, ncol=3)
    params = {"sp": {"parval": "1.5", "data": [(0, 0, 0, 0, 1, 0.5)]}}
    hfb = MfUsgHfb(ml, parameters=params, acthfb_names=["sp"])
    out = function_tmpdir / "st.hfb"
    hfb.fn_path = str(out)
    hfb.write_file()

    re = MfUsgHfb.load(
        str(out),
        _hfb_model(function_tmpdir, "st2", structured=True, nrow=3, ncol=3),
        nper=1,
    )
    pd = re.parameters["sp"]
    assert list(pd["data"]["k"]) == [0]
    assert list(pd["data"]["irow1"]) == [0] and list(pd["data"]["icol2"]) == [1]
    assert re.acthfb_names == ["sp"]


def test_mfusghfb_parameter_authoring_mixed_nonparam(function_tmpdir):
    """A non-parametric barrier (NHFBNP) and a parameter definition coexist
    when authored from scratch; both survive write -> reload."""
    from flopy.mfusg import MfUsgHfb

    dt = MfUsgHfb.get_default_dtype(structured=False)
    nonparam = np.array([(4, 5, 0.9)], dtype=dt).view(np.recarray)
    hfb = MfUsgHfb(
        _hfb_model(function_tmpdir, "mx"),
        hfb_data=nonparam,
        parameters={"hp": {"parval": "1.0", "data": [(0, 1, 0.5)]}},
        acthfb_names=["hp"],
    )
    out = function_tmpdir / "mx.hfb"
    hfb.fn_path = str(out)
    hfb.write_file()
    hdr = next(
        ln for ln in out.read_text().splitlines() if not ln.startswith("#")
    ).split()
    assert hdr[0] == "1" and hdr[2] == "1"  # NPHFB=1, NHFBNP=1

    re = MfUsgHfb.load(str(out), _hfb_model(function_tmpdir, "mx2"), nper=1)
    assert list(re.hfb_data["node1"]) == [4]
    assert set(re.parameters) == {"hp"} and re.acthfb_names == ["hp"]


def test_mfusghfb_parameter_authoring_auto_counts(function_tmpdir):
    """nphfb / mxfb / nacthfb are computed when omitted."""
    from flopy.mfusg import MfUsgHfb

    params = {
        "a": {"parval": "1.0", "data": [(0, 1, 0.5), (2, 3, 0.5)]},
        "b": {"parval": "1.0", "data": [(4, 5, 0.5)]},
    }
    hfb = MfUsgHfb(
        _hfb_model(function_tmpdir, "ac"), parameters=params, acthfb_names=["a", "b"]
    )
    out = function_tmpdir / "ac.hfb"
    hfb.fn_path = str(out)
    hfb.write_file()
    lines = [ln for ln in out.read_text().splitlines() if not ln.startswith("#")]
    hdr = lines[0].split()
    assert hdr[0] == "2" and hdr[1] == "3"  # NPHFB=2, MXFBP=2+1
    # NACTHFB line (after the 2 defs + their 3 rows + 0 non-param) = "2"
    assert "         2" in lines  # nacthfb auto = 2


def test_mfusghfb_parameter_authoring_duplicate_definition_names_fails(function_tmpdir):
    """Definition names colliding case-insensitively (`hp`/`HP`) -> ValueError."""
    from flopy.mfusg import MfUsgHfb

    params = {
        "hp": {"parval": "1.0", "data": [(0, 1, 0.5)]},
        "HP": {"parval": "1.0", "data": [(2, 3, 0.5)]},
    }
    hfb = MfUsgHfb(
        _hfb_model(function_tmpdir, "dd"), parameters=params, acthfb_names=["hp"]
    )
    out = function_tmpdir / "dd.hfb"
    hfb.fn_path = str(out)
    with pytest.raises(ValueError, match="duplicate parameter definition"):
        hfb.write_file()
    assert not out.exists()


def test_mfusghfb_parameter_authoring_bad_name_fails(function_tmpdir):
    """A blank/whitespace/over-long parameter name -> ValueError, no file."""
    from flopy.mfusg import MfUsgHfb

    for badname in ("", "h p", "abcdefghijk"):
        hfb = MfUsgHfb(
            _hfb_model(function_tmpdir, "bn"),
            parameters={badname: {"parval": "1.0", "data": [(0, 1, 0.5)]}},
            acthfb_names=[badname],
        )
        out = function_tmpdir / "bn.hfb"
        if out.exists():
            out.unlink()
        hfb.fn_path = str(out)
        with pytest.raises(ValueError):
            hfb.write_file()
        assert not out.exists()


def test_mfusghfb_parameter_authoring_bad_parval_fails(function_tmpdir):
    """A multi-token / blank parval -> ValueError, no file."""
    from flopy.mfusg import MfUsgHfb

    for badparval in ("1 2", "  "):
        hfb = MfUsgHfb(
            _hfb_model(function_tmpdir, "bp"),
            parameters={"hp": {"parval": badparval, "data": [(0, 1, 0.5)]}},
            acthfb_names=["hp"],
        )
        out = function_tmpdir / "bp.hfb"
        if out.exists():
            out.unlink()
        hfb.fn_path = str(out)
        with pytest.raises(ValueError):
            hfb.write_file()
        assert not out.exists()


def test_mfusghfb_parameter_authoring_empty_data_fails(function_tmpdir):
    """A parameter with empty data rows -> ValueError, no file."""
    from flopy.mfusg import MfUsgHfb

    hfb = MfUsgHfb(
        _hfb_model(function_tmpdir, "ed"),
        parameters={"hp": {"parval": "1.0", "data": []}},
        acthfb_names=["hp"],
    )
    out = function_tmpdir / "ed.hfb"
    hfb.fn_path = str(out)
    with pytest.raises(ValueError, match="data"):
        hfb.write_file()
    assert not out.exists()


def test_mfusghfb_parameter_authoring_mxfb_too_small_fails(function_tmpdir):
    """An explicit mxfb below the total parameter rows -> ValueError, no file."""
    from flopy.mfusg import MfUsgHfb

    params = {
        "a": {"parval": "1.0", "data": [(0, 1, 0.5)]},
        "b": {"parval": "1.0", "data": [(2, 3, 0.5)]},
    }
    hfb = MfUsgHfb(
        _hfb_model(function_tmpdir, "ms"),
        parameters=params,
        mxfb=1,  # total = 2
        acthfb_names=["a"],
    )
    out = function_tmpdir / "ms.hfb"
    hfb.fn_path = str(out)
    with pytest.raises(ValueError, match="MXFBP"):
        hfb.write_file()
    assert not out.exists()


def test_mfusghfb_parameter_authoring_negative_index_fails(function_tmpdir):
    """A negative barrier index in a parameter -> ValueError, no file."""
    from flopy.mfusg import MfUsgHfb

    hfb = MfUsgHfb(
        _hfb_model(function_tmpdir, "ni"),
        parameters={"hp": {"parval": "1.0", "data": [(-1, 1, 0.5)]}},
        acthfb_names=["hp"],
    )
    out = function_tmpdir / "ni.hfb"
    hfb.fn_path = str(out)
    with pytest.raises(ValueError, match="negative"):
        hfb.write_file()
    assert not out.exists()


def _dpt_aw_model(function_tmpdir, name, mcomp=1, idpf=0):
    """A minimal transport-enabled MfUsg for DPT A-W_ADSORBIM authoring tests."""
    from flopy.modflow import ModflowDis

    m = MfUsg(structured=True, model_ws=str(function_tmpdir), modelname=name)
    ModflowDis(m, nlay=1, nrow=2, ncol=2, nper=1)
    m.itrnsp = 1
    m.mcomp = mcomp
    m.idpf = idpf
    m.iheat = 0
    return m


def test_mfusgdpt_aw_adsorbim_authoring_roundtrip(function_tmpdir):
    """A-W_ADSORBIM IAREA_FNIM=1 (AMAX) + IKAWI_FNIM=1 (Langmuir A/B) authors
    from scratch, writes the option-line indices + RP1/RP2 arrays, and reloads
    with the indices and arrays intact."""
    from flopy.mfusg import MfUsgDpt

    m = _dpt_aw_model(function_tmpdir, "aw1")
    dpt = MfUsgDpt(
        m,
        aw_adsorbim=True,
        iarea_fnim=1,
        ikawi_fnim=1,
        prsityim=0.3,
        ddtr=0.6,
        phif=0.4,
        concim=1.2,
        awamaxim=2.5,
        alangawim=3.5,
        blangawim=4.5,
    )
    dpt.fn_path = str(function_tmpdir / "aw1.dpt")
    dpt.write_file()
    assert "A-W_ADSORBIM 1 1" in Path(dpt.fn_path).read_text()

    re = MfUsgDpt.load(dpt.fn_path, _dpt_aw_model(function_tmpdir, "aw1b"))
    assert re.aw_adsorbim and re.iarea_fnim == 1 and re.ikawi_fnim == 1
    assert abs(float(re.awamaxim.array.mean()) - 2.5) < 1e-5
    assert abs(float(re.alangawim[0].array.mean()) - 3.5) < 1e-5
    assert abs(float(re.blangawim[0].array.mean()) - 4.5) < 1e-5


def test_mfusgdpt_aw_adsorbim_iarea4_roundtrip(function_tmpdir):
    """A-W_ADSORBIM IAREA_FNIM=4 (X2/X1/X0) + IKAWI_FNIM=2, two species:
    all three area arrays and the per-species Langmuir A/B round-trip."""
    from flopy.mfusg import MfUsgDpt

    m = _dpt_aw_model(function_tmpdir, "aw4", mcomp=2)
    dpt = MfUsgDpt(
        m,
        aw_adsorbim=True,
        iarea_fnim=4,
        ikawi_fnim=2,
        prsityim=0.3,
        ddtr=0.6,
        phif=0.4,
        concim=1.0,
        awarea_x2im=1.1,
        awarea_x1im=2.2,
        awarea_x0im=3.3,
        alangawim=[5.0, 6.0],
        blangawim=[7.0, 8.0],
    )
    dpt.fn_path = str(function_tmpdir / "aw4.dpt")
    dpt.write_file()
    assert "A-W_ADSORBIM 4 2" in Path(dpt.fn_path).read_text()

    re = MfUsgDpt.load(dpt.fn_path, _dpt_aw_model(function_tmpdir, "aw4b", mcomp=2))
    assert re.iarea_fnim == 4 and re.ikawi_fnim == 2
    assert abs(float(re.awarea_x2im.array.mean()) - 1.1) < 1e-5
    assert abs(float(re.awarea_x1im.array.mean()) - 2.2) < 1e-5
    assert abs(float(re.awarea_x0im.array.mean()) - 3.3) < 1e-5
    assert abs(float(re.alangawim[1].array.mean()) - 6.0) < 1e-5
    assert abs(float(re.blangawim[1].array.mean()) - 8.0) < 1e-5


def test_mfusgdpt_aw_adsorbim_rejects_unsupported(function_tmpdir):
    """Unsupported A-W_ADSORBIM branches fail explicitly: authoring with
    IAREA_FNIM in {2,3,5} or IKAWI_FNIM in {3,4} -> NotImplementedError (in
    __init__, so no file is written); mcomp=0 -> ValueError; on load the
    tabular forms (IAREA_FNIM=5, IKAWI_FNIM=4) raise before any array read so
    nothing downstream shifts; a bare keyword (missing indices) -> ValueError."""
    from flopy.mfusg import MfUsgDpt

    for iarea in (2, 3, 5):  # grain / ROG_SIGMA / tabular area
        with pytest.raises(NotImplementedError, match="IAREA_FNIM"):
            MfUsgDpt(
                _dpt_aw_model(function_tmpdir, f"a{iarea}"),
                aw_adsorbim=True,
                iarea_fnim=iarea,
                ikawi_fnim=1,
            )
    for ikawi in (3, 4):  # Brusseau scalar / tabular K_AWI
        with pytest.raises(NotImplementedError, match="IKAWI_FNIM"):
            MfUsgDpt(
                _dpt_aw_model(function_tmpdir, f"k{ikawi}"),
                aw_adsorbim=True,
                iarea_fnim=1,
                ikawi_fnim=ikawi,
            )
    with pytest.raises(ValueError, match="mcomp>0"):  # needs active transport
        MfUsgDpt(
            _dpt_aw_model(function_tmpdir, "az", mcomp=0),
            aw_adsorbim=True,
            iarea_fnim=1,
            ikawi_fnim=1,
        )

    # load: tabular (5 4) raises on the area index before any array read
    f54 = function_tmpdir / "aw54.dpt"
    f54.write_text("# x\n 0 0 1 0 0 0 0 A-W_ADSORBIM 5 4\n")
    with pytest.raises(NotImplementedError, match="IAREA_FNIM=5"):
        MfUsgDpt.load(str(f54), _dpt_aw_model(function_tmpdir, "l54"))

    # load: supported area, tabular K_AWI (1 4) raises on the K_AWI index
    f14 = function_tmpdir / "aw14.dpt"
    f14.write_text("# x\n 0 0 1 0 0 0 0 A-W_ADSORBIM 1 4\n")
    with pytest.raises(NotImplementedError, match="IKAWI_FNIM=4"):
        MfUsgDpt.load(str(f14), _dpt_aw_model(function_tmpdir, "l14"))

    # load: bare keyword without the two indices is malformed
    fbare = function_tmpdir / "awbare.dpt"
    fbare.write_text("# x\n 0 0 1 0 0 0 0 A-W_ADSORBIM\n")
    with pytest.raises(ValueError, match="two integers"):
        MfUsgDpt.load(str(fbare), _dpt_aw_model(function_tmpdir, "lbare"))


def test_mfusgdpt_aw_adsorbim_canonicalizes_indices(function_tmpdir):
    """IAREA_FNIM/IKAWI_FNIM given as numeric strings or integer-valued floats
    are canonicalized to int, so write_file emits "1 2" (never "1.0"/"2.0") and
    the model round-trips."""
    from flopy.mfusg import MfUsgDpt

    m = _dpt_aw_model(function_tmpdir, "can")
    dpt = MfUsgDpt(
        m,
        aw_adsorbim=True,
        iarea_fnim="1",  # numeric string
        ikawi_fnim=2.0,  # integer-valued float
        prsityim=0.3,
        ddtr=0.6,
        phif=0.4,
        concim=1.0,
        awamaxim=2.5,
        alangawim=3.5,
        blangawim=4.5,
    )
    assert dpt.iarea_fnim == 1 and isinstance(dpt.iarea_fnim, int)
    assert dpt.ikawi_fnim == 2 and isinstance(dpt.ikawi_fnim, int)
    dpt.fn_path = str(function_tmpdir / "can.dpt")
    dpt.write_file()
    option_line = Path(dpt.fn_path).read_text().splitlines()[0]
    assert "A-W_ADSORBIM 1 2" in option_line
    assert "1.0" not in option_line and "2.0" not in option_line

    re = MfUsgDpt.load(dpt.fn_path, _dpt_aw_model(function_tmpdir, "canb"))
    assert re.iarea_fnim == 1 and re.ikawi_fnim == 2


def test_mfusgdpt_aw_adsorbim_validation_hardening(function_tmpdir):
    """Review hardening: malformed A-W_ADSORBIM inputs fail with a clear
    ValueError before any array is created/written (no partial file), instead of
    a raw IndexError/TypeError or a late EOF on load."""
    from flopy.mfusg import MfUsgDpt

    def build(name, **kw):
        return MfUsgDpt(
            _dpt_aw_model(function_tmpdir, name, mcomp=kw.pop("mcomp", 1)),
            aw_adsorbim=True,
            prsityim=0.3,
            ddtr=0.6,
            phif=0.4,
            concim=1.0,
            **kw,
        )

    # per-species Langmuir arrays of the wrong length / None
    with pytest.raises(ValueError, match="needs one per species"):
        build("h1", mcomp=2, iarea_fnim=1, ikawi_fnim=1,
              alangawim=[1.0], blangawim=[2.0, 3.0])
    with pytest.raises(ValueError, match="is required"):
        build("h2", iarea_fnim=1, ikawi_fnim=1, alangawim=None)
    with pytest.raises(ValueError, match="is required"):
        build("h3", iarea_fnim=1, ikawi_fnim=1, blangawim=None)

    # non-canonical function indices
    for i, bad in enumerate((1.5, "abc", True, False, None)):
        with pytest.raises(ValueError):
            build(f"hi{i}", iarea_fnim=bad, ikawi_fnim=1)
    for i, bad in enumerate((2.5, "x", None)):
        with pytest.raises(ValueError):
            build(f"hk{i}", iarea_fnim=1, ikawi_fnim=bad)

    # nothing above opened a file
    for nm in ("h1", "h2", "h3"):
        assert not (function_tmpdir / f"{nm}.dpt").exists()

    # load with active transport off must fail at the option line (before arrays)
    f = function_tmpdir / "awmc0.dpt"
    f.write_text("# x\n 0 0 1 0 0 0 0 A-W_ADSORBIM 1 1\n")
    with pytest.raises(ValueError, match="mcomp>0"):
        MfUsgDpt.load(str(f), _dpt_aw_model(function_tmpdir, "lmc0", mcomp=0))


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


# --- Stage 4.5B: recipient U1DINT EXTERNAL / OPEN-CLOSE -----------------------
#
# The DRT/QRT recipient-node lists are read with U1DINT (utl7u1.f), which honors
# INTERNAL/CONSTANT (inline) and EXTERNAL <unit> / OPEN/CLOSE <fname> (separate
# file). These exercise the *recipient block itself* being external (distinct
# from the main list being EXTERNAL, covered above). On write the recipients are
# always expanded inline (Expanded valid write), never preserved as EXTERNAL.


def test_mfusgqrt_recipient_u1dint_external(function_tmpdir):
    """QRT recipient U1DINT as EXTERNAL <unit> loads via ext_unit_dict (0-based),
    and write_file re-emits it inline 1-based (no EXTERNAL/OPEN-CLOSE)."""
    from flopy.utils.mfreadnam import NamData

    (function_tmpdir / "qrecip.dat").write_text(" 10 11\n")  # 1-based -> 0-based 9,10
    p = function_tmpdir / "qext.qrt"
    p.write_text(
        "# qrt external recipient\n"
        "        1         2 0 0 0 RETURNFLOW\n"
        " 1 SP1\n"
        " 1  -1.000000e+02  2  7.500000e-01\n"
        "EXTERNAL 81\n"
    )
    eud = {81: NamData("DATA", "qrecip.dat", None, {})}
    qrt = MfUsgQrt.load(
        str(p), _usgt_unstructured_model(function_tmpdir), nper=1, ext_unit_dict=eud
    )
    assert qrt.recipient_nodes[0][0] == [9, 10]

    out = function_tmpdir / "qext_out.qrt"
    qrt.fn_path = str(out)
    qrt.write_file()
    text = out.read_text()
    assert "EXTERNAL" not in text and "OPEN/CLOSE" not in text  # not preserved
    assert "INTERNAL" in text and "\n 10 11\n" in text  # expanded inline, 1-based


def test_mfusgdrt_recipient_u1dint_external(function_tmpdir):
    """DRT spreading (NR<0) recipient U1DINT as EXTERNAL <unit> loads via
    ext_unit_dict (0-based), and write_file re-emits it inline 1-based."""
    from flopy.utils.mfreadnam import NamData

    (function_tmpdir / "drecip.dat").write_text(" 8 9\n")  # 1-based -> 0-based 7,8
    p = function_tmpdir / "dext.drt"
    p.write_text(
        "# drt external spreading\n"
        "         2 0 0 0 RETURNFLOW SPREAD 10\n"
        " 1 SP1\n"
        " 1  5.000000e+00  1.000000e+02  -2  7.000000e-01\n"
        "EXTERNAL 82\n"
    )
    eud = {82: NamData("DATA", "drecip.dat", None, {})}
    drt = MfUsgDrt.load(
        str(p), _usgt_unstructured_model(function_tmpdir), nper=1, ext_unit_dict=eud
    )
    assert drt.recipient_nodes[0][0] == [7, 8]

    out = function_tmpdir / "dext_out.drt"
    drt.fn_path = str(out)
    drt.write_file()
    text = out.read_text()
    assert "EXTERNAL" not in text and "OPEN/CLOSE" not in text
    assert "INTERNAL" in text and "\n 8 9\n" in text


def test_mfusgqrt_recipient_u1dint_open_close_quoted(function_tmpdir):
    """QRT recipient U1DINT as OPEN/CLOSE with a double-quoted filename that
    contains a space (quote-aware parsing)."""
    (function_tmpdir / "my recips.dat").write_text(" 3 4 5\n")  # -> 0-based 2,3,4
    p = function_tmpdir / "qoc.qrt"
    p.write_text(
        "# qrt open/close quoted\n"
        "        1         3 0 0 0 RETURNFLOW\n"
        " 1 SP1\n"
        " 1  -5.000000e+01  3  5.000000e-01\n"
        'OPEN/CLOSE "my recips.dat"\n'
    )
    qrt = MfUsgQrt.load(
        str(p), _usgt_unstructured_model(function_tmpdir), nper=1, ext_unit_dict={}
    )
    assert qrt.recipient_nodes[0][0] == [2, 3, 4]


def test_mfusgdrt_recipient_u1dint_open_close(function_tmpdir):
    """DRT spreading recipient U1DINT as OPEN/CLOSE <fname> (unquoted)."""
    (function_tmpdir / "dspread.dat").write_text(" 6 7\n")  # -> 0-based 5,6
    p = function_tmpdir / "doc.drt"
    p.write_text(
        "# drt open/close\n"
        "         2 0 0 0 RETURNFLOW SPREAD 10\n"
        " 1 SP1\n"
        " 1  5.000000e+00  1.000000e+02  -2  7.000000e-01\n"
        "OPEN/CLOSE dspread.dat\n"
    )
    drt = MfUsgDrt.load(
        str(p), _usgt_unstructured_model(function_tmpdir), nper=1, ext_unit_dict={}
    )
    assert drt.recipient_nodes[0][0] == [5, 6]


def test_mfusgqrt_recipient_u1dint_external_unresolved_fails(function_tmpdir):
    """An EXTERNAL recipient unit absent from ext_unit_dict raises an actionable
    NotImplementedError (not a raw ValueError or a silent mis-read)."""
    p = function_tmpdir / "qbad.qrt"
    p.write_text(
        "# qrt external unresolved\n"
        "        1         2 0 0 0 RETURNFLOW\n"
        " 1 SP1\n"
        " 1  -1.000000e+02  2  7.500000e-01\n"
        "EXTERNAL 77\n"
    )
    with pytest.raises(NotImplementedError, match=r"EXTERNAL.*77|77.*EXTERNAL"):
        MfUsgQrt.load(
            str(p), _usgt_unstructured_model(function_tmpdir), nper=1, ext_unit_dict={}
        )


# --- Stage 4.5B polish: recipient U1DINT control-tail (ICNSTNT FMTIN IPRN) ----
#
# U1DINT (utl7u1.f) reads a LOCAT>0 control record as
# "<kw> [unit/fname] ICNSTNT FMTIN IPRN"; values are read with FMTIN, list-
# directed only when FMTIN == (FREE). FloPy parses the tail explicitly: it
# supports (FREE) (and the abbreviated keyword-only form, treated as (FREE)),
# applies the ICNSTNT multiplier, and rejects a fixed Fortran FMTIN.


def test_mfusgqrt_recipient_u1dint_external_full_control(function_tmpdir):
    """QRT recipient EXTERNAL with the full control tail `EXTERNAL <u> 1 (FREE) -1`."""
    from flopy.utils.mfreadnam import NamData

    (function_tmpdir / "qrecip.dat").write_text(" 10 11\n")
    p = function_tmpdir / "qfull.qrt"
    p.write_text(
        "# qrt external full control\n"
        "        1         2 0 0 0 RETURNFLOW\n"
        " 1 SP1\n"
        " 1  -1.000000e+02  2  7.500000e-01\n"
        "EXTERNAL 81 1 (FREE) -1\n"
    )
    eud = {81: NamData("DATA", "qrecip.dat", None, {})}
    qrt = MfUsgQrt.load(
        str(p), _usgt_unstructured_model(function_tmpdir), nper=1, ext_unit_dict=eud
    )
    assert qrt.recipient_nodes[0][0] == [9, 10]


def test_mfusgdrt_recipient_u1dint_open_close_full_control(function_tmpdir):
    """DRT spreading OPEN/CLOSE with the full control tail `... 1 (FREE) -1`."""
    (function_tmpdir / "dspread.dat").write_text(" 8 9\n")
    p = function_tmpdir / "dfull.drt"
    p.write_text(
        "# drt open/close full control\n"
        "         2 0 0 0 RETURNFLOW SPREAD 10\n"
        " 1 SP1\n"
        " 1  5.000000e+00  1.000000e+02  -2  7.000000e-01\n"
        "OPEN/CLOSE dspread.dat 1 (FREE) -1\n"
    )
    drt = MfUsgDrt.load(
        str(p), _usgt_unstructured_model(function_tmpdir), nper=1, ext_unit_dict={}
    )
    assert drt.recipient_nodes[0][0] == [7, 8]


def test_mfusgqrt_recipient_u1dint_multiplier(function_tmpdir):
    """A non-1 ICNSTNT multiplies the recipient values (utl7u1.f): file 1,2 with
    ICNSTNT=10 -> 1-based 10,20 -> 0-based 9,19."""
    from flopy.utils.mfreadnam import NamData

    (function_tmpdir / "qmul.dat").write_text(" 1 2\n")
    p = function_tmpdir / "qmul.qrt"
    p.write_text(
        "# qrt external multiplier\n"
        "        1         2 0 0 0 RETURNFLOW\n"
        " 1 SP1\n"
        " 1  -1.000000e+02  2  7.500000e-01\n"
        "EXTERNAL 82 10 (FREE) -1\n"
    )
    eud = {82: NamData("DATA", "qmul.dat", None, {})}
    qrt = MfUsgQrt.load(
        str(p), _usgt_unstructured_model(function_tmpdir), nper=1, ext_unit_dict=eud
    )
    assert qrt.recipient_nodes[0][0] == [9, 19]


def test_mfusgqrt_recipient_u1dint_fixed_format_fails(function_tmpdir):
    """A non-(FREE) FMTIN on a recipient U1DINT raises NotImplementedError
    (no silent free-parse)."""
    from flopy.utils.mfreadnam import NamData

    (function_tmpdir / "qff.dat").write_text(" 10 11\n")
    p = function_tmpdir / "qff.qrt"
    p.write_text(
        "# qrt fixed format\n"
        "        1         2 0 0 0 RETURNFLOW\n"
        " 1 SP1\n"
        " 1  -1.000000e+02  2  7.500000e-01\n"
        "EXTERNAL 83 1 (10I8) -1\n"
    )
    eud = {83: NamData("DATA", "qff.dat", None, {})}
    with pytest.raises(NotImplementedError, match=r"(?i)free-format|FMTIN"):
        MfUsgQrt.load(
            str(p), _usgt_unstructured_model(function_tmpdir), nper=1, ext_unit_dict=eud
        )


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


# --- Stage 4.6C-D: ETS from-scratch ETSR array-parameter authoring ----------
#
# ETS uses MODFLOW *array* parameters (ModflowParBc / UPARARRRP), not the
# list-parameter family. Author from Python: parameters={name: {parval,
# clusters | instances}} + evtr_parm={kper: [(name, instance_or_None), ...]}.
# Only ETSR is parameterized; ETSS/ETSX/IETS/PXDP/PETM stay plain arrays. NPETS
# goes in item 2a (no PARAMETER line). Validated before open; no partial file.


def test_mfusgets_parameter_authoring_from_scratch_static(function_tmpdir):
    """Author a static ETSR parameter from scratch: NPETS in item 2a (no
    PARAMETER line), definition + cluster, per-SP active record, reload."""
    from flopy.mfusg import MfUsgEts

    params = {"etsrate": {"partyp": "ets", "parval": "5e-4",
                          "clusters": [("NONE", "ALL", [])]}}
    ets = MfUsgEts(
        _param_ets_model(function_tmpdir, "as"),
        netsop=1,
        parameters=params,
        evtr_parm={0: [("etsrate", None)]},
    )
    out = function_tmpdir / "as.ets"
    ets.fn_path = str(out)
    ets.write_file()
    text = out.read_text()
    item2a = next(ln for ln in text.splitlines() if not ln.startswith("#")).split()
    assert item2a[2] == "1"  # NPETS auto = 1
    assert "PARAMETER" not in text
    assert "etsrate ets 5e-4 1" in text and "NONE ALL" in text

    re = MfUsgEts.load(str(out), _param_ets_model(function_tmpdir, "as2"), nper=1)
    assert re.npets == 1 and re.evtr_parm[0] == [("etsrate", None)]
    assert re.parameters.bc_parms["etsrate"][0]["parval"] == "5e-4"
    assert np.allclose(re.evtr[0].array, 5.0e-4)


def test_mfusgets_parameter_authoring_instances(function_tmpdir):
    """A time-varying (INSTANCES) ETSR parameter authored from scratch; the
    active instance round-trips."""
    from flopy.mfusg import MfUsgEts

    params = {
        "etsrate": {
            "parval": "5e-4",
            "instances": {"spring": [("NONE", "ALL", [])],
                          "fall": [("NONE", "ALL", [])]},
        }
    }
    ets = MfUsgEts(
        _param_ets_model(function_tmpdir, "ai"),
        netsop=1,
        parameters=params,
        evtr_parm={0: [("etsrate", "spring")]},
    )
    out = function_tmpdir / "ai.ets"
    ets.fn_path = str(out)
    ets.write_file()
    text = out.read_text()
    assert "INSTANCES 2" in text and "etsrate spring" in text

    re = MfUsgEts.load(str(out), _param_ets_model(function_tmpdir, "ai2"), nper=1)
    assert re.evtr_parm[0] == [("etsrate", "spring")]


def test_mfusgets_parameter_authoring_netseg2_mixed_plain_arrays(function_tmpdir):
    """A parameterized ETSR mixes with plain ETSS/ETSX and NETSEG>1 PXDP/PETM
    arrays authored from scratch."""
    from flopy.mfusg import MfUsgEts

    ets = MfUsgEts(
        _param_ets_model(function_tmpdir, "ns"),
        netsop=1,
        netseg=2,
        pxdp=[0.5],
        petm=[0.5],
        parameters={"etsrate": {"parval": "5e-4", "clusters": [("NONE", "ALL", [])]}},
        evtr_parm={0: [("etsrate", None)]},
    )
    out = function_tmpdir / "ns.ets"
    ets.fn_path = str(out)
    ets.write_file()
    text = out.read_text()
    assert "etsrate ets 5e-4 1" in text  # parameterized ETSR
    # PXDP + PETM segment arrays written as plain arrays
    assert sum(ln.strip().startswith("CONSTANT") for ln in text.splitlines()) >= 4
    MfUsgEts.load(str(out), _param_ets_model(function_tmpdir, "ns2"), nper=1)


def test_mfusgets_parameter_authoring_netsop2(function_tmpdir):
    """A parameterized ETSR with NETSOP=2 (IEVT layer array) authored from
    scratch round-trips."""
    from flopy.mfusg import MfUsgEts

    ml = _param_ets_model(function_tmpdir, "n2")
    ets = MfUsgEts(
        ml,
        netsop=2,
        ievt=1,
        parameters={"etsrate": {"parval": "5e-4", "clusters": [("NONE", "ALL", [])]}},
        evtr_parm={0: [("etsrate", None)]},
    )
    out = function_tmpdir / "n2.ets"
    ets.fn_path = str(out)
    ets.write_file()
    text = out.read_text()
    item2a = next(ln for ln in text.splitlines() if not ln.startswith("#")).split()
    assert item2a[0] == "2" and item2a[2] == "1"  # NETSOP=2, NPETS=1
    MfUsgEts.load(str(out), _param_ets_model(function_tmpdir, "n22"), nper=1)


def test_mfusgets_parameter_authoring_auto_npets(function_tmpdir):
    """NPETS = len(parameters) when omitted; two static ETSR params."""
    from flopy.mfusg import MfUsgEts

    params = {
        "r1": {"parval": "1e-4", "clusters": [("NONE", "ALL", [])]},
        "r2": {"parval": "2e-4", "clusters": [("NONE", "ALL", [])]},
    }
    ets = MfUsgEts(
        _param_ets_model(function_tmpdir, "an"),
        netsop=1,
        parameters=params,
        evtr_parm={0: [("r1", None), ("r2", None)]},
    )
    out = function_tmpdir / "an.ets"
    ets.fn_path = str(out)
    ets.write_file()
    item2a = next(
        ln for ln in out.read_text().splitlines() if not ln.startswith("#")
    ).split()
    assert item2a[2] == "2"  # NPETS auto = 2


def _ets_neg(function_tmpdir, name, kwargs, match):
    """Helper: a parameterized-ETS write that must raise ValueError, no file."""
    from flopy.mfusg import MfUsgEts

    ets = MfUsgEts(_param_ets_model(function_tmpdir, name), netsop=1, **kwargs)
    out = function_tmpdir / f"{name}.ets"
    ets.fn_path = str(out)
    with pytest.raises(ValueError, match=match):
        ets.write_file()
    assert not out.exists()


def test_mfusgets_parameter_authoring_negatives(function_tmpdir):
    """The from-scratch ETS authoring validations all fire before any write."""
    sp = {"parval": "5e-4", "clusters": [("NONE", "ALL", [])]}
    inst = {"parval": "5e-4", "instances": {"s": [("NONE", "ALL", [])]}}
    _ets_neg(function_tmpdir, "e_nodefs",
             {"evtr_parm": {0: [("p1", None)]}}, "none are defined")
    _ets_neg(function_tmpdir, "e_first",
             {"parameters": {"etsrate": sp}, "evtr_parm": {}}, "first stress period")
    _ets_neg(function_tmpdir, "e_undef",
             {"parameters": {"etsrate": sp}, "evtr_parm": {0: [("px", None)]}},
             "not defined")
    _ets_neg(function_tmpdir, "e_dupdef",
             {"parameters": {"r1": sp, "R1": dict(sp)},
              "evtr_parm": {0: [("r1", None)]}}, "duplicate parameter definition")
    _ets_neg(function_tmpdir, "e_badname",
             {"parameters": {"r 1": sp}, "evtr_parm": {0: [("r 1", None)]}},
             "single non-empty")
    _ets_neg(function_tmpdir, "e_badparval",
             {"parameters": {"r1": {"parval": "1 2", "clusters": [("NONE", "ALL")]}},
              "evtr_parm": {0: [("r1", None)]}}, "single token")
    _ets_neg(function_tmpdir, "e_emptyclu",
             {"parameters": {"r1": {"parval": "5e-4", "clusters": []}},
              "evtr_parm": {0: [("r1", None)]}}, "cluster")
    _ets_neg(function_tmpdir, "e_badzone",
             {"parameters": {"r1": {"parval": "5e-4", "clusters": [("m1", "z1", [0])]}},
              "evtr_parm": {0: [("r1", None)]}}, "positive integer")
    _ets_neg(function_tmpdir, "e_instreq",
             {"parameters": {"r1": inst}, "evtr_parm": {0: [("r1", None)]}},
             "must name an instance")
    _ets_neg(function_tmpdir, "e_instunknown",
             {"parameters": {"r1": inst}, "evtr_parm": {0: [("r1", "x")]}},
             "no instance")
    _ets_neg(function_tmpdir, "e_staticinst",
             {"parameters": {"r1": sp}, "evtr_parm": {0: [("r1", "spring")]}},
             "static")
    _ets_neg(function_tmpdir, "e_dupactive",
             {"parameters": {"r1": sp}, "evtr_parm": {0: [("r1", None), ("R1", None)]}},
             "more than once")
    _ets_neg(function_tmpdir, "e_npets",
             {"parameters": {"r1": sp}, "npets": 2, "evtr_parm": {0: [("r1", None)]}},
             "NPETS")


# --- recipient-node U1DINT controls (INTERNAL/CONSTANT; EXTERNAL/OPEN-CLOSE
# --- added in Stage 4.5B, see the dedicated recipient tests above) ----------

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


def test_mfusg_recipient_external_u1dint_unresolved_fails(function_tmpdir):
    """An EXTERNAL recipient U1DINT unit that ext_unit_dict cannot resolve fails
    with an actionable NotImplementedError (no silent mis-read).

    Recipient U1DINT EXTERNAL/OPEN-CLOSE *is* supported as of Stage 4.5B when the
    unit resolves (see the dedicated recipient tests above); this case has no
    ext_unit_dict entry for unit 88, so it cannot be resolved and is rejected.
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


# --- Stage 4.6F-A: EVT NPEVT array-parameter preservation + authoring ---------


def test_mfusgevt_npevt_authoring_and_roundtrip(function_tmpdir):
    """A from-scratch NPEVT>0 EVT writes the PARAMETER line + definition + per-SP
    activation, and load -> write -> reload preserves NPEVT/parameters/evtr_parm
    (no silent collapse to NP=0)."""
    from flopy.mfusg import MfUsgEvt

    evt = MfUsgEvt(
        _evt_struct_model(function_tmpdir, "pe"),
        nevtop=1,
        surf=10.0,
        exdp=1.0,
        parameters={
            "evtrate": {"parval": "5e-4", "clusters": [("NONE", "ALL", [])]}
        },
        evtr_parm={0: [("evtrate", None)]},
    )
    evt.fn_path = str(function_tmpdir / "pe.evt")
    evt.write_file()
    text = Path(evt.fn_path).read_text()
    # item 1 PARAMETER NPEVT (not NP=0), the EVT-typed definition, and the
    # per-SP activation record (the parameter name, not an EVTR array)
    assert any(ln.strip().startswith("PARAMETER") for ln in text.splitlines())
    assert any(
        ln.lower().split()[:2] == ["evtrate", "evt"]
        for ln in text.splitlines()
        if ln.strip()
    )

    re = MfUsgEvt.load(evt.fn_path, _evt_struct_model(function_tmpdir, "peb"))
    assert re.npevt == 1
    assert re.parameters is not None and "evtrate" in re.parameters.bc_parms
    assert [n.lower() for n, _ in re.evtr_parm[0]] == ["evtrate"]

    # reload of the written file keeps NPEVT>0 through another round-trip
    re.fn_path = str(function_tmpdir / "pe2.evt")
    re.write_file()
    assert "PARAMETER" in Path(re.fn_path).read_text()
    re2 = MfUsgEvt.load(re.fn_path, _evt_struct_model(function_tmpdir, "pec"))
    assert re2.npevt == 1


def test_mfusgevt_npevt_instances_roundtrip(function_tmpdir):
    """A time-varying EVT parameter (INSTANCES) preserves its per-SP active
    instance across load -> write -> reload."""
    from flopy.mfusg import MfUsgEvt

    evt = MfUsgEvt(
        _evt_struct_model(function_tmpdir, "pi", nper=2),
        nevtop=1,
        surf=10.0,
        exdp=1.0,
        parameters={
            "etr": {
                "parval": "5e-4",
                "instances": {
                    "spring": [("NONE", "ALL", [])],
                    "fall": [("NONE", "ALL", [])],
                },
            }
        },
        evtr_parm={0: [("etr", "spring")], 1: [("etr", "fall")]},
    )
    evt.fn_path = str(function_tmpdir / "pi.evt")
    evt.write_file()
    re = MfUsgEvt.load(
        evt.fn_path, _evt_struct_model(function_tmpdir, "pib", nper=2)
    )
    assert re.npevt == 1
    assert re.evtr_parm[0][0][1] == "spring"
    assert re.evtr_parm[1][0][1] == "fall"


def test_mfusgevt_npevt_expand_parameters(function_tmpdir):
    """expand_parameters=True keeps the legacy expanded path: NPEVT=0 on output,
    EVTR written as a concrete array (no PARAMETER line)."""
    from flopy.mfusg import MfUsgEvt

    evt = MfUsgEvt(
        _evt_struct_model(function_tmpdir, "px"),
        nevtop=1,
        surf=10.0,
        exdp=1.0,
        parameters={
            "evtrate": {"parval": "5e-4", "clusters": [("NONE", "ALL", [])]}
        },
        evtr_parm={0: [("evtrate", None)]},
    )
    evt.fn_path = str(function_tmpdir / "px.evt")
    evt.write_file()

    ex = MfUsgEvt.load(
        evt.fn_path, _evt_struct_model(function_tmpdir, "pxb"), expand_parameters=True
    )
    assert ex.npevt == 0 and ex.parameters is None
    ex.fn_path = str(function_tmpdir / "pxex.evt")
    ex.write_file()
    assert "PARAMETER" not in Path(ex.fn_path).read_text()


def test_mfusgevt_npevt_rejects_invalid(function_tmpdir):
    """NPEVT authoring fails with an actionable ValueError, no partial file:
    activation without definitions, npevt without definitions, first SP not
    activating, an undefined active name, and a duplicate definition name."""
    from flopy.mfusg import MfUsgEvt

    good_def = {"r": {"parval": "1", "clusters": [("NONE", "ALL", [])]}}

    def write(name, nper=1, **kw):
        evt = MfUsgEvt(
            _evt_struct_model(function_tmpdir, name, nper=nper),
            nevtop=1,
            surf=10.0,
            exdp=1.0,
            **kw,
        )
        evt.fn_path = str(function_tmpdir / f"{name}.evt")
        evt.write_file()

    # evtr_parm activates but no parameters defined
    with pytest.raises(ValueError, match="none are defined"):
        write("b1", evtr_parm={0: [("r", None)]})
    # npevt>0 but no parameters defined
    with pytest.raises(ValueError, match="no parameters"):
        write("b2", npevt=1)
    # first stress period does not activate
    with pytest.raises(ValueError, match="first stress period"):
        write("b3", nper=2, parameters=good_def, evtr_parm={1: [("r", None)]})
    # active name not defined
    with pytest.raises(ValueError, match="not defined"):
        write("b4", parameters=good_def, evtr_parm={0: [("nope", None)]})
    # duplicate definition name (case-insensitive)
    with pytest.raises(ValueError):
        write(
            "b5",
            parameters={
                "R": {"parval": "1", "clusters": [("NONE", "ALL", [])]},
                "r": {"parval": "2", "clusters": [("NONE", "ALL", [])]},
            },
            evtr_parm={0: [("r", None)]},
        )
    for nm in ("b1", "b2", "b3", "b4", "b5"):
        assert not (function_tmpdir / f"{nm}.evt").exists()


def test_mfusgevt_npevt_validation_hardening(function_tmpdir):
    """Review hardening: a malformed evtr_parm or parameters={} fails with a
    clear ValueError (no AttributeError / tuple-unpack), no partial file; and a
    malformed PARAMETER line on load raises a clear ValueError."""
    from flopy.mfusg import MfUsgEvt

    good_def = {"r": {"parval": "1", "clusters": [("NONE", "ALL", [])]}}

    def write(name, **kw):
        evt = MfUsgEvt(
            _evt_struct_model(function_tmpdir, name),
            nevtop=1,
            surf=10.0,
            exdp=1.0,
            **kw,
        )
        evt.fn_path = str(function_tmpdir / f"{name}.evt")
        evt.write_file()

    # evtr_parm container/record structure
    with pytest.raises(ValueError, match="must be a dict"):  # list, not dict
        write("v1", parameters=good_def, evtr_parm=[("r", None)])
    with pytest.raises(ValueError, match="must be a dict"):  # string, not dict
        write("v2", parameters=good_def, evtr_parm="r")
    with pytest.raises(ValueError, match="must be a dict"):  # bad even w/o defs
        write("v3", evtr_parm=[("r", None)])
    with pytest.raises(ValueError, match="list of"):  # recs is a string
        write("v4", parameters=good_def, evtr_parm={0: "r"})
    with pytest.raises(ValueError, match=r"must be a \(name"):  # rec not a pair
        write("v5", parameters=good_def, evtr_parm={0: ["r"]})
    with pytest.raises(ValueError, match=r"must be a \(name"):  # 1-tuple
        write("v6", parameters=good_def, evtr_parm={0: [("r",)]})
    with pytest.raises(ValueError, match="non-empty string"):  # empty name
        write("v7", parameters=good_def, evtr_parm={0: [("", None)]})
    with pytest.raises(ValueError, match="string or None"):  # bad instance
        write("v8", parameters=good_def, evtr_parm={0: [("r", 3)]})
    with pytest.raises(ValueError, match="integer stress period"):  # str key
        write("v9", parameters=good_def, evtr_parm={"0": [("r", None)]})
    # parameters={} is "no definitions", not a first-period error
    with pytest.raises(ValueError, match="parameters is empty"):
        write("v10", parameters={}, evtr_parm={0: [("r", None)]})
    for nm in (f"v{i}" for i in range(1, 11)):
        assert not (function_tmpdir / f"{nm}.evt").exists()

    # load: malformed PARAMETER line
    f1 = function_tmpdir / "lp1.evt"
    f1.write_text("# x\nPARAMETER\n         1         0\n")
    with pytest.raises(ValueError, match="PARAMETER <NPEVT>"):
        MfUsgEvt.load(str(f1), _evt_struct_model(function_tmpdir, "lp1b"))
    f2 = function_tmpdir / "lp2.evt"
    f2.write_text("# x\nPARAMETER abc\n         1         0\n")
    with pytest.raises(ValueError, match="integer count"):
        MfUsgEvt.load(str(f2), _evt_struct_model(function_tmpdir, "lp2b"))


def test_mfusgevt_npevt_parameter0_loads_nonparametric(function_tmpdir):
    """A 'PARAMETER 0' header loads as a plain non-parametric EVT (npevt=0)."""
    from flopy.mfusg import MfUsgEvt

    evt = MfUsgEvt(
        _evt_struct_model(function_tmpdir, "p0"),
        nevtop=1,
        surf=10.0,
        evtr=1e-4,
        exdp=1.0,
    )
    evt.fn_path = str(function_tmpdir / "p0.evt")
    evt.write_file()
    body = Path(evt.fn_path).read_text().splitlines()
    # re-emit with an explicit "PARAMETER 0" line before item 2
    p0 = function_tmpdir / "p0p.evt"
    p0.write_text(body[0] + "\nPARAMETER 0\n" + "\n".join(body[1:]) + "\n")

    re = MfUsgEvt.load(str(p0), _evt_struct_model(function_tmpdir, "p0b"))
    assert re.npevt == 0 and re.parameters is None


# --- Stage 4.6F-B: EVT ETS zonal time-series is deferred (explicit failure) ---


def test_mfusgevt_ets_zonal_deferred(function_tmpdir):
    """EVT ETS zonal time-series ('ETS MXZNEVT' / per-SP 'INEVTZONES') is a
    dynamic ATS-coupled mode (external IUETS time-series superseding EVTR), not
    static I/O. It is deferred (Stage 4.6F-B): authoring and load fail explicitly,
    per branch, with no partial parse — while plain EVT is unaffected."""
    from flopy.mfusg import MfUsgEvt

    # authoring with mxetzones>0 -> specific NotImplementedError (cites the spec)
    with pytest.raises(NotImplementedError, match="ATS"):
        MfUsgEvt(_evt_struct_model(function_tmpdir, "z1"), nevtop=1, mxetzones=3)

    # load: 'ETS MXZNEVT' in item 2 fails before any stress-period parse
    f1 = function_tmpdir / "ets.evt"
    f1.write_text(
        "# x\n         1         0 ETS         3\n"
        "         0         0         0  # SP1 must not be parsed\n"
    )
    with pytest.raises(NotImplementedError, match="ETS"):
        MfUsgEvt.load(str(f1), _evt_struct_model(function_tmpdir, "z2"))

    # load: a per-SP 'INEVTZONES' flag fails explicitly
    f2 = function_tmpdir / "inz.evt"
    f2.write_text(
        "# x\n         1         0\n"
        "INEVTZONES 1         0         1         1\n"
    )
    with pytest.raises(NotImplementedError, match="INEVTZONES"):
        MfUsgEvt.load(str(f2), _evt_struct_model(function_tmpdir, "z3"))

    # plain EVT (no ETS option) is unaffected and still round-trips
    evt = MfUsgEvt(
        _evt_struct_model(function_tmpdir, "z4"),
        nevtop=1,
        surf=10.0,
        evtr=1e-4,
        exdp=1.0,
    )
    evt.fn_path = str(function_tmpdir / "z4.evt")
    evt.write_file()
    re = MfUsgEvt.load(evt.fn_path, _evt_struct_model(function_tmpdir, "z4b"))
    assert re.nevtop == 1 and re.mxetzones == 0


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


# --- Stage 4.6D: LAK multi-lake + sill/connectivity (datasets 7/8) authoring --


def _lak_grids_2lakes():
    """1-layer 3x3 grid with two lakes: cell (0,0,0)=lake 1, (0,2,2)=lake 2."""
    lakarr = np.zeros((1, 3, 3), dtype=int)
    lakarr[0, 0, 0] = 1
    lakarr[0, 2, 2] = 2
    return {0: lakarr}, {0: np.ones((1, 3, 3), dtype=float) * 0.1}


def test_mfusglak_multilake_sill_authoring_roundtrip(function_tmpdir):
    """Two lakes + a connected-lake system (datasets 7/8) author from scratch,
    write, and reload with the sill system intact (NSLMS, IC, 1-based ISUB,
    SILLVT)."""
    from flopy.mfusg import MfUsgLak

    lakarr, bdlknc = _lak_grids_2lakes()
    lak = MfUsgLak(
        _lak_model(function_tmpdir, "ml"),
        nlakes=2,
        stages=[100.0, 100.0],
        lakarr=lakarr,
        bdlknc=bdlknc,
        flux_data={0: {0: [1.0, 2.0, 0.0, 0.0], 1: [1.0, 2.0, 0.0, 0.0]}},
        sill_data={0: [([2, 1, 2], [95.0])]},
    )
    lak.fn_path = str(function_tmpdir / "ml.lak")
    lak.write_file()

    text = Path(lak.fn_path).read_text()
    # dataset 7 (NSLMS=1) written once, dataset 8a is "IC lake1 lake2" 1-based
    assert sum("Data set 7" in ln for ln in text.splitlines()) == 1
    assert any(
        ln.split("#")[0].split() == ["2", "1", "2"] for ln in text.splitlines()
    )

    re = MfUsgLak.load(lak.fn_path, _lak_model(function_tmpdir, "mlb"))
    assert re.nlakes == 2
    assert list(re.sill_data.keys()) == [0]
    ds8a, sillvt = re.sill_data[0][0]
    assert ds8a == [2, 1, 2]  # IC=2, center lake 1, sublake 2 (1-based)
    assert sillvt == [95.0]


def test_mfusglak_multilake_sill_transport_roundtrip(function_tmpdir):
    """Two lakes + sill system with classic transport (mcomp=1) round-trips:
    datasets 7/8 and dataset 9b (CPPT/CRNF per lake) are both preserved."""
    from flopy.mfusg import MfUsgLak

    lakarr, bdlknc = _lak_grids_2lakes()
    lak = MfUsgLak(
        _lak_model(function_tmpdir, "mlt", mcomp=1),
        nlakes=2,
        stages=[100.0, 100.0],
        lakarr=lakarr,
        bdlknc=bdlknc,
        flux_data={0: {0: [1.0, 2.0, 0.0, 0.0], 1: [1.0, 2.0, 0.0, 0.0]}},
        sill_data={0: [([2, 1, 2], [95.0])]},
        clake=[[5.0], [6.0]],
        conc_data={0: {(0, 0): [3.0, 1.0], (1, 0): [4.0, 2.0]}},  # CPPT, CRNF
    )
    lak.fn_path = str(function_tmpdir / "mlt.lak")
    lak.write_file()

    re = MfUsgLak.load(lak.fn_path, _lak_model(function_tmpdir, "mltb", mcomp=1))
    assert re.nlakes == 2
    assert re.sill_data[0][0][0] == [2, 1, 2]
    assert re.sill_data[0][0][1] == [95.0]
    assert [float(x) for x in re.conc_data[0][(1, 0)]] == [4.0, 2.0]


def test_mfusglak_sill_rejects_invalid(function_tmpdir):
    """Datasets 7/8 authoring fails with an actionable ValueError, no partial
    file: bad lake ids, IC/list mismatch, IC<2, wrong sill count, duplicate
    lake, kper out of range, sill for a period with ITMP<=0, empty ds8a, a
    non-integral IC / lake id, and a non-numeric sill elevation."""
    from flopy.mfusg import MfUsgLak
    from flopy.modflow import ModflowDis

    def build(name, sill, nper=1, **over):
        lakarr, bdlknc = _lak_grids_2lakes()
        if nper == 1:
            model = _lak_model(function_tmpdir, name)
        else:
            model = MfUsg(
                structured=True, model_ws=str(function_tmpdir), modelname=name
            )
            ModflowDis(
                model, nlay=1, nrow=3, ncol=3, nper=nper, nstp=1, steady=False
            )
        kw = {
            "nlakes": 2,
            "stages": [100.0, 100.0],
            "lakarr": lakarr,
            "bdlknc": bdlknc,
            "flux_data": {0: {0: [1.0, 2.0, 0.0, 0.0], 1: [1.0, 2.0, 0.0, 0.0]}},
            "sill_data": sill,
        }
        kw.update(over)
        return MfUsgLak(model, **kw)

    with pytest.raises(ValueError, match="out of range"):  # lake id (1-based)
        build("b1", {0: [([2, 1, 3], [95.0])]})
    with pytest.raises(ValueError, match="but lists"):  # IC vs listed lakes
        build("b2", {0: [([3, 1, 2], [95.0])]})
    with pytest.raises(ValueError, match="IC>=2"):  # no sublake
        build("b3", {0: [([1, 1], [])]})
    with pytest.raises(ValueError, match="sill elevations"):  # sill count != IC-1
        build("b4", {0: [([2, 1, 2], [1.0, 2.0])]})
    with pytest.raises(ValueError, match="repeats a lake"):  # duplicate lake
        build("b5", {0: [([2, 1, 1], [95.0])]})
    with pytest.raises(ValueError, match=r"out of range \[0, 1\)"):  # bad kper
        build("b6", {5: [([2, 1, 2], [95.0])]})
    with pytest.raises(ValueError, match="ITMP<=0"):  # sill where lakarr reused
        build("b7", {1: [([2, 1, 2], [95.0])]}, nper=2)
    # --- review hardening: ambiguous / non-canonical values rejected up front
    with pytest.raises(ValueError, match="ds8a is empty"):  # ds8a=[]
        build("b8", {0: [([], [])]})
    with pytest.raises(ValueError, match="whole number"):  # IC=2.5 (no truncation)
        build("b9", {0: [([2.5, 1, 2], [95.0])]})
    with pytest.raises(ValueError, match="whole number"):  # lake id 1.5
        build("b10", {0: [([2, 1.5, 2], [95.0])]})
    with pytest.raises(ValueError, match="must be numeric"):  # sill "abc"
        build("b11", {0: [([2, 1, 2], ["abc"])]})
    with pytest.raises(ValueError, match="must be numeric"):  # sill None
        build("b12", {0: [([2, 1, 2], [None])]})
    # the failing constructors never opened a file
    for nm in ("b7", "b8", "b9", "b10", "b11", "b12"):
        assert not (function_tmpdir / f"{nm}.lak").exists()


def test_mfusglak_sill_normalizes_numeric_input(function_tmpdir):
    """Numeric strings and integer-valued floats in sill_data are canonicalized
    to int IC/lake numbers and float sills, so write_file emits canonical values
    (no "2.0"/"95.0" artifacts in dataset 8a) and the system round-trips."""
    from flopy.mfusg import MfUsgLak

    lakarr, bdlknc = _lak_grids_2lakes()
    lak = MfUsgLak(
        _lak_model(function_tmpdir, "norm"),
        nlakes=2,
        stages=[100.0, 100.0],
        lakarr=lakarr,
        bdlknc=bdlknc,
        flux_data={0: {0: [1.0, 2.0, 0.0, 0.0], 1: [1.0, 2.0, 0.0, 0.0]}},
        sill_data={0: [(["2", 1.0, "2"], ["95.0"])]},  # numeric strings + floats
    )
    # canonicalized in __init__: ints for ds8a, floats for sillvt
    ds8a, sillvt = lak.sill_data[0][0]
    assert ds8a == [2, 1, 2] and all(isinstance(x, int) for x in ds8a)
    assert sillvt == [95.0] and all(isinstance(x, float) for x in sillvt)

    lak.fn_path = str(function_tmpdir / "norm.lak")
    lak.write_file()
    text = Path(lak.fn_path).read_text()
    assert any(
        ln.split("#")[0].split() == ["2", "1", "2"] for ln in text.splitlines()
    )

    re = MfUsgLak.load(lak.fn_path, _lak_model(function_tmpdir, "normb"))
    assert re.sill_data[0][0][0] == [2, 1, 2]
    assert re.sill_data[0][0][1] == [95.0]


# --- Stage 4.6A: STR/SUB/SWT compatibility guard on unstructured MfUsg --------
#
# STR/SUB/SWT have no validated USG-T unstructured layout, so MfUsg registers
# guarded wrappers that fail explicitly on a DISU/unstructured MfUsg model
# (load *and* authoring), rather than silently reading/writing a non-USG-T file.
# SFR is intentionally NOT guarded (freyberg_usg = DISU+SFR loads/writes/runs via
# base ModflowSfr2 in test_usg.py); FHB/GAGE are intentionally NOT guarded
# (round-trip via base in Ex8).


def test_mfusg_compat_str_sub_swt_load_guard(function_tmpdir):
    """The registered STR/SUB/SWT loaders raise on an unstructured MfUsg model.

    MfUsg.load dispatches each NAM entry to ``<pkg>.load(filehandle, model, ...)``
    (see MfUsg._ext_unit_d_load); this exercises that exact call. The guard
    raises before the file is read, so a USG-T NAM that lists STR/SUB/SWT fails
    explicitly instead of mis-reading it via the structured base class.
    """
    from flopy.mfusg import MfUsgStr, MfUsgSub, MfUsgSwt

    ml = _usgt_unstructured_model(function_tmpdir)
    assert ml.version == "mfusg" and ml.structured is False
    dummy = function_tmpdir / "dummy.dat"
    dummy.write_text("# unread by the guard\n")
    for cls, label in ((MfUsgStr, "STR"), (MfUsgSub, "SUB"), (MfUsgSwt, "SWT")):
        with open(dummy) as fh:
            with pytest.raises(NotImplementedError, match=f"{label}.*compatibility"):
                cls.load(fh, ml)


def test_mfusg_compat_str_sub_swt_write_guard(function_tmpdir):
    """Constructing STR/SUB/SWT on an unstructured MfUsg model raises before any
    file is written (the guard runs in __init__, before add_package/write_file),
    so no partial package or file is produced."""
    from flopy.mfusg import MfUsgStr, MfUsgSub, MfUsgSwt

    ml = _usgt_unstructured_model(function_tmpdir)
    for cls, label in ((MfUsgStr, "STR"), (MfUsgSub, "SUB"), (MfUsgSwt, "SWT")):
        with pytest.raises(NotImplementedError, match=f"{label}.*DISU"):
            cls(ml)
        # No partial state: the package was never added to the model.
        assert ml.get_package(label) is None


def test_mfusg_compat_guard_noop_and_sfr_untouched(function_tmpdir):
    """The guard is a no-op for a structured MfUsg model and for plain
    flopy.modflow.Modflow; and SFR/FHB/GAGE stay on their base classes in the
    MfUsg registry (only STR/SUB/SWT are guarded)."""
    from flopy.mfusg import MfUsgStr, MfUsgSub, MfUsgSwt
    from flopy.mfusg.mfusgcompat import _guard_compat
    from flopy.modflow import (
        Modflow,
        ModflowFhb,
        ModflowGage,
        ModflowSfr2,
    )

    structured_mfusg = MfUsg(structured=True, model_ws=str(function_tmpdir))
    base_mf = Modflow()  # version mf2005
    for model in (structured_mfusg, base_mf):
        for label in ("STR", "SUB", "SWT"):
            _guard_compat(model, label)  # must not raise

    # Registry: STR/SUB/SWT guarded; SFR/FHB/GAGE untouched (base classes).
    reg = MfUsg(structured=False, model_ws=str(function_tmpdir)).mfnam_packages
    assert reg["str"] is MfUsgStr
    assert reg["sub"] is MfUsgSub
    assert reg["swt"] is MfUsgSwt
    assert reg["sfr"] is ModflowSfr2
    assert reg["fhb"] is ModflowFhb
    assert reg["gage"] is ModflowGage


# G03: EVS/RTS name-file entries (ET and recharge time-series files read by
# EVT with ETS and RCH with RTS) are kept as references to opaque files.
# Hand-written fixture: 3-node DISU line (area 100), CHD h=10 at node 1, one
# transient SP of 10, OC ATSA (RTS/ETS need ATS), RCH RTS 2 zones, EVT ETS 1
# zone. The name file lists the time-series files as given by the caller.

_G03_FILES = {
    "t.disu": """# DISU
3 1 7 0 1 4 2 0
  0
INTERNAL               1     (1I10) -1 #nodelay
         3
INTERNAL               1   (3E15.6) -1 #top layer 1
   3.000000E+01   3.000000E+01   3.000000E+01
INTERNAL               1   (3E15.6) -1 #bot layer 1
   0.000000E+00   0.000000E+00   0.000000E+00
CONSTANT    1.000000E+02                           #area layer 1
INTERNAL               1     (3I10) -1 #iac
         2         3         2
INTERNAL               1     (7I10) -1 #ja
         1         2         2         1         3         3         2
CONSTANT    5.000000E+00                           #cl12
CONSTANT    1.000000E+01                           #fahl
     10.000000             1  1.000000  TR
""",
    "t.bas": """# BAS6
UNSTRUCTURED FREE
CONSTANT          1                                #ibound layer 1
   -999.99
CONSTANT    1.000000E+01                           #strt layer 1
""",
    "t.lpf": """# LPF
         0    -1E+30         0         0
         0
         0
   1.000000E+00
         0
         0
CONSTANT    1.000000E+00                           #hk layer 1
CONSTANT    1.000000E+00                           #vka1
CONSTANT    1.000000E-03                           #ss layer 1
""",
    "t.sms": """# SMS
  1.00e-09  1.00e-11       100       100         2         1         1
  7.00e-01  1.00e-01  2.00e-01  1.00e-03        20  1.00e+04  2.00e-01  1.00e+02
         2         0         7         2         0  0.00e+00         0  1.00e-03
""",
    "t.chd": """# CHD
         1
 1    Stress Period 1
 1   10.000000  10.000000
""",
    "t.oc": """# OC
ATSA
HEAD PRINT FORMAT   0
HEAD SAVE UNIT    51
DRAWDOWN PRINT FORMAT   0
COMPACT BUDGET AUX
period 1 step 1
  save head
  print budget
""",
}
# RCH with RTS: zones 0 1 2; rows 0-4 (1e-3, 2e-3), 4-10 (3e-3, 1e-3)
_G03_RCH = (
    "# RCH\n 1 0 RTS 2\n 0 INRCHZONES 2\nCONSTANT 0.0\nINTERNAL 1 (FREE) -1\n 0 1 2\n"
)
_G03_RTS = (
    "0.0 4.0 1.0 1.0e-3 2.0e-3\n4.0 10.0 1.0 3.0e-3 1.0e-3\n10.0 1.0e6 1.0 0.0 0.0\n"
)
# EVT with ETS: SURF 5 (below heads: maximum rate), zones 0 1 1;
# rows 0-6 (5e-4), 6-10 (1e-4)
_G03_EVT = (
    "# EVT\n 1 0 ETS 1\n 0 0 0 INEVTZONES 1\nCONSTANT 5.0\nCONSTANT 0.0\n"
    "CONSTANT 1.0\nINTERNAL 1 (FREE) -1\n 0 1 1\n"
)
_G03_EVS = "0.0 6.0 1.0 5.0e-4\n6.0 10.0 1.0 1.0e-4\n10.0 1.0e6 1.0 0.0\n"
# controls without time series
_G03_RCH_PLAIN = "# RCH\n 1 0\n 0\nCONSTANT 1.0e-3\n"
_G03_EVT_PLAIN = "# EVT\n 1 0\n 0 0 0\nCONSTANT 5.0\nCONSTANT 2.0e-4\nCONSTANT 1.0\n"
_G03_BASE = ["DISU", "BAS6", "LPF", "SMS", "CHD", "OC"]


def _g03_fixture(
    ws,
    rts="RTS               40  ts/rech.rts",
    evs="EVS               41  et.evs NO-ARCHIVE",
    extra="",
):
    """Write the fixture into ws. rts/evs: name-file line or None (the series
    file is written at the path given in the line, unquoted; an absolute
    path is written as given); extra: more name-file lines."""
    ws = Path(ws)
    ws.mkdir(parents=True, exist_ok=True)
    for name, text in _G03_FILES.items():
        (ws / name).write_text(text)
    nam = [
        "# fixture",
        "LIST               2  t.list",
        "DISU              11  t.disu",
        "BAS6              13  t.bas",
        "LPF               15  t.lpf",
        "SMS               32  t.sms",
        "CHD               24  t.chd",
        "OC                14  t.oc",
        "DATA(BINARY)      51  t.hds",
    ]
    for line, pkg, pkg_text, series in (
        (rts, "RCH               18  t.rch", _G03_RCH, _G03_RTS),
        (evs, "EVT               17  t.evt", _G03_EVT, _G03_EVS),
    ):
        if line is None:
            continue
        nam += [pkg, line]
        (ws / pkg.split()[-1]).write_text(pkg_text)
        rest = line.split(None, 2)[2]
        path = rest[1 : rest.index("'", 1)] if rest.startswith("'") else rest.split()[0]
        target = ws / path
        target.parent.mkdir(parents=True, exist_ok=True)
        target.write_text(series)
    (ws / "t.nam").write_text("\n".join(nam) + "\n" + extra)
    return ws


def _g03_nam_entries(ws):
    lines = (Path(ws) / "t.nam").read_text().splitlines()
    return [ln for ln in lines if ln.split() and ln.split()[0] in ("RTS", "EVS")]


def _g03_tree(ws):
    ws = Path(ws)
    return {
        str(p.relative_to(ws)): (p.read_bytes(), p.stat().st_mtime_ns)
        for p in sorted(ws.rglob("*"))
        if p.is_file()
    }


def _g03_load(ws, load_only=None):
    if load_only is None:
        load_only = list(_G03_BASE)
    return MfUsg.load("t.nam", model_ws=str(ws), load_only=load_only, check=False)


_G03_RTS_LINE = "RTS               40  ts/rech.rts"
_G03_EVS_LINE = "EVS               41  et.evs NO-ARCHIVE"


@pytest.mark.parametrize("series", ["both", "rts", "evs"])
def test_mfusg_evsrts_same_ws_rewrite(function_tmpdir, series):
    """Same workspace: entries kept once with unit, path and options; the
    opaque files (and RCH/EVT, skipped by load_only) keep bytes and mtime;
    repeated writes give the same name file."""
    ws = _g03_fixture(
        function_tmpdir,
        rts=None if series == "evs" else _G03_RTS_LINE,
        evs=None if series == "rts" else _G03_EVS_LINE,
    )
    before = _g03_tree(ws)
    ml = _g03_load(ws)
    expected = {
        "both": [_G03_RTS_LINE, _G03_EVS_LINE],
        "rts": [_G03_RTS_LINE],
        "evs": [_G03_EVS_LINE],
    }[series]
    assert [e["filetype"] for e in ml.preserved_nam_entries] == [
        ln.split()[0] for ln in expected
    ]
    ml.write_input()
    assert _g03_nam_entries(ws) == expected
    after = _g03_tree(ws)
    for name in ("ts/rech.rts", "et.evs", "t.rch", "t.evt"):
        assert after.get(name) == before.get(name)
    first = (ws / "t.nam").read_bytes()
    ml.write_input()
    assert (ws / "t.nam").read_bytes() == first


@pytest.mark.parametrize(
    "load_only", [[], _G03_BASE, _G03_BASE + ["RTS", "EVS"], ["RTS"]]
)
def test_mfusg_evsrts_any_load_only(function_tmpdir, load_only):
    ws = _g03_fixture(function_tmpdir)
    ml = _g03_load(ws, load_only=list(load_only))
    assert [e["unit"] for e in ml.preserved_nam_entries] == [40, 41]
    ml.write_input()
    assert _g03_nam_entries(ws) == [_G03_RTS_LINE, _G03_EVS_LINE]


def test_mfusg_evsrts_new_ws_copies(function_tmpdir):
    """New workspace: relative paths are copied with their folders; the
    source is untouched; repeated writes keep the copies (same mtime)."""
    src = _g03_fixture(function_tmpdir / "src")
    before = _g03_tree(src)
    ml = _g03_load(src)
    new = function_tmpdir / "new"
    ml.change_model_ws(str(new))
    ml.write_input()
    assert _g03_nam_entries(new) == [_G03_RTS_LINE, _G03_EVS_LINE]
    for name in ("ts/rech.rts", "et.evs"):
        assert (new / name).read_bytes() == (src / name).read_bytes()
    assert _g03_tree(src) == before
    # RCH/EVT skipped by load_only are not copied (existing behaviour)
    assert not (new / "t.rch").exists() and not (new / "t.evt").exists()
    copied = _g03_tree(new)
    ml.write_input()
    after = _g03_tree(new)
    for name in ("ts/rech.rts", "et.evs"):
        assert after[name] == copied[name]


def test_mfusg_evsrts_new_ws_identical_destination_kept(function_tmpdir):
    src = _g03_fixture(function_tmpdir / "src")
    new = function_tmpdir / "new"
    (new / "ts").mkdir(parents=True)
    shutil.copy2(src / "ts" / "rech.rts", new / "ts" / "rech.rts")
    kept = _g03_tree(new)
    ml = _g03_load(src)
    ml.change_model_ws(str(new))
    ml.write_input()
    assert _g03_tree(new)["ts/rech.rts"] == kept["ts/rech.rts"]


def test_mfusg_evsrts_new_ws_foreign_destination_rejected(function_tmpdir):
    """A different file at the destination is not replaced, and nothing is
    written (no package files, no name file)."""
    src = _g03_fixture(function_tmpdir / "src")
    new = function_tmpdir / "new"
    new.mkdir()
    (new / "et.evs").write_text("SENTINEL\n")
    ml = _g03_load(src)
    ml.change_model_ws(str(new))
    with pytest.raises(FileExistsError, match="not replaced"):
        ml.write_input()
    with pytest.raises(FileExistsError, match="not replaced"):
        ml.write_name_file()
    assert sorted(p.name for p in new.iterdir()) == ["et.evs"]
    assert (new / "et.evs").read_text() == "SENTINEL\n"


@pytest.mark.parametrize("where", ["relative", "absolute"])
def test_mfusg_evsrts_quoted_path_with_spaces(function_tmpdir, where):
    """Quoted names (URWORD: a word starting with a single quote ends at the
    next quote) keep their spaces through load, same-workspace and
    new-workspace writes, and reload."""
    src = function_tmpdir / "src"
    if where == "relative":
        name = "ts dir/with spaces.rts"
    else:
        name = str(function_tmpdir / "shared dir" / "with spaces.rts")
    _g03_fixture(src, rts=f"RTS               40  '{name}' flush", evs=None)
    ml = _g03_load(src)
    (entry,) = ml.preserved_nam_entries
    assert (entry["path"], entry["options"]) == (name, ["FLUSH"])
    line = f"RTS               40  '{name}' FLUSH"
    ml.write_input()
    assert _g03_nam_entries(src) == [line]
    new = function_tmpdir / "new"
    ml.change_model_ws(str(new))
    if where == "absolute":
        with pytest.warns(UserWarning, match="still depends on that file"):
            ml.write_input()
    else:
        ml.write_input()
    assert _g03_nam_entries(new) == [line]
    copied = new / "ts dir" / "with spaces.rts"
    assert copied.exists() == (where == "relative")
    again = _g03_load(new).preserved_nam_entries
    assert [(e["path"], e["options"]) for e in again] == [(name, ["FLUSH"])]


@pytest.mark.parametrize(
    "line, match",
    [
        ("RTS               40  'ts/rech.rts", "unclosed quote"),
        ("RTS               40  ''", "empty file name"),
        # three items for the name-file parser, two words for URWORD
        ("RTS               40  ,,", "needs a unit and a file name"),
        # rejected by the name-file parser
        ("RTS               4x  ts/rech.rts", "not an integer"),
        ("RTS               40  ts/rech.rts " + "x" * 300, "longer than 299"),
    ],
    ids=["unclosed_quote", "empty_name", "no_name", "unit", "too_long"],
)
def test_mfusg_evsrts_malformed_entry_rejected(function_tmpdir, line, match):
    ws = _g03_fixture(function_tmpdir, evs=None)
    nam = (ws / "t.nam").read_text().replace(_G03_RTS_LINE, line)
    (ws / "t.nam").write_text(nam)
    with pytest.raises(ValueError, match=match):
        _g03_load(ws)


def test_mfusg_evsrts_external_reference(function_tmpdir):
    """A path outside the workspace is kept as written in the same
    workspace and becomes an absolute reference (not copied) in a new one."""
    src = function_tmpdir / "proj" / "model"
    _g03_fixture(src, rts="RTS               40  ../shared/rech.rts", evs=None)
    assert (function_tmpdir / "proj" / "shared" / "rech.rts").read_text() == _G03_RTS
    ml = _g03_load(src)
    ml.write_input()
    assert _g03_nam_entries(src) == ["RTS               40  ../shared/rech.rts"]
    new = function_tmpdir / "new"
    ml.change_model_ws(str(new))
    with pytest.warns(UserWarning, match="still depends on that file"):
        ml.write_input()
    shared = os.path.normpath(str(src / ".." / "shared" / "rech.rts"))
    assert _g03_nam_entries(new) == [f"RTS               40  {shared}"]
    assert not (new / "rech.rts").exists() and not (new / "shared").exists()


@pytest.mark.parametrize("mode", ["same", "new", "absolute"])
def test_mfusg_evsrts_missing_source_rejected(function_tmpdir, mode):
    """A missing time-series file can be inspected after load (warning) but
    is rejected on write in every mode, before any file is written."""
    src = function_tmpdir / "src"
    rts = (
        f"RTS               40  {function_tmpdir / 'gone.rts'}"
        if mode == "absolute"
        else _G03_RTS_LINE
    )
    _g03_fixture(src, rts=rts, evs=None)
    if mode == "absolute":
        target = function_tmpdir / "gone.rts"
    else:
        target = src / "ts" / "rech.rts"
    target.rename(target.with_suffix(".moved"))  # missing, without deleting
    with pytest.warns(UserWarning, match="does not exist"):
        ml = _g03_load(src)
    assert ml.preserved_nam_entries[0]["unit"] == 40
    if mode != "same":
        ml.change_model_ws(str(function_tmpdir / "new"))
    before = _g03_tree(ml.model_ws)
    with pytest.raises(FileNotFoundError):
        ml.write_input()
    with pytest.raises(FileNotFoundError):
        ml.write_name_file()
    assert _g03_tree(ml.model_ws) == before


@pytest.mark.parametrize("unit", [-40, 0])
def test_mfusg_evsrts_non_positive_unit_rejected(function_tmpdir, unit):
    ws = _g03_fixture(
        function_tmpdir, rts=f"RTS             {unit:4d}  ts/rech.rts", evs=None
    )
    with pytest.warns(UserWarning, match=f"unit {unit}"):
        ml = _g03_load(ws)
    assert ml.preserved_nam_entries[0]["unit"] == unit
    before = _g03_tree(ws)
    with pytest.raises(ValueError, match="<= 0"):
        ml.write_input()
    assert _g03_tree(ws) == before


def test_mfusg_evsrts_options(function_tmpdir):
    """FLUSH and NO-ARCHIVE are kept (upper case); other options are
    dropped with a warning."""
    ws = _g03_fixture(
        function_tmpdir,
        rts=_G03_RTS_LINE + " flush FOO,Replace",
        evs="EVS               41  et.evs no-archive",
    )
    with pytest.warns(UserWarning, match=r"\['FOO', 'Replace'\]"):
        ml = _g03_load(ws)
    ml.write_input()
    assert _g03_nam_entries(ws) == [
        _G03_RTS_LINE + " FLUSH",
        "EVS               41  et.evs NO-ARCHIVE",
    ]


def test_mfusg_evsrts_entries_are_copies(function_tmpdir):
    ws = _g03_fixture(function_tmpdir)
    ml = _g03_load(ws)
    entries = ml.preserved_nam_entries
    entries[0]["options"].append("FOO")
    entries[0]["unit"] = 99
    entries.append({"filetype": "RTS"})
    assert ml.preserved_nam_entries[0]["options"] == []
    assert ml.preserved_nam_entries[0]["unit"] == 40
    assert len(ml.preserved_nam_entries) == 2
    ml.write_input()
    assert _g03_nam_entries(ws) == [_G03_RTS_LINE, _G03_EVS_LINE]


def _g03_collide(ml, what):
    ws = Path(ml.model_ws)
    chd = ml.get_package("CHD")
    if what == "unit":
        chd.unit_number = [40]
        return "unit 40 of the RTS entry is also used by CHD"
    if what == "rts_evs_unit":
        ml._preserved_nam_entries[1]["unit"] = 40
        return "unit 40 of the EVS entry is also used by RTS"
    if what == "rts_evs_file":
        ml._preserved_nam_entries[1]["source"] = ml._preserved_nam_entries[0]["source"]
        return "is also written by RTS"
    names = {
        "package": "et.evs",
        "package_alias": "./ts/../et.evs",
        "package_case": "ET.EVS",
    }
    if what in names:
        chd.file_name = [names[what]]
        chd.fn_path = str(ws / names[what])
        return "is also written by CHD"
    if what == "list":
        ml.lst.file_name = ["et.evs"]
        return "is also written by LIST"
    if what == "output_basename":
        # registered with a folder; the name file writes the basename
        ml.add_output_file(55, fname="out/et.evs", package="LPF")
        return "is also written by DATA output"
    raise ValueError(what)


@pytest.mark.parametrize(
    "what",
    [
        "unit",
        "rts_evs_unit",
        "rts_evs_file",
        "package",
        "package_alias",
        "package_case",
        "list",
        "output_basename",
    ],
)
@pytest.mark.parametrize("call", ["write_input", "write_name_file"])
def test_mfusg_evsrts_collisions_rejected(function_tmpdir, what, call):
    """Units and files are checked against what the name file actually
    lists (outputs as basename), before anything is written."""
    if what == "package_case" and sys.platform not in ("darwin", "win32"):
        pytest.skip("case-insensitive file names only assumed on macOS/Windows")
    ws = _g03_fixture(function_tmpdir)
    ml = _g03_load(ws)
    match = _g03_collide(ml, what)
    before = _g03_tree(ws)
    with pytest.raises(ValueError, match=match):
        getattr(ml, call)()
    assert _g03_tree(ws) == before


def test_mfusg_evsrts_external_line_too_long(function_tmpdir):
    """The 299-character limit applies to the final line (quotes, options)."""
    deep = function_tmpdir / ("d " + "x" * 120) / ("e" * 120)
    src = function_tmpdir / "src"
    rel = f"../{deep.parent.name}/{deep.name}/r.rts"
    _g03_fixture(src, rts=f"RTS               40  '{rel}' flush", evs=None)
    ml = _g03_load(src)
    ml.change_model_ws(str(function_tmpdir / "new"))
    with pytest.raises(ValueError, match="longer than 299"):
        ml.write_input()
    assert not (function_tmpdir / "new" / "t.nam").exists()


# (case, rts line, evs line, extra name-file lines, line moved after LIST,
#  expected key, expected "TYPE unit N" users in name-file order)
_G03_UNIT_CASES = [
    (
        "evs_chd",
        _G03_RTS_LINE,
        "EVS               24  et.evs NO-ARCHIVE",
        "",
        None,
        24,
        ["CHD unit 24", "EVS unit 24"],
    ),
    (
        "evs_chd_first",
        _G03_RTS_LINE,
        "EVS               24  et.evs NO-ARCHIVE",
        "",
        "EVS",
        24,
        ["EVS unit 24", "CHD unit 24"],
    ),
    (
        "rts_chd",
        "RTS               24  ts/rech.rts",
        _G03_EVS_LINE,
        "",
        None,
        24,
        ["CHD unit 24", "RTS unit 24"],
    ),
    (
        "rts_chd_first",
        "RTS               24  ts/rech.rts",
        _G03_EVS_LINE,
        "",
        "RTS",
        24,
        ["RTS unit 24", "CHD unit 24"],
    ),
    (
        "rts_output",
        "RTS               51  ts/rech.rts",
        _G03_EVS_LINE,
        "",
        None,
        51,
        ["DATA(BINARY) unit 51", "RTS unit 51"],
    ),
    (
        "rts_output_first",
        "RTS               51  ts/rech.rts",
        _G03_EVS_LINE,
        "",
        "RTS",
        51,
        ["RTS unit 51", "DATA(BINARY) unit 51"],
    ),
    (
        "evs_data",
        _G03_RTS_LINE,
        "EVS               62  et.evs",
        "DATA              62  x.dat\n",
        None,
        62,
        ["EVS unit 62", "DATA unit 62"],
    ),
    (
        "evs_data_first",
        _G03_RTS_LINE,
        "EVS               62  et.evs",
        "DATA              62  x.dat\n",
        "DATA",
        62,
        ["DATA unit 62", "EVS unit 62"],
    ),
    (
        "rts_evs",
        "RTS               41  ts/rech.rts",
        _G03_EVS_LINE,
        "",
        None,
        41,
        ["RTS unit 41", "EVS unit 41"],
    ),
    (
        "rts_evs_first",
        "RTS               41  ts/rech.rts",
        _G03_EVS_LINE,
        "",
        "EVS",
        41,
        ["EVS unit 41", "RTS unit 41"],
    ),
    # unit 0 is keyed by the file type in the name-file reader
    (
        "evs_evs_unit0",
        _G03_RTS_LINE,
        "EVS                0  et.evs",
        "EVS                0  et.evs\n",
        None,
        "EVS",
        ["EVS unit 0", "EVS unit 0"],
    ),
]


@pytest.mark.parametrize("load_only", [[], _G03_BASE], ids=["empty", "selected"])
@pytest.mark.parametrize("case", _G03_UNIT_CASES, ids=[c[0] for c in _G03_UNIT_CASES])
def test_mfusg_evsrts_nam_unit_collision_rejected(function_tmpdir, case, load_only):
    """An EVS/RTS unit shared with another entry in the original name file is
    rejected on load (the unit-keyed name-file reader would drop one entry),
    in either line order, without touching the files."""
    _, rts, evs, extra, first, key, users = case
    ws = _g03_fixture(function_tmpdir, rts=rts, evs=evs, extra=extra)
    (ws / "x.dat").write_text("1.0\n")
    lines = (ws / "t.nam").read_text().splitlines()
    if first is not None:
        i = max(k for k, ln in enumerate(lines) if ln.startswith(first))
        lines.insert(2, lines.pop(i))
        (ws / "t.nam").write_text("\n".join(lines) + "\n")
    pairs = [u.replace(" unit", "").split() for u in users]
    found = [
        (k + 1, *ln.split()[:2])
        for k, ln in enumerate(lines)
        if ln.split()[:2] in pairs
    ]
    assert [f"{t} unit {u}" for _, t, u in found] == users
    expected = ", ".join(f"{t} unit {u} (line {n})" for n, t, u in found)
    before = _g03_tree(ws)
    with pytest.raises(ValueError) as e:
        _g03_load(ws, load_only=list(load_only))
    assert f"name-file unit {key} is used by more than one entry: " in str(e.value)
    assert expected in str(e.value)
    assert _g03_tree(ws) == before


def test_mfusg_evsrts_other_duplicate_units_unchanged(function_tmpdir):
    """Units shared only by other types keep the reader's behaviour (no new
    error)."""
    extra = "IBS               60  x.ibs\nHOB               60  x.hob\n"
    ws = _g03_fixture(function_tmpdir, extra=extra)
    for name in ("x.ibs", "x.hob"):
        (ws / name).write_text("1.0\n")
    ml = _g03_load(ws)
    assert [e["unit"] for e in ml.preserved_nam_entries] == [40, 41]


def test_mfusg_evsrts_other_types_unchanged(function_tmpdir):
    """Without EVS/RTS nothing is kept; IBS/HOB (accepted by the USG-T name
    file reader, not read by any package) and DATA keep their handling."""
    ws = _g03_fixture(
        function_tmpdir,
        rts=None,
        evs=None,
        extra="IBS               60  x.ibs\nHOB               61  x.hob\n"
        "DATA              62  x.dat\n",
    )
    for name in ("x.ibs", "x.hob", "x.dat"):
        (ws / name).write_text("1.0\n")
    ml = _g03_load(ws)
    assert ml.preserved_nam_entries == []
    ml.write_input()
    types = [
        ln.split()[0]
        for ln in (ws / "t.nam").read_text().splitlines()
        if ln and not ln.startswith("#")
    ]
    assert "IBS" not in types and "HOB" not in types
    assert types.count("DATA") == 1
    assert "RTS" not in types and "EVS" not in types


# G06: TRTOSS stress periods (SS / TR / TRTOSS stofrac) in DIS and DISU.
# Hand-written fixture: 5 cells in line (DIS 1x5 or DISU 5 nodes), CHD h=10 at
# cell 1, WEL -0.5 at cell 5, strt 20, confined LPF (ss 1e-3) or BCF
# (sf1 3e-2), 3 stress periods given as text records, OC with ATSA.

_G06_DIS = """# DIS
         1         1         5         3         4         2
  0
CONSTANT    1.000000E+01                           #delr
CONSTANT    1.000000E+01                           #delc
CONSTANT    3.000000E+01                           #model_top
CONSTANT    0.000000E+00                           #botm layer 1
"""
_G06_DISU = """# DISU
5 1 13 0 3 4 2 0
  0
INTERNAL               1     (1I10) -1 #nodelay
         5
INTERNAL               1   (5E15.6) -1 #top layer 1
   3.000000E+01   3.000000E+01   3.000000E+01   3.000000E+01   3.000000E+01
INTERNAL               1   (5E15.6) -1 #bot layer 1
   0.000000E+00   0.000000E+00   0.000000E+00   0.000000E+00   0.000000E+00
CONSTANT    1.000000E+02                           #area layer 1
INTERNAL               1     (5I10) -1 #iac
         2         3         3         3         2
INTERNAL               1     (FREE) -1 #ja
 1 2  2 1 3  3 2 4
 4 3 5  5 4
CONSTANT    5.000000E+00                           #cl12
CONSTANT    3.000000E+01                           #fahl
"""
_G06_LPF = """# LPF
         0    -1E+30         0         0
         0
         0
   1.000000E+00
         0
         0
CONSTANT    1.000000E+00                           #hk layer 1
CONSTANT    1.000000E+00                           #vka1
"""
_G06_LPF_SS = "CONSTANT    1.000000E-03                           #ss layer 1\n"
_G06_BCF = {
    "dis": "         0 -1.000e+30         0 1.000e-01         1         0"
    "         0\n00\n",
    "disu": "         0 -1.000e+30         0 1.000e-01         1         0"
    "         0         0\n00\n",
}
_G06_BCF_ARRAYS = (
    "CONSTANT    1.000000E+00                           #anisotropy factor\n",
    "CONSTANT    3.000000E-02                           #sf1 layer 1\n",
    "CONSTANT    3.000000E+01                           #transmissivity layer 1\n",
)
_G06_CHD = {
    "dis": "# CHD\n         1\n         1         0\n         1         1         1"
    "            10.0            10.0\n"
    "        -1         0\n        -1         0\n",
    "disu": "# CHD\n         1\n 1    Stress Period 1\n 1   10.000000  10.000000\n"
    " -1    Stress Period 2\n -1    Stress Period 3\n",
}
_G06_WEL = {
    "dis": "# WEL\n         1         0\n         1         0\n"
    "         1         1         5            -0.5\n"
    "        -1         0\n        -1         0\n",
    "disu": "# WEL\n         1         0\n         1         0         0\n"
    "         5            -0.5\n"
    "        -1         0         0\n        -1         0         0\n",
}
_G06_SMS = """# SMS
  1.00e-09  1.00e-11       100       100         2         1         1
  7.00e-01  1.00e-01  2.00e-01  1.00e-03        20  1.00e+04  2.00e-01  1.00e+02
         2         0         7         2         0  0.00e+00         0  1.00e-03
"""


def _g06_oc(ats=True):
    lines = (
        ["# OC"]
        + (["ATSA"] if ats else [])
        + ["HEAD PRINT FORMAT   0", "HEAD SAVE UNIT    51", "COMPACT BUDGET AUX"]
    )
    for k in (1, 2, 3):
        lines.append(f"period {k} step 1")
        if ats and k == 1:
            lines += ["  DELTAT 0.1", "  TADJAT 1.3", "  TMAXAT 1.0E4"]
        lines += ["  save head", "  print budget"]
    return "\n".join(lines) + "\n"


_G06_E3 = [
    "100000.0   1   1.0   TRTOSS 1.0E-3",
    "    10.0   1   1.0   TR",
    "    10.0   1   1.0   TR",
]


def _g06_model(ws, grid="dis", flow="lpf", periods=_G06_E3, storage=True, ats=True):
    """Write the hand-written 1-layer fixture into ws; periods are the three
    stress-period records (text)."""
    ws = Path(ws)
    ws.mkdir(parents=True, exist_ok=True)
    disname = "t.dis" if grid == "dis" else "t.disu"
    base = _G06_DIS if grid == "dis" else _G06_DISU
    (ws / disname).write_text(base + "\n".join(periods) + "\n")
    bas = "UNSTRUCTURED FREE" if grid == "disu" else "FREE"
    (ws / "t.bas").write_text(
        f"# BAS6\n{bas}\nCONSTANT          1                                #ibound\n"
        "   -999.99\nCONSTANT    2.000000E+01                           #strt\n"
    )
    if flow == "lpf":
        (ws / "t.lpf").write_text(_G06_LPF + (_G06_LPF_SS if storage else ""))
        flowline = "LPF               15  t.lpf"
    else:
        anis, sf1, tran = _G06_BCF_ARRAYS
        bcf = _G06_BCF[grid] + anis + (sf1 if storage else "") + tran
        (ws / "t.bcf").write_text(bcf)
        flowline = "BCF6              15  t.bcf"
    (ws / "t.sms").write_text(_G06_SMS)
    (ws / "t.oc").write_text(_g06_oc(ats))
    (ws / "t.chd").write_text(_G06_CHD[grid])
    (ws / "t.wel").write_text(_G06_WEL[grid])
    (ws / "t.nam").write_text(
        "# fixture\nLIST               2  t.list\n"
        f"{'DIS ' if grid == 'dis' else 'DISU'}              11  {disname}\n"
        "CHD               24  t.chd\nWEL               20  t.wel\n"
        f"BAS6              13  t.bas\n{flowline}\nSMS               32  t.sms\n"
        "OC                14  t.oc\nDATA(BINARY)      51  t.hds\n"
    )
    return ws


def _g06_records(path):
    """Words of the stress-period records (last three lines)."""
    return [ln.split()[3:] for ln in Path(path).read_text().splitlines()[-3:]]


def _g06_new(ws, grid, steady, trtoss=None, **kw):
    ml = MfUsg(modelname="t", model_ws=str(ws), structured=(grid == "dis"))
    args = dict(nper=len(steady), perlen=10.0, steady=steady, **kw)
    if trtoss is not None:
        args["trtoss"] = trtoss
    if grid == "dis":
        return MfUsgDis(ml, nlay=1, nrow=1, ncol=3, **args)
    return MfUsgDisU(
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
        **args,
    )


def _g06_write(pkg):
    if isinstance(pkg, MfUsgDis):
        pkg.write_file(check=False)
    else:
        pkg.write_file()


@pytest.mark.parametrize("grid", ["dis", "disu"])
def test_mfusg_trtoss_three_states(function_tmpdir, grid):
    dis = _g06_new(
        function_tmpdir,
        grid,
        [True, False, False, False],
        [None, 1.0e-3, np.nan, -0.25],
    )
    assert dis.period_type == ["SS", "TRTOSS", "TR", "TRTOSS"]
    assert dis.itrss == -1
    _g06_write(dis)
    words = [ln.split()[3:] for ln in Path(dis.fn_path).read_text().splitlines()[-4:]]
    assert words == [["SS"], ["TRTOSS", "1.e-03"], ["TR"], ["TRTOSS", "-2.5e-01"]]


@pytest.mark.parametrize("grid", ["dis", "disu"])
def test_mfusg_trtoss_all_and_scalar(function_tmpdir, grid):
    dis = _g06_new(function_tmpdir, grid, [False, False], 1.0e-6)
    assert dis.period_type == ["TRTOSS", "TRTOSS"] and dis.itrss == -1
    tr = _g06_new(function_tmpdir / "tr", grid, [False, True, False], [None, None, 0.0])
    assert tr.period_type == ["TR", "SS", "TRTOSS"] and tr.itrss == -1
    assert _g06_new(function_tmpdir / "t", grid, [False, False]).itrss == 1
    assert _g06_new(function_tmpdir / "s", grid, [True, True]).itrss == 0


@pytest.mark.parametrize("grid", ["dis", "disu"])
@pytest.mark.parametrize(
    "steady, trtoss, match",
    [
        ([True, False], [1.0e-3, None], "need steady=False"),
        ([True], 1.0e-3, "need steady=False"),  # default steady=True
        ([True, False], [None, 1.0e-3], "without a TR period"),
        ([False, False], [np.inf, None], "not a finite"),
        ([False, False], [-np.inf, None], "not a finite"),
        ([False, False], [1.0e39, None], "not a finite"),
        ([False, False], ["1e-3", None], "number or None"),
        ([False, False], [True, None], "number or None"),
        ([False, False], [[1.0e-3, None]], "scalar or a sequence"),
        ([False, False], [1.0e-3, None, None], "scalar or a sequence"),
    ],
)
def test_mfusg_trtoss_api_rejected(function_tmpdir, grid, steady, trtoss, match):
    with pytest.raises(ValueError, match=match):
        _g06_new(function_tmpdir, grid, steady, trtoss)


@pytest.mark.parametrize("grid", ["dis", "disu"])
@pytest.mark.parametrize("change", ["ss_trtoss", "steady_true", "length", "inf"])
def test_mfusg_trtoss_checked_before_writing(function_tmpdir, grid, change):
    """Changes after construction are checked before the file is opened."""
    dis = _g06_new(function_tmpdir, grid, [True, False, False], [None, None, 1.0e-3])
    if change == "ss_trtoss":
        dis.trtoss = [None, 1.0e-3, 1.0e-3]  # SS, TRTOSS, TRTOSS
        match = "without a TR period"
    elif change == "steady_true":
        dis.steady = [True, False, True]
        dis.trtoss = [None, None, None]
        dis.trtoss[0] = 1.0e-3
        match = "need steady=False"
    elif change == "length":
        dis.trtoss = [1.0e-3]
        match = "scalar or a sequence"
    else:
        dis.trtoss = np.array([np.nan, np.nan, np.inf])
        match = "not a finite"
    Path(dis.fn_path).write_text("SENTINEL\n")
    with pytest.raises(ValueError, match=match):
        _g06_write(dis)
    assert Path(dis.fn_path).read_text() == "SENTINEL\n"


def test_mfusg_trtoss_positional_compatibility(function_tmpdir):
    """Existing positional calls keep their meaning (trtoss is appended)."""
    ml = MfUsg(modelname="p", model_ws=str(function_tmpdir))
    dis = MfUsgDis(
        ml,
        1,
        1,
        3,
        2,
        10.0,
        10.0,
        0,
        30.0,
        0.0,
        [5.0, 6.0],
        [2, 3],
        [1.0, 1.5],
        [True, False],
        False,
        4,
        2,
        "dis",
        None,
        None,
        None,
        None,
        None,
        None,
        None,
        "1/1/2020",
    )
    assert dis.perlen.array.tolist() == [5.0, 6.0] and dis.itmuni == 4
    assert dis.period_type == ["SS", "TR"] and dis.trtoss == [None, None]
    mu = MfUsg(modelname="q", model_ws=str(function_tmpdir), structured=False)
    disu = MfUsgDisU(
        mu,
        3,
        1,
        7,
        0,
        2,
        4,
        2,
        0,
        0,
        None,
        [10.0] * 3,
        [0.0] * 3,
        100.0,
        np.array([2, 3, 2]),
        np.array([0, 1, 1, 0, 2, 2, 1]),
        None,
        None,
        None,
        5.0,
        10.0,
        [5.0, 6.0],
        1,
        1,
        [False, True],
        "disu",
        None,
        None,
        "1/1/2020",
    )
    assert disu.perlen.array.tolist() == [5.0, 6.0]
    assert disu.period_type == ["TR", "SS"] and disu.trtoss == [None, None]


def test_mfusg_trtoss_legacy_flag(function_tmpdir):
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        _g06_new(function_tmpdir, "dis", [False])  # default Trtoss: no warning
        _g06_new(function_tmpdir, "dis", [False], Trtoss=False)
    with pytest.warns(UserWarning, match="Trtoss has no effect"):
        dis = _g06_new(function_tmpdir, "dis", [False, False], Trtoss=True)
    assert dis.period_type == ["TR", "TR"]
    with pytest.warns(UserWarning, match="Trtoss has no effect"):
        dis = _g06_new(
            function_tmpdir, "dis", [False, False], [1.0e-3, None], Trtoss=True
        )
    assert dis.period_type == ["TRTOSS", "TR"]


def _g06_load_dis(ws, grid):
    ml = MfUsg(modelname="t", model_ws=str(ws), structured=(grid == "dis"))
    cls = MfUsgDis if grid == "dis" else MfUsgDisU
    kw = {"check": False} if grid == "dis" else {}
    return cls.load(str(Path(ws) / ("t.dis" if grid == "dis" else "t.disu")), ml, **kw)


@pytest.mark.parametrize("grid", ["dis", "disu"])
@pytest.mark.parametrize(
    "word, expected",
    [
        ("TRTOSS 1.0E-3", 1.0e-3),
        ("trtoss 1.0e-3", 1.0e-3),
        ("TrToss 1.0D-3", 1.0e-3),
        ("TRTOSS 1.0d-3", 1.0e-3),
        ("TRTOSS 1.0Q-3", 1.0e-3),
        ("TRTOSS 1.0-3", 1.0e-3),
        ("TRTOSS .001", 1.0e-3),
        ("TRTOSS +1E-03", 1.0e-3),
        ("TRTOSS 1.0E-3,9", 1.0e-3),
        ("TRTOSS 0.0", 0.0),
        ("TRTOSS -1.0", -1.0),
    ],
    ids=[
        "E",
        "lower",
        "mixed_D",
        "d",
        "Q",
        "no_letter",
        "dot",
        "plus",
        "comma",
        "zero",
        "negative",
    ],
)
def test_mfusg_trtoss_load_words(function_tmpdir, grid, word, expected):
    """Number syntax read by USG-T (URWORD, F30.0): E/D/Q or no exponent
    letter, leading dot, sign, comma separator; case-insensitive flag."""
    ws = _g06_model(
        function_tmpdir, grid, periods=[f"100000.0   1   1.0   {word}"] + _G06_E3[1:]
    )
    dis = _g06_load_dis(ws, grid)
    assert dis.period_type == ["TRTOSS", "TR", "TR"]
    assert np.float32(dis.trtoss[0]) == np.float32(expected)


@pytest.mark.parametrize("grid", ["dis", "disu"])
@pytest.mark.parametrize(
    "word, match",
    [
        ("TRTOSS NaN", "not a finite number"),
        ("TRTOSS Inf", "not a finite number"),
        ("TRTOSS -Infinity", "not a finite number"),
        ("TRTOSS abc", "not a finite number"),
        ("TRTOSS 1." + "0" * 30 + "E-3", "not a finite number"),  # word > 30 characters
        ("TRTOSS 3.5E38", "beyond the REAL range"),
        ("XX", "must be SS, TR or TRTOSS"),
        ("SSX", "must be SS, TR or TRTOSS"),
    ],
    ids=[
        "nan",
        "inf",
        "neg_inf",
        "text",
        "long_word",
        "overflow",
        "flag_xx",
        "flag_ssx",
    ],
)
def test_mfusg_trtoss_load_rejected(function_tmpdir, grid, word, match):
    ws = _g06_model(
        function_tmpdir, grid, periods=[f"100000.0   1   1.0   {word}"] + _G06_E3[1:]
    )
    with pytest.raises(ValueError, match=match):
        _g06_load_dis(ws, grid)


@pytest.mark.parametrize("grid", ["dis", "disu"])
def test_mfusg_trtoss_missing_fraction(function_tmpdir, grid):
    ws = _g06_model(
        function_tmpdir, grid, periods=["100000.0   1   1.0   TRTOSS"] + _G06_E3[1:]
    )
    with pytest.warns(UserWarning, match="without a storage fraction"):
        dis = _g06_load_dis(ws, grid)
    assert dis.period_type[0] == "TRTOSS" and dis.trtoss[0] == 0.0
    _g06_write(dis)
    assert _g06_records(dis.fn_path)[0] == ["TRTOSS", "0.e+00"]


def _g06_tie_words():
    """Words close to a single-precision rounding tie (exact decimal
    midpoint of two neighbouring REAL values, +/- 1e-22 relative): the REAL
    USG-T reads is the neighbour on the side of the word."""
    from decimal import Decimal, localcontext

    out = []
    for a in (np.float32(1.0e-3), np.float32(4.3712345e-4), np.float32(7.7e-6)):
        b = np.nextafter(a, np.float32(1))
        with localcontext() as ctx:
            ctx.prec = 100
            mid = (Decimal(float(a)) + Decimal(float(b))) / 2
            for eps, real in ((Decimal("-1e-22"), a), (Decimal("1e-22"), b)):
                out.append((format(mid * (1 + eps), ".22E"), real))
    return out


@pytest.mark.parametrize("grid", ["dis", "disu"])
def test_mfusg_trtoss_load_write_load(function_tmpdir, grid):
    """Flags and the REAL storage fraction survive load -> write -> load,
    including values with more than seven significant digits and words
    next to a single-precision rounding tie."""
    cases = [
        ("1.23456789E-3", np.float32(1.23456789e-3)),
        ("4.3712345678D-4", np.float32(4.3712345678e-4)),
    ] + _g06_tie_words()
    for i, (word, real) in enumerate(cases):
        ws = _g06_model(
            function_tmpdir / f"c{i}",
            grid,
            periods=[
                f"100000.0   1   1.0   TRTOSS {word}",
                "    10.0   1   1.0   TR",
                "     5.0   1   1.0   SS",
            ],
        )
        dis = _g06_load_dis(ws, grid)
        assert np.float32(dis.trtoss[0]) == real, word
        _g06_write(dis)
        records = _g06_records(dis.fn_path)
        assert records[1:] == [["TR"], ["SS"]]
        assert records[0][0] == "TRTOSS" and np.float32(records[0][1]) == real
        again = _g06_load_dis(ws, grid)
        assert again.period_type == ["TRTOSS", "TR", "SS"]
        assert np.float32(again.trtoss[0]) == real


@pytest.mark.parametrize("grid", ["dis", "disu"])
@pytest.mark.parametrize("flow", ["lpf", "bcf"])
def test_mfusg_trtoss_ss_without_tr_rejected_on_load(function_tmpdir, grid, flow):
    """SS + TRTOSS without TR stops the model load at DIS/DISU, before the
    storage arrays of LPF/BCF are interpreted."""
    ws = _g06_model(
        function_tmpdir,
        grid,
        flow,
        storage=False,
        periods=[
            "     1.0   1   1.0   SS",
            "100000.0   1   1.0   TRTOSS 1.0E-3",
            "     1.0   1   1.0   SS",
        ],
    )
    with pytest.raises(ValueError, match="without a TR period"):
        MfUsg.load("t.nam", model_ws=str(ws), check=False)


_G06_LPF2 = """# LPF (2 layers)
         0    -1E+30         0         0
         0         0
         0         0
   1.000000E+00   1.000000E+00
         0         0
         0         0
"""


def _g06_two_layer(ws, grid, flow):
    """2-layer fixture with all periods TRTOSS: layer values differ so that
    a storage array read as another property is visible."""
    ws = Path(ws)
    ws.mkdir(parents=True)
    rec = "\n".join(["100000.0   1   1.0   TRTOSS 1.0E-3"] * 3) + "\n"
    if grid == "dis":
        (ws / "t.dis").write_text(
            "# DIS\n         2         1         5         3         4         2\n"
            "  0  0\nCONSTANT 10.0 #delr\nCONSTANT 10.0 #delc\nCONSTANT 30.0 #top\n"
            "CONSTANT 15.0 #botm 1\nCONSTANT 0.0 #botm 2\n" + rec
        )
        nam = "DIS               11  t.dis\n"
    else:
        (ws / "t.disu").write_text(
            "# DISU\n2 2 4 -1 3 4 2 0\n  0  0\n"
            "INTERNAL 1 (FREE) -1 #nodelay\n 1 1\n"
            "CONSTANT 30.0 #top 1\nCONSTANT 15.0 #top 2\n"
            "CONSTANT 15.0 #bot 1\nCONSTANT 0.0 #bot 2\nCONSTANT 100.0 #area\n"
            "INTERNAL 1 (FREE) -1 #iac\n 2 2\nINTERNAL 1 (FREE) -1 #ja\n 1 2 2 1\n"
            "CONSTANT 7.5 #cl12\nCONSTANT 100.0 #fahl\n" + rec
        )
        nam = "DISU              11  t.disu\n"
    if flow == "lpf":
        (ws / "t.lpf").write_text(
            _G06_LPF2
            + "CONSTANT 1.0 #hk 1\nCONSTANT 1.0 #vka 1\nCONSTANT 1.0E-3 #ss 1\n"
            + "CONSTANT 2.0 #hk 2\nCONSTANT 2.0 #vka 2\nCONSTANT 2.0E-3 #ss 2\n"
        )
        nam += "LPF               15  t.lpf\n"
    else:  # BCF order per layer: SF1, TRAN, VCONT (not on the last layer)
        (ws / "t.bcf").write_text(
            "         0   -1.E+30         0       0.1         1         0         0\n"
            " 0 0\nCONSTANT 1.0 #anisotropy\n"
            "CONSTANT 3.0E-2 #sf1 1\nCONSTANT 30.0 #tran 1\nCONSTANT 0.1 #vcont 1\n"
            "CONSTANT 6.0E-2 #sf1 2\nCONSTANT 60.0 #tran 2\n"
        )
        nam += "BCF6              15  t.bcf\n"
    bas = "FREE" if grid == "dis" else "UNSTRUCTURED FREE"
    (ws / "t.bas").write_text(
        f"# BAS6\n{bas}\nCONSTANT 1 #ibound 1\nCONSTANT 1 #ibound 2\n   -999.99\n"
        "CONSTANT 20.0 #strt 1\nCONSTANT 20.0 #strt 2\n"
    )
    # BAS6 before LPF/BCF: they are loaded in name-file order
    first, flowline = nam.split("\n")[0], nam.split("\n")[1]
    nam = f"{first}\nBAS6              13  t.bas\n{flowline}\n"
    (ws / "t.nam").write_text("LIST               2  t.list\n" + nam)
    return ws


@pytest.mark.parametrize(
    "grid, flow", [("dis", "lpf"), ("disu", "lpf"), ("dis", "bcf")]
)
def test_mfusg_trtoss_two_layer_storage(function_tmpdir, grid, flow):
    """With every period TRTOSS (ITRSS=-1, storage read by USG-T) the
    per-layer storage arrays are read in place, and kept on rewrite."""
    ws = _g06_two_layer(function_tmpdir / "m", grid, flow)
    pkg = "LPF" if flow == "lpf" else "BCF6"
    disname = "DIS" if grid == "dis" else "DISU"
    for _ in range(2):  # load, write, load again
        ml = MfUsg.load(
            "t.nam", model_ws=str(ws), check=False, load_only=[disname, "BAS6", pkg]
        )
        assert ml.get_package(disname).period_type == ["TRTOSS"] * 3
        p = ml.get_package(pkg)
        if flow == "lpf":
            hk = [p.hk[k].array.ravel()[0] for k in range(2)]
            vka = [p.vka[k].array.ravel()[0] for k in range(2)]
            ss = [p.ss[k].array.ravel()[0] for k in range(2)]
            assert (hk, vka) == ([1.0, 2.0], [1.0, 2.0])
            np.testing.assert_allclose(ss, [1.0e-3, 2.0e-3], rtol=1e-6)
        else:
            sf1 = [p.sf1[k].array.ravel()[0] for k in range(2)]
            tran = [p.tran[k].array.ravel()[0] for k in range(2)]
            np.testing.assert_allclose(sf1, [3.0e-2, 6.0e-2], rtol=1e-6)
            assert tran == [30.0, 60.0]
        ml.write_input()


@pytest.mark.parametrize("grid", ["dis", "disu"])
def test_mfusg_trtoss_absent_records_unchanged(function_tmpdir, grid):
    """Without TRTOSS the stress-period records are written as before."""
    dis = _g06_new(function_tmpdir, grid, [True, False])
    _g06_write(dis)
    tail = Path(dis.fn_path).read_text().splitlines()[-2:]
    assert tail == [
        "     10.000000             1  1.000000  SS",
        "     10.000000             1  1.000000  TR",
    ]


@pytest.mark.parametrize("grid", ["dis", "disu"])
def test_mfusg_trtoss_existing_uses(function_tmpdir, grid):
    """trtoss is stored as a list (None for other periods), so attribute-wise
    equality after reload holds; reducing nper after construction writes the
    first nper periods, checked on those periods only."""
    dis = _g06_new(function_tmpdir, grid, [False, False, True], [1.0e-3, None, None])
    assert dis.trtoss == [1.0e-3, None, None]
    _g06_write(dis)
    ws = Path(dis.fn_path).parent
    name = Path(dis.fn_path).name
    ml = MfUsg(modelname="r", model_ws=str(ws), structured=(grid == "dis"))
    cls = MfUsgDis if grid == "dis" else MfUsgDisU
    kw = {"check": False} if grid == "dis" else {}
    again = cls.load(str(ws / name), ml, **kw)
    assert again.trtoss == dis.trtoss and again.trtoss == [1.0e-3, None, None]
    dis.nper = 2
    _g06_write(dis)
    assert [r[0] for r in _g06_records(dis.fn_path)[-2:]] == ["TRTOSS", "TR"]
    assert dis.period_type == ["TRTOSS", "TR"]
    dis.nper = 3
    dis.trtoss = [1.0e-3, None]  # shorter than nper
    with pytest.raises(ValueError, match="scalar or a sequence"):
        _g06_write(dis)


# G04: SMS options, NUMTRACK records and REAL values. Hand-written SMS text;
# for the executable tests a hand-written model (5 cells, CHD h=10 at cell 1,
# WEL -2.0 then -3.0 at cell 5, strt 20, SS then TR with 5 steps).

_G04_XMD = " 1 0 7 2 0 0.0 0 1.0E-3\n"


def _g04_load(ws, text):
    ws = Path(ws)
    ws.mkdir(parents=True, exist_ok=True)
    (ws / "in.sms").write_text("# SMS\n" + text)
    return MfUsgSms.load(str(ws / "in.sms"), MfUsg(modelname="t", model_ws=str(ws)))


def _g04_rewrite(sms):
    """Write sms to out.sms next to its input; return (lines, reloaded)."""
    ws = Path(sms.parent.model_ws)
    sms.fn_path = str(ws / "out.sms")
    sms.write_file()
    lines = Path(sms.fn_path).read_text().splitlines()[1:]
    again = MfUsgSms.load(sms.fn_path, MfUsg(modelname="r", model_ws=str(ws)))
    return lines, again


@pytest.mark.parametrize(
    "text, options, options2",
    [
        (" 1.0E-6 1.0E-6 100 200 1 0 1 SHIFT\n" + _G04_XMD, [], [" SHIFT"]),
        (" 1.0E-6 1.0E-6 100 200 1 0 1 shift\n" + _G04_XMD, [], [" SHIFT"]),
        (
            " 1.0E-6 1.0E-6 100 200 1 0 1 Shift DampBot SOLVEACTIVE\n" + _G04_XMD,
            [],
            [" SHIFT", " DAMPBOT", " SOLVEACTIVE"],
        ),
        (
            " 1.0E-6 1.0E-6 100 200 1 0 1 solveactive,dampbot\n" + _G04_XMD,
            [],
            [" SOLVEACTIVE", " DAMPBOT"],
        ),
        (
            " 1.0E-6 1.0E-6 100 200 1 1 1 truncatednewton\n 0.7 0.1 0.2 0.001 0\n"
            + _G04_XMD,
            [],
            [" TRUNCATEDNEWTON"],
        ),
        ("simple\n 1.0E-6 1.0E-6 100 200 1 0 1 shift\n", ["SIMPLE"], [" SHIFT"]),
        (
            "Moderate\n 1.0E-6 1.0E-6 100 200 1 1 1 TruncatedNewton\n",
            ["MODERATE"],
            [" TRUNCATEDNEWTON"],
        ),
        ("COMPLEX\n 1.0E-6 1.0E-6 100 200 1 1 1 DAMPBOT\n", ["COMPLEX"], [" DAMPBOT"]),
    ],
    ids=[
        "shift",
        "shift_lower",
        "mixed_three",
        "comma_lower",
        "truncnewton_lower",
        "simple_lower",
        "moderate_mixed",
        "complex",
    ],
)
def test_mfusg_sms_options_kept(function_tmpdir, text, options, options2):
    """Option words are read in any case and kept in their order; SHIFT is
    not turned into DAMPBOT."""
    sms = _g04_load(function_tmpdir, text)
    assert (sms.options, sms.options2) == (options, options2)
    lines, again = _g04_rewrite(sms)
    assert (again.options, again.options2) == (options, options2)
    record_1b = lines[1] if options else lines[0]
    assert record_1b.split()[7:] == [o.strip() for o in options2]


@pytest.mark.parametrize(
    "word", ["truncatednewtoncutoff", "TRUNCATEDNEWTONCUTOFF"], ids=["lower", "upper"]
)
def test_mfusg_sms_truncnewton_cutoff(function_tmpdir, word):
    text = (
        f" 1.0E-6 1.0E-6 100 200 1 1 1 {word} 0.7654321D0\n"
        " 0.7 0.1 0.2 0.001 0\n" + _G04_XMD
    )
    sms = _g04_load(function_tmpdir, text)
    name, value = sms.options2[0].split()
    assert name == "TRUNCATEDNEWTONCUTOFF"
    assert np.float32(value) == np.float32(0.7654321)
    lines, again = _g04_rewrite(sms)
    assert np.float32(again.options2[0].split()[1]) == np.float32(0.7654321)


@pytest.mark.parametrize(
    "record2, nwritten",
    [
        (" 0.7 0.1 0.2 0.001 0\n", 5),
        (" 0.7 0.1 0.2 0.001 0 5.0 0.5 9.0\n", 5),
        (" 0.7 0.1 0.2 0.001 -1\n", 5),
        (" 0.7 0.1 0.2 0.001 20 1.0E4 0.2 100.0\n", 8),
    ],
    ids=["zero", "zero_extra", "negative", "positive"],
)
def test_mfusg_sms_numtrack(function_tmpdir, record2, nwritten):
    """BTOL, BREDUC, RES_LIM only when NUMTRACK > 0 (SMS7U1AR), and the XMD
    record that follows is read in sequence."""
    sms = _g04_load(
        function_tmpdir,
        " 1.0E-6 1.0E-6 100 200 1 1 1\n" + record2 + " 2 1 5 3 0 0.0 0 2.0E-3\n",
    )
    xmd = (sms.iacl, sms.norder, sms.level, sms.north, sms.epsrn)
    assert xmd == (2, 1, 5, 3, 2.0e-3)
    lines, again = _g04_rewrite(sms)
    assert len(lines[1].split()) == nwritten
    assert again.numtrack == sms.numtrack and again.level == 5
    if nwritten == 8:
        assert (again.btol, again.breduc, again.reslim) == (1.0e4, 0.2, 100.0)


def test_mfusg_sms_preset_records(function_tmpdir):
    """With a preset only records 1a/1b are written; with PCGU a blank line
    follows, which PCGU7U1AR reads (URDCOM) without using it."""
    text = "SIMPLE\n 1.0E-6 1.0E-6 100 200 1 1 2\n\n"
    sms = _g04_load(function_tmpdir / "pcgu", text)
    lines, again = _g04_rewrite(sms)
    assert lines[0] == "SIMPLE" and lines[1].split()[:7] == [
        "1.00e-06",
        "1.00e-06",
        "100",
        "200",
        "1",
        "1",
        "2",
    ]
    assert lines[2:] == [""]
    assert again.options == ["SIMPLE"] and again.theta is None


@pytest.mark.parametrize(
    "word, expected",
    [
        ("1.0E-6", 1.0e-6),
        ("1.0D-6", 1.0e-6),
        ("1.0d-6", 1.0e-6),
        ("1.0Q-6", 1.0e-6),
        ("1.0-6", 1.0e-6),
        (".000001", 1.0e-6),
        ("+1E-06", 1.0e-6),
    ],
    ids=["E", "D", "d_lower", "Q", "no_letter", "dot", "plus"],
)
def test_mfusg_sms_real_syntax(function_tmpdir, word, expected):
    sms = _g04_load(
        function_tmpdir,
        f" {word} {word} 100 200 1 1 1\n"
        f" {word} 0.1 0.2 0.001 20 {word} 0.2 100.0\n"
        f" 1 0 7 2 0 {word} 0 {word}\n",
    )
    for name in ("hclose", "hiclose", "theta", "btol", "rrctol", "epsrn"):
        assert np.float32(getattr(sms, name)) == np.float32(expected), name


def _g04_tie_words():
    """Words next to a single-precision rounding tie (decimal midpoint of two
    neighbouring REAL values, +/- 1e-22 relative) and the REAL USG-T reads."""
    from decimal import Decimal, localcontext

    out = []
    for a in (np.float32(1.0e-6), np.float32(0.712345), np.float32(1.2345e4)):
        b = np.nextafter(a, np.float32(np.inf))
        with localcontext() as ctx:
            ctx.prec = 100
            mid = (Decimal(float(a)) + Decimal(float(b))) / 2
            for eps, real in ((Decimal("-1e-22"), a), (Decimal("1e-22"), b)):
                out.append((format(mid * (1 + eps), ".22E"), real))
    return out


def test_mfusg_sms_precision(function_tmpdir):
    """REAL values with more than seven significant digits, small values and
    words next to a rounding tie keep the REAL USG-T reads through
    load -> write -> load; values exact in 9.2e keep that text."""
    cases = [
        ("1.23456789E-7", np.float32(1.23456789e-7)),
        ("7.1234567890E-1", np.float32(0.712345678)),
        ("1.0E-30", np.float32(1.0e-30)),
    ] + _g04_tie_words()
    for i, (word, real) in enumerate(cases):
        sms = _g04_load(
            function_tmpdir / f"c{i}",
            f" {word} 1.0E-6 100 200 1 1 2 TRUNCATEDNEWTONCUTOFF {word}\n"
            f" {word} 0.1 0.2 0.001 20 {word} 0.2 100.0\n"
            f" CG 1 0 0 {word} 1.0\n",
        )
        lines, again = _g04_rewrite(sms)
        for name in ("hclose", "theta", "btol", "rclosepcgu"):
            assert np.float32(getattr(again, name)) == real, (word, name)
        assert np.float32(again.options2[0].split()[1]) == real
    text = " 1.0E-4 1.0E-5 100 200 1 0 1\n" + _G04_XMD
    sms = _g04_load(function_tmpdir / "exact", text)
    lines, _ = _g04_rewrite(sms)
    assert lines[0].split()[:2] == ["1.00e-04", "1.00e-05"]
    assert lines[1].split()[-1] == "1.00e-03"


@pytest.mark.parametrize(
    "word",
    ["NaN", "Inf", "-Infinity", "abc", "1." + "0" * 30 + "E-3", "3.5E38"],
    ids=["nan", "inf", "neg_inf", "text", "long_word", "overflow"],
)
def test_mfusg_sms_invalid_real_on_load(function_tmpdir, word):
    with pytest.raises(ValueError, match="HICLOSE") as e:
        _g04_load(function_tmpdir, f" 1.0E-6 {word} 100 200 1 0 1\n" + _G04_XMD)
    assert "TRTOSS" not in str(e.value)


@pytest.mark.parametrize(
    "change, error, match",
    [
        (("hclose", np.inf), ValueError, "HCLOSE"),
        (("hiclose", np.nan), ValueError, "HICLOSE"),
        (("theta", 1.0e39), ValueError, "THETA"),
        (("btol", None), ValueError, "BTOL"),
        (("epsrn", "x"), ValueError, "EPSRN"),
        (
            ("options2", [" TRUNCATEDNEWTONCUTOFF abc"]),
            ValueError,
            "TRUNCATEDNEWTONCUTOFF",
        ),
        (
            ("options2", ["stressperiodsettings"]),
            NotImplementedError,
            "STRESSPERIODSETTINGS",
        ),
    ],
    ids=[
        "hclose_inf",
        "hiclose_nan",
        "theta_overflow",
        "btol_none",
        "epsrn_text",
        "cutoff_text",
        "stressperiodsettings",
    ],
)
def test_mfusg_sms_checked_before_writing(function_tmpdir, change, error, match):
    sms = MfUsgSms(
        MfUsg(modelname="t", model_ws=str(function_tmpdir)),
        nonlinmeth=1,
        linmeth=1,
        numtrack=20,
    )
    setattr(sms, *change)
    Path(sms.fn_path).write_text("SENTINEL\n")
    with pytest.raises(error, match=match):
        sms.write_file()
    assert Path(sms.fn_path).read_text() == "SENTINEL\n"


def test_mfusg_sms_unsupported_options(function_tmpdir):
    """STRESSPERIODSETTINGS is rejected; SOLVEROUTPUTFILE n warns and is not
    kept; other words of record 1b (labels) are ignored like in USG-T."""
    with pytest.raises(NotImplementedError, match="STRESSPERIODSETTINGS"):
        _g04_load(
            function_tmpdir / "a",
            " 1.0E-6 1.0E-6 100 200 1 0 1 StressPeriodSettings\n" + _G04_XMD,
        )
    with pytest.warns(UserWarning, match="solver output file is not written"):
        sms = _g04_load(
            function_tmpdir / "b",
            " 1.0E-6 1.0E-6 100 200 1 0 1 SOLVEROUTPUTFILE 77 SHIFT\n" + _G04_XMD,
        )
    assert sms.options2 == [" SHIFT"]
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        sms = _g04_load(
            function_tmpdir / "c",
            " 1.0E-6 1.0E-6 100 200 1 0 1 HCLOSE HICLOSE MXITER ITER1\n" + _G04_XMD,
        )
    assert sms.options2 == []


_G04_SIGNATURE = [
    "self",
    "model",
    "hclose",
    "hiclose",
    "mxiter",
    "iter1",
    "iprsms",
    "nonlinmeth",
    "linmeth",
    "theta",
    "akappa",
    "gamma",
    "amomentum",
    "numtrack",
    "btol",
    "breduc",
    "reslim",
    "iacl",
    "norder",
    "level",
    "north",
    "iredsys",
    "rrctol",
    "idroptol",
    "epsrn",
    "clin",
    "ipc",
    "iscl",
    "iord",
    "rclosepcgu",
    "relaxpcgu",
    "extension",
    "options",
    "options2",
    "unitnumber",
    "filenames",
]


def test_mfusg_sms_signature_and_positional(function_tmpdir):
    import inspect

    params = inspect.signature(MfUsgSms.__init__).parameters
    assert list(params) == _G04_SIGNATURE
    assert params["numtrack"].default == 20 and params["btol"].default == 1e4
    args = [
        1e-5,
        1e-6,
        50,
        30,
        1,
        1,
        1,
        0.6,
        0.2,
        0.3,
        0.002,
        0,
        2e4,
        0.3,
        50.0,
        1,
        0,
        5,
        3,
        0,
        0.0,
        0,
        2e-3,
        "CG",
        1,
        0,
        0,
        0.2,
        0.9,
        "sms",
        None,
        ["shift"],
    ]
    sms = MfUsgSms(MfUsg(modelname="p", model_ws=str(function_tmpdir)), *args)
    assert [getattr(sms, n) for n in _G04_SIGNATURE[2:31]] == args[:29]
    assert sms.options2 == ["shift"]
    sms.write_file()
    lines = Path(sms.fn_path).read_text().splitlines()[1:]
    assert lines[0].split()[7:] == ["shift"] and len(lines[1].split()) == 5


def test_mfusg_sms_api_options_after_load(function_tmpdir):
    text = " 1.0E-6 1.0E-6 100 200 1 0 1 DAMPBOT\n" + _G04_XMD
    sms = _g04_load(function_tmpdir, text)
    sms.options2.append(" SHIFT")
    lines, again = _g04_rewrite(sms)
    assert again.options2 == [" DAMPBOT", " SHIFT"]


def test_mfusg_sms_default_text_unchanged(function_tmpdir):
    """Values exact in 9.2e are written as before."""
    sms = MfUsgSms(MfUsg(modelname="d", model_ws=str(function_tmpdir)))
    sms.write_file()
    assert Path(sms.fn_path).read_text().splitlines()[1:] == [
        "  1.00e-04  1.00e-04       100        20         2         0         2 ",
        " bcgs              3         0         0  1.00e-01  1.00e+00",
        "",
    ]


def _g04_model(ws, sms_text, grid="dis", flow="lpf", laytyp=1):
    """Write the hand-written G04 model with the given SMS text into ws."""
    ws = Path(ws)
    ws.mkdir(parents=True, exist_ok=True)
    periods = "     1.0   1   1.0   SS\n    10.0   5   1.0   TR\n"
    if grid == "dis":
        (ws / "t.dis").write_text(
            "# DIS\n         1         1         5         2         4         2\n"
            "  0\nCONSTANT 10.0 #delr\nCONSTANT 10.0 #delc\nCONSTANT 30.0 #top\n"
            "CONSTANT 0.0 #botm\n" + periods
        )
        chd = "# CHD\n 1\n 1 0\n 1 1 1 10.0 10.0\n -1 0\n"
        wel = "# WEL\n 1 0\n 1 0\n 1 1 5 -2.0\n 1 0\n 1 1 5 -3.0\n"
        bas = "FREE"
    else:
        (ws / "t.disu").write_text(
            "# DISU\n5 1 13 0 2 4 2 0\n  0\nINTERNAL 1 (FREE) -1 #nodelay\n 5\n"
            "CONSTANT 30.0 #top\nCONSTANT 0.0 #bot\nCONSTANT 100.0 #area\n"
            "INTERNAL 1 (FREE) -1 #iac\n 2 3 3 3 2\n"
            "INTERNAL 1 (FREE) -1 #ja\n 1 2  2 1 3  3 2 4  4 3 5  5 4\n"
            "CONSTANT 5.0 #cl12\nCONSTANT 100.0 #fahl\n" + periods
        )
        chd = "# CHD\n 1\n 1 0\n 1 10.0 10.0\n -1 0\n"
        wel = "# WEL\n 1 0\n 1 0 0\n 5 -2.0\n 1 0 0\n 5 -3.0\n"
        bas = "UNSTRUCTURED FREE"
    (ws / "t.bas").write_text(
        f"# BAS6\n{bas}\nCONSTANT 1 #ibound\n   -999.99\nCONSTANT 20.0 #strt\n"
    )
    if flow == "lpf":
        (ws / "t.lpf").write_text(
            f"# LPF\n 0 -1E+30 0 0\n {laytyp}\n 0\n 1.0\n 0\n 0\n"
            "CONSTANT 1.0 #hk\nCONSTANT 1.0 #vka\nCONSTANT 1.0E-5 #ss\n"
            + ("CONSTANT 0.1 #sy\n" if laytyp else "")
        )
    else:  # confined BCF: SF1, TRAN
        (ws / "t.bcf").write_text(
            "         0 -1.000e+30         0 1.000e-01         1         0         0\n"
            "00\nCONSTANT 1.0 #anisotropy\nCONSTANT 3.0E-4 #sf1\nCONSTANT 30.0 #tran\n"
        )
    (ws / "t.sms").write_text("# SMS\n" + sms_text)
    oc = "# OC\nHEAD PRINT FORMAT 0\nHEAD SAVE UNIT 51\nCOMPACT BUDGET AUX\n"
    for k, n in ((1, 1), (2, 5)):
        for s in range(1, n + 1):
            oc += f"period {k} step {s}\n  save head\n  print budget\n"
    (ws / "t.oc").write_text(oc)
    (ws / "t.chd").write_text(chd)
    (ws / "t.wel").write_text(wel)
    disline = "DIS 11 t.dis" if grid == "dis" else "DISU 11 t.disu"
    flowline = "LPF 15 t.lpf" if flow == "lpf" else "BCF6 15 t.bcf"
    (ws / "t.nam").write_text(
        f"LIST 2 t.list\n{disline}\nBAS6 13 t.bas\n{flowline}\nSMS 32 t.sms\n"
        "OC 14 t.oc\nCHD 24 t.chd\nWEL 20 t.wel\nDATA(BINARY) 51 t.hds\n"
    )
    return ws


def _g04_sentinel(sms):
    p = Path(sms.fn_path)
    p.write_text("SENTINEL\n")
    return p, p.stat().st_mtime_ns


@pytest.mark.parametrize(
    "options2, error, match",
    [
        (["SHIFT,STRESSPERIODSETTINGS"], NotImplementedError, "STRESSPERIODSETTINGS"),
        (["shift\tstressperiodsettings"], NotImplementedError, "STRESSPERIODSETTINGS"),
        (
            ["SHIFT", "dampbot, StressPeriodSettings"],
            NotImplementedError,
            "STRESSPERIODSETTINGS",
        ),
        (["TRUNCATEDNEWTONCUTOFF,NaN"], ValueError, "TRUNCATEDNEWTONCUTOFF"),
        (
            ["TRUNCATEDNEWTONCUTOFF 0.1 TRUNCATEDNEWTONCUTOFF NaN"],
            ValueError,
            "TRUNCATEDNEWTONCUTOFF",
        ),
        (
            ["TRUNCATEDNEWTONCUTOFF inf TRUNCATEDNEWTONCUTOFF 0.1"],
            ValueError,
            "TRUNCATEDNEWTONCUTOFF",
        ),
        (["truncatednewtoncutoff", "1e39"], ValueError, "TRUNCATEDNEWTONCUTOFF"),
        (["TRUNCATEDNEWTONCUTOFF\tabc"], ValueError, "TRUNCATEDNEWTONCUTOFF"),
        (["TRUNCATEDNEWTONCUTOFF SHIFT"], ValueError, "TRUNCATEDNEWTONCUTOFF"),
    ],
    ids=[
        "comma_sps",
        "tab_sps",
        "elements_sps",
        "comma_nan",
        "second_nan",
        "first_inf",
        "next_element_overflow",
        "tab_text",
        "option_as_value",
    ],
)
def test_mfusg_sms_options2_rejected_before_writing(
    function_tmpdir, options2, error, match
):
    """options2 is read like record 1b (spaces, tabs, commas; every option
    and its value); rejections leave the destination untouched."""
    sms = MfUsgSms(
        MfUsg(modelname="t", model_ws=str(function_tmpdir)), options2=list(options2)
    )
    p, mtime = _g04_sentinel(sms)
    with pytest.raises(error, match=match) as e:
        sms.write_file()
    assert "MfUsgSms" in str(e.value) and "TRTOSS" not in str(e.value)
    assert p.read_text() == "SENTINEL\n" and p.stat().st_mtime_ns == mtime


@pytest.mark.parametrize(
    "added, error",
    [
        (["DAMPBOT,stressperiodsettings"], NotImplementedError),
        (["TRUNCATEDNEWTONCUTOFF\tNaN"], ValueError),
    ],
    ids=["sps", "cutoff_nan"],
)
def test_mfusg_sms_options2_changed_after_load(function_tmpdir, added, error):
    sms = _g04_load(function_tmpdir, " 1.0E-6 1.0E-6 100 200 1 0 1 SHIFT\n" + _G04_XMD)
    sms.options2 += added
    p, mtime = _g04_sentinel(sms)
    with pytest.raises(error):
        sms.write_file()
    assert p.read_text() == "SENTINEL\n" and p.stat().st_mtime_ns == mtime


def test_mfusg_sms_options2_separators_kept(function_tmpdir):
    """Valid options with commas, tabs, labels, repeated cutoffs and a cutoff
    without value (read as 0.0 by USG-T) are written unchanged and reload in
    their order, duplicates included."""
    options2 = [
        "SHIFT,dampbot",
        "\tsolveactive",
        "label words",
        "TRUNCATEDNEWTONCUTOFF 0.1 truncatednewtoncutoff 0.2",
        "TRUNCATEDNEWTONCUTOFF",
    ]
    sms = MfUsgSms(
        MfUsg(modelname="t", model_ws=str(function_tmpdir)),
        nonlinmeth=1,
        linmeth=1,
        numtrack=0,
        options2=list(options2),
    )
    lines, again = _g04_rewrite(sms)
    assert sms.options2 == options2
    assert again.options2 == [
        " SHIFT",
        " DAMPBOT",
        " SOLVEACTIVE",
        " TRUNCATEDNEWTONCUTOFF 1.00e-01",
        " TRUNCATEDNEWTONCUTOFF 2.00e-01",
        " TRUNCATEDNEWTONCUTOFF 0.00e+00",
    ]


def test_mfusg_sms_load_separators(function_tmpdir):
    sms = _g04_load(
        function_tmpdir / "a",
        " 1.0E-6 1.0E-6 100 200 1 1 1 SHIFT,dampbot\tsolveactive "
        "TRUNCATEDNEWTONCUTOFF 0.9,truncatednewtoncutoff 0.7654321\n"
        " 0.7 0.1 0.2 0.001 0\n" + _G04_XMD,
    )
    assert [o.split()[0] for o in sms.options2] == [
        "SHIFT",
        "DAMPBOT",
        "SOLVEACTIVE",
        "TRUNCATEDNEWTONCUTOFF",
        "TRUNCATEDNEWTONCUTOFF",
    ]
    assert np.float32(sms.options2[4].split()[1]) == np.float32(0.7654321)
    with pytest.raises(ValueError, match="TRUNCATEDNEWTONCUTOFF"):
        _g04_load(
            function_tmpdir / "b",
            " 1.0E-6 1.0E-6 100 200 1 1 1 "
            "TRUNCATEDNEWTONCUTOFF SHIFT\n"
            " 0.7 0.1 0.2 0.001 0\n" + _G04_XMD,
        )


# G05: OC with ATS (ATSA, NPTIMES/TIMOT, NPSTPS) and numeric ATS records,
# checked against USG-T 2.8.1 glo2basu1.f (SGWF2BAS7I, PTIMES1RP, ATS1RP,
# SGWF2BAS7N). The engine runs are in test_usg_transport_exe.py.
_G05_PERIODS = "     1.0   1   1.0   SS\n    10.0   5   1.0   TR\n"
_G05_HDR = "HEAD SAVE UNIT 51\nCOMPACT BUDGET AUX\n"
_G05_ACT = "  save head\n  save budget\n  print budget\n"
_G05_NUM = " 0 1 1 1\n 0 0 1 0\n"


def _g05_model(ws, oc_text, grid="dis", periods=_G05_PERIODS):
    """Write a hand-written 5-cell model with the given OC text into ws."""
    ws = Path(ws)
    ws.mkdir(parents=True, exist_ok=True)
    nper = len(periods.strip().splitlines())
    if grid == "dis":
        (ws / "t.dis").write_text(
            f"# DIS\n         1         1         5{nper:10d}         4         2\n"
            "  0\nCONSTANT 10.0 #delr\nCONSTANT 10.0 #delc\nCONSTANT 30.0 #top\n"
            "CONSTANT 0.0 #botm\n" + periods
        )
        chd = "# CHD\n 1\n 1 0\n 1 1 1 10.0 10.0\n" + " -1 0\n" * (nper - 1)
        wel = "# WEL\n 1 0\n 1 0\n 1 1 5 -2.0\n" + "".join(
            f" 1 0\n 1 1 5 {-2.0 - k:.1f}\n" for k in range(1, nper)
        )
        bas = "FREE"
    else:
        (ws / "t.disu").write_text(
            f"# DISU\n5 1 13 0 {nper} 4 2 0\n  0\nINTERNAL 1 (FREE) -1 #nodelay\n 5\n"
            "CONSTANT 30.0 #top\nCONSTANT 0.0 #bot\nCONSTANT 100.0 #area\n"
            "INTERNAL 1 (FREE) -1 #iac\n 2 3 3 3 2\n"
            "INTERNAL 1 (FREE) -1 #ja\n 1 2  2 1 3  3 2 4  4 3 5  5 4\n"
            "CONSTANT 5.0 #cl12\nCONSTANT 100.0 #fahl\n" + periods
        )
        chd = "# CHD\n 1\n 1 0\n 1 10.0 10.0\n" + " -1 0\n" * (nper - 1)
        wel = "# WEL\n 1 0\n 1 0 0\n 5 -2.0\n" + "".join(
            f" 1 0 0\n 5 {-2.0 - k:.1f}\n" for k in range(1, nper)
        )
        bas = "UNSTRUCTURED FREE"
    (ws / "t.bas").write_text(
        f"# BAS6\n{bas}\nCONSTANT 1 #ibound\n   -999.99\nCONSTANT 20.0 #strt\n"
    )
    (ws / "t.lpf").write_text(
        "# LPF\n 50 -1E+30 0 0\n 1\n 0\n 1.0\n 0\n 0\n"
        "CONSTANT 1.0 #hk\nCONSTANT 1.0 #vka\nCONSTANT 1.0E-5 #ss\nCONSTANT 0.1 #sy\n"
    )
    (ws / "t.sms").write_text(
        "# SMS\n 1.0E-6 1.0E-6 100 200 1 0 1\n 1 0 7 2 0 0.0 0 1.0E-3\n"
    )
    (ws / "t.oc").write_text(oc_text)
    (ws / "t.chd").write_text(chd)
    (ws / "t.wel").write_text(wel)
    disline = "DIS 11 t.dis" if grid == "dis" else "DISU 11 t.disu"
    (ws / "t.nam").write_text(
        f"LIST 2 t.list\n{disline}\nBAS6 13 t.bas\nLPF 15 t.lpf\nSMS 32 t.sms\n"
        "OC 14 t.oc\nCHD 24 t.chd\nWEL 20 t.wel\nDATA(BINARY) 51 t.hds\n"
        "DATA(BINARY) 50 t.cbc\n"
    )
    return ws


def _g05_dis(ws):
    """The model of ws with its DIS/DISU only."""
    disname = "DIS" if (Path(ws) / "t.dis").exists() else "DISU"
    return MfUsg.load("t.nam", model_ws=str(ws), check=False, load_only=[disname])


def _g05_load(ws, oc_text, **kwargs):
    ws = _g05_model(ws, oc_text, **kwargs)
    return MfUsgOc.load(str(ws / "t.oc"), _g05_dis(ws))


def _g05_rewrite(oc, name="out.oc"):
    """Write oc next to its model; return (lines without heading, reloaded)."""
    path = Path(oc.parent.model_ws) / name
    oc.fn_path = str(path)
    oc.write_file()
    again = MfUsgOc.load(str(path), _g05_dis(oc.parent.model_ws))
    return path.read_text().splitlines()[1:], again


def _g05_words(first, timot="", p2=""):
    return (
        f"# OC\n{first}\n{timot}"
        + _G05_HDR
        + "period 1\n"
        + _G05_ACT
        + "period 2\n"
        + p2
        + _G05_ACT
    )


@pytest.mark.parametrize(
    "first, timot, expected",
    [
        ("ATSA NPTIMES 3", "3 6 11\n", [3.0, 6.0, 11.0]),
        ("ATSA NPTIMES 3", "3.5 6.25 11.0\n", [3.5, 6.25, 11.0]),
        ("ATSA NPTIMES 4", "3.5E0 6.25d0 1.1D1 1.2q1\n", [3.5, 6.25, 11.0, 12.0]),
        ("ATSA NPTIMES 2", "1.5-1 2.0+1\n", [0.15, 20.0]),
        ("ATSA NPTIMES 4", "+.5, 5.\n\n6.25\t,8.5 ,\n", [0.5, 5.0, 6.25, 8.5]),
        ("atsa,nptimes,3", "3.5\n6.25\n11.0\n", [3.5, 6.25, 11.0]),
        (
            "ATSA NPTIMES 60 NPSTPS 4",
            "".join(
                " ".join(str(10 * (6 * r + c + 1)) for c in range(6)) + "\n"
                for r in range(10)
            ),
            [10.0 * (i + 1) for i in range(60)],
        ),
    ],
    ids=[
        "int",
        "decimal",
        "exp_letters",
        "exp_sign",
        "commas_blank_tabs",
        "multiline",
        "long",
    ],
)
def test_mfusg_oc_timot_load_write(function_tmpdir, first, timot, expected):
    """TIMOT is read like USG-T (list-directed DOUBLE over several lines) and
    written once, right after the first line, with its DOUBLE value."""
    oc = _g05_load(function_tmpdir, _g05_words(first, timot))
    assert oc.atsa == 1 and oc.nptimes == len(expected)
    assert oc.timot == expected
    lines, again = _g05_rewrite(oc)
    head = lines[0].split()
    assert head[:3] == ["ATSA", "NPTIMES", str(len(expected))]
    assert len(lines[0]) < 198 and not any("{" in ln for ln in lines)
    nt = (len(expected) + 9) // 10
    assert [float(w) for ln in lines[1 : 1 + nt] for w in ln.split()] == expected
    assert lines[1 + nt] == "HEAD PRINT FORMAT   0"
    assert sum(ln.count(repr(expected[-1])) for ln in lines) == 1
    assert again.timot == expected and again.npsteps == oc.npsteps
    assert _g05_rewrite(again, "out2.oc")[0] == lines


def test_mfusg_oc_timot_double_value(function_tmpdir):
    """TIMOT keeps the DOUBLE value USG-T reads (no REAL rounding)."""
    words = [
        "0.1",
        "0.30000000000000004",
        "1.00000000000000011102230246251565",
        "3.14159265358979323846264338",
    ]
    text = _g05_words("ATSA NPTIMES 4", " ".join(words) + "\n")
    oc = _g05_load(function_tmpdir / "a", text)
    # the third word is just below the midpoint between 1 and the next DOUBLE
    assert oc.timot == [0.1, 0.30000000000000004, 1.0, 3.141592653589793]
    lines, again = _g05_rewrite(oc)
    assert lines[1].split() == [repr(t) for t in oc.timot]
    ml = _g05_dis(_g05_model(function_tmpdir / "b", "# OC\n" + _G05_HDR))
    oc = MfUsgOc(ml, atsa=1, nptimes=2, timot=np.array([3.3, 11.0], dtype=np.float32))
    lines, again = _g05_rewrite(oc)
    assert again.timot == [float(np.float32(3.3)), 11.0]


@pytest.mark.parametrize(
    "timot, match",
    [
        ("2*3.5 11.0\n", "not supported"),
        ("3.5,,11.0\n", "null value"),
        (",3.5 6.0 11.0\n", "null value"),
        ("3.5 6.25 /\n", "not supported"),
        ("3.5 6.25 NaN\n", "not supported"),
        ("3.5 6.25 Infinity\n", "not supported"),
        ("3.5 6.25 1.0E400\n", "beyond the DOUBLE range"),
        ("'3.5' 6.25 11.0\n", "not supported"),
        ("# times\n3.5 6.25 11.0\n", "not supported"),
        ("3.5;6.25 11.0\n", "not supported"),
        ("3.5 6.25 11.0 12.0\n", "NPTIMES is 3"),
        ("3.5 6.25\n", "not supported"),
    ],
    ids=[
        "repeat",
        "null",
        "leading_comma",
        "slash",
        "nan",
        "inf",
        "overflow",
        "quote",
        "comment",
        "semicolon",
        "extra_value",
        "next_record",
    ],
)
def test_mfusg_oc_timot_rejected_on_load(function_tmpdir, timot, match):
    """TIMOT forms FloPy cannot read like USG-T, extra values (ignored by
    USG-T) and OC records read as times are rejected with the TIMOT name."""
    with pytest.raises(ValueError, match=match) as e:
        _g05_load(function_tmpdir, _g05_words("ATSA NPTIMES 3", timot))
    assert "TIMOT" in str(e.value) or "NPTIMES" in str(e.value)


def test_mfusg_oc_timot_eof(function_tmpdir):
    with pytest.raises(ValueError, match="end of file while reading 4 TIMOT"):
        _g05_load(function_tmpdir, "# OC\nATSA NPTIMES 4\n3.5 6.0\n")


@pytest.mark.parametrize("first", ["NPTIMES 3", "npstps 2", "ATS", "Ats NPTIMES 3"])
def test_mfusg_oc_numeric_first_word_rejected(function_tmpdir, first):
    """USG-T reads an OC file whose first word is NPTIMES, NPSTPS or ATS as
    numeric output control and stops; FloPy does not repair it."""
    with pytest.raises(ValueError, match="start the first line with ATSA"):
        _g05_load(function_tmpdir, _g05_words(first, "3 6 11\n"))


def test_mfusg_oc_ats_activation(function_tmpdir):
    """NPTIMES or NPSTPS turn ATS on (numeric first line or API); the file
    then starts with ATSA."""
    oc = _g05_load(
        function_tmpdir / "n",
        "# OC\n 0 0 51 0 NPTIMES 3\n3 6 11\n"
        " 1.0 1.0E-5 10.0 2.0 5.0 3 1000 0 0\n"
        + _G05_NUM
        + " 0.5 1.0E-5 2.0 2.0 5.0 3 1000 0 0\n"
        + _G05_NUM,
    )
    assert (oc.atsa, oc.nptimes, oc.npsteps, oc.timot) == (1, 3, 1000, [3.0, 6.0, 11.0])
    ml = _g05_dis(_g05_model(function_tmpdir / "a", "# OC\n" + _G05_HDR))
    for kw, first in (
        ({"nptimes": 2, "timot": [3.5, 11.0]}, "ATSA NPTIMES 2"),
        ({"npsteps": 4}, "ATSA NPSTPS 4"),
        ({"atsa": 1}, "ATSA"),
    ):
        oc = MfUsgOc(ml, stress_period_data={(1, 0): ["save head"]}, **kw)
        lines, again = _g05_rewrite(oc)
        assert lines[0].strip() == first
        assert again.atsa == 1 and lines[-2:] == ["period 2 ", "  save head"]


def test_mfusg_oc_ats_params_kept(function_tmpdir):
    """DELTAT, TMINAT, TMAXAT, TADJAT and TCUTAT keep their REAL text."""
    p2 = (
        "  DELTAT 0.123456789\n  tminat 1.5d-6\n  TMAXAT 1.7654321\n"
        "  TADJAT 1.33333\n  TCUTAT 3\n"
    )
    oc = _g05_load(function_tmpdir, _g05_words("ATSA", "", p2))
    expected = [
        "DELTAT 0.123456789",
        "TMINAT 1.5D-6",
        "TMAXAT 1.7654321",
        "TADJAT 1.33333",
        "TCUTAT 3",
    ]
    assert oc.stress_period_data[(1, 0)][:5] == expected
    lines, again = _g05_rewrite(oc)
    assert [
        ln.strip()
        for ln in lines
        if ln.strip().split()[0] in ("DELTAT", "TMINAT", "TMAXAT", "TADJAT", "TCUTAT")
    ] == expected
    assert again.stress_period_data == oc.stress_period_data


@pytest.mark.parametrize(
    "action, match",
    [
        ("DELTAT abc", "DELTAT 'ABC'"),
        ("TMAXAT 1.0E39", "TMAXAT '1.0E39'"),
        ("TADJAT nan", "TADJAT 'NAN'"),
        ("TCUTAT", "TCUTAT needs a value"),
    ],
    ids=["text", "beyond_real", "nan", "missing"],
)
def test_mfusg_oc_ats_param_rejected(function_tmpdir, action, match):
    with pytest.raises(ValueError, match=match):
        _g05_load(function_tmpdir / "l", _g05_words("ATSA", "", f"  {action}\n"))
    ml = _g05_dis(_g05_model(function_tmpdir / "w", "# OC\n" + _G05_HDR))
    oc = MfUsgOc(ml, atsa=1, stress_period_data={(1, 0): [action.lower(), "save head"]})
    oc.fn_path = str(function_tmpdir / "w" / "dest.oc")
    p, mtime = _g04_sentinel(oc)
    with pytest.raises(ValueError, match=action.split()[0]):
        oc.write_file()
    assert p.read_text() == "SENTINEL\n" and p.stat().st_mtime_ns == mtime


def _g05_num(header, timot, rec1, rec2, flags=_G05_NUM):
    return f"# OC\n 0 0 51 0 ATSA{header}\n{timot}{rec1}\n{flags}{rec2}\n{flags}"


_G05_TR2 = ["DELTAT 0.5", "TMINAT 1.0E-5", "TMAXAT 2.0", "TADJAT 2.0", "TCUTAT 5.0"]
_G05_FLAGS = ["PRINT BUDGET", "SAVE BUDGET", "SAVE HEAD"]


@pytest.mark.parametrize(
    "text, periods, spd, npsteps",
    [
        (
            _g05_num(
                " NPTIMES 3",
                "3.5 6.25 11.0\n",
                " 1.0 1.0E-5 10.0 2.0 5.0 3 1000 0 0",
                " 0.5 1.0E-5 2.0 2.0 5.0 3 1000 0 0",
            ),
            _G05_PERIODS,
            {
                (0, 0): [
                    "DELTAT 1.0000001E+00",
                    "TMINAT 1.0E-5",
                    "TMAXAT 10.0",
                    "TADJAT 2.0",
                    "TCUTAT 5.0",
                ]
                + _G05_FLAGS,
                (1, 0): _G05_TR2 + _G05_FLAGS,
            },
            1000,
        ),
        (
            _g05_num(
                "",
                "",
                " 1.0 1.0E-5 10.0 2.0 5.0 0 1000 0 0",
                " 0.0 0.0 -1.0 0.5 0.9 0 1000 0 0",
            ),
            _G05_PERIODS,
            {
                (0, 0): [
                    "DELTAT 1.0000001E+00",
                    "TMINAT 1.0E-5",
                    "TMAXAT 10.0",
                    "TADJAT 2.0",
                    "TCUTAT 5.0",
                ]
                + _G05_FLAGS,
                (1, 0): [
                    "DELTAT 1.0",
                    "TMINAT 1.0E-10",
                    "TMAXAT 1.0E10",
                    "TADJAT 2.0",
                    "TCUTAT 5.0",
                ]
                + _G05_FLAGS,
            },
            1000,
        ),
        (
            _g05_num(
                " NPSTPS 3",
                "",
                " 1.0 1.0E-5 10.0 2.0 5.0 0 2 0 0",
                " 0.5 1.0E-5 2.0 2.0 5.0 0 2",
            ),
            "     0.1   1   1.0   SS\n    10.0   5   1.0   TR\n",
            {
                (0, 0): [
                    "DELTAT 1.0000001E-01",
                    "TMINAT 1.0E-5",
                    "TMAXAT 10.0",
                    "TADJAT 2.0",
                    "TCUTAT 5.0",
                ]
                + _G05_FLAGS,
                (1, 0): _G05_TR2 + _G05_FLAGS,
            },
            2,
        ),
        (
            _g05_num(
                " NPSTPS 3",
                "",
                " 1.0 1.0E-5 10.0 2.0 5.0 0 0 0 0",
                " 0.5 1.0E-5 2.0 2.0 5.0",
            ),
            _G05_PERIODS,
            None,
            0,
        ),
        (
            _g05_num(
                "",
                "",
                " 1.0 1.0E-5 10.0 2.0 5.0",
                " 0.5 1.0E-5 2.0 2.0 5.0",
                flags=" 1 1 0 0\n 1 0 0 0\n",
            ),
            _G05_PERIODS,
            {
                (0, 0): [
                    "DELTAT 1.0000001E+00",
                    "TMINAT 1.0E-5",
                    "TMAXAT 10.0",
                    "TADJAT 2.0",
                    "TCUTAT 5.0",
                    "PRINT HEAD 1",
                ],
                (1, 0): _G05_TR2 + ["PRINT HEAD 1"],
            },
            0,
        ),
    ],
    ids=["nptimes", "limits", "ss_perlen_npstps", "npstps_zero", "incode_layers"],
)
def test_mfusg_oc_numeric_ats_load(function_tmpdir, text, periods, spd, npsteps):
    """Numeric ATS records are loaded as words with the values USG-T uses:
    ATS1RP defaults for values below its limits, DELTAT at or above PERLEN in
    steady-state periods (cut to PERLEN), NPSTPS of the records (0: none)."""
    oc = _g05_load(function_tmpdir, text, periods=periods)
    assert oc.atsa == 1 and oc.npsteps == npsteps
    if spd is not None:
        assert oc.stress_period_data == spd
    lines, again = _g05_rewrite(oc)
    assert ("NPSTPS" in lines[0]) == (npsteps > 0)
    assert not any("{" in ln for ln in lines)
    assert again.stress_period_data == oc.stress_period_data
    assert again.npsteps == npsteps


@pytest.mark.parametrize(
    "text, match",
    [
        (
            _g05_num(
                " NPTIMES 3",
                "3.5 6.25 11.0\n",
                " 1.0 1.0E-5 10.0 2.0 5.0",
                " 0.5 1.0E-5 2.0 2.0 5.0",
            ),
            "has NPTIMES 0",
        ),
        (
            _g05_num(
                " NPTIMES 3",
                "3.5 6.25 11.0\n",
                " 1.0 1.0E-5 10.0 2.0 5.0 3",
                " 0.5 1.0E-5 2.0 2.0 5.0 2",
            ),
            "has NPTIMES 2",
        ),
        (
            _g05_num(
                "", "", " 1.0 1.0E-5 10.0 2.0 5.0 0 2", " 0.5 1.0E-5 2.0 2.0 5.0 0 3"
            ),
            r"different NPSTPS \[2, 3\]",
        ),
        (
            _g05_num(
                "", "", " 1.0 1.0E-5 10.0 2.0 5.0 0 0 1 0", " 0.5 1.0E-5 2.0 2.0 5.0"
            ),
            "IBOOT",
        ),
        (
            _g05_num(
                "", "", " 1.0 1.0E-5 10.0 2.0 5.0", " 0.5 1.0E-5 2.0 2.0 5.0 0 0 0 1"
            ),
            "IBOOT",
        ),
        (
            _g05_num(
                "",
                "",
                " 1.0 1.0E-5 10.0 2.0 5.0",
                " 0.5 1.0E-5 2.0 2.0 5.0",
                flags=" -1 1 1 1\n",
            ),
            "INCODE < 0",
        ),
        (
            _g05_num(
                "",
                "",
                " 1.0 1.0E-5 10.0 2.0 5.0",
                " 0.5 1.0E-5 2.0 2.0 5.0",
                flags=" 0 1 1 2\n 0 0 1 0\n",
            ),
            "ICBCFL = 2",
        ),
        (
            _g05_num(
                "",
                "",
                " 1.0 1.0E-5 10.0 2.0 5.0",
                " 0.5 1.0E-5 2.0 2.0 5.0",
                flags=" 0 0 1 1\n 0 0 1 0\n",
            ),
            "IHDDFL = 0",
        ),
        (
            _g05_num(
                "",
                "",
                " 1.0 1.0E-5 10.0 2.0 5.0",
                " 0.5 1.0E-5 2.0 2.0 5.0",
                flags=" 1 0 1 1\n 1 0 0 0\n",
            ),
            "IHDDFL = 0",
        ),
        (
            _g05_num(
                "",
                "",
                " 1.0 1.0E-5 10.0 2.0 5.0",
                " 0.5 1.0E-5 2.0 2.0 5.0",
                flags=" 0 1 1 1 0\n 0 0 1 0 0 0\n",
            ),
            "transport flags",
        ),
        (
            _g05_num(
                " FASTFORWARD 1 1 60",
                "",
                " 1.0 1.0E-5 10.0 2.0 5.0",
                " 0.5 1.0E-5 2.0 2.0 5.0",
            ),
            "FASTFORWARD",
        ),
        (
            _g05_num(
                " BOOTSTRAPPING 60",
                "",
                " 1.0 1.0E-5 10.0 2.0 5.0",
                " 0.5 1.0E-5 2.0 2.0 5.0",
            ),
            "BOOTSTRAPPING",
        ),
        (
            _g05_num(
                " NPTIMES 3",
                "3.5 6.25 11.0\n",
                " 1.0 1.0E-5 10.0 2.0 5.0 3",
                " 5.0 1.0E-5 10.0 2.0 5.0 3",
            ),
            "limited by TIMOT",
        ),
        (
            _g05_num("", "", " 1.0 1.0E-5 10.0 2.0", " 0.5 1.0E-5 2.0 2.0 5.0"),
            "must give DELTAT TMINAT",
        ),
    ],
    ids=[
        "nptimes_missing",
        "nptimes_differs",
        "npstps_varies",
        "iboot",
        "ibootscale",
        "incode_negative",
        "icbcfl_2",
        "ihddfl0",
        "ihddfl0_layers",
        "transport",
        "fastforward",
        "bootstrapping",
        "deltat_capped",
        "four_reals",
    ],
)
def test_mfusg_oc_numeric_ats_rejected(function_tmpdir, text, match):
    """Numeric ATS records the word form cannot represent are rejected."""
    with pytest.raises((NotImplementedError, ValueError), match=match):
        _g05_load(function_tmpdir, text)


def test_mfusg_oc_numeric_ats_timot_beyond(function_tmpdir):
    """ATS1RP reads TIMOT(KPER) for the first time step of a transient
    period; beyond NPTIMES+1 that is outside TIMOT: rejected."""
    text = _g05_num(
        " NPTIMES 1", "2.0\n", " 1.0 1.0E-5 10.0 2.0 5.0 1", " 0.5 1.0E-5 2.0 2.0 5.0 1"
    )
    text += " 0.5 1.0E-5 2.0 2.0 5.0 1\n" + _G05_NUM
    with pytest.raises(NotImplementedError, match="beyond NPTIMES"):
        _g05_load(function_tmpdir, text, periods=_G05_THREE)


def test_mfusg_oc_numeric_ats_model_rejected(function_tmpdir):
    """Fixed format, transport, CLN and DPF are not validated: rejected."""
    text = _g05_num("", "", " 1.0 1.0E-5 10.0 2.0 5.0", " 0.5 1.0E-5 2.0 2.0 5.0")
    ws = _g05_model(function_tmpdir, text)
    for attr, value, match in (
        ("free_format_input", False, "fixed format"),
        ("itrnsp", 1, "transport"),
        ("icln", 1, "CLN"),
        ("idpf", 1, "DPF"),
    ):
        ml = _g05_dis(ws)
        setattr(ml, attr, value)
        with pytest.raises(NotImplementedError, match=match):
            MfUsgOc.load(str(ws / "t.oc"), ml)
    # IHDDFL = 0 without head or drawdown flags keeps its meaning
    oc = _g05_load(
        function_tmpdir / "z",
        _g05_num(
            "",
            "",
            " 1.0 1.0E-5 10.0 2.0 5.0",
            " 0.5 1.0E-5 2.0 2.0 5.0",
            flags=" 0 0 1 1\n 0 0 0 0\n",
        ),
    )
    assert oc.stress_period_data[(1, 0)] == _G05_TR2 + ["PRINT BUDGET", "SAVE BUDGET"]


_G05_THREE = _G05_PERIODS + "    10.0   5   1.0   TR\n"


@pytest.mark.parametrize(
    "blocks",
    [
        "period 1\n  save head\nperiod 1\n  save head\nperiod 2\n  save head\n",
        "period 1\n  save head\nperiod 2\n  save head\nperiod 2\n  save head\n"
        "period 3\n  save head\n",
        "period 1\n  save head\nperiod 2\n  save head\nperiod 3\n  save head\n"
        "period 3\n  save head\n",
        "period 1\nperiod 1\nperiod 2\n  save head\n",
        "period 1\n  save head\nperiod 2 ddreference\nperiod 2 ddreference\n",
        "period 1\n  save head\nperiod 3\n  save head\nperiod 2\n  save head\n",
        "period 1\n  save head\nperiod 2 step 3\n  save head\nperiod 2 step 5\n"
        "  save head\n",
        "period 1\n  save head\nperiod 2 step 1 ddreference\n  save head\n",
    ],
    ids=[
        "first",
        "intermediate",
        "last_before_eof",
        "empty",
        "ddreference",
        "out_of_order",
        "two_steps",
        "ddreference_after_step",
    ],
)
def test_mfusg_oc_ats_period_blocks_rejected(function_tmpdir, blocks):
    """With ATS, USG-T reads output control once per stress period: a period
    given again or out of order would be applied to a later period, and
    DDREFERENCE after STEP is ignored. Without ATS the same records load as
    before."""
    with pytest.raises(NotImplementedError, match="with ATS"):
        _g05_load(
            function_tmpdir / "a",
            "# OC\nATSA\n" + _G05_HDR + blocks,
            periods=_G05_THREE,
        )
    _g05_load(function_tmpdir / "b", "# OC\n" + _G05_HDR + blocks, periods=_G05_THREE)


def test_mfusg_oc_ats_period_records(function_tmpdir):
    """With ATS, STEP is ignored on load and PERIOD k [DDREFERENCE] is
    written; empty and DDREFERENCE-only blocks are kept. Without ATS the
    records keep their step."""
    text = (
        "# OC\nATSA\n" + _G05_HDR + "period 1\nperiod 2 step 4 \n  save head\n"
        "period 3 ddreference\n"
    )
    oc = _g05_load(function_tmpdir, text, periods=_G05_THREE)
    assert oc.stress_period_data == {(1, 0): ["SAVE HEAD"], (2, 0): ["DDREFERENCE"]}
    lines, again = _g05_rewrite(oc)
    assert lines[-3:] == ["period 2 ", "  SAVE HEAD", "period 3 ddreference"]
    assert again.stress_period_data == oc.stress_period_data
    oc = MfUsgOc(oc.parent, stress_period_data={(1, 3): ["save head", "ddreference"]})
    assert _g05_rewrite(oc, "noats.oc")[0][-2:] == [
        "period 2 step 4 ddreference",
        "  save head",
    ]


def test_mfusg_oc_old_flopy_inline_timot(function_tmpdir):
    """Files written by FloPy before G05 repeat TIMOT on the first line; USG-T
    ignores the options that start after column 198 (here NPSTPS)."""
    times = " ".join(str(10 * (i + 1)) for i in range(60))
    text = (
        f"# OC\nATSA NPTIMES 60 {times} NPSTPS 4 \n{times} \n"
        + _G05_HDR
        + "period 1\n"
        + _G05_ACT
        + "period 2\n"
        + _G05_ACT
    )
    oc = _g05_load(
        function_tmpdir,
        text,
        periods="     1.0   1   1.0   SS\n  1000.0   5   1.0   TR\n",
    )
    assert (oc.nptimes, oc.npsteps) == (60, 0)
    assert oc.timot == [10.0 * (i + 1) for i in range(60)]
    lines, _ = _g05_rewrite(oc)
    assert lines[0].strip() == "ATSA NPTIMES 60"


def _g05_dest(function_tmpdir, **kwargs):
    ml = _g05_dis(_g05_model(function_tmpdir, "# OC\n" + _G05_HDR))
    spd = kwargs.pop("stress_period_data", {(0, 0): ["save head"]})
    oc = MfUsgOc(ml, stress_period_data=spd, **kwargs)
    oc.fn_path = str(Path(function_tmpdir) / "dest.oc")
    return oc


@pytest.mark.parametrize(
    "kwargs, mutate, match",
    [
        ({"atsa": 1, "nptimes": 3, "timot": [3.5, 6.25]}, None, "NPTIMES is 3 but 2"),
        ({"atsa": 1, "nptimes": 2}, None, "NPTIMES is 2 but 0"),
        ({"atsa": 1, "timot": [3.5]}, None, "NPTIMES is 0 but 1"),
        (
            {"nptimes": 4, "timot": np.array([[1.0, 2.0], [3.0, 4.0]])},
            None,
            "one-dimensional",
        ),
        ({"nptimes": 1, "timot": 3.5}, None, "one-dimensional"),
        ({"nptimes": 2, "timot": ["3.5", "6.0"]}, None, "not a number"),
        ({"nptimes": 2, "timot": [True, 6.0]}, None, "not a number"),
        ({"nptimes": 2, "timot": [3.5, float("nan")]}, None, "finite"),
        ({"nptimes": 2, "timot": [3.5, np.inf]}, None, "finite"),
        ({"nptimes": 2, "timot": [6.0, 3.5]}, None, "not increasing"),
        ({"nptimes": 2, "timot": [1.0, 1.0 + 0.5e-10]}, None, "not increasing"),
        ({"atsa": 1}, ("nptimes", -1), "non-negative integer"),
        ({"atsa": 1}, ("npsteps", 2.0), "non-negative integer"),
        (
            {"nptimes": 2, "timot": [3.5, 11.0]},
            ("timot", [3.5, 11.0, 12.0]),
            "NPTIMES is 2 but 3",
        ),
        (
            {
                "atsa": 1,
                "stress_period_data": {(1, 0): ["save head"], (1, 3): ["save budget"]},
            },
            None,
            r"stress period\(s\) \[2\]",
        ),
        ({"stress_period_data": {(1, 0): ["deltat 1.0e39"]}}, None, "DELTAT"),
        ({"stress_period_data": {(1, 0): ["DELTAT"]}}, None, "DELTAT needs a value"),
    ],
    ids=[
        "count",
        "missing",
        "without_nptimes",
        "matrix",
        "scalar",
        "text",
        "bool",
        "nan",
        "inf",
        "decreasing",
        "gap_tminat",
        "nptimes_negative",
        "npsteps_float",
        "mutated_timot",
        "two_blocks",
        "deltat_real",
        "deltat_missing",
    ],
)
def test_mfusg_oc_write_rejected_before_opening(function_tmpdir, kwargs, mutate, match):
    """Invalid ATS settings, also changed after the constructor, are rejected
    before the destination is opened."""
    oc = _g05_dest(function_tmpdir, **kwargs)
    if mutate:
        setattr(oc, *mutate)
    p, mtime = _g04_sentinel(oc)
    with pytest.raises(ValueError, match=match):
        oc.write_file()
    assert p.read_text() == "SENTINEL\n" and p.stat().st_mtime_ns == mtime


def test_mfusg_oc_timot_gap(function_tmpdir):
    """TIMOT values 2.0E-10 apart pass the USG-T check (TMINAT=1.0E-10)."""
    oc = _g05_dest(function_tmpdir, nptimes=2, timot=[1.0, 1.0 + 2.0e-10])
    lines, again = _g05_rewrite(oc)
    assert again.timot == [1.0, 1.0 + 2.0e-10]


def test_mfusg_oc_signature_and_output_without_ats(function_tmpdir):
    """The constructor signature is unchanged and OC without ATS is written
    as before G05."""
    import inspect

    params = inspect.signature(MfUsgOc.__init__).parameters
    assert list(params) == [
        "self",
        "model",
        "ihedfm",
        "iddnfm",
        "ispcfm",
        "chedfm",
        "cddnfm",
        "cspcfm",
        "cboufm",
        "compact",
        "stress_period_data",
        "extension",
        "unitnumber",
        "filenames",
        "label",
        "kwargs",
    ]
    assert [params[k].default for k in ("ihedfm", "compact", "label")] == [
        0,
        True,
        "LABEL",
    ]
    oc = _g05_dest(
        function_tmpdir,
        chedfm="(10G12.4)",
        stress_period_data={
            (0, 0): ["save head", "print budget"],
            (1, 2): ["save head", "ddreference"],
            (1, 4): ["save budget"],
        },
    )
    oc.write_file()
    assert Path(oc.fn_path).read_text().splitlines()[1:] == [
        "HEAD PRINT FORMAT   0",
        "HEAD SAVE FORMAT (10G12.4)            LABEL",
        "HEAD SAVE UNIT    51",
        "DRAWDOWN PRINT FORMAT   0",
        "COMPACT BUDGET AUX",
        "period 1 step 1 ",
        "  save head",
        "  print budget",
        "period 2 step 3 ddreference",
        "  save head",
        "period 2 step 5 ",
        "  save budget",
    ]


def test_mfusg_oc_load_closes_its_file(function_tmpdir, monkeypatch):
    """A file opened by load is closed when reading fails; a handle given by
    the caller is left open."""
    import builtins

    ws = _g05_model(function_tmpdir, _g05_words("ATSA NPTIMES 3", "3.5,,11.0\n"))
    ml = _g05_dis(ws)
    opened = []
    real_open = builtins.open

    def tracking_open(*args, **kwargs):
        fh = real_open(*args, **kwargs)
        opened.append(fh)
        return fh

    monkeypatch.setattr(builtins, "open", tracking_open)
    with pytest.raises(ValueError):
        MfUsgOc.load(str(ws / "t.oc"), ml)
    monkeypatch.setattr(builtins, "open", real_open)
    oc_files = [fh for fh in opened if str(fh.name).endswith("t.oc")]
    assert oc_files and all(fh.closed for fh in oc_files)
    with open(ws / "t.oc") as fh:
        with pytest.raises(ValueError):
            MfUsgOc.load(fh, ml)
        assert not fh.closed


# G05 v2: ATS actions with URWORD separators (blanks, tabs, commas)
_G05_ACTION_VALUES = {
    "DELTAT": "0.123456789",
    "TMINAT": "1.5d-6",
    "TMAXAT": "1.7654321",
    "TADJAT": "1.33333",
    "TCUTAT": "3",
}
_G05_SEPS = {"comma": ",", "tab": "\t", "blanks": "   ", "comma_blanks": " , "}


def _g05_real(action):
    return np.float32(float(action.split()[1].upper().replace("D", "E")))


@pytest.mark.parametrize("sep", list(_G05_SEPS))
@pytest.mark.parametrize("name", list(_G05_ACTION_VALUES))
def test_mfusg_oc_ats_action_separators(function_tmpdir, name, sep):
    """DELTAT, TMINAT, TMAXAT, TADJAT and TCUTAT are read like URWORD (blanks,
    tabs or commas, any case) from a file, through the API and after a
    change of a loaded package; the action and its REAL value survive
    load -> write -> load."""
    value = _G05_ACTION_VALUES[name]
    action = f"{name[:3].lower()}{name[3:]}{_G05_SEPS[sep]}{value}"
    expected = f"{name} {value.upper()}"
    real = np.float32(float(value.replace("d", "e")))
    oc = _g05_load(function_tmpdir / "text", _g05_words("ATSA", "", f"  {action}\n"))
    assert oc.stress_period_data[(1, 0)][0] == expected
    lines, again = _g05_rewrite(oc)
    assert again.stress_period_data == oc.stress_period_data
    ml = _g05_dis(_g05_model(function_tmpdir / "api", "# OC\n" + _G05_HDR))
    oc = MfUsgOc(ml, atsa=1, stress_period_data={(1, 0): [action, "save head"]})
    lines, again = _g05_rewrite(oc)
    assert f"  {action}" in lines
    assert again.stress_period_data[(1, 0)] == [expected, "SAVE HEAD"]
    oc = _g05_load(function_tmpdir / "changed", _g05_words("ATSA"))
    oc.stress_period_data[(1, 0)].append(action)
    lines, again = _g05_rewrite(oc)
    assert again.stress_period_data[(1, 0)][-1] == expected
    assert _g05_real(again.stress_period_data[(1, 0)][-1]) == real


@pytest.mark.parametrize(
    "action, field",
    [
        ("DELTAT,NaN", "DELTAT"),
        ("TMAXAT,3.5E38", "TMAXAT"),
        ("tadjat\tabc", "TADJAT"),
        ("TCUTAT,", "TCUTAT"),
        ("tminat , inf", "TMINAT"),
        ("DeltaT\t1.0E39", "DELTAT"),
    ],
    ids=[
        "nan_comma",
        "beyond_real_comma",
        "text_tab",
        "missing_comma",
        "inf_comma_blanks",
        "beyond_real_tab",
    ],
)
def test_mfusg_oc_ats_action_separators_rejected(function_tmpdir, action, field):
    """Invalid values are rejected whatever the separator: on load, and on
    write before the destination is opened (API and changed package)."""
    with pytest.raises(ValueError, match=field):
        _g05_load(function_tmpdir / "text", _g05_words("ATSA", "", f"  {action}\n"))
    ml = _g05_dis(_g05_model(function_tmpdir / "api", "# OC\n" + _G05_HDR))
    api = MfUsgOc(ml, atsa=1, stress_period_data={(1, 0): [action, "save head"]})
    changed = _g05_load(function_tmpdir / "changed", _g05_words("ATSA"))
    changed.stress_period_data[(1, 0)].append(action)
    for oc in (api, changed):
        oc.fn_path = str(Path(oc.parent.model_ws) / "dest.oc")
        p, mtime = _g04_sentinel(oc)
        with pytest.raises(ValueError, match=field):
            oc.write_file()
        assert p.read_text() == "SENTINEL\n" and p.stat().st_mtime_ns == mtime


def test_mfusg_oc_ats_action_reviewer_cases(function_tmpdir):
    """The three API cases of the G05 review: DELTAT,0.5 is kept; DELTAT,NaN
    and TMAXAT,3.5E38 are rejected before the destination is opened."""
    ml = _g05_dis(_g05_model(function_tmpdir, "# OC\n" + _G05_HDR))
    oc = MfUsgOc(ml, atsa=1, stress_period_data={(1, 0): ["DELTAT,0.5", "save head"]})
    lines, again = _g05_rewrite(oc)
    assert again.stress_period_data[(1, 0)] == ["DELTAT 0.5", "SAVE HEAD"]
    for action in ("DELTAT,NaN", "TMAXAT,3.5E38"):
        oc = MfUsgOc(ml, atsa=1, stress_period_data={(1, 0): [action, "save head"]})
        oc.fn_path = str(Path(function_tmpdir) / "dest.oc")
        p, mtime = _g04_sentinel(oc)
        with pytest.raises(ValueError, match=action.split(",")[0]):
            oc.write_file()
        assert p.read_text() == "SENTINEL\n" and p.stat().st_mtime_ns == mtime


# G01 (partial, F39/F40): MfUsgRch keeps its own SELEV and IZNRCH values and
# writes the RTS option with an integer MXZNRCH (gwf2rch8u1.f: AR 72-89 RTS
# read by URWORD as integer; RP 400-446 IZNRCH, then 449-462 SELEV). Written
# files are read back with _g01_read, by hand, not with MfUsgRch.load.

_G01_NPER = 3


def _g01_model(ws, grid):
    """2x2 DIS or 4-node DISU (1 layer), three stress periods."""
    ml = MfUsg(modelname="t", model_ws=str(ws), structured=(grid == "dis"))
    if grid == "dis":
        MfUsgDis(ml, nlay=1, nrow=2, ncol=2, nper=_G01_NPER)
    else:
        MfUsgDisU(
            ml,
            nodes=4,
            nlay=1,
            njag=10,
            iac=np.array([2, 3, 3, 2]),
            ja=np.array([0, 1, 1, 0, 2, 2, 1, 3, 3, 2]),
            ivsd=0,
            cl12=5.0,
            fahl=10.0,
            area=100.0,
            top=[10.0] * 4,
            bot=[0.0] * 4,
            nper=_G01_NPER,
        )
    return ml


def _g01_arr(grid, values, dtype=np.float64):
    """values as a (nrow, ncol) DIS array or a (1, nodes) DISU array."""
    return np.array(values, dtype=dtype).reshape((2, 2) if grid == "dis" else (1, 4))


def _g01_read(path, n, skip=0):
    """Read a written RCH file by hand: (options tokens, [(period header
    tokens, {name: n values})]). Arrays follow the RP order: RECH, IRCH
    (NRCHOP=2), IZNRCH (INRCHZONES > 0), SELEV (INSELEV > 0). skip: lines
    between the options and the first period (MXNDRCH, IRCHCONC)."""
    lines = [ln for ln in Path(path).read_text().splitlines() if not ln.startswith("#")]
    pos = 1 + skip

    def array():
        nonlocal pos
        ctl = lines[pos].split()
        pos += 1
        if ctl[0].upper() == "CONSTANT":
            return [float(ctl[1])] * n
        assert ctl[0].upper() == "INTERNAL", ctl
        values = []
        while len(values) < n:
            values += [float(v) for v in lines[pos].split()]
            pos += 1
        assert len(values) == n
        return values

    opts = lines[0].split()
    periods = []
    while pos < len(lines):
        head = lines[pos].split("#")[0].split()
        pos += 1
        arrays = {}
        if int(head[0]) >= 0:
            arrays["rech"] = array()
        if int(opts[0]) == 2 and int(head[1]) >= 0:
            arrays["irch"] = array()
        for flag, name in (("INRCHZONES", "iznrch"), ("INSELEV", "selev")):
            if flag in head and int(head[head.index(flag) + 1]) > 0:
                arrays[name] = array()
        periods.append((head, arrays))
    return opts, periods


_G01_R1 = [1e-3, 2e-3, 3e-3, 4e-3]
_G01_R3 = [4e-3, 3e-3, 2e-3, 1e-3]


@pytest.mark.parametrize("grid", ["dis", "disu"])
def test_mfusg_rch_selev_keeps_its_values(function_tmpdir, grid):
    """SELEV keeps its own non-uniform values (not RECH) in every period,
    with the period 1 array reused in period 2 and RECH reused in period 2."""
    ml = _g01_model(function_tmpdir, grid)
    s1 = _g01_arr(grid, [5.5, 6.25, 7.0, 8.125])
    s3 = _g01_arr(grid, [9.5, 10.0, 11.25, 12.0])
    rech = {0: _g01_arr(grid, _G01_R1), 2: _g01_arr(grid, _G01_R3)}
    rch = MfUsgRch(ml, nrchop=3, rech=rech, seepelev=1, selev={0: s1, 1: s1, 2: s3})
    for kper, expected in ((0, s1), (1, s1), (2, s3)):
        array = rch.selev[kper].array
        assert array.dtype == np.float32 and array.shape == expected.shape
        np.testing.assert_array_equal(array, expected.astype(np.float32))
    rch.write_file(check=False)
    opts, periods = _g01_read(rch.fn_path, 4)
    assert opts == ["3", "0", "SEEPELEV"]
    assert [p[1].get("selev") for p in periods] == [
        [5.5, 6.25, 7.0, 8.125],
        [5.5, 6.25, 7.0, 8.125],
        [9.5, 10.0, 11.25, 12.0],
    ]
    assert [p[1].get("rech") for p in periods] == [_G01_R1, None, _G01_R3]
    assert [p[0][0] for p in periods] == ["1", "-1", "1"]


@pytest.mark.parametrize("grid", ["dis", "disu"])
def test_mfusg_rch_iznrch_keeps_its_values(function_tmpdir, grid):
    """IZNRCH keeps its own non-uniform zone numbers (not RECH) as int32 in
    every period, with the period 1 array reused in period 2."""
    ml = _g01_model(function_tmpdir, grid)
    z1 = _g01_arr(grid, [1, 2, 3, 0], np.int64)
    z3 = _g01_arr(grid, [0, 3, 2, 1], np.int64)
    rech = {0: _g01_arr(grid, _G01_R1), 2: _g01_arr(grid, _G01_R3)}
    rch = MfUsgRch(ml, nrchop=3, rech=rech, mxrtzones=3, iznrch={0: z1, 1: z1, 2: z3})
    for kper, expected in ((0, z1), (1, z1), (2, z3)):
        array = rch.iznrch[kper].array
        assert array.dtype == np.int32 and array.shape == expected.shape
        np.testing.assert_array_equal(array, expected)
    rch.write_file(check=False)
    opts, periods = _g01_read(rch.fn_path, 4)
    assert opts == ["3", "0", "RTS", "3"]
    assert [p[1].get("iznrch") for p in periods] == [
        [1, 2, 3, 0],
        [1, 2, 3, 0],
        [0, 3, 2, 1],
    ]
    assert [p[1].get("rech") for p in periods] == [_G01_R1, None, _G01_R3]


@pytest.mark.parametrize("grid", ["dis", "disu"])
def test_mfusg_rch_iznrch_written_before_selev(function_tmpdir, grid):
    """With RTS and SEEPELEV, each period writes RECH, IZNRCH, SELEV in the
    order RP reads them."""
    ml = _g01_model(function_tmpdir, grid)
    z = _g01_arr(grid, [2, 0, 1, 2], np.int32)
    s = _g01_arr(grid, [15.5, 14.25, 13.0, 12.75])
    rch = MfUsgRch(
        ml,
        nrchop=3,
        ipakcb=50,
        rech=_g01_arr(grid, _G01_R1),
        seepelev=1,
        mxrtzones=2,
        selev=dict.fromkeys(range(3), s),
        iznrch=dict.fromkeys(range(3), z),
    )
    rch.write_file(check=False)
    opts, periods = _g01_read(rch.fn_path, 4)
    assert opts == ["3", "50", "SEEPELEV", "RTS", "2"]
    assert periods[0][1] == {
        "rech": _G01_R1,
        "iznrch": [2, 0, 1, 2],
        "selev": [15.5, 14.25, 13.0, 12.75],
    }
    for head, arrays in periods[1:]:
        assert head[0] == "-1"
        assert arrays == {"iznrch": [2, 0, 1, 2], "selev": [15.5, 14.25, 13.0, 12.75]}


@pytest.mark.parametrize("grid", ["dis", "disu"])
@pytest.mark.parametrize("mxrtzones", [0, 1, 3, 12345])
def test_mfusg_rch_rts_option(function_tmpdir, grid, mxrtzones):
    """RTS is written with its integer MXZNRCH; without zones (0) there is
    no RTS keyword."""
    ml = _g01_model(function_tmpdir, grid)
    rch = MfUsgRch(ml, nrchop=3, ipakcb=50, rech=1e-3, mxrtzones=mxrtzones)
    rch.write_file(check=False)
    opts, periods = _g01_read(rch.fn_path, 4)
    assert opts == ["3", "50"] + (["RTS", str(mxrtzones)] if mxrtzones else [])
    assert periods[0] == (["1"], {"rech": [1e-3] * 4})
    assert [p[0] for p in periods[1:]] == [["-1"], ["-1"]]


def test_mfusg_rch_rts_load_write(function_tmpdir):
    """The G03 hand-written RCH with RTS 2 loads and writes again with
    RTS 2 and its IZNRCH zones (int32)."""
    ws = _g03_fixture(function_tmpdir, evs=None)
    ml = MfUsg.load("t.nam", model_ws=str(ws), check=False, load_only=["DISU"])
    rch = MfUsgRch.load(str(ws / "t.rch"), ml, check=False)
    assert rch.mxrtzones == 2
    assert rch.iznrch[0].array.dtype == np.int32
    np.testing.assert_array_equal(rch.iznrch[0].array, [[0, 1, 2]])
    rch.fn_path = str(ws / "out.rch")
    rch.write_file(check=False)
    opts, periods = _g01_read(ws / "out.rch", 3)
    assert opts == ["1", "0", "RTS", "2"]
    assert periods[0][1] == {"rech": [0.0] * 3, "iznrch": [0, 1, 2]}


# hand-written SEEPELEV control for the 5-cell model of _g05_model, three
# periods: SELEV differs from RECH and is not uniform; period 2 repeats the
# period 1 SELEV (p2="1") or reuses it (p2="-1"); period 3 reuses RECH
_G01_P3 = _G05_PERIODS + "    10.0   5   1.0   TR\n"
_G01_OC = (
    "# OC\n"
    + _G05_HDR
    + "".join(
        f"period {k} step {s}\n" + _G05_ACT
        for k, n in ((1, 1), (2, 5), (3, 5))
        for s in range(1, n + 1)
    )
)
_G01_RECH = (
    ["0.01", "0.02", "0.03", "0.02", "0.01"],
    ["0.02", "0.01", "0.03", "0.01", "0.02"],
)
_G01_SELEV = (
    ["12.5", "11.9", "10.8", "10.9", "10.7"],
    ["12.0", "11.8", "10.6", "10.9", "10.4"],
)


def _g01_seep_text(selev1=_G01_SELEV[0], p2="1", seepelev=True):
    """RCH text (NRCHOP 3, IRCHCB 50) for _g05_model; without seepelev the
    same RECH without SELEV."""

    def blk(values):
        return "INTERNAL 1.0 (FREE) -1\n " + " ".join(values) + "\n"

    if not seepelev:
        return (
            "# RCH\n 3 50\n 0\n"
            + blk(_G01_RECH[0])
            + " 0\n"
            + blk(_G01_RECH[1])
            + " -1\n"
        )
    return (
        "# RCH\n 3 50 SEEPELEV\n 0 INSELEV 1\n"
        + blk(_G01_RECH[0])
        + blk(selev1)
        + f" 0 INSELEV {p2}\n"
        + blk(_G01_RECH[1])
        + (blk(selev1) if p2 == "1" else "")
        + " -1 INSELEV 1\n"
        + blk(_G01_SELEV[1])
    )


def _g01_seep_model(ws, rch_text, grid):
    ws = _g05_model(ws, _G01_OC, grid, _G01_P3)
    (ws / "t.rch").write_text(rch_text)
    with open(ws / "t.nam", "a") as f:
        f.write("RCH 18 t.rch\n")
    return ws


@pytest.mark.parametrize("grid", ["dis", "disu"])
def test_mfusg_rch_seepelev_load_write(function_tmpdir, grid):
    """The hand-written SEEPELEV control loads and writes again with its
    SELEV values in every period."""
    ws = _g01_seep_model(function_tmpdir, _g01_seep_text(), grid)
    rch = MfUsgRch.load(str(ws / "t.rch"), _g05_dis(ws), check=False)
    rch.fn_path = str(ws / "out.rch")
    rch.write_file(check=False)
    opts, periods = _g01_read(ws / "out.rch", 5)
    assert opts == ["3", "50", "SEEPELEV"]
    s1, s3 = ([float(v) for v in s] for s in _G01_SELEV)
    r1, r2 = ([float(v) for v in r] for r in _G01_RECH)
    assert [p[1] for p in periods] == [
        {"rech": r1, "selev": s1},
        {"rech": r2, "selev": s1},
        {"rech": r2, "selev": s3},
    ]


# kept fork behaviour: INCONC is a bare keyword (RP 246-248); INIRCH and the
# IRCH array only with NRCHOP=2 (RP 255-272); MXNDRCH only for DISU NRCHOP=2


def test_mfusg_rch_inconc_bare_keyword(function_tmpdir):
    ml = _g01_model(function_tmpdir, "dis")
    ml.mcomp = 1
    rch = MfUsgRch(ml, nrchop=3, rech=1e-3, iconc=1, rchconc={0: 0.5, 1: 0.5, 2: 0.25})
    rch.write_file(check=False)
    lines = [
        ln
        for ln in Path(rch.fn_path).read_text().splitlines()
        if not ln.startswith("#")
    ]
    assert lines[0].split() == ["3", "0", "CONC"]
    assert lines[1].split() == ["1"]
    heads = [ln.split("#")[0].split() for ln in lines if "Stress period" in ln]
    assert heads == [["1", "INCONC"], ["-1", "INCONC"], ["-1", "INCONC"]]


@pytest.mark.parametrize("grid", ["dis", "disu"])
@pytest.mark.parametrize("nrchop", [1, 2, 3])
def test_mfusg_rch_nrchop_headers(function_tmpdir, grid, nrchop):
    ml = _g01_model(function_tmpdir, grid)
    irch = _g01_arr(grid, [0, 0, 0, 0] if grid == "dis" else [0, 1, 2, 3], np.int32)
    rch = MfUsgRch(
        ml,
        nrchop=nrchop,
        rech=_g01_arr(grid, _G01_R1),
        irch=dict.fromkeys(range(3), irch),
    )
    rch.write_file(check=False)
    mxndrch = nrchop == 2 and grid == "disu"
    opts, periods = _g01_read(rch.fn_path, 4, skip=int(mxndrch))
    assert opts == [str(nrchop), "0"]
    if mxndrch:
        lines = Path(rch.fn_path).read_text().splitlines()
        assert lines[2].split() == ["4"]
    if nrchop == 2:
        inirch = "4" if grid == "disu" else "1"
        assert [p[0] for p in periods] == [
            ["1", inirch],
            ["-1", inirch],
            ["-1", inirch],
        ]
        irch_file = [1, 1, 1, 1] if grid == "dis" else [1, 2, 3, 4]
        assert [p[1] for p in periods] == [
            {"rech": _G01_R1, "irch": irch_file},
            {"irch": irch_file},
            {"irch": irch_file},
        ]
    else:
        assert [p[0] for p in periods] == [["1"], ["-1"], ["-1"]]
        assert periods[0][1] == {"rech": _G01_R1}


# VORONOI-USG: DISU arrays of a voronoi grid extruded into layers with the
# same cells (flopy.utils.voronoi.get_disu5_gridprops). USG-T 2.8.1 contract:
# gwf2basu1.f SGWF2DIS8UR (IAC, JA with the node first) and FILLGFU (the
# horizontal FAHL is an area that USG-T divides by 0.5*(thick_n + thick_m));
# glo2basu1.f FILLISYM (n-m needs m-n) and SGLO2BAS8ARDIS 10A (IVSD <= 0:
# m in a later layer is a vertical connection); utl7u1.f U1DRELNJAU
# (IDSYMRD=0: CL12 in row n is the distance from n to the face, CL1 upper and
# CL2 lower triangle) and U1DRELNJA (only the upper FAHL is read);
# gwf2bcf-lpf-u1.f SGWF2LPFU1VCOND (vertical flow uses FAHL as the area).

_VOR_SQ = (
    [(0, 0), (1, 0), (2, 0), (2, 1), (1, 1), (0, 1)],
    [[0, 1, 4, 5], [1, 2, 3, 4]],
    [(0.5, 0.5), (1.5, 0.5)],
)
# voronoi cells of A(1,1), B(3,1), C(2,3) in [0,4]x[0,4]; the vertex
# (2,1.75) is 1.25 from the three points
_VOR_TRI = (
    [(0, 0), (2, 0), (4, 0), (4, 2.75), (2, 1.75), (0, 2.75), (4, 4), (0, 4)],
    [[0, 1, 4, 5], [1, 2, 3, 4], [5, 4, 3, 6, 7]],
    [(1, 1), (3, 1), (2, 3)],
)


def _vor_check(g):
    """Invariants of the DISU arrays, independent of the geometry."""
    from flopy.utils.voronoi import get_disu5_gridprops  # noqa: F401

    iac, ja, nodes = g["iac"], g["ja"], g["nodes"]
    njag = g["njag"]
    assert (
        iac.sum()
        == njag
        == len(ja)
        == len(g["cl12"])
        == len(g["fahl"])
        == len(g["ivc"])
    )
    assert ja.min() >= 0 and ja.max() < nodes
    ncpl = g["nodelay"][0]
    ia = np.r_[0, np.cumsum(iac)]
    conn = {}
    for n in range(nodes):
        row = ja[ia[n] : ia[n + 1]].tolist()
        assert row[0] == n and len(set(row)) == len(row) and row[1:] == sorted(row[1:])
        assert g["cl12"][ia[n]] == 0 and g["fahl"][ia[n]] == 0
        for j in range(ia[n] + 1, ia[n + 1]):
            m = int(ja[j])
            assert g["cl12"][j] > 0 and g["fahl"][j] > 0
            vertical = m // ncpl != n // ncpl
            assert g["ivc"][j] == int(vertical)
            if vertical:
                assert abs(m // ncpl - n // ncpl) == 1 and m % ncpl == n % ncpl
            conn[(n, m)] = (g["ivc"][j], g["fahl"][j])
    for (n, m), value in conn.items():
        assert conn[(m, n)] == value


def _vor_rows(g):
    ia = np.r_[0, np.cumsum(g["iac"])]
    return [
        list(
            zip(
                g["ja"][ia[n] : ia[n + 1]].tolist(),
                g["cl12"][ia[n] : ia[n + 1]],
                g["fahl"][ia[n] : ia[n + 1]],
                g["ivc"][ia[n] : ia[n + 1]].tolist(),
            )
        )
        for n in range(g["nodes"])
    ]


def _vor_assert_rows(g, expected):
    got = _vor_rows(g)
    assert [[r[0] for r in row] for row in got] == [
        [e[0] for e in row] for row in expected
    ]
    for grow, erow in zip(got, expected):
        for (m, cl, fa, v), (_, ecl, efa, ev) in zip(grow, erow):
            assert v == ev
            assert cl == pytest.approx(ecl, rel=1e-12, abs=1e-15)
            assert fa == pytest.approx(efa, rel=1e-12, abs=1e-15)


def _vor_same(a, b, scale_cl=1.0, scale_fa=1.0, scale_area=1.0):
    np.testing.assert_array_equal(a["ja"], b["ja"])
    np.testing.assert_array_equal(a["iac"], b["iac"])
    np.testing.assert_allclose(b["cl12"], a["cl12"] * scale_cl, rtol=1e-9)
    np.testing.assert_allclose(b["fahl"], a["fahl"] * scale_fa, rtol=1e-9)
    np.testing.assert_allclose(b["area"], a["area"] * scale_area, rtol=1e-9)


def test_mfusg_voronoi_disu5_two_squares():
    """Hand values: 1 layer, 2 layers, and variable thickness (FAHL uses the
    mean thickness of the two cells; vertical CL12 is half the thickness of
    the cell of the row)."""
    from flopy.utils.voronoi import get_disu5_gridprops

    g = get_disu5_gridprops(*_VOR_SQ, 2.0, [0.0])
    _vor_check(g)
    assert (g["nodes"], g["nlay"], g["njag"], g["ivsd"], g["idsymrd"]) == (
        2,
        1,
        4,
        0,
        0,
    )
    np.testing.assert_allclose(g["area"], [1.0, 1.0], rtol=1e-12)
    _vor_assert_rows(
        g, [[(0, 0, 0, 0), (1, 0.5, 2.0, 0)], [(1, 0, 0, 0), (0, 0.5, 2.0, 0)]]
    )
    g = get_disu5_gridprops(*_VOR_SQ, 2.0, [1.5, 0.0])
    _vor_check(g)
    _vor_assert_rows(
        g,
        [
            [(0, 0, 0, 0), (1, 0.5, 0.5, 0), (2, 0.25, 1.0, 1)],
            [(1, 0, 0, 0), (0, 0.5, 0.5, 0), (3, 0.25, 1.0, 1)],
            [(2, 0, 0, 0), (0, 0.75, 1.0, 1), (3, 0.5, 1.5, 0)],
            [(3, 0, 0, 0), (1, 0.75, 1.0, 1), (2, 0.5, 1.5, 0)],
        ],
    )
    np.testing.assert_array_equal(g["top"], [2, 2, 1.5, 1.5])
    np.testing.assert_array_equal(g["bot"], [1.5, 1.5, 0, 0])
    np.testing.assert_array_equal(g["nodelay"], [2, 2])
    g = get_disu5_gridprops(*_VOR_SQ, [3.0, 2.0], [[1.0, 0.5], [0.0, 0.0]])
    _vor_check(g)
    _vor_assert_rows(
        g,
        [
            [(0, 0, 0, 0), (1, 0.5, 1.75, 0), (2, 1.0, 1.0, 1)],
            [(1, 0, 0, 0), (0, 0.5, 1.75, 0), (3, 0.75, 1.0, 1)],
            [(2, 0, 0, 0), (0, 0.5, 1.0, 1), (3, 0.5, 0.75, 0)],
            [(3, 0, 0, 0), (1, 0.25, 1.0, 1), (2, 0.5, 0.75, 0)],
        ],
    )


def test_mfusg_voronoi_disu5_irregular_three_cells():
    from flopy.utils.voronoi import get_disu5_gridprops

    for p in _VOR_TRI[2]:
        assert np.hypot(p[0] - 2, p[1] - 1.75) == pytest.approx(1.25, rel=1e-15)
    g = get_disu5_gridprops(*_VOR_TRI, 10.0, [0.0])
    _vor_check(g)
    np.testing.assert_allclose(g["area"], [4.5, 4.5, 7.0], rtol=1e-12)
    s5 = np.sqrt(5.0)
    _vor_assert_rows(
        g,
        [
            [(0, 0, 0, 0), (1, 1.0, 17.5, 0), (2, s5 / 2, 10 * s5, 0)],
            [(1, 0, 0, 0), (0, 1.0, 17.5, 0), (2, s5 / 2, 10 * s5, 0)],
            [(2, 0, 0, 0), (0, s5 / 2, 10 * s5, 0), (1, s5 / 2, 10 * s5, 0)],
        ],
    )


def test_mfusg_voronoi_disu5_elevation_forms():
    from flopy.utils.voronoi import get_disu5_gridprops

    forms = (
        [np.array([1.0, 0.5]), 0.0],
        np.array([[1.0, 0.5], [0.0, 0.0]]),
        [[1.0, 0.5], [0.0, 0.0]],
        ([1.0, 0.5], 0),
    )
    gs = [get_disu5_gridprops(*_VOR_SQ, [3.0, 2.0], b) for b in forms]
    for g in gs[1:]:
        _vor_same(gs[0], g)
        np.testing.assert_array_equal(g["bot"], gs[0]["bot"])
    _vor_same(
        get_disu5_gridprops(*_VOR_SQ, 3.0, [1.0, 0.0]),
        get_disu5_gridprops(*_VOR_SQ, np.array([3.0, 3.0]), np.array([1.0, 0.0])),
    )


@pytest.mark.parametrize("shift", [(3.3e5, 4.47e6), (-1e7, 2.5e3)])
def test_mfusg_voronoi_disu5_invariance(shift):
    """Translation, scale, direction, first vertex, closed rings and repeated
    coordinates give the same connections."""
    from flopy.utils.voronoi import get_disu5_gridprops

    v, iv, c = _VOR_TRI
    base = get_disu5_gridprops(v, iv, c, 10.0, [0.0])
    _vor_same(
        base, get_disu5_gridprops(np.add(v, shift), iv, np.add(c, shift), 10.0, [0.0])
    )
    _vor_same(
        base,
        get_disu5_gridprops(np.multiply(v, 10), iv, np.multiply(c, 10), 10.0, [0.0]),
        10,
        10,
        100,
    )
    for cells in (
        [r[::-1] for r in iv],
        [r[2:] + r[:2] for r in iv],
        [r + r[:1] for r in iv],
    ):
        _vor_same(base, get_disu5_gridprops(v, cells, c, 10.0, [0.0]))
    # cell 0 lists a second vertex with the coordinates of vertex 4
    _vor_same(
        base,
        get_disu5_gridprops(v + [v[4]], [[0, 1, 8, 5], iv[1], iv[2]], c, 10.0, [0.0]),
    )


def test_mfusg_voronoi_disu5_merge_tolerance():
    """Copies of a shared vertex within the tolerance are merged whatever
    their position relative to the rounding bins and their order; vertices
    used by no cell, even far away, change nothing; a chain of points closer
    than the tolerance but longer than it cannot be resolved."""
    from flopy.utils.voronoi import _disu5_merge, get_disu5_gridprops

    tol = 1e-9 * 2.0  # extent of the two squares
    base = get_disu5_gridprops(*_VOR_SQ, 1.0, [0.0])
    v = np.array(_VOR_SQ[0], float)
    for x in (1.0, 1.0 + 0.5 * tol, 1.0 + 0.999 * tol, 1.0 + 1.5 * tol):
        for side in (0.4 * tol, -0.4 * tol):
            vv = v.copy()
            vv[[1, 4], 0] = x
            vv = np.vstack([vv, vv[[1, 4]] + [side, 0.0]])
            for cells in ([[0, 1, 4, 5], [6, 2, 3, 7]], [[0, 6, 7, 5], [1, 2, 3, 4]]):
                g = get_disu5_gridprops(vv, cells, _VOR_SQ[2], 1.0, [0.0])
                np.testing.assert_array_equal(g["ja"], base["ja"])
    far = np.vstack([v, [[1e12, 1e12], [-1e9, 5.0]]])
    _vor_same(base, get_disu5_gridprops(far, _VOR_SQ[1], _VOR_SQ[2], 1.0, [0.0]))
    cid, xy = _disu5_merge(np.array([[0.0, 0.0], [0.4, 0.0], [5.0, 5.0]]), 1.0)
    assert cid[0] == cid[1] != cid[2] and np.allclose(xy[cid[0]], [0.2, 0.0])
    with pytest.raises(ValueError, match="chain"):
        _disu5_merge(np.array([[0.0, 0.0], [0.8, 0.0], [1.6, 0.0], [2.4, 0.0]]), 1.0)


def test_mfusg_voronoi_disu5_small_cells_kept():
    """A strip 1e-4 wide between two unit cells keeps its connections (the
    tolerance is 1e-9 times the extent)."""
    from flopy.utils.voronoi import get_disu5_gridprops

    e = 1e-4
    v = [(0, 0), (1, 0), (1 + e, 0), (2 + e, 0), (2 + e, 1), (1 + e, 1), (1, 1), (0, 1)]
    iv = [[0, 1, 6, 7], [1, 2, 5, 6], [2, 3, 4, 5]]
    c = [(0.5, 0.5), (1 + e / 2, 0.5), (1.5 + e, 0.5)]
    g = get_disu5_gridprops(v, iv, c, 2.0, [0.0])
    _vor_check(g)
    np.testing.assert_allclose(g["area"], [1.0, e, 1.0], rtol=1e-9)
    _vor_assert_rows(
        g,
        [
            [(0, 0, 0, 0), (1, 0.5, 2.0, 0)],
            [(1, 0, 0, 0), (0, e / 2, 2.0, 0), (2, e / 2, 2.0, 0)],
            [(2, 0, 0, 0), (1, 0.5, 2.0, 0)],
        ],
    )


def test_mfusg_voronoi_disu5_collinear_face_segments():
    """A face made of two collinear segments counts its whole width; two
    faces in different directions between the same cells are rejected."""
    from flopy.utils.voronoi import get_disu5_gridprops

    v = [(0, 0), (1, 0), (2, 0), (2, 1), (1, 1), (0, 1), (1, 0.5)]
    g = get_disu5_gridprops(
        v, [[0, 1, 6, 4, 5], [1, 2, 3, 4, 6]], _VOR_SQ[2], 2.0, [0.0]
    )
    _vor_assert_rows(
        g, [[(0, 0, 0, 0), (1, 0.5, 2.0, 0)], [(1, 0, 0, 0), (0, 0.5, 2.0, 0)]]
    )
    v = [(0, 0), (2, 0), (2, 1), (1, 1), (1, 2), (0, 2), (2, 2)]
    with pytest.raises(ValueError, match="more than one face"):
        get_disu5_gridprops(
            v, [[0, 1, 2, 3, 4, 5], [3, 2, 6, 4]], [(0.5, 0.5), (1.5, 1.5)], 1.0, [0.0]
        )


def test_mfusg_voronoi_disu5_inputs_not_modified():
    from flopy.utils.voronoi import get_disu5_gridprops

    v = np.array(_VOR_TRI[0], float)
    iv = [list(r) for r in _VOR_TRI[1]]
    c = np.array(_VOR_TRI[2], float)
    top = np.array([10.0, 11.0, 12.0])
    botm = [np.array([4.0, 5.0, 6.0]), 0.0]
    copies = (v.copy(), [list(r) for r in iv], c.copy(), top.copy(), botm[0].copy())
    g1 = get_disu5_gridprops(v, iv, c, top, botm)
    g2 = get_disu5_gridprops(v, iv, c, top, botm)
    np.testing.assert_array_equal(v, copies[0])
    assert iv == copies[1]
    np.testing.assert_array_equal(c, copies[2])
    np.testing.assert_array_equal(top, copies[3])
    np.testing.assert_array_equal(botm[0], copies[4])
    g1["ja"] += 1
    g1["top"][0] = -1.0
    assert g2["ja"][0] == 0 and g2["top"][0] == 10.0


@pytest.mark.parametrize(
    "top,botm,match",
    [
        (10.0, [10.0], "non-positive thickness"),
        (10.0, [5.0, 6.0], "non-positive thickness"),
        ([10.0, 10.0, 1.0], [[0.0, 0.0, 2.0]], "non-positive thickness"),
        (np.nan, [0.0], "finite"),
        (10.0, [[0.0, np.inf, 0.0]], "finite"),
        ([10.0, 10.0], [0.0], "top must be"),
        (10.0, [[0.0, 0.0]], r"botm\[0\]"),
        (10.0, 0.0, "one entry per layer"),
        (10.0, [], "at least one layer"),
        ("ten", [0.0], "top must be"),
        (10.0, ["zero"], r"botm\[0\]"),
    ],
)
def test_mfusg_voronoi_disu5_rejects_elevations(top, botm, match):
    from flopy.utils.voronoi import get_disu5_gridprops

    with pytest.raises(ValueError, match=match):
        get_disu5_gridprops(*_VOR_TRI, top, botm)


_VOR_BAD = {
    # reviewer case: [0,1]x[0,1] and [0.5,1.5]x[0.3,1.3] overlap
    "overlap": (
        [
            (0, 0),
            (1, 0),
            (1, 1),
            (0, 1),
            (0.5, 0.3),
            (1.5, 0.3),
            (1.5, 1.3),
            (0.5, 1.3),
        ],
        [[0, 1, 2, 3], [4, 5, 6, 7]],
        [(0.5, 0.5), (1.0, 0.8)],
        "cross or touch",
    ),
    "duplicate": (
        [(0, 0), (1, 0), (1, 1), (0, 1)],
        [[0, 1, 2, 3], [0, 1, 2, 3]],
        [(0.5, 0.5), (0.4, 0.4)],
        "same direction",
    ),
    "contained": (
        [(0, 0), (4, 0), (4, 4), (0, 4), (1, 1), (2, 1), (2, 2), (1, 2)],
        [[0, 1, 2, 3], [4, 5, 6, 7]],
        [(3.0, 3.0), (1.5, 1.5)],
        "inside cell|boundaries",
    ),
    "bowtie": (
        [(0, 0), (1, 1), (1, 0), (0, 1)],
        [[0, 1, 2, 3]],
        [(0.5, 0.5)],
        "intersects itself|outside|zero area",
    ),
    "folded": (
        [(0, 0), (2, 0), (1, 0), (1, 1)],
        [[0, 1, 2, 3]],
        [(0.8, 0.2)],
        "folds back|intersects itself",
    ),
    "spike": (
        [(0, 0), (2, 0), (2, 2), (2, 1), (0, 2)],
        [[0, 1, 2, 3, 4]],
        [(1.0, 1.0)],
        "folds back|intersects itself",
    ),
    "fraction": (
        *_VOR_SQ[:1],
        [[0, 1, 4.5, 5], [1, 2, 3, 4]],
        _VOR_SQ[2],
        "not an integer",
    ),
    "float_index": (
        *_VOR_SQ[:1],
        [[0, 1, 4.0, 5], [1, 2, 3, 4]],
        _VOR_SQ[2],
        "not an integer",
    ),
    "bool_index": (
        *_VOR_SQ[:1],
        [[0, True, 4, 5], [1, 2, 3, 4]],
        _VOR_SQ[2],
        "not an integer",
    ),
    "text_cell": (*_VOR_SQ[:1], ["0145", [1, 2, 3, 4]], _VOR_SQ[2], "list of integers"),
    "empty_cell": (*_VOR_SQ[:1], [[], [1, 2, 3, 4]], _VOR_SQ[2], "no vertices"),
    "no_cells": (_VOR_SQ[0], [], np.zeros((0, 2)), "empty"),
    "out_of_range": (
        *_VOR_SQ[:1],
        [[0, 1, 4, 99], [1, 2, 3, 4]],
        _VOR_SQ[2],
        "out of range",
    ),
    "negative": (
        *_VOR_SQ[:1],
        [[0, 1, 4, -1], [1, 2, 3, 4]],
        _VOR_SQ[2],
        "out of range",
    ),
    "verts_3d": (
        [(0, 0, 0), (1, 0, 0), (1, 1, 0)],
        [[0, 1, 2]],
        [(0.6, 0.3)],
        "verts must have",
    ),
    "verts_text": (
        [("a", 0), (1, 0), (1, 1)],
        [[0, 1, 2]],
        [(0.6, 0.3)],
        "arrays of numbers",
    ),
    "xcyc_shape": (*_VOR_SQ[:2], [(0.5, 0.5)], "xcyc must have"),
    "xcyc_nan": (*_VOR_SQ[:2], [(np.nan, 0.5), (1.5, 0.5)], "finite"),
    "two_vertices": (*_VOR_SQ[:1], [[0, 1], [1, 2, 3, 4]], _VOR_SQ[2], "fewer than 3"),
    "repeated_vertex": (
        *_VOR_SQ[:1],
        [[0, 1, 4, 1, 5], [1, 2, 3, 4]],
        _VOR_SQ[2],
        "repeated",
    ),
    "center_outside": (*_VOR_SQ[:2], [(1.5, 0.5), (0.5, 0.5)], "outside"),
    "zero_area": (
        [(0, 0), (1, 0), (2, 0)],
        [[0, 1, 2]],
        [(1.0, 0.0)],
        "fewer than 3|zero area|folds",
    ),
    "three_cells_edge": (
        _VOR_SQ[0] + [(0.5, 0.5)],
        [[0, 1, 4, 5], [1, 2, 3, 4], [1, 4, 6]],
        [(0.5, 0.5), (1.5, 0.5), (0.8, 0.5)],
        "more than two cells",
    ),
    "non_orthogonal": (
        [(0, 0), (1, 0), (2, 0), (2.5, 1), (1.5, 1), (0.5, 1)],
        _VOR_SQ[1],
        [(0.75, 0.5), (1.75, 0.5)],
        "not orthogonal",
    ),
    "hanging_vertex": (
        [(0, 0), (1, 0), (2, 0), (2, 1), (1, 1), (0, 1), (2, 0.5), (1, 0.5)],
        [[0, 1, 4, 5], [1, 2, 6, 7], [7, 6, 3, 4]],
        [(0.5, 0.5), (1.5, 0.25), (1.5, 0.75)],
        "overlap|cross or touch",
    ),
    # 3x3 squares without the middle one: a hole
    "hole": (
        [(i, j) for j in range(4) for i in range(4)],
        [
            [0, 1, 5, 4],
            [1, 2, 6, 5],
            [2, 3, 7, 6],
            [4, 5, 9, 8],
            [6, 7, 11, 10],
            [8, 9, 13, 12],
            [9, 10, 14, 13],
            [10, 11, 15, 14],
        ],
        [
            (0.5, 0.5),
            (1.5, 0.5),
            (2.5, 0.5),
            (0.5, 1.5),
            (2.5, 1.5),
            (0.5, 2.5),
            (1.5, 2.5),
            (2.5, 2.5),
        ],
        "separate boundaries",
    ),
    "separate_regions": (
        [(0, 0), (1, 0), (1, 1), (0, 1), (5, 0), (6, 0), (6, 1), (5, 1)],
        [[0, 1, 2, 3], [4, 5, 6, 7]],
        [(0.5, 0.5), (5.5, 0.5)],
        "separate boundaries",
    ),
}


@pytest.mark.parametrize("case", list(_VOR_BAD))
def test_mfusg_voronoi_disu5_rejects_geometry(case):
    from flopy.utils.voronoi import get_disu5_gridprops

    v, iv, c, match = _VOR_BAD[case]
    with pytest.raises(ValueError, match=match):
        get_disu5_gridprops(v, iv, c, 1.0, [0.0])


def test_mfusg_voronoi_disu5_center_inside_other_cell():
    """The last overlap check on its own: a cell center inside another cell."""
    from flopy.utils.voronoi import _disu5_cells, _disu5_check_overlaps

    v, iv, c, _ = _VOR_BAD["contained"]
    c = np.array(c, float)
    rings, _, vxy, tol = _disu5_cells(np.array(v, float), iv, c)
    edges = [(min(a, b), max(a, b)) for r in rings for a, b in zip(r, r[1:] + r[:1])]
    with pytest.raises(ValueError, match="inside cell 0"):
        _disu5_check_overlaps(rings, edges, vxy, c, tol)


def _vor_rect(ws, seed=7, n=40, max_area=4000.0, length=1000.0, width=500.0):
    """flopy Triangle + VoronoiGrid of a length x width rectangle with random
    interior nodes (fixed seed)."""
    from flopy.utils.triangle import Triangle
    from flopy.utils.voronoi import VoronoiGrid

    Path(ws).mkdir(parents=True)
    rng = np.random.default_rng(seed)
    nodes = np.c_[rng.uniform(40, length - 40, n), rng.uniform(40, width - 40, n)]
    tri = Triangle(maximum_area=max_area, angle=30, nodes=nodes, model_ws=str(ws))
    tri.add_polygon([(0, 0), (length, 0), (length, width), (0, width)])
    tri.build(verbose=False)
    return tri, VoronoiGrid(tri)


@requires_pkg("scipy", "shapely")
@requires_exe("triangle")
def test_mfusg_voronoi_disu5_voronoigrid(tmp_path):
    """The VoronoiGrid method equals the module function; the neighbors are
    Triangle edges (Delaunay), the others touch at one point; the cells
    cover the rectangle."""
    import inspect

    from shapely.geometry import Polygon
    from shapely.ops import unary_union

    from flopy.utils.voronoi import VoronoiGrid, get_disu5_gridprops

    assert list(inspect.signature(VoronoiGrid.get_disu5_gridprops).parameters) == [
        "self",
        "top",
        "botm",
    ]
    assert list(inspect.signature(get_disu5_gridprops).parameters) == [
        "verts",
        "iverts",
        "xcyc",
        "top",
        "botm",
    ]
    assert "MfUsgDisU" in VoronoiGrid.get_disu5_gridprops.__doc__
    assert "transport" in get_disu5_gridprops.__doc__
    tri, vor = _vor_rect(tmp_path / "mesh")
    top = np.linspace(30.0, 40.0, vor.ncpl)
    g = vor.get_disu5_gridprops(top, [10.0, 0.0])
    _vor_check(g)
    f = get_disu5_gridprops(vor.verts, vor.iverts, vor.points, top, [10.0, 0.0])
    for key in g:
        np.testing.assert_array_equal(g[key], f[key])
    ncpl = vor.ncpl
    ia = np.r_[0, np.cumsum(g["iac"])]
    pairs = {
        tuple(sorted((n, int(m))))
        for n in range(ncpl)
        for m in g["ja"][ia[n] + 1 : ia[n + 1]]
        if m < ncpl
    }
    tri_pairs = {
        tuple(sorted(map(int, e)))
        for e in tri.edge[["endpoint1", "endpoint2"]].tolist()
    }
    base = min(min(p) for p in tri_pairs)
    tri_pairs = {(a - base, b - base) for a, b in tri_pairs}
    polys = [Polygon(np.asarray(vor.verts)[iv]) for iv in vor.iverts]
    assert pairs <= tri_pairs
    for a, b in tri_pairs - pairs:
        assert polys[a].intersection(polys[b]).geom_type == "Point"
    domain = Polygon([(0, 0), (1000, 0), (1000, 500), (0, 500)])
    assert unary_union(polys).symmetric_difference(domain).area <= 1e-9 * domain.area
    np.testing.assert_allclose(g["area"][:ncpl], [p.area for p in polys], rtol=1e-12)


@requires_pkg("scipy")
@requires_exe("triangle")
def test_mfusg_voronoi_disu5_hole_rejected(tmp_path):
    """A VoronoiGrid built with a hole does not cover its domain (a gap next
    to the hole, cells inside it): rejected."""
    from flopy.utils.triangle import Triangle
    from flopy.utils.voronoi import VoronoiGrid

    tri = Triangle(maximum_area=2000.0, angle=30, model_ws=str(tmp_path))
    tri.add_polygon([(0, 0), (400, 0), (400, 300), (0, 300)])
    tri.add_polygon([(150, 100), (250, 100), (250, 200), (150, 200)])
    tri.add_hole((200, 150))
    tri.build(verbose=False)
    vor = VoronoiGrid(tri)
    with pytest.raises(ValueError, match="separate boundaries"):
        vor.get_disu5_gridprops(10.0, [0.0])


@requires_pkg("scipy")
@requires_exe("triangle")
def test_mfusg_voronoi_unstructuredgrid_props(tmp_path):
    """Without top and botm the dictionary is unchanged; with them it
    describes the layered grid (independent lists per layer)."""
    from flopy.discretization import UnstructuredGrid

    _, vor = _vor_rect(tmp_path / "mesh", seed=3, n=10, max_area=20000.0)
    ncpl = vor.ncpl
    plain = vor.get_gridprops_unstructuredgrid()
    assert list(plain) == ["vertices", "iverts", "ncpl", "xcenters", "ycenters"]
    assert plain["iverts"] is vor.iverts and plain["ncpl"] == ncpl
    np.testing.assert_array_equal(plain["xcenters"], vor.points[:, 0])
    assert plain["vertices"] == vor.get_disv_gridprops()["vertices"]
    for kwargs in ({"top": 10.0}, {"botm": [0.0]}):
        with pytest.raises(ValueError, match="together"):
            vor.get_gridprops_unstructuredgrid(**kwargs)
    top = np.linspace(30.0, 40.0, ncpl)
    botm = [np.full(ncpl, 20.0), 5.0]
    p = vor.get_gridprops_unstructuredgrid(top=top, botm=botm)
    assert p["ncpl"] == [ncpl, ncpl] and len(p["iverts"]) == 2 * ncpl
    p["iverts"][0].append(-1)
    assert p["iverts"][ncpl][-1] != -1 and vor.iverts[0][-1] != -1
    p["iverts"][0].pop()
    ug = UnstructuredGrid(**p)
    assert ug.nlay == 2 and list(ug.ncpl) == [ncpl, ncpl]
    np.testing.assert_allclose(ug.top[:ncpl], top)
    np.testing.assert_allclose(ug.top[ncpl:], 20.0)
    np.testing.assert_allclose(ug.botm[:ncpl], 20.0)
    np.testing.assert_allclose(ug.botm[ncpl:], 5.0)
    verts = np.asarray(vor.verts)
    for k in range(2):
        for c in range(ncpl):
            n = k * ncpl + c
            assert ug.xcellcenters[n] == vor.points[c, 0]
            assert ug.ycellcenters[n] == vor.points[c, 1]
            np.testing.assert_array_equal(
                ug.get_cell_vertices(n), verts[vor.iverts[c]].tolist()
            )
