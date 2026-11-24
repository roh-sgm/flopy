"""MfUsgRch: IRCHCONC per component, INCONC given once, DISU NRCHOP=2 with
IRCH reused or changed, and option keywords in any case.

USG-T 2.8.1 (gwf2rch8u1.f; 2.7.0 runs the same files):
- AR :91-97 CONC or CONCENTRATION turns on the recharge concentration;
  :112-120 reads MCOMPT values IRCHCONC (list-directed) and sizes RCHCONC
  with their sum; gwt2bndsu1.f uses a component only when its flag is 1.
- RP :227-253 reads the INRCHZONES / INSELEV / INCONC words in any case
  (URWORD, utl7u1.f :89-96); without INCONC the concentrations are reused
  (:466); with INCONC one array per component with IRCHCONC not 0 (:471-486).
- RP :264-272 DISU NRCHOP=2: INIRCH is the number of IRCH nodes, < 0 reuses
  IRCH (:338-342).
- gwt2bndsu1.f :3508, :3521 (:3586, :3599): the transport of a period that
  reuses IRCH adds no recharge mass, so FloPy does not write that form.
- Structured NRCHOP=2 (RP :338-360): INIRCH >= 0 reads IRCH for all the
  NROW*NCOL cells, but transport applies the recharge concentration only to
  the first INIRCH cells (gwt2bndsu1.f :3521), so the INIRCH read is kept,
  CONC writes NROW*NCOL by default and warns below it.
- RTS with NRCHOP=2: the RTS records read during a period update the
  recharge of the first INIRCH cells only and none when the period reuses
  IRCH (glo2basu1.f ATS1AJ :3342, :3353), so the flow depends on INIRCH: RTS
  follows the CONC rule (NROW*NCOL by default, warning below it, a period
  that reuses IRCH is not written).
  Period 1 must read IRCH (RP :338-342 would reuse an IRCH never read).
- Structured RTS (AR :103, :107; RP :438-439): IZNRCH has NROW*NCOL
  positions, so the zone of an IRCH cell below layer 1 is never read and
  that cell gets no RTS recharge (probe with both engines).
"""

import builtins
import re
import warnings
from pathlib import Path

import numpy as np
import pytest

from autotest.usgt_exe import USGT_EXE, requires_usgt
from flopy.mbase import run_model
from flopy.mfusg import (
    MfUsg,
    MfUsgBas,
    MfUsgBct,
    MfUsgDisU,
    MfUsgLpf,
    MfUsgOc,
    MfUsgRch,
    MfUsgSms,
)
from flopy.modflow import ModflowChd, ModflowDis
from flopy.utils import Util2d

_KEYWORDS = {
    "CONC",
    "CONCENTRATION",
    "RTS",
    "SEEPELEV",
    "INCONC",
    "INSELEV",
    "INRCHZONES",
}


# ------------------------------------------------------------- fixtures
def _adjacency(n, nlay):
    """iac/ja (zero-based) of nlay layers of n nodes in a line."""
    nbr = {i: [] for i in range(n * nlay)}
    for k in range(nlay):
        for j in range(n - 1):
            a, b = k * n + j, k * n + j + 1
            nbr[a].append(b)
            nbr[b].append(a)
    for k in range(nlay - 1):
        for j in range(n):
            a, b = k * n + j, (k + 1) * n + j
            nbr[a].append(b)
            nbr[b].append(a)
    iac, ja = [], []
    for i in range(n * nlay):
        iac.append(1 + len(nbr[i]))
        ja += [i, *sorted(nbr[i])]
    return np.array(iac), np.array(ja)


def _ml(ws, nper=2, mcomp=2, nlay=1, n=5, idisp=0):
    """DISU model of nlay layers of n nodes (BCT with mcomp components when
    mcomp > 0)."""
    ws = Path(ws)
    ws.mkdir(parents=True, exist_ok=True)
    ml = MfUsg(modelname="t", model_ws=str(ws), structured=False, exe_name=USGT_EXE)
    iac, ja = _adjacency(n, nlay)
    dz = 30.0 / nlay
    MfUsgDisU(
        ml,
        nodes=n * nlay,
        nlay=nlay,
        njag=len(ja),
        iac=iac,
        ja=ja,
        ivsd=0,
        nodelay=[n] * nlay,
        cl12=5.0,
        fahl=10.0,
        top=[30.0 - k * dz for k in range(nlay)],
        bot=[30.0 - (k + 1) * dz for k in range(nlay)],
        area=100.0,
        nper=nper,
        perlen=10.0,
        nstp=2,
        steady=False,
    )
    if mcomp:
        MfUsgBct(
            ml,
            itrnsp=1,
            mcomp=mcomp,
            idisp=idisp,
            dl=0.0,
            dt=0.0,
            diffnc=0.0,
            prsity=0.2,
            conc=0.0,
        )
    return ml


def _blk(values):
    return "INTERNAL 1.0 (FREE) -1\n " + " ".join(str(v) for v in values) + "\n"


def _iblk(values):
    return "INTERNAL 1 (FREE) -1\n " + " ".join(str(v) for v in values) + "\n"


_R5 = [1.0e-3, 2.0e-3, 3.0e-3, 4.0e-3, 5.0e-3]
_C5 = [5.0, 6.0, 7.0, 8.0, 9.0]
_D5 = [0.5, 0.25, 0.75, 1.5, 2.5]


def _controls(text):
    """Control records of an RCH file: tokens before '#', 'ARRAY' for each
    array (an INTERNAL header and its one data line, or a CONSTANT)."""
    lines = text.splitlines()
    if lines and lines[0].startswith("#"):
        lines = lines[1:]
    out, skip = [], False
    for line in lines:
        if skip:
            skip = False
            continue
        tokens = line.split("#")[0].split()
        if not tokens:
            continue
        if tokens[0].upper() == "INTERNAL":
            out.append(["ARRAY"])
            skip = True
        elif tokens[0].upper() == "CONSTANT":
            out.append(["ARRAY"])
        else:
            out.append(tokens)
    return out


def _values(rch):
    """Arrays of an MfUsgRch by stress period (keys as stored)."""
    out = {"rech": {k: v.array.tolist() for k, v in rch.rech.transient_2ds.items()}}
    if rch.irch is not None and rch.nrchop == 2:
        out["irch"] = {k: v.array.tolist() for k, v in rch.irch.transient_2ds.items()}
    if rch.rchconc is not None:
        out["rchconc"] = {
            k: [u3d[i].array.tolist() for i in range(u3d.shape[0])]
            for k, u3d in rch.rchconc.transient_3ds.items()
        }
    return out


def _write(rch, ws, name="out.rch"):
    rch.fn_path = str(Path(ws) / name)
    rch.write_file(check=False)
    return (Path(ws) / name).read_text()


def _reload(path, ml):
    """MfUsgRch.load into ml; a previous RCH package is removed first (one
    RCH package per model)."""
    if ml.get_package("RCH") is not None:
        ml.remove_package("RCH")
    return MfUsgRch.load(str(path), ml, check=False)


def _load(ws, text, ml, name="t.rch"):
    path = Path(ws) / name
    path.write_text(text)
    return _reload(path, ml)


def _load_write_reload(ws, text, ml):
    """(loaded, written text, reloaded) of a hand-written RCH file."""
    rch = _load(ws, text, ml)
    out = _write(rch, ws)
    again = _reload(Path(ws) / "out.rch", ml)
    return rch, out, again


# ------------------------------------------------------------ IRCHCONC
_H_F41 = (
    "# RCH\n 3 50 CONC\n 1 0\n 1 INCONC\n" + _blk(_R5) + _blk(_C5) + " 1\n" + _blk(_R5)
)


def test_rch_irchconc_1_0_hand_file_load_write(function_tmpdir):
    """A hand-written file with MCOMP=2 and IRCHCONC 1 0 has one
    concentration array; FloPy keeps the flags and writes the same file."""
    ml = _ml(function_tmpdir, nper=2)
    rch, out, again = _load_write_reload(function_tmpdir, _H_F41, ml)
    assert rch.irchconc == [1, 0]
    assert _controls(out) == _controls(_H_F41)
    assert again.irchconc == [1, 0]
    assert _values(again) == _values(rch)
    assert sorted(rch.rchconc.transient_3ds) == [0]
    np.testing.assert_allclose(rch.rchconc.transient_3ds[0][0].array.ravel(), _C5)
    np.testing.assert_array_equal(rch.rchconc.transient_3ds[0][1].array, 0.0)


def test_rch_irchconc_0_1_api(function_tmpdir):
    """irchconc=[0, 1]: IRCHCONC 0 1 and only the component 2 array."""
    ml = _ml(function_tmpdir, nper=1)
    rch = MfUsgRch(
        ml,
        nrchop=3,
        ipakcb=50,
        rech=1e-3,
        iconc=1,
        irchconc=[0, 1],
        rchconc={0: [0.0, 7.0]},
    )
    out = _write(rch, function_tmpdir)
    assert _controls(out) == [
        ["3", "50", "CONC"],
        ["0", "1"],
        ["1", "INCONC"],
        ["ARRAY"],
        ["ARRAY"],
    ]
    again = _reload(function_tmpdir / "out.rch", ml)
    assert again.irchconc == [0, 1]
    np.testing.assert_array_equal(again.rchconc.transient_3ds[0][1].array, 7.0)
    np.testing.assert_array_equal(again.rchconc.transient_3ds[0][0].array, 0.0)


def test_rch_irchconc_default_all_components(function_tmpdir):
    """Without irchconc every component is read (IRCHCONC 1 1), as before."""
    ml = _ml(function_tmpdir, nper=1)
    rch = MfUsgRch(ml, nrchop=3, ipakcb=50, rech=1e-3, iconc=1, rchconc={0: [3.0, 4.0]})
    assert list(rch.irchconc) == [1, 1]
    out = _write(rch, function_tmpdir)
    assert _controls(out)[1:3] == [["1", "1"], ["1", "INCONC"]]
    assert len(_controls(out)) == 6  # header, flags, period, RECH, 2 arrays


_BAD_CONC = {
    "flag_2": {"iconc": 1, "irchconc": [2, 0], "rchconc": 1.0},
    "flag_negative": {"iconc": 1, "irchconc": [-1, 1], "rchconc": 1.0},
    "flags_short": {"iconc": 1, "irchconc": [1], "rchconc": 1.0},
    "flags_long": {"iconc": 1, "irchconc": [1, 0, 1], "rchconc": 1.0},
    "flag_float": {"iconc": 1, "irchconc": [1.0, 0], "rchconc": 1.0},
    "flag_bool": {"iconc": 1, "irchconc": [True, 0], "rchconc": 1.0},
    "flags_text": {"iconc": 1, "irchconc": "10", "rchconc": 1.0},
    "irchconc_without_conc": {"iconc": 0, "irchconc": [1, 0]},
    "rchconc_without_conc": {"iconc": 0, "rchconc": 1.0},
    "no_period1_conc": {"iconc": 1, "irchconc": [1, 0]},
    "conc_from_period2": {"iconc": 1, "rchconc": {1: [1.0, 0.0]}},
}


@pytest.mark.parametrize("case", list(_BAD_CONC))
def test_rch_conc_rejected_before_writing(function_tmpdir, case):
    """Flags other than 0/1, a wrong number of flags, concentration data
    without CONC and a period 1 without concentrations are rejected by the
    constructor (ValueError); the package is not added to the model."""
    ml = _ml(function_tmpdir, nper=2)
    with pytest.raises(ValueError, match=r"irchconc|rchconc|IRCHCONC|RCHCONC"):
        MfUsgRch(ml, nrchop=3, ipakcb=50, rech=1e-3, **_BAD_CONC[case])
    assert ml.get_package("RCH") is None


def test_rch_conc_all_flags_zero_needs_no_concentration(function_tmpdir):
    """IRCHCONC 0 0: no concentration array is read, so none is needed."""
    ml = _ml(function_tmpdir, nper=2)
    rch = MfUsgRch(ml, nrchop=3, ipakcb=50, rech=1e-3, iconc=1, irchconc=[0, 0])
    out = _write(rch, function_tmpdir)
    assert _controls(out) == [["3", "50", "CONC"], ["0", "0"], ["1"], ["ARRAY"], ["-1"]]


def test_rch_conc_needs_bct_counts(function_tmpdir):
    """CONC reads MCOMPT flags of the BCT package: without BCT counts the
    constructor raises ValueError."""
    ml = _ml(function_tmpdir, nper=1, mcomp=0)
    with pytest.raises(ValueError, match="BCT"):
        MfUsgRch(ml, nrchop=3, rech=1e-3, iconc=1, rchconc=1.0)
    assert ml.get_package("RCH") is None


@pytest.mark.parametrize("case", ["hand_counts", "bct_removed"])
def test_rch_conc_counts_set_by_hand(function_tmpdir, case):
    """Counts set by hand on a model without BCT (itrnsp 0) are used; a BCT
    package that was removed leaves stale counts (ValueError)."""
    if case == "hand_counts":
        ml = _ml(function_tmpdir, nper=1, mcomp=0)
        ml.mcomp = 2
        rch = MfUsgRch(ml, nrchop=3, ipakcb=50, rech=1e-3, iconc=1, rchconc=1.0)
        assert _controls(_write(rch, function_tmpdir))[1] == ["1", "1"]
    else:
        ml = _ml(function_tmpdir, nper=1)
        ml.remove_package("BCT")
        with pytest.raises(ValueError, match="BCT"):
            MfUsgRch(ml, nrchop=3, rech=1e-3, iconc=1, rchconc=1.0)
        assert ml.get_package("RCH") is None


@pytest.mark.parametrize("case", ["bct_replaced", "flags_changed", "conc_turned_on"])
def test_rch_conc_checked_again_at_write(function_tmpdir, case):
    """Counts or attributes changed after construction are checked before
    the file is opened (ValueError, no file)."""
    ml = _ml(function_tmpdir, nper=1)
    rch = MfUsgRch(ml, nrchop=3, ipakcb=50, rech=1e-3, iconc=1, rchconc=1.0)
    if case == "bct_replaced":
        ml.remove_package("BCT")
        MfUsgBct(ml, itrnsp=1, mcomp=3, idisp=0, prsity=0.2, conc=0.0)
    elif case == "flags_changed":
        rch.irchconc = [1, 2]
    else:
        ml.remove_package("RCH")
        rch = MfUsgRch(ml, nrchop=3, ipakcb=50, rech=1e-3)
        rch.iconc = 1
        rch.irchconc = [1, 0]
    rch.fn_path = str(function_tmpdir / "out.rch")
    with pytest.raises(ValueError, match=r"MCOMPT|irchconc|IRCHCONC|RCHCONC"):
        rch.write_file(check=False)
    assert not (function_tmpdir / "out.rch").exists()


def test_rch_load_closes_file_when_counts_missing(function_tmpdir, monkeypatch):
    """A CONC file loaded without BCT counts raises ValueError and the file
    opened by load is closed."""
    ml = _ml(function_tmpdir, nper=2, mcomp=0)
    path = function_tmpdir / "t.rch"
    path.write_text(_H_F41)
    opened = []
    real_open = builtins.open

    def spy(*args, **kwargs):
        fh = real_open(*args, **kwargs)
        opened.append(fh)
        return fh

    monkeypatch.setattr(builtins, "open", spy)
    with pytest.raises(ValueError, match="BCT"):
        MfUsgRch.load(str(path), ml, check=False)
    monkeypatch.undo()
    mine = [fh for fh in opened if Path(fh.name) == path]
    assert mine and all(fh.closed for fh in opened)


# --------------------------------------------------------------- INCONC
_H_ONCE = (
    "# RCH\n 3 50 CONC\n 1 1\n 1 INCONC\n"
    + _blk(_R5)
    + _blk(_C5)
    + _blk(_D5)
    + " 1\n"
    + _blk(_R5)
    + " 1\n"
    + _blk(_R5)
)


def test_rch_inconc_given_once_load_write(function_tmpdir):
    """INCONC only in period 1 (USG-T reuses the concentrations): FloPy
    keeps the concentrations of period 1 only and writes INCONC once."""
    ml = _ml(function_tmpdir, nper=3)
    rch, out, again = _load_write_reload(function_tmpdir, _H_ONCE, ml)
    assert sorted(rch.rchconc.transient_3ds) == [0]
    assert _controls(out) == _controls(_H_ONCE)
    assert _values(again) == _values(rch)


def test_rch_inconc_once_api(function_tmpdir):
    """rchconc given once (a scalar, period 1) with three periods: INCONC in
    period 1 only, and the file loads again."""
    ml = _ml(function_tmpdir, nper=3)
    rch = MfUsgRch(ml, nrchop=3, ipakcb=50, rech=1e-3, iconc=1, rchconc=2.0)
    out = _write(rch, function_tmpdir)
    heads = [c for c in _controls(out) if c != ["ARRAY"]]
    assert heads == [["3", "50", "CONC"], ["1", "1"], ["1", "INCONC"], ["-1"], ["-1"]]
    again = _reload(function_tmpdir / "out.rch", ml)
    assert sorted(again.rchconc.transient_3ds) == [0]
    np.testing.assert_array_equal(again.rchconc.transient_3ds[0][1].array, 2.0)


# ---------------------------------------------------- DISU NRCHOP=2 IRCH
_H_IRCH = (
    "# RCH\n 2 50\n 5\n 1 5\n"
    + _blk(_R5)
    + _iblk([1, 2, 3, 4, 5])
    + " 1 -1\n"
    + _blk(_R5)
    + " 1 5\n"
    + _blk(_R5)
    + _iblk([5, 4, 3, 2, 1])
)


def test_rch_disu_nrchop2_irch_reuse_and_change_load_write(function_tmpdir):
    """IRCH in period 1, reused (INIRCH -1) in period 2 and changed in
    period 3: FloPy keeps the reuse and writes the same records."""
    ml = _ml(function_tmpdir, nper=3, mcomp=0)
    rch, out, again = _load_write_reload(function_tmpdir, _H_IRCH, ml)
    assert sorted(rch.irch.transient_2ds) == [0, 2]
    assert _controls(out) == _controls(_H_IRCH)
    assert _values(again) == _values(rch)
    assert again.irch.transient_2ds[2].array.ravel().tolist() == [4, 3, 2, 1, 0]


def test_rch_disu_nrchop2_irch_period1_only_api(function_tmpdir):
    """irch given for period 1 only: later periods write INIRCH -1 and no
    IRCH array (before: the node count without the array)."""
    ml = _ml(function_tmpdir, nper=3, mcomp=0)
    rch = MfUsgRch(
        ml,
        nrchop=2,
        ipakcb=50,
        rech=np.full((1, 5), 1e-3),
        irch=np.array([[0, 1, 2, 3, 4]]),
    )
    out = _write(rch, function_tmpdir)
    heads = [c for c in _controls(out) if c != ["ARRAY"]]
    assert heads == [["2", "50"], ["5"], ["1", "5"], ["-1", "-1"], ["-1", "-1"]]
    again = _reload(function_tmpdir / "out.rch", ml)
    assert sorted(again.irch.transient_2ds) == [0]


_H_IRCH_SIZE = (
    "# RCH\n 2 50\n 5\n 1 5\n"
    + _blk(_R5)
    + _iblk([1, 2, 3, 4, 5])
    + " 1 4\n"
    + _blk(_R5[:4])
    + _iblk([2, 3, 4, 5])
)


def test_rch_disu_nrchop2_node_count_changes_load_write(function_tmpdir):
    """A period with fewer IRCH nodes (INIRCH 4 <= MXNDRCH 5) is kept."""
    ml = _ml(function_tmpdir, nper=2, mcomp=0)
    rch, out, again = _load_write_reload(function_tmpdir, _H_IRCH_SIZE, ml)
    assert _controls(out) == _controls(_H_IRCH_SIZE)
    assert _values(again) == _values(rch)


def _u2d(ml, values, dtype=np.int32):
    arr = np.array([values], dtype=dtype)
    return Util2d(ml, arr.shape, dtype, arr, "irch")


@pytest.mark.parametrize(
    "case", ["no_irch_period1", "rech_size", "count_changes_rech_reused"]
)
def test_rch_disu_nrchop2_rejected_before_writing(function_tmpdir, case):
    """IRCH missing in period 1, a RECH array whose size is not the IRCH node
    count, and a new node count while RECH is reused are rejected before the
    file is opened."""
    ml = _ml(function_tmpdir, nper=2, mcomp=0)
    rech = np.full((1, 5), 1e-3)
    irch = {0: np.array([[0, 1, 2, 3, 4]])}
    if case == "no_irch_period1":
        irch = {1: np.array([[0, 1, 2, 3, 4]])}
    elif case == "rech_size":
        rech = np.full((1, 4), 1e-3)
    else:
        irch = {0: _u2d(ml, [0, 1, 2, 3, 4]), 1: _u2d(ml, [0, 1, 2, 3])}
    rch = MfUsgRch(ml, nrchop=2, ipakcb=50, rech=rech, irch=irch)
    rch.fn_path = str(function_tmpdir / "out.rch")
    with pytest.raises(ValueError, match=r"IRCH|irch"):
        rch.write_file(check=False)
    assert not (function_tmpdir / "out.rch").exists()


_H_CONC_REUSE = (
    "# RCH\n 2 50 CONC\n 5\n 1 0\n 1 5 INCONC\n"
    + _blk(_R5)
    + _iblk([1, 2, 3, 4, 5])
    + _blk(_C5)
    + " -1 -1\n"
)


@pytest.mark.parametrize("source", ["load", "api"])
def test_rch_conc_nrchop2_irch_reuse_not_written(function_tmpdir, source):
    """CONC with a period that reuses IRCH: USG-T adds no recharge mass in
    that period (gwt2bndsu1.f :3508), so FloPy raises NotImplementedError
    before opening the file; the file still loads."""
    ml = _ml(function_tmpdir, nper=2)
    if source == "load":
        rch = _load(function_tmpdir, _H_CONC_REUSE, ml)
    else:
        rch = MfUsgRch(
            ml,
            nrchop=2,
            ipakcb=50,
            rech=np.full((1, 5), 1e-3),
            irch=np.array([[0, 1, 2, 3, 4]]),
            iconc=1,
            irchconc=[1, 0],
            rchconc=1.0,
        )
    rch.fn_path = str(function_tmpdir / "out.rch")
    with pytest.raises(NotImplementedError, match="irch"):
        rch.write_file(check=False)
    assert not (function_tmpdir / "out.rch").exists()


_H_RTS_LAYER2 = (
    "# RCH\n 2 50 RTS 1\n 10\n 1 5 INRCHZONES 1\n"
    + _blk([0.0] * 5)
    + _iblk([6, 7, 8, 9, 10])
    + _iblk([1] * 5)
)


@pytest.mark.parametrize("case", ["api_layer2", "load_mxndrch10", "api_permuted"])
def test_rch_rts_nrchop2_irch_nodes(function_tmpdir, case):
    """RTS with DISU NRCHOP=2: USG-T takes the zone of node IRCH(NN) from
    IZNRCH, read in positions 1..NIRCH only. A node above NIRCH is not
    written (NotImplementedError; the zone is never read, also with MXNDRCH
    above the node); nodes permuted within 1..NIRCH are written."""
    ml = _ml(function_tmpdir, nper=1, mcomp=0, nlay=2)
    if case == "load_mxndrch10":
        rch = _load(function_tmpdir, _H_RTS_LAYER2, ml)
    else:
        nodes = [5, 6, 7, 8, 9] if case == "api_layer2" else [4, 3, 2, 1, 0]
        rch = MfUsgRch(
            ml,
            nrchop=2,
            ipakcb=50,
            rech=np.zeros((1, 5)),
            irch=np.array([nodes]),
            mxrtzones=1,
            iznrch=np.ones((1, 5), dtype=int),
        )
    rch.fn_path = str(function_tmpdir / "out.rch")
    if case == "api_permuted":
        rch.write_file(check=False)
        again = _reload(rch.fn_path, ml)
        assert again.irch.transient_2ds[0].array.ravel().tolist() == [4, 3, 2, 1, 0]
    else:
        with pytest.raises(NotImplementedError, match="IZNRCH"):
            rch.write_file(check=False)
        assert not (function_tmpdir / "out.rch").exists()


# ------------------------------------------------------------- keywords
_KW_CASES = {
    "conc_lower": (
        " 3 50 conc\n 1 1\n 0 inconc\n",
        " 3 50 CONC\n 1 1\n 0 INCONC\n",
        _blk(_R5) + _blk(_C5) + _blk(_D5),
    ),
    "concentration_mixed": (
        " 3 50 Concentration\n 1 1\n 0 InConc\n",
        " 3 50 CONC\n 1 1\n 0 INCONC\n",
        _blk(_R5) + _blk(_C5) + _blk(_D5),
    ),
    "inconc_in_comment": (
        " 3 50 CONC\n 1 1\n 0 # inconc\n",
        " 3 50 CONC\n 1 1\n 0 INCONC\n",
        _blk(_R5) + _blk(_C5) + _blk(_D5),
    ),
    "seepelev_lower": (
        " 3 50 seepelev\n 0 inselev 1\n",
        " 3 50 SEEPELEV\n 0 INSELEV 1\n",
        _blk(_R5) + _blk(_C5),
    ),
    "rts_lower": (
        " 1 0 rts 2\n 0 inrchzones 2\n",
        " 1 0 RTS 2\n 0 INRCHZONES 2\n",
        _blk(_R5) + _iblk([0, 1, 2, 1, 0]),
    ),
    "not_whole_words": (
        " 3 50 XCONC CONCX\n 1 XINCONC\n",
        " 3 50\n 1\n",
        _blk(_R5),
    ),
}


@pytest.mark.parametrize("case", list(_KW_CASES))
def test_rch_keywords_any_case(function_tmpdir, case):
    """Option and period keywords are whole words in any case, also after
    '#' (URWORD scans the whole line); a word that only contains one is not
    a keyword."""
    lower, upper, arrays = _KW_CASES[case]
    ml = _ml(function_tmpdir, nper=1)
    a = _load(function_tmpdir, "# RCH\n" + lower + arrays, ml, "a.rch")
    b = _load(function_tmpdir, "# RCH\n" + upper + arrays, ml, "b.rch")
    assert (a.iconc, a.seepelev, a.mxrtzones) == (b.iconc, b.seepelev, b.mxrtzones)
    assert _values(a) == _values(b)
    assert _write(a, function_tmpdir, "a_out.rch") == _write(
        b, function_tmpdir, "b_out.rch"
    )


def test_rch_written_comments_have_no_keywords(function_tmpdir):
    """The comments FloPy writes on the option line and the period lines
    contain no RCH keyword as a whole word (USG-T reads words after '#')."""
    ml = _ml(function_tmpdir, nper=3)
    rch = MfUsgRch(
        ml,
        nrchop=2,
        ipakcb=50,
        rech={0: np.full((1, 5), 1e-3), 2: np.full((1, 5), 2e-3)},
        irch={
            0: np.array([[0, 1, 2, 3, 4]]),
            1: np.array([[0, 1, 2, 3, 4]]),
            2: np.array([[0, 1, 2, 3, 4]]),
        },
        seepelev=1,
        selev={0: np.full((1, 5), 9.0)},
        mxrtzones=2,
        iznrch={0: np.array([[0, 1, 2, 1, 0]])},
        iconc=1,
        irchconc=[1, 0],
        rchconc=1.0,
    )
    lines = _write(rch, function_tmpdir).splitlines()
    records = [lines[1]] + [ln for ln in lines if "Stress period" in ln]
    assert len(records) == 4
    for line in records:
        comment = line.split("#", 1)[1] if "#" in line else ""
        assert not _KEYWORDS & set(comment.upper().replace(",", " ").split()), line


# ------------------------------------------------------- engine controls
def _exe_model(ws, nper, mcomp, nlay=1):
    """Flow (and transport with mcomp > 0) model written by FloPy, without
    RCH; the RCH file is added to the name file as t.rch. DISU with IDISP=1
    and zero dispersion: the local 2.8.1 and 2.7.0 builds stop with a
    segmentation fault for MCOMP=2 with IDISP=0, with or without RCH."""
    ml = _ml(ws, nper=nper, mcomp=mcomp, nlay=nlay, idisp=1)
    MfUsgBas(ml, ibound=1, strt=10.0)
    MfUsgLpf(ml, laytyp=0, hk=1.0, vka=1.0, ss=1.0e-3, ipakcb=50)
    MfUsgSms(ml, linmeth=1)
    ModflowChd(ml, stress_period_data={0: [[0, 10.0, 10.0]]})
    act = ["save head", "save budget", "print budget"]
    if mcomp:
        act.append("save concentration")
    MfUsgOc(ml, stress_period_data={(k, s): act for k in range(nper) for s in range(2)})
    ml.write_input()
    with open(Path(ws) / "t.nam", "a") as f:
        f.write("RCH               18  t.rch\n")
    return ml


def _exe_run(ws):
    ws = Path(ws)
    ok, buff = run_model(USGT_EXE, "t.nam", model_ws=str(ws), silent=True, report=True)
    (ws / "engine_stdout.txt").write_text("\n".join(str(b) for b in buff))
    assert ok and any("normal termination" in str(b).lower() for b in buff)
    lst = (ws / "t.list").read_text(errors="replace")
    out = {
        "hds": (ws / "t.hds").read_bytes(),
        "cbc": (ws / "t.cbc").read_bytes(),
        "rch": re.findall(r"^\s+RECHARGE =\s+(\S+)\s+RECHARGE =\s+(\S+)", lst, re.M),
        "rch_mass": re.findall(r"^\s+RCH MASS FLUX\S* =\s+(\S+)", lst, re.M),
    }
    if (ws / "t.con").exists():
        out["con"] = (ws / "t.con").read_bytes()
    assert out["rch"]
    return out


def _con_species(path, mcomp):
    """Concentrations of an unstructured binary file, one list per species
    (records alternate by species within each output time)."""
    b = Path(path).read_bytes()
    recs, i = [], 0
    while i < len(b):
        n0, n1 = np.frombuffer(b, "<i4", 2, i + 32)
        n = int(n1 - n0 + 1)
        recs.append(np.frombuffer(b, "<f4", n, i + 44).copy())
        i += 44 + 4 * n
    return [recs[k::mcomp] for k in range(mcomp)]


def _f41_text(flags, conc):
    return (
        f"# RCH\n 3 50 CONC\n {flags}\n 1 INCONC\n" + _blk(_R5) + _blk(conc) + " -1\n"
    )


@requires_usgt
@pytest.mark.parametrize("source", ["load_write", "api"])
def test_usgt_exe_rch_irchconc_matches_hand(function_tmpdir, source):
    """MCOMP=2 with IRCHCONC 1 0: FloPy's file runs like the hand-written
    one (heads, budget, concentrations); component 2 gets no recharge mass.
    Changing the component 1 concentration, or reading component 2 instead
    (0 1), changes the concentrations."""
    ctrl_ws = function_tmpdir / "control"
    _exe_model(ctrl_ws, nper=2, mcomp=2)
    (ctrl_ws / "t.rch").write_text(_f41_text("1 0", _C5))
    control = _exe_run(ctrl_ws)
    sp1, sp2 = _con_species(ctrl_ws / "t.con", 2)
    assert max(float(r.max()) for r in sp1) > 0.0
    assert all(np.all(r == 0.0) for r in sp2)

    ws = function_tmpdir / source
    ml = _exe_model(ws, nper=2, mcomp=2)
    if source == "load_write":
        rch = _load(ws, _f41_text("1 0", _C5), ml)
    else:
        rch = MfUsgRch(
            ml,
            nrchop=3,
            ipakcb=50,
            rech=np.array([_R5]),
            iconc=1,
            irchconc=[1, 0],
            rchconc={0: [np.array([_C5]), 0.0]},
        )
    _write(rch, ws, "t.rch")
    assert _exe_run(ws) == control

    for name, text in (
        ("conc_changed", _f41_text("1 0", _D5)),
        ("flags_0_1", _f41_text("0 1", _C5)),
    ):
        other_ws = function_tmpdir / name
        _exe_model(other_ws, nper=2, mcomp=2)
        (other_ws / "t.rch").write_text(text)
        other = _exe_run(other_ws)
        assert other["hds"] == control["hds"]
        assert other["con"] != control["con"]
    o1, o2 = _con_species(function_tmpdir / "flags_0_1" / "t.con", 2)
    assert all(np.all(r == 0.0) for r in o1)
    assert max(float(r.max()) for r in o2) > 0.0


def _irch_text(p2, p3):
    """DISU 2 layers x 5 nodes, NRCHOP=2: IRCH on layer 2 in period 1, then
    the period 2 and 3 records given."""
    return "# RCH\n 2 50\n 5\n 1 5\n" + _blk(_R5) + _iblk([6, 7, 8, 9, 10]) + p2 + p3


_IRCH_P3 = " -1 5\n" + _iblk([1, 2, 3, 4, 5])


@requires_usgt
@pytest.mark.parametrize("source", ["load_write", "api"])
def test_usgt_exe_rch_disu_nrchop2_irch_reuse_and_change(function_tmpdir, source):
    """DISU NRCHOP=2: IRCH on layer 2, reused in period 2 (INIRCH -1) and
    moved to layer 1 in period 3. FloPy's file runs like the hand-written
    one; reusing IRCH in period 3 instead changes the heads."""
    ctrl_ws = function_tmpdir / "control"
    _exe_model(ctrl_ws, nper=3, mcomp=0, nlay=2)
    (ctrl_ws / "t.rch").write_text(_irch_text(" -1 -1\n", _IRCH_P3))
    control = _exe_run(ctrl_ws)

    ws = function_tmpdir / source
    ml = _exe_model(ws, nper=3, mcomp=0, nlay=2)
    if source == "load_write":
        rch = _load(ws, _irch_text(" -1 -1\n", _IRCH_P3), ml)
    else:
        rch = MfUsgRch(
            ml,
            nrchop=2,
            ipakcb=50,
            rech=np.array([_R5]),
            irch={0: np.array([[5, 6, 7, 8, 9]]), 2: np.array([[0, 1, 2, 3, 4]])},
        )
    text = _write(rch, ws, "t.rch")
    assert [c for c in _controls(text) if c != ["ARRAY"]][3] == (
        ["-1", "-1"] if source == "api" else ["1", "-1"]
    )
    assert _exe_run(ws) == control

    other_ws = function_tmpdir / "p3_reused"
    _exe_model(other_ws, nper=3, mcomp=0, nlay=2)
    (other_ws / "t.rch").write_text(_irch_text(" -1 -1\n", " -1 -1\n"))
    assert _exe_run(other_ws)["hds"] != control["hds"]


# ------------------------------------------- structured NRCHOP=2 INIRCH
def _dis(ws, nper=2, mcomp=1, nrow=2, ncol=3, nlay=1):
    """Structured model of nlay layers (BCT with mcomp components when
    mcomp > 0, IDISP=0)."""
    ws = Path(ws)
    ws.mkdir(parents=True, exist_ok=True)
    ml = MfUsg(modelname="t", model_ws=str(ws), structured=True, exe_name=USGT_EXE)
    ModflowDis(
        ml,
        nlay=nlay,
        nrow=nrow,
        ncol=ncol,
        nper=nper,
        perlen=10.0,
        nstp=2,
        steady=False,
        delr=10.0,
        delc=10.0,
        top=30.0,
        botm=[30.0 * (nlay - k - 1) / nlay for k in range(nlay)],
    )
    if mcomp:
        MfUsgBct(ml, itrnsp=1, mcomp=mcomp, idisp=0, prsity=0.2, conc=0.0)
    return ml


_R6 = [1.0e-3, 2.0e-3, 3.0e-3, 4.0e-3, 5.0e-3, 6.0e-3]
_C6 = [5.0, 6.0, 7.0, 8.0, 9.0, 10.0]


def _s_conc(n1, n2):
    """Structured 2 x 3, CONC (1 component), IRCH read in both periods with
    INIRCH n1 and n2."""
    return (
        f"# RCH\n 2 50 CONC\n 1\n 1 {n1} INCONC\n"
        + _blk(_R6)
        + "CONSTANT 1\n"
        + _blk(_C6)
        + f" -1 {n2}\nCONSTANT 1\n"
    )


_S_FLOW = "# RCH\n 2 50\n 1 6\n" + _blk(_R6) + _iblk([1] * 6) + " -1 -1\n"


def _inirch_of(text):
    """INIRCH of the stress period records of a file written by FloPy."""
    lines = [ln for ln in text.splitlines() if "Stress period" in ln]
    return [ln.split("#")[0].split()[1] for ln in lines]


def _write_quiet(rch, ws, name="out.rch"):
    """_write with the INIRCH warning turned into an error."""
    with warnings.catch_warnings():
        warnings.filterwarnings("error", message=".*INIRCH")
        return _write(rch, ws, name)


_LOW = r"first INIRCH cells.*without mass.*same count.*inirch_ncpl=True"


_KEPT = {
    "conc_1": (1, _s_conc(1, 1), ["1", "1"]),
    "conc_ncpl": (1, _s_conc(6, 6), ["6", "6"]),
    "conc_mixed": (1, _s_conc(6, 3), ["6", "3"]),
    "flow_reuse": (0, _S_FLOW, ["6", "-1"]),
}


@pytest.mark.parametrize("case", list(_KEPT))
def test_rch_structured_inirch_kept_on_load_write(function_tmpdir, case):
    """Load keeps the INIRCH of every period and the written file has the
    same INIRCH (and data); with CONC an INIRCH below NROW*NCOL warns, with
    the periods named."""
    mcomp, text, inirch = _KEPT[case]
    ml = _dis(function_tmpdir, mcomp=mcomp)
    rch = _load(function_tmpdir, text, ml)
    if case in ("conc_1", "conc_mixed"):
        periods = "1 (INIRCH 1), 2 (INIRCH 1)" if case == "conc_1" else "2 (INIRCH 3)"
        with pytest.warns(UserWarning, match=_LOW) as record:
            out = _write(rch, function_tmpdir)
        assert periods in str(record[0].message)
    else:
        out = _write_quiet(rch, function_tmpdir)
    assert _inirch_of(out) == inirch
    assert rch.inirch == {k: int(v) for k, v in enumerate(inirch)}
    assert rch.inirch_ncpl is False
    again = _reload(function_tmpdir / "out.rch", ml)
    assert again.inirch == rch.inirch
    assert _values(again) == _values(rch)


@pytest.mark.parametrize("conc", [1, 0])
def test_rch_structured_inirch_api_default(function_tmpdir, conc):
    """API without inirch: with CONC the periods that read irch write
    NROW*NCOL (no warning); without CONC they write 1 as before."""
    ml = _dis(function_tmpdir, mcomp=conc)
    kwargs = {"iconc": 1, "rchconc": 1.0} if conc else {}
    rch = MfUsgRch(
        ml,
        nrchop=2,
        ipakcb=50,
        rech=np.array([_R6[:3], _R6[3:]]),
        irch={0: 0, 1: 0},
        **kwargs,
    )
    out = _write_quiet(rch, function_tmpdir)
    assert _inirch_of(out) == (["6", "6"] if conc else ["1", "1"])


def test_rch_structured_inirch_api_low_warns(function_tmpdir):
    """API with CONC and inirch below NROW*NCOL: written as given, warns."""
    ml = _dis(function_tmpdir)
    rch = MfUsgRch(
        ml,
        nrchop=2,
        rech=1e-3,
        irch={0: 0, 1: 0},
        iconc=1,
        rchconc=1.0,
        inirch={0: 6, 1: 2},
    )
    with pytest.warns(UserWarning, match=_LOW):
        out = _write(rch, function_tmpdir)
    assert _inirch_of(out) == ["6", "2"]


@pytest.mark.parametrize("source", ["constructor", "loaded"])
def test_rch_structured_inirch_ncpl_override(function_tmpdir, source):
    """inirch_ncpl=True writes NROW*NCOL in every period that reads irch,
    without a warning, in the constructor or set on a loaded package."""
    ml = _dis(function_tmpdir)
    if source == "constructor":
        rch = MfUsgRch(
            ml,
            nrchop=2,
            rech=1e-3,
            irch={0: 0, 1: 0},
            iconc=1,
            rchconc=1.0,
            inirch={0: 1, 1: 1},
            inirch_ncpl=True,
        )
    else:
        rch = _load(function_tmpdir, _s_conc(1, 3), ml)
        rch.inirch_ncpl = True
    out = _write_quiet(rch, function_tmpdir)
    assert _inirch_of(out) == ["6", "6"]
    assert _reload(function_tmpdir / "out.rch", ml).inirch == {0: 6, 1: 6}


_S_CONC_REUSE = (
    "# RCH\n 2 50 CONC\n 1\n 1 6 INCONC\n"
    + _blk(_R6)
    + "CONSTANT 1\n"
    + _blk(_C6)
    + " -1 -1\n"
)


@pytest.mark.parametrize("source", ["load", "api"])
def test_rch_structured_conc_irch_reuse_not_written(function_tmpdir, source):
    """Structured CONC with a period that reuses IRCH: USG-T adds no
    recharge mass in that period (gwt2bndsu1.f :3508; probe P7 with both
    engines), so FloPy raises NotImplementedError before opening the file,
    also with inirch_ncpl; the file still loads."""
    ml = _dis(function_tmpdir)
    if source == "load":
        rch = _load(function_tmpdir, _S_CONC_REUSE, ml)
    else:
        rch = MfUsgRch(ml, nrchop=2, rech=1e-3, irch=0, iconc=1, rchconc=1.0)
    rch.fn_path = str(function_tmpdir / "out.rch")
    for override in (False, True):
        rch.inirch_ncpl = override
        with pytest.raises(NotImplementedError, match="every stress period"):
            rch.write_file(check=False)
        assert not (function_tmpdir / "out.rch").exists()
    if source == "load":
        assert rch.inirch == {0: 6, 1: -1}


_BAD_INIRCH = {
    "disu": ("disu", {"inirch": {0: 5}}, "structured"),
    "nrchop3": ("dis", {"nrchop": 3, "inirch_ncpl": True}, "structured"),
    "ncpl_not_bool": ("dis", {"inirch_ncpl": "yes"}, "True or False"),
    "not_int": ("dis", {"inirch": {0: 1.5}}, "not an integer"),
    "negative_with_irch": ("dis", {"inirch": {0: -1}}, "does not read it"),
    "positive_without_irch": ("dis", {"inirch": {1: 3}}, "needs irch data"),
    "above_ncpl": ("dis", {"inirch": {0: 7, 1: 6}}, "above NROW"),
}


@pytest.mark.parametrize("case", list(_BAD_INIRCH))
def test_rch_structured_inirch_rejected(function_tmpdir, case):
    """Wrong inirch/inirch_ncpl: ValueError in the constructor or, for the
    checks that need the irch data, at write time before opening the file;
    also when set on the package after it was built."""
    grid, kwargs, match = _BAD_INIRCH[case]
    if grid == "disu":
        ml = _ml(function_tmpdir, nper=1, mcomp=1)
        base = {
            "nrchop": 2,
            "rech": np.full((1, 5), 1e-3),
            "irch": np.array([[0, 1, 2, 3, 4]]),
        }
    else:
        ml = _dis(function_tmpdir)
        base = {
            "nrchop": 2,
            "rech": 1e-3,
            "irch": {0: 0} if "positive" in case else {0: 0, 1: 0},
        }
    base.update(iconc=1, rchconc=1.0)
    args = {**base, **kwargs}
    if case in ("negative_with_irch", "positive_without_irch", "above_ncpl"):
        if case == "positive_without_irch":
            args.update(iconc=0, rchconc=None)
        rch = MfUsgRch(ml, **args)
        rch.fn_path = str(function_tmpdir / "out.rch")
        with pytest.raises(ValueError, match=match):
            rch.write_file(check=False)
        assert not (function_tmpdir / "out.rch").exists()
    else:
        with pytest.raises(ValueError, match=match):
            MfUsgRch(ml, **args)
        rch = MfUsgRch(ml, **base)
        name = next(iter(kwargs.keys() - {"nrchop"}))
        if "nrchop" in kwargs:
            rch.nrchop = 3
        setattr(rch, name, kwargs[name])
        rch.fn_path = str(function_tmpdir / "out.rch")
        with pytest.raises(ValueError, match=match):
            rch.write_file(check=False)
        assert not (function_tmpdir / "out.rch").exists()


def _s5_text(inirch):
    """Probe P2 file: structured 1 x 5, NRCHOP=2, IRCH 1, CONC, one period."""
    return (
        f"# RCH\n 2 50 CONC\n 1\n 0 {inirch} INCONC\n"
        + _blk([1.0e-3] * 5)
        + "CONSTANT 1\n"
        + _blk([1.0] * 5)
    )


def _dis_exe_model(ws):
    """Probe P2 model: structured 1 x 5, one period of 2 steps, CHD in cell
    1, MCOMP=1 with IDISP=0; RCH added to the name file as t.rch."""
    ml = _dis(ws, nper=1, mcomp=1, nrow=1, ncol=5)
    MfUsgBas(ml, ibound=1, strt=10.0)
    MfUsgLpf(ml, laytyp=0, hk=1.0, vka=1.0, ss=1.0e-3, ipakcb=50)
    MfUsgSms(ml, linmeth=1)
    ModflowChd(ml, stress_period_data={0: [[0, 0, 0, 10.0, 10.0]]})
    act = ["save head", "save budget", "print budget", "save concentration"]
    MfUsgOc(ml, stress_period_data={(0, s): act for s in range(2)})
    ml.write_input()
    with open(Path(ws) / "t.nam", "a") as f:
        f.write("RCH               18  t.rch\n")
    return ml


@requires_usgt
@pytest.mark.parametrize(
    "source", ["load_write_1", "load_write_ncpl", "api", "override"]
)
def test_usgt_exe_rch_structured_inirch(function_tmpdir, source):
    """Structured 1 x 5 with CONC: the hand files with INIRCH 1 and 5 run
    with the same flow, but only INIRCH 5 adds the full recharge mass.
    Loading and writing either file keeps its result; the API file and the
    INIRCH 1 file written with inirch_ncpl=True run like INIRCH 5."""
    runs = {}
    for name, inirch in (("hand_1", 1), ("hand_ncpl", 5)):
        ws = function_tmpdir / name
        _dis_exe_model(ws)
        (ws / "t.rch").write_text(_s5_text(inirch))
        runs[name] = _exe_run(ws)
    assert runs["hand_1"]["hds"] == runs["hand_ncpl"]["hds"]
    assert runs["hand_1"]["cbc"] == runs["hand_ncpl"]["cbc"]
    assert runs["hand_1"]["rch"] == runs["hand_ncpl"]["rch"]
    # conc 1.0: full mass equals the recharge volume, INIRCH 1 (CHD cell) none
    volume = [cum for cum, _ in runs["hand_ncpl"]["rch"]]
    assert runs["hand_ncpl"]["rch_mass"] == volume
    assert max(float(v) for v in volume) > 0.0
    assert all(float(m) == 0.0 for m in runs["hand_1"]["rch_mass"])

    ws = function_tmpdir / source
    ml = _dis_exe_model(ws)
    if source == "api":
        rch = MfUsgRch(
            ml,
            nrchop=2,
            ipakcb=50,
            rech=np.full((1, 5), 1.0e-3),
            irch=0,
            iconc=1,
            rchconc=1.0,
        )
    else:
        rch = _load(ws, _s5_text(5 if source == "load_write_ncpl" else 1), ml)
        rch.inirch_ncpl = source == "override"
    if source == "load_write_1":
        with pytest.warns(UserWarning, match=_LOW):
            text = _write(rch, ws, "t.rch")
    else:
        text = _write_quiet(rch, ws, "t.rch")
    expected = "hand_1" if source == "load_write_1" else "hand_ncpl"
    assert _inirch_of(text) == ["1" if expected == "hand_1" else "5"]
    assert _exe_run(ws) == runs[expected]


# ------------------------------------- structured NRCHOP=2 period 1, RTS
_MISSING_P1 = r"irch has no data for stress period 1: USG-T reads IRCH in period 1"
_S_P1_FLOW = "# RCH\n 2 50\n 1 -1\n" + _blk(_R6) + " -1 6\n" + _iblk([1] * 6)
_S_P1_CONC = (
    "# RCH\n 2 50 CONC\n 1\n 1 -1 INCONC\n"
    + _blk(_R6)
    + _blk(_C6)
    + " -1 6\n"
    + _iblk([1] * 6)
)


@pytest.mark.parametrize("case", ["api_flow", "api_conc", "hand_flow", "hand_conc"])
def test_rch_structured_irch_missing_period1(function_tmpdir, case):
    """Structured NRCHOP=2 without IRCH in period 1 (no irch data for period
    1, or a hand file with INIRCH -1 in period 1): USG-T would reuse an IRCH
    it never read and run without recharge, so write_file raises ValueError
    before opening the file, with or without CONC and inirch_ncpl; the hand
    file still loads."""
    conc = case.endswith("conc")
    ml = _dis(function_tmpdir, mcomp=1 if conc else 0)
    if case.startswith("api"):
        kwargs = {"iconc": 1, "rchconc": 1.0} if conc else {}
        rch = MfUsgRch(ml, nrchop=2, ipakcb=50, rech=1e-3, irch={1: 0}, **kwargs)
    else:
        rch = _load(function_tmpdir, _S_P1_CONC if conc else _S_P1_FLOW, ml)
        assert rch.inirch == {0: -1, 1: 6}
        assert list(rch.irch.transient_2ds) == [1]
    rch.fn_path = str(function_tmpdir / "out.rch")
    for override in (False, True):
        rch.inirch_ncpl = override
        with pytest.raises(ValueError, match=_MISSING_P1):
            rch.write_file(check=False)
        assert not (function_tmpdir / "out.rch").exists()


_S_RTS_LAYER2 = (
    "# RCH\n 2 50 RTS 1\n 0 6 INRCHZONES 1\n"
    + _blk([0.0] * 6)
    + _iblk([1, 2, 2, 1, 1, 1])
    + _iblk([1] * 6)
)


@pytest.mark.parametrize("case", ["api_layer2", "load_layer2", "api_layer1"])
def test_rch_structured_rts_irch_layer(function_tmpdir, case):
    """Structured RTS with NRCHOP=2: IZNRCH has NROW*NCOL positions and
    USG-T takes the zone of cell IRCH(NN), so a cell below layer 1 gets no
    RTS recharge (probe H3, both engines). An irch layer above 1 is not
    written (NotImplementedError before opening the file; the file still
    loads); irch in layer 1 is written."""
    ml = _dis(function_tmpdir, nper=1, mcomp=0, nlay=2)
    if case == "load_layer2":
        rch = _load(function_tmpdir, _S_RTS_LAYER2, ml)
    else:
        layers = [[0, 1, 1], [0, 0, 0]] if case == "api_layer2" else [[0] * 3] * 2
        rch = MfUsgRch(
            ml,
            nrchop=2,
            ipakcb=50,
            rech=0.0,
            irch=np.array(layers),
            mxrtzones=1,
            iznrch=np.ones((2, 3), dtype=int),
        )
    rch.fn_path = str(function_tmpdir / "out.rch")
    if case == "api_layer1":
        rch.write_file(check=False)
        again = _reload(rch.fn_path, ml)
        assert again.irch.transient_2ds[0].array.tolist() == [[0] * 3] * 2
    else:
        with pytest.raises(NotImplementedError, match=r"layer 2: with RTS.*IZNRCH"):
            rch.write_file(check=False)
        assert not (function_tmpdir / "out.rch").exists()


# --------------------------------------------- RTS with NRCHOP=2: INIRCH
# one line per row: USG-T reads a structured array row by row (U2DINT)
_ZN6 = "INTERNAL 1 (FREE) -1\n 1 2 1\n 2 1 2\n"
_C6_ROWS = "INTERNAL 1.0 (FREE) -1\n 5.0 6.0 7.0\n 8.0 9.0 10.0\n"
_ZN23 = np.array([[1, 2, 1], [2, 1, 2]])
_LOW_RTS = r"RTS records read during a period.*first INIRCH cells.*:3353.*inirch_ncpl"
_REUSE_RTS = r"with RTS, USG-T applies none of the RTS records.*every stress period"


def _s_rts(n1, n2, conc=False):
    """Structured 2 x 3, RTS 2 (and CONC with 1 component), INIRCH n1 in
    period 1 and n2 in period 2 (< 0: IRCH reused). RECH is read in both
    periods: load keeps a reused RECH as data, and with RTS a RECH read
    again applies the RTS record to every cell (gwf2rch8u1.f :429-441)."""
    out = (
        "# RCH\n 2 50"
        + (" CONC" if conc else "")
        + " RTS 2\n"
        + (" 1\n" if conc else "")
    )
    out += f" 0 {n1} INRCHZONES 2" + (" INCONC" if conc else "") + "\n"
    out += "CONSTANT 0.0\nCONSTANT 1\n" + _ZN6 + (_C6_ROWS if conc else "")
    out += f" 0 {n2}\nCONSTANT 0.0\n"
    return out + ("CONSTANT 1\n" if n2 >= 0 else "")


def _rts_api(ml, irch, conc=False, **kwargs):
    extra = {"iconc": 1, "rchconc": 1.0} if conc else {}
    return MfUsgRch(
        ml,
        nrchop=2,
        ipakcb=50,
        rech={0: 0.0, 1: 0.0},
        irch=irch,
        mxrtzones=2,
        iznrch={0: _ZN23},
        iniznrch={0: 2},
        **extra,
        **kwargs,
    )


@pytest.mark.parametrize("case", ["rts", "conc_rts", "none"])
def test_rch_rts_structured_inirch_api_default(function_tmpdir, case):
    """API without inirch: with RTS (alone or with CONC) the periods that
    read irch write NROW*NCOL without a warning; without RTS and CONC they
    write 1 as before."""
    ml = _dis(function_tmpdir, mcomp=1 if case == "conc_rts" else 0)
    if case == "none":
        rch = MfUsgRch(ml, nrchop=2, ipakcb=50, rech=0.0, irch={0: 0, 1: 0})
    else:
        rch = _rts_api(ml, {0: 0, 1: 0}, conc=case == "conc_rts")
    out = _write_quiet(rch, function_tmpdir)
    assert _inirch_of(out) == (["1", "1"] if case == "none" else ["6", "6"])


@pytest.mark.parametrize(
    "case", ["load_rts", "load_conc_rts", "api_rts", "api_conc_rts"]
)
def test_rch_rts_structured_inirch_low(function_tmpdir, case):
    """RTS with INIRCH below NROW*NCOL (loaded file or API inirch): written
    as given with one warning that names the options (RTS, or CONC and RTS)
    and the RTS flow effect; inirch_ncpl=True writes NROW*NCOL without a
    warning."""
    conc = "conc" in case
    ml = _dis(function_tmpdir, mcomp=1 if conc else 0)
    if case.startswith("load"):
        rch = _load(function_tmpdir, _s_rts(1, 3, conc), ml)
        assert rch.inirch == {0: 1, 1: 3}
    else:
        rch = _rts_api(ml, {0: 0, 1: 0}, conc, inirch={0: 1, 1: 3})
    with warnings.catch_warnings(record=True) as record:
        warnings.simplefilter("always")
        out = _write(rch, function_tmpdir)
    msgs = [str(w.message) for w in record if "INIRCH below" in str(w.message)]
    assert len(msgs) == 1
    head = "CONC and RTS with INIRCH" if conc else "MfUsgRch: RTS with INIRCH"
    assert head in msgs[0]
    assert "1 (INIRCH 1), 2 (INIRCH 3)" in msgs[0]
    assert re.search(_LOW_RTS, msgs[0])
    assert bool(re.search(_LOW, msgs[0])) is conc
    assert _inirch_of(out) == ["1", "3"]
    assert _reload(function_tmpdir / "out.rch", ml).inirch == {0: 1, 1: 3}
    rch.inirch_ncpl = True
    assert _inirch_of(_write_quiet(rch, function_tmpdir)) == ["6", "6"]


@pytest.mark.parametrize(
    "case", ["s_load", "s_api", "s_conc_load", "disu_load", "disu_api"]
)
def test_rch_rts_irch_reuse_not_written(function_tmpdir, case):
    """RTS with NRCHOP=2 and a period that reuses IRCH: USG-T applies none
    of the RTS records read after the start of that period (glo2basu1.f
    :3342, :3353), so FloPy raises NotImplementedError before opening the
    file (structured, also with inirch_ncpl, and DISU); the file still
    loads. With CONC the message gives both causes."""
    if case.startswith("disu"):
        ml = _ml(function_tmpdir, nper=2, mcomp=0)
        if case == "disu_load":
            text = (
                "# RCH\n 2 50 RTS 1\n 5\n 0 5 INRCHZONES 1\n"
                + _blk([0.0] * 5)
                + _iblk([1, 2, 3, 4, 5])
                + _iblk([1] * 5)
                + " -1 -1\n"
            )
            rch = _load(function_tmpdir, text, ml)
        else:
            rch = MfUsgRch(
                ml,
                nrchop=2,
                ipakcb=50,
                rech={0: np.zeros((1, 5))},
                irch={0: np.array([[0, 1, 2, 3, 4]])},
                mxrtzones=1,
                iznrch={0: np.ones((1, 5), dtype=int)},
                iniznrch={0: 1},
            )
    else:
        conc = "conc" in case
        ml = _dis(function_tmpdir, mcomp=1 if conc else 0)
        if case.endswith("load"):
            rch = _load(function_tmpdir, _s_rts(6, -1, conc), ml)
            assert rch.inirch == {0: 6, 1: -1}
        else:
            rch = _rts_api(ml, 0)
    rch.fn_path = str(function_tmpdir / "out.rch")
    for override in (False, True) if case.startswith("s_") else (False,):
        rch.inirch_ncpl = override
        with pytest.raises(NotImplementedError, match=_REUSE_RTS) as err:
            rch.write_file(check=False)
        assert ("with CONC" in str(err.value)) is ("conc" in case)
        assert "stress period 2 (INIRCH -1)" in str(err.value)
        assert not (function_tmpdir / "out.rch").exists()


_RTS4 = (
    "0.0 4.0 1.0 1.0e-3 2.0e-3\n4.0 10.0 1.0 3.0e-3 1.0e-3\n"
    "10.0 14.0 1.0 2.0e-3 2.0e-3\n14.0 1.0e6 1.0 5.0e-3 4.0e-3\n"
)


def _rts_exe_model(ws):
    """Structured 2 x 3, one layer, two periods with ATSA, CHD in cell
    (1, 1), flow only; RCH (t.rch) and RTS (rech.rts, 4 records, two of
    them read inside the periods) added to the name file."""
    ml = _dis(ws, nper=2, mcomp=0)
    MfUsgBas(ml, ibound=1, strt=10.0)
    MfUsgLpf(ml, laytyp=0, hk=1.0, vka=1.0, ss=1.0e-3, ipakcb=50)
    MfUsgSms(ml, linmeth=1)
    ModflowChd(ml, stress_period_data={0: [[0, 0, 0, 10.0, 10.0]]})
    MfUsgOc(ml, stress_period_data={(0, 0): ["save head"]})
    ml.write_input()
    ws = Path(ws)
    nam = (ws / "t.nam").read_text().splitlines()
    unit = next(ln.split()[1] for ln in nam if ln.split()[-1] == "t.hds")
    (ws / "t.oc").write_text(
        f"# OC\nATSA\nHEAD PRINT FORMAT   0\nHEAD SAVE UNIT    {unit}\n"
        "DRAWDOWN PRINT FORMAT   0\nCOMPACT BUDGET AUX\n"
        + "".join(
            f"period {k} step 1\n  save head\n  save budget\n  print budget\n"
            for k in (1, 2)
        )
    )
    (ws / "rech.rts").write_text(_RTS4)
    with open(ws / "t.nam", "a") as f:
        f.write("RCH               18  t.rch\nRTS               40  rech.rts\n")
    return ml


@requires_usgt
@pytest.mark.parametrize(
    "source", ["api", "load_write_1", "load_write_ncpl", "override"]
)
def test_usgt_exe_rch_rts_structured_inirch(function_tmpdir, source):
    """Structured RTS with ATSA and RTS records read inside the periods:
    the hand files with INIRCH 1 and NROW*NCOL give different flow (the RTS
    records read during a period update the first INIRCH cells only).
    FloPy's API file runs like INIRCH NROW*NCOL; a loaded file keeps its
    result; inirch_ncpl=True on the INIRCH 1 file runs like NROW*NCOL."""
    runs = {}
    for name, n in (("hand_1", 1), ("hand_ncpl", 6)):
        ws = function_tmpdir / name
        _rts_exe_model(ws)
        (ws / "t.rch").write_text(_s_rts(n, n))
        runs[name] = _exe_run(ws)
    assert runs["hand_1"]["rch"] != runs["hand_ncpl"]["rch"]
    assert runs["hand_1"]["hds"] != runs["hand_ncpl"]["hds"]

    ws = function_tmpdir / source
    ml = _rts_exe_model(ws)
    if source == "api":
        rch = _rts_api(ml, {0: 0, 1: 0})
    else:
        rch = _load(
            ws, _s_rts(6, 6) if source == "load_write_ncpl" else _s_rts(1, 1), ml
        )
        rch.inirch_ncpl = source == "override"
    if source == "load_write_1":
        with pytest.warns(UserWarning, match=_LOW_RTS):
            text = _write(rch, ws, "t.rch")
    else:
        text = _write_quiet(rch, ws, "t.rch")
    expected = "hand_1" if source == "load_write_1" else "hand_ncpl"
    assert _inirch_of(text) == (["1", "1"] if expected == "hand_1" else ["6", "6"])
    assert _exe_run(ws) == runs[expected]
