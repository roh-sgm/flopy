"""MfUsgRch RTS zones (INRCHZONES, IZNRCH) and seepage elevations (INSELEV,
SELEV) per stress period.

USG-T RP (gwf2rch8u1.f, 2.8.1 and 2.7.0): INRCHZONES n > 0 reads IZNRCH and
sets the number of RTS zones to n; the RTS array is allocated in period 1, so
period 1 reads zones and the count stays fixed; n <= 0 or no keyword reuses.
INSELEV n > 0 reads SELEV and turns seepage on, 0 or no keyword turns it off,
n < 0 keeps the previous state (a -1 after a 0 does not turn it on again).
Arrays are compared by value, not by their text.
"""

import io
import re
from pathlib import Path

import numpy as np
import pytest

from autotest.test_usg_transport import (
    _G01_RECH,
    _G01_SELEV,
    _G05_ACT,
    _G05_HDR,
    _G05_PERIODS,
    _g01_read,
    _g03_fixture,
    _g05_dis,
    _g05_model,
)
from autotest.usgt_exe import USGT_EXE, requires_usgt
from flopy.mbase import run_model
from flopy.mfusg import MfUsgRch
from flopy.utils import HeadFile, HeadUFile, Transient2d

# ------------------------------------------------------------- fixtures
# SEEPELEV: the 5-cell model of _g05_model with 4 periods (SS, TR x3)
_P4 = _G05_PERIODS + "    10.0   5   1.0   TR\n" * 2
_OC4 = (
    "# OC\n"
    + _G05_HDR
    + "".join(
        f"period {k} step {s}\n" + _G05_ACT
        for k, n in ((1, 1), (2, 5), (3, 5), (4, 5))
        for s in range(1, n + 1)
    )
)
_S1, _S2 = _G01_SELEV


def _blk(values):
    return "INTERNAL 1.0 (FREE) -1\n " + " ".join(values) + "\n"


def _seep_text(flags, arrays):
    """flags: INSELEV token of each period (None: no keyword); arrays: SELEV
    values of each period or None. RECH: periods 1 and 2, then reused."""
    out = "# RCH\n 3 50 SEEPELEV\n"
    for k, (flag, arr) in enumerate(zip(flags, arrays)):
        out += (" 0" if k < 2 else " -1") + (f" INSELEV {flag}" if flag else "")
        out += "\n" + (_blk(_G01_RECH[k]) if k < 2 else "")
        out += _blk(arr) if arr is not None else ""
    return out


_Q1 = (["1", "0", "-1", "1"], [_S1, None, None, _S2])
_Q4 = (["1", "-1", "-1", "1"], [_S1, None, None, _S2])
_Z1 = (["0", "1", "-1", "-1"], [None, _S1, None, None])


def _seep_model(ws, text, grid="dis"):
    ws = _g05_model(ws, _OC4, grid, _P4)
    (ws / "t.rch").write_text(text)
    with open(ws / "t.nam", "a") as f:
        f.write("RCH 18 t.rch\n")
    return ws


def _arr(values):
    return np.array([float(v) for v in values]).reshape(1, 5)


_A10, _A20 = np.full((1, 5), 10.0), np.full((1, 5), 20.0)
_Z0 = np.zeros((1, 5), int)


# RTS: the 3-node G03 model with 3 periods (ATS) and 2 zones
_RTS2 = (
    "0.0 4.0 1.0 1.0e-3 2.0e-3\n4.0 12.0 1.0 3.0e-3 1.0e-3\n"
    "12.0 22.0 1.0 0.5e-3 4.0e-3\n22.0 1.0e6 1.0 2.0e-3 0.0\n"
)
_ZONES = ("0 1 2", "0 2 1", None)
_T8_SEEP = (("1", "10.2 10.05 10.01"), ("-1", None), ("0", None))


def _rts_text(sp=("2", "2", None), zones=_ZONES, seep=None):
    """sp: INRCHZONES token of each period (None: no keyword)."""
    out = "# RCH\n 1 0 RTS 2" + (" SEEPELEV" if seep else "") + "\n"
    for k in range(3):
        out += (" 0" if k == 0 else " -1") + (f" INRCHZONES {sp[k]}" if sp[k] else "")
        out += f" INSELEV {seep[k][0]}\n" if seep else "\n"
        out += "CONSTANT 0.0\n" if k == 0 else ""
        if sp[k] and int(sp[k]) > 0:
            out += f"INTERNAL 1 (FREE) -1\n {zones[k]}\n"
        if seep and seep[k][1]:
            out += f"INTERNAL 1.0 (FREE) -1\n {seep[k][1]}\n"
    return out


def _rts_model(ws, text):
    ws = _g03_fixture(ws, evs=None)
    disu = (ws / "t.disu").read_text().replace("3 1 7 0 1 4 2 0", "3 1 7 0 3 4 2 0")
    period = "     10.000000             1  1.000000  TR\n"
    (ws / "t.disu").write_text(disu.replace(period, period * 3))
    (ws / "t.chd").write_text((ws / "t.chd").read_text() + " -1\n -1\n")
    (ws / "t.oc").write_text(
        (ws / "t.oc").read_text()
        + "".join(f"period {k} step 1\n  save head\n  print budget\n" for k in (2, 3))
    )
    (ws / "t.rch").write_text(text)
    (ws / "ts" / "rech.rts").write_text(_RTS2)
    return ws


def _heads(path_or_text):
    text = path_or_text
    if isinstance(path_or_text, Path):
        text = path_or_text.read_text()
    return [
        ln.split("#")[0].split() for ln in text.splitlines() if "Stress period" in ln
    ]


def _rewrite(ws):
    rch = MfUsgRch.load(str(ws / "t.rch"), _g05_dis(ws), check=False)
    rch.fn_path = str(ws / "t.rch")
    rch.write_file(check=False)
    return rch


def _seep_rch(ws, **kwargs):
    return MfUsgRch(
        _g05_dis(ws),
        nrchop=3,
        ipakcb=50,
        rech={0: _arr(_G01_RECH[0]), 1: _arr(_G01_RECH[1])},
        seepelev=1,
        **kwargs,
    )


# ------------------------------------------------------- load and write


def test_rch_inselev_1_0_m1_1_load_write(function_tmpdir):
    """INSELEV 1 / 0 / -1 / 1: the flags are kept, SELEV only in periods 1
    and 4, and the period 4 array is the one of the file."""
    ws = _seep_model(function_tmpdir, _seep_text(*_Q1))
    rch = _rewrite(ws)
    assert sorted(rch.selev.transient_2ds) == [0, 3]
    np.testing.assert_array_equal(rch.selev[2].array, _arr(_S1).astype(np.float32))
    opts, periods = _g01_read(ws / "t.rch", 5)
    assert [p[0][1:] for p in periods] == [["INSELEV", f] for f in _Q1[0]]
    assert [p[1].get("selev") for p in periods] == [
        [float(v) for v in _S1],
        None,
        None,
        [float(v) for v in _S2],
    ]
    # the loader keeps RECH of every period (written again, same values)
    r1, r2 = ([float(v) for v in r] for r in _G01_RECH)
    assert [p[1].get("rech") for p in periods] in (
        [r1, r2, None, None],
        [r1, r2, r2, r2],
    )


@pytest.mark.parametrize("first", ["0", None, "-1"], ids=["zero", "absent", "m1"])
def test_rch_inselev_off_first_then_on(function_tmpdir, first):
    """Seepage off in period 1 (0, no keyword or -1) and SELEV read in
    period 2: period 1 is written as read (no keyword as 0)."""
    flags = [first, "1", "-1", "-1"]
    ws = _seep_model(function_tmpdir, _seep_text(flags, _Z1[1]))
    rch = _rewrite(ws)
    assert sorted(rch.selev.transient_2ds) == [1]
    _, periods = _g01_read(ws / "t.rch", 5)
    assert [p[0][1:] for p in periods] == [
        ["INSELEV", first or "0"],
        ["INSELEV", "1"],
        ["INSELEV", "-1"],
        ["INSELEV", "-1"],
    ]
    assert periods[1][1]["selev"] == [float(v) for v in _S1]


def test_rch_inselev_api_defaults(function_tmpdir):
    """Without inselev: 1 in the periods with SELEV, -1 in the others (also
    period 1 without SELEV: seepage stays off)."""
    ws = _seep_model(function_tmpdir, _seep_text(*_Z1))
    rch = _seep_rch(ws, selev={1: _arr(_S1), 3: _arr(_S2)})
    rch.fn_path = str(ws / "t.rch")
    rch.write_file(check=False)
    assert [h[1:] for h in _heads(ws / "t.rch")] == [
        ["INSELEV", "-1"],
        ["INSELEV", "1"],
        ["INSELEV", "-1"],
        ["INSELEV", "1"],
    ]


def test_rch_inselev_list_is_per_period(function_tmpdir):
    """inselev as one int per period writes like the same dict."""
    ws = _seep_model(function_tmpdir, _seep_text(*_Q1))
    texts = []
    for flags in ([1, 0, -1, 1], {0: 1, 1: 0, 2: -1, 3: 1}):
        rch = _seep_rch(ws, selev={0: _arr(_S1), 3: _arr(_S2)}, inselev=flags)
        rch.fn_path = str(ws / "t.rch")
        rch.write_file(check=False)
        texts.append((ws / "t.rch").read_text())
    assert texts[0] == texts[1]
    assert [h[1:] for h in _heads(ws / "t.rch")] == [["INSELEV", f] for f in _Q1[0]]


@pytest.mark.parametrize("flag", ["-1", "0"])
def test_rch_inrchzones_reuse_load_write(function_tmpdir, flag):
    """INRCHZONES -1 or 0 in period 3 is kept and reads no array; no
    keyword stays without keyword; period 2 zones are its own."""
    ws = _rts_model(function_tmpdir, _rts_text(sp=("2", "2", flag)))
    rch = _rewrite(ws)
    assert sorted(rch.iznrch.transient_2ds) == [0, 1]
    np.testing.assert_array_equal(rch.iznrch[2].array, [[0, 2, 1]])
    _, periods = _g01_read(ws / "t.rch", 3)
    assert [p[0][1:] for p in periods] == [
        ["INRCHZONES", "2"],
        ["INRCHZONES", "2"],
        ["INRCHZONES", flag],
    ]
    assert [p[1].get("iznrch") for p in periods] == [[0, 1, 2], [0, 2, 1], None]
    ws2 = _rts_model(function_tmpdir / "absent", _rts_text())
    _rewrite(ws2)
    assert [h[1:] for h in _heads(ws2 / "t.rch")] == [
        ["INRCHZONES", "2"],
        ["INRCHZONES", "2"],
        [],
    ]


def test_rch_iznrch_given_once(function_tmpdir):
    """iznrch only for period 1: INRCHZONES mxrtzones and the zones in
    period 1, no keyword and no array later."""
    ws = _rts_model(function_tmpdir, _rts_text())
    rch = MfUsgRch(
        _g05_dis(ws), nrchop=1, rech=0.0, mxrtzones=2, iznrch=np.array([[0, 1, 2]])
    )
    rch.fn_path = str(ws / "t.rch")
    rch.write_file(check=False)
    _, periods = _g01_read(ws / "t.rch", 3)
    assert [p[0][1:] for p in periods] == [["INRCHZONES", "2"], [], []]
    assert [p[1].get("iznrch") for p in periods] == [[0, 1, 2], None, None]


def test_rch_rts_and_seepelev_order(function_tmpdir):
    """RTS and SEEPELEV together: RECH, IZNRCH, SELEV in RP order and the
    record after each period read in place."""
    ws = _rts_model(function_tmpdir, _rts_text(seep=_T8_SEEP))
    rch = _rewrite(ws)
    assert sorted(rch.iznrch.transient_2ds) == [0, 1]
    assert sorted(rch.selev.transient_2ds) == [0]
    _, periods = _g01_read(ws / "t.rch", 3)
    assert [p[0][1:] for p in periods] == [
        ["INRCHZONES", "2", "INSELEV", "1"],
        ["INRCHZONES", "2", "INSELEV", "-1"],
        ["INSELEV", "0"],
    ]
    assert periods[0][1] == {
        "rech": [0.0] * 3,
        "iznrch": [0, 1, 2],
        "selev": [10.2, 10.05, 10.01],
    }
    assert periods[1][1].get("iznrch") == [0, 2, 1] and "selev" not in periods[1][1]
    assert periods[1][1].get("rech", [0.0] * 3) == [0.0] * 3


# ------------------------------------------------------------- checks


def _rejects(match, ws, **kwargs):
    ml = _g05_dis(ws)
    before = (list(ml.packagelist), list(ml.output_units), list(ml.output_fnames))
    with pytest.raises(ValueError, match=match):
        MfUsgRch(ml, **{"nrchop": 3, "ipakcb": 50, "rech": 1e-3, **kwargs})
    assert (list(ml.packagelist), list(ml.output_units), list(ml.output_fnames)) == (
        before
    )


@pytest.mark.parametrize(
    "kwargs, match",
    [
        ({"iznrch": np.zeros((1, 5), int)}, "RTS"),
        ({"iniznrch": {0: 1}}, "RTS"),
        ({"selev": np.ones((1, 5))}, "SEEPELEV"),
        ({"inselev": {0: 1}}, "SEEPELEV"),
        ({"mxrtzones": 2}, "period 1"),
        ({"mxrtzones": 2, "iznrch": {1: np.zeros((1, 5), int)}}, "period 1"),
        (
            {"mxrtzones": 2, "iznrch": np.zeros((1, 5), int), "iniznrch": {0: 0}},
            "does not read",
        ),
        (
            {"mxrtzones": 2, "iznrch": np.zeros((1, 5), int), "iniznrch": {1: 2}},
            "needs iznrch",
        ),
        (
            {
                "mxrtzones": 2,
                "iznrch": {0: np.zeros((1, 5), int), 1: np.ones((1, 5), int)},
                "iniznrch": {0: 3, 1: 2},
            },
            "fixed",
        ),
        (
            {
                "mxrtzones": 2,
                "iznrch": {0: np.zeros((1, 5), int), 1: np.ones((1, 5), int)},
                "iniznrch": {1: 1},
            },
            "fixed",
        ),
        (
            {"seepelev": 1, "selev": {1: np.ones((1, 5))}, "inselev": {1: -1}},
            "does not read",
        ),
        (
            {"seepelev": 1, "selev": {0: np.ones((1, 5))}, "inselev": {2: 1}},
            "needs selev",
        ),
        ({"seepelev": 1, "inselev": {0: 1}}, "needs selev"),
        ({"seepelev": 1, "selev": np.ones((1, 4))}, "shape"),
        ({"mxrtzones": 2, "iznrch": np.ones((1, 6), int)}, "shape"),
        ({"seepelev": 1, "selev": {5: np.ones((1, 5))}}, "stress period"),
        ({"seepelev": 1, "selev": np.ones((1, 5)), "inselev": {0: 1.0}}, "integer"),
        ({"seepelev": 1, "selev": np.ones((1, 5)), "inselev": {0: True}}, "integer"),
        (
            {"mxrtzones": 2, "iznrch": np.ones((1, 5), int), "iniznrch": {0: 2.5}},
            "integer",
        ),
        ({"mxrtzones": 2.0, "iznrch": np.ones((1, 5), int)}, "mxrtzones"),
        ({"seepelev": 2}, "seepelev"),
        ({"seepelev": 1, "selev": np.ones((1, 5)), "inselev": [1, -1]}, "values"),
        (
            {"seepelev": 1, "selev": np.ones((1, 5)), "inselev": {-1: 1}},
            "stress period",
        ),
        ({"seepelev": 1, "selev": {0: _A10, 0.9: _A20}}, "key 0.9"),
        ({"seepelev": 1, "selev": {0: _A10, -0.5: _A20}}, "key -0.5"),
        ({"seepelev": 1, "selev": {True: _A10}}, "key True"),
        ({"seepelev": 1, "selev": {"0": _A10}}, "key '0'"),
        ({"seepelev": 1, "selev": {0: _A10, 4: _A20}}, "key 4"),
        ({"mxrtzones": 2, "iznrch": {0: _Z0, 0.5: _Z0}}, "key 0.5"),
        ({"mxrtzones": 2, "iznrch": {np.float64(0): _Z0}}, "key"),
    ],
    ids=[
        "iznrch_no_rts",
        "iniznrch_no_rts",
        "selev_no_seepelev",
        "inselev_no_seepelev",
        "rts_no_zones",
        "rts_zones_not_period1",
        "zones_flag0",
        "flag_no_zones",
        "count_later_differs_from_period1",
        "count_change",
        "selev_flag_m1",
        "inselev_no_selev",
        "inselev_selev_none",
        "selev_shape",
        "iznrch_shape",
        "selev_period_range",
        "inselev_float",
        "inselev_bool",
        "iniznrch_float",
        "mxrtzones_float",
        "seepelev_2",
        "inselev_list_length",
        "inselev_negative_period",
        "selev_key_collision_0_09",
        "selev_key_negative_fraction",
        "selev_key_bool",
        "selev_key_text",
        "selev_key_nper",
        "iznrch_key_fraction",
        "iznrch_key_numpy_float",
    ],
)
def test_rch_rts_seepelev_rejected(function_tmpdir, kwargs, match):
    """Rejected before the package is added to the model (no package, no
    cbc output unit)."""
    ws = _seep_model(function_tmpdir, _seep_text(*_Q1), "disu")
    _rejects(match, ws, **kwargs)


def test_rch_rts_period1_count_differs_from_header(function_tmpdir):
    """The period 1 INRCHZONES sets the number of RTS zones (the header
    count is not used): kept on load and write, used as the default of later
    periods with zones, also after a change of iniznrch."""
    ws = _rts_model(function_tmpdir, _rts_text().replace("RTS 2", "RTS 1"))
    rch = _rewrite(ws)
    assert (rch.mxrtzones, rch.iniznrch) == (1, {0: 2, 1: 2})
    opts, periods = _g01_read(ws / "t.rch", 3)
    assert opts == ["1", "0", "RTS", "1"]
    assert [p[0][1:] for p in periods] == [["INRCHZONES", "2"]] * 2 + [[]]
    zones = {0: np.array([[0, 1, 2]]), 1: np.array([[0, 2, 1]])}
    rch = MfUsgRch(
        _g05_dis(ws), nrchop=1, rech=0.0, mxrtzones=1, iznrch=zones, iniznrch={0: 2}
    )
    rch.fn_path = str(ws / "t.rch")
    rch.write_file(check=False)
    assert [h[1:] for h in _heads(ws / "t.rch")] == [["INRCHZONES", "2"]] * 2 + [[]]
    rch = MfUsgRch(_g05_dis(ws), nrchop=1, rech=0.0, mxrtzones=2, iznrch=zones)
    rch.iniznrch = {0: 3}
    rch.write_file(check=False)
    assert [h[1:] for h in _heads(ws / "t.rch")] == [["INRCHZONES", "3"]] * 2 + [[]]
    rch.iniznrch = {0: 3, 1: 2}
    with pytest.raises(ValueError, match="fixed in period 1 to 3"):
        rch.write_file(check=False)


def test_rch_numpy_integer_keys(function_tmpdir):
    """numpy integer keys are stress periods."""
    ws = _seep_model(function_tmpdir, _seep_text(*_Q1))
    rch = _seep_rch(ws, selev={np.int64(0): _arr(_S1), np.int32(3): _arr(_S2)})
    assert sorted(rch.selev.transient_2ds) == [0, 3]
    rch.fn_path = str(ws / "t.rch")
    rch.write_file(check=False)
    assert [h[1:] for h in _heads(ws / "t.rch")] == [
        ["INSELEV", f] for f in ("1", "-1", "-1", "1")
    ]


@pytest.mark.parametrize(
    "key", [0.5, "1", True, 4], ids=["fraction", "text", "bool", "nper"]
)
def test_rch_array_keys_checked_at_write(function_tmpdir, key):
    """A key put into selev.transient_2ds after construction is checked
    before any output."""
    ws = _seep_model(function_tmpdir, _seep_text(*_Q1))
    rch = _seep_rch(ws, selev={0: _arr(_S1), 3: _arr(_S2)})
    rch.fn_path = str(ws / "t.rch")
    old = (ws / "t.rch").read_text()
    rch.selev.transient_2ds[key] = rch.selev.transient_2ds[3]
    with pytest.raises(ValueError, match="stress period key"):
        rch.write_file(check=False)
    assert (ws / "t.rch").read_text() == old


def _t2d_state(t2d):
    return (
        id(t2d.transient_2ds),
        {
            k: (id(u), u.locat, u.format.fortran, u.cnstnt, u.filename, u.array.copy())
            for k, u in t2d.transient_2ds.items()
        },
    )


def _same_state(a, b):
    assert a[0] == b[0] and a[1].keys() == b[1].keys()
    for k in a[1]:
        assert a[1][k][:5] == b[1][k][:5]
        np.testing.assert_array_equal(a[1][k][5], b[1][k][5])


@pytest.mark.parametrize("stage", ["early", "late", "accepted"])
@pytest.mark.parametrize("name", ["selev", "iznrch"])
def test_rch_given_arrays_unchanged(function_tmpdir, name, stage):
    """A Transient2d (locat, format, factor) or a dict given as selev /
    iznrch is not changed by a rejected (early: option missing; late:
    flag contradiction) or an accepted construction; the model is not
    changed by a rejection."""
    ws = _seep_model(function_tmpdir, _seep_text(*_Q1), "disu")
    ml = _g05_dis(ws)
    dtype = np.float32 if name == "selev" else np.int32
    t2d = Transient2d(
        ml,
        (1, 5),
        dtype,
        {0: 2.0, 2: 3.0},
        name="given",
        fmtin="(5G15.6)",
        cnstnt=2.0,
        locat=77,
    )
    raw = {0: np.full((1, 5), 2, dtype), 2: np.full((1, 5), 3, dtype)}
    raw_ids = {k: id(v) for k, v in raw.items()}
    option = {"selev": {"seepelev": 1}, "iznrch": {"mxrtzones": 2}}[name]
    flags = {"selev": {"inselev": {0: 0}}, "iznrch": {"iniznrch": {0: 0}}}[name]
    kwargs = {"early": {}, "late": {**option, **flags}, "accepted": option}[stage]
    before = (list(ml.packagelist), list(ml.output_units), list(ml.output_fnames))
    for value in (t2d, raw):
        state = _t2d_state(t2d)
        if stage == "accepted":
            MfUsgRch(ml, nrchop=3, rech=1e-3, **{name: value}, **kwargs)
        else:
            with pytest.raises(ValueError):
                MfUsgRch(ml, nrchop=3, ipakcb=50, rech=1e-3, **{name: value}, **kwargs)
            assert (
                list(ml.packagelist),
                list(ml.output_units),
                list(ml.output_fnames),
            ) == before
        _same_state(state, _t2d_state(t2d))
        assert {k: id(v) for k, v in raw.items()} == raw_ids
        assert [float(raw[k].ravel()[0]) for k in (0, 2)] == [2.0, 3.0]


def test_rch_checked_again_after_changes(function_tmpdir):
    """write_file checks the changed attributes before any output: the
    existing file, the .chk file and a caller handle are not touched."""
    ws = _seep_model(function_tmpdir, _seep_text(*_Q1))
    rch = _seep_rch(ws, selev={0: _arr(_S1), 3: _arr(_S2)})
    rch.fn_path = str(ws / "t.rch")
    old = (ws / "t.rch").read_text()
    changes = [
        ("inselev", {1: 1}, "needs selev"),
        ("inselev", {0: 0}, "does not read"),
        ("inselev", {0: 1.5}, "integer"),
        ("seepelev", 0, "SEEPELEV"),
        ("iniznrch", {0: 1}, "RTS"),
    ]
    for name, value, match in changes:
        saved = getattr(rch, name)
        setattr(rch, name, value)
        handle = io.StringIO()
        with pytest.raises(ValueError, match=match):
            rch.write_file(check=False, f=handle)
        assert handle.getvalue() == "" and not handle.closed
        with pytest.raises(ValueError, match=match):
            rch.write_file(check=True)
        assert (ws / "t.rch").read_text() == old
        assert not list(Path(ws).glob("*.chk"))
        setattr(rch, name, saved)
    rch.selev[1] = _arr(_S2)  # a new array for period 2
    rch.inselev = {1: -1}
    with pytest.raises(ValueError, match="does not read"):
        rch.write_file(check=False)
    rch.inselev = None
    handle = io.StringIO()
    rch.write_file(check=False, f=handle)
    assert not handle.closed
    assert [h[1:] for h in _heads(handle.getvalue())] == [
        ["INSELEV", "1"],
        ["INSELEV", "1"],
        ["INSELEV", "-1"],
        ["INSELEV", "1"],
    ]


# ----------------------------------------------------------------- engine


def _run(ws):
    ok, buff = run_model(USGT_EXE, "t.nam", model_ws=str(ws), silent=True, report=True)
    lst = (Path(ws) / "t.list").read_text()
    out = {
        "normal": ok and any("normal termination" in b.lower() for b in buff),
        "rch": re.findall(r"^\s+RECHARGE =\s+(\S+)\s+RECHARGE =\s+(\S+)", lst, re.M),
        "pct": re.findall(r"PERCENT DISCREPANCY =\s+(\S+)", lst),
    }
    assert out["normal"] and "NEGATIVE TIME STEP" not in lst
    cls = HeadFile if (Path(ws) / "t.dis").exists() else HeadUFile
    with cls(str(Path(ws) / "t.hds")) as hf:
        out["times"] = hf.get_times()
        out["heads"] = [
            np.asarray(hf.get_data(idx=i), float).ravel()
            for i in range(len(out["times"]))
        ]
    assert out["times"] and out["rch"]
    assert all(np.all(np.isfinite(h)) for h in out["heads"])
    assert np.all(np.isfinite([float(v) for r in out["rch"] for v in r]))
    return out


def _same(a, b):
    assert (a["rch"], a["pct"], a["times"]) == (b["rch"], b["pct"], b["times"])
    for x, y in zip(a["heads"], b["heads"]):
        np.testing.assert_array_equal(x, y)


def _differ(a, b):
    assert a["rch"] != b["rch"] or any(
        not np.array_equal(x, y) for x, y in zip(a["heads"], b["heads"])
    )


@requires_usgt
@pytest.mark.parametrize("source", ["load_write", "api"])
def test_usgt_exe_rch_rts_two_zones_periods(function_tmpdir, source):
    """Two RTS zones in three periods: zones read in period 1, other zones
    in period 2 (same count), reused in period 3. FloPy's file runs like
    the control; without the period 2 zones the heads change."""
    control = _run(_rts_model(function_tmpdir / "control", _rts_text()))
    other = _run(_rts_model(function_tmpdir / "other", _rts_text(sp=("2", None, None))))
    _differ(other, control)
    ws = _rts_model(function_tmpdir / source, _rts_text())
    if source == "load_write":
        _rewrite(ws)
    else:
        rch = MfUsgRch(
            _g05_dis(ws),
            nrchop=1,
            rech=0.0,
            mxrtzones=2,
            iznrch={0: np.array([[0, 1, 2]]), 1: np.array([[0, 2, 1]])},
        )
        rch.fn_path = str(ws / "t.rch")
        rch.write_file(check=False)
    _same(_run(ws), control)


@requires_usgt
@pytest.mark.parametrize("source", ["load_write", "api"])
def test_usgt_exe_rch_rts_header_count_differs(function_tmpdir, source):
    """Header RTS 1 with INRCHZONES 2 runs like RTS 2 (USG-T uses the
    period 1 count); FloPy keeps it and its file runs like the control."""
    text = _rts_text(zones=("0 2 1", "0 1 2", None), sp=("2", "2", "-1"))
    control = _run(
        _rts_model(function_tmpdir / "control", text.replace("RTS 2", "RTS 1"))
    )
    _same(_run(_rts_model(function_tmpdir / "header2", text)), control)
    ws = _rts_model(function_tmpdir / source, text.replace("RTS 2", "RTS 1"))
    if source == "load_write":
        _rewrite(ws)
    else:
        rch = MfUsgRch(
            _g05_dis(ws),
            nrchop=1,
            rech=0.0,
            mxrtzones=1,
            iznrch={0: np.array([[0, 2, 1]]), 1: np.array([[0, 1, 2]])},
            iniznrch={0: 2},
        )
        rch.fn_path = str(ws / "t.rch")
        rch.write_file(check=False)
    assert (ws / "t.rch").read_text().splitlines()[1].split() == ["1", "0", "RTS", "1"]
    _same(_run(ws), control)


@requires_usgt
@pytest.mark.parametrize("flag", ["-1", "0"])
def test_usgt_exe_rch_inrchzones_reuse(function_tmpdir, flag):
    """INRCHZONES -1 / 0 in period 3 reuses the zones (= no keyword); the
    rewritten file runs like it."""
    control = _run(_rts_model(function_tmpdir / "control", _rts_text()))
    ws = _rts_model(function_tmpdir / "w", _rts_text(sp=("2", "2", flag)))
    _same(_run(ws), control)
    _rewrite(ws)
    _same(_run(ws), control)


@requires_usgt
@pytest.mark.parametrize("source", ["load_write", "api"])
def test_usgt_exe_rch_inselev_1_0_m1_1(function_tmpdir, source):
    """INSELEV 1 / 0 / -1 / 1 (the -1 after a 0 keeps seepage off): FloPy's
    file runs like the control; 1 / -1 / -1 / 1 gives other results."""
    control = _run(_seep_model(function_tmpdir / "control", _seep_text(*_Q1)))
    _differ(_run(_seep_model(function_tmpdir / "q4", _seep_text(*_Q4))), control)
    ws = _seep_model(function_tmpdir / source, _seep_text(*_Q1))
    if source == "load_write":
        _rewrite(ws)
    else:
        rch = _seep_rch(ws, selev={0: _arr(_S1), 3: _arr(_S2)}, inselev=[1, 0, -1, 1])
        rch.fn_path = str(ws / "t.rch")
        rch.write_file(check=False)
    _same(_run(ws), control)


@requires_usgt
@pytest.mark.parametrize("source", ["load_write", "api"])
def test_usgt_exe_rch_seepelev_on_in_period2(function_tmpdir, source):
    """Seepage off in period 1 and on from period 2 (INSELEV 0 / 1 / -1 /
    -1): FloPy's file (api: SELEV only for period 2, default flags) runs
    like the control; never on gives other results."""
    control = _run(_seep_model(function_tmpdir / "control", _seep_text(*_Z1)))
    never = _seep_text(["0"] * 4, [None] * 4)
    _differ(_run(_seep_model(function_tmpdir / "never", never)), control)
    ws = _seep_model(function_tmpdir / source, _seep_text(*_Z1))
    if source == "load_write":
        _rewrite(ws)
    else:
        rch = _seep_rch(ws, selev={1: _arr(_S1)})
        rch.fn_path = str(ws / "t.rch")
        rch.write_file(check=False)
    _same(_run(ws), control)


@requires_usgt
def test_usgt_exe_rch_rts_seepelev_together(function_tmpdir):
    """RTS zones and SEEPELEV in the same file (INSELEV 1 / -1 / 0): the
    rewritten file runs like the control; SEEPELEV changes the results."""
    control = _run(_rts_model(function_tmpdir / "control", _rts_text(seep=_T8_SEEP)))
    _differ(_run(_rts_model(function_tmpdir / "plain", _rts_text())), control)
    ws = _rts_model(function_tmpdir / "w", _rts_text(seep=_T8_SEEP))
    _rewrite(ws)
    _same(_run(ws), control)
