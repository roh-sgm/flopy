"""MfUsgOc, output control by words: solver words per stress period, REAL
text of HCLOSE and BTOL, COMPACT BUDGET with and without AUX, and records
USG-T does not accept.

Checked against USG-T 2.8.1 glo2basu1.f: SGWF2BAS7N (4826-4898 solver words,
4916-4917 unknown word) and SGWF2BAS7J (1450-1467 COMPACT BUDGET [AUX],
1497-1506 skipped and unknown words). USG-T 2.7.0 reads HCLOSE, BTOL and
MXITER only (4826-4846) and stops at the other solver words (4864-4865).

The OC files are hand-written. The engine runs need USGT_EXE.
"""

import inspect
import re
from pathlib import Path

import numpy as np
import pytest

from autotest.test_usg_transport import _g04_sentinel, _g05_dis, _g05_model
from autotest.test_usg_transport_exe import _g05x_run
from autotest.usgt_exe import requires_usgt
from flopy.mfusg import MfUsgOc
from flopy.utils import CellBudgetFile

_HDR = "# OC\nHEAD PRINT FORMAT 0\nHEAD SAVE UNIT 51\nDRAWDOWN PRINT FORMAT 0\n"
_ACT = "  save head\n  save budget\n  print budget\n"

# (1, 4): stress period 2, time step 5 of the model of _g05_model
_P2 = (1, 4)


def _text(p2="", compact="COMPACT BUDGET AUX\n", header=""):
    """Hand-written OC by words; p2 goes first in the last time step."""
    return (
        _HDR
        + compact
        + header
        + "period 1 step 1\n"
        + _ACT
        + "period 2 step 5\n"
        + p2
        + _ACT
    )


def _load(ws, text):
    ws = _g05_model(ws, text)
    return MfUsgOc.load(str(ws / "t.oc"), _g05_dis(ws))


def _rewrite(oc, name="out.oc"):
    """Write oc next to its model; return (lines without heading, reloaded)."""
    path = Path(oc.parent.model_ws) / name
    oc.fn_path = str(path)
    oc.write_file()
    again = MfUsgOc.load(str(path), _g05_dis(oc.parent.model_ws))
    return path.read_text().splitlines()[1:], again


def _real(action):
    return np.float32(float(action.split()[1].upper().replace("D", "E")))


# solver words of SGWF2BAS7N; NUMTRACK is an INTEGER, the others are REAL
_VALUES = {
    "HCLOSE": "1.23456789E-5",
    "BTOL": "1.2345678",
    "THETA": "0.85",
    "AKAPPA": "1.0E-4",
    "GAMMA": "5.0d-2",
    "AMOMENTUM": ".01",
    "NUMTRACK": "3",
    "BREDUC": "0.25",
    "RESLIM": "100.",
    "TRUNCATEDNEWTONCUTOFF": "0.5",
}
_SEPS = {"blank": " ", "comma": ",", "tab": "\t", "comma_blanks": " , "}


@pytest.mark.parametrize("sep", list(_SEPS))
@pytest.mark.parametrize("name", list(_VALUES))
def test_oc_solver_words_kept(function_tmpdir, name, sep):
    """The solver words are read like URWORD (blanks, tabs or commas, any
    case) and keep their value text from a file, through the API and after a
    change of a loaded package; the action survives load -> write -> load."""
    value = _VALUES[name]
    action = f"{name[:3].lower()}{name[3:]}{_SEPS[sep]}{value}"
    expected = f"{name} {value.upper()}"
    oc = _load(function_tmpdir / "text", _text(p2=f"  {action}\n"))
    assert oc.stress_period_data[_P2][0] == expected
    lines, again = _rewrite(oc)
    assert f"  {expected}" in lines
    assert again.stress_period_data == oc.stress_period_data
    ml = _g05_dis(_g05_model(function_tmpdir / "api", _text()))
    oc = MfUsgOc(ml, stress_period_data={_P2: [action, "save head"]})
    lines, again = _rewrite(oc)
    assert f"  {action}" in lines
    assert again.stress_period_data[_P2] == [expected, "SAVE HEAD"]
    oc = _load(function_tmpdir / "changed", _text())
    oc.stress_period_data[_P2].append(action)
    lines, again = _rewrite(oc)
    assert again.stress_period_data[_P2][-1] == expected
    if name != "NUMTRACK":
        assert _real(again.stress_period_data[_P2][-1]) == _real(f"{name} {value}")


def test_oc_hclose_btol_full_precision(function_tmpdir):
    """HCLOSE and BTOL keep the REAL value of the file (they were rounded to
    five significant digits)."""
    p2 = "  HCLOSE 1.23456789E-5\n  BTOL 1.2345678\n"
    oc = _load(function_tmpdir, _text(p2=p2))
    lines, again = _rewrite(oc)
    for acts in (oc.stress_period_data[_P2], again.stress_period_data[_P2]):
        assert _real(acts[0]) == np.float32(1.23456789e-5)
        assert _real(acts[1]) == np.float32(1.2345678)
    assert "  HCLOSE 1.23456789E-5" in lines and "  BTOL 1.2345678" in lines


# a complete hand-written file in the form FloPy writes (one blank between
# words); also the control of the engine run
_HAND = (
    _HDR + "COMPACT BUDGET\n"
    "PERIOD 1 STEP 1\n"
    "  HCLOSE 1.23456789E-5\n"
    "  BTOL 1.2345678\n"
    "  MXITER 90\n"
    "  THETA 0.85\n"
    "  AKAPPA 1.0E-4\n"
    "  GAMMA 5.0D-2\n"
    "  AMOMENTUM .01\n"
    "  NUMTRACK 0\n"
    "  BREDUC 0.25\n"
    "  RESLIM 100.\n"
    "  TRUNCATEDNEWTONCUTOFF 0.5\n"
    "  SAVE HEAD\n"
    "  SAVE BUDGET\n"
    "  PRINT BUDGET\n"
    "PERIOD 2 STEP 5\n"
    "  THETA 0.9\n"
    "  TRUNCATEDNEWTONCUTOFF 0\n"
    "  SAVE HEAD\n"
    "  SAVE BUDGET\n"
    "  PRINT BUDGET\n"
)
# what USG-T echoes for _HAND, in order
_HAND_ECHO = [
    ("HCLOSE", 1.23456789e-5),
    ("BTOL", 1.2345678),
    ("MXITER", 90),
    ("THETA", 0.85),
    ("AKAPPA", 1.0e-4),
    ("GAMMA", 0.05),
    ("AMOMENTUM", 0.01),
    ("NUMTRACK", 0),
    ("BREDUC", 0.25),
    ("RESLIM", 100.0),
    ("TRUNCATEDNEWTONCUTOFF", 0.5),
    ("THETA", 0.9),
    ("TRUNCATEDNEWTONCUTOFF", 0.0),
]


def test_oc_rewrite_matches_hand_written_file(function_tmpdir):
    """load -> write of a hand-written file gives the same records, word by
    word: solver words, value text and COMPACT BUDGET without AUX."""
    oc = _load(function_tmpdir, _HAND)
    lines, again = _rewrite(oc)
    assert [ln.upper().split() for ln in lines] == [
        ln.split() for ln in _HAND.splitlines()[1:]
    ]
    assert again.stress_period_data == oc.stress_period_data


@pytest.mark.parametrize(
    "action, match",
    [
        ("THETA abc", "THETA value 'ABC': FloPy requires"),
        ("AKAPPA,NaN", "AKAPPA value 'NAN': FloPy requires"),
        ("gamma\tinf", "GAMMA value 'INF': FloPy requires"),
        ("RESLIM 1.0E39", "RESLIM value '1.0E39': FloPy requires"),
        ("TRUNCATEDNEWTONCUTOFF 0.5.1", "TRUNCATEDNEWTONCUTOFF value '0.5.1': FloPy"),
        ("BREDUC", "BREDUC needs a value"),
        ("amomentum ,", "AMOMENTUM needs a value"),
        ("NUMTRACK 1.5", "NUMTRACK value '1.5': FloPy requires"),
        ("NUMTRACK 99999999999", "NUMTRACK value '99999999999': FloPy requires"),
        ("NUMTRACK " + "0" * 30 + "3", "NUMTRACK value '0{30}3': FloPy requires"),
        ("NUMTRACK", "NUMTRACK needs a value"),
        ("HCLOSE 1e", "HCLOSE value '1E': FloPy requires"),
        ("btol", "BTOL needs a value"),
    ],
    ids=[
        "theta_text",
        "akappa_nan",
        "gamma_inf",
        "reslim_beyond_real",
        "cutoff_text",
        "breduc_missing",
        "amomentum_missing_comma",
        "numtrack_real",
        "numtrack_beyond_integer",
        "numtrack_31_characters",
        "numtrack_missing",
        "hclose_text",
        "btol_missing",
    ],
)
def test_oc_solver_words_rejected(function_tmpdir, action, match):
    """FloPy requires an explicit finite REAL value after a solver word (an
    INTEGER after NUMTRACK): other values, and words without a value, are
    rejected on load, and on write before the destination is opened (API and
    changed package). This is a guard of FloPy: USG-T reads some of them (a
    missing value as 0, '1e' as 1.0)."""
    with pytest.raises(ValueError, match=match):
        _load(function_tmpdir / "text", _text(p2=f"  {action}\n"))
    ml = _g05_dis(_g05_model(function_tmpdir / "api", _text()))
    api = MfUsgOc(ml, stress_period_data={_P2: [action.lower(), "save head"]})
    changed = _load(function_tmpdir / "changed", _text())
    changed.stress_period_data[_P2].append(action)
    for oc in (api, changed):
        oc.fn_path = str(Path(oc.parent.model_ws) / "dest.oc")
        p, mtime = _g04_sentinel(oc)
        with pytest.raises(ValueError, match=match):
            oc.write_file()
        assert p.read_text() == "SENTINEL\n" and p.stat().st_mtime_ns == mtime


def test_oc_numtrack_30_characters_kept(function_tmpdir):
    """A NUMTRACK value of 30 characters, the longest word URWORD converts,
    is kept."""
    action = "NUMTRACK " + "0" * 29 + "3"
    oc = _load(function_tmpdir, _text(p2=f"  {action}\n"))
    assert oc.stress_period_data[_P2][0] == action
    lines, again = _rewrite(oc)
    assert f"  {action}" in lines
    assert again.stress_period_data == oc.stress_period_data


@pytest.mark.parametrize(
    "first", ["period 1 step 1", "period 2 step 4"], ids=["period_1", "period_2"]
)
@pytest.mark.parametrize("name", list(_VALUES))
def test_oc_solver_word_before_first_period_rejected(function_tmpdir, name, first):
    """USG-T reads the solver words after a PERIOD record only (it stops at
    one before the first PERIOD): load raises instead of dropping the word or
    moving it to the first time step."""
    text = _HDR + f"COMPACT BUDGET AUX\n  {name} {_VALUES[name]}\n{first}\n" + _ACT
    with pytest.raises(ValueError, match=f"starts with '{name}'"):
        _load(function_tmpdir, text)


def _check_warnings(ws, actions):
    ml = _g05_dis(_g05_model(ws, _text()))
    oc = MfUsgOc(ml, stress_period_data={(0, 0): actions})
    return [str(d) for d in oc.check(verbose=False).summary_array["desc"]]


def test_oc_check_recognizes_solver_words(function_tmpdir):
    """check() takes the solver words with a value as valid actions and
    warns when the value is missing."""
    actions = [f"{name.lower()} {value}" for name, value in _VALUES.items()]
    assert _check_warnings(function_tmpdir / "ok", actions + ["save head"]) == []
    for name in _VALUES:
        descs = _check_warnings(function_tmpdir / name, ["save head", name.lower()])
        assert any(name in d and "requires a value" in d for d in descs), (name, descs)


@pytest.mark.parametrize(
    "line, written",
    [
        ("COMPACT BUDGET", "COMPACT BUDGET"),
        ("compact budget", "COMPACT BUDGET"),
        ("COMPACT,BUDGET", "COMPACT BUDGET"),
        ("COMPACT BUDGET FILES", "COMPACT BUDGET"),
        ("COMPACT BUDGET AUX", "COMPACT BUDGET AUX"),
        ("compact budget auxiliary", "COMPACT BUDGET AUX"),
        ("Compact,Budget,Aux", "COMPACT BUDGET AUX"),
        ("COMPACT\tBUDGET\tAUXILIARY", "COMPACT BUDGET AUX"),
        ("COMPACT BUDGET AUX\nCOMPACT BUDGET", "COMPACT BUDGET AUX"),
        ("COMPACT BUDGET\nCOMPACT BUDGET AUX", "COMPACT BUDGET AUX"),
    ],
    ids=[
        "no_aux",
        "no_aux_lower",
        "no_aux_comma",
        "no_aux_files",
        "aux",
        "auxiliary_lower",
        "aux_comma",
        "auxiliary_tab",
        "aux_then_no_aux",
        "no_aux_then_aux",
    ],
)
def test_oc_compact_budget_aux_roundtrip(function_tmpdir, line, written):
    """COMPACT BUDGET is written with AUX only when the file has AUX or
    AUXILIARY as its third word (in any COMPACT record: USG-T does not turn
    IAUXSV off again); the form survives load -> write -> load."""
    oc = _load(function_tmpdir, _text(compact=line + "\n"))
    lines, again = _rewrite(oc)
    assert [ln for ln in lines if ln.startswith("COMPACT")] == [written]
    lines, _ = _rewrite(again, "again.oc")
    assert [ln for ln in lines if ln.startswith("COMPACT")] == [written]


def test_oc_compact_aux_default_and_keyword(function_tmpdir):
    """COMPACT BUDGET AUX stays the default; the compact_aux keyword writes
    COMPACT BUDGET without AUX. The constructor signature does not change."""
    ml = _g05_dis(_g05_model(function_tmpdir, _text()))
    spd = {(0, 0): ["save head", "save budget"]}
    oc = MfUsgOc(ml, stress_period_data=spd, compact_aux=False)
    lines, again = _rewrite(oc)
    assert "COMPACT BUDGET" in lines and "COMPACT BUDGET AUX" not in lines
    assert oc.compact_aux is False and again.compact_aux is False
    oc = MfUsgOc(ml, stress_period_data=spd)
    lines, again = _rewrite(oc, "default.oc")
    assert "COMPACT BUDGET AUX" in lines
    assert oc.compact_aux is True and again.compact_aux is True
    oc = MfUsgOc(ml, stress_period_data=spd, compact=False, compact_aux=False)
    lines, _ = _rewrite(oc, "none.oc")
    assert not [ln for ln in lines if "COMPACT" in ln]
    assert "compact_aux" not in inspect.signature(MfUsgOc.__init__).parameters


@pytest.mark.parametrize(
    "text, match",
    [
        (_text(p2="  FROBNICATE 1\n"), "FROBNICATE"),
        (_text(p2="  5 6\n"), "'5'"),
        (_text(p2=",,,\n"), "''"),
        (_text(p2="  ATSA\n"), "ATSA"),
        (_text(p2="  FASTFORWARD 1 1 60\n"), "FASTFORWARD"),
        (_text(header="FOO BAR\n"), "FOO"),
        (_text(header="12\n"), "'12'"),
        (_text(header="NPSTPS 4\n"), "NPSTPS"),
    ],
    ids=[
        "period_word",
        "period_integer",
        "period_commas",
        "period_atsa",
        "period_fastforward",
        "header_word",
        "header_integer",
        "header_option",
    ],
)
def test_oc_unknown_word_rejected_on_load(function_tmpdir, text, match):
    """A record whose first word USG-T does not accept (it stops with ERROR
    READING OUTPUT CONTROL INPUT DATA) is not dropped: load raises."""
    with pytest.raises(ValueError, match=match):
        _load(function_tmpdir, text)


def test_oc_skipped_header_lines_still_load(function_tmpdir):
    """Lines that start with ATSA, FASTFORWARD or FASTFORWARDC before the
    first PERIOD are skipped like USG-T does (they are options only on the
    first record)."""
    plain = _load(function_tmpdir / "plain", _text())
    header = "ATSA\nFASTFORWARD 1 1 60\nfastforwardc,1,1,61\n"
    oc = _load(function_tmpdir / "skipped", _text(header=header))
    assert oc.atsa == 0
    assert oc.stress_period_data == plain.stress_period_data


# engine runs: the hand-written file against the file FloPy writes

_ECHO = re.compile(
    r"^ {4}(HCLOSE|BTOL|MXITER|THETA|AKAPPA|GAMMA|AMOMENTUM|NUMTRACK|BREDUC|RESLIM"
    r"|TRUNCATEDNEWTONCUTOFF) *(\S+) *$",
    re.M,
)
# a GHB with an auxiliary variable: USG-T WEL saves its auxiliary data whatever
# IAUXSV is (gwf2wel7u1.f 505-515, 604-606), GHB only with AUX (300-310, 379-380)
_GHB_AUX = "# GHB\n 1 50 AUX IFACE\n 1 0\n 1 1 3 15.0 0.5 2\n -1 0\n"


def _run(ws, oc_text, change=None, ghb=None):
    """Run the 5-cell model with oc_text (loaded and written by FloPy when
    change is given; change(oc) may modify the package first)."""
    ws = _g05_model(ws, oc_text)
    if ghb is not None:
        (ws / "t.ghb").write_text(ghb)
        with open(ws / "t.nam", "a") as f:
            f.write("GHB 23 t.ghb\n")
    if change is not None:
        oc = MfUsgOc.load(str(ws / "t.oc"), _g05_dis(ws))
        change(oc)
        oc.fn_path = str(ws / "t.oc")
        oc.write_file()
    out = _g05x_run(ws)
    lst = (ws / "t.list").read_text()
    out["echo"] = _ECHO.findall(lst)
    out["aux_echo"] = lst.count("AUXILIARY DATA WILL BE SAVED IN CELL-BY-CELL")
    major, minor = re.search(r"USG-TRANSPORT VERSION\s+(\d+)\.(\d+)", lst).groups()
    out["version"] = (int(major), int(minor))
    out["ghb"] = []
    if out["normal"] and ghb is not None:
        cbc = CellBudgetFile(str(ws / "t.cbc"))
        out["ghb"] = [
            (r.dtype.names, r.tolist()) for r in cbc.get_data(text="HEAD DEP BOUNDS")
        ]
    return out


def _same(a, b):
    for key in ("normal", "failed", "pct", "budget_at", "delt", "hds", "cbc"):
        assert a[key] == b[key], key
    for key in ("echo", "aux_echo", "ghb"):
        assert a[key] == b[key], key
    np.testing.assert_array_equal(a["heads"], b["heads"])


def _unchanged(oc):
    pass


def _replace(old, new):
    def change(oc):
        for acts in oc.stress_period_data.values():
            acts[:] = [new if a == old else a for a in acts]

    return change


@requires_usgt
def test_usgt_exe_oc_solver_words_rewrite_matches_control(function_tmpdir):
    """USG-T 2.8 accepts the solver words of a hand-written OC file and
    echoes their values; the file FloPy writes after loading it runs the
    same; a changed value changes the echo."""
    plain = _run(function_tmpdir / "plain", _text())
    assert plain["normal"] and plain["echo"] == []
    if plain["version"] < (2, 8):
        pytest.skip(
            "USG-T before 2.8 reads HCLOSE, BTOL and MXITER only and stops at "
            "the other solver words"
        )
    control = _run(function_tmpdir / "control", _HAND)
    assert control["normal"] and not control["failed"] and control["hds"]
    assert [name for name, _ in control["echo"]] == [name for name, _ in _HAND_ECHO]
    for (_, got), (name, value) in zip(control["echo"], _HAND_ECHO):
        assert float(got) == pytest.approx(value, rel=5e-6, abs=0), name
    _same(_run(function_tmpdir / "rewrite", _HAND, _unchanged), control)
    changed = _run(
        function_tmpdir / "changed", _HAND, _replace("THETA 0.85", "THETA 0.6")
    )
    assert changed["normal"]
    assert [float(v) for n, v in changed["echo"] if n == "THETA"] == [
        pytest.approx(0.6, rel=5e-6),
        pytest.approx(0.9, rel=5e-6),
    ]
    cutoff = _run(
        function_tmpdir / "cutoff",
        _HAND,
        _replace("TRUNCATEDNEWTONCUTOFF 0.5", "TRUNCATEDNEWTONCUTOFF 0.75"),
    )
    assert cutoff["normal"]
    assert [float(v) for n, v in cutoff["echo"] if n == "TRUNCATEDNEWTONCUTOFF"] == [
        0.75,
        0.0,
    ]


@requires_usgt
def test_usgt_exe_oc_hclose_btol_precision(function_tmpdir):
    """The HCLOSE and BTOL values USG-T echoes for the file FloPy writes are
    those of the hand-written file (six digits in the list file)."""
    text = _text(p2="  HCLOSE 1.23456E-5\n  BTOL 1.2345678\n")
    control = _run(function_tmpdir / "control", text)
    assert control["normal"] and not control["failed"] and control["hds"]
    assert [name for name, _ in control["echo"]] == ["HCLOSE", "BTOL"]
    assert float(control["echo"][0][1]) == pytest.approx(1.23456e-5, rel=1e-9)
    assert float(control["echo"][1][1]) == pytest.approx(1.23457, rel=1e-9)
    _same(_run(function_tmpdir / "rewrite", text, _unchanged), control)
    changed = _run(
        function_tmpdir / "changed",
        text,
        _replace("HCLOSE 1.23456E-5", "HCLOSE 2.5E-5"),
    )
    assert changed["normal"]
    assert float(changed["echo"][0][1]) == pytest.approx(2.5e-5, rel=1e-9)


@requires_usgt
def test_usgt_exe_oc_compact_budget_aux(function_tmpdir):
    """COMPACT BUDGET without AUX: the file FloPy writes saves no auxiliary
    data in the budget file, like the hand-written file; with AUX (file or
    compact_aux) the GHB records carry IFACE."""
    no_aux = _text(compact="COMPACT BUDGET\n")
    control = _run(function_tmpdir / "control", no_aux, ghb=_GHB_AUX)
    assert control["normal"] and control["aux_echo"] == 0 and control["ghb"]
    assert all("IFACE" not in names for names, _ in control["ghb"])
    _same(_run(function_tmpdir / "rewrite", no_aux, _unchanged, _GHB_AUX), control)
    with_aux = _run(function_tmpdir / "control_aux", _text(), ghb=_GHB_AUX)
    assert with_aux["normal"] and with_aux["aux_echo"] == 1
    assert all("IFACE" in names for names, _ in with_aux["ghb"])
    _same(
        _run(function_tmpdir / "rewrite_aux", _text(), _unchanged, _GHB_AUX), with_aux
    )

    def set_aux(oc):
        oc.compact_aux = True

    _same(_run(function_tmpdir / "changed", no_aux, set_aux, _GHB_AUX), with_aux)
