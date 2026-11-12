"""Tests of the run and output-reading helpers of the USG-Transport conduit
tutorial (mfusg_transport_tutorial03_Conduit.py).

The helpers are loaded from the tutorial itself: only its import statements
and its function definitions are executed (selected with ast), so the model
cells, which need the executables, are not run. The concentration files are
written here with struct, independently of the FloPy writers, in the
structured HeadFile layout: kstp, kper (int32), pertim, totim (real),
text (16 bytes), ncol, nrow, ilay (int32), followed by nrow * ncol reals.
"""

import ast
import os
import struct
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest

TUTORIAL = (
    Path(__file__).parents[1]
    / ".docs"
    / "Notebooks"
    / "mfusg_transport_tutorial03_Conduit.py"
)
HELPERS = {"new_run_root", "write_and_run", "output_calendar", "read_last_conc"}

NLAY, NROW, NCOL = 2, 3, 4
CALENDAR = [(1, 1, 10.0), (2, 1, 20.0), (3, 1, 30.0)]  # (kstp, kper, totim)
NORMAL = " Normal termination of simulation"


@pytest.fixture(scope="module")
def namespace():
    tree = ast.parse(TUTORIAL.read_text())
    body = [
        node
        for node in tree.body
        if isinstance(node, (ast.Import, ast.ImportFrom))
        or (isinstance(node, ast.FunctionDef) and node.name in HELPERS)
    ]
    missing = HELPERS - {n.name for n in body if isinstance(n, ast.FunctionDef)}
    assert not missing, f"helpers not found in {TUTORIAL.name}: {sorted(missing)}"
    namespace = {}
    exec(compile(ast.Module(body, []), str(TUTORIAL), "exec"), namespace)
    return namespace


@pytest.fixture
def tut(namespace):
    return SimpleNamespace(**{name: namespace[name] for name in HELPERS})


@pytest.fixture
def opened(namespace, monkeypatch):
    """HeadFile objects created by read_last_conc."""
    readers = []

    class TrackedHeadFile(namespace["HeadFile"]):
        def __init__(self, *args, **kwargs):
            super().__init__(*args, **kwargs)
            readers.append(self)

    monkeypatch.setitem(namespace, "HeadFile", TrackedHeadFile)
    return readers


def value(step, ilay, offset=0.0, ncol=NCOL):
    cells = np.arange(NROW * ncol, dtype=float).reshape(NROW, ncol)
    return 100.0 * step + 10.0 * ilay + 0.1 * cells + offset


def records(
    calendar=CALENDAR,
    layers=range(1, NLAY + 1),
    label="CONC",
    ncol=NCOL,
    precision="single",
    offset=0.0,
):
    real = "f" if precision == "single" else "d"
    out = b""
    for step, (kstp, kper, totim) in enumerate(calendar, 1):
        for ilay in layers:
            out += struct.pack(
                f"<2i2{real}16s3i",
                kstp,
                kper,
                totim,
                totim,
                f"{label:>16}".encode(),
                ncol,
                NROW,
                ilay,
            )
            vals = value(step, ilay, offset, ncol).ravel()
            out += struct.pack(f"<{NROW * ncol}{real}", *vals)
    return out


def write(path, data):
    path.write_bytes(data)
    return str(path)


def read(tut, fname, nlay=NLAY, calendar=CALENDAR):
    return tut.read_last_conc(fname, nlay, NROW, NCOL, calendar)


# new_run_root


def test_new_run_root(tut, tmp_path):
    root1, ws1 = tut.new_run_root(tmp_path)
    root2, ws2 = tut.new_run_root(tmp_path)
    assert root1 != root2
    for root, ws in ((root1, ws1), (root2, ws2)):
        assert Path(root).parent == tmp_path
        assert Path(root).name.startswith("Ex3_CLN_Conduit_")
        assert sorted(ws) == ["dispersion", "no_dispersion"]
        assert {Path(p).parent for p in ws.values()} == {Path(root)}
        assert all(os.path.isdir(p) and not os.listdir(p) for p in ws.values())
    assert len(set(ws1.values()) | set(ws2.values())) == 4


# write_and_run


class Model:
    name = "Ex3_CLN_Conduit"
    exe_name = "/path/to/usgt"

    def __init__(self, model_ws, success=True, buff=None):
        self.model_ws = str(model_ws)
        self.success = success
        self.buff = ["  USG-TRANSPORT VERSION 9.9.9", NORMAL] if buff is None else buff
        self.calls = []

    def write_input(self):
        self.calls.append("write")
        Path(self.model_ws, f"{self.name}.nam").write_text("")

    def run_model(self, **kwargs):
        self.calls.append(("run", kwargs))
        return self.success, self.buff


def test_write_and_run(tut, tmp_path, capsys):
    mf = Model(tmp_path)
    tut.write_and_run(mf, "no_dispersion")
    assert mf.calls == ["write", ("run", {"silent": True, "report": True})]
    out = capsys.readouterr().out
    assert out == (
        "no_dispersion: USG-Transport version 9.9.9; "
        f"Normal termination of simulation; {tmp_path}\n"
    )


def test_write_and_run_version_not_reported(tut, tmp_path, capsys):
    tut.write_and_run(Model(tmp_path, buff=[NORMAL]), "dispersion")
    assert "USG-Transport version not reported;" in capsys.readouterr().out


@pytest.mark.parametrize(
    "success, buff, nnormal",
    [
        (False, [NORMAL], 1),
        (True, [" ABNORMAL TERMINATION OF SIMULATION"], 0),
        (True, [" Failure to converge"], 0),
        (True, [], 0),
        (True, [NORMAL, NORMAL], 2),
        (True, [" Normal termination of simulation, see listing"], 0),
    ],
    ids=["no_success", "abnormal", "not_converged", "no_output", "two", "inexact"],
)
def test_write_and_run_failure(tut, tmp_path, success, buff, nnormal):
    mf = Model(tmp_path, success, buff + ["last line"])
    with pytest.raises(RuntimeError) as err:
        tut.write_and_run(mf, "dispersion")
    msg = str(err.value)
    assert msg.startswith(
        f"dispersion: Ex3_CLN_Conduit did not terminate normally in {tmp_path} "
        f"(executable /path/to/usgt, {nnormal} normal termination lines):"
    )
    assert msg.endswith("last line")
    assert [c if c == "write" else c[0] for c in mf.calls] == ["write", "run"]


def test_write_and_run_workspace_not_empty(tut, tmp_path):
    (tmp_path / "Ex3_CLN_Conduit.con").write_bytes(b"old output")
    mf = Model(tmp_path)
    with pytest.raises(RuntimeError, match=r"no_dispersion: workspace .* is not empty"):
        tut.write_and_run(mf, "no_dispersion")
    assert mf.calls == []


# output_calendar


def test_output_calendar_of_tutorial_dis(tut):
    import flopy

    m = flopy.modflow.Modflow()
    dis = flopy.modflow.ModflowDis(m, 2, 3, 4, nper=1, perlen=3000.0, nstp=100)
    calendar = tut.output_calendar(dis)
    assert len(calendar) == 100
    assert calendar[0] == (1, 1, 30.0) and calendar[-1] == (100, 1, 3000.0)
    assert [c[0] for c in calendar] == list(range(1, 101))


def test_output_calendar_periods(tut):
    import flopy

    m = flopy.modflow.Modflow()
    dis = flopy.modflow.ModflowDis(
        m, 1, 1, 1, nper=2, perlen=[1.0, 7.0], nstp=[1, 3], tsmult=[1.0, 2.0]
    )
    assert tut.output_calendar(dis) == [
        (1, 1, 1.0),
        (1, 2, 2.0),
        (2, 2, 4.0),
        (3, 2, 8.0),
    ]


# read_last_conc


@pytest.mark.parametrize("precision", ["single", "double"])
def test_read_last_conc(tut, tmp_path, opened, precision):
    fname = write(tmp_path / "m.con", records(precision=precision))
    data = read(tut, fname)
    assert data.shape == (NLAY, NROW, NCOL)
    for k in range(NLAY):
        np.testing.assert_allclose(data[k], value(3, k + 1), rtol=1e-6)
    assert len(opened) == 1 and opened[0].file.closed


def test_read_last_conc_time_within_stored_precision(tut, tmp_path):
    last = float(np.nextafter(np.float32(30.0), np.float32(31.0)))
    calendar = CALENDAR[:2] + [(3, 1, last)]
    read(tut, write(tmp_path / "m.con", records(calendar)))


def test_read_last_conc_missing(tut, tmp_path):
    with pytest.raises(FileNotFoundError, match=r"m\.con: concentration file not"):
        read(tut, str(tmp_path / "m.con"))


def test_read_last_conc_truncated_header(tut, tmp_path):
    data = records()
    fname = write(tmp_path / "m.con", data[: len(data) - 4 * NROW * NCOL - 20])
    with pytest.raises(ValueError, match=r"m\.con: truncated record header"):
        read(tut, fname)


def shifted(calendar, dt=10.0):
    return [(kstp, kper, totim + dt) for kstp, kper, totim in calendar]


# (file contents, error message) for files that read_last_conc must reject
# after opening them
BAD = {
    "truncated_payload": (records()[:-2], r"bytes; expected \d+ \(truncated"),
    "other_label_appended": (
        records() + records(CALENDAR[-1:], label="HEAD"),
        r"bytes; expected \d+ \(truncated record or other data",
    ),
    "label_conc01": (records(label="CONC01"), r"text 'conc' selects \['CONC01'\]"),
    "labels_conc_and_conc01": (
        records() + records(CALENDAR[-1:], label="CONC01"),
        r"selects \['CONC', 'CONC01'\]",
    ),
    "time_step_missing": (records(CALENDAR[:2]), "4 CONC records; expected 6"),
    # a shorter run written over the file of a longer, earlier run
    "stale_time_step_after_new": (
        records() + records([(4, 1, 40.0)], offset=5000.0),
        "8 CONC records; expected 6",
    ),
    "calendar_shifted": (
        records(shifted(CALENDAR)),
        r"record 1 has .* \(1, 1, 1, 4, 3, 20\.0\); expected \(1, 1, 1, 4, 3, 10\.0\)",
    ),
    "time_off_by_1e-3": (
        records(CALENDAR[:2] + [(3, 1, 30.001)]),
        r"record 5 has .*30\.00099",
    ),
    "time_going_back": (
        records([CALENDAR[0], CALENDAR[2], CALENDAR[1]]),
        r"record 3 has \([^)]*\) \(3, 1, 1, 4, 3, 30\.0\); "
        r"expected \(2, 1, 1, 4, 3, 20\.0\)",
    ),
    "layer_duplicated_replaces_other": (
        records(layers=(1, 1)),
        r"record 2 has \([^)]*\) \(1, 1, 1, 4, 3, 10\.0\); "
        r"expected \(1, 1, 2, 4, 3, 10\.0\)",
    ),
    "layer_missing": (records(layers=(1,)), r"record 2 has \([^)]*\) \(2, 1, 1,"),
    "layers_swapped": (records(layers=(2, 1)), r"record 1 has \([^)]*\) \(1, 1, 2,"),
    "ncol_differs": (
        records(ncol=NCOL - 1),
        r"record 1 has \([^)]*\) \(1, 1, 1, 3, 3,",
    ),
}


@pytest.mark.parametrize("case", sorted(BAD))
def test_read_last_conc_rejects(tut, tmp_path, opened, case):
    data, match = BAD[case]
    fname = write(tmp_path / "m.con", data)
    with pytest.raises(ValueError, match=match):
        read(tut, fname)
    assert len(opened) == 1 and opened[0].file.closed


@pytest.mark.parametrize(
    "nlay, calendar, match",
    [
        (1, CALENDAR, r"record 2 has \([^)]*\) \(1, 1, 2,"),
        (3, CALENDAR, r"record 3 has \([^)]*\) \(2, 1, 1,"),
        (NLAY, CALENDAR[:2], "6 CONC records; expected 4"),
    ],
    ids=["fewer_layers", "more_layers", "fewer_time_steps"],
)
def test_read_last_conc_other_model(tut, tmp_path, opened, nlay, calendar, match):
    fname = write(tmp_path / "m.con", records())
    with pytest.raises(ValueError, match=match):
        read(tut, fname, nlay, calendar)
    assert len(opened) == 1 and opened[0].file.closed


@pytest.mark.parametrize("bad", [np.nan, np.inf])
def test_read_last_conc_non_finite(tut, tmp_path, opened, bad):
    data = bytearray(records())
    data[-4:] = struct.pack("<f", bad)
    fname = write(tmp_path / "m.con", bytes(data))
    with pytest.raises(ValueError, match=r"last time step has shape .* non-finite"):
        read(tut, fname)
    assert len(opened) == 1 and opened[0].file.closed


def test_read_last_conc_cannot_tell_complete_stale_records(tut, tmp_path):
    # complete records of an earlier run with the same calendar look valid:
    # only the new, empty workspace of each part (write_and_run) excludes them
    old = records(offset=5000.0)
    fname = write(tmp_path / "m.con", old)
    np.testing.assert_allclose(read(tut, fname)[0], value(3, 1, 5000.0), rtol=1e-6)
