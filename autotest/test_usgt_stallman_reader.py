"""Tests of the output-reading helpers of the USG-Transport Stallman heat
tutorial (mfusg_transport_tutorial06_StallmanHeatTransport.py).

The helpers are loaded from the tutorial itself: only its import statements
and its function definitions are executed (selected with ast), so the model
cells, which need the executables, are not run. The binary files are written
here with struct, independently of the FloPy writers, in the USG-T layout:
kstp, kper (int32), pertim, totim (real), text (16 bytes), nstrt, nend,
ilay (int32), followed by nend - nstrt + 1 reals.
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
    / "mfusg_transport_tutorial06_StallmanHeatTransport.py"
)
HELPERS = {
    "new_run_root",
    "write_and_run",
    "check_records",
    "solute_label",
    "read_profile",
}

# (kper, kstp, totim): the tutorial plots kstp 6 of stress period 601
STEPS = [(1, 1, 1.0), (601, 5, 2.0), (601, 6, 3.0)]
LAST_STEP = (5, 600)  # zero-based (kstp, kper)
LAYERS = ((1, 1), (2, 1), (3, 1), (4, 1))  # (ilay, number of nodes)
BASE = {"CONC": 100.0, "CONC01": 100.0, "CONC02": 200.0, "TMPR": 300.0}


@pytest.fixture(scope="module")
def tut():
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
    return SimpleNamespace(**{name: namespace[name] for name in HELPERS})


def value(label, step, ilay, node, offset=0.0):
    return BASE.get(label, 900.0) + 10 * ilay + step + 0.1 * node + offset


def records(labels, precision="single", layers=LAYERS, steps=STEPS, offset=0.0):
    real = "f" if precision == "single" else "d"
    out = b""
    for step, (kper, kstp, totim) in enumerate(steps, 1):
        for label in labels:
            nstrt = 1
            for ilay, n in layers:
                nend = nstrt + n - 1
                out += struct.pack(
                    f"<2i2{real}16s3i",
                    kstp,
                    kper,
                    totim,
                    totim,
                    f"{label:>16}".encode(),
                    nstrt,
                    nend,
                    ilay,
                )
                vals = [value(label, step, ilay, j, offset) for j in range(n)]
                out += struct.pack(f"<{n}{real}", *vals)
                nstrt = nend + 1
    return out


def expected(label, step=3, layers=LAYERS):
    return [value(label, step, ilay, j) for ilay, n in layers for j in range(n)]


def write(path, data):
    path.write_bytes(data)
    return str(path)


@pytest.mark.parametrize(
    "labels, solute",
    [
        (["CONC", "TMPR"], "CONC"),
        (["CONC01", "TMPR"], "CONC01"),
        (["CONC"], "CONC"),
        (["CONC01", "CONC02", "TMPR"], "CONC01"),
    ],
)
def test_solute_label_and_profiles(tut, tmp_path, labels, solute):
    fname = write(tmp_path / "model.con", records(labels))
    assert tut.solute_label(fname) == solute
    conc = tut.read_profile(fname, solute, LAST_STEP, 4)
    assert conc.shape == (4,)
    np.testing.assert_allclose(conc, expected(solute), rtol=1e-6)
    if "TMPR" in labels:
        tmpr = tut.read_profile(fname, "TMPR", LAST_STEP, 4)
        np.testing.assert_allclose(tmpr, expected("TMPR"), rtol=1e-6)


@pytest.mark.parametrize(
    "labels, match",
    [
        (["CONC", "CONC01", "TMPR"], r"found \['CONC', 'CONC01'\]"),
        (["TMPR"], r"found \[\]"),
        (["CONC02", "TMPR"], r"found \[\]"),
    ],
)
def test_solute_label_not_unique(tut, tmp_path, labels, match):
    fname = write(tmp_path / "model.con", records(labels))
    with pytest.raises(ValueError, match=match):
        tut.solute_label(fname)


@pytest.mark.parametrize(
    "labels, label, other",
    [
        (["CONC", "CONCX"], "CONC", "CONCX"),
        (["CONC01", "CONC010"], "CONC01", "CONC010"),
    ],
)
def test_accidental_prefix_match(tut, tmp_path, labels, label, other):
    fname = write(tmp_path / "model.con", records(labels))
    assert tut.solute_label(fname) == label
    with pytest.raises(ValueError, match=f"also selects \\['{other}'\\]"):
        tut.read_profile(fname, label, LAST_STEP, 4)


@pytest.mark.parametrize("precision", ["single", "double"])
def test_precision(tut, tmp_path, precision):
    fname = write(tmp_path / "model.con", records(["CONC01", "TMPR"], precision))
    profile = tut.read_profile(fname, "CONC01", LAST_STEP, 4)
    np.testing.assert_allclose(
        profile, expected("CONC01"), rtol=1e-12 if precision == "double" else 1e-6
    )


def test_irregular_layers(tut, tmp_path):
    layers = ((1, 3), (2, 2), (3, 1))
    fname = write(tmp_path / "model.con", records(["CONC01", "CONC02"], layers=layers))
    profile = tut.read_profile(fname, "CONC02", LAST_STEP, 6)
    np.testing.assert_allclose(profile, expected("CONC02", layers=layers), rtol=1e-6)
    with pytest.raises(ValueError, match="has 6 values"):
        tut.read_profile(fname, "CONC02", LAST_STEP, 3)


def test_explicit_time_step(tut, tmp_path):
    fname = write(tmp_path / "model.con", records(["CONC01", "TMPR"]))
    # zero-based (4, 600) is kstp 5 of stress period 601 (step 2 in STEPS)
    profile = tut.read_profile(fname, "TMPR", (4, 600), 4)
    np.testing.assert_allclose(profile, expected("TMPR", step=2), rtol=1e-6)
    with pytest.raises(ValueError, match="kstpkper not found"):
        tut.read_profile(fname, "TMPR", (6, 600), 4)


def test_non_finite_profile(tut, tmp_path):
    data = bytearray(records(["CONC"]))
    data[-4:] = struct.pack("<f", np.nan)
    fname = write(tmp_path / "model.con", bytes(data))
    with pytest.raises(ValueError, match="non-finite"):
        tut.read_profile(fname, "CONC", LAST_STEP, 4)


def test_truncated_header(tut, tmp_path):
    data = records(["CONC01", "TMPR"])
    fname = write(tmp_path / "model.con", data[: len(data) - 4 - 30])
    with pytest.raises(ValueError, match=r"model\.con: truncated record header"):
        tut.solute_label(fname)


def test_truncated_last_payload_of_other_label(tut, tmp_path):
    # the last record is TMPR; the CONC01 records are complete
    data = records(["CONC01", "TMPR"])
    fname = write(tmp_path / "model.con", data[:-2])
    with pytest.raises(ValueError, match=r"last record \(TMPR.*is truncated"):
        tut.read_profile(fname, "CONC01", LAST_STEP, 4)
    with pytest.raises(ValueError, match="is truncated"):
        tut.solute_label(fname)


@pytest.mark.parametrize("label", ["TMPR", "CONC"])
def test_stale_records_going_back_in_time(tut, tmp_path, label):
    # a shorter run written over a longer, earlier heat and solute run
    old = records(["CONC01", "TMPR"], offset=5000.0)
    new = records([label])
    fname = write(tmp_path / "model.con", new + old[len(new) :])
    with pytest.raises(ValueError, match=r"totim goes back from 3\.0 to"):
        tut.read_profile(fname, label, LAST_STEP, 4)


def test_equal_times_within_step_accepted(tut, tmp_path):
    fname = write(tmp_path / "model.con", records(["CONC01", "CONC02", "TMPR"]))
    assert tut.check_records(fname) == {"CONC01", "CONC02", "TMPR"}


def test_new_run_root(tut, tmp_path):
    root1, ws1 = tut.new_run_root(tmp_path)
    root2, ws2 = tut.new_run_root(tmp_path)
    assert root1 != root2
    for root, ws in ((root1, ws1), (root2, ws2)):
        assert Path(root).parent == tmp_path
        assert sorted(ws) == ["heat", "heat_solute", "solute"]
        assert {Path(p).parent for p in ws.values()} == {Path(root)}
        assert all(os.path.isdir(p) and not os.listdir(p) for p in ws.values())
    assert len(set(ws1.values()) | set(ws2.values())) == 6


def test_write_and_run(tut):
    calls = []

    class Model:
        name = "Ex6_Stallman"
        model_ws = "ws"

        def __init__(self, success):
            self.success = success

        def write_input(self):
            calls.append("write")

        def run_model(self, **kwargs):
            calls.append("run")
            return self.success, ["line 1", "stopped"]

    tut.write_and_run(Model(True))
    assert calls == ["write", "run"]
    with pytest.raises(RuntimeError, match="did not terminate normally in ws"):
        tut.write_and_run(Model(False))
