"""Tests of the conduit tutorial (tutorial03) checks of the CI script
.github/scripts/usgt_engine.py (snapshot-t03, report-t03) and regression
tests of its earlier commands (snapshot, report, fetch, install-usgt).

The script is loaded from its file. Run folders, listing files, concentration
files and executed notebooks are written here with struct and json; the
"executables" are inert files that are hashed, never run. Most tests use a
small grid and calendar (patched into the module); one end-to-end test runs
the script in a separate Python process with the tutorial's real dimensions.
"""

import ast
import importlib.util
import json
import ntpath
import os
import struct
import subprocess
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).parents[1]
SCRIPT = ROOT / ".github" / "scripts" / "usgt_engine.py"
TUTORIAL = ROOT / ".docs" / "Notebooks" / "mfusg_transport_tutorial03_Conduit.py"
NORMAL = "Normal termination of simulation"
SHAPE = (2, 2, 3)  # nlay, nrow, ncol
CALENDAR = ((1, 1, 1.0), (2, 1, 2.0), (3, 1, 3.0))
PARTS = ("no_dispersion", "dispersion")


@pytest.fixture
def ue(monkeypatch):
    spec = importlib.util.spec_from_file_location("usgt_engine_t03", SCRIPT)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    monkeypatch.setattr(mod, "T03_SHAPE", SHAPE)
    monkeypatch.setattr(mod, "T03_CALENDAR", CALENDAR)
    monkeypatch.setattr(sys, "argv", [str(SCRIPT)])
    return mod


def con_records(nlay=2, nrow=2, ncol=3, calendar=CALENDAR, label="CONC"):
    """[kstp, kper, pertim, totim, label, ncol, nrow, ilay, values] records."""
    return [
        [k, p, t, t, label, ncol, nrow, lay, [0.5] * (nrow * ncol)]
        for k, p, t in calendar
        for lay in range(1, nlay + 1)
    ]


def write_con(path, records, real="f"):
    with open(path, "wb") as f:
        for k, p, pt, t, label, nc, nr, lay, values in records:
            text = label.encode().rjust(6).ljust(16)
            f.write(struct.pack(f"<2i2{real}16s3i", k, p, pt, t, text, nc, nr, lay))
            f.write(struct.pack(f"<{len(values)}{real}", *values))


def make_part(ws, version="2.8.1", records=None, real="f", skip=()):
    ws.mkdir(parents=True)
    for ext in (
        "bas", "bcf", "bct", "chd", "cln", "dis", "nam", "oc", "sms", "wel",
        "hds", "cbc", "BCT.cbc", "clncb", "clncn", "clnhd",
    ):  # fmt: skip
        if f"Ex3_CLN_Conduit.{ext}" not in skip:
            (ws / f"Ex3_CLN_Conduit.{ext}").write_bytes(b"x")
    if "Ex3_CLN_Conduit.list" not in skip:
        lines = [f"  VERSION USG-TRANSPORT VERSION {v}\n" for v in version.split(",")]
        (ws / "Ex3_CLN_Conduit.list").write_text("".join(lines) + " Elapsed\n")
    if "Ex3_CLN_Conduit.con" not in skip:
        write_con(ws / "Ex3_CLN_Conduit.con", records or con_records(), real)


def make_root(runs, name="Ex3_CLN_Conduit_abc12345", parts=PARTS, **kw):
    root = runs / name
    for part in parts:
        make_part(root / part, **kw)
    return root


def lines_for(root, version="2.8.1", term=NORMAL):
    """stdout of the tutorial cells: run folder line, then one cell per part."""
    return [
        [f"Model output files are written in {root}\n"],
        [f"no_dispersion: USG-Transport version {version}; {term}; "
         f"{os.path.join(root, 'no_dispersion')}\n"],
        [f"dispersion: USG-Transport version {version}; {term}; "
         f"{os.path.join(root, 'dispersion')}\n"],
    ]  # fmt: skip


def write_notebook(runs, cells, name="mfusg_transport_tutorial03_Conduit.ipynb"):
    """cells: one list of stdout text fragments (or output dicts) per cell."""
    nb = {"cells": [], "nbformat": 4, "nbformat_minor": 5, "metadata": {}}
    for texts in cells:
        outs = [
            t if isinstance(t, dict) else
            {"output_type": "stream", "name": "stdout", "text": [t]}
            for t in texts
        ]  # fmt: skip
        nb["cells"].append({"cell_type": "code", "source": "", "outputs": outs})
    path = runs / name
    path.write_text(json.dumps(nb), encoding="utf-8")
    return path


def make_install(tmp_path, monkeypatch, env=True):
    exe = tmp_path / "bin" / "mfusgt"
    exe.parent.mkdir()
    exe.write_bytes(b"inert fixture, never run")
    import hashlib

    rec = {
        "command": "install-usgt",
        "errors": [],
        "exe": str(exe),
        "sha256": hashlib.sha256(exe.read_bytes()).hexdigest(),
        "zip_sha256": "0" * 64,
        "USGT_EXE": str(exe),
    }
    path = tmp_path / "install-usgt.json"
    path.write_text(json.dumps(rec))
    if env:
        monkeypatch.setenv("USGT_EXE", str(exe))
    else:
        monkeypatch.delenv("USGT_EXE", raising=False)
    return path, rec


class Case:
    """runs folder, snapshot, run folders, notebook, install and report."""

    def __init__(self, ue, tmp_path, monkeypatch, runs_name="Notebooks"):
        self.ue, self.tmp, self.mp = ue, tmp_path, monkeypatch
        self.runs = tmp_path / runs_name
        self.runs.mkdir()
        self.out = tmp_path / "usgt-engine"
        self.install, self.install_rec = make_install(tmp_path, monkeypatch)
        self.export = tmp_path / "usgt-t03-evidence"

    def snapshot(self):
        rc = self.ue.main(
            ["snapshot-t03", "--runs", str(self.runs), "--out", str(self.out)]
        )
        self.snap = json.loads((self.out / "snapshot-t03.json").read_text())
        return rc

    def touch(self, *paths, dt=1.0):
        t = self.snap["time"] + dt
        for p in paths:
            os.utime(p, (t, t))

    def report(self, export=True):
        argv = ["report-t03", "--runs", str(self.runs), "--out", str(self.out)]
        argv += ["--snapshot", str(self.out / "snapshot-t03.json")]
        argv += ["--install", str(self.install)]
        if export:
            argv += ["--export", str(self.export)]
        rc = self.ue.main(argv)
        self.rec = json.loads((self.out / "report-t03.json").read_text())
        return rc

    def valid(self, name="Ex3_CLN_Conduit_abc12345", **kw):
        root = make_root(self.runs, name, **kw)
        nb = write_notebook(self.runs, lines_for(root.resolve()))
        self.touch(root, nb)
        return root, nb


@pytest.fixture
def case(ue, tmp_path, monkeypatch):
    c = Case(ue, tmp_path, monkeypatch)
    assert c.snapshot() == 0
    return c


def errors_text(rec):
    return "\n".join(rec["errors"])


# T01: both parts valid, exported evidence
def test_t01_valid_run(case):
    root, nb = case.valid()
    assert case.report() == 0, case.rec["errors"]
    rec = case.rec
    assert rec["selected_root"] == str(root.resolve())
    (att,) = rec["attempts"]
    assert att["selected"] and att["complete"] and att["errors"] == []
    for part, p in zip(PARTS, att["parts"]):
        assert p["part"] == part and p["ws"] == str(root.resolve() / part)
        assert p["termination"] == NORMAL
        assert p["termination_source"] == "stdout of the executed notebook"
        assert p["version_observed"] == {"stdout": "2.8.1", "list": ["2.8.1"]}
        assert p["list"]["normal_termination_lines"] == 0
        assert p["con"]["labels"] == ["CONC"] and p["con"]["precision"] == "single"
        assert p["con"]["first"] == [1, 1, 1, 1.0] and p["con"]["last"] == [
            3,
            1,
            2,
            3.0,
        ]
        assert p["files_missing"] == p["files_extra"] == []
    eng = rec["engine_identity"]
    assert eng["basis"].endswith("not observed per part")
    assert eng["sha256_now"] == eng["sha256"] == case.install_rec["sha256"]
    files = {f["path"]: f for f in rec["exported"]["files"]}
    assert set(files) == {
        "mfusg_transport_tutorial03_Conduit.ipynb",
        "no_dispersion/Ex3_CLN_Conduit.list",
        "dispersion/Ex3_CLN_Conduit.list",
    }
    assert (case.export / "dispersion/Ex3_CLN_Conduit.list").read_bytes() == (
        root / "dispersion/Ex3_CLN_Conduit.list"
    ).read_bytes()
    assert not list(case.export.rglob("*.con"))


# T02: the notebook names a run folder that existed before the snapshot
def test_t02_old_root_named(ue, tmp_path, monkeypatch):
    c = Case(ue, tmp_path, monkeypatch)
    root = make_root(c.runs)
    c.snapshot()
    nb = write_notebook(c.runs, lines_for(root.resolve()))
    c.touch(nb)
    assert c.report() == 1
    assert "not a new run folder" in errors_text(c.rec)
    assert c.rec["ignored_old_roots"] == [str(root.resolve())]
    assert c.rec["selected_root"] is None


# T03: an earlier incomplete attempt is recorded, not counted
def test_t03_retry_incomplete_attempt(case):
    bad = make_root(case.runs, "Ex3_CLN_Conduit_first000", parts=("no_dispersion",))
    case.touch(bad)
    root, _ = case.valid("Ex3_CLN_Conduit_second00")
    assert case.report() == 0, case.rec["errors"]
    att = {a["root"]: a for a in case.rec["attempts"]}
    a = att[str(bad.resolve())]
    assert not a["selected"] and not a["complete"] and a["errors"]
    p = a["parts"][0]
    assert p["stdout_line"] is None and p["termination"] is None
    assert p["version_observed"]["stdout"] is None
    assert p["stdout"].startswith("not available")
    assert att[str(root.resolve())]["selected"]


# T04: a complete earlier attempt with a wrong version counts
def test_t04_retry_complete_invalid(case):
    bad = make_root(case.runs, "Ex3_CLN_Conduit_first000", version="2.7.0")
    case.touch(bad)
    case.valid("Ex3_CLN_Conduit_second00")
    assert case.report() == 1
    assert "versions ['2.7.0']" in errors_text(case.rec)


# T05: the selected run is incomplete; another complete run never replaces it
def test_t05_no_substitution(case):
    other = make_root(case.runs, "Ex3_CLN_Conduit_other000")
    sel = make_root(case.runs, "Ex3_CLN_Conduit_select00", parts=("no_dispersion",))
    nb = write_notebook(case.runs, lines_for(sel.resolve()))
    case.touch(other, sel, nb)
    assert case.report() == 1
    assert case.rec["selected_root"] == str(sel.resolve())
    assert f"{sel.resolve() / 'dispersion'}: not a folder" in errors_text(case.rec)
    att = {a["root"]: a for a in case.rec["attempts"]}
    assert att[str(other.resolve())]["errors"] == []
    assert not att[str(other.resolve())]["selected"]


# T06: notebook missing, older than the snapshot, unchanged since it
@pytest.mark.parametrize("kind", ["missing", "older", "unchanged"])
def test_t06_notebook_provenance(ue, tmp_path, monkeypatch, kind):
    c = Case(ue, tmp_path, monkeypatch)
    root = make_root(c.runs)
    if kind == "unchanged":
        nb = write_notebook(c.runs, lines_for(root.resolve()))
    c.snapshot()
    if kind == "older":
        nb = write_notebook(c.runs, lines_for(root.resolve()))
        c.touch(nb, dt=-10)
    if kind == "unchanged":
        c.touch(nb)
    # the run folder must be new for this test: move the old one aside
    new = make_root(c.runs, "Ex3_CLN_Conduit_new00000")
    c.touch(new)
    assert c.report() == 1
    expected = {
        "missing": "missing (not written or not executed)",
        "older": "older than the snapshot",
        "unchanged": "unchanged since the snapshot",
    }[kind]
    assert expected in errors_text(c.rec)
    assert c.rec["notebook"]["exists"] == (kind != "missing")


# T07: zero or two run folder lines
@pytest.mark.parametrize("n", [0, 2])
def test_t07_root_line_count(case, n):
    root = make_root(case.runs)
    cells = lines_for(root.resolve())
    cells = cells[1:] if n == 0 else [cells[0]] + cells
    nb = write_notebook(case.runs, cells)
    case.touch(root, nb)
    assert case.report() == 1
    assert f"{n} run folder lines, expected 1" in errors_text(case.rec)
    assert case.rec["selected_root"] is None


# T08: part line missing, duplicated, out of order or malformed
@pytest.mark.parametrize("kind", ["missing", "duplicated", "reversed", "malformed"])
def test_t08_part_lines(case, kind):
    root = make_root(case.runs)
    cells = lines_for(root.resolve())
    if kind == "missing":
        cells = cells[:2]
    elif kind == "duplicated":
        cells = cells + [cells[2]]
    elif kind == "reversed":
        cells = [cells[0], cells[2], cells[1]]
    else:
        cells[1] = ["no_dispersion: USG-Transport version 2.8.1 without fields\n"]
    nb = write_notebook(case.runs, cells)
    case.touch(root, nb)
    assert case.report() == 1
    msg = {
        "missing": "dispersion: 0 part lines, expected 1",
        "duplicated": "dispersion: 2 part lines, expected 1",
        "reversed": "part lines out of order",
        "malformed": "malformed part line",
    }[kind]
    assert msg in errors_text(case.rec)


# T09: workspace of a part line in another run folder, part or tutorial
@pytest.mark.parametrize(
    "ws", ["Ex3_CLN_Conduit_other000/dispersion", "{root}/no_dispersion", "Ex6_S/heat"]
)
def test_t09_foreign_workspace(case, ws):
    root = make_root(case.runs)
    cells = lines_for(root.resolve())
    ws = str(case.runs / ws.format(root=root.resolve()))
    cells[2] = [f"dispersion: USG-Transport version 2.8.1; {NORMAL}; {ws}\n"]
    nb = write_notebook(case.runs, cells)
    case.touch(root, nb)
    assert case.report() == 1
    assert "dispersion: workspace" in errors_text(case.rec)


# T10: stdout of the Stallman tutorial in the notebook
@pytest.mark.parametrize("mixed", [False, True])
def test_t10_stallman_stdout(case, mixed):
    root = make_root(case.runs)
    cells = lines_for(root.resolve())
    ex6 = [f"Model output files are written in {case.runs / 'Ex6_Stallman_x'}\n"]
    cells = [ex6] + (cells if mixed else cells[1:])
    nb = write_notebook(case.runs, cells)
    case.touch(root, nb)
    assert case.report() == 1
    msg = "2 run folder lines" if mixed else "is not named Ex3_CLN_Conduit_*"
    assert msg in errors_text(case.rec)


# T11: a normal termination text in the .list does not replace stdout
def test_t11_termination_only_in_list(case):
    root = make_root(case.runs)
    for part in PARTS:
        with open(root / part / "Ex3_CLN_Conduit.list", "a") as f:
            f.write(f" {NORMAL}\n")
    nb = write_notebook(case.runs, lines_for(root.resolve(), term="Stopped"))
    case.touch(root, nb)
    assert case.report() == 1
    assert "no_dispersion: termination 'Stopped'" in errors_text(case.rec)
    p = case.rec["attempts"][0]["parts"][0]
    assert p["list"]["normal_termination_lines"] == 1
    assert p["termination"] == "Stopped"


# T12: versions in stdout and .list
@pytest.mark.parametrize(
    "stdout, listed, msg",
    [
        ("2.7.0", "2.8.1", "version '2.7.0', expected 2.8.1"),
        ("2.8.1", "2.7.0", "versions ['2.7.0']"),
        ("2.8.1", "", "versions []"),
        ("2.8.1", "2.8.1,2.7.0", "versions ['2.7.0', '2.8.1']"),
        ("not reported", "2.8.1", "version 'not reported'"),
    ],
)
def test_t12_versions(case, stdout, listed, msg):
    root = make_root(case.runs, version=listed)
    if not listed:
        for part in PARTS:
            (root / part / "Ex3_CLN_Conduit.list").write_text("no version\n")
    nb = write_notebook(case.runs, lines_for(root.resolve(), version=stdout))
    case.touch(root, nb)
    assert case.report() == 1
    assert msg in errors_text(case.rec)


def recs_drop(recs, i):
    return recs[:i] + recs[i + 1 :]


def recs_set(recs, i, j, v):
    recs[i][j] = v
    return recs


# T13 and T14: concentration files that must be rejected
@pytest.mark.parametrize(
    "kind",
    [
        "missing", "header_truncated", "payload_truncated", "conc01", "conc_tmpr",
        "layer_missing", "layers_swapped", "step_missing", "totim_1e-3",
        "ncol_nrow_swapped", "nlay3", "extra_record",
    ],
)  # fmt: skip
def test_t13_t14_con_rejected(case, kind):
    recs = con_records()
    skip = ()
    if kind == "missing":
        skip = ("Ex3_CLN_Conduit.con",)
    elif kind == "conc01":
        recs = con_records(label="CONC01")
    elif kind == "conc_tmpr":
        recs = recs_set(recs, 3, 4, "TMPR")
    elif kind == "layer_missing":
        recs = recs_drop(recs, 3)
    elif kind == "layers_swapped":
        recs[2], recs[3] = recs[3], recs[2]
    elif kind == "step_missing":
        recs = con_records(calendar=(CALENDAR[0], CALENDAR[2], CALENDAR[2]))
    elif kind == "totim_1e-3":
        recs = recs_set(recs, 5, 3, 3.0 + 1e-3)
    elif kind == "ncol_nrow_swapped":
        recs = con_records(nrow=3, ncol=2)
    elif kind == "nlay3":
        recs = con_records(nlay=3)
    elif kind == "extra_record":
        recs = recs + recs[-1:]
    root = make_root(case.runs, records=recs, skip=skip)
    con = root / "dispersion" / "Ex3_CLN_Conduit.con"
    if kind in ("header_truncated", "payload_truncated"):
        data = con.read_bytes()
        # header: keep 10 bytes of the last record header; payload: drop 4 bytes
        cut = 6 * 4 + 44 - 10 if kind == "header_truncated" else 4
        con.write_bytes(data[: len(data) - cut])
    nb = write_notebook(case.runs, lines_for(root.resolve()))
    case.touch(root, nb)
    assert case.report() == 1
    text = errors_text(case.rec)
    assert "Ex3_CLN_Conduit.con" in text
    p = case.rec["attempts"][0]["parts"][1]
    if kind != "missing":
        # the partial diagnosis keeps the identity of the rejected file
        assert p["con"]["sha256"] and p["con"]["size"] == con.stat().st_size


# T15: double precision file
def test_t15_double_precision(case):
    root, _ = case.valid(real="d")
    assert case.report() == 0, case.rec["errors"]
    assert case.rec["attempts"][0]["parts"][0]["con"]["precision"] == "double"


# T16: non-finite values at the last time step
@pytest.mark.parametrize("value", [float("nan"), float("inf")])
def test_t16_non_finite(case, value):
    recs = con_records()
    recs[-1][8] = [value] * 6
    root = make_root(case.runs, records=recs)
    nb = write_notebook(case.runs, lines_for(root.resolve()))
    case.touch(root, nb)
    assert case.report() == 1
    assert "non-finite values at the last time step" in errors_text(case.rec)


def test_con_structured_dims_checked_before_reading(ue, tmp_path):
    missing = tmp_path / "not_there.con"
    for dims in [(0, 2, 3), (2, -1, 3), (2, 2, 3.0)]:
        with pytest.raises(ValueError, match="invalid dimensions"):
            ue.con_structured(missing, *dims, CALENDAR)
    with pytest.raises(ValueError, match="empty calendar"):
        ue.con_structured(missing, 2, 2, 3, ())


# T17: Windows paths, checked as strings (not on a Windows file system)
def test_t17_windows_path_strings(ue):
    root = r"D:\a\flopy\flopy\flopy\.docs\Notebooks\Ex3_CLN_Conduit_dhw4ciib"
    head = f"USG-Transport version 2.8.1; {NORMAL};"
    lines = [
        (0, 1, f"Model output files are written in {root}"),
        (1, 1, f"no_dispersion: {head} {root}\\no_dispersion"),
        (2, 1, f"dispersion: {head} {root.upper()}\\dispersion"),
    ]
    rec, errors = ue.t03_stdout(lines, ntpath)
    assert errors == [] and rec["root"] == root
    lines[2] = (
        2,
        1,
        f"dispersion: USG-Transport version 2.8.1; {NORMAL}; {root}x\\dispersion",
    )
    assert "dispersion: workspace" in "".join(ue.t03_stdout(lines, ntpath)[1])


# T18: paths with spaces, on the file system
def test_t18_paths_with_spaces(ue, tmp_path, monkeypatch):
    c = Case(ue, tmp_path, monkeypatch, runs_name="Note books")
    c.snapshot()
    root, _ = c.valid("Ex3_CLN_Conduit_a b c d")
    assert c.report() == 0, c.rec["errors"]
    assert c.rec["selected_root"] == str(root.resolve())


# T19: links (out of and within scope) and entries outside the runs folder
LINK_CASES = [
    "root_link",
    "part_link",
    "list_link",
    "part_link_inside",
    "list_link_inside",
    "outside",
    "name",
]


@pytest.mark.parametrize("kind", LINK_CASES)
def test_t19_links_and_scope(case, kind):
    if not hasattr(os, "symlink"):
        pytest.skip("no symbolic links")
    real = make_root(case.tmp / "elsewhere", "Ex3_CLN_Conduit_real0000")
    root = case.runs / "Ex3_CLN_Conduit_abc12345"
    if kind == "root_link":
        root.symlink_to(real, target_is_directory=True)
    elif kind == "part_link":
        make_root(case.runs, parts=("no_dispersion",))
        (root / "dispersion").symlink_to(real / "dispersion", target_is_directory=True)
    elif kind == "list_link":
        make_root(case.runs, skip=("Ex3_CLN_Conduit.list",))
        for part in PARTS:
            (root / part / "Ex3_CLN_Conduit.list").symlink_to(
                real / part / "Ex3_CLN_Conduit.list"
            )
    elif kind == "part_link_inside":
        make_root(case.runs, parts=("no_dispersion",))
        (root / "dispersion").symlink_to(root / "no_dispersion")
    elif kind == "list_link_inside":
        make_root(case.runs)
        lst = root / "dispersion" / "Ex3_CLN_Conduit.list"
        lst.rename(lst.with_name("copy.list"))
        lst.symlink_to(lst.with_name("copy.list"))
    elif kind == "outside":
        root = real
        make_root(case.runs, "Ex3_CLN_Conduit_inside00")
        case.touch(case.runs / "Ex3_CLN_Conduit_inside00")
    else:
        root = make_root(case.runs, "Ex3_other")
    nb = write_notebook(case.runs, lines_for(root))
    case.touch(nb)
    if not root.is_symlink() and root.exists():
        case.touch(root)
    assert case.report() == 1
    text = errors_text(case.rec)
    msg = {
        "root_link": "not a run folder directly in",
        "part_link": "dispersion: not a folder directly in",
        "list_link": "Ex3_CLN_Conduit.list: not a regular file directly in",
        "part_link_inside": "dispersion: not a folder directly in",
        "list_link_inside": "Ex3_CLN_Conduit.list: not a regular file directly in",
        "outside": "not a new run folder directly in",
        "name": "is not named Ex3_CLN_Conduit_*",
    }[kind]
    assert msg in text
    exported = {f["path"] for f in case.rec["exported"]["files"]}
    lists = {p for p in exported if p.endswith(".list")}
    if kind in ("part_link", "part_link_inside"):
        assert lists == {"no_dispersion/Ex3_CLN_Conduit.list"}
    elif kind == "list_link_inside":
        assert lists == {"no_dispersion/Ex3_CLN_Conduit.list"}
    else:
        assert lists == set()


# T20: executable identity inferred from the installation and USGT_EXE
IDENTITY_CASES = [
    "install_errors",
    "env_differs",
    "env_unset",
    "hash_differs",
    "install_missing",
    "exe_link",
]


@pytest.mark.parametrize("kind", IDENTITY_CASES)
def test_t20_engine_identity(case, kind):
    case.valid()
    rec = dict(case.install_rec)
    exe = Path(rec["exe"])
    if kind == "install_errors":
        rec["errors"] = ["sha256 differs"]
    elif kind == "env_differs":
        case.mp.setenv("USGT_EXE", str(exe.with_name("mfusg_gsi")))
    elif kind == "env_unset":
        case.mp.delenv("USGT_EXE")
    elif kind == "hash_differs":
        exe.write_bytes(b"another inert file")
    elif kind == "install_missing":
        case.install = case.tmp / "no-install.json"
    else:
        target = exe.with_name("target")
        exe.rename(target)
        exe.symlink_to(target)
    if kind != "install_missing":
        case.install.write_text(json.dumps(rec))
    assert case.report() == 1
    msg = {
        "install_errors": "install-usgt failed",
        "env_differs": "is not the installed executable",
        "env_unset": "USGT_EXE is not set",
        "hash_differs": "sha256 differs from install-usgt.json",
        "install_missing": "no-install.json",
        "exe_link": "is not a regular file",
    }[kind]
    assert msg in errors_text(case.rec)
    assert case.rec["engine_identity"]["basis"].endswith("not observed per part")
    # the part evidence is kept although the engine check failed
    assert case.rec["attempts"][0]["errors"] == []


# T21: snapshot failed, missing, of another folder or of another command
@pytest.mark.parametrize("kind", ["missing", "errors", "other_runs", "command"])
def test_t21_snapshot(case, kind):
    case.valid()
    path = case.out / "snapshot-t03.json"
    snap = json.loads(path.read_text())
    if kind == "missing":
        path.rename(path.with_name("moved.json"))
    else:
        if kind == "errors":
            snap["errors"] = ["FileNotFoundError: runs"]
        elif kind == "other_runs":
            snap["runs"] = str(case.tmp)
        else:
            snap["command"] = "snapshot"
        path.write_text(json.dumps(snap))
    assert case.report() == 1
    msg = "snapshot-t03.json" if kind == "missing" else "snapshot-t03 failed"
    assert msg in errors_text(case.rec)


def test_snapshot_t03_failure_writes_json(ue, tmp_path):
    out = tmp_path / "out"
    rc = ue.main(
        ["snapshot-t03", "--runs", str(tmp_path / "nowhere"), "--out", str(out)]
    )
    rec = json.loads((out / "snapshot-t03.json").read_text())
    assert rc == 1 and "FileNotFoundError" in rec["errors"][0]


# stdout fragments, several outputs per cell, other streams and cell limits
def test_stdout_fragments_and_cells(ue):
    def out(text, name="stdout"):
        return {"output_type": "stream", "name": name, "text": text}

    first = [
        out(["Model output ", "files are"]),
        out("warning\n", "stderr"),
        out(" written in /r/Ex3_CLN_Conduit_x\n"),
        {"output_type": "display_data", "data": {"text/plain": "no_dispersion: x"}},
    ]
    nb = {
        "cells": [
            {"outputs": first},
            {"outputs": [out("dispersion: ")]},
            {"outputs": [out("USG-Transport version")]},
            {"cell_type": "markdown"},
        ]
    }
    assert ue.notebook_stdout_lines(nb) == [
        (0, 1, "Model output files are written in /r/Ex3_CLN_Conduit_x"),
        (1, 1, "dispersion: "),
        (2, 1, "USG-Transport version"),
    ]


def test_malformed_notebook(case):
    root, nb = case.valid()
    nb.write_text("{not json")
    case.touch(nb)
    assert case.report() == 1
    assert "JSONDecodeError" in errors_text(case.rec)
    # the run folder is still checked and recorded
    assert case.rec["attempts"][0]["complete"]


def test_export_folder_must_be_new(case):
    case.valid()
    case.export.mkdir()
    (case.export / "kept.txt").write_text("earlier staging")
    assert case.report() == 1
    assert "export folder already exists" in errors_text(case.rec)
    assert [p.name for p in case.export.iterdir()] == ["kept.txt"]


def test_failed_run_exports_available_evidence(case):
    root = make_root(case.runs, parts=("no_dispersion",))
    nb = write_notebook(case.runs, lines_for(root.resolve())[:2])
    case.touch(root, nb)
    assert case.report() == 1
    files = {f["path"] for f in case.rec["exported"]["files"]}
    assert files == {
        "mfusg_transport_tutorial03_Conduit.ipynb",
        "no_dispersion/Ex3_CLN_Conduit.list",
    }


def test_calendar_matches_tutorial():
    """T03_CALENDAR and T03_SHAPE match the DIS of the tutorial."""
    spec = importlib.util.spec_from_file_location("usgt_engine_cal", SCRIPT)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    tree = ast.parse(TUTORIAL.read_text())
    consts = {
        t.id: n.value.value
        for n in tree.body
        if isinstance(n, ast.Assign) and isinstance(n.value, ast.Constant)
        for t in n.targets
        if isinstance(t, ast.Name)
    }
    (dis,) = [
        n for n in ast.walk(tree)
        if isinstance(n, ast.Call) and getattr(n.func, "id", None) == "ModflowDis"
    ]  # fmt: skip
    kw = {k.arg: ast.literal_eval(k.value) for k in dis.keywords if k.arg in
          ("nper", "perlen", "nstp", "tsmult")}  # fmt: skip
    assert [a.id for a in dis.args] == ["mf", "nlay", "nrow", "ncol"]
    assert kw == {"nper": 1, "perlen": 3000.0, "nstp": 100}
    assert mod.T03_SHAPE == (consts["nlay"], consts["nrow"], consts["ncol"])
    nstp, perlen = kw["nstp"], kw["perlen"]
    assert mod.T03_CALENDAR == tuple(
        (k, 1, perlen * k / nstp) for k in range(1, nstp + 1)
    )
    assert mod.T03_MODEL in TUTORIAL.read_text()


# T22, T23, T24: the Stallman tutorial commands are unchanged
def write_con_usgt(path, labels, nodes=(2, 3)):
    """USG-T STREAM layout: kstp, kper, pertim, totim, text, nstrt, nend, ilay."""
    with open(path, "wb") as f:
        for label in labels:
            start = 1
            for lay, n in enumerate(nodes, 1):
                text = label.encode().rjust(6).ljust(16)
                head = (1, 1, 1.0, 1.0, text, start, start + n - 1, lay)
                f.write(struct.pack("<2i2f16s3i", *head))
                f.write(struct.pack(f"<{n}f", *[1.0] * n))
                start += n


def make_stallman(runs, name="Ex6_Stallman_abc12345"):
    root = runs / name
    for sim, labels in (("heat", ["TMPR"]), ("solute", ["CONC"]),
                        ("heat_solute", ["CONC01", "TMPR"])):  # fmt: skip
        (root / sim).mkdir(parents=True)
        (root / sim / "Ex6_Stallman.list").write_text(" USG-TRANSPORT VERSION 2.8.1\n")
        write_con_usgt(root / sim / "Ex6_Stallman.con", labels)
    return root


def test_t22_t23_stallman_regression(ue, tmp_path):
    runs, out = tmp_path / "Notebooks", tmp_path / "out"
    runs.mkdir()
    assert ue.main(["snapshot", "--runs", str(runs), "--out", str(out)]) == 0
    snap = json.loads((out / "snapshot.json").read_text())
    assert set(snap) == {"command", "argv", "errors", "python", "python_version",
                         "runs", "time", "roots", "notebook_mtime"}  # fmt: skip
    root = make_stallman(runs)
    make_root(runs)  # a conduit run folder in the same folder (T24)
    write_notebook(
        runs,
        [[f"Model output files are written in {root.resolve()}\n"]],
        name="mfusg_transport_tutorial06_StallmanHeatTransport.ipynb",
    )
    argv = ["report", "--runs", str(runs), "--snapshot", str(out / "snapshot.json")]
    assert ue.main(argv + ["--out", str(out)]) == 0
    rec = json.loads((out / "report.json").read_text())
    assert set(rec) == {"command", "argv", "errors", "python", "python_version",
                        "expected_version", "snapshot", "attempts",
                        "ignored_old_roots", "notebook_roots"}  # fmt: skip
    (att,) = rec["attempts"]
    assert att["root"] == str(root.resolve()) and att["selected"] and att["complete"]
    assert [r["labels"] for r in att["runs"]] == [
        ["TMPR"],
        ["CONC"],
        ["CONC01", "TMPR"],
    ]


def test_t22_stallman_selected_run_must_pass(ue, tmp_path):
    runs, out = tmp_path / "Notebooks", tmp_path / "out"
    runs.mkdir()
    assert ue.main(["snapshot", "--runs", str(runs), "--out", str(out)]) == 0
    root = make_stallman(runs)
    (root / "solute" / "Ex6_Stallman.list").write_text(" USG-TRANSPORT VERSION 2.7.0\n")
    write_notebook(
        runs,
        [[f"Model output files are written in {root.resolve()}\n"]],
        name="mfusg_transport_tutorial06_StallmanHeatTransport.ipynb",
    )
    argv = ["report", "--runs", str(runs), "--snapshot", str(out / "snapshot.json")]
    assert ue.main(argv + ["--out", str(out)]) == 1
    rec = json.loads((out / "report.json").read_text())
    assert any("version 2.7.0, expected 2.8.1" in e for e in rec["errors"])


def test_t24_conduit_report_ignores_stallman(case):
    stallman = make_stallman(case.runs)
    case.touch(stallman)
    root, _ = case.valid()
    assert case.report() == 0, case.rec["errors"]
    assert [a["root"] for a in case.rec["attempts"]] == [str(root.resolve())]


def test_earlier_commands_without_network(ue, tmp_path):
    dest = tmp_path / "dl"
    dest.mkdir()
    (dest / "mfusgt-v2.8.1-linux.zip").write_bytes(b"not the pinned asset")
    out = tmp_path / "out"
    argv = [
        "fetch",
        "--asset",
        "mfusgt-v2.8.1",
        "--dest",
        str(dest),
        "--ostag",
        "linux",
    ]
    assert ue.main(argv + ["--out", str(out)]) == 1
    rec = json.loads((out / "fetch-mfusgt-v2.8.1.json").read_text())
    assert rec["errors"][0].startswith("FileExistsError") and rec["attempts"] == []
    argv = ["install-usgt", "--downloads", str(dest), "--bindir", str(tmp_path / "b")]
    assert ue.main(argv + ["--ostag", "linux", "--out", str(out)]) == 1
    rec = json.loads((out / "install-usgt.json").read_text())
    assert rec["errors"] == [
        f"{dest / 'mfusgt-v2.8.1-linux.zip'}: sha256 differs from the pinned value"
    ]
    assert not (tmp_path / "b").exists()


def test_end_to_end_real_dimensions(tmp_path, monkeypatch):
    """The script run as a program (another Python process, no executable
    of a model), with the tutorial's real grid and calendar."""
    install, rec = make_install(tmp_path, monkeypatch)
    runs, out = tmp_path / "Notebooks", tmp_path / "usgt-engine"
    runs.mkdir()
    env = os.environ | {"USGT_EXE": rec["exe"], "PYTHONDONTWRITEBYTECODE": "1"}

    def run(*args):
        cmd = [sys.executable, "-I", str(SCRIPT), *args, "--out", str(out)]
        return subprocess.run(cmd, env=env, capture_output=True, text=True, check=False)

    assert run("snapshot-t03", "--runs", str(runs)).returncode == 0
    snap = json.loads((out / "snapshot-t03.json").read_text())
    calendar = [(k, 1, 30.0 * k) for k in range(1, 101)]
    root = make_root(runs, records=con_records(2, 100, 100, calendar))
    nb = write_notebook(runs, lines_for(root.resolve()))
    for p in (root, nb):
        os.utime(p, (snap["time"] + 1, snap["time"] + 1))
    args = ["report-t03", "--runs", str(runs), "--install", str(install)]
    args += ["--snapshot", str(out / "snapshot-t03.json")]
    args += ["--export", str(tmp_path / "usgt-t03-evidence")]
    proc = run(*args)
    rep = json.loads((out / "report-t03.json").read_text())
    assert proc.returncode == 0, rep["errors"]
    assert json.loads(proc.stdout) == rep
    con = rep["attempts"][0]["parts"][1]["con"]
    assert (con["size"], con["n_records"], con["last"]) == (
        8008800,
        200,
        [100, 1, 2, 3000.0],
    )
    # one broken part: return code 1 and the JSON still written
    (root / "dispersion" / "Ex3_CLN_Conduit.con").write_bytes(b"\0" * 44)
    args[-1] = str(tmp_path / "usgt-t03-evidence-2")
    proc = run(*args)
    rep = json.loads((out / "report-t03.json").read_text())
    assert proc.returncode == 1 and "8008800" in "\n".join(rep["errors"])
