"""Install and check the USG-Transport executable used in CI.

fetch: download one pinned release asset (URL, size and sha256 below) for
    this platform into --dest, and fail unless the size and sha256 match.
install-usgt: check the layout of the pinned mfusgt archive, extract only
    the mfusgt executable into an empty --bindir, and write its absolute
    path as USGT_EXE to GITHUB_ENV. Fails on any unexpected archive member,
    or if a file of the same name is in one of the --others directories.
snapshot: record the Ex6_Stallman_* run folders and the executed tutorial
    notebook found in --runs before the notebooks are run.
report: check the run of the Stallman tutorial named in the executed
    notebook, which must be newer than the snapshot: every listing file must
    report exactly the expected version and every .con file must hold
    exactly the expected record labels. Other new run folders (attempts)
    are checked and recorded too.
snapshot-t03: record the Ex3_CLN_Conduit_* run folders and the executed
    conduit tutorial notebook (tutorial03) found in --runs before the
    notebooks are run.
report-t03: check both parts (no_dispersion, dispersion) of the run of the
    conduit tutorial named in its executed notebook, which must be new
    with respect to the snapshot: the stdout lines of the notebook must name
    each part once, in order, with the expected version, normal termination
    and workspace; every listing file must report exactly the expected
    version and every .con file must hold exactly the expected CONC records.
    Other new run folders (attempts) are checked and recorded too. The
    executable identity is inferred from install-usgt.json and USGT_EXE; it
    is not observed per part. The executed notebook and the listing files
    of the selected run are copied to a new --export folder.

Every command writes a JSON record to --out, also when it fails. Only the
standard library is used.
"""

import argparse
import hashlib
import json
import math
import os
import platform
import re
import stat
import struct
import sys
import time
import traceback
import urllib.request
import zipfile
from pathlib import Path

USGT_VERSION = "2.8.1"
VERSION_RE = re.compile(r"USG-TRANSPORT VERSION\s+(\S+)")
LABELS = {"heat": {"TMPR"}, "solute": {"CONC"}, "heat_solute": {"CONC01", "TMPR"}}
NOTEBOOK = "mfusg_transport_tutorial06_StallmanHeatTransport.ipynb"
ROOT_RE = re.compile(r"Model output files are written in (\S+)")

# release assets, from the GitHub releases API (2026-11-10)
_EXE = "https://github.com/MODFLOW-ORG/executables/releases/download/29.0"
_USGT = "https://github.com/MODFLOW-ORG/mfusgt/releases/download/v2.8.1"
PINS = {
    "executables-29.0": {
        "prefix": "modflow_executables-29.0",
        "linux": (
            f"{_EXE}/linux.zip",
            45010855,
            "ab17192a3531dec61d4f4eab3da10f618e85d1f792d5edc71306a3e06e0e22c3",
        ),
        "macarm": (
            f"{_EXE}/macarm.zip",
            23469698,
            "13dc4841f28f47552f1999dde088e8ce1908b97f3b65b937ea603be01ac2dd30",
        ),
        "win64": (
            f"{_EXE}/win64.zip",
            47860985,
            "13d1491da8417df85b1084f03bf987e493120857072543ea6f1dd33c0ed1f81e",
        ),
    },
    "mfusgt-v2.8.1": {
        "prefix": "mfusgt-v2.8.1",
        "linux": (
            f"{_USGT}/linux.zip",
            3463978,
            "6a19614f8b4257bc29efea0de1dcfd912e7fbd35ccee4f8dd4204473909804ed",
        ),
        "macarm": (
            f"{_USGT}/macarm.zip",
            2589976,
            "20f8ae1b1d333effdf46796c7c7b6d96b4877aae6cacb1409f0c1d0ea69ae0b5",
        ),
        "win64": (
            f"{_USGT}/win64.zip",
            3896957,
            "90cfb5983e6dc5e670ddfd71da036276181dd12e4690823a0e86c7b1a88f40a4",
        ),
    },
}


def get_ostag():
    """Same tags as flopy/utils/get_modflow.py get_ostag()."""
    if sys.platform.startswith("linux"):
        return "linux"
    if sys.platform.startswith("win"):
        return "win" + ("64" if sys.maxsize > 2**32 else "32")
    if sys.platform.startswith("darwin"):
        return "mac" + ("arm" if platform.processor() == "arm" else "")
    raise ValueError(f"platform {sys.platform!r} not supported")


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for block in iter(lambda: f.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def asset_path(dest, asset, ostag):
    return Path(dest) / f"{PINS[asset]['prefix']}-{ostag}.zip"


def fetch(args, record):
    ostag = args.ostag or get_ostag()
    if ostag not in PINS[args.asset]:
        raise ValueError(f"no pinned {args.asset} asset for ostag {ostag!r}")
    url, size, digest = PINS[args.asset][ostag]
    dst = asset_path(args.dest, args.asset, ostag)
    part = dst.with_name(dst.name + ".part")
    record.update(asset=args.asset, ostag=ostag, url=url, path=str(dst))
    record.update(expected_size=size, expected_sha256=digest, attempts=[])
    if dst.exists() or part.exists():
        raise FileExistsError(f"{dst} or {part} already exists")
    for attempt in range(1, args.tries + 1):
        try:
            with urllib.request.urlopen(url, timeout=120) as r, open(part, "xb") as f:
                record["attempts"].append({"attempt": attempt, "final_url": r.url})
                for block in iter(lambda: r.read(1 << 20), b""):
                    f.write(block)
            break
        except OSError as err:
            record["attempts"].append({"attempt": attempt, "error": repr(err)})
            if part.exists():
                part.rename(part.with_name(f"{part.name}.failed{attempt}"))
            if attempt == args.tries:
                raise
            time.sleep(5 * attempt)
    record.update(size=part.stat().st_size, sha256=sha256(part))
    if record["size"] != size or record["sha256"] != digest:
        record["errors"].append(f"{url}: size or sha256 differ from the pinned values")
        return
    part.rename(dst)


def check_usgt_zip(zf, ostag):
    """Members of the mfusgt archive as recorded, the name of the executable
    member, and the layout errors. The release workflow of mfusgt v2.8.1
    zips the files of bin/ without folders (zip -j, Compress-Archive bin\\*):
    mfusgt and mfusgtdbl, with .exe on Windows. Anything else is an error."""
    ext = ".exe" if ostag.startswith("win") else ""
    exe, allowed = f"mfusgt{ext}", {f"mfusgt{ext}", f"mfusgtdbl{ext}"}
    members, errors, seen = [], [], set()
    for info in zf.infolist():
        mode = info.external_attr >> 16
        members.append(
            {
                "name": info.filename,
                "size": info.file_size,
                "crc": info.CRC,
                "mode": oct(mode),
                "flag_bits": info.flag_bits,
            }
        )
        name = info.filename
        if name.casefold() in seen:
            errors.append(f"duplicate member {name!r}")
        seen.add(name.casefold())
        if info.is_dir() or "/" in name or "\\" in name or ":" in name:
            errors.append(f"member {name!r} is not a plain file name")
        elif stat.S_ISLNK(mode):
            errors.append(f"member {name!r} is a symbolic link")
        elif info.flag_bits & 0x1:
            errors.append(f"member {name!r} is encrypted")
        elif name not in allowed:
            errors.append(f"unexpected member {name!r}; expected {sorted(allowed)}")
    if exe not in {m["name"] for m in members}:
        errors.append(f"no member named {exe!r}")
    return members, exe, errors


def locate_exe(bindir, ostag):
    """Absolute path of the mfusgt executable in bindir, without searching
    PATH or the current directory: it must be the only entry of bindir
    whose name starts with mfusgt (any case), have the exact expected name,
    be a regular file (not a link) and resolve inside bindir."""
    bindir = Path(bindir).resolve(strict=True)
    name = "mfusgt.exe" if ostag.startswith("win") else "mfusgt"
    errors = []
    variants = sorted(
        p.name for p in bindir.iterdir() if p.name.casefold().startswith("mfusgt")
    )
    if variants != [name]:
        errors.append(f"{bindir}: expected only {name!r}, found {variants}")
        return None, errors
    exe = bindir / name
    st = os.lstat(exe)
    if not stat.S_ISREG(st.st_mode):
        errors.append(f"{exe} is not a regular file")
    elif exe.resolve(strict=True).parent != bindir:
        errors.append(f"{exe} resolves outside {bindir}")
    elif not ostag.startswith("win") and not os.access(exe, os.X_OK):
        errors.append(f"{exe} is not executable")
    return (None if errors else str(exe)), errors


def install_usgt(args, record):
    ostag = args.ostag or get_ostag()
    zpath = asset_path(args.downloads, "mfusgt-v2.8.1", ostag)
    bindir = Path(args.bindir).expanduser()
    record.update(ostag=ostag, zip=str(zpath), bindir=str(bindir))
    record["zip_sha256"] = sha256(zpath)
    if record["zip_sha256"] != PINS["mfusgt-v2.8.1"][ostag][2]:
        record["errors"].append(f"{zpath}: sha256 differs from the pinned value")
        return
    with zipfile.ZipFile(zpath) as zf:
        record["members"], name, errors = check_usgt_zip(zf, ostag)
        record["errors"] += errors
        if errors:
            return
        bindir.mkdir(parents=True, exist_ok=True)
        if any(bindir.iterdir()):
            record["errors"].append(f"{bindir} is not empty")
            return
        with zf.open(name) as src, open(bindir / name, "xb") as dst:
            for block in iter(lambda: src.read(1 << 20), b""):
                dst.write(block)
    if not ostag.startswith("win"):
        os.chmod(bindir / name, 0o755)
    exe, errors = locate_exe(bindir, ostag)
    record["errors"] += errors
    if exe is None:
        return
    record.update(exe=exe, size=os.path.getsize(exe), sha256=sha256(exe))
    for other in args.others:
        other = Path(other).expanduser()
        if other.is_dir() and (other / name).exists():
            record["errors"].append(f"{other} also holds {name}")
    record["on_path"] = {
        n: [
            str(Path(d) / n)
            for d in os.environ.get("PATH", "").split(os.pathsep)
            if d and (Path(d) / n).exists()
        ]
        for n in (name, "mfusg_gsi" + Path(name).suffix)
    }
    if not record["errors"] and "GITHUB_ENV" in os.environ:
        with open(os.environ["GITHUB_ENV"], "a") as f:
            f.write(f"USGT_EXE={exe}\n")
        record["USGT_EXE"] = exe


def con_labels(path):
    """Labels of a USG-T STREAM concentration file, reading it as single or
    double precision; raises ValueError if neither reads it to the end."""
    data = Path(path).read_bytes()
    for real, size in (("f", 4), ("d", 8)):
        head = struct.Struct(f"<2i2{real}16s3i")
        pos, labels = 0, set()
        try:
            while pos < len(data):
                kstp, kper, _, _, text, nstrt, nend, _ = head.unpack_from(data, pos)
                n = nend - nstrt + 1
                if kstp < 1 or kper < 1 or n < 1 or not text.isascii():
                    raise ValueError
                labels.add(text.decode().strip())
                pos += head.size + n * size
            if pos == len(data):
                return sorted(labels)
        except (struct.error, ValueError, UnicodeDecodeError):
            pass
    raise ValueError(f"{path}: not a single or double precision STREAM file")


def list_versions(path):
    """(line number, version) of every 'USG-TRANSPORT VERSION <v>' in path.
    The version is written from one string in the program (in the source,
    in the listing header and on the BAS line), at no fixed line."""
    found = []
    with open(path, errors="replace") as f:
        for i, line in enumerate(f, 1):
            found += [(i, m.group(1)) for m in VERSION_RE.finditer(line)]
    return found


def check_root(root):
    """Check the three simulations of one run folder; return the run
    records, the errors and whether every simulation wrote its files."""
    runs, errors, complete = [], [], True
    for name, expected in LABELS.items():
        ws = root / name
        run = {"ws": str(ws)}
        runs.append(run)
        lists = sorted(ws.glob("*.list")) if ws.is_dir() else []
        cons = sorted(ws.glob("*.con")) if ws.is_dir() else []
        if len(lists) != 1 or len(cons) != 1:
            errors.append(f"{ws}: {len(lists)} .list, {len(cons)} .con")
            complete = False
            continue
        try:
            found = list_versions(lists[0])
            run["versions"] = found[:10] + (["..."] if len(found) > 10 else [])
            versions = {v for _, v in found}
            if not versions:
                errors.append(f"{lists[0]}: no USG-TRANSPORT VERSION line")
            elif len(versions) > 1:
                errors.append(f"{lists[0]}: ambiguous versions {sorted(versions)}")
            elif versions != {USGT_VERSION}:
                errors.append(
                    f"{lists[0]}: version {versions.pop()}, expected {USGT_VERSION}"
                )
            run["con_sha256"] = sha256(cons[0])
            run["labels"] = con_labels(cons[0])
        except (OSError, ValueError) as err:
            errors.append(f"{ws}: {err!r}")
            continue
        if set(run["labels"]) != expected:
            errors.append(
                f"{cons[0]}: labels {run['labels']}, expected {sorted(expected)}"
            )
    return runs, errors, complete


def run_roots(runs):
    return sorted(p for p in Path(runs).glob("Ex6_Stallman_*") if p.is_dir())


def snapshot(args, record):
    nb = Path(args.runs) / NOTEBOOK
    record.update(
        runs=str(Path(args.runs).resolve()),
        time=time.time(),
        roots=[str(p.resolve()) for p in run_roots(args.runs)],
        notebook_mtime=nb.stat().st_mtime if nb.exists() else None,
    )


def report(args, record):
    snap = json.loads(Path(args.snapshot).read_text())
    runs = Path(args.runs).resolve()
    record.update(expected_version=USGT_VERSION, snapshot=snap, attempts=[])
    if snap["errors"] or Path(snap["runs"]) != runs:
        record["errors"].append("snapshot failed or was taken in another folder")
    old = set(snap["roots"])
    new = [p.resolve() for p in run_roots(runs) if str(p.resolve()) not in old]
    record["ignored_old_roots"] = sorted(old)
    # the executed notebook names the folder of its own (last, passing) run;
    # jupytext 1.19.6 executes it in the folder of the .ipynb it writes
    nb = runs / NOTEBOOK
    selected = None
    if not nb.exists() or nb.stat().st_mtime < snap["time"]:
        record["errors"].append(f"{nb}: missing or older than the snapshot")
    else:
        cells = json.loads(nb.read_text(encoding="utf-8"))["cells"]
        text = "".join(
            "".join(o.get("text", "")) for c in cells for o in c.get("outputs", [])
        )
        named = ROOT_RE.findall(text)
        record["notebook_roots"] = named
        if len(named) != 1:
            record["errors"].append(f"{nb}: {len(named)} run folders named")
        else:
            selected = (runs / named[0]).resolve()
            if selected not in new:
                record["errors"].append(f"{selected}: not a new run folder")
    if not new:
        record["errors"].append(f"no new Ex6_Stallman_* run folder in {runs}")
    for root in new:
        runs_, errors, complete = check_root(root)
        record["attempts"].append(
            {
                "root": str(root),
                "selected": root == selected,
                "complete": complete,
                "mtime": root.stat().st_mtime,
                "runs": runs_,
                "errors": errors,
            }
        )
        if root.stat().st_mtime < snap["time"]:
            record["errors"].append(f"{root}: older than the snapshot")
        if root == selected or complete:
            record["errors"] += errors


# conduit tutorial (mfusg_transport_tutorial03_Conduit.py)
T03_NOTEBOOK = "mfusg_transport_tutorial03_Conduit.ipynb"
T03_PREFIX = "Ex3_CLN_Conduit_"
T03_MODEL = "Ex3_CLN_Conduit"
T03_PARTS = ("no_dispersion", "dispersion")
T03_EXTS = (
    "bas", "bcf", "bct", "chd", "cln", "dis", "nam", "oc", "sms", "wel",
    "list", "con", "hds", "cbc", "BCT.cbc", "clncb", "clncn", "clnhd",
)  # fmt: skip
T03_FILES = tuple(f"{T03_MODEL}.{ext}" for ext in T03_EXTS)
# (nlay, nrow, ncol) and (kstp, kper, totim) of every time step of the DIS
# of the tutorial: nper=1, perlen=3000.0, nstp=100, tsmult 1 (the default)
T03_SHAPE = (2, 100, 100)
T03_CALENDAR = tuple((kstp, 1, 30.0 * kstp) for kstp in range(1, 101))
T03_ROOT_LINE = re.compile(r"Model output files are written in (.+)")
T03_PART_LINE = re.compile(
    r"(no_dispersion|dispersion): USG-Transport version (.+?); (.+?); (.+)"
)
T03_NORMAL = "normal termination of simulation"
T03_BASIS = (
    "inferred from install-usgt.json and USGT_EXE of the report step; "
    "not observed per part"
)


def norm_path(path, pathmod=os.path):
    """Path string normalized for comparison on the platform of pathmod."""
    return pathmod.normcase(pathmod.normpath(path))


def notebook_stdout_lines(nb):
    """(cell, line number, text) of every stdout line of a notebook: the
    stdout stream outputs of one cell are joined in order (their text may
    be one string or a list of fragments) before splitting into lines; a
    line never continues into another cell."""
    lines = []
    for c, cell in enumerate(nb["cells"]):
        text = ""
        for out in cell.get("outputs", []):
            if out.get("output_type") == "stream" and out.get("name") == "stdout":
                part = out.get("text", "")
                text += part if isinstance(part, str) else "".join(part)
        lines += [(c, i, line) for i, line in enumerate(text.splitlines(), 1)]
    return lines


def t03_stdout(lines, pathmod=os.path):
    """Run folder and part lines named in the stdout lines of the executed
    conduit tutorial, and the errors: exactly one run folder line and one
    line per part, in order, with the expected version, normal termination
    and the workspace <run folder>/<part>."""
    errors = []
    roots = [(c, i, m[1]) for c, i, s in lines if (m := T03_ROOT_LINE.fullmatch(s))]
    found = {part: [] for part in T03_PARTS}
    for c, i, s in lines:
        if not s.startswith(tuple(f"{part}:" for part in T03_PARTS)):
            continue
        m = T03_PART_LINE.fullmatch(s)
        if m is None:
            errors.append(f"cell {c} line {i}: malformed part line {s!r}")
            continue
        found[m[1]].append(
            {"part": m[1], "version": m[2], "termination": m[3], "ws": m[4]}
            | {"cell": c, "line": i, "text": s}
        )
    record = {"root_lines": [r[2] for r in roots], "root": None, "parts": {}}
    if len(roots) != 1:
        errors.append(f"{len(roots)} run folder lines, expected 1")
    else:
        root = roots[0][2]
        if not pathmod.basename(pathmod.normpath(root)).startswith(T03_PREFIX):
            errors.append(f"run folder {root!r} is not named {T03_PREFIX}*")
        else:
            record["root"] = root
    for part in T03_PARTS:
        record["parts"][part] = found[part]
        if len(found[part]) != 1:
            errors.append(f"{part}: {len(found[part])} part lines, expected 1")
            continue
        line = found[part][0]
        if line["version"] != USGT_VERSION:
            errors.append(
                f"{part}: version {line['version']!r}, expected {USGT_VERSION}"
            )
        if line["termination"].lower() != T03_NORMAL:
            errors.append(f"{part}: termination {line['termination']!r}")
        if record["root"] is not None:
            expected = norm_path(pathmod.join(record["root"], part), pathmod)
            if norm_path(line["ws"], pathmod) != expected:
                errors.append(
                    f"{part}: workspace {line['ws']!r} is not in the run folder"
                )
            if (line["cell"], line["line"]) < roots[0][:2]:
                errors.append(f"{part}: part line before the run folder line")
    nd, d = (found[part] for part in T03_PARTS)
    if len(nd) == len(d) == 1 and (nd[0]["cell"], nd[0]["line"]) > (
        d[0]["cell"],
        d[0]["line"],
    ):
        errors.append("part lines out of order")
    return record, errors


def con_structured(path, nlay, nrow, ncol, calendar):
    """Check a structured concentration file (MODFLOW head file layout:
    kstp, kper, pertim, totim, text, ncol, nrow, ilay, then nrow * ncol
    reals): one CONC record per layer and time step of calendar, in order,
    times within 4 units of the stored precision, nothing else in the file
    and finite values at the last time step. The size is checked against
    the expected records before any record is read. Raises ValueError."""
    dims = (nlay, nrow, ncol)
    if not all(type(n) is int and n > 0 for n in dims) or not calendar:
        raise ValueError(f"invalid dimensions {dims} or empty calendar")
    size = os.path.getsize(path)
    nrec = len(calendar) * nlay
    layouts = []
    for precision, real in (("single", "f"), ("double", "d")):
        head = struct.Struct(f"<2i2{real}16s3i")
        layouts.append(
            (precision, real, head, head.size + nrow * ncol * struct.calcsize(real))
        )
    sizes = [nrec * layout[3] for layout in layouts]
    if size not in sizes:
        raise ValueError(
            f"{path}: {size} bytes; expected {sizes[0]} (single) or {sizes[1]} "
            f"(double) for {len(calendar)} time steps x {nlay} layers"
        )
    precision, real, head, rec_size = layouts[sizes.index(size)]
    eps = 2.0**-23 if real == "f" else 2.0**-52
    expected = [(k, p, lay, t) for k, p, t in calendar for lay in range(1, nlay + 1)]
    digest, first, last = hashlib.sha256(), None, None
    with open(path, "rb") as f:
        for i, (kstp, kper, ilay, totim) in enumerate(expected):
            block = f.read(rec_size)
            if len(block) != rec_size:
                raise ValueError(f"{path}: record {i + 1} is truncated")
            digest.update(block)
            k1, k2, _, t, text, nc, nr, il = head.unpack_from(block)
            try:
                label = text.decode("ascii").strip()
            except UnicodeDecodeError:
                label = repr(text)
            got = (k1, k2, il, nc, nr)
            if label != "CONC" or got != (kstp, kper, ilay, ncol, nrow):
                raise ValueError(
                    f"{path}: record {i + 1} has label {label!r} and (kstp, "
                    f"kper, ilay, ncol, nrow) {got}; expected 'CONC' and "
                    f"{(kstp, kper, ilay, ncol, nrow)}"
                )
            if not math.isclose(t, totim, rel_tol=4 * eps, abs_tol=0.0):
                raise ValueError(f"{path}: record {i + 1} totim {t}, expected {totim}")
            if i >= nrec - nlay:
                values = struct.unpack_from(f"<{nrow * ncol}{real}", block, head.size)
                if not all(map(math.isfinite, values)):
                    raise ValueError(f"{path}: non-finite values at the last time step")
            first = first or [k1, k2, il, t]
            last = [k1, k2, il, t]
        if f.read(1):
            raise ValueError(f"{path}: data after the last expected record")
    return {
        "path": str(path),
        "sha256": digest.hexdigest(),
        "size": size,
        "precision": precision,
        "n_records": nrec,
        "labels": ["CONC"],
        "nlay": nlay,
        "nrow": nrow,
        "ncol": ncol,
        "n_steps": len(calendar),
        "first": first,
        "last": last,
        "final_finite": True,
    }


def inside(path, parent):
    """Whether path is a real entry (not a link) whose resolved location is
    directly inside the resolved parent."""
    path = Path(path)
    return not path.is_symlink() and path.resolve().parent == Path(parent).resolve()


def t03_roots(runs):
    """Run folders of the conduit tutorial in runs (resolved), and the
    entries named like one that are not real folders directly in runs."""
    roots, rejected = [], []
    for p in sorted(Path(runs).glob(f"{T03_PREFIX}*")):
        if p.is_dir() and inside(p, runs):
            roots.append(p.resolve())
        else:
            rejected.append(str(p))
    return roots, rejected


def t03_file(ws, name):
    """ws/name if it is a regular file (not a link) directly inside ws."""
    p = Path(ws) / name
    try:
        mode = os.lstat(p).st_mode
    except FileNotFoundError:
        return None
    return p if stat.S_ISREG(mode) and inside(p, ws) else None


def check_root_t03(root, stdout=None):
    """Check both parts of one run folder of the conduit tutorial; return
    the part records, the errors and whether every part wrote its .list and
    .con files. stdout holds the part lines of the executed notebook for
    the selected run folder; other run folders have no attributable stdout,
    so their version and termination are recorded as not available."""
    parts, errors, complete = [], [], True
    nlay, nrow, ncol = T03_SHAPE
    for part in T03_PARTS:
        ws = Path(root) / part
        rec = {"part": part, "ws": str(ws), "stdout_line": None}
        rec.update(termination=None, termination_source=None, errors=[])
        rec["version_observed"] = {"stdout": None, "list": None}
        rec.update(files_missing=[], files_extra=[], files_rejected=[])
        rec.update(list=None, con=None)
        parts.append(rec)
        err = rec["errors"]
        if stdout is None or part not in stdout:
            rec["stdout"] = "not available: no stdout of this part of this run folder"
        else:
            line = stdout[part]
            rec.update(stdout_line=line["text"], termination=line["termination"])
            rec["termination_source"] = "stdout of the executed notebook"
            rec["version_observed"]["stdout"] = line["version"]
        if not (ws.is_dir() and inside(ws, root)):
            err.append(f"{ws}: not a folder directly in {root}")
            complete = False
            errors += err
            continue
        names = {p.name for p in ws.iterdir()}
        rec["files_missing"] = sorted(set(T03_FILES) - names)
        rec["files_extra"] = sorted(names - set(T03_FILES))
        if rec["files_missing"]:
            err.append(f"{ws}: missing {rec['files_missing']}")
        files = {}
        for name in sorted(set(T03_FILES) & names):
            files[name] = t03_file(ws, name)
            if files[name] is None:
                rec["files_rejected"].append(name)
                err.append(f"{ws / name}: not a regular file directly in {ws}")
        lst, con = files.get(f"{T03_MODEL}.list"), files.get(f"{T03_MODEL}.con")
        if lst is None or con is None:
            complete = False
        if lst is not None:
            try:
                found = list_versions(lst)
                text = lst.read_text(errors="replace").lower()
            except OSError as e:
                err.append(f"{lst}: {e!r}")
            else:
                rec["list"] = {
                    "path": str(lst),
                    "sha256": sha256(lst),
                    "size": lst.stat().st_size,
                    "versions": found[:10] + (["..."] if len(found) > 10 else []),
                    "normal_termination_lines": text.count(T03_NORMAL),
                }
                versions = sorted({v for _, v in found})
                rec["version_observed"]["list"] = versions
                if versions != [USGT_VERSION]:
                    err.append(f"{lst}: versions {versions}, expected {USGT_VERSION}")
        if con is not None:
            try:
                rec["con"] = con_structured(con, nlay, nrow, ncol, T03_CALENDAR)
            except (OSError, ValueError) as e:
                rec["con"] = {"path": str(con), "sha256": sha256(con)}
                rec["con"]["size"] = con.stat().st_size
                err.append(str(e))
        errors += err
    return parts, errors, complete


def engine_identity(install, errors):
    """Executable identity inferred from install-usgt.json (keys exe,
    sha256, zip_sha256, errors, USGT_EXE) and USGT_EXE of this process; the
    sha256 of the installed file is computed again without running it.
    This does not show which executable ran each part."""
    rec = {"basis": T03_BASIS, "install": str(install)}
    try:
        inst = json.loads(Path(install).read_text())
    except (OSError, ValueError) as e:
        errors.append(f"{install}: {e!r}")
        return rec
    rec.update({k: inst.get(k) for k in ("exe", "sha256", "zip_sha256")})
    rec.update(install_errors=inst.get("errors"), USGT_EXE_install=inst.get("USGT_EXE"))
    rec["USGT_EXE_env"] = env = os.environ.get("USGT_EXE")
    exe = inst.get("exe")
    if inst.get("command") != "install-usgt" or inst.get("errors") or not exe:
        errors.append(f"{install}: install-usgt failed or names no executable")
        return rec
    if "USGT_EXE" in inst and norm_path(inst["USGT_EXE"]) != norm_path(exe):
        errors.append(f"{install}: USGT_EXE {inst['USGT_EXE']!r} is not exe {exe!r}")
    if not env:
        errors.append("USGT_EXE is not set in the report step")
    elif norm_path(env) != norm_path(exe):
        errors.append(f"USGT_EXE {env!r} is not the installed executable {exe!r}")
    try:
        mode = os.lstat(exe).st_mode
        if not stat.S_ISREG(mode):
            errors.append(f"{exe} is not a regular file")
        else:
            rec["sha256_now"] = sha256(exe)
            if rec["sha256_now"] != inst.get("sha256"):
                errors.append(f"{exe}: sha256 differs from install-usgt.json")
    except OSError as e:
        errors.append(f"{exe}: {e!r}")
    return rec


def export_t03(dest, notebook, selected, errors):
    """Copy the executed notebook and the listing files of the selected run
    folder (only files already checked to be in scope) to the new folder
    dest; never write into an existing folder."""
    dest = Path(dest)
    rec = {"dir": str(dest), "files": []}
    if os.path.lexists(dest):
        errors.append(f"{dest}: export folder already exists")
        return rec
    dest.mkdir(parents=True)
    sources = [(notebook, T03_NOTEBOOK)] if notebook is not None else []
    for part in T03_PARTS if selected is not None else ():
        ws = Path(selected) / part
        lst = t03_file(ws, f"{T03_MODEL}.list") if inside(ws, selected) else None
        if lst is not None:
            sources.append((lst, f"{part}/{T03_MODEL}.list"))
    for src, rel in sources:
        dst = dest / rel
        dst.parent.mkdir(exist_ok=True)
        with open(src, "rb") as a, open(dst, "xb") as b:
            for block in iter(lambda: a.read(1 << 20), b""):
                b.write(block)
        entry = {"source": str(src), "path": rel, "sha256": sha256(dst)}
        entry["size"] = dst.stat().st_size
        rec["files"].append(entry)
        if entry["sha256"] != sha256(src):
            errors.append(f"{dst}: copy differs from {src}")
    return rec


def snapshot_t03(args, record):
    t = time.time()
    runs = Path(args.runs).resolve(strict=True)
    roots, rejected = t03_roots(runs)
    nb = runs / T03_NOTEBOOK
    notebook = None
    if os.path.lexists(nb):
        st = os.lstat(nb)
        notebook = {"path": str(nb), "mtime": st.st_mtime, "sha256": None}
        if stat.S_ISREG(st.st_mode):
            notebook["sha256"] = sha256(nb)
    record.update(runs=str(runs), time=t, roots=[str(p) for p in roots])
    record.update(rejected_roots=rejected, notebook=notebook)


def report_t03(args, record):
    errors = record["errors"]
    runs = Path(args.runs).resolve()
    record.update(expected_version=USGT_VERSION, snapshot=None, notebook=None)
    record.update(engine_identity=None, ignored_old_roots=[], rejected_roots=[])
    record.update(selected_root=None, attempts=[], exported=None)
    snap, old, snap_time = {}, set(), None
    try:
        snap = json.loads(Path(args.snapshot).read_text())
        record["snapshot"] = snap
        old, snap_time = set(snap["roots"]), snap["time"]
        if (
            snap.get("command") != "snapshot-t03"
            or snap["errors"]
            or Path(snap["runs"]) != runs
        ):
            errors.append("snapshot-t03 failed or was taken in another folder")
    except (OSError, ValueError, KeyError, TypeError) as e:
        errors.append(f"snapshot {args.snapshot}: {e!r}")
    record["ignored_old_roots"] = sorted(old)
    roots, rejected = t03_roots(runs)
    record["rejected_roots"] = rejected
    errors += [f"{p}: not a run folder directly in {runs}" for p in rejected]
    new = [p for p in roots if str(p) not in old]
    if not new:
        errors.append(f"no new {T03_PREFIX}* run folder in {runs}")

    # the executed notebook names the run folder and the workspace of each
    # part of its own run; jupytext writes it next to the tutorial
    nb = runs / T03_NOTEBOOK
    notebook = record["notebook"] = {"path": str(nb), "exists": os.path.lexists(nb)}
    stdout, nb_ok = None, None
    if not notebook["exists"]:
        errors.append(f"{nb}: missing (not written or not executed)")
    elif t03_file(runs, T03_NOTEBOOK) is None:
        errors.append(f"{nb}: not a regular file directly in {runs}")
    else:
        nb_ok = nb
        st = nb.stat()
        notebook.update(mtime=st.st_mtime, size=st.st_size, sha256=sha256(nb))
        notebook["newer_than_snapshot"] = (
            snap_time is not None and st.st_mtime >= snap_time
        )
        if not notebook["newer_than_snapshot"]:
            errors.append(f"{nb}: older than the snapshot")
        if (snap.get("notebook") or {}).get("sha256") == notebook["sha256"]:
            errors.append(f"{nb}: unchanged since the snapshot")
        try:
            lines = notebook_stdout_lines(json.loads(nb.read_text(encoding="utf-8")))
        except (OSError, ValueError, KeyError, TypeError, AttributeError) as e:
            errors.append(f"{nb}: {e!r}")
        else:
            parsed, errs = t03_stdout(lines)
            notebook.update(parsed)
            errors += [f"{nb}: {e}" for e in errs]
            stdout = {p: v[0] for p, v in parsed["parts"].items() if len(v) == 1}

    selected = None
    named = notebook.get("root")
    if named is not None:
        path = Path(named) if Path(named).is_absolute() else runs / named
        if path.is_dir() and inside(path, runs) and path.resolve() in new:
            selected = path.resolve()
        else:
            errors.append(f"{named}: not a new run folder directly in {runs}")
    record["selected_root"] = None if selected is None else str(selected)
    for root in new:
        try:
            parts, errs, complete = check_root_t03(
                root, stdout if root == selected else None
            )
        except OSError as e:
            parts, errs, complete = [], [f"{root}: {e!r}"], False
        record["attempts"].append(
            {
                "root": str(root),
                "selected": root == selected,
                "complete": complete,
                "mtime": root.stat().st_mtime,
                "parts": parts,
                "errors": errs,
            }
        )
        if snap_time is None or root.stat().st_mtime < snap_time:
            errors.append(f"{root}: not newer than the snapshot")
        if root == selected or complete:
            errors += errs
    record["engine_identity"] = engine_identity(args.install, errors)
    if args.export:
        record["exported"] = export_t03(args.export, nb_ok, selected, errors)


COMMANDS = {
    "fetch": fetch,
    "install-usgt": install_usgt,
    "snapshot": snapshot,
    "report": report,
    "snapshot-t03": snapshot_t03,
    "report-t03": report_t03,
}


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    sub = parser.add_subparsers(dest="command", required=True)
    p = sub.add_parser("fetch")
    p.add_argument("--asset", required=True, choices=sorted(PINS))
    p.add_argument("--dest", required=True)
    p.add_argument("--ostag")
    p.add_argument("--tries", type=int, default=3)
    p = sub.add_parser("install-usgt")
    p.add_argument("--downloads", required=True)
    p.add_argument("--bindir", required=True)
    p.add_argument("--others", nargs="*", default=[])
    p.add_argument("--ostag")
    for name in ("snapshot", "report"):
        p = sub.add_parser(name)
        p.add_argument("--runs", required=True)
        if name == "report":
            p.add_argument("--snapshot", required=True)
    p = sub.add_parser("snapshot-t03")
    p.add_argument("--runs", required=True)
    p = sub.add_parser("report-t03")
    p.add_argument("--runs", required=True)
    p.add_argument("--snapshot", required=True)
    p.add_argument("--install", required=True)
    p.add_argument("--export")
    for p in sub.choices.values():
        p.add_argument("--out", required=True)
    args = parser.parse_args(argv)
    record = {
        "command": args.command,
        "argv": sys.argv[1:],
        "errors": [],
        "python": sys.executable,
        "python_version": sys.version.split()[0],
    }
    try:
        COMMANDS[args.command](args, record)
    except Exception as err:
        record["errors"].append(f"{type(err).__name__}: {err}")
        record["traceback"] = traceback.format_exc()
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    tag = f"-{args.asset}" if args.command == "fetch" else ""
    (out / f"{args.command}{tag}.json").write_text(json.dumps(record, indent=2) + "\n")
    print(json.dumps(record, indent=2))
    return 1 if record["errors"] else 0


if __name__ == "__main__":
    sys.exit(main())
