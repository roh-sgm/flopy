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

Every command writes a JSON record to --out, also when it fails. Only the
standard library is used.
"""

import argparse
import hashlib
import json
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


COMMANDS = {
    "fetch": fetch,
    "install-usgt": install_usgt,
    "snapshot": snapshot,
    "report": report,
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
