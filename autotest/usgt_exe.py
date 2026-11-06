"""USG-Transport executable for the optional run tests in test_usg_transport.py
and test_usg_transport_exe.py.

``USGT_EXE`` (a name on PATH or a path) is used when set, without falling back
to another executable. Otherwise the first of ``mfusgt`` (MODFLOW executables
30.0 and later) and ``mfusg_gsi`` (earlier releases) found on PATH is used. The
version of the selected executable is not inferred from its name.
"""

import os
import shutil

import pytest

USGT_EXE_NAMES = ("mfusgt", "mfusg_gsi")


def select_usgt_exe(environ=None, which=shutil.which):
    """Return (absolute path, None), or (None, reason) when no executable is found."""
    environ = os.environ if environ is None else environ
    if "USGT_EXE" in environ:
        override = environ["USGT_EXE"]
        path = which(override) if override else None
        if path is None:
            return None, f"USGT_EXE={override!r} is not an executable (no fallback)"
        return os.path.abspath(path), None
    for name in USGT_EXE_NAMES:
        path = which(name)
        if path is not None:
            return os.path.abspath(path), None
    return None, (
        "missing executable: none of " + ", ".join(USGT_EXE_NAMES) + " on PATH "
        "and USGT_EXE is not set"
    )


USGT_EXE, USGT_EXE_REASON = select_usgt_exe()
requires_usgt = pytest.mark.skipif(USGT_EXE is None, reason=USGT_EXE_REASON or "")
