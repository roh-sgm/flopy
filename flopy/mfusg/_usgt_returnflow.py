"""
Shared helpers for USG-Transport return-flow list packages (DRT8, QRT8).

Both ``gwf2drt8u.f`` and ``gwf2QRT8u.f`` read variable-length lists of
return-flow recipient nodes with ``U1DINT`` -- one array block per source
cell, written after the main stress-period list -- plus an optional per-cell
integer concentration-change type (``IDCHNGTYP`` / ``IQCHNGTYP``) read at the
end of each data line when the ``CHANGEC`` option is active.

These helpers parse and emit the common ``U1DINT`` ``INTERNAL`` / ``CONSTANT``
control records. ``EXTERNAL`` and ``OPEN/CLOSE`` recipient lists are not
supported and raise ``NotImplementedError`` rather than being silently
mis-read (USG-T design rule: fail explicitly on unsupported syntax).
"""

from __future__ import annotations


def read_u1dint_list(f, n):
    """Read ``n`` integers written as a U1DINT array control block.

    Parameters
    ----------
    f : file-like
        Open text handle positioned at the control-record line.
    n : int
        Number of integers to read.

    Returns
    -------
    list of int
        The ``n`` integer values exactly as written in the file (1-based node
        ids for return-flow recipients).
    """
    if n <= 0:
        return []
    cntrl = f.readline()
    tokens = cntrl.split()
    if not tokens:
        raise ValueError("Expected a U1DINT control record, got a blank line.")
    kw = tokens[0].upper()
    if kw == "CONSTANT":
        return [int(tokens[1])] * n
    if kw in ("EXTERNAL", "OPEN/CLOSE"):
        raise NotImplementedError(
            f"U1DINT '{kw}' recipient-node lists are not supported; "
            "use INTERNAL or CONSTANT."
        )
    if kw != "INTERNAL":
        raise NotImplementedError(
            "Only free-format INTERNAL/CONSTANT U1DINT control records are "
            f"supported for recipient-node lists (got: {tokens[0]!r})."
        )
    # INTERNAL <icnstnt> <FMTIN> <IPRN>; values follow on subsequent lines.
    icnstnt = int(tokens[1]) if len(tokens) > 1 else 0
    mult = icnstnt if icnstnt != 0 else 1
    vals: list[int] = []
    while len(vals) < n:
        line = f.readline()
        if not line:
            raise ValueError("Unexpected EOF reading U1DINT values.")
        vals.extend(int(float(t)) for t in line.split())
    vals = vals[:n]
    if mult != 1:
        vals = [v * mult for v in vals]
    return vals


def write_u1dint_list(f, values):
    """Write integers as a U1DINT ``INTERNAL (FREE)`` control block.

    ``values`` are written verbatim (already 1-based for node ids). ``IPRN`` is
    set to ``-1`` so USG-T does not echo the list to the listing file.
    """
    f.write("INTERNAL  1  (FREE)  -1\n")
    f.write(" " + " ".join(str(int(v)) for v in values) + "\n")
