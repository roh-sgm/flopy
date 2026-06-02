"""
Shared helpers for USG-Transport return-flow list packages (DRT8, QRT8).

Both ``gwf2drt8u.f`` and ``gwf2QRT8u.f`` read variable-length lists of
return-flow recipient nodes with ``U1DINT`` -- one array block per source
cell, written after the main stress-period list -- plus an optional per-cell
integer concentration-change type (``IDCHNGTYP`` / ``IQCHNGTYP``) read at the
end of each data line when the ``CHANGEC`` option is active.

These helpers parse and emit the common ``U1DINT`` control records. Reading
supports ``INTERNAL`` / ``CONSTANT`` (data inline in the package stream) and,
as of Stage 4.5B, ``EXTERNAL <unit>`` and ``OPEN/CLOSE <fname>`` (data in a
separate file, resolved via ``ext_unit_dict`` / ``model.model_ws``), matching
how USG-T's ``U1DINT`` (utl7u1.f) reads the recipient-node arrays in both
``gwf2drt8u.f`` and ``gwf2QRT8u.f``. An ``EXTERNAL`` unit that cannot be
resolved raises an actionable ``NotImplementedError`` rather than being
silently mis-read.

The ``LOCAT>0`` control records carry an ``ICNSTNT FMTIN IPRN`` tail
(``EXTERNAL``/``OPEN/CLOSE`` after their unit/filename). Reading is **free-format
only** (Stage 4.5B polish): ``FMTIN`` must be ``(FREE)`` -- a fixed Fortran
format (e.g. ``(10I8)``) raises an actionable ``NotImplementedError`` instead of
being free-parsed silently. The abbreviated keyword-only forms (``EXTERNAL
<unit>`` / ``OPEN/CLOSE <fname>`` with no tail) are accepted as a read
convenience and treated as ``(FREE)``. The ``ICNSTNT`` multiplier is applied
when non-zero (matching utl7u1.f). Writing always emits an inline
``INTERNAL (FREE)`` block (``Expanded valid write`` -- EXTERNAL/OPEN-CLOSE and
fixed formats are never produced on output).
"""

from __future__ import annotations

import os

from ._usgt_list import _resolve_external_filename, parse_open_close


def _read_ints_free(source, n):
    """Read ``n`` integers free-format from ``source``, spanning lines."""
    vals: list[int] = []
    while len(vals) < n:
        line = source.readline()
        if not line:
            raise ValueError("Unexpected EOF reading U1DINT values.")
        vals.extend(int(float(t)) for t in line.split())
    return vals[:n]


def _parse_u1dint_tail(tail):
    """Parse a U1DINT free-format control tail ``[ICNSTNT [FMTIN [IPRN]]]``.

    ``tail`` is the token list after the keyword (and, for EXTERNAL/OPEN-CLOSE,
    after the unit/filename). Returns ``(icnstnt, fmtin)``. ``fmtin`` is ``None``
    for the abbreviated keyword-only form (no tail) -- FloPy then assumes
    ``(FREE)`` as a read convenience. ``IPRN`` (listing-print flag) is ignored.
    """
    icnstnt = int(tail[0]) if len(tail) >= 1 else 0
    fmtin = tail[1] if len(tail) >= 2 else None
    return icnstnt, fmtin


def _require_free_format(fmtin, package):
    """Reject a non-``(FREE)`` FMTIN: recipient U1DINT is read list-directed only."""
    if fmtin is not None and fmtin.upper() != "(FREE)":
        raise NotImplementedError(
            f"{package}: recipient U1DINT lists are read free-format only; "
            f"FMTIN {fmtin!r} (a fixed Fortran format) is not supported yet. "
            "Write the list with FMTIN (FREE), or inline the nodes "
            "(INTERNAL/CONSTANT)."
        )


def read_u1dint_list(f, n, model=None, ext_unit_dict=None, package="U1DINT"):
    """Read ``n`` integers written as a U1DINT array control block.

    Parameters
    ----------
    f : file-like
        Open text handle positioned at the control-record line.
    n : int
        Number of integers to read.
    model : MfUsg, optional
        Used to resolve relative ``OPEN/CLOSE`` paths and ``EXTERNAL`` units.
    ext_unit_dict : dict, optional
        NAM external-unit dictionary (for ``EXTERNAL``).
    package : str
        Package name, for error messages.

    Returns
    -------
    list of int
        The ``n`` integer values exactly as written (1-based node ids for
        return-flow recipients), with the U1DINT ``ICNSTNT`` multiplier applied
        when it is non-zero (matching utl7u1.f).
    """
    if n <= 0:
        return []
    cntrl = f.readline()
    tokens = cntrl.split()
    if not tokens:
        raise ValueError("Expected a U1DINT control record, got a blank line.")
    kw = tokens[0].upper()

    if kw == "CONSTANT":
        # CONSTANT <icnstnt>: LOCAT=0, every value = ICNSTNT (no list follows).
        return [int(tokens[1])] * n

    # LOCAT>0 forms (INTERNAL / EXTERNAL / OPEN/CLOSE): an ICNSTNT [FMTIN [IPRN]]
    # tail follows the keyword (and the unit/filename). Resolve the source file
    # (None => inline in f) and the tail, then validate FMTIN *before* opening
    # any external file so a bad format never leaks a file descriptor.
    if kw == "INTERNAL":
        fname = None  # values follow inline in f
        icnstnt, fmtin = _parse_u1dint_tail(tokens[1:])
    elif kw == "EXTERNAL":
        unit = int(tokens[1])
        fname = _resolve_external_filename(model, ext_unit_dict, unit)
        if fname is None:
            raise NotImplementedError(
                f"{package}: U1DINT EXTERNAL unit {unit} could not be resolved "
                "(no matching ext_unit_dict entry). Pass ext_unit_dict, or "
                "inline the recipient nodes (INTERNAL/CONSTANT)."
            )
        icnstnt, fmtin = _parse_u1dint_tail(tokens[2:])
    elif kw == "OPEN/CLOSE":
        fname, rest = parse_open_close(cntrl)
        if model is not None and not os.path.isabs(fname):
            fname = os.path.join(model.model_ws, fname)
        icnstnt, fmtin = _parse_u1dint_tail(rest)
    else:
        raise NotImplementedError(
            "Only free-format INTERNAL / CONSTANT / EXTERNAL / OPEN-CLOSE "
            "U1DINT control records are supported for recipient-node lists "
            f"(got: {tokens[0]!r})."
        )

    _require_free_format(fmtin, package)

    if fname is None:
        vals = _read_ints_free(f, n)
    else:
        # Close the opened file even if the read raises (e.g. unexpected EOF).
        source = open(fname, "r")
        try:
            vals = _read_ints_free(source, n)
        finally:
            source.close()

    mult = icnstnt if icnstnt != 0 else 1
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
