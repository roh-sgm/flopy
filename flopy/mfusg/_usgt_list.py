"""
Shared reader for MODFLOW list-control records at the start of a USG-T
stress-period list (the main lists of SGB, QRT, and DRT).

The USG-T Fortran list readers (`ULSTRDU`, `SGWF2QRT8LR`, `SGWF2DRT8LR`)
accept, before the data rows, the standard MODFLOW list controls::

    EXTERNAL <unit>      read the rows from another NAM-declared unit's file
    OPEN/CLOSE <fname>   read the rows from a named file
    SFAC <value>         scale the package's primary data column

`begin_list_block` consumes any of these and hands back the source the rows
should be read from, the scale factor, and the first already-read data line, so
the caller's row loop is identical whether or not controls were present.

Writes emit expanded inline rows (``Expanded valid write``): the in-memory
values are already final (SFAC applied on read), so no control records are
produced on output.
"""

import os


def _open_close_filename(line):
    """Extract the filename from an ``OPEN/CLOSE`` control line.

    Supports an unquoted token or a single/double-quoted name (which may itself
    contain spaces), mirroring the MODFLOW ``URWORD`` quoted-word convention and
    FloPy's quote-stripping for ``OPEN/CLOSE`` array records.
    """
    parts = line.split(None, 1)  # ["OPEN/CLOSE", "<remainder>"]
    rest = parts[1].strip() if len(parts) > 1 else ""
    if rest[:1] in ("'", '"'):
        quote = rest[0]
        end = rest.find(quote, 1)
        if end != -1:
            return rest[1:end]
        return rest[1:].strip()  # unmatched quote: best effort
    tokens = rest.split()
    return tokens[0] if tokens else ""


def _resolve_external_filename(model, ext_unit_dict, unit):
    """Best-effort filename for an EXTERNAL unit via ext_unit_dict."""
    if ext_unit_dict and unit in ext_unit_dict:
        fname = getattr(ext_unit_dict[unit], "filename", None)
        if fname:
            if model is not None and not os.path.isabs(fname):
                return os.path.join(model.model_ws, fname)
            return fname
    return None


def begin_list_block(f, model=None, ext_unit_dict=None, package="list"):
    """Consume leading SFAC / EXTERNAL / OPEN-CLOSE controls before list rows.

    Parameters
    ----------
    f : file-like
        Handle positioned at the line after the stress-period header.
    model : MfUsg, optional
        Used to resolve relative ``OPEN/CLOSE`` paths and ``EXTERNAL`` units.
    ext_unit_dict : dict, optional
        NAM external-unit dictionary (for ``EXTERNAL``).
    package : str
        Package name, for error messages.

    Returns
    -------
    (source, sfac, first_line, to_close)
        ``source`` is the handle the *remaining* rows (and, for QRT/DRT, the
        recipient ``U1DINT`` blocks) should be read from. ``sfac`` is the scale
        factor (1.0 if none). ``first_line`` is the first data row, already read
        from ``source``. ``to_close`` is a handle to close when finished
        (for ``EXTERNAL``/``OPEN/CLOSE``) or ``None``.
    """
    sfac = 1.0
    source = f
    to_close = None

    line = f.readline()
    tok = line.split()
    kw = tok[0].upper() if tok else ""

    if kw == "EXTERNAL":
        unit = int(tok[1])
        fname = _resolve_external_filename(model, ext_unit_dict, unit)
        if fname is None:
            raise NotImplementedError(
                f"{package}: EXTERNAL unit {unit} could not be resolved "
                "(no ext_unit_dict entry); EXTERNAL list input is unsupported "
                "without it."
            )
        source = open(fname, "r")
        to_close = source
        line = source.readline()
        tok = line.split()
        kw = tok[0].upper() if tok else ""
    elif kw == "OPEN/CLOSE":
        fname = _open_close_filename(line)
        if model is not None and not os.path.isabs(fname):
            fname = os.path.join(model.model_ws, fname)
        source = open(fname, "r")
        to_close = source
        line = source.readline()
        tok = line.split()
        kw = tok[0].upper() if tok else ""

    if kw == "SFAC":
        sfac = float(tok[1])
        line = source.readline()

    return source, sfac, line, to_close
