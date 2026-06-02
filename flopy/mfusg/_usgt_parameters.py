"""Shared helpers to preserve MODFLOW-style named parameters in USG-T packages.

USG-T inherits MODFLOW's array- and list-parameter machinery (``parutl7.f``:
``UPARARRRP`` / ``UPARARRSUB2`` for *array* parameters such as ETS/EVT/RCH;
``UPARLSTRP`` / ``UPARLSTSUB`` for *list* parameters such as SGB/DRT/QRT).

FloPy's :class:`flopy.modflow.ModflowParBc` already *parses* both forms into a
structured ``bc_parms`` dictionary on load. Historically USG-T packages then
expanded the parameters to plain arrays and wrote ``NP*=0`` ("Expanded valid
write"). These helpers add the missing *write* side for the array-parameter
form, so a loaded parameterized package can be written back with its parameter
syntax intact (load -> write -> reload) instead of being silently expanded.

Array-parameter ``bc_parms`` layout (from ``ModflowParBc.loadarray``)::

    bc_parms[name] = [
        {"partyp": str, "parval": str, "nclu": int, "timevarying": bool},
        {instance_name: [[mltarr, zonarr, [izone, ...]], ...nclu clusters]},
    ]

Current per-package status (all reuse the helpers below):

* **ETS** — array parameters (ETSR rate) are *preserved* (Stage 4.4A); the
  array-parameter write path lives here.
* **HFB/DRT** — list parameters are *preserved* and *authorable from scratch*
  (Stage 4.4B/4.4D + 4.6C-A/4.6B), including per-stress-period activations
  (DRT) / a global active list (HFB).
* **QRT** — list parameters are *structurally* preserved (Stage 4.4E); active
  parameters round-trip but are not execution-guaranteed (Fortran scales the
  wrong field).
* **SGB** — *definition* preservation + from-scratch *definition authoring* only
  (Stage 4.4C / 4.6C-B); **active** SGB parameters are unsupported because USG-T
  2.7 reads them with ``PARTYP='SGB'`` but activates them as ``PTYP='G'`` (a type
  conflict that aborts the run), so any activation raises ``NotImplementedError``.

The shared name/``parval`` validators (``check_parameter_name`` /
``check_parval``) and the list/array read/write helpers below are reused by the
from-scratch authoring paths; parameter ``INSTANCES`` remain unsupported
everywhere.

Notes
-----
* Surface array packages (ETS/EVT/RCH) read clusters with ``ILFLG=0`` in the
  Fortran, i.e. the cluster line has **no layer field** -- it is
  ``MLTARR ZONARR [zone ...]``. ``ModflowParBc.loadarray`` parses exactly this
  form, and :func:`write_array_parameter_defs` mirrors it.
* The optional per-activation print flag (``IPF`` in ``UPARARRSUB2``) is a
  listing-output control, not data; it is not preserved.
"""

import numbers


def check_parameter_name(name, what, prefix):
    """Validate a MODFLOW parameter name (definition key or activation name).

    ``UPARLSTRP`` (``parutl7.f``) reads ``PARNAM`` as a single ``URWORD`` word
    into a ``CHARACTER*10`` buffer and upper-cases it, so the name must be a
    non-empty, whitespace-free token of at most 10 characters. Raises
    ``ValueError`` otherwise. ``what`` describes the field and ``prefix`` is the
    caller's error prefix (e.g. ``"MfUsgHfb.write_file"``), so packages get a
    consistent but self-identifying message. Shared by the from-scratch
    parameter-authoring paths (DRT/HFB; QRT/SGB to follow).
    """
    if not isinstance(name, str):
        raise ValueError(
            f"{prefix}: {what} must be a string, got {type(name).__name__}."
        )
    if not name or any(c.isspace() for c in name):
        raise ValueError(
            f"{prefix}: {what} {name!r} must be a single non-empty token with no "
            "whitespace (the Fortran reads PARNAM as one word)."
        )
    if len(name) > 10:
        raise ValueError(
            f"{prefix}: {what} {name!r} exceeds the 10-character Fortran PARNAM "
            "limit (CHARACTER*10), which would truncate or collide with another "
            "name."
        )


def check_parval(parval, name, prefix):
    """Validate a parameter value: a real number (``int``/``float``/``np.integer``/
    ``np.floating`` -- any ``numbers.Real``, excluding ``bool``), or a single
    whitespace-free string token. ``UPARLSTRP`` reads ``PARVAL`` as one numeric
    value via ``URWORD``, so a missing/blank value or a multi-token string (e.g.
    ``"1 2"``) raises ``ValueError``. Returns ``parval`` unchanged."""
    if parval is None:
        raise ValueError(f"{prefix}: parameter '{name}' is missing 'parval'.")
    if isinstance(parval, bool) or not isinstance(parval, (numbers.Real, str)):
        raise ValueError(
            f"{prefix}: parameter '{name}' parval must be a number or a "
            f"single-token string; got {parval!r}."
        )
    if isinstance(parval, str):
        if not parval.strip():
            raise ValueError(f"{prefix}: parameter '{name}' is missing 'parval'.")
        if any(c.isspace() for c in parval.strip()):
            raise ValueError(
                f"{prefix}: parameter '{name}' parval {parval!r} must be a single "
                "token (the Fortran reads one PARVAL value)."
            )
    return parval


def _build_array_clusters(name, clusters, prefix):
    """Validate + canonicalize an array-parameter instance's clusters.

    Each cluster is ``(MLTARR, ZONARR[, zones])`` -- ``ZONARR='ALL'`` means the
    whole grid (an empty zone list). Returns the
    ``[[mltarr, zonarr, [izone, ...]], ...]`` form ``ModflowParBc`` uses. Raises
    ``ValueError`` for an empty cluster list or non-positive-integer zones.
    """
    if not clusters:
        raise ValueError(
            f"{prefix}: parameter '{name}' needs at least one cluster "
            "(MLTARR ZONARR [zones...])."
        )
    out = []
    for clu in clusters:
        if len(clu) < 2:
            raise ValueError(
                f"{prefix}: parameter '{name}' cluster {clu!r} needs at least "
                "MLTARR and ZONARR."
            )
        mltarr, zonarr = str(clu[0]), str(clu[1])
        zones = list(clu[2]) if len(clu) > 2 and clu[2] is not None else []
        if zonarr.lower() == "all":
            zones = []
        else:
            izones = []
            for z in zones:
                if not float(z).is_integer() or int(z) <= 0:
                    raise ValueError(
                        f"{prefix}: parameter '{name}' cluster zone {z!r} must be "
                        "a positive integer (or use ZONARR='ALL' with no zones)."
                    )
                izones.append(int(z))
            zones = izones
        out.append([mltarr, zonarr, zones])
    return out


def build_array_parameter_bc_parms(parameters, partyp, prefix):
    """Build a validated ``ModflowParBc.bc_parms`` dict from an ergonomic
    from-scratch ``parameters`` mapping (UPARARRRP array parameters, e.g. ETS).

    ``parameters`` is ``{name: {"parval": ..., "clusters": [...]}}`` for a static
    parameter, or ``{name: {"parval": ..., "instances": {inst: [clusters]}}}`` for
    a time-varying (``INSTANCES``) one. ``partyp`` defaults to / is validated
    against the expected type (e.g. ``"ets"``). ``nclu`` is computed from the
    clusters (and validated if given). Names/instances obey the Fortran ``PARNAM``
    limit (single token, <=10 chars, unique case-insensitively); ``parval`` is a
    number or single token. Raises ``ValueError`` before any file is written.

    Returns the ``bc_parms`` dict (keys lower-cased) ready for
    ``ModflowParBc(bc_parms)``.
    """
    bc_parms = {}
    for name, pdef in parameters.items():
        check_parameter_name(name, "parameter name", prefix)
        if not isinstance(pdef, dict):
            raise ValueError(
                f"{prefix}: parameter '{name}' must be a dict, got "
                f"{type(pdef).__name__}."
            )
        ptyp = pdef.get("partyp", partyp)
        if str(ptyp).lower() != partyp.lower():
            raise ValueError(
                f"{prefix}: parameter '{name}' partyp must be '{partyp}'; got {ptyp!r}."
            )
        parval = check_parval(pdef.get("parval"), name, prefix)
        if "instances" in pdef and pdef["instances"] is not None:
            timevarying = True
            insts = pdef["instances"]
            if not insts:
                raise ValueError(
                    f"{prefix}: parameter '{name}' has an empty 'instances' map."
                )
            pinst = {}
            nclu = None
            for instnam, clusters in insts.items():
                check_parameter_name(instnam, f"instance name of '{name}'", prefix)
                bcinst = _build_array_clusters(name, clusters, prefix)
                if nclu is None:
                    nclu = len(bcinst)
                elif len(bcinst) != nclu:
                    raise ValueError(
                        f"{prefix}: parameter '{name}' instances must all have the "
                        f"same number of clusters (nclu={nclu})."
                    )
                if instnam.lower() in pinst:
                    raise ValueError(
                        f"{prefix}: parameter '{name}' has duplicate instance name "
                        f"'{instnam}' (case-insensitive)."
                    )
                pinst[instnam.lower()] = bcinst
        else:
            timevarying = False
            bcinst = _build_array_clusters(name, pdef.get("clusters"), prefix)
            nclu = len(bcinst)
            pinst = {"static": bcinst}
        if "nclu" in pdef and pdef["nclu"] != nclu:
            raise ValueError(
                f"{prefix}: parameter '{name}' declares nclu={pdef['nclu']} but "
                f"carries {nclu} clusters."
            )
        if name.lower() in bc_parms:
            raise ValueError(
                f"{prefix}: duplicate parameter definition name '{name}' "
                "(case-insensitive); the Fortran upper-cases PARNAM."
            )
        bc_parms[name.lower()] = [
            {
                "partyp": str(ptyp).lower(),
                "parval": parval,
                "nclu": nclu,
                "timevarying": timevarying,
            },
            pinst,
        ]
    return bc_parms


def write_array_parameter_defs(f, pak_parms):
    """Write array-parameter definitions (``UPARARRRP`` grammar) to ``f``.

    Mirrors :meth:`flopy.modflow.ModflowParBc.loadarray` so the block
    round-trips byte-for-structure.

    Parameters
    ----------
    f : file handle
        Open text handle positioned where the definition block belongs (after
        item 2a / the optional ESFACTOR record for ETS).
    pak_parms : flopy.modflow.ModflowParBc
        Parsed array parameters (``pak_parms.bc_parms``).
    """
    for name, (pdict, pinst) in pak_parms.bc_parms.items():
        nclu = pdict["nclu"]
        timevarying = pdict["timevarying"]
        line = f"{name} {pdict['partyp']} {pdict['parval']} {nclu}"
        if timevarying:
            line += f" INSTANCES {len(pinst)}"
        f.write(line + "\n")
        for instnam, clusters in pinst.items():
            if timevarying:
                f.write(f"{instnam}\n")
            for mltarr, zonarr, izones in clusters:
                parts = [str(mltarr), str(zonarr)]
                parts.extend(str(int(iz)) for iz in izones)
                f.write(" ".join(parts) + "\n")


def read_active_array_parameters(f, count, pak_parms):
    """Read ``count`` active-parameter records (``UPARARRSUB2`` grammar).

    Each record is ``PNAME [INSTANCE] [IPF]``. The instance token is only
    consumed for parameters that were defined as time-varying (``INSTANCES``);
    a trailing print flag is ignored.

    Parameters
    ----------
    f : file handle
    count : int
        Number of active-parameter lines to read (the per-period ``INETSR``).
    pak_parms : flopy.modflow.ModflowParBc
        Parsed parameter definitions, used to tell whether a named parameter is
        time-varying.

    Returns
    -------
    list of (name, instance_or_None)
        File order is preserved.
    """
    records = []
    for _ in range(count):
        t = f.readline().split("#")[0].strip().split()
        name = t[0]
        pdef = pak_parms.bc_parms.get(name.lower())
        timevarying = bool(pdef[0]["timevarying"]) if pdef else False
        instance = t[1] if (timevarying and len(t) > 1) else None
        records.append((name, instance))
    return records


def write_active_array_parameters(f, records):
    """Write active-parameter records, one ``PNAME [INSTANCE]`` per line."""
    for name, instance in records:
        if instance and str(instance).lower() != "static":
            f.write(f"{name} {instance}\n")
        else:
            f.write(f"{name}\n")


# ---------------------------------------------------------------------------
# List parameters (UPARLSTAL / UPARLSTRP / UPARLSTSUB) -- HFB (4.4B), SGB (4.4C).
# ---------------------------------------------------------------------------
#
# A list parameter owns NLST list rows; its value scales a designated column on
# expansion. The definition header (UPARLSTRP) is
# ``PARNAM PARTYP PARVAL NLST [INSTANCES n]``; the NLST rows that follow are
# package-specific (e.g. HFB barrier rows ``LAYER IROW1 ICOL1 IROW2 ICOL2
# FACTOR`` / ``NODE1 NODE2 FACTOR``; SGB rows ``NODE GRADIENT [aux ...]``), so
# the row reader/writer stays in the package. Activation (UPARLSTSUB) reads one
# parameter name per active parameter (plus an instance name when the parameter
# is time-varying). The helpers below carry parameter names only; callers that
# do not support INSTANCES reject NUMINST>0 themselves.
#
# Packages whose header carries the parameter count inline (SGB/DRT/QRT, via
# ``UPARLSTAL``) use ``read_list_parameter_count`` / ``write_list_parameter_count``
# for the leading ``PARAMETER NP MXL`` record. HFB has no such line (its NPHFB is
# in item 1).


def read_list_parameter_count(line):
    """Parse a leading ``PARAMETER NP MXL`` record (``UPARLSTAL`` grammar).

    Returns ``(np, mxl)``; ``(0, 0)`` when the line is not a PARAMETER record
    (so the same line is then parsed as the package header by the caller).
    """
    t = line.split("#")[0].strip().split()
    if t and t[0].upper() == "PARAMETER":
        np_ = int(t[1])
        mxl = int(t[2]) if len(t) > 2 else 0
        return np_, mxl
    return 0, 0


def write_list_parameter_count(f, np_, mxl):
    """Write a leading ``PARAMETER NP MXL`` record."""
    f.write(f"PARAMETER {np_} {mxl}\n")


def read_list_parameter_header(line):
    """Parse a list-parameter definition header (``UPARLSTRP`` grammar).

    ``PARNAM PARTYP PARVAL NLST [INSTANCES n]``.

    Returns
    -------
    (name, partyp, parval, nlst, numinst)
        ``parval`` is kept as the original string (value text preserved);
        ``numinst`` is 0 when there is no ``INSTANCES`` token.
    """
    t = line.split("#")[0].strip().split()
    name = t[0]
    partyp = t[1]
    parval = t[2]
    nlst = int(t[3])
    numinst = 0
    if len(t) > 4 and t[4].upper() == "INSTANCES":
        numinst = int(t[5])
    return name, partyp, parval, nlst, numinst


def write_list_parameter_header(f, name, partyp, parval, nlst):
    """Write a list-parameter definition header (no ``INSTANCES``)."""
    f.write(f"{name} {partyp} {parval} {nlst}\n")


def read_active_list_parameters(f, count):
    """Read ``count`` active list-parameter names (``UPARLSTSUB`` activation).

    Returns a list of names in file order (HFB activations carry only the name).
    """
    names = []
    for _ in range(count):
        t = f.readline().split("#")[0].strip().split()
        names.append(t[0])
    return names


def write_active_list_parameters(f, names):
    """Write active list-parameter names, one per line."""
    for name in names:
        f.write(f"{name}\n")


def resolve_list_parameter(parameters, name):
    """Resolve a list-parameter definition by name, case-insensitively.

    The Fortran upper-cases both the defined and the activated parameter name
    (``UPARLSTSUB`` / ``SGWF2DRT8LS`` / ``SGWF2QRT8LS``), so an activation can
    differ from its definition only in case. Returns the matching definition
    (``parameters[key]``) or ``None``.
    """
    if not parameters:
        return None
    if name in parameters:
        return parameters[name]
    lname = name.lower()
    for key, pdef in parameters.items():
        if key.lower() == lname:
            return pdef
    return None
