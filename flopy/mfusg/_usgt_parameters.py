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

Only the **array-parameter** write path is implemented here (used by ETS for
the ETSR rate array). The analogous **list-parameter** write path (SGB/DRT/QRT)
is documented in ``USGT_STAGE4_04_PARAMETERS.md`` and is intentionally not yet
implemented; those packages still fail explicitly on ``NP*>0``.

Notes
-----
* Surface array packages (ETS/EVT/RCH) read clusters with ``ILFLG=0`` in the
  Fortran, i.e. the cluster line has **no layer field** -- it is
  ``MLTARR ZONARR [zone ...]``. ``ModflowParBc.loadarray`` parses exactly this
  form, and :func:`write_array_parameter_defs` mirrors it.
* The optional per-activation print flag (``IPF`` in ``UPARARRSUB2``) is a
  listing-output control, not data; it is not preserved.
"""


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
