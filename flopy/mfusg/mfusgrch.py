"""
mfusgrch module.  Contains the MfUsgRch class. Note that the user can access
the MfUsgRch class as `flopy.mfusg.MfUsgRch`.

"""

import copy
import io
import warnings

import numpy as np

from ..modflow.mfparbc import ModflowParBc as mfparbc
from ..modflow.mfrch import ModflowRch
from ..utils import Transient2d, Transient3d, Util2d, read1d
from ..utils.flopy_io import line_parse
from ..utils.utils_def import (
    get_pak_vals_shape,
    get_unitnumber_from_ext_unit_dict,
)
from .mfusg import MfUsg


def _is_int(value):
    return isinstance(value, (int, np.integer)) and not isinstance(value, bool)


def _words(line):
    """Upper-case words of a whole line: URWORD splits at blanks, commas and
    tabs and reads words in any case (utl7u1.f :89-96), comments included."""
    return line.replace(",", " ").replace("\t", " ").upper().split()


def _period_flags(value, nper, name):
    """INRCHZONES/INSELEV values as {kper: int}; a list or array has one
    value per stress period."""
    if value is None:
        return {}
    if isinstance(value, dict):
        flags = dict(value)
    elif isinstance(value, (list, tuple, np.ndarray)):
        if len(value) != nper:
            raise ValueError(
                f"{name} has {len(value)} values for {nper} stress periods"
            )
        flags = dict(enumerate(value))
    else:
        raise ValueError(f"{name} must be a dict {{kper: int}} or one int per period")
    for kper, flag in flags.items():
        if not _is_int(kper) or not 0 <= kper < nper:
            raise ValueError(f"{name}: stress period {kper!r} not in 0..{nper - 1}")
        if not _is_int(flag):
            raise ValueError(f"{name}[{kper}] = {flag!r} is not an integer")
    return {int(k): int(v) for k, v in flags.items()}


def _array_keys(value, nper, name):
    """Check the stress-period keys of selev/iznrch as given: a dict or a
    Transient2d needs integer keys in 0..nper-1 (no float, bool or text keys,
    which Transient2d would convert); other values are period 1."""
    if isinstance(value, Transient2d):
        keys = value.transient_2ds
    elif isinstance(value, dict):
        keys = value
    else:
        return
    for kper in keys:
        if not _is_int(kper) or not 0 <= kper < nper:
            raise ValueError(
                f"{name}: stress period key {kper!r} is not an integer in 0..{nper - 1}"
            )


def _as_transient2d(model, value, dtype, name):
    """value as a new Transient2d; a given Transient2d is not changed
    (Transient2d(value=Transient2d) rebuilds the Util2d of its dict)."""
    if value is None:
        return None
    if isinstance(value, Transient2d):
        value = copy.copy(value)
        value.transient_2ds = dict(value.transient_2ds)
        shape = value.shape
    else:
        shape = get_pak_vals_shape(model, value)
    return Transient2d(model, shape, dtype, value, name=name)


def _rts_seepelev_options(nper, mxrtzones, seepelev, iznrch, selev, iniznrch, inselev):
    """Checks that need no array conversion; returns the INRCHZONES and
    INSELEV values as {kper: int}."""
    if not _is_int(mxrtzones) or mxrtzones < 0:
        raise ValueError(f"mxrtzones = {mxrtzones!r} must be an integer >= 0")
    if not _is_int(seepelev) or seepelev not in (0, 1):
        raise ValueError(f"seepelev = {seepelev!r} must be 0 or 1")
    izn = _period_flags(iniznrch, nper, "iniznrch")
    sel = _period_flags(inselev, nper, "inselev")
    if mxrtzones == 0 and (iznrch is not None or izn):
        raise ValueError("iznrch/iniznrch need the RTS option (mxrtzones > 0)")
    if not seepelev and (selev is not None or sel):
        raise ValueError("selev/inselev need the SEEPELEV option (seepelev=1)")
    _array_keys(iznrch, nper, "iznrch")
    _array_keys(selev, nper, "selev")
    return izn, sel


def _rts_seepelev_plan(
    nper, shape, mxrtzones, seepelev, iznrch, selev, iniznrch, inselev
):
    """Check the RTS zones (INRCHZONES, IZNRCH) and seepage elevations
    (INSELEV, SELEV) and return the [INRCHZONES, INSELEV] values of every
    stress period, None for a keyword that is not written.

    USG-T RP (gwf2rch8u1.f): INRCHZONES n > 0 reads IZNRCH and sets the
    number of RTS zones to n (the RTS header count is not used then); the
    RTS array is allocated with the period 1 count, so period 1 must read
    zones and later positive values must keep that count. n <= 0 or no
    keyword reuses the zones. INSELEV n > 0 reads SELEV and turns seepage
    on; 0 or no keyword turns it off; n < 0 keeps the previous state and
    SELEV. An array is read only with a value > 0, so a period with an
    array needs a value > 0 and the other way round. Without explicit
    values: a period with an array writes INRCHZONES with the period 1
    count (mxrtzones in period 1) / INSELEV 1, any other period writes no
    INRCHZONES / INSELEV -1."""
    izn, sel = _rts_seepelev_options(
        nper, mxrtzones, seepelev, iznrch, selev, iniznrch, inselev
    )
    defined = {}
    for name, t2d in (("iznrch", iznrch), ("selev", selev)):
        defined[name] = set()
        if t2d is None:
            continue
        for kper, u2d in t2d.transient_2ds.items():
            if tuple(u2d.shape) != tuple(shape):
                raise ValueError(
                    f"{name} period {kper + 1} has shape {u2d.shape}, RECH {shape}"
                )
            defined[name].add(kper)
    plan = []
    count = None  # number of RTS zones set in period 1
    for kper in range(nper):
        n = s = None
        if mxrtzones > 0:
            new = kper in defined["iznrch"]
            default = (mxrtzones if count is None else count) if new else None
            n = izn.get(kper, default)
            if n is not None and n > 0:
                if not new:
                    raise ValueError(
                        f"INRCHZONES {n} in period {kper + 1} needs iznrch data"
                    )
                if count is None:
                    count = n
                elif n != count:
                    raise ValueError(
                        f"INRCHZONES {n} in period {kper + 1}: the number of RTS "
                        f"zones is fixed in period 1 to {count}"
                    )
            elif new:
                raise ValueError(
                    f"iznrch is given for period {kper + 1} but INRCHZONES {n} "
                    "does not read it"
                )
            elif kper == 0:
                raise ValueError("RTS needs iznrch and INRCHZONES > 0 in period 1")
        if seepelev and (selev is not None or sel):
            new = kper in defined["selev"]
            s = sel.get(kper, 1 if new else -1)
            if s > 0 and not new:
                raise ValueError(f"INSELEV {s} in period {kper + 1} needs selev data")
            if s <= 0 and new:
                raise ValueError(
                    f"selev is given for period {kper + 1} but INSELEV {s} "
                    "does not read it"
                )
        plan.append([n, s])
    return plan


def _irchconc_flags(irchconc, mcompt):
    """IRCHCONC as MCOMPT values 0 or 1 (default: 1 for every component).

    USG-T RP reads an array for a flag not 0, but BCT applies a component
    only when its flag is 1 and AR sizes RCHCONC with the sum of the flags
    (gwf2rch8u1.f :116-120, gwt2bndsu1.f :3513-3533), so other values are
    rejected."""
    if irchconc is None:
        return [1] * mcompt
    if not isinstance(irchconc, (list, tuple, np.ndarray)):
        raise ValueError(
            f"irchconc = {irchconc!r} must be a list of MCOMPT = {mcompt} values 0 or 1"
        )
    flags = list(irchconc)
    if len(flags) != mcompt:
        raise ValueError(
            f"irchconc has {len(flags)} values; IRCHCONC needs MCOMPT = {mcompt}"
        )
    for icomp, flag in enumerate(flags):
        if not _is_int(flag) or flag not in (0, 1):
            raise ValueError(
                f"irchconc[{icomp}] = {flag!r}: IRCHCONC must be 0 or 1 (USG-T "
                "uses the recharge concentration of a component only when its "
                "flag is 1)"
            )
    return [int(flag) for flag in flags]


def _conc_flags(iconc, mcompt, irchconc, rchconc):
    """Check the CONC data (IRCHCONC, RCHCONC as a Transient3d) and return
    the IRCHCONC flags, None without the CONC option."""
    if not iconc:
        if irchconc is not None or rchconc is not None:
            raise ValueError("irchconc/rchconc need the CONC option (iconc=1)")
        return None
    flags = _irchconc_flags(irchconc, mcompt)
    if rchconc is not None and rchconc.shape[0] != mcompt:
        raise ValueError(
            f"rchconc has {rchconc.shape[0]} components; RCHCONC needs MCOMPT = "
            f"{mcompt}"
        )
    if 1 in flags and (rchconc is None or 0 not in rchconc.transient_3ds):
        raise ValueError(
            "IRCHCONC reads recharge concentrations but rchconc has no data for "
            "stress period 1: USG-T allocates RCHCONC without initial values "
            "(gwf2rch8u1.f :120) and would use it unread"
        )
    return flags


def _conc_mcompt(model, who):
    """MCOMPT for the CONC option: the BCT counts (model._transport_mcompt),
    or counts set by hand on a model without BCT and itrnsp 0 when MCOMPT
    is at least 1."""
    if (
        model._transport_context()[0] == "none"
        and _is_int(model.mcomp)
        and _is_int(model.iheat)
        and model.mcompt >= 1
    ):
        return model.mcompt
    return model._transport_mcompt(who)


# with RTS, RP applies the current RTS record to every IRCH cell when it
# reads RECH, and load keeps RECH as data of every period (as upstream)
_RTS_RECH_NOTE = (
    "RECH is written as data in every period that has it (load keeps it in "
    "every period), so a file that reused RECH (INRECH < 0) in those periods "
    "gets another RTS recharge when written again: USG-T applies the current "
    "RTS record to every IRCH cell when it reads RECH (gwf2rch8u1.f :429-442)."
)


def _irch_reuse_warning(reuse, conc, rts):
    """Warn about the stress periods ([kper, INIRCH]) that reuse IRCH with
    CONC and/or RTS."""
    periods = ", ".join(f"{k + 1} (INIRCH {n})" for k, n in reuse)
    opts = []
    why = []
    if conc:
        opts.append("CONC")
        why.append(
            "With CONC, USG-T adds no recharge mass in a period that reuses "
            "IRCH (gwt2bndsu1.f :3508)."
        )
    if rts:
        opts.append("RTS")
        why.append(
            "With RTS, USG-T applies none of the RTS records read after the "
            "start of a period that reuses IRCH (glo2basu1.f :3342, :3353)."
        )
        why.append(_RTS_RECH_NOTE)
    warnings.warn(
        f"MfUsgRch: {' and '.join(opts)} with IRCH reused in stress period "
        f"{periods}: "
        + " ".join(why)
        + " Set inirch_ncpl=True to write the current IRCH in those periods, "
        "or remove the negative inirch of those periods."
    )


def _disu_irch_plan(
    model, rech, irch, rchconc, flags, nper, mxrtzones, inirch_ncpl, values
):
    """[INIRCH, IRCH file entry] of every stress period for DISU NRCHOP=2,
    and the [kper, INIRCH] of the periods that reuse IRCH with CONC or RTS.

    USG-T RP: INIRCH >= 0 is the number of IRCH nodes and reads IRCH; < 0
    reuses IRCH and its node count (gwf2rch8u1.f :272, :338-342), which also
    sizes RECH and the concentrations of the period (:318, :482). values
    holds the INIRCH < 0 read or given (inirch). A period without irch data
    writes, with CONC or RTS, the current IRCH with its node count, unless
    values has a value for it and inirch_ncpl is off (a reuse); without
    CONC and RTS it reuses IRCH (the value, or -1). A value in a period
    with irch data raises ValueError. Period 1 must read IRCH."""
    plan = []
    reuse = []
    nirch = None
    for kper in range(nper):
        u2d = irch.transient_2ds.get(kper)
        changed = False
        if u2d is not None:
            if kper in values:
                raise ValueError(
                    f"irch is given for stress period {kper + 1} but INIRCH "
                    f"{values[kper]} does not read it"
                )
            changed = nirch is not None and u2d.array.size != nirch
            nirch = u2d.array.size
            # RTS takes the zone of node IRCH(NN) from IZNRCH, read only in
            # positions 1..NIRCH (gwf2rch8u1.f :412, :438-439)
            if mxrtzones and u2d.array.size and u2d.array.max() + 1 > nirch:
                raise NotImplementedError(
                    f"irch of stress period {kper + 1} has node "
                    f"{u2d.array.max() + 1} above the {nirch} IRCH nodes: with "
                    "RTS, USG-T looks up the zone of that node in IZNRCH "
                    "positions it never reads and adds no RTS recharge there"
                )
            entry = Util2d(
                model, u2d.shape, np.int32, u2d.array + 1, u2d.name
            ).get_file_entry()
            plan.append([nirch, entry])
        elif kper == 0:
            raise ValueError(
                "irch has no data for stress period 1: USG-T reads IRCH in "
                "period 1 (NRCHOP=2)"
            )
        elif (flags is not None or mxrtzones) and (inirch_ncpl or kper not in values):
            plan.append([nirch, entry])
        else:
            n = values.get(kper, -1)
            if flags is not None or mxrtzones:
                reuse.append([kper, n])
            plan.append([n, ""])
        inrech = rech.get_kper_entry(kper)[0]
        if inrech >= 0 and rech[kper].array.size != nirch:
            raise ValueError(
                f"rech of stress period {kper + 1} has {rech[kper].array.size} "
                f"values; IRCH has {nirch} nodes"
            )
        if inrech < 0 and changed:
            raise ValueError(
                f"IRCH of stress period {kper + 1} has {nirch} nodes but RECH is "
                "reused with the previous node count; give rech for that period"
            )
        if flags is not None and rchconc is not None and kper in rchconc.transient_3ds:
            u3d = rchconc.transient_3ds[kper]
            for icomp, flag in enumerate(flags):
                if flag and u3d[icomp].array.size != nirch:
                    raise ValueError(
                        f"rchconc component {icomp + 1} of stress period "
                        f"{kper + 1} has {u3d[icomp].array.size} values; IRCH "
                        f"has {nirch} nodes"
                    )
    return plan, reuse


def _inirch_options(model, nrchop, nper, inirch, inirch_ncpl):
    """Check inirch (NRCHOP=2 only; values < 0 only on an unstructured grid,
    where INIRCH is the size of the irch array) and inirch_ncpl (NRCHOP=2
    only) and return the INIRCH values as {kper: int}."""
    if inirch_ncpl not in (False, True):
        raise ValueError(f"inirch_ncpl = {inirch_ncpl!r} must be True or False")
    values = _period_flags(inirch, nper, "inirch")
    if inirch is not None and nrchop != 2:
        raise ValueError("inirch needs nrchop=2 (structured grid or DISU)")
    if not model.structured:
        given = {k: n for k, n in values.items() if n >= 0}
        if given:
            raise ValueError(
                f"inirch {given} on an unstructured grid: INIRCH is the size "
                "of the irch array of the period; inirch only takes values < 0 "
                "(reuse irch)"
            )
    if inirch_ncpl and nrchop != 2:
        raise ValueError("inirch_ncpl needs nrchop=2 (structured grid or DISU)")
    return values


def _structured_inirch(irch, values, nper, ncpl, conc, inirch_ncpl, mxrtzones):
    """[INIRCH, kper of the irch written] of every stress period for a
    structured grid with NRCHOP=2, the [kper, INIRCH] of the periods that
    read IRCH with CONC or RTS and an INIRCH below NROW*NCOL, and the
    [kper, INIRCH] of the periods that reuse IRCH with CONC or RTS.

    USG-T RP: INIRCH >= 0 reads IRCH for the NROW*NCOL cells, < 0 reuses it
    (gwf2rch8u1.f :338-360). Without CONC and RTS the result does not depend
    on INIRCH. Transport applies the recharge concentration to the first
    INIRCH cells only (gwt2bndsu1.f :3521) and to none for INIRCH <= 0
    (:3508); the RTS records read during a period update the recharge of
    the first INIRCH cells only (glo2basu1.f ATS1AJ :3342, :3353), none for
    INIRCH < 0. A period with irch data writes NROW*NCOL with inirch_ncpl,
    else the INIRCH read (values), else NROW*NCOL with CONC or RTS and 1
    without; a period without data writes, with CONC or RTS, the current
    irch with NROW*NCOL when it has no INIRCH value or with inirch_ncpl,
    else the INIRCH read or given (< 0, a reuse), or -1 without a value;
    period 1 must read IRCH."""
    full = conc or bool(mxrtzones)
    plan = []
    low = []
    reuse = []
    src = None
    for kper in range(nper):
        n = values.get(kper)
        if kper in irch.transient_2ds:
            src = kper
            # RTS takes the zone of cell IRCH(NN) from IZNRCH, which has
            # NROW*NCOL positions (gwf2rch8u1.f :103, :107, :356, :438-439)
            layer = irch.transient_2ds[kper].array.max() + 1
            if mxrtzones and layer > 1:
                raise NotImplementedError(
                    f"irch of stress period {kper + 1} has layer {layer}: with "
                    "RTS, USG-T looks up the zone of a cell below layer 1 in "
                    "IZNRCH positions it never reads and adds no RTS recharge "
                    "there"
                )
            if n is not None and n < 0:
                raise ValueError(
                    f"irch is given for stress period {kper + 1} but INIRCH "
                    f"{n} does not read it"
                )
            if inirch_ncpl:
                n = ncpl
            elif n is None:
                n = ncpl if full else 1
            # IRCH has NROW*NCOL positions (gwf2rch8u1.f :102, :152)
            if full and n > ncpl:
                raise ValueError(
                    f"INIRCH {n} in stress period {kper + 1} is above NROW*NCOL "
                    f"= {ncpl}: USG-T would read IRCH positions that do not "
                    "exist"
                )
            if full and n < ncpl:
                low.append([kper, n])
        else:
            if kper == 0:
                raise ValueError(
                    "irch has no data for stress period 1: USG-T reads IRCH in "
                    "period 1 (NRCHOP=2)"
                )
            if n is not None and n >= 0:
                raise ValueError(
                    f"INIRCH {n} in stress period {kper + 1} needs irch data"
                )
            if full and (inirch_ncpl or n is None):
                n = ncpl
            else:
                n = -1 if n is None else n
                if full:
                    reuse.append([kper, n])
        plan.append([n, src])
    return plan, low, reuse


class MfUsgRch(ModflowRch):
    """
    MFUSG TRANSPORT Recharge Package Class.

    Parameters
    ----------
    model : model object
        The model object (of type :class:`flopy.modflow.mf.Modflow`) to which
        this package will be added.
    ipakcb : int, optional
        Toggles whether cell-by-cell budget data should be saved. If None or zero,
        budget data will not be saved (default is None).
    nrchop : int
        is the recharge option code.
        1: Recharge to top grid layer only
        2: Recharge to layer defined in irch
        3: Recharge to highest active cell (default is 3).
    rech : float or filename or ndarray or dict keyed on kper (zero-based)
        Recharge flux (default is 1.e-3, which is used for all stress periods)
    irch : int or filename or ndarray or dict keyed on kper (zero-based)
        Layer (for an unstructured grid) or node (for an unstructured grid) to
        which recharge is applied in each vertical column (only used when
        nrchop=2). Default is 0, which is used for all stress periods.
    extension : string
        Filename extension (default is 'rch')
    unitnumber : int
        File unit number (default is None).
    filenames : str or list of str
        Filenames to use for the package and the output files. If
        filenames=None the package name will be created using the model name
        and package extension and the cbc output name will be created using
        the model name and .cbc extension (for example, modflowtest.cbc),
        if ipakcb is a number greater than zero. If a single string is passed
        the package will be set to the string and cbc output names will be
        created using the model name and .cbc extension, if ipakcb is a
        number greater than zero. To define the names for all package files
        (input and output) the length of the list of strings should be 2.
        Default is None.
    inselev : dict {kper: int} or one int per stress period, optional
        INSELEV of each stress period (SEEPELEV option): > 0 reads selev of
        the period and turns seepage on, 0 turns it off, < 0 keeps the
        previous state and SELEV. Periods without a value write 1 when selev
        has an array for the period and -1 otherwise (default is None).
    iniznrch : dict {kper: int} or one int per stress period, optional
        INRCHZONES of each stress period (RTS option): > 0 reads iznrch of
        the period; the period 1 value is the number of RTS zones (it may
        differ from mxrtzones, which USG-T does not use then) and later
        values > 0 must keep it; <= 0 reuses the zones. Periods without a
        value write that number (mxrtzones in period 1) when iznrch has an
        array for the period and no keyword otherwise. Period 1 must read
        zones (default is None). With nrchop=2 on an unstructured grid USG-T
        takes the zone of IRCH node n from position n of IZNRCH, not from the
        position of the entry; a node above the number of IRCH nodes is not
        written (NotImplementedError: that zone is never read). On a
        structured grid IZNRCH has NROW*NCOL positions, so with nrchop=2 an
        irch layer above 1 is not written (NotImplementedError: USG-T adds no
        RTS recharge to that cell). With nrchop=2, a period without irch
        data writes the current irch again with its full count; a period
        with a negative inirch (read or given) is written as a reuse with a
        warning (USG-T applies none of the RTS records read after the start
        of that period, glo2basu1.f :3342, :3353), and inirch_ncpl=True
        writes the current irch there instead. A loaded file that also
        reuses RECH (INRECH < 0) in that period does not keep its result
        when written again (see inirch).
    iconc : int
        1 writes the CONC option: recharge concentrations of the components
        selected by irchconc (needs the component counts of a BCT package, or
        counts set by hand on a model without BCT; default is 0).
    rchconc : float or ndarray or dict keyed on kper (zero-based), optional
        Recharge concentrations, one array per component (MCOMPT), only with
        iconc=1. A period with data writes INCONC and the arrays of the
        components with flag 1; other periods reuse them. Period 1 needs data
        when a flag is 1 (default is None). On an unstructured grid with
        nrchop=2, a period whose irch changes the number of nodes and has no
        rchconc data reuses the previous concentrations by position: USG-T
        gives entry n of the new IRCH the concentration kept in position n
        (an engine behavior; FloPy does not check it).
    irchconc : list of int, optional
        IRCHCONC, one value 0 or 1 per component (MCOMPT); only with iconc=1.
        Default is 1 for every component.
    inirch : dict {kper: int} or one int per stress period, optional
        INIRCH of each stress period, nrchop=2 only: >= 0 reads irch of the
        period, < 0 reuses it. Load keeps the values read: every period on a
        structured grid, the values < 0 (reuse) on an unstructured grid, so
        a loaded flow-only file keeps its reuse when CONC or RTS is added
        later (with the warning below). On an unstructured grid INIRCH is
        the size of the irch array, so inirch only takes values < 0 (a value
        >= 0, also 0, raises ValueError); give it as a dict, as the list
        form needs a value for period 1, which reads irch. On a structured
        grid a period with irch data and without a value writes NROW*NCOL
        with iconc=1 or RTS (mxrtzones > 0) and 1 otherwise. A period
        without irch data and without a value writes, with CONC or RTS, the
        current irch again with its full count (NROW*NCOL, or its node count
        on an unstructured grid), and -1 (a reuse) without them. Without
        CONC and RTS the USG-T result does not depend on the value. With
        CONC, transport applies the recharge concentration only to the first
        INIRCH cells (gwt2bndsu1.f :3521) and to none for INIRCH <= 0
        (:3508); with RTS, the RTS records read during a period update the
        recharge of the first INIRCH cells only (glo2basu1.f :3342, :3353),
        none for INIRCH < 0, so the flow depends on INIRCH. With CONC or RTS
        a value below NROW*NCOL gives a warning at write time, and so does a
        negative value (a reuse: USG-T adds no recharge mass in that period
        with CONC and applies no RTS record read during it with RTS);
        inirch_ncpl=True writes the current irch with its full count there,
        and so does inirch=None on a loaded package (on a structured grid
        it also drops the low INIRCH values read). Rebuilding a loaded
        package with MfUsgRch(...) without passing inirch treats it as a new
        model. A key of inirch at or above the number of stress periods (for
        example after nper is reduced on a loaded package) raises ValueError
        at write time; a stress period added later, without irch data and
        without a value, is written as in a new model. With CONC or RTS the
        current irch is written as an INTERNAL array in every period without
        irch data, so the file grows with the number of stress periods; a
        negative inirch for those periods writes a reuse instead, with the
        consequences above. A negative value in a period with irch data, or
        in period 1, raises ValueError at write time, and so does a period 1
        without irch data (also a loaded file with INIRCH < 0 in period 1).
        A loaded file with RTS and a stress period that reuses RECH (INRECH
        < 0) with an INIRCH below NROW*NCOL, or that reuses IRCH (INIRCH <
        0, also on an unstructured grid), does not keep its result when
        written again: like upstream FloPy, load keeps RECH as data of every
        period and write_file reads it again in each one (INRECH 1), so
        USG-T applies the current RTS record to every IRCH cell at the start
        of that period (gwf2rch8u1.f :429-442) and resets the recharge of
        the cells the RTS records did not update; the warnings say so
        (default is None).
    inirch_ncpl : bool
        nrchop=2 only. On a structured grid, True writes INIRCH = NROW*NCOL
        in every stress period that reads irch, whatever inirch holds. With
        CONC or RTS, True also writes the current irch with its full count
        (NROW*NCOL, or its node count on an unstructured grid) in every
        period with a negative inirch (read or given); a period without irch
        data and without a value already does so. It exists to correct a
        file whose INIRCH would apply the recharge concentration or the RTS
        records to part of the cells, or to none in a period that reuses
        irch; it can be set on a loaded package before write_file (default
        is False: write inirch and the reuse). With RTS, a loaded file that
        reused RECH in those periods still does not keep its result (see
        inirch).

    Attributes
    ----------

    Methods
    -------

    See Also
    --------

    Notes
    -----
    Parameters are supported in Flopy only when reading in existing models.
    Parameter values are converted to native values in Flopy and the
    connection to "parameters" is thus nonexistent.

    Examples
    --------

    """

    def __init__(
        self,
        model,
        nrchop=3,
        ipakcb=None,
        rech=1e-3,
        irch=0,
        seepelev=0,
        mxrtzones=0,
        iconc=0,
        selev=None,
        iznrch=None,
        rchconc=None,
        unitnumber=None,
        filenames=None,
        inselev=None,
        iniznrch=None,
        irchconc=None,
        inirch=None,
        inirch_ncpl=False,
    ):
        """Package constructor."""
        msg = (
            "Model object must be of type flopy.mfusg.MfUsg\n"
            f"but received type: {type(model)}."
        )
        assert isinstance(model, MfUsg), msg

        # SELEV/IZNRCH are checked before the package is added to the model:
        # options and keys first, then the arrays converted without changing
        # the given objects
        nper = model.nrow_ncol_nlay_nper[3]
        _rts_seepelev_options(
            nper, mxrtzones, seepelev, iznrch, selev, iniznrch, inselev
        )
        selev = _as_transient2d(model, selev, np.float32, "rech_selev")
        iznrch = _as_transient2d(model, iznrch, np.int32, "rech_izn")
        _rts_seepelev_plan(
            nper,
            get_pak_vals_shape(model, rech),
            mxrtzones,
            seepelev,
            iznrch,
            selev,
            iniznrch,
            inselev,
        )
        # IRCHCONC/RCHCONC: MCOMPT from the component counts (ValueError
        # without them), checked before the package is added to the model
        mcompt = _conc_mcompt(model, "MfUsgRch") if iconc else None
        if rchconc is not None and iconc:
            rchconc = Transient3d(
                model,
                (mcompt,) + get_pak_vals_shape(model, rech),
                np.float32,
                rchconc,
                name="rech_conc",
            )
        irchconc = _conc_flags(iconc, mcompt, irchconc, rchconc)
        _inirch_options(model, nrchop, nper, inirch, inirch_ncpl)

        # call base package constructor
        super().__init__(
            model,
            nrchop=nrchop,
            ipakcb=ipakcb,
            rech=rech,
            irch=irch,
            unitnumber=unitnumber,
            filenames=filenames,
        )

        self.mxrtzones = mxrtzones
        self.seepelev = seepelev
        self.iconc = iconc
        self.selev = selev
        self.iznrch = iznrch
        self.inselev = inselev
        self.iniznrch = iniznrch

        self.irchconc = irchconc
        self.rchconc = rchconc
        self._built_mcompt = mcompt
        self.inirch = inirch
        self.inirch_ncpl = inirch_ncpl

    def write_file(self, check=True, f=None):
        """
        Write the package file.

        Parameters
        ----------
        check : boolean
            Check package data for common errors. (default True)

        Returns
        -------
        None

        """
        nrow, ncol, nlay, nper = self.parent.nrow_ncol_nlay_nper
        # checked again (attributes may have changed) before any output
        plan = _rts_seepelev_plan(
            nper,
            self.rech.shape,
            self.mxrtzones,
            self.seepelev,
            self.iznrch,
            self.selev,
            self.iniznrch,
            self.inselev,
        )
        mcompt = None
        if self.iconc:
            mcompt = _conc_mcompt(self.parent, "MfUsgRch")
            if self._built_mcompt is not None and mcompt != self._built_mcompt:
                raise ValueError(
                    f"MfUsgRch: MCOMPT changed from {self._built_mcompt} to "
                    f"{mcompt} after this package was built"
                )
        flags = _conc_flags(self.iconc, mcompt, self.irchconc, self.rchconc)
        values = _inirch_options(
            self.parent, self.nrchop, nper, self.inirch, self.inirch_ncpl
        )
        irch_plan = None
        reuse = []
        if self.nrchop == 2 and not self.parent.structured:
            irch_plan, reuse = _disu_irch_plan(
                self.parent,
                self.rech,
                self.irch,
                self.rchconc,
                flags,
                nper,
                self.mxrtzones,
                self.inirch_ncpl,
                values,
            )
        inirch_plan = None
        if self.nrchop == 2 and self.parent.structured:
            inirch_plan, low, reuse = _structured_inirch(
                self.irch,
                values,
                nper,
                nrow * ncol,
                flags is not None,
                self.inirch_ncpl,
                self.mxrtzones,
            )
            if low:
                periods = ", ".join(f"{k + 1} (INIRCH {n})" for k, n in low)
                opts = []
                why = []
                if flags is not None:
                    opts.append("CONC")
                    why.append(
                        "USG-T applies the recharge concentration only to the "
                        "first INIRCH cells (gwt2bndsu1.f :3521), so the "
                        "recharge of the other cells enters without mass and "
                        "dilutes the concentrations, and the mass budget does "
                        "not show it because it uses the same count (:3586, "
                        ":3599)."
                    )
                if self.mxrtzones:
                    opts.append("RTS")
                    why.append(
                        "USG-T applies the RTS records read during a period "
                        "only to the recharge of the first INIRCH cells "
                        "(glo2basu1.f :3353), so the other cells keep their "
                        "previous recharge and the flow is not the one of "
                        "INIRCH = NROW*NCOL."
                    )
                    why.append(_RTS_RECH_NOTE)
                warnings.warn(
                    f"MfUsgRch: {' and '.join(opts)} with INIRCH below NROW*NCOL "
                    f"= {nrow * ncol} in stress period {periods}: "
                    + " ".join(why)
                    + " Set inirch_ncpl=True to write INIRCH = NROW*NCOL."
                )
        if reuse:
            _irch_reuse_warning(reuse, flags is not None, self.mxrtzones)
        # allows turning off package checks when writing files at model level
        if check:
            self.check(
                f=f"{self.name[0]}.chk",
                verbose=self.parent.verbose,
                level=1,
            )
        # the text is built in memory; the file is opened at the end
        f_rch = io.StringIO()
        f_rch.write(f"{self.heading}\n")
        f_rch.write(f"{self.nrchop:10d}{self.ipakcb:10d}")
        if self.seepelev:
            f_rch.write(" SEEPELEV")
        if self.iconc:
            f_rch.write(" CONC")
        if self.mxrtzones:
            f_rch.write(f" RTS {self.mxrtzones:4d}")
        f_rch.write("\n")

        if self.nrchop == 2:
            if self.parent.structured:
                irch = {}
                for kper, u2d in self.irch.transient_2ds.items():
                    irch[kper] = u2d.array + 1
                irch = Transient2d(
                    self.parent,
                    self.irch.shape,
                    self.irch.dtype,
                    irch,
                    self.irch.name,
                )
            else:
                mxndrch = np.max(
                    [u2d.array.size for kper, u2d in self.irch.transient_2ds.items()]
                )
                f_rch.write(f"{mxndrch:10d}\n")

        if self.iconc:
            for flag in flags:
                f_rch.write(f"{flag:10d}")
            f_rch.write("\n")

        for kper in range(nper):
            inrech, file_entry_rech = self.rech.get_kper_entry(kper)
            if self.nrchop == 2:
                if self.parent.structured:
                    inirch, src = inirch_plan[kper]
                    file_entry_irch = irch.get_kper_entry(src)[1]
                else:
                    inirch, file_entry_irch = irch_plan[kper]
            else:
                inirch = -1

            if self.nrchop == 2:
                f_rch.write(f"{inrech:10d}{inirch:10d} ")
            else:
                f_rch.write(f"{inrech:10d} ")

            inrchzones, inselev = plan[kper]
            if inrchzones is not None:
                f_rch.write(f" INRCHZONES {inrchzones}")
            if inselev is not None:
                f_rch.write(f" INSELEV {inselev}")
            # INCONC only for a period with data; RP reuses them otherwise
            inconc = self.rchconc is not None and kper in self.rchconc.transient_3ds
            if inconc:
                f_rch.write(" INCONC")
            f_rch.write(f" # Stress period {kper + 1}\n")

            if inrech >= 0:
                f_rch.write(file_entry_rech)
            if self.nrchop == 2:
                if inirch >= 0:
                    f_rch.write(file_entry_irch)

            # arrays only for INRCHZONES/INSELEV > 0 (RP reads them then)
            if inrchzones is not None and inrchzones > 0:
                f_rch.write(self.iznrch.transient_2ds[kper].get_file_entry())
            if inselev is not None and inselev > 0:
                f_rch.write(self.selev.transient_2ds[kper].get_file_entry())
            # one array per component with IRCHCONC 1, in component order
            if inconc:
                u3d = self.rchconc.transient_3ds[kper]
                for icomp, flag in enumerate(flags):
                    if flag:
                        f_rch.write(u3d[icomp].get_file_entry())
        # a handle given by the caller is written and left open
        if f is not None:
            f.write(f_rch.getvalue())
        else:
            with open(self.fn_path, "w") as f_out:
                f_out.write(f_rch.getvalue())

    @classmethod
    def load(cls, f, model, nper=None, ext_unit_dict=None, check=True):
        """
        Load an existing package.

        Parameters
        ----------
        f : filename or file handle
            File to load.
        model : model object
            The model object (of type :class:`flopy.modflow.mf.Modflow`) to
            which this package will be added.
        nper : int
            The number of stress periods.  If nper is None, then nper will be
            obtained from the model object. (default is None).
        ext_unit_dict : dictionary, optional
            If the arrays in the file are specified using EXTERNAL,
            or older style array control records, then `f` should be a file
            handle.  In this case ext_unit_dict is required, which can be
            constructed using the function
            :class:`flopy.utils.mfreadnam.parsenamefile`.
        check : boolean
            Check package data for common errors. (default True)

        Returns
        -------
        rch : MfUsgRch object
            MfUsgRch object.

        Examples
        --------

        >>> import flopy
        >>> m = flopy.modflow.Modflow()
        >>> rch = flopy.modflow.MfUsgRch.load('test.rch', m)

        """
        if model.verbose:
            print("loading rch package file...")

        openfile = not hasattr(f, "read")
        if openfile:
            filename = f
            f = open(filename, "r")

        # dataset 0 -- header
        while True:
            line = f.readline()
            if line[0] != "#":
                break
        npar = 0

        if "parameter" in line.lower():
            raw = line.strip().split()
            npar = int(raw[1])
            if npar > 0:
                if model.verbose:
                    print(f"   Parameters detected. Number of parameters = {npar}")
            line = f.readline()
        # dataset 2
        t = line_parse(line)
        nrchop = int(t[0])
        ipakcb = int(t[1])

        # item 2 - options
        # whole words in any case on the whole line, as URWORD reads them
        # (also after '#')
        words = _words(line)
        seepelev = 0
        if "SEEPELEV" in words:
            seepelev = 1

        iconc = 0
        if "CONC" in words or "CONCENTRATION" in words:
            iconc = 1

        mxrtzones = 0
        if "RTS" in words:
            idx = words.index("RTS")
            mxrtzones = int(words[idx + 1])

        # dataset 2b for mfusg
        if not model.structured and nrchop == 2:
            line = f.readline()
            t = line_parse(line)
            mxndrch = int(t[0])

        irchconc = None
        if iconc:
            try:
                mcomp = _conc_mcompt(model, "MfUsgRch.load")
            except ValueError:
                if openfile:
                    f.close()
                raise
            irchconc = np.empty((mcomp), dtype=np.int32)
            if model.verbose:
                print(f"   loading IRCHCONC[{mcomp}] ...")
            irchconc = read1d(f, irchconc)

        if model.verbose:
            print(f"   irchconc: {irchconc}")

        # dataset 3 and 4 - parameters data
        pak_parms = None
        if npar > 0:
            pak_parms = mfparbc.loadarray(f, npar, model.verbose)

        if nper is None:
            nrow, ncol, nlay, nper = model.get_nrow_ncol_nlay_nper()
        else:
            nrow, ncol, nlay, _ = model.get_nrow_ncol_nlay_nper()

        u2d_shape = (nrow, ncol)

        # read data for every stress period
        current_rech = []
        rech = {}
        irch = None
        if nrchop == 2:
            irch = {}
        current_irch = []
        # NRCHOP=2: the INIRCH read in every period
        inirch_read = {}

        # arrays only of the periods that read them (flag > 0); flags as read,
        # INSELEV 0 for a period without the keyword (RP sets ISELEV = 0)
        iznrch = {}
        iniznrch = {}
        selev = {}
        inselev = {}

        rchconc = None
        if iconc:
            rchconc = {}
        inconc = [0] * nper
        current_rchconc = []

        for iper in range(nper):
            line = f.readline()
            t = line_parse(line)
            inrech = int(t[0])

            if nrchop == 2:
                inirch = int(t[1])
                inirch_read[iper] = inirch
                if (not model.structured) and (inirch >= 0):
                    u2d_shape = (1, inirch)
            elif not model.structured:
                u2d_shape = (1, ncol[0])

            words = _words(line)
            inselev[iper] = 0
            if "INSELEV" in words:
                idx = words.index("INSELEV")
                inselev[iper] = int(words[idx + 1])

            if "INRCHZONES" in words:
                idx = words.index("INRCHZONES")
                iniznrch[iper] = int(words[idx + 1])

            if "INCONC" in words:
                inconc[iper] = 1

            if inrech >= 0:
                if npar == 0:
                    if model.verbose:
                        print(f"   loading rech stress period {iper + 1:3d}...")
                    t = Util2d.load(
                        f,
                        model,
                        u2d_shape,
                        np.float32,
                        "rech",
                        ext_unit_dict,
                    )
                else:
                    parm_dict = {}
                    for ipar in range(inrech):
                        line = f.readline()
                        t = line.strip().split()
                        pname = t[0].lower()
                        try:
                            c = t[1].lower()
                            instance_dict = pak_parms.bc_parms[pname][1]
                            if c in instance_dict:
                                iname = c
                            else:
                                iname = "static"
                        except:
                            iname = "static"
                        parm_dict[pname] = iname
                    t = mfparbc.parameter_bcfill(model, u2d_shape, parm_dict, pak_parms)

                current_rech = t
            rech[iper] = current_rech

            if nrchop == 2:
                if inirch >= 0:
                    if model.verbose:
                        print(f"   loading irch stress period {iper + 1:3d}...")
                    t = Util2d.load(
                        f, model, u2d_shape, np.int32, "irch", ext_unit_dict
                    )
                    current_irch = Util2d(
                        model, u2d_shape, np.int32, t.array - 1, "irch"
                    )
                # a period that reuses IRCH (INIRCH < 0) has no data
                if inirch >= 0:
                    irch[iper] = current_irch

            # Item 9 for mfusg transport
            if mxrtzones:
                if iniznrch.get(iper, 0) > 0:
                    if model.verbose:
                        print(f"   loading iznrch stress period {iper + 1:3d}...")
                    iznrch[iper] = Util2d.load(
                        f, model, u2d_shape, np.int32, "iznrch", ext_unit_dict
                    )

            # Item 10 for mfusg transport
            if seepelev:
                if inselev[iper] > 0:
                    if model.verbose:
                        print(f"   loading selev stress period {iper + 1:3d}...")
                    selev[iper] = Util2d.load(
                        f, model, u2d_shape, np.float32, "selev", ext_unit_dict
                    )

            # Item 11 for mfusg transport
            if iconc:
                if inconc[iper]:
                    if model.verbose:
                        print(f"   loading rch conc stress period {iper + 1:3d}...")

                    # components with IRCHCONC 0 are not read
                    current_rchconc = [0.0] * mcomp
                    for icomp in range(mcomp):
                        if irchconc[icomp]:
                            if model.verbose:
                                print(
                                    f"   loading rch conc stress period {iper + 1:3d}"
                                    f" component {icomp + 1}..."
                                )
                            current_rchconc[icomp] = Util2d.load(
                                f,
                                model,
                                u2d_shape,
                                np.float32,
                                "rchconc",
                                ext_unit_dict,
                            )
                    # only the periods with INCONC (the others reuse them)
                    rchconc[iper] = current_rchconc

        if openfile:
            f.close()

        unitnumber, filenames = get_unitnumber_from_ext_unit_dict(
            model, cls, ext_unit_dict, ipakcb
        )
        # inirch: every period on a structured grid; on DISU only the reuse
        # (INIRCH < 0), the other counts come from the irch arrays
        inirch_kept = None
        if nrchop == 2 and model.structured:
            inirch_kept = inirch_read
        elif nrchop == 2:
            inirch_kept = {k: n for k, n in inirch_read.items() if n < 0} or None

        # create recharge package instance
        rch = cls(
            model,
            nrchop=nrchop,
            ipakcb=ipakcb,
            rech=rech,
            irch=irch,
            seepelev=seepelev,
            mxrtzones=mxrtzones,
            iconc=iconc,
            selev=selev or None,
            iznrch=iznrch or None,
            rchconc=rchconc or None,
            unitnumber=unitnumber,
            filenames=filenames,
            # without SEEPELEV/RTS, RP ignores INSELEV/INRCHZONES
            inselev=inselev if seepelev else None,
            iniznrch=iniznrch if mxrtzones else None,
            irchconc=irchconc,
            inirch=inirch_kept,
        )
        if check:
            rch.check(
                f=f"{rch.name[0]}.chk",
                verbose=rch.parent.verbose,
                level=0,
            )
        return rch
