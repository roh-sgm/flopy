"""gridgen2gsf module — build a GSF from Gridgen / DISV geometry.

Non-interactive helper **inspired by** the ``GRIDGEN2GSF`` utility
(``gridgen2gsf.f90``) for turning an already-built DISV-style / Gridgen geometry
into a MODFLOW-USG Grid Specification File.

Scope and references:

* Primary format spec: gwutil_a section 2.17 (``MODFLOW-USG Grid Specification
  File``).
* ``gridgen2gsf.f90`` is consulted only as a reference for the two GSF vertex
  layouts it emits (vertex-parsimonious and non-parsimonious).
* This helper does **not** parse GRIDGEN2GSF interactive prompts or its
  definition / quadtree input files, and does **not** reproduce the Fortran's
  grid construction (refinement, thresholds, rotation, offsets, quadtree
  structure). It starts from a geometry the caller already has: a
  ``disv_gridprops`` dict, a flopy ``Gridgen`` object, or an ``UnstructuredGrid``.

It is a thin, **separate** front-end: it accepts those geometry sources and
delegates GSF authoring (vertex layout, validation, 0-based/1-based ids,
``to_grid``) to the :class:`flopy.mfusg.MfUsgGsf` package, which it does not
modify. For the parsimonious/shared mode it first drops vertices not used by any
(surviving) cell and compacts the ids — the GRIDGEN2GSF vertex-parsimonious
behaviour.
"""

from __future__ import annotations

import numpy as np

from .mfusgsf import MfUsgGsf


def gridgen_to_gsf(
    model,
    source,
    top=1.0,
    botm=0.0,
    vertex_mode="shared",
    skip_degenerate=False,
    header="UNSTRUCTURED",
    extra_header=None,
):
    """Build an :class:`~flopy.mfusg.MfUsgGsf` from a Gridgen / DISV geometry.

    Parameters
    ----------
    model : flopy.mfusg.MfUsg
        Parent model the GSF package is attached to.
    source : dict or flopy.utils.Gridgen or flopy.discretization.UnstructuredGrid
        Geometry to convert:

        * a MODFLOW 6 ``disv_gridprops`` dict (keys ``"vertices"`` and
          ``"cell2d"``);
        * a flopy ``Gridgen`` object — its ``get_gridprops_disv()`` is used;
        * an ``UnstructuredGrid``.
    top, botm : float or array-like, optional
        The GSF top/bottom surfaces: a scalar applied to every vertex, or a
        per-vertex array. Defaults ``1.0`` / ``0.0`` (a unit-thickness slab, as
        in the teaching notebook). Real elevations must be supplied explicitly —
        a Gridgen / DISV source carries no unambiguous per-vertex z (the
        elevation is collapsed to per-cell top/botm).
    vertex_mode : str, optional
        ``"shared"`` / ``"parsimonious"`` (default) reuses vertex ids between
        neighbouring cells; ``"cell"`` / ``"nonparsimonious"`` gives every cell
        its own unique top/bottom vertices (8 per quad), never shared. Both write
        each node's top vertices then its bottom vertices, so
        ``to_grid()`` / ``from_gridspec(..., split_vertices=True)`` recovers
        top/botm.
    skip_degenerate : bool, optional
        DISV / Gridgen sources only: drop cells with fewer than 3 unique
        vertices and renumber the remaining nodes (default ``False`` -> raise).
    header : str, optional
        ``"UNSTRUCTURED"`` (default) or ``"UNSTRUCTURED GWF"``.
    extra_header : sequence of int, optional
        The line-2 ``IZ IC`` flags (default ``(1, 1)``; spec 2.17 requires 1 1).

    Returns
    -------
    flopy.mfusg.MfUsgGsf
        A semantic GSF package attached to ``model``. Call ``write_file()`` to
        write it and ``to_grid()`` to recover an ``UnstructuredGrid``.

    Notes
    -----
    Node ids are contiguous and 0-based internally / 1-based in the file; cells
    with a duplicated closing vertex are handled automatically. All validation
    is performed by ``MfUsgGsf``. In the shared/parsimonious mode, vertices not
    used by any surviving cell are dropped and the ids compacted before
    authoring (per-vertex ``top``/``botm`` arrays are remapped alongside);
    scalar ``top``/``botm`` pass through unchanged. The cell mode emits unique
    per-cell vertices, so unused source vertices never appear and no compaction
    is needed.
    """
    from ..discretization.unstructuredgrid import UnstructuredGrid

    # A flopy Gridgen object -> use its DISV grid properties.
    if hasattr(source, "get_gridprops_disv"):
        source = source.get_gridprops_disv()

    if isinstance(source, dict):
        if "vertices" not in source or "cell2d" not in source:
            raise ValueError(
                "disv_gridprops must contain 'vertices' and 'cell2d' keys."
            )
        if MfUsgGsf._canon_vertex_mode(vertex_mode) == "shared":
            source, top, botm = _compact_shared_disv(source, top, botm, skip_degenerate)
        return MfUsgGsf.from_disv_gridprops(
            model,
            source,
            top=top,
            botm=botm,
            skip_degenerate=skip_degenerate,
            vertex_mode=vertex_mode,
            header=header,
            extra_header=extra_header,
        )

    if isinstance(source, UnstructuredGrid):
        nverts = np.asarray(source.verts).shape[0]
        top_zverts = np.broadcast_to(np.asarray(top, dtype=float), (nverts,))
        bot_zverts = np.broadcast_to(np.asarray(botm, dtype=float), (nverts,))
        return MfUsgGsf.from_grid(
            model,
            source,
            top_zverts=top_zverts,
            bot_zverts=bot_zverts,
            vertex_mode=vertex_mode,
            header=header,
            extra_header=extra_header,
        )

    raise TypeError(
        "gridgen_to_gsf 'source' must be a disv_gridprops dict, a flopy Gridgen "
        f"object, or an UnstructuredGrid; got {type(source).__name__}."
    )


def _compact_shared_disv(disv_gridprops, top, botm, skip_degenerate):
    """Drop unused vertices and compact ids for the shared/parsimonious mode.

    Mirrors the GRIDGEN2GSF vertex-parsimonious pass: only vertices referenced by
    a surviving cell are kept, and their ids are remapped to ``0..n_used-1``. A
    cell dropped by ``skip_degenerate`` does not retain its exclusive vertices.
    Per-vertex ``top``/``botm`` arrays are remapped to the kept subset; scalars
    are returned unchanged. Returns ``(disv_gridprops, top, botm)`` — the inputs
    unchanged when every vertex is already used.
    """
    verts = disv_gridprops["vertices"]
    cell2d = disv_gridprops["cell2d"]
    nvert = len(verts)

    used = set()
    surviving = []  # (record, stripped 0-based vertex ids) in original order
    for rec in cell2d:
        ncvert = int(rec[3])
        ids = MfUsgGsf._strip_closing([int(v) for v in rec[4 : 4 + ncvert]])
        if skip_degenerate and len(set(ids)) < 3:
            # This cell will be dropped by from_disv_gridprops; its exclusive
            # vertices must not be retained.
            continue
        used.update(ids)
        surviving.append((rec, ids))

    if len(used) == nvert:
        return disv_gridprops, top, botm  # nothing unused; no remap needed

    used_sorted = sorted(used)
    remap = {old: new for new, old in enumerate(used_sorted)}

    new_vertices = [
        (new, verts[old][1], verts[old][2]) for new, old in enumerate(used_sorted)
    ]
    new_cell2d = [
        [rec[0], rec[1], rec[2], len(ids), *[remap[v] for v in ids]]
        for rec, ids in surviving
    ]

    new_disv = dict(disv_gridprops)
    new_disv["vertices"] = new_vertices
    new_disv["cell2d"] = new_cell2d
    if "nvert" in new_disv:
        new_disv["nvert"] = len(new_vertices)

    new_top = top if np.ndim(top) == 0 else np.asarray(top)[used_sorted]
    new_botm = botm if np.ndim(botm) == 0 else np.asarray(botm)[used_sorted]
    return new_disv, new_top, new_botm
