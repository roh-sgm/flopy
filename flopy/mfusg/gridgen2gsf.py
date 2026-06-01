"""gridgen2gsf module — build a GSF from Gridgen / DISV geometry.

Non-interactive Python equivalent of the ``GRIDGEN2GSF`` utility
(``gridgen2gsf.f90``): convert a quadtree / DISV-style geometry into a
MODFLOW-USG Grid Specification File. The GSF format spec is gwutil_a section
2.17; the Fortran is used only as a reference for the two geometry layouts it
writes (vertex-parsimonious and non-parsimonious), not for its interactive flow.

This is a thin, **separate** front-end: it accepts several geometry sources and
delegates the GSF authoring (vertex layout, validation, 0-based/1-based id
handling, ``to_grid``) to the :class:`flopy.mfusg.MfUsgGsf` package. It does not
modify ``MfUsgGsf``.
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
    is performed by ``MfUsgGsf``.
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
