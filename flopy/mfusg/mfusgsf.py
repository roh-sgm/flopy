"""mfusgsf module. Contains the MfUsgGsf class.

Grid Specification File (GSF) — USG-Transport / MODFLOW-USG.

Format (unstructured)::

    UNSTRUCTURED [GWF]                 <- header (1 or 2 tokens)
    NNODES  NLAY  [extra ...]          <- node count, layer count, gridgen flags
    NVERTS                             <- number of vertices
    X(1)  Y(1)  Z(1)                   <- one line per vertex (NVERTS lines)
    ...
    NODENO XC YC ZC LAY NVERT V1 .. VN <- one line per node (NNODES lines)

``NODENO`` and the per-cell vertex ids ``V1..VN`` are 1-based in the file.
The GSF is a grid-specification file produced by gridgen-style tools and
consumed by post-processors (FloPy's ``UnstructuredGrid.from_gridspec``); it is
**not** read by the MODFLOW-USG / USG-T solver, so there is no Fortran reader to
audit. ``from_gridspec`` only reads ``NNODES`` from line 2 and derives the layer
count from the per-node ``LAY`` column, so the trailing line-2 integers (``1 1``
in observed gridgen output) are preserved verbatim on load and default to
``(1, 1)`` when authoring.

This class supports two mutually-exclusive modes:

* **Semantic** — build from Python with ``vertices`` + ``node_data`` (the
  primary authoring workflow; no existing ``.gsf`` required) or recover it from
  a file with ``MfUsgGsf.load(..., parse=True)``. Node and vertex references are
  0-based in the Python API and 1-based in the written file; layers are likewise
  0-based internally and 1-based in the file.
* **Raw round-trip** — ``MfUsgGsf.load`` defaults to ``parse=False`` and keeps
  the file ``lines`` verbatim for byte-faithful rewrite (what ``MfUsg.load``
  uses, and the fallback when ``parse=True`` meets an unrecognised file).

``to_grid()`` returns a fully-parsed ``UnstructuredGrid`` via
``UnstructuredGrid.from_gridspec`` in either mode.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np

from ..pakbase import Package
from .mfusg import MfUsg


def _fmt(v):
    """Format a float so it reloads exactly (repr round-trips Python floats)."""
    return repr(float(v))


class MfUsgGsf(Package):
    """MODFLOW-USG Grid Specification File (GSF) Package.

    Parameters
    ----------
    model : flopy.mfusg.MfUsg
        Parent model.
    lines : list of str, optional
        Raw lines of the GSF file (raw round-trip mode). Stored verbatim and
        written unchanged.
    vertices : array-like, optional
        Semantic mode: ``(nverts, 3)`` array of vertex ``(x, y, z)`` coordinates,
        indexed by 0-based vertex id.
    node_data : list of dict, optional
        Semantic mode: one dict per node with keys ``"xc"``, ``"yc"``,
        ``"layer"`` (0-based), ``"vertices"`` (0-based vertex ids), and optional
        ``"node"`` (0-based, defaults to list position) and ``"zc"`` (defaults
        to 0.0). ``vertices`` and ``node_data`` must be supplied together.
    header : str, optional
        Header line for semantic mode (default ``"UNSTRUCTURED"``; may be
        ``"UNSTRUCTURED GWF"``).
    nlay : int, optional
        Number of layers for semantic mode. Defaults to ``max(layer) + 1``.
    extra_header : sequence of int, optional
        Trailing integers on line 2 after ``NNODES NLAY`` (gridgen flags that
        the grid consumer ignores). Defaults to ``(1, 1)``.
    extension : str, optional
        File extension (default "gsf").
    unitnumber : int, optional
        File unit number (default 110).
    filenames : str or list of str, optional
        Package filename override.

    Notes
    -----
    Provide either raw ``lines`` or the semantic ``vertices`` + ``node_data``,
    never both — supplying both raises ``ValueError`` (the modes are mutually
    exclusive, and ``write_file`` selects the semantic branch on
    ``node_data is not None``, so it never silently prefers one). With none of
    them the package writes an empty file.
    """

    def __init__(
        self,
        model,
        lines=None,
        vertices=None,
        node_data=None,
        header="UNSTRUCTURED",
        nlay=None,
        extra_header=None,
        extension="gsf",
        unitnumber=None,
        filenames=None,
    ):
        msg = (
            "Model object must be of type flopy.mfusg.MfUsg\n"
            f"but received type: {type(model)}."
        )
        assert isinstance(model, MfUsg), msg

        # Mutually exclusive modes, judged by explicit presence (is not None) so
        # validation matches the test write_file uses to pick a branch.
        raw_provided = lines is not None
        semantic_provided = vertices is not None or node_data is not None
        if raw_provided and semantic_provided:
            raise ValueError(
                "MfUsgGsf accepts either raw 'lines' or semantic "
                "'vertices'+'node_data', not both."
            )
        if semantic_provided and (vertices is None or node_data is None):
            raise ValueError("Semantic GSF requires both 'vertices' and 'node_data'.")

        if unitnumber is None:
            unitnumber = MfUsgGsf._defaultunit()

        filenames = self._prepare_filenames(filenames, num=1)

        super().__init__(
            model,
            extension=extension,
            name=self._ftype(),
            unit_number=unitnumber,
            filenames=filenames,
        )
        self._generate_heading()

        if semantic_provided:
            verts, nodes, nnodes, nlay = self._normalize_semantic(
                vertices, node_data, nlay
            )
            self.vertices = verts
            self.node_data = nodes
            self.nnodes = nnodes
            self.nlay = nlay
            self.header = header
            self.extra_header = (
                tuple(int(x) for x in extra_header)
                if extra_header is not None
                else (1, 1)
            )
            self.lines = None
        else:
            self.lines = list(lines) if lines else []
            self.vertices = None
            self.node_data = None
            self.nnodes = 0
            self.nlay = None
            self.header = None
            self.extra_header = None

        self.parent.add_package(self)

    # ------------------------------------------------------------------
    # Authoring / validation
    # ------------------------------------------------------------------
    @staticmethod
    def _normalize_semantic(vertices, node_data, nlay):
        """Validate semantic GSF inputs; return (verts, nodes, nnodes, nlay)."""
        verts = np.asarray(vertices, dtype=float)
        if verts.ndim != 2 or verts.shape[1] != 3:
            raise ValueError(
                "GSF 'vertices' must be an (nverts, 3) array of (x, y, z); "
                f"got shape {verts.shape}."
            )
        nverts = verts.shape[0]

        nodes = []
        for i, rec in enumerate(node_data):
            node = int(rec.get("node", i))
            if node < 0:
                raise ValueError(f"GSF node id must be 0-based (>= 0); got {node}.")
            if "xc" not in rec or "yc" not in rec:
                raise ValueError("GSF node_data requires 'xc' and 'yc'.")
            if "vertices" not in rec:
                raise ValueError("GSF node_data requires 'vertices'.")
            layer = int(rec.get("layer", 0))
            if layer < 0:
                raise ValueError(f"GSF layer must be 0-based (>= 0); got {layer}.")
            vids = [int(v) for v in rec["vertices"]]
            for v in vids:
                if v < 0 or v >= nverts:
                    raise ValueError(
                        f"GSF node {node} references vertex {v}, out of range "
                        f"[0, {nverts - 1}] (vertices are 0-based)."
                    )
            nodes.append(
                {
                    "node": node,
                    "xc": float(rec["xc"]),
                    "yc": float(rec["yc"]),
                    "zc": float(rec.get("zc", 0.0)),
                    "layer": layer,
                    "vertices": vids,
                }
            )

        nnodes = len(nodes)
        if nlay is None:
            nlay = max((r["layer"] for r in nodes), default=-1) + 1
        return verts, nodes, nnodes, int(nlay)

    @classmethod
    def from_grid(
        cls, model, grid, zverts=None, header="UNSTRUCTURED", extra_header=None
    ):
        """Build a semantic GSF from an ``UnstructuredGrid``.

        Per-vertex Z elevations are **not** retained by ``UnstructuredGrid``
        (``from_gridspec`` collapses them into per-cell top/botm), so they must
        be supplied via ``zverts`` (one value per grid vertex); otherwise this
        raises rather than inventing elevations. The per-node ``zc`` is taken as
        the cell midpoint ``(top + botm) / 2`` from the grid, and the layer of
        each cell is derived from ``grid.ncpl``.
        """
        verts_xy = np.asarray(grid.verts, dtype=float)
        nverts = verts_xy.shape[0]
        if zverts is None:
            raise ValueError(
                "from_grid requires per-vertex z elevations via 'zverts' "
                f"(length {nverts}); UnstructuredGrid does not retain vertex z "
                "(from_gridspec collapses it into per-cell top/botm)."
            )
        zverts = np.asarray(zverts, dtype=float).ravel()
        if len(zverts) != nverts:
            raise ValueError(
                f"'zverts' has length {len(zverts)}; expected {nverts} (one per "
                "grid vertex)."
            )
        vertices = [(verts_xy[i, 0], verts_xy[i, 1], zverts[i]) for i in range(nverts)]

        top = np.asarray(grid.top, dtype=float).ravel()
        bot = np.asarray(grid.botm, dtype=float).ravel()
        xc = np.asarray(grid.xcellcenters, dtype=float).ravel()
        yc = np.asarray(grid.ycellcenters, dtype=float).ravel()
        ncpl = np.asarray(grid.ncpl).ravel()
        layers = np.repeat(np.arange(len(ncpl)), ncpl)  # 0-based layer per cell

        node_data = [
            {
                "node": c,
                "xc": float(xc[c]),
                "yc": float(yc[c]),
                "zc": float((top[c] + bot[c]) / 2.0),
                "layer": int(layers[c]),
                "vertices": [int(v) for v in iv],
            }
            for c, iv in enumerate(grid.iverts)
        ]
        return cls(
            model,
            vertices=vertices,
            node_data=node_data,
            header=header,
            nlay=len(ncpl),
            extra_header=extra_header,
        )

    # ------------------------------------------------------------------
    # write_file
    # ------------------------------------------------------------------
    def write_file(self, check=False):
        """Write the GSF file (semantic if authored, else raw lines)."""
        with open(self.fn_path, "w") as f:
            if self.node_data is not None:
                self._write_semantic(f)
            else:
                for line in self.lines:
                    f.write(line if line.endswith("\n") else line + "\n")

    def _write_semantic(self, f):
        """Write semantic vertex/node data as a valid GSF (1-based file ids)."""
        f.write(f"{self.header}\n")
        extra = " ".join(str(int(x)) for x in self.extra_header)
        line2 = f"{self.nnodes} {self.nlay}"
        if extra:
            line2 += " " + extra
        f.write(line2 + "\n")
        f.write(f"{len(self.vertices)}\n")
        for x, y, z in self.vertices:
            f.write(f"{_fmt(x)} {_fmt(y)} {_fmt(z)}\n")
        for rec in self.node_data:
            verts1 = " ".join(str(v + 1) for v in rec["vertices"])
            f.write(
                f"{rec['node'] + 1} {_fmt(rec['xc'])} {_fmt(rec['yc'])} "
                f"{_fmt(rec['zc'])} {rec['layer'] + 1} "
                f"{len(rec['vertices'])} {verts1}\n"
            )

    # ------------------------------------------------------------------
    # to_grid
    # ------------------------------------------------------------------
    def to_grid(self):
        """Return an ``UnstructuredGrid`` parsed from this GSF.

        Requires the GSF file to exist on disk; if it was created
        programmatically it is written first.

        Returns
        -------
        flopy.discretization.UnstructuredGrid
        """
        from ..discretization.unstructuredgrid import UnstructuredGrid

        path = self.fn_path
        if not Path(path).exists():
            self.write_file()
        return UnstructuredGrid.from_gridspec(path, split_vertices=True)

    # ------------------------------------------------------------------
    # load
    # ------------------------------------------------------------------
    @classmethod
    def load(cls, f, model, nper=None, ext_unit_dict=None, check=False, parse=False):
        """Load an existing GSF file.

        By default (``parse=False``) the file ``lines`` are kept verbatim for a
        byte-faithful round-trip. With ``parse=True`` the file is parsed into
        semantic ``vertices`` + ``node_data``; if it cannot be parsed (or is
        inconsistent) the loader falls back to the raw round-trip instead of
        writing partial data.
        """
        msg = (
            "Model object must be of type flopy.mfusg.MfUsg\n"
            f"but received type: {type(model)}."
        )
        assert isinstance(model, MfUsg), msg

        if model.verbose:
            print("loading gsf package file...")

        if hasattr(f, "read"):
            fh = f
        else:
            fh = open(f, "r")

        try:
            lines = fh.readlines()
        finally:
            if not hasattr(f, "read"):
                fh.close()

        unitnumber = None
        filenames = [None]
        if ext_unit_dict is not None:
            unitnumber_ext, fname_ext = model.get_ext_dict_attr(
                ext_unit_dict, filetype=cls._ftype()
            )
            if unitnumber_ext is not None:
                unitnumber = unitnumber_ext
            if fname_ext is not None:
                filenames = [fname_ext]

        if parse:
            try:
                kw = cls._parse_semantic(lines)
                return cls(model, **kw, unitnumber=unitnumber, filenames=filenames)
            except (ValueError, IndexError) as e:
                if model.verbose:
                    print(f"  GSF semantic parse failed ({e}); keeping raw round-trip.")

        return cls(
            model,
            lines=lines,
            unitnumber=unitnumber,
            filenames=filenames,
        )

    @staticmethod
    def _parse_semantic(lines):
        """Parse raw GSF lines into semantic 0-based constructor kwargs.

        Raises ``ValueError``/``IndexError`` on any structural problem so
        ``load`` can fall back to the raw round-trip.
        """
        filtered = [
            ln for ln in lines if ln.strip() and not ln.lstrip().startswith("#")
        ]

        head_toks = filtered[0].split()
        if not head_toks or head_toks[0].upper() != "UNSTRUCTURED":
            raise ValueError(
                f"GSF must start with 'UNSTRUCTURED'; got {filtered[0].strip()!r}."
            )
        header = "UNSTRUCTURED" + (
            " GWF" if any(t.upper() == "GWF" for t in head_toks[1:]) else ""
        )

        line2 = filtered[1].split()
        nnodes = int(line2[0])
        nlay = int(line2[1]) if len(line2) > 1 else None
        extra_header = [int(t) for t in line2[2:]]

        nverts = int(filtered[2].split()[0])

        idx = 3
        vertices = []
        for _ in range(nverts):
            p = filtered[idx].split()
            idx += 1
            vertices.append((float(p[0]), float(p[1]), float(p[2])))

        node_data = []
        for _ in range(nnodes):
            p = filtered[idx].split()
            idx += 1
            nv = int(p[5])
            vids = [int(p[6 + k]) - 1 for k in range(nv)]
            node_data.append(
                {
                    "node": int(p[0]) - 1,
                    "xc": float(p[1]),
                    "yc": float(p[2]),
                    "zc": float(p[3]),
                    "layer": int(float(p[4])) - 1,
                    "vertices": vids,
                }
            )

        return {
            "vertices": vertices,
            "node_data": node_data,
            "header": header,
            "nlay": nlay,
            "extra_header": extra_header,
        }

    @staticmethod
    def _ftype():
        return "GSF"

    @staticmethod
    def _defaultunit():
        return 110
