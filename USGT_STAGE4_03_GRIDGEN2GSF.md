# Stage 4.3 - gridgen2gsf utility

Goal: a clean, non-interactive Python utility that builds a GSF file from
Gridgen / DISV-style geometry — the equivalent of the `GRIDGEN2GSF` program,
without its interactive prompts. This follows Stage 4.2 (GSF semantic support)
and is GSF-adjacent; it is numbered 03 as requested and is independent of the
master-plan's other Stage 4.3 package work.

## References

- Primary GSF spec: gwutil_a section 2.17 (`MODFLOW-USG Grid Specification
  File`).
- Secondary (geometry patterns only): `gridgen2gsf.f90` (`GRIDGEN2GSF`) — its
  two output layouts (vertex-parsimonious and non-parsimonious). The interactive
  flow and quadtree-file reading are intentionally **not** reproduced.

## Design (chosen)

A **separate module** `flopy/mfusg/gridgen2gsf.py` exposing one non-interactive
function:

```python
gridgen_to_gsf(model, source, top=1.0, botm=0.0, vertex_mode="shared",
               skip_degenerate=False, header="UNSTRUCTURED", extra_header=None)
    -> MfUsgGsf
```

It does **not** modify `MfUsgGsf` (Stage 4.2 stays closed). It is a thin
front-end that delegates the GSF authoring — vertex layout, validation,
0-based/1-based ids, `to_grid()` — to `MfUsgGsf`'s public authoring API
(`from_disv_gridprops` / `from_grid`), so there is no duplicated geometry logic.

`source` may be:

- a MODFLOW 6 `disv_gridprops` dict (`"vertices"` + `"cell2d"`),
- a flopy `Gridgen` object (its `get_gridprops_disv()` is used), or
- an `UnstructuredGrid`.

`top` / `botm` are the GSF top/bottom **surfaces** (scalar applied to all
vertices, or a per-vertex array). They default to `1.0` / `0.0` — a
unit-thickness slab, as in the teaching notebook — because a Gridgen / DISV
source has no unambiguous per-vertex z. For an `UnstructuredGrid` source the
scalars are broadcast to per-vertex `top_zverts` / `bot_zverts`.

### Modes

- `vertex_mode="shared"` / `"parsimonious"` (default): neighbouring cells reuse
  vertex ids (one shared `2*nverts` vertex set).
- `vertex_mode="cell"` / `"nonparsimonious"`: every cell owns unique top/bottom
  vertices (8 per quad), never shared. Caller polygon order is kept within each
  half; this is a non-shared generalization and is **not** byte-equivalent to
  the GRIDGEN2GSF quadtree winding. Only the top-half/bottom-half split (what
  `split_vertices=True` needs) is guaranteed.

### Validation (inherited from `MfUsgGsf`)

- node ids contiguous and ordered `0..nnodes-1`;
- `IZ IC` = `(1, 1)` (spec 2.17);
- degenerate cells (< 3 unique vertices) raise unless `skip_degenerate=True`
  (then dropped + nodes renumbered);
- 0-based internal node/vertex ids, 1-based in the file;
- a per-cell closing-duplicate vertex is dropped automatically.

## Tests (`autotest/test_usg_transport.py`)

- `test_gridgen_to_gsf_disv_modes` — parsimonious shares vertex ids between two
  adjacent quads; non-parsimonious does not (8 vertices/quad); both recover
  top/botm through `to_grid()` (`split_vertices=True`).
- `test_gridgen_to_gsf_source_types` — a Gridgen-like object
  (`get_gridprops_disv`) and an `UnstructuredGrid` source.
- `test_gridgen_to_gsf_validation` — invalid `vertex_mode`, invalid source type,
  malformed `disv_gridprops`, and degenerate-cell raise / `skip_degenerate`.

## Validation

```bash
python -m pytest autotest/test_usg_transport.py -k mfusggsf -q
python -m pytest autotest/test_usg_transport.py -q
python -m pytest autotest/test_usg_transport_exe.py -q
USGT_EXE=/Users/roh.sgm/Documents/GitHub_Projects/GW-Software-Compiled-via-Claude/gfortran/usgt_2.7/usgt_270_arm python -m pytest autotest/test_usg_transport.py autotest/test_usg_transport_exe.py -q
git diff --check
git status --short
```

Result: `-k mfusggsf` 17 passed; focused **126 passed**; exe **3 passed**;
combined **129 passed** under the USG-T 2.7 ARM binary; `git diff --check` clean.

## Status

Done. New module + 3 tests; `MfUsgGsf` unchanged. GSF is not solver input, so
there is no executable smoke test (validated through `to_grid` / `from_gridspec`
and the GSF write path).
