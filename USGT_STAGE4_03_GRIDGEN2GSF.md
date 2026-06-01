# Stage 4.3 - gridgen2gsf utility

Goal: a clean, non-interactive Python helper **inspired by** the `GRIDGEN2GSF`
program for building a GSF file from an already-built Gridgen / DISV-style
geometry. This follows Stage 4.2 (GSF semantic support) and is GSF-adjacent; it
is numbered 03 as requested and is independent of the master-plan's other Stage
4.3 package work.

## References

- Primary GSF spec: gwutil_a section 2.17 (`MODFLOW-USG Grid Specification
  File`).
- Secondary (geometry patterns only): `gridgen2gsf.f90` (`GRIDGEN2GSF`) — its
  two output layouts (vertex-parsimonious and non-parsimonious).

Scope: this helper is **not** a full port of `GRIDGEN2GSF`. It does not parse
the Fortran's interactive prompts or its definition / quadtree input files, and
it does not reproduce the Fortran's grid construction (refinement, thresholds,
rotation, offsets, quadtree structure). It starts from a geometry the caller
already has.

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
  vertex ids. Vertices not used by any surviving cell are dropped and the ids
  compacted before authoring (the GRIDGEN2GSF vertex-parsimonious pass);
  per-vertex `top`/`botm` arrays are remapped alongside, scalars pass through.
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
- `test_gridgen_to_gsf_parsimonious_compacts` — an unused source vertex is
  dropped and ids compacted; per-vertex `top`/`botm` arrays are remapped; cell
  mode also ignores the unused vertex (8 per quad).
- `test_gridgen_to_gsf_skip_degenerate_compacts` — vertices used only by a
  skipped degenerate cell are not retained.

## Validation

```bash
python -m pytest autotest/test_usg_transport.py -k gridgen_to_gsf -q   # focal
python -m pytest autotest/test_usg_transport.py -k mfusggsf -q          # GSF package
python -m pytest autotest/test_usg_transport.py -q
python -m pytest autotest/test_usg_transport_exe.py -q
USGT_EXE=/Users/roh.sgm/Documents/GitHub_Projects/GW-Software-Compiled-via-Claude/gfortran/usgt_2.7/usgt_270_arm python -m pytest autotest/test_usg_transport.py autotest/test_usg_transport_exe.py -q
git diff --check
git status --short
```

`-k mfusggsf` covers the GSF package; the utility's own tests are under
`-k gridgen_to_gsf` (the focal filter).

## Status

Done.

Follow-up (resolved):

1. **Parsimonious now compacts** unused vertices for DISV / Gridgen sources:
   before delegating, `gridgen_to_gsf` drops vertices not referenced by any
   surviving cell and remaps the `cell2d` ids (helper `_compact_shared_disv`).
   `skip_degenerate=True` drops a degenerate cell's exclusive vertices too;
   per-vertex `top`/`botm` arrays are remapped, scalars unchanged. The cell mode
   needs no compaction (it emits unique per-cell vertices).
2. **Docstring no longer overpromises**: "inspired by GRIDGEN2GSF", with the
   non-reproduction scope (interactive/definition files, refinement, thresholds,
   rotation, offsets, quadtree) stated explicitly.
3. Focal validation filter documented as `-k gridgen_to_gsf`.

Five tests (`..._disv_modes`, `..._source_types`, `..._validation`,
`..._parsimonious_compacts`, `..._skip_degenerate_compacts`); `MfUsgGsf`
unchanged. GSF is not solver input, so there is no executable smoke test
(validated through `to_grid` / `from_gridspec` and the GSF write path).
`-k gridgen_to_gsf` **5 passed**; `-k mfusggsf` **17 passed**; focused
**128 passed**; exe **3 passed**; combined **131 passed** under the USG-T 2.7
ARM binary; `git diff --check` clean.
