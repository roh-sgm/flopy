# Stage 4.2 - GSF Semantic Support

Goal: turn `MfUsgGsf` from a raw text wrapper plus `to_grid()` helper into a
semantic GSF package that can be authored from Python data and written as valid
USG-T/MODFLOW-USG grid-specification input.

## Current State

`MfUsgGsf` currently:

- stores raw GSF lines,
- writes them back verbatim,
- can call `UnstructuredGrid.from_gridspec()` via `to_grid()`,
- does not expose semantic vertex/node data,
- cannot author a GSF from arrays without manually building text lines.

Roadmap status: `Raw/text round-trip`.

## Source References

Primary sources:

- existing `flopy/mfusg/mfusgsf.py`,
- `flopy.discretization.UnstructuredGrid.from_gridspec`,
- USG-T examples that include GSF files,
- any Fortran reader path that consumes GSF/NAM entries if present.

GSF may be more of a grid-spec file than a package with heavy Fortran parsing,
so the implementation should lean on FloPy grid abstractions where possible.

## Design Requirements

- Provide semantic attributes for:
  - `structured` / `unstructured` header,
  - `nnodes`,
  - `nlay`,
  - vertices `(x, y, z)`,
  - per-node records: node id, center coordinates, layer, vertex ids.
- Internal node and vertex references should be 0-based where exposed through
  Python APIs.
- File output must be 1-based.
- Preserve raw round-trip support for existing files.
- Support construction from:
  - explicit arrays/lists, and
  - an existing `UnstructuredGrid` where practical.
- Avoid duplicating large grid parsing code if FloPy already provides a stable
  utility.

## Suggested API Shape

The agent should choose the final API after inspecting `UnstructuredGrid`, but a
reasonable target is:

```python
MfUsgGsf(
    model,
    vertices=[(x, y, z), ...],
    node_data=[
        {
            "node": 0,
            "xc": ...,
            "yc": ...,
            "zc": ...,
            "layer": 0,
            "vertices": [0, 1, 2],
        },
    ],
)
```

Optional helpers:

- `MfUsgGsf.from_grid(model, grid)`
- `gsf.to_grid()`
- `gsf.to_arrays()`

## Required Tests

Add focused tests in `autotest/test_usg_transport.py` or a new targeted file if
that is cleaner:

- semantic load of a minimal triangular GSF,
- from-scratch authoring of the same minimal GSF,
- reload confirms 0-based node/vertex references,
- written file uses 1-based ids,
- `to_grid()` still returns `UnstructuredGrid`,
- raw text round-trip remains protected,
- invalid vertex references fail clearly.

## Documentation Updates

Update:

- `USGT_roadmap.md`
- `USGT_PACKAGE_BACKLOG.md`
- `USGT_improvements.md`

Promote GSF only if semantic authoring and reload are tested. If the final
choice is to keep GSF raw/text, explain why and keep status unchanged.

## Validation

```bash
python -m pytest autotest/test_usg_transport.py -q
python -m pytest autotest/test_usg_transport_exe.py -q
USGT_EXE=/Users/roh.sgm/Documents/GitHub_Projects/GW-Software-Compiled-via-Claude/gfortran/usgt_2.7/usgt_270_arm python -m pytest autotest/test_usg_transport.py autotest/test_usg_transport_exe.py -q
git diff --check
git status --short
```

## Agent Prompt

```text
Goal: implement semantic GSF authoring/load/write support for USG-T 2.7 in FloPy.

Read first:
- USGT_STAGE4_MASTER_PLAN.md
- USGT_STAGE4_02_GSF_SEMANTIC.md
- USGT_roadmap.md
- USGT_PACKAGE_BACKLOG.md
- flopy/mfusg/mfusgsf.py
- flopy/discretization/unstructuredgrid.py
- autotest/test_usg_transport.py

Task:
Turn `MfUsgGsf` into a semantic grid-spec package. It must still preserve raw round-trip behavior, but should also support construction from Python arrays/lists and write valid GSF with 1-based file ids. Exposed Python node/vertex references should be 0-based.

Required tests:
- semantic minimal GSF load,
- from-scratch authoring + reload,
- 0-based internal / 1-based file ids,
- `to_grid()` still works,
- raw text round-trip still works,
- invalid references fail clearly.

Update docs:
- USGT_roadmap.md
- USGT_PACKAGE_BACKLOG.md
- USGT_improvements.md

Validation:
python -m pytest autotest/test_usg_transport.py -q
python -m pytest autotest/test_usg_transport_exe.py -q
USGT_EXE=/Users/roh.sgm/Documents/GitHub_Projects/GW-Software-Compiled-via-Claude/gfortran/usgt_2.7/usgt_270_arm python -m pytest autotest/test_usg_transport.py autotest/test_usg_transport_exe.py -q
git diff --check
git status --short

Constraints:
- Push only to origin, never upstream.
- Do not touch MF6-TID or unrelated packages.
- Do not promote GSF to Full unless semantic constructor, writer, load/reload, and tests are complete.
```

## Stage 4.2 resolution — done

`MfUsgGsf` was promoted from `Raw/text round-trip` to **`Full (authoring)`**.

Key finding: **GSF is not read by the USG-T solver.** A grep of the USG-T 2.7
source found no `.gsf` / GRIDSPEC reader — GSF is produced by gridgen-style
tools and consumed by post-processors (FloPy's
`UnstructuredGrid.from_gridspec`). So there is no Fortran reader to audit and an
executable smoke test is not meaningful; correctness is validated via
`from_gridspec`/`to_grid` and semantic write/reload.

Implemented:

- Semantic constructor `MfUsgGsf(model, vertices=[(x,y,z),...], node_data=[{node,
  xc, yc, zc, layer, vertices}], header=..., nlay=..., extra_header=...)`.
- `parse=True` loader (`_parse_semantic`) recovering the same structure, with
  raw fallback on any unparseable/inconsistent file.
- Semantic writer producing valid GSF with **1-based** node/vertex/layer ids
  from **0-based** Python data.
- `from_grid(model, grid, zverts)` — builds from an `UnstructuredGrid`; requires
  per-vertex `zverts` because the grid does not retain vertex z (collapsed into
  per-cell top/botm), failing explicitly rather than inventing elevations.
- `to_grid()` unchanged (delegates to `from_gridspec`).
- Mutually-exclusive modes: raw `lines` vs semantic `vertices`+`node_data`
  (judged by explicit presence, matching the `write_file` branch test); raw
  round-trip kept as the `load` default.
- Line-2 trailing gridgen integers preserved verbatim on load, default `(1, 1)`
  on authoring (the grid consumer ignores them).

Tests added in `autotest/test_usg_transport.py`:
`test_mfusggsf_semantic_load`, `test_mfusggsf_authoring_from_scratch`,
`test_mfusggsf_write_reload_roundtrip`,
`test_mfusggsf_rejects_invalid_and_mixed_modes`, `test_mfusggsf_from_grid`
(the three existing raw `lines`/`to_grid` tests remain green).

Validation: focused **114 passed**, exe **3 passed**, combined **117 passed**
under the USG-T 2.7 ARM binary; `git diff --check` clean.

## Stage 4.2 review follow-up — resolved

Incorporated the teaching notebook's practical GSF authoring pattern
(`Advanced3_geopandas.ipynb`, cells 29/30/36/55) and closed four findings. GSF
stays `Full (authoring)`.

1. **`UNSTRUCTURED GWF` header (P2).** `UnstructuredGrid.from_gridspec` had an
   operator-precedence bug (`not (A) or (B)`) that rejected the valid two-token
   header. Fixed to accept exactly `UNSTRUCTURED` / `UNSTRUCTURED GWF` and reject
   anything else; the freyberg GSF and comment-prefixed headers still parse.
2. **Trailing content under `parse=True` (P2).** `_parse_semantic` now raises if
   non-comment lines remain after the `NVERTS + NNODES` records, so `load`
   falls back to the raw round-trip rather than rewriting the file without the
   extra lines.
3. **Notebook top/bottom convention (P2).** `from_grid` gained
   `top_zverts` + `bot_zverts`: the writer doubles the vertices (all top, then
   all bottom) and each node lists its top ids followed by the matching bottom
   ids (top id + `totalverts`), so `from_gridspec(..., split_vertices=True)`
   reconstructs correct top/botm. Single-surface `zverts` is kept as an explicit
   legacy mode; a per-cell closing-duplicate vertex is dropped automatically.
4. **DISV mapping (P2/P3).** Added
   `from_disv_gridprops(model, disv_gridprops, top, botm, skip_degenerate=False)`
   building a single-layer GSF from a MODFLOW 6 DISV 2D template, with
   automatic closing-vertex removal and (optional) degenerate-cell skipping plus
   consecutive node renumbering.

Hardening in the semantic constructor: the header must be `UNSTRUCTURED` /
`UNSTRUCTURED GWF`; node ids must be unique; `nlay` must be `>= max(layer) + 1`.

Tests added (`autotest/test_usg_transport.py`):
`test_mfusggsf_unstructured_gwf_header`,
`test_mfusggsf_parse_rejects_trailing_content`,
`test_mfusggsf_from_grid_top_bottom`, `test_mfusggsf_from_disv_gridprops`,
`test_mfusggsf_hardening_rejects`.

Validation: focused **119 passed**, exe **3 passed**, combined **122 passed**
under the USG-T 2.7 ARM binary; `git diff --check` clean.
