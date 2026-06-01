# USG-Transport improvements on this fork

This fork of [modflowpy/flopy](https://github.com/modflowpy/flopy) carries
USG-Transport additions on `develop`, aimed at improving how
`MfUsg.load(...).write_input()` round-trips existing USG-T models and
extending post-processing utilities for transport output. Nothing else is
changed.

Tested against **USG-T 1.8** and **USG-T 2.7** binaries on both structured and
unstructured grids.

These are incremental improvements on top of the USG-T work already in
upstream flopy. They live here while testing continues.

## What's added

### Stage 4 — programmatic authoring (in progress)

Guided by `USGT_STAGE4_MASTER_PLAN.md`. Closing authoring gaps so USG-T packages
can be built from Python/numpy without first loading an existing model.

- **Stage 4.1 — TIB semantic support (2026-05-31):** promoted `MfUsgTib` from
  `Raw/text round-trip` to **`Full (authoring)`**, reversing the Stage 3 Card 4
  "keep raw/text" decision. Audited `GWF2TIB1RP` in `glo2basu1.f` and added a
  semantic `stress_period_data` constructor plus a `parse=True` loader covering
  the full grammar: flow blocks (`NIB0` inactivate / `NIB1` activate / `NIBM1`
  prescribed-head, with `HEAD`/`AVHEAD`/bare) and transport blocks
  (`NICB0`/`NICB1`/`NICBM1`, with multi-component `CONC`/`AVCONC`/bare); `NIB0`/
  `NICB0` node lists via `U1DINT`; 3- vs 6-int header keyed to BCT presence;
  0-based internal / 1-based file. `MfUsgTib.load` still **defaults to byte-exact
  raw round-trip** (what `MfUsg.load` uses); `parse=True` returns semantic data
  and falls back to raw on `EXTERNAL`/`OPEN-CLOSE` `U1DINT` rather than writing
  partial data. Authoring validates transport-requires-BCT, `CONC` length =
  MCOMP, and 0-based nodes. Six new focused tests (from-scratch non-transport +
  transport authoring, semantic load, write/reload byte-stable, raw fallback on
  unsupported syntax, explicit authoring rejections) join the existing raw
  round-trip test, plus a USG-T 2.7 executable smoke test
  (`test_usgt_exe_tib_prescribed_head_from_scratch`) that proves a FloPy-authored
  TIB `NIBM1`+`HEAD` record actually bends the head solution when run. Focused
  suite **107 passed**, exe suite **3 passed**, combined **110 passed** under the
  USG-T 2.7 ARM binary.
  - **Review follow-up (resolved):** made the three input modes
    (`stress_period_data`/`blocks`/`raw_body`) mutually exclusive — more than one
    mode now raises `ValueError` instead of `write_file` silently preferring
    `raw_body`; and `parse=True` now raises on a premature EOF (fewer headers
    than `nper`) so `load` falls back to the raw body rather than padding the
    file with synthetic no-op stress periods. Two tests added; focused suite
    **109 passed**, combined **112 passed** under the ARM binary. TIB stays
    `Full (authoring)`.
  - **Polish (resolved):** the mode-exclusivity check counted modes with a
    truthiness test while `write_file` branched on `is not None`, so an
    explicitly-empty mode (`raw_body=""`, `stress_period_data={}`, `blocks={}`)
    could slip past validation and still steer the writer into a silent branch.
    Validation now uses explicit presence (`is not None`) to match the writer:
    a single explicitly-empty mode is valid, but mixing it with another mode
    raises. Regression cases added; focused suite **109 passed**, combined
    **112 passed** under the ARM binary.
- **Stage 4.2 — GSF semantic support (2026-05-31):** promoted `MfUsgGsf` from
  `Raw/text round-trip` to **`Full (authoring)`**. Added a semantic
  `vertices` + `node_data` constructor, a `parse=True` loader, a semantic writer,
  and `from_grid(model, grid, zverts)`, covering the GSF grid spec: header
  (`UNSTRUCTURED [GWF]`), `nnodes`/`nlay`, vertex `(x, y, z)`, and per-node
  `(node, xc, yc, zc, layer, vertices)`. Nodes/vertices/layers are 0-based in the
  Python API and 1-based in the file. GSF is **not** read by the USG-T solver
  (no Fortran reader — it is consumed by `UnstructuredGrid.from_gridspec`), so
  correctness is validated through `from_gridspec`/`to_grid` and semantic
  write/reload rather than an executable smoke test. `from_grid` requires
  per-vertex `zverts` because `UnstructuredGrid` does not retain vertex z
  (collapsed into per-cell top/botm) — it fails explicitly rather than inventing
  elevations. `MfUsgGsf.load` still **defaults to the byte-faithful raw `lines`
  round-trip** (what `MfUsg.load` uses), with raw fallback when `parse=True`
  meets an unparseable file; raw `lines` and semantic data are mutually
  exclusive. Five new tests (semantic load, from-scratch authoring, write/reload,
  invalid-ref + mixed-mode rejection, `from_grid`) join the three existing raw
  tests. Focused suite **114 passed**, exe suite **3 passed**, combined
  **117 passed** under the USG-T 2.7 ARM binary.
  - **Review follow-up (resolved):** (1) fixed an operator-precedence bug in
    `UnstructuredGrid.from_gridspec` that wrongly rejected the valid
    `UNSTRUCTURED GWF` header (now accepts exactly `UNSTRUCTURED` /
    `UNSTRUCTURED GWF`, rejects others); (2) `parse=True` now rejects trailing
    non-comment content after the node records, so such files fall back to the
    raw round-trip instead of being rewritten without the extra lines;
    (3) generalised `from_grid` to the USG-T **top/bottom doubled-vertex**
    convention (`top_zverts` + `bot_zverts`) that `from_gridspec(...,
    split_vertices=True)` reconstructs into correct top/botm — the teaching-
    notebook pattern — keeping single-surface `zverts` as an explicit legacy
    mode and auto-dropping per-cell closing-duplicate vertices; (4) added
    `from_disv_gridprops(model, disv_gridprops, top, botm, skip_degenerate=...)`
    mapping a MODFLOW 6 DISV 2D template to a single-layer GSF (closing-vertex
    removal, degenerate-cell skip + node renumbering). Hardening: header,
    unique-node-id, and `nlay >= max(layer)+1` validation in the semantic
    constructor. Five follow-up tests added; focused suite **119 passed**, exe
    **3 passed**, combined **122 passed** under the ARM binary. GSF stays
    `Full (authoring)`.
  - **Vertex-mode follow-up (resolved):** using `gridgen2gsf.f90` as a secondary
    reference for geometry generation (primary spec stays gwutil_a 2.17), added a
    `vertex_mode` option to `from_grid`/`from_disv_gridprops` emulating the two
    GRIDGEN2GSF layouts: `"shared"`/`"parsimonious"` (default — neighbouring
    cells reuse vertex ids) and `"cell"`/`"nonparsimonious"` (every cell owns
    unique top/bottom vertices, 8 per quad, never shared). Both write each node's
    top vertices then bottom vertices so `from_gridspec(split_vertices=True)`
    recovers top/botm; `"cell"` keeps the caller's polygon order within each half
    (not byte-equivalent to the gridgen quadtree `1,4,3,2 / 5,8,7,6` winding for
    non-quads). Also documented the line-2 `IZ IC` flags per spec 2.17 (both must
    be 1). One test added (`test_mfusggsf_vertex_modes`); focused suite
    **120 passed**, exe **3 passed**, combined **123 passed** under the ARM
    binary. GSF stays `Full (authoring)`.
  - **Final follow-up (resolved):** tightened the GSF semantic contract against
    spec 2.17. (1) `_parse_semantic` now reuses `_normalize_header`, so only
    `UNSTRUCTURED` / `UNSTRUCTURED GWF` parse and e.g. `UNSTRUCTURED EXTRA GWF`
    falls back to raw. (2) `IZ`/`IC` are validated: authoring defaults to
    `(1, 1)`, accepts an omitted value as the spec's assumed `1 1`, and rejects
    anything else or a bad length (e.g. `(1,1,9)`); on load, `nnode nlay` or
    `nnode nlay 1 1` parse semantically while `0 1`, `1 0`, or an odd token count
    fall back to raw. (3) node ids must be `0..nnodes-1` contiguous and ordered —
    authoring raises on a gap/duplicate/reorder, and a file whose `inode` column
    is not `1..nnode` in order falls back to raw. (4) clarified that
    `vertex_mode="cell"` is a non-shared *generalization* (caller polygon order
    within each half) and is **not** byte-equivalent to the GRIDGEN2GSF quadtree
    `1,4,3,2 / 5,8,7,6` winding; only the top/bottom-half split is guaranteed.
    Three tests added (`test_mfusggsf_parse_header_strict`,
    `test_mfusggsf_iz_ic_flags`, `test_mfusggsf_inode_validation`); focused suite
    **123 passed**, exe **3 passed**, combined **126 passed** under the ARM
    binary. GSF stays `Full (authoring)`.
- **Stage 4.3 — gridgen2gsf utility (2026-05-31):** added a clean,
  non-interactive `gridgen_to_gsf(model, source, top=, botm=, vertex_mode=, ...)`
  in a **separate** module `flopy/mfusg/gridgen2gsf.py` (exported from
  `flopy.mfusg`) — the Python equivalent of the `GRIDGEN2GSF` program without its
  interactive flow. It builds a `MfUsgGsf` from a MODFLOW 6 `disv_gridprops`
  dict, a flopy `Gridgen` object (`get_gridprops_disv()`), or an
  `UnstructuredGrid`, delegating all GSF authoring (vertex layout, validation,
  0-based/1-based ids, `to_grid`) to `MfUsgGsf`'s public API — so `MfUsgGsf` is
  **unchanged** (Stage 4.2 stays closed) and there is no duplicated geometry
  logic. Supports `vertex_mode="shared"`/`"parsimonious"` and
  `"cell"`/`"nonparsimonious"`; `top`/`botm` are the GSF top/bottom surfaces
  (default unit slab `1`/`0`, broadcast per-vertex for grid sources). Three tests
  added (DISV shared-vs-cell + `to_grid` top/botm; Gridgen-like + UnstructuredGrid
  sources; invalid mode/source/geometry + degenerate handling). Primary spec
  gwutil_a 2.17; `gridgen2gsf.f90` used only as a geometry reference. Focused
  suite **126 passed**, exe **3 passed**, combined **129 passed** under the
  USG-T 2.7 ARM binary.
  - **Review follow-up (resolved):** (1) the shared/parsimonious mode now
    **compacts** vertices for DISV/Gridgen sources — `gridgen_to_gsf` drops
    vertices not used by any surviving cell and remaps the `cell2d` ids before
    delegating (helper `_compact_shared_disv`), mirroring the GRIDGEN2GSF
    vertex-parsimonious pass; `skip_degenerate=True` also drops a degenerate
    cell's exclusive vertices, and per-vertex `top`/`botm` arrays are remapped
    (scalars unchanged); the cell mode needs no compaction. (2) the docstring no
    longer claims to be a "Python equivalent" — it is a helper *inspired by*
    GRIDGEN2GSF, and explicitly does not parse the Fortran's interactive /
    definition / quadtree files or reproduce its grid construction (refinement,
    thresholds, rotation, offsets, quadtree). (3) the utility's focal test filter
    is `-k gridgen_to_gsf`. Two tests added
    (`test_gridgen_to_gsf_parsimonious_compacts`,
    `test_gridgen_to_gsf_skip_degenerate_compacts`); `MfUsgGsf` unchanged.
    `-k gridgen_to_gsf` **5 passed**, focused suite **128 passed**, exe
    **3 passed**, combined **131 passed** under the ARM binary.
- **Stage 4 — OC Fullness Card A (2026-05-31):** audited `MfUsgOc` against
  USG-T 2.7 Output Control (`glo2basu1.f` SGWF2BAS7I/J/N) and hardened it.
  Implemented the `BOOTSTRAPPING` header — `__init__`/`write_file` now author and
  preserve it (previously `load` read the units but they were dropped on
  round-trip), written on the **first OC line** (the only place USG-T's setup
  parses it; the per-record reader rejects a standalone `BOOTSTRAPPING` line).
  Fixed `load` to keep **layer-qualified** `PRINT`/`SAVE HEAD`/`DRAWDOWN`/`CONC`
  actions (it previously truncated to two tokens) and to recapture
  **`DDREFERENCE`** from the period line. Added six tests (BOOTSTRAPPING header,
  per-SP `BOOTSTRAP`/`NOBOOTSTRAP`/`BOOTSTRAPSCALE`/`NOBOOTSTRAPSCALE`,
  CONC/BUDGET output blocks, `SAVE IBOUND`, layer-qualified, `DDREFERENCE`).
  **Decision: kept `✅ (intentionally not Full)`** with explicit gaps —
  `SAVE IBOUND` is commented out in USG-T 2.7 (solver rejects it; FloPy preserves
  the keyword); `FASTFORWARD`/`FASTFORWARDC` separate-line placement and
  `BOOTSTRAPPING` execution are not exe-verified; numeric-format OC rewrites as
  words. EVT/MDT/LAK untouched. See `USGT_STAGE4_OC_FULLNESS.md`. `-k mfusgoc`
  **7 passed**, focused suite **134 passed**, exe **3 passed**, combined
  **137 passed** under the USG-T 2.7 ARM binary.
  - **Polish follow-up (resolved):** (1) `write_file` now emits a
    DDREFERENCE-only period line (`period k step s ddreference`); previously such
    a period was dropped and the flag leaked into the next written period.
    (2) `check()` recognises the valid USG-T OC actions — the single-word
    toggles `BOOTSTRAP`/`NOBOOTSTRAP`/`BOOTSTRAPSCALE`/`NOBOOTSTRAPSCALE` and
    `DDREFERENCE`, and the keyword-value params `DELTAT`/`TMINAT`/`TMAXAT`/
    `TADJAT`/`TCUTAT`/`HCLOSE`/`BTOL`/`MXITER` **when a value is given** (a bare
    parametric keyword is flagged, since `SGWF2BAS7N` reads a value after it) —
    instead of emitting false "ignored" warnings, while still flagging genuinely
    unknown actions. (3) `check()` emits a specific warning that `SAVE IBOUND` is
    preserved by FloPy but rejected by USG-T 2.7 (write/load still preserve it).
    Four tests added; `-k mfusgoc` **11 passed**, focused suite **138 passed**,
    exe **3 passed**, combined **141 passed** under the ARM binary. OC stays
    `✅ (intentionally not Full)`.
- **Stage 4 — EVT Fullness Card B (2026-06-01):** audited `MfUsgEvt` against
  `gwf2evt8u1.f` and hardened it. Fixed a real round-trip bug — `load` read
  `IETFACTOR` but never passed it to the constructor, so it reset to 0 (and the
  `ETFACTOR` array was lost on the next write) — and made `write_file` emit the
  dataset-1 3-integer line (`NEVTOP IEVTCB IETFACTOR`) whenever transport is
  active and the `ETFACTOR` array when `ietfactor>0`. Added validation
  (`NEVTOP` 1..3, `ETFACTOR` length == MCOMP, `NEVTOP=2` structured `IEVT` layer
  range) and made the unsupported `ETS` zonal time-series fail explicitly
  (`NotImplementedError`). Six synthetic tests cover NEVTOP=1/2/3 (structured +
  unstructured `NEVTOP=2` with `MXNDEVT`; `IEVT` 0-based internal / 1-based
  file), per-SP reuse, transport `IETFACTOR` 0/<0/>0 with per-`MCOMP`
  `ETFACTOR`, and the negative cases; plus a USG-T 2.7 executable smoke
  (`test_usgt_exe_evt_from_scratch`). **Decision: kept `✅ (intentionally not
  Full)`** — gaps: ETS zonal time-series (unsupported) and `NPEVT` parameters
  (expanded-valid-write; note ETS array parameters are now *preserved* as of
  Stage 4.4A, but EVT's `NPEVT` still expands). OC/MDT/LAK untouched. See
  `USGT_STAGE4_EVT_FULLNESS.md`. `-k mfusgevt` **7 passed**, focused suite
  **144 passed**, exe **4 passed**, combined **148 passed** under the USG-T 2.7
  ARM binary.
  - **Review follow-up (resolved):** (1) a scalar `ETFACTOR` (e.g.
    `etfactor=2.5` with MCOMP=1) passed the length check but crashed
    `write_file` at `self.etfactor[icomp]`; `__init__` now normalizes
    `self.etfactor` to a 1-D array, so scalar (MCOMP=1) and array (MCOMP>1) both
    author/load/reload and an incompatible length still raises. (2) unstructured
    `NEVTOP=2` `IEVT` is now validated as a 0-based node index in `[0, NODES-1]`
    (file 1-based) via `node_count`; `< 0` or `>= NODES` raises (previously only
    the structured layer index was checked). Two tests added; `-k mfusgevt`
    **9 passed**, focused suite **146 passed**, exe **4 passed**, combined
    **150 passed** under the ARM binary. EVT stays `✅ (intentionally not Full)`.
- **Stage 4 — MDT Fullness Card C (2026-06-01):** audited `MfUsgMdt` against
  `gwt2mdtu1.for` and hardened it. Fixed three real bugs: `load` dropped
  `FRAHK`/`FRADARCY` (it upper-cased the header line but matched the lowercase
  keywords, so e.g. the Ex7 Multispecies `FRAHK` was lost on round-trip);
  `write_file(f=handle)` raised `NameError` (`f_obj` was only assigned when
  `f is None`, and it also closed a caller-supplied handle); and a stray
  `print(f"line={line}")` debug line. Added from-scratch authoring + round-trip
  for every main branch (header options parsed/written only when `IDPF==0`; base
  arrays with `VOLFRACMD` dropped when `IDPF!=0`; per-species `KDMD`/`DECAYMD`/
  `YIELDMD`/`DIFFMD`; `AIOLD1MD`/`AIOLD2MD` under `TSHIFTMD>0`; multi-species),
  plus validation (FRAHK⊕FRADARCY, IDPF-options, `MULTIFILE_MD` needs
  `imdtcf>0`, per-`MCOMP` list lengths → clear `ValueError`). Eight synthetic
  tests; the three Ex7 real-model round-trip/run tests are kept (no separate
  from-scratch MDT exe smoke — the Ex7 runs already exercise MDT execution).
  **Decision: kept `✅ (intentionally not Full)`** — gaps: the species loop uses
  `MCOMP` (chained-decay `NTCOMP>MCOMP` unverified) and the AI1/AI2 output
  binaries are authored but not read. OC/EVT/LAK untouched. See
  `USGT_STAGE4_MDT_FULLNESS.md`. `-k mfusgmdt` **8 passed**, focused suite
  **154 passed**, exe **4 passed**, combined **158 passed** under the USG-T 2.7
  ARM binary.
  - **Review follow-up (resolved):** fixed a TSHIFTMD misalignment. USG-T reads
    AIOLD1MD/AIOLD2MD only when `TSHIFTMD > 1e-10`, but the writer used
    `tshiftmd > 0` for the AIOLD arrays and a fixed `{:9.2f}` keyword format that
    rounded a small valid value (e.g. `1e-6`) to `0.00` — so the solver read
    TSHIFTMD as 0.0, skipped AIOLD, and the file desynced. A module constant
    `MDT_TSHIFT_THRESHOLD = 1e-10` is now used consistently (constructor, writer,
    IDPF-options check, load), and the keyword is written with a general format
    (`{:15.7g}`) so small valid values survive. One test added
    (`test_mfusgmdt_tshiftmd_threshold`); `-k mfusgmdt` **9 passed**, focused
    suite **155 passed**, exe **4 passed**, combined **159 passed** under the ARM
    binary. MDT stays `✅ (intentionally not Full)`.

### Stage 3 — upstream-readiness pass (2026-05-31)

Guided by `USGT_STAGE3_COMPLETION_PLAN.md`. Card-by-card hardening toward an
honest, upstream-ready USG-T 2.7 story.

- **Card 1 — status audit:** added a per-package decision table to
  `USGT_roadmap.md`; relabeled `TIB` and `GSF` as `Raw/text round-trip` (no
  semantic constructor by design) rather than a bare checkmark.
- **Card 2 — parameter strategy (ETS/HFB/SGB/QRT/DRT):** original decision was
  to keep `Expanded valid write` (ETS array parameters expand to concrete arrays
  and write `NPETS=0`) and explicit `NotImplementedError` for list parameters
  (`NPSGB`/`NPQRT`/`NPDRT` via `UPARLSTAL`, `NPHFB`). Fortran audit: `parutl7.f`
  (`UPARLSTAL`/`UPARLSTRP`/`UPARLSTSUB` list params) and `mfparbc` (ETS/EVT
  array params). Nothing writes incomplete parametric syntax. → **Superseded for
  ETS by Stage 4.4A** (below): ETS array parameters are now *preserved* on
  write/reload; the list-parameter packages still keep their explicit failures.
- **Stage 4.4A — ETS array-parameter preservation (executed):** ETS no longer
  expands by default — `MfUsgEts.load` preserves the parsed parameter
  definitions (`self.parameters`) and per-period activation records
  (`self.evtr_parm`), and `write_file` re-emits `NPETS>0` in item 2a, the
  definition blocks, and the activation records (incl. `INSTANCES`), while
  ETSS/ETSX/IETS/PXDP/PETM stay plain arrays (the mix the Fortran allows). USG-T
  reads `NPETS` from item 2a (`UPARARRAL` called with `IN=-1`), so no
  `PARAMETER` line is written. New shared helper
  `flopy/mfusg/_usgt_parameters.py` adds the array-parameter **write** side and
  reuses `ModflowParBc` as the parser (no second parser; no ad-hoc strings).
  Opt-in `expand_parameters=True` keeps the legacy expanded path; from-scratch
  parameter *authoring* (`npets>0` without loaded defs) still raises
  `NotImplementedError`, so ETS stays ⚠️ Partial / not `Full`. Tests: `-k
  mfusgets` **10 passed** (4 new preservation/instances tests; the old
  expand test now opts in via `expand_parameters=True`); focused **169**, exe
  **4**, combined **173** (ARM). The list-parameter write path was then
  implemented for HFB (Stage 4.4B, below); SGB/DRT/QRT remain pending. See
  `USGT_STAGE4_04_PARAMETERS_ETS.md`.
- **Stage 4.4B — HFB list-parameter preservation (executed):** HFB uses MODFLOW
  *list* parameters (`UPARLSTRP`/`UPARLSTSUB`), a different grammar from ETS:
  each parameter owns `NLST` barrier rows and its value scales `HYDCHR`.
  `MfUsgHfb.load` now preserves parameterized files — it stores the definitions
  (`self.parameters = {name: {"partyp","parval","nlst","data"}}`, barrier rows
  0-based), the non-parametric barriers (`self.hfb_data`), and the
  active-parameter names (`self.acthfb_names`); `write_file` re-emits the
  definition blocks, the non-parametric barriers, and `NACTHFB` + active names
  (barrier rows 1-based on file). The shared helper
  `flopy/mfusg/_usgt_parameters.py` gains list-parameter header/activation
  read+write functions (reusable by SGB/DRT/QRT later); the package keeps its own
  barrier-row reader/writer. Not promoted to `Full`: from-scratch parameter
  authoring, `TRANSIENT_HFB`+`NPHFB>0` (the Fortran would redefine params each SP
  under `ITERP=1`), and parameter `INSTANCES` (the Fortran aborts) all raise
  `NotImplementedError`; non-parametric HFB is unchanged. Tests: `-k mfusghfb`
  **8 passed** (5 new: unstructured/structured round-trip, params+non-param mix,
  from-scratch fail, transient+params fail; 3 non-parametric regressions kept);
  focused **173**, exe **4**, combined **177** (ARM). See
  `USGT_STAGE4_04_PARAMETERS_HFB.md`.
- **Stage 4.4B follow-up — HFB list controls (executed):** review found that a
  valid HFB file may begin a barrier list with `SFAC`/`OPEN/CLOSE`/`EXTERNAL`
  (per `SGWF2HFB7RL`/`SGWF2HFB7RLU`), which the previous loader crashed on with a
  raw `ValueError`. `MfUsgHfb._read_hfb_rows` now consumes these via the shared
  `_usgt_list.begin_list_block` — the same helper SGB/QRT/DRT use, **reused
  unchanged** — for both the per-parameter `NLST` rows and the non-parametric
  `NHFBNP` rows: `SFAC` scales `HYDCHR` per block, `OPEN/CLOSE` resolves against
  `model_ws` (quoted names supported), `EXTERNAL` resolves via `ext_unit_dict`
  (unresolvable → `NotImplementedError`, not `ValueError`). The writer still emits
  expanded inline rows (round-trip is semantically exact: `FACTOR*SFAC` baked in,
  `PARVAL` preserved). `write_file` now also requires `nacthfb ==
  len(acthfb_names)` (else `ValueError`). Tests: `-k mfusghfb` **15 passed** (7
  new: SFAC non-param/param/mixed, OPEN/CLOSE plain+quoted, EXTERNAL positive,
  EXTERNAL-without-dict fail, nacthfb-mismatch fail); SGB/QRT/DRT `_usgt_list`
  regressions green; focused **180**, exe **4**, combined **184** (ARM). See
  `USGT_STAGE4_04_PARAMETERS_HFB.md`.
- **Stage 4.4B polish — HFB parameter writer validation (executed):** the writer
  only guarded `parameters is None`, so a `NPHFB>0` header could still be written
  with no body from `parameters={}` or an inconsistent dict (a partial, invalid
  file). `write_file` now validates the preserved parameter state up front (new
  `_validate_parameter_write`, before the file is opened): from-scratch authoring
  (`parameters` `None`/empty) raises `NotImplementedError`; inconsistent
  definitions raise `ValueError` — `len(parameters) != NPHFB`, a def missing
  `partyp`/`parval`/`nlst`/`data`, `len(data) != nlst`, `nacthfb !=
  len(acthfb_names)`, or an active name not defined (case-insensitive). Files
  read by `load` satisfy these invariants, so valid round-trips are unaffected.
  Tests: `-k mfusghfb` **19 passed** (4 new negative tests, incl. an assertion
  that no partial file is written); focused **184**, exe **4**, combined **188**
  (ARM). See `USGT_STAGE4_04_PARAMETERS_HFB.md`.
- **Stage 4.4C — SGB list-parameter preservation (executed):** SGB `NPSGB>0`
  list parameters now load → write → reload with their syntax intact.
  `MfUsgSgb.load` consumes the leading `PARAMETER NPSGB MXS` record (`UPARLSTAL`),
  stores the per-parameter definitions (`UPARLSTRP`: name/partyp/parval/nlst +
  `NLST` `NODE GRADIENT [aux]` rows, 0-based) in `self.parameters` and the per-SP
  active-parameter names (`UPARLSTSUB`) in `self.active_params`, and keeps `MXS`;
  `write_file` re-emits all of it after validating the state up front. SFAC is
  inert on the gradient (Fortran ISCLOC=2 scales a dummy column), AUX is
  preserved, and `ITMP<0` reuse is kept (written expanded, matching the existing
  non-parametric behavior). The shared helper `_usgt_parameters.py` gains
  `read_list_parameter_count` / `write_list_parameter_count` (the `PARAMETER NP
  MXL` record); the existing list-parameter header/activation helpers are reused
  unchanged; SGB keeps its own `_read_sgb_rows`/`_write_sgb_rows`. Not promoted to
  `Full` on the parametric axis: parameter `INSTANCES` (Fortran-supported) and
  from-scratch parameter authoring raise `NotImplementedError`; inconsistent
  state raises `ValueError` (no partial file). Non-parametric SGB unchanged.
  Tests: `-k mfusgsgb` **14 passed** (9 new: round-trip, mixed non-param+active,
  reuse+active, AUX, SFAC-inert, OPEN/CLOSE, INSTANCES fail, from-scratch fail,
  inconsistent fail; the old parameters-fail test superseded); `-k "mfusgsgb or
  usgt_list"` **17**; focused **192**, exe **4**, combined **196** (ARM). DRT/QRT
  remain explicit-fail. See `USGT_STAGE4_04_PARAMETERS_SGB.md`.
- **Card 3 — DPT `A-W_ADSORBIM`:** decision is **explicitly unsupported**
  (deferred). Fortran audit of `dpt2aw_adsorb.f` (`AW_ADSORBIM1AL`) shows the
  option triggers a cascade of conditional arrays (zone map, tabular area
  functions, Langmuir isotherm arrays); too large and rare to model in v1.
  `MfUsgDpt.load` raises `NotImplementedError` at the option line before any
  extra read; test covers both the bare keyword and the `IAREA_FNIM IKAWI_FNIM`
  form. Spec recorded in roadmap Gap §6.
- **Card 4 — TIB:** decision is **keep raw/text round-trip for v1** (no semantic
  constructor). Fortran grammar (`GWF2TIB1RP`: `NIB0/NIB1/NIBM1
  [NICB0/NICB1/NICBM1]` + `U1DINT` node lists + node/head + transport blocks)
  recorded in the backlog. The raw round-trip preserves multi-node `U1DINT`
  continuation lines (tested) and round-trips Model A bit-for-bit; semantic
  authoring deferred until a target model needs it. → **Superseded by Stage 4.1**
  (semantic authoring + `parse=True` load now implemented; raw round-trip kept as
  the `load` default).
- **Card 5 — LAK → executed as Stage 4 LAK Fullness Card D; decision is still
  keep `✅` not Full.** From-scratch authoring is now real and tested for the main
  branches (no-transport, classic transport `CPPT`/`CRNF`, `TRANSPORTBOUNDARY`
  with MCOMP>1, `TABLEINPUT`), and six authoring bugs were fixed: `conc_data`
  mis-assignment (`{0: sill_data}` → `{0: conc_data}`); `write_file` crash when
  `flux_data` is None (now required); `conc_data` accessed when `mcomp==0`;
  `TRANSPORTBOUNDARY` dataset-9b written per-component instead of one
  `CLAKE(1:NSOL)` line per lake; load stored conc as strings (now float); and the
  `transportboundary` flag not being synced to the header keyword (so the file
  used the boundary layout without the keyword and would not reload). Added
  validation (flux_data required; TRANSPORTBOUNDARY needs transport; clake
  nlakes×mcomp; transport needs conc_data). `Ex8_Lake` real-model round-trip/run
  kept. Stays `✅` not Full — gaps: sill/connectivity (ds 7/8) + multi-lake
  systems round-trip but aren't authored from scratch; TABLEINPUT bathymetry
  table contents are external; GAGE coupling separate. See
  `USGT_STAGE4_LAK_FULLNESS.md`. **Review follow-up:** hardened three authoring
  inputs that `__init__` accepted but `write_file()` then crashed on with a raw
  `IndexError`/`KeyError`/`TypeError` — TABLEINPUT now requires exactly one
  `tab_file`/`tab_unit` per lake; `conc_data` is validated per written period
  (every `(lake, component)` present; classic 2 values for `WTHDRW>=0` / 3 for
  `WTHDRW<0`; TRANSPORTBOUNDARY a single value); `flux_data` requires one
  dataset-9a entry per lake — all clear `ValueError`s. `-k mfusglak` 10 passed
  (5 original + 5 follow-up).
- **Card 6 — base-class compatibility (SFR/STR/GAGE/FHB/SUB/SWT):** all six
  classified **Compatibility-only** in the roadmap (base MODFLOW-2005 classes;
  CLN is the project's coupling). Usage scan: only **FHB and GAGE** appear in a
  target/example model (the Ex8 lake model), where they load via the base class
  and round-trip is exercised; SFR/STR/SUB/SWT are used by no target model and
  are out of scope. No USG-T-specific semantics added; a dedicated per-package
  card would be opened only if a real model requires it.
- **Card 7 — QRT/DRT recipient `U1DINT` controls:** the Fortran `U1DINT`
  technically accepts `EXTERNAL`/`OPEN/CLOSE` for recipient-node lists, but
  these are short inline blocks in practice. Decision: support `INTERNAL` and
  `CONSTANT` (both tested — a CONSTANT spreading list expands to the repeated
  node); keep `EXTERNAL`/`OPEN/CLOSE` recipient lists as a documented
  `NotImplementedError` (rare-within-rare). Main-list controls
  (`SFAC`/`OPEN-CLOSE`/`EXTERNAL`) remain fully handled by `_usgt_list`.
- **Card 8 — plain-`✅` hardening:** promoted to `✅ Full` (authoring +
  round-trip + Fortran-verified): **CLN** (PROCESSCCF/ISHAPE/GENERAL_SEC),
  **DPF** (FRAHK/IUZONTABIM/SC2IM/immobile Richards), **DIS**/**DISU**
  (foundational; DISU large-grid `free_format_npl` formatting protected by
  `test_usg.py`), and **PCB** (new node/iSpec/conc authoring + round-trip
  test). Kept honest `✅` (not Full) with strengthened coverage/notes:
  **OC** (new `ATSA` authoring round-trip; `BOOTSTRAPPING` parsed on load),
  **EVT** (new transport `IETFACTOR`/`ETFACTOR` authoring test), and **MDT**
  (round-trip validated via the three Ex7 Matrix-Diffusion real models; field
  order not yet independently Fortran-audited).
- **Card 9 — executable end-to-end validation:** new opt-in suite
  `autotest/test_usg_transport_exe.py`, gated by `@requires_exe` on the
  `USGT_EXE` environment variable (default `mfusg_gsi`), so it skips cleanly
  when the executable is absent and is excluded from the default focused
  command. It authors models **from scratch**, runs them under USG-T 2.7, and
  checks outputs: a 1-D CHD flow model reproduces the analytical head gradient
  `[8, 6.5, 5, 3.5, 2]` with a closed `.list` budget (note: `linmeth=1`/PCGU —
  the default `linmeth=2`/XMD does not converge this trivial system), and a
  BCT+PCB transport model runs, emits a `.con`, and closes the species-isolated
  mass budget (`MfusgTransportListBudget`). The real-model `Ex1..Ex9` run tests
  (incl. Ex7 multi-species) share the **same** `USGT_EXE` contract (Stage 3
  follow-up), so the whole executable tier is selected by one env var.
- **Card 10 — upstream infrastructure track:** wrote `USGT_UPSTREAM_INFRA.md`.
  Verified the current reality (`get-modflow` installs `MODFLOW-ORG/executables`
  bundles; USG-Transport `mfusg_gsi` is a separate GSI build not in any
  FloPy-installable distribution; `@requires_exe` resolves a name on PATH or an
  absolute path). Proposes
  the path to satisfy the Nov-2024 criteria: a tagged USG-T source release with
  its own build CI + smoke tests, executable autotests living with that source,
  and only then an *optional* FloPy dependency via `USGT_EXE`/`get-modflow`.
  Execution is outside this fork.

All ten Stage 3 cards are complete; see `USGT_STAGE3_COMPLETION_PLAN.md`.

### Phase 2 hardening (2026-05-30)

Post-implementation critical review (`USGT_PHASE2_REVIEW.md`) found five
correctness/robustness gaps in the Priority-1 work; all are closed with
regression tests that fail before the fix and pass after. A follow-up
re-review (`USGT_PHASE2_REREVIEW.md`) added a polish pass — quote-aware
`OPEN/CLOSE` filenames (incl. single-quoted names with spaces) in
`_usgt_list`, positive `EXTERNAL` tests via `ext_unit_dict` for SGB/QRT/DRT,
a DRT omitted-recipient zero test, and review-doc cleanup. Suite:
`autotest/test_usg_transport.py` → **95 passed**.

| File | Fix |
|---|---|
| `flopy/mfusg/mfusgdrt.py` | `MfUsgDrt` is unstructured-only: structured construction raises `NotImplementedError` (use `ModflowDrt`); the dead structured `write_file` delegation was removed; structured `load` still delegates to `ModflowDrt.load`. (Previously structured authoring crashed with `AttributeError`.) |
| `flopy/mfusg/_usgt_list.py` (new) | Shared reader `begin_list_block` for the SGB/QRT/DRT main lists: consumes `SFAC` (scales Q for QRT, COND for DRT; inert on SGB gradient per Fortran ISCLOC), `OPEN/CLOSE`, and `EXTERNAL` (via `ext_unit_dict`) before row parsing; unresolved `EXTERNAL` raises `NotImplementedError`. Writes emit expanded inline rows (`Expanded valid write`). (Previously a leading `SFAC` raised a raw `ValueError`.) |
| `flopy/mfusg/mfusgbas.py` | `IHM` now takes an optional `IUIHM`: bare `IHM` or `IHM <option>` → `iuihm=0`; `IHM <int>` → that unit. (Previously a bare `IHM` raised `IndexError`.) |
| `flopy/mfusg/mfusgqrt.py`, `mfusgdrt.py` | `write_file` validates `recipient_nodes` against the stress list when `RETURNFLOW` is active: exactly one recipient list per record (omitted ⇒ all-zero); a length mismatch raises `ValueError` instead of silently dropping/shifting return-flow metadata. |
| `flopy/mfusg/_tabrich.py` | `node_count` falls back to DIS grid dimensions for `structured=False` models without a DISU package, and raises a clear `ValueError` if neither DISU nor DIS is present (was an opaque `AttributeError`). |

Status impact: SGB/QRT/DRT list-control *input* is classified
`Expanded valid write` (controls are read and expanded to inline rows on
output). The from-scratch authoring paths remain `Full semantic`.

### USG-T 2.7 Priority-1 packages (2026-05-30)

The four highest-value coverage gaps against the USG-T 2.7 Fortran source are
now closed (parser + writer + from-scratch authoring tests + round-trip tests,
all in `autotest/test_usg_transport.py`):

| Package | File | What it does |
|---|---|---|
| `MfUsgSgb` | `flopy/mfusg/mfusgsgb.py` | **New.** Specified Gradient Boundary (`glo2sgbu1.f`). Node-based `(node, gradient)` list, AUX transport concentrations, `ITMP/-1` reuse. Registered as `"sgb"`, so `MfUsg.load()` no longer silently skips SGB. `NPSGB>0` fails explicitly. |
| `MfUsgQrt` | `flopy/mfusg/mfusgqrt.py` | **New.** Sink with Return Flow (`gwf2QRT8u.f`). Per-sink `(node, q, rfprop)` plus variable-length recipient-node lists (`NodQRT` via `U1DINT`), `CHANGEC`/`IQCHNGTYP` transport, AUX, reuse. `AUTOFLOWREDUCE` preserved; `NPQRT>0` and `TRANSIENTQ` fail explicitly. Registered as `"qrt"`. |
| `MfUsgDrt` | `flopy/mfusg/mfusgdrt.py` | **New** (replaces base `ModflowDrt` in the registry). DRT8 (`gwf2drt8u.f`): EL+COND, `RETURNFLOW` single recipient (`NR>0`) or `SPREAD` multi-node (`NR<0`, `U1DINT` block), `CHANGEC`/`IDCHNGTYP` transport, AUX, reuse. `NPDRT>0` fails explicitly; structured grids delegate to base `ModflowDrt`. |
| `MfUsgBcf` / `MfUsgLpf` TABRICH | `flopy/mfusg/mfusgbcf.py`, `mfusglpf.py`, `_tabrich.py` | TABRICH items 1c (`IUZONTAB` zone map) and 1d (`RETCRVS`, shape `(nuzones, nutabrows, 3)` = capillary head / saturation / relative permeability) are now authored/loaded/written via a shared helper. For LPF the per-layer Richards arrays are skipped under TABRICH (matching `ITABRICH/=0`) and a token-index/`int` parse bug was fixed. Incomplete TABRICH writes fail explicitly. |

Shared helpers: `flopy/mfusg/_usgt_returnflow.py` (DRT/QRT recipient-node
`U1DINT` lists) and `flopy/mfusg/_tabrich.py` (BCF/LPF 1c/1d).

### Priority-2 closeout (2026-05-30)

| File | What changed |
|---|---|
| `flopy/mfusg/mfusgbas.py` | `RICHARDS_HP` (Richards mode with pressure-head initial values; implies effective Richards mode for BCF/LPF) and `IHM [IUIHM]` (integrated-hydrologic-model coupling flag + debug unit) are now authored/loaded/written. The option-line cleaner keeps `_` so `RICHARDS_HP` is one token. Gap §1 resolved. |
| `flopy/mfusg/mfusgdpt.py` | `MfUsgDpt.load` raises `NotImplementedError` for the immobile air-water adsorption option `A-W_ADSORBIM` instead of silently shifting subsequent reads. Gap §6. |
| `autotest/test_usg_transport.py` | Added authoring tests for BAS `RICHARDS_HP`/`IHM`, ETS `NETSOP=2` and `IESFACTOR`, explicit-failure tests for HFB `NPHFB>0` (load) and DPT `A-W_ADSORBIM` (load). ETS `npets>0`, HFB transient `IHFBRD=>0/0/-1`, and TIB multi-node `U1DINT` raw round-trip were already covered. |

ETS, HFB, and DPT remain classified Partial on purpose: parameter-syntax
preservation (ETS), parameterized barriers (HFB `NPHFB>0`), and the DPT
air-water adsorption sub-mode are intentionally out of scope and fail
explicitly rather than producing incomplete or mis-parsed files.

### Priority-3/4/5 review (2026-05-30)

- **BCT / DDF from-scratch authoring tests** added (previously these were only
  exercised through real-model round-trips): minimal BCT (1 species, IDISP=1),
  BCT IDISP=2 (full dispersion tensor), BCT multi-species (MCOMP=2), and a DDF
  NONLINEAR density-table round-trip. This justifies their `✅ Full` status.
- **CHD/RIV/GHB/DRN, WEL, CLN, DPF, TVM** already carry from-scratch authoring
  and/or round-trip tests; reviewed and left as-is.
- **Scope decisions (Priority 4):** `GSF` text round-trip + `to_grid()` is
  sufficient (tested). `LAK` now carries from-scratch `TABLEINPUT` /
  `TRANSPORTBOUNDARY` / classic-transport authoring + round-trip tests (Stage 4
  Card D; six bugs fixed) plus the `Ex8_Lake` round-trip, but stays `✅` not Full
  (sill/connectivity + multi-lake authoring, TABLEINPUT table contents, and GAGE
  remain out of scope). `SFR/STR/GAGE/FHB/SUB/SWT` are **compatibility-only**
  (base MODFLOW-2005 classes; USG-T unstructured records not validated; CLN is
  the preferred coupling) — documented, not silently "supported".
- **Post-processing (Priority 5):** `MfusgTransportListBudget` old/new format
  and multi-species isolation are tested. Real-model run-validation stays a
  manual, out-of-CI tier (see the Validation section below); the `Ex*`
  load+write round-trips are the in-CI real-model regression.

### Round-trip fixes (original five)

| Branch | What it does |
|---|---|
| `usgt/rch-transport-fix` | In `MfUsgRch.write_file`, `"# Stress period {kper + 1}"` becomes an f-string so the comment substitutes per SP, and `INRECH` is followed by `INIRCH` only when `NRCHOP == 2` (matching the read path). Both details were harmless at runtime but left the written RCH file not round-trippable. |
| `usgt/cln-load-none-unit` | In `MfUsgCln.__init__`, treat `None` entries in `unitnumber` as 0 instead of calling `int(None)`. This comes up when a CLN-declared output unit is not declared in the NAM's `ext_unit_dict`. |
| `usgt/tib-package` | Adds a new `MfUsgTib` class for the Transient Ibound package. **Stage 4.1: semantic** `stress_period_data` authoring (flow `NIB0/NIB1/NIBM1` + transport `NICB0/NICB1/NICBM1` blocks) and a `parse=True` loader, on top of the byte-exact raw round-trip that remains the `load` default (and the `parse=True` fallback for `EXTERNAL`/`OPEN-CLOSE` `U1DINT`). |
| `usgt/bas-preserve-unstructured` | `MfUsgBas.write_file` re-emits the `UNSTRUCTURED` keyword when `parent.structured is False`. The load side already reads the token; adding it to write closes the round-trip. |
| `usgt/nam-rebase-output-paths` | `BaseModel._reset_external` stores the basename of output files on `change_model_ws`, and `Modflow.write_name_file` preserves subdirectories for external input files. Previously the NAM writer also reduced inputs to basenames, breaking valid `DATA` paths in subfolders. |

### New packages and utilities

| What | File | What it does |
|---|---|---|
| `MfUsgChd` | `flopy/mfusg/mfusgchd.py` | CHD package for unstructured USG-T grids. Node-based (replaces k/i/j), supports AUX concentration variables. Internal `node` values are 0-based and file I/O is 1-based. Full `load` and `write_file` for the USG-T format (`NACT    Stress Period N` headers, `-1` reuse). |
| `MfUsgRiv` | `flopy/mfusg/mfusgriv.py` | RIV package for unstructured USG-T grids. Internal `node` values are 0-based and file I/O is 1-based. Supports AUX concentration and a trailing reach-ID column (`irch`) that is written positionally without being declared as AUX. `irch` is auto-detected from the first data row when loading. |
| `MfUsgEts` | `flopy/mfusg/mfusgets.py` | Segmented Evapotranspiration (ETS) package for USG-T 2.7. Supports `NETSEG > 1` with per-SP `PXDP`/`PETM` segment arrays and the `IESFACTOR` transport flag. Parameterized ETS files load by expanding parameters to concrete arrays and then write as valid non-parametric `NPETS=0`; preserving parameter syntax is not yet supported. |
| `MfUsgGhb` | `flopy/mfusg/mfusgghb.py` | GHB package for unstructured USG-T grids. Internal `node` values are 0-based and file I/O is 1-based. Supports AUX concentration variables, stores `ipakcb` (CBC unit). Full `load` and `write_file`. Load registry now maps `"ghb"` to `MfUsgGhb`. |
| `MfUsgDrn` | `flopy/mfusg/mfusgdrn.py` | DRN package for unstructured USG-T grids. Same pattern as GHB: internal `node` values are 0-based, file I/O is 1-based, AUX support, `ipakcb`. Full `load` and `write_file`. Load registry maps `"drn"` to `MfUsgDrn`. |
| `MfUsgTvm` | `flopy/mfusg/mfusgtvm.py` | TVM2 (Time-Variant Materials) package. Semantic implementation: global interpolation controls plus nper+1 stress-period boundary records, with 0-based internal nodes and 1-based file I/O. Missing boundaries emit all-zero headers on write. |
| `MfUsgGsf` | `flopy/mfusg/mfusgsf.py` | Grid Specification File wrapper. Text round-trip (stores raw lines). `to_grid()` delegates to `UnstructuredGrid.from_gridspec()` for full geometric parsing. Load registry maps `"gsf"` to `MfUsgGsf`. |
| `MfusgTransportListBudget` | `flopy/utils/mflistfile.py` | Reads transport species budget from a USG-T listing file for a single species. Handles both **old** USG-T format (transport blocks use `VOLUMETRIC BUDGET`, same keyword as flow) and **new** format (transport blocks use `MASS BUDGET`). Instantiate once per species: `MfusgTransportListBudget("model.lst", species=2)`. Returns the same recarrays / DataFrames as `MfusgListBudget`. |

Authoring note: semantic packages are expected to support direct construction
from Python/numpy inputs, not only `load()` + `write_file()` round-trips. The
USG-T boundary package tests now include from-scratch CHD/RIV/GHB/DRN creation
with 0-based internal nodes and 1-based file output.

### Authoring-focused fixes (2026-05-29)

| File | Fix |
|---|---|
| `flopy/mfusg/mfusgbas.py` | `MfUsgBas.load` now preserves the `UNSTRUCTURED` option on the constructed package. Programmatic tests cover write/load of `PRINTFV`, `CONVERGE`, `UNSTRUCTURED`, `FREE`, `PRINTTIME`, `SHOWPROGRESS`, `RICHARDS`, `DPIN`, `DPOUT`, `DPIO`, `SY-ALL`, and `STOPERROR`. |
| `flopy/mfusg/mfusgwel.py` | `options=None` is normalized before AUX auto-registration, so programmatic WEL authoring with AUX fields no longer crashes. Tests cover WEL rates assigned to GWF and CLN nodes, `ITMP NP ITMPCLN` headers, AUX concentrations, and 0-based internal node storage after reload. CLN connectivity stays in the CLN package. |
| `flopy/mfusg/mfusgcln.py`, `flopy/mfusg/cln_dtypes.py` | CLN now supports direct authoring and reload of `PROCESSCCF`/`ICLNGWCB`, `GENERAL_SEC`, and the 9-field `ISHAPE` node-property format used by rectangular/general conduit shapes. |
| `flopy/mfusg/mfusgdpf.py` | DPF now sets `model.idpf=1` when constructed programmatically, writes/loads optional `FRAHK`, supports `IUZONTABIM` for TABRICH models, writes `SC2IM` only for convertible layers, and supports immobile Richards arrays (`alphaIM`, `betaIM`, `srIM`, `brookIM`, optional `bPIM`). |

### Density-coupled round-trip fixes (2026-05-19)

Four bugs found during round-trip testing against a real-world Vistas-generated
density-coupled USG-T 2.7 model (BCT IDISP=2, DDF, 1382 stress periods):

| File | Fix |
|---|---|
| `flopy/mfusg/mfusgddf.py` | `ithickav` default: `default_val=1` → `default_val=0`. USG-T treats the absent ITHICKAV field as 0 (arithmetic averaging). Loading it as 1 changed transmissivity in the density layer and caused transport divergence. |
| `flopy/mfusg/mfusgchd.py` | `shead`/`ehead` dtype: `np.float32` → `np.float64`. Float32 round-trip loss (~0.04 mm per CHD node) accumulated over 469 density-coupled SPs into concentration differences that marginally failed outer-loop convergence (ΔC > CICLOSE=1e-8). |
| `flopy/mfusg/mfusgriv.py` | `stage` dtype: `np.float32` → `np.float64`, same reason. |
| `flopy/mfusg/mfusgbas.py` | `CONVERGE` option: (1) `converge=converge` was missing from the `cls(...)` call in `load`, so `self.converge` was always `False`; (2) `write_file` never emitted `CONVERGE` even when `self.converge=True`. The `CONVERGE` keyword in BAS6 tells USG-T to use the coupled flow–transport convergence criterion in the outer nonlinear loop — without it the model failed to accept time step 2 of the first pumping stress period after 250 iterations. |

### Minor fixes

| File | Fix |
|---|---|
| `flopy/mfusg/mfusgbct.py` | Removed a stray `print()` left from development. Fixed file-handle management in `write_file`: the file is now closed only when opened internally (`close_on_exit` flag), so callers that pass an open handle are not surprised. `ICBUND` now loads as `np.int32`, so round-trip writes integer constants instead of `1.000000E+00` values that USG-T rejects. |

### Fortran source audit, new implementations, and bug fixes (2026-05-20)

All packages cross-referenced against the USG-Transport 2.7.0 Fortran source code
(`USGT_V_2-7-0_Source_Code/`). See `USGT_roadmap.md` for the full coverage table.

#### New: `MfUsgHfb` — HFB with TRANSIENT_HFB support

| File | What it does |
|---|---|
| `flopy/mfusg/mfusghfb.py` | HFB6 (Hydraulic Flow Barrier) implementation for unstructured USG-T grids. Extends `ModflowHfb` with: (a) node-based `(node1, node2, hydchr)` unstructured format, (b) `TRANSIENT_HFB` keyword and the Fortran `IHFBRD` reuse/read flag semantics. Non-parametric static and transient HFB load/write are semantic; parameterized HFB (`NPHFB > 0`) now fails explicitly until preservation/expansion is implemented. |

#### New: `MfUsgTvm` — TVM2 full semantic implementation

`mfusgtvm.py` was rewritten from a verbatim text round-tripper into a complete
semantic implementation verified against `tvmu2.f`:

- Parameters: `itvmprint`, `tvmlogbasehk/vka/ss/sy`, `tvmddftr`, `tvmlogbasepor`
- `stress_period_data`: `dict[int, dict[str, np.recarray]]` — boundary index 0..nper,
  properties `hk/vka/ss/sy/ddftr/por`, nodes 0-based
- Transport-aware: with BCT → 7 global + 6 SP fields; without BCT → 6 + 5 fields
- nper+1 boundary blocks confirmed against `TVMU2AR` + `TVMU2RP` call structure
- Fortran fixed-format `(I10,F10.0)` per record confirmed
- 20 autotests in `autotest_local/test_mfusg_tvm.py` (all pass)

#### LPF Richards fix

| File | Fix |
|---|---|
| `flopy/mfusg/mfusglpf.py` | `Util2d.__eq__` returns `False` for any non-`Util2d` argument, so `if self.laytyp == 5:` was always `False`. This silently disabled Richards-equation array initialization (ALPHA/BETA/SR/BROOK) for any model with LAYTYP=5, even though the load path read them correctly. Fixed: replaced with `np.any(self.laytyp.array == 5)`. |

#### Bugs found in Fortran source audit — all fixed

| File | Bug | Severity | Fortran reference |
|---|---|---|---|
| `mfusgrch.py:396` | `t.index("INIZNRCH")` — `INIZNRCH` is an internal Fortran variable, not the file keyword. Raises `ValueError` on any model using RTS recharge zones. | **Critical** | `gwf2rch8u1.f` keyword is `INRCHZONES` |
| `mfusgdpt.py:355,507` | DLIM written/loaded when only `idpf` is True. Fortran requires `IDPF≠0 AND IDISPIM≠0`. Models with `idpf=1, idispim=0` wrote an extra array, shifting all subsequent reads. | High | `gwt2dptu1.f` line 281: `IF(IDPF.NE.0.AND.IDISPIM.NE.0)THEN` |
| `mfusgddf.py` | ISHARP (sharp-interface model flag) missing entirely from `__init__`, `write_file`, and `load`. ISHARP is the 6th numeric field after IMPHDD. Loading any file with ISHARP≠0 silently lost the flag; round-trip discarded it. | High | `density.f` line 64: `CALL URWORD(...,ISHARP,...)` |
| `mfusgdpf.py:191` | `f_obj` unbound `NameError` when `write_file(f=<open handle>)` called directly. The `else: f_obj = f` branch was missing. | Low | — |

## Install

```bash
conda create -n usgt -y python=3.12 numpy pandas matplotlib jupyter shapely
conda activate usgt
git clone --branch develop https://github.com/roh-sgm/flopy.git
cd flopy
pip install -e .
```

## Quick check

```python
import flopy
print(flopy.__version__)                              # 3.11.0.dev0 (or newer)
print('MfUsgTib' in dir(flopy.mfusg))                # True
print('MfUsgChd' in dir(flopy.mfusg))                # True
print('MfUsgRiv' in dir(flopy.mfusg))                # True
print('MfUsgGhb' in dir(flopy.mfusg))                # True
print('MfUsgDrn' in dir(flopy.mfusg))                # True
print('MfUsgTvm' in dir(flopy.mfusg))                # True
print('MfUsgGsf' in dir(flopy.mfusg))                # True
print('MfUsgEts' in dir(flopy.mfusg))                # True
print('MfusgTransportListBudget' in dir(flopy.utils)) # True

# Multi-species transport budget reader
from flopy.utils import MfusgTransportListBudget
lst_s1 = MfusgTransportListBudget("model.lst", species=1)
lst_s2 = MfusgTransportListBudget("model.lst", species=2)
inc, cum = lst_s1.get_budget()
df_inc, df_cum = lst_s2.get_dataframes(start_datetime="2000-01-01")
```

All other flopy APIs behave exactly as upstream; the patches only change the
text that `MfUsg.load(...).write_input()` produces for USG-T models.

## Validation

Exercised end-to-end against two real-world unstructured USG-T models built
with Groundwater Vistas, using their Vistas-generated input packages as the
reference:

**Model A — BCT + CLN + TIB transport** (≈ 112 k nodes, 493 stress periods, USG-T 1.8):
FloPy-rewritten input produces a listing file matching the reference bit-for-bit on every
budget column — flow and transport (`max|diff| = 0`, 20 flow + 16 mass components).
Wall-clock runtime equivalent to the reference.

**Model B — BCT IDISP=2 + DDF density-coupled transport** (≈ 19 k nodes, 1382 stress periods,
4322 time steps, USG-T 2.7 ARM):

| Output | Time steps | Max difference |
|--------|-----------|----------------|
| LST budget | 4322 / 4322 ✓ | 0.000e+00 |
| HDS heads | 4322 / 4322 ✓ | 0.000e+00 m |
| CON concentrations | 4322 / 4322 ✓ | 0.000e+00 |
| CBB (7 record types) | 4322 / 4322 ✓ | 0.000e+00 |

**Both models: bit-for-bit identical to the Vistas reference on all outputs.**

Tested binaries:
- USG-T 2.7.0 (ARM) — Model B full validation
- USG-T 1.8 (ARM and x86) — Model A full validation

## Upstream

Held on the fork while testing continues. No PRs to `modflowpy/flopy` have
been opened yet.
