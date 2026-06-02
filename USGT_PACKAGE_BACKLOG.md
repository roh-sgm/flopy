# USG-T Package Backlog

Purpose: turn the USG-Transport support work into executable package-level
tasks for agents, with an explicit bias toward **programmatic authoring** as
well as round-trip preservation.

This file is a working backlog, not a proof of completion. A package is only
done when an agent has verified the Fortran I/O, implemented semantic load/write
where appropriate, added from-scratch authoring tests, added round-trip tests,
and updated the coverage docs.

Primary source of truth:

`/Users/roh.sgm/Documents/SGM Co/GDM/01_Software/GSI/USG-Transport 2.7.0 (2026.03.21)/USGT_V_2-7-0_Source_Code`

Use the USG-T 2.7 I/O manual as secondary reference. If manual and Fortran
disagree, follow the Fortran.

---

## Closeout status (2026-05-30)

All priority tiers below have been worked through. Commits on `develop`:
`b0de152f` (checkpoint of prior P0/P2/P3 work), `ad7de360` (P1 packages),
`036d714d` (P2 closeout), `6429dcc6` (P3/P4/P5 closeout).

- **P0 — guardrails:** accepted (NAM external paths; CHD/RIV/GHB/DRN/WEL 0-based
  nodes), tested.
- **P1 — highest-value gaps:** **done.** `MfUsgSgb`, `MfUsgQrt`, `MfUsgDrt`
  (new/registered), and BCF/LPF `TABRICH` 1c/1d. Each: parser + writer +
  authoring + round-trip tests, verified against Fortran.
- **P2 — partial packages:** ETS (NETSEG/NETSOP/IESFACTOR authoring; params
  expand or fail explicitly), HFB (static/transient/IHFBRD; `NPHFB>0` fails),
  BAS (`RICHARDS_HP`, `IHM` implemented), DPT (`A-W_ADSORBIM` explicit fail),
  TIB (raw/text in P2; promoted to **Full (authoring)** in Stage 4.1).
- **P3 — review:** BCT (1-species/IDISP=2/multi-species) and DDF (NONLINEAR
  table) from-scratch authoring tests added; other Full packages reviewed.
- **P4 — base-class packages:** GSF text round-trip in P4, promoted to **Full
  (authoring)** in Stage 4.2; LAK from-scratch authoring + round-trip tested in
  Stage 4 Card D (six bugs fixed; kept `✅` not Full, gaps listed) plus Ex8;
  SFR/STR/GAGE/FHB/SUB/SWT compatibility-only (documented).
- **P5 — post-processing:** transport list-budget tested (old/new/multi-species);
  real-model run-validation is a documented manual tier (Ex* loads in CI).

Test status (current, after Stage 3 + follow-up): `autotest/test_usg_transport.py`
→ **101 passed** when a USG-T executable resolves (`Ex1..Ex9` run), or
**85 passed, 16 skipped** with no executable (the `Ex*` run tests skip). With
the opt-in `autotest/test_usg_transport_exe.py`, both suites together →
**103 passed** under `USGT_EXE` (USG-T 2.7), or **85 passed, 18 skipped** with
no resolvable executable. (Historical counts below: 80 after the initial
backlog pass, 95 after Phase 2 + polish.)

## Critical review addendum (2026-05-30)

Follow-up review found that several closeout claims needed a Phase 2 hardening
pass before they should be treated as upstream-ready. See
`USGT_PHASE2_REVIEW.md` for the detailed reviewer brief and the
per-finding resolution table.

**Phase 2 status: all five reopened items are RESOLVED**, plus a follow-up
polish pass (`USGT_PHASE2_REREVIEW.md`: quoted `OPEN/CLOSE`, positive
`EXTERNAL` tests, DRT zero-recipient test, doc cleanup). On `develop`;
`python -m pytest autotest/test_usg_transport.py -q` → **95 passed**.
Summary of closures:

- **P0 - DRT structured authoring:** RESOLVED. `MfUsgDrt` is now unstructured
  USG-T DRT8 only; structured construction raises `NotImplementedError` (use
  `ModflowDrt`), `load` delegates structured files to the base class.
- **P0 - SGB/QRT/DRT list controls:** RESOLVED (Expanded valid write). Shared
  `_usgt_list.begin_list_block` handles `SFAC`, `OPEN/CLOSE`, and `EXTERNAL`
  (via `ext_unit_dict`) on the main lists; unresolved `EXTERNAL` raises
  `NotImplementedError`; writes emit expanded inline rows.
- **P1 - BAS `IHM`:** RESOLVED. `IHM` parses an optional `IUIHM`; bare `IHM`
  (or `IHM` before another option) yields `iuihm=0`. No `IndexError`.
- **P1 - QRT/DRT recipient validation:** RESOLVED. `write_file` requires one
  `recipient_nodes` list per record when `RETURNFLOW` is active (omitted ⇒
  all-zero); mismatches raise `ValueError`.
- **P2 - TABRICH node count contract:** RESOLVED. `_tabrich.node_count` falls
  back to DIS grid dimensions for `structured=False` models without DISU and
  raises a clear `ValueError` if neither is present.

Resulting status labels: SGB/QRT/DRT list-control input is **Expanded valid
write** (controls read, expanded inline on output). The from-scratch authoring
paths (inline rows, recipients, TABRICH) remain `Full semantic` where parser +
writer + authoring + round-trip tests all exist.

---

## Non-Negotiable Design Rules

1. Internal FloPy data is 0-based. USG-T files are 1-based.
2. Every semantic package must be buildable from Python/numpy inputs without
   first loading an existing file.
3. Round-trip tests and programmatic authoring tests are separate test cases.
4. Do not mark a package `Full semantic` unless it has parser, writer,
   authoring tests, and round-trip tests.
5. If MODFLOW parameter syntax is read but not preserved, either write expanded
   non-parametric arrays intentionally or fail explicitly. Do not silently write
   incomplete parameter syntax.
6. WEL applies rates to GWF and CLN nodes. CLN owns CLN connectivity, geometry,
   CLN-CLN links, and CLN-GWF links.
7. NAM/external-file behavior is global risk. Preserve input subdirectories;
   rebase outputs only when that is the intended FloPy behavior.
8. Prefer small synthetic tests in default CI. Real models belong in a slow or
   optional validation suite.

## Package Status Vocabulary

- `Full semantic`: semantic constructor, load, write, authoring tests, round-trip
  tests, verified against Fortran.
- `Expanded valid write`: can read a richer syntax but writes a valid expanded
  form instead of preserving original syntax.
- `Raw/text round-trip`: preserves existing text but does not claim a semantic
  API.
- `Patch candidate`: code may exist in this branch, but still needs critical
  review before being treated as accepted design.

## Agent Gate

Every agent should complete these gates in order:

1. Write a short spec from the Fortran source: datasets, optional flags,
   conditional reads, reuse rules, indexing, parameters, and external arrays.
2. Compare the current FloPy implementation to that spec.
3. Decide the package class: `Full semantic`, `Expanded valid write`, or
   `Raw/text round-trip`.
4. Implement only the scoped package or cross-cutting behavior assigned.
5. Add one authoring test and one round-trip test at minimum.
6. Run targeted tests and `autotest/test_usg_transport.py` if the package is in
   that suite.
7. Update `USGT_roadmap.md`, `USGT_improvements.md`, and this backlog if the
   status changes.

---

## Priority 0 - Cross-Cutting Guardrails

### NAM and external files

Status: **Accepted** — committed and covered by
`test_modflow_name_file_preserves_input_external_paths`; input `DATA`
subdirectories are preserved and outputs rebased to basenames. (Checkpoint
`b0de152f`.)

Why it matters: upstream-ready package support can still fail if name-file
rewriting breaks external input paths.

Tasks:

- Review `Modflow.write_name_file` and `_reset_external` behavior.
- Preserve input `DATA` and `DATA(BINARY)` subdirectories.
- Rebase output files to basenames only where outputs are intentionally moved.
- Confirm binary and text external input cases.
- Confirm model workspace changes do not strand external arrays.

Required tests:

- Input external arrays in subdirectories remain in NAM as subdirectory paths.
- External output files are rebased as output basenames.
- Both text and binary `DATA` entries are covered.

Acceptance:

- No global behavior change outside NAM/external handling without regression
  tests and docs.

### Shared USG-T list packages

Status: **Accepted** — CHD/RIV/GHB/DRN/WEL keep `node` 0-based internally and
write 1-based; AUX concentrations and `-1` reuse are covered by authoring and
round-trip tests. (Checkpoint `b0de152f`; SGB/QRT/DRT follow the same pattern.)

Why it matters: many USG-T packages use node-based records. One indexing error
silently shifts the physical target cell.

Tasks:

- Keep `node` internal values 0-based.
- Write file nodes 1-based.
- Confirm `-1` stress-period reuse remains compatible with `MfList`.
- Keep AUX concentration fields round-trippable and authorable.

Required tests:

- Programmatic creation with 0-based nodes.
- Written file inspection for 1-based nodes.
- Reload returns 0-based nodes.
- AUX fields survive write/load.

Acceptance:

- The package API never asks users to provide 1-based node ids.

---

## Priority 1 - Highest Value Gaps

### BCF6 / LPF - TABRICH

FloPy classes: `MfUsgBcf`, `MfUsgLpf`

Fortran: `gwf2bcf-lpf-u1.f`

Status: **DONE** (2026-05-30). `IUZONTAB` + `RETCRVS` authored/loaded/written
for BCF and LPF via shared `_tabrich.py`; LPF Richards arrays skipped under
TABRICH; incomplete TABRICH write fails explicitly. Authoring + round-trip +
shape-validation tests in `autotest/test_usg_transport.py`.

Problem:

- Current support recognizes `TABRICH`, `NUZONES`, and `NUTABROWS`, but does not
  semantically read/write all TABRICH datasets.
- Missing datasets are `IUZONTAB` and `RETCRVS`.
- A model authored with `tabrich=True` can produce an incomplete BCF/LPF file.

Fortran facts to verify:

- Read location of `IUZONTAB`.
- Shape/order of `RETCRVS`.
- Differences between BCF and LPF paths, if any.
- Interaction with `LAYCON=5` / `LAYTYP=5`.
- Interaction with `BUBBLEPT`, `FULLYDRY`, and `ALTSTO`.

Implementation tasks:

- Add constructor args for `iuzontab` and `retcrvs`.
- Choose a stable Python representation for `RETCRVS`, likely ndarray with
  explicit axes `(nuzones, nutabrows, 3)`.
- Implement write order exactly as Fortran reads it.
- Implement load order exactly as Fortran writes it.
- Decide whether BCF and LPF should share helpers.
- Validate shapes at construction time.

Required tests:

- BCF authoring with `TABRICH`, `IUZONTAB`, and `RETCRVS`.
- LPF authoring with `TABRICH`, `IUZONTAB`, and `RETCRVS`.
- Round-trip BCF TABRICH file.
- Round-trip LPF TABRICH file.
- Negative shape validation for malformed `RETCRVS`.

Acceptance:

- A new BCF/LPF TABRICH model can be authored from scratch and written as a
  complete USG-T input file.

### DRT - USG-T 2.7 extensions

FloPy class: currently base `ModflowDrt`

Fortran: `gwf2drt8u.f`

Status: **DONE** (2026-05-30). `MfUsgDrt` (subclass of `ModflowDrt`) added and
registered. Node-based EL+COND, RETURNFLOW single/spread recipients,
`CHANGEC`/`IDCHNGTYP`, AUX, reuse; structured delegates to base. Authoring +
round-trip + reuse tests added.

**DONE** (Stage 4.4D): `NPDRT>0` list parameters are now **preserved**
(load → write → reload, including activations and per-row recipients). Fortran
audit confirmed DRT is internally consistent — definition and activation both
use `PARTYP='DRT'` (`gwf2drt8u.f:135` / `:1108`), unlike SGB's `'SGB'`/`'G'`
mismatch — so active DRT parameters are valid. `load` keeps `NPDRT`/`MXL`
(item 1) and the per-parameter definitions (`UPARLSTRP` header + `NLST` drain
rows + RETURNFLOW recipients/spreading via the shared `_read_drain_rows`) in
`self.parameters`, and the per-SP active names (`SGWF2DRT8LS`) in
`self.active_params`; the parameter value scales `COND` (`IPVL=5`). `write_file`
re-emits all of it after validating up front (no partial file). Not `Full` on
the parametric axis: from-scratch authoring and `INSTANCES` raise
`NotImplementedError`; `MXL`/consistency → `ValueError`. See
`USGT_STAGE4_04_PARAMETERS_DRT.md`.

**DONE** (Stage 4.4D review follow-up): (P1) `MXADRT` now sized to the active
total per period (`non-parametric + Σ active-parameter NLST`); the writer used
the non-parametric count only, which would abort the Fortran (`NDRTCL >
MXADRT`). (P2) activated SPREAD (`NR<0`) parameter recipients documented as
structural-round-trip only, not execution-guaranteed (USG-T copies `DRTF` but
not `NodDRT` on activation). (P3) roadmap Gap §7 de-staled. 11 parameter tests
(`-k mfusgdrt` 22 passed), incl. a multi-active-parameter `MXADRT=4` case.

Problem:

- USG-T DRT8 extends base DRT with transport and flow-reduction controls.
- The base FloPy class does not expose all USG-T fields.

Fortran facts to verify:

- Header items and package options.
- Stress-period record layout.
- AUX concentration behavior.
- `IQCHANGEC` / `IQCHNGTYP` semantics.
- `MXSPREADNDS` and spreading return-flow records.
- Node indexing for unstructured grids.
- Reuse behavior for stress periods.

Implementation tasks:

- Decide whether to subclass base DRT or create `MfUsgDrt`.
- Add USG-T load/write registry mapping if new class is used.
- Add dtype definitions for all stress-period variants.
- Keep source and return nodes 0-based internally.
- Write all file nodes 1-based.
- Add explicit unsupported-feature failures if a sub-mode is not implemented in
  v1.

Required tests:

- Programmatic DRT with transport AUX concentrations.
- Programmatic DRT with `IQCHANGEC` / `IQCHNGTYP`.
- Programmatic DRT with `MXSPREADNDS`.
- Round-trip minimal USG-T DRT file.
- Stress-period reuse.

Acceptance:

- A DRT package using USG-T transport extensions can be authored without hand
  editing the package file.

### SGB - Specified Gradient Boundary

FloPy class: missing

Fortran: `glo2sgbu1.f`

Status: implemented; non-parametric Full (authoring) + parameter-**definition**-
preserving (`NPSGB`, no activations) as of Stage 4.4C + review follow-up. Active
SGB parameters are unsupported (USG-T 2.7 Fortran type conflict).

Problem:

- **DONE** (2026-05-30): `MfUsgSgb` added (`flopy/mfusg/mfusgsgb.py`) and
  registered as `"sgb"`. Node-based `(node, gradient)`, AUX, `ITMP/-1` reuse,
  0-based internal / 1-based file, `NPSGB>0` explicit failure. Authoring +
  round-trip + NAM-registry + parameter-failure tests added.
- **DONE** (Stage 4.4C): `NPSGB>0` parameter **definitions** are preserved
  (load → write → reload). `load` keeps the `PARAMETER NPSGB MXS` record
  (`UPARLSTAL`) and the per-parameter definitions (`UPARLSTRP`:
  name/partyp/parval/nlst + `NLST` `NODE GRADIENT [aux]` rows, 0-based) in
  `self.parameters`; `MXS` preserved. SFAC inert on the gradient (Fortran
  ISCLOC=2). Uses the shared list-parameter helpers in `_usgt_parameters.py`.
- **DONE** (Stage 4.4C review follow-up): **active SGB parameters are
  unsupported** and now fail explicitly. Re-audit showed USG-T 2.7 defines SGB
  parameters as `PARTYP='SGB'` (`UPARLSTRP`, glo2sgbu1.f:97) but activates them
  as `PTYP='G'` (`UPARLSTSUB`, glo2sgbu1.f:185); a parameter has one type, so any
  activation trips "Parameter type conflict" (`parutl7.f:684/800`) and aborts the
  run. So per-SP `NP>0` and `INSTANCES` raise `NotImplementedError` (load/write).
  Writer hardened: `MXS<=0` or `MXS <` total definition rows → `ValueError`;
  inconsistent definitions → `ValueError`; no partial file. SGB is
  **definition-preserving only**. 10 parameter tests (`-k mfusgsgb` 15 passed).
  See `USGT_STAGE4_04_PARAMETERS_SGB.md`.

Fortran facts to verify:

- NAM/CUNIT package key.
- Header items.
- Stress-period block structure.
- Whether data is node-based, face-based, or geometry-dependent.
- Reuse semantics.
- Any transport concentration coupling.
- Output units, if any.

Implementation tasks:

- Create `MfUsgSgb`.
- Add registry mapping in `MfUsg`.
- Define semantic dtype(s).
- Implement constructor, load, write.
- Decide whether v1 supports all options or fails explicitly for rare modes.

Required tests:

- Programmatic minimal SGB package.
- Round-trip minimal SGB file.
- NAM load registry test.
- 0-based internal / 1-based file indexing, if node-based.

Acceptance:

- `MfUsg.load()` no longer skips SGB models, and a minimal SGB can be authored
  from scratch.

### QRT - Sink with Return Flow

FloPy class: missing

Fortran: `gwf2QRT8u.f`

Status: implemented; non-parametric Full (authoring) + parameter
structural-preserving (`NPQRT`, not execution-guaranteed) as of Stage 4.4E.

Problem:

- **DONE** (2026-05-30): `MfUsgQrt` added (`flopy/mfusg/mfusgqrt.py`) and
  registered as `"qrt"`. Per-sink `(node, q, rfprop)` + recipient-node lists,
  `CHANGEC`/`IQCHNGTYP`, AUX, reuse; `TRANSIENTQ` explicit failure (lifted in
  Stage 4.5A, below); shares `_usgt_returnflow.py` with DRT. Authoring (minimal,
  return-flow concentration, multi-recipient) + round-trip + registry tests.
- **DONE** (Stage 4.4E): `NPQRT>0` list parameters are **structurally preserved**
  (load → write → reload of definitions + recipient blocks + per-SP activations;
  `MXAQRT` = non-parametric + active rows; `MXRTCELLS` reflects definition
  recipients). Audit: QRT is type-consistent (def + activation both
  `PARTYP='QRT'`, `gwf2QRT8u.f:183`/`:1118`), so active params are type-valid —
  **but not execution-guaranteed**: the parameter value scales `QRTF(5)=NumRT`,
  not `Q` (`IPVL=5`; a Fortran bug, `SFAC` scales `Q` at `ISCLOC=4`), and
  `NodQRT` is not copied on activation (like DRT). FloPy does not apply the
  parameter value to `Q`. From-scratch parameter authoring and `INSTANCES` →
  `NotImplementedError`; `MXL`/consistency → `ValueError`. Reuses the shared
  list-parameter helpers. See `USGT_STAGE4_04_PARAMETERS_QRT.md`. **This
  completes the Stage 4.4 list-parameter family (HFB/SGB/DRT/QRT).**
- **DONE** (Stage 4.4E review follow-up, applies to QRT and DRT): active
  parameter names are resolved **case-insensitively** when sizing `MXAQRT`/
  `MXADRT` (shared `resolve_list_parameter` helper), matching the validation and
  the Fortran `UPCASE`; previously a `qp`/`QP` (or `dp`/`DP`) case mismatch
  dropped the active `NLST` from the header. Duplicate activation of a parameter
  in one stress period now raises `ValueError`. 4 new tests; `-k mfusgqrt` 20,
  `-k mfusgdrt` 24 passed.
- **DONE** (Stage 4.5A): the inline **`TRANSIENTQ`** transient extraction-flow
  time series is now supported (replacing the prior explicit failure). Fortran
  audit (`gwf2QRT8u.f`): the block is read in `GWF2QRT8U1AR` after the parameter
  definitions and before stress periods (times control + `NBDQTIM` times, values
  control + exactly `MXAQRT` rows of `(node, NBDQTIM values)`); `GWF2QRT8U1AD`
  overrides `QRTF(4)=Q` only (recipients untouched); `NBDQTIM<0` selects
  staircasing; `IQRTN` is an informational node tag (series applied
  positionally). FloPy adds `transientq_times`/`values`/`nodes`/`staircase` +
  `CNSTM` multipliers, preserved raw on load → write → reload and authorable from
  scratch; 0-based internal / 1-based file. Explicit failures: `TRANSIENTQ` +
  `NPQRT>0` (the Fortran reads `BDQV` past its `MXAQRT` allocation when `MXL>0`),
  external-unit data, dimension mismatch. 6 new tests; `-k mfusgqrt` 26 passed.
  See `USGT_STAGE4_05_QRT_TRANSIENTQ.md`.
- **DONE** (Stage 4.5A review follow-up): closed two TRANSIENTQ contract gaps.
  (1) **External-unit data now fails explicitly** — `_read_transientq_block`
  validates each control line's `IQRTUN` against the QRT package's inline unit
  (resolved from the NAM via `get_ext_dict_attr(..., pop_key=False)`, else
  `_defaultunit()`, matching `write_file`'s `self.unit_number[0]`); a non-inline
  unit raises `NotImplementedError` before any data is read as inline. (2)
  **`TRANSIENTQ` must be the last item-1 option** — `_parse_header` rejects any
  trailing token with `ValueError` (the Fortran branch does not loop back). Also
  added non-negative 0-based integer validation for `transientq_nodes` on write.
  5 new tests (external times/values, trailing option, negative node, non-unit
  `CNSTM` round-trip); `-k mfusgqrt` **31 passed**, focused **229**, combined
  **233** (ARM). Only `mfusgqrt.py` touched in code.

Fortran facts to verify:

- NAM/CUNIT package key and expected file extension.
- Header items.
- Stress-period record layout.
- Transport concentration fields.
- `IQCHANGEC` / `IQCHNGTYP` behavior.
- `MXSPREADNDS` behavior.
- Reuse rules.

Implementation tasks:

- Create `MfUsgQrt`.
- Add registry mapping in `MfUsg`.
- Decide whether to share helpers with `MfUsgDrt`.
- Define semantic dtype(s) for extraction, return, transport, and spread nodes.
- Implement constructor, load, write.

Required tests:

- Programmatic minimal QRT.
- Programmatic QRT with return-flow concentration.
- Programmatic QRT with spread nodes.
- Round-trip minimal QRT.
- NAM load registry test.

Acceptance:

- A QRT-using USG-T model can be loaded/written without silently dropping the
  package, and minimal QRT can be authored from Python.

---

## Priority 2 - Important Partial Packages

### ETS - parameters

FloPy class: `MfUsgEts`

Fortran: `gwf2ets8u1.f`

Status: **parameter-preserving for the ETSR array parameter (Stage 4.4A,
executed)**; semantic for non-parametric ETS. **Authoring tests added**
(2026-05-30): NETSEG=1, NETSEG>1, NETSOP=2, IESFACTOR, and explicit `npets>0`
write failure.

Resolution (Stage 4.4A): `MfUsgEts.load` now *preserves* ETSR array parameters
by default — it keeps the parsed definitions (`self.parameters`, a
`ModflowParBc`) and the per-period activation records (`self.evtr_parm`), and
`write_file` re-emits `NPETS>0` in item 2a, the definition blocks, and the
activation records (with `INSTANCES`), while ETSS/ETSX/IETS/PXDP/PETM stay plain
arrays. USG-T reads `NPETS` from item 2a (`UPARARRAL` is called with `IN=-1`),
so there is no `PARAMETER` line. The shared write helper is
`flopy/mfusg/_usgt_parameters.py` (reuses `ModflowParBc`; no second parser).
Opt-in `expand_parameters=True` keeps the legacy expanded path (`NPETS=0`).
Still not `Full`: from-scratch parameter *authoring* (`npets>0` without loaded
defs) raises `NotImplementedError`. See `USGT_STAGE4_04_PARAMETERS_ETS.md`.

Required tests (all green, `-k mfusgets` 10 passed):

- `NETSEG=1` authoring. (Done.)
- `NETSEG>1` authoring. (Done.)
- `NETSOP=2` authoring. (Done.)
- `IESFACTOR` authoring. (Done.)
- Parameterized load **preserves** (`npets`, defs, per-SP records). (Done.)
- Write **preserves** parameter syntax + mixes plain ETSS/ETSX arrays. (Done.)
- Reload of the preserved file round-trips. (Done.)
- Time-varying parameter (`INSTANCES`) round-trips. (Done.)
- Expanded fallback (`expand_parameters=True`) writes `NPETS=0`. (Done —
  regression.)
- Programmatic `npets>0` (no defs) fails explicitly. (Done.)

Acceptance:

- Users can author non-parametric ETS from scratch. (Yes.)
- Loaded ETSR parameters are preserved on write/reload; from-scratch parameter
  authoring is an explicit, documented `NotImplementedError`. (Yes.)

### HFB - parameterized barriers

FloPy class: `MfUsgHfb`

Fortran: `gwf2hfb7u1.f`

Status: **parameter-preserving for HFB list parameters (Stage 4.4B, executed)**;
semantic for non-parametric static/transient HFB. **Tests cover**: static,
structured static, transient `IHFBRD=>0/0/-1`, parameterized load → write →
reload (unstructured + structured), parameters mixed with non-parametric
barriers, and the explicit-failure cases.

Resolution (Stage 4.4B): HFB uses MODFLOW **list** parameters
(`UPARLSTRP`/`UPARLSTSUB`) — each parameter owns `NLST` barrier rows and its value
scales `HYDCHR`. `MfUsgHfb.load` now *preserves* them: it stores the definitions
(`self.parameters = {name: {"partyp","parval","nlst","data"}}`, barrier rows
0-based), the non-parametric barriers (`self.hfb_data`), and the active-parameter
names (`self.acthfb_names`); `write_file` re-emits the definition blocks, the
non-parametric barriers, and `NACTHFB` + the active names. Shared list-parameter
helpers live in `flopy/mfusg/_usgt_parameters.py`. See
`USGT_STAGE4_04_PARAMETERS_HFB.md`.

Not `Full`: from-scratch parameter authoring (`NPHFB>0` without loaded defs),
`TRANSIENT_HFB`+`NPHFB>0` (the Fortran redefines params each SP under `ITERP=1`),
and parameter `INSTANCES` (the Fortran aborts) all raise `NotImplementedError`.

List-control follow-up (executed): re-audit of `SGWF2HFB7RL`/`SGWF2HFB7RLU`
showed barrier lists may begin with `SFAC`/`OPEN/CLOSE`/`EXTERNAL` (both the
per-parameter `NLST` rows and the non-parametric `NHFBNP` rows). `MfUsgHfb` now
consumes these via the shared `_usgt_list.begin_list_block` (reused unchanged):
`SFAC` scales `HYDCHR` per block, `OPEN/CLOSE` resolves against `model_ws` (quoted
names supported), `EXTERNAL` resolves via `ext_unit_dict` (unresolvable →
`NotImplementedError`). The writer still emits expanded inline rows.

Writer-validation polish (executed): `write_file` validates the preserved
parameter state before opening the file, so a `NPHFB>0` header is never written
without a complete, consistent body. From-scratch authoring (`parameters` `None`
or empty `{}`) raises `NotImplementedError`; inconsistent definitions raise
`ValueError` — `len(parameters) != NPHFB`, a def missing
`partyp`/`parval`/`nlst`/`data`, `len(data) != nlst`, `nacthfb !=
len(acthfb_names)`, or an active name not defined (case-insensitive). Files read
by `load` satisfy these invariants, so valid round-trips are unaffected.

Stage 4.4 final polish (executed): `write_file` now also rejects activating the
same parameter more than once (case-insensitive) with `ValueError` — the Fortran
`SGWF2HFB7SUB` aborts "already activated". This matches the duplicate-active
guard added to DRT/QRT, so all four list-parameter packages are consistent.
1 new test (`-k mfusghfb` 20 passed).

Required tests (all green, `-k mfusghfb` 20 passed):

- Static non-parametric authoring. (Done.)
- Transient non-parametric authoring; `IHFBRD=-1`, `0`, `>0`. (Done.)
- Structured and unstructured. (Done.)
- Parameterized load → write → reload preserves defs + activations. (Done.)
- Parameters mixed with non-parametric barriers. (Done.)
- 0-based internal / 1-based file for parameter barrier rows. (Done.)
- From-scratch parametric authoring → `NotImplementedError`. (Done.)
- `TRANSIENT_HFB`+`NPHFB>0` → `NotImplementedError`. (Done.)

Acceptance:

- Parameterized HFB is preserved (load → write → reload), with the unsupported
  variants failing explicitly. No silent partial writes.

### BAS6 - niche options

FloPy class: `MfUsgBas`

Fortran: `glo2basu1.f`

Status: **DONE** (2026-05-30). `RICHARDS_HP` and `IHM [IUIHM]` implemented and
tested; option-line cleaner keeps `_`. Common options already covered.

Problem:

- `RICHARDS_HP` and `IHM [IUIHM]` are not yet explicitly modeled.

Tasks:

- Audit exact Fortran behavior for `RICHARDS_HP`.
- Audit exact Fortran behavior for `IHM [IUIHM]`.
- Decide whether these options need semantic fields or raw option preservation.
- Keep `FREE`, `UNSTRUCTURED`, `DPIN/DPOUT/DPIO`, `SY-ALL`, and
  `SHOWPROGRESS` behavior stable.

Required tests:

- Programmatic BAS with `RICHARDS_HP`, if implemented.
- Programmatic BAS with `IHM [IUIHM]`, if implemented.
- Round-trip of unknown-but-preserved options if raw preservation is chosen.

Acceptance:

- BAS options can be authored from scratch without relying on previous load.

### DPT - immobile-domain adsorption validation

FloPy class: `MfUsgDpt`

Fortran: `gwt2dptu1.f`, `dpt2aw_adsorb.f`

Status: **Reviewed** (2026-05-30). DLIM requires `IDPF/=0 AND IDISPIM/=0`; the
immobile air-water adsorption option `A-W_ADSORBIM` now fails explicitly on
load instead of silently shifting reads.

Problem:

- `DLIM` condition was corrected, but air-water interface adsorption for the
  immobile domain needs audit. (Audited: `A-W_ADSORBIM` reads extra function
  indices/arrays via `dpt2aw_adsorb.f`; unsupported -> explicit load failure.)

Tasks:

- Audit Fortran conditions for DPT adsorption arrays.
- Verify interaction with BCT `A-W_ADSORB`.
- Verify interaction with DPF and `IDISPIM`.
- Add explicit failures for unsupported adsorption sub-modes.

Required tests:

- DPT authoring without DPF.
- DPT authoring with DPF and `IDISPIM=0`.
- DPT authoring with DPF and `IDISPIM!=0`.
- Adsorption-path round-trip or explicit unsupported failure.

Acceptance:

- DPT conditional array layout cannot shift after authoring.

### TIB - semantic parser decision

FloPy class: `MfUsgTib`

Fortran: `glo2basu1.f`

Status: ✅ **Full (authoring)** — **done in Stage 4.1** (reverses the Stage 3
Card 4 "keep raw/text for v1" decision). `MfUsgTib` now has a semantic
`stress_period_data` constructor and a `parse=True` loader that cover the full
`GWF2TIB1RP` grammar below, including transport concentration blocks.
`MfUsgTib.load` still **defaults to byte-exact raw round-trip** (multi-node
`U1DINT` continuation lines preserved; this is what `MfUsg.load` uses) and the
`parse=True` path falls back to raw on unsupported `U1DINT`
(`EXTERNAL`/`OPEN-CLOSE`) rather than writing partial data.

Fortran-derived grammar (`GWF2TIB1RP` in `glo2basu1.f`), per stress period:

```
NIB0 NIB1 NIBM1 [NICB0 NICB1 NICBM1]      # 3 ints; 6 when BCT/transport active
[U1DINT list of NIB0 node ids]            # cells set IBOUND=0 (inactivated)
[NIB1  lines: ICELL [HEAD v | AVHEAD]]    # cells set IBOUND=1
[NIBM1 lines: ICELL [HEAD v | AVHEAD]]    # cells set IBOUND=-1 (CHD-like)
[transport: U1DINT NICB0 list]            # cells set ICBUND=0 (CONC=CINACT)
[transport: NICB1 lines: ICELL [CONC c.. | AVCONC]]   # ICBUND=1
[transport: NICBM1 lines: ICELL [CONC c.. | AVCONC]]  # ICBUND=-1
```

A count `<= 0` skips that block (IBOUND/ICBUND persist; no re-zeroing). The
header is 6 ints whenever the model has an active BCT (mirrors `INBCT`),
3 otherwise. Internal nodes 0-based; file ids 1-based; `CONC` has MCOMP values.

Semantic API (per-SP dict keyed by 0-based `kper`): `ib0`/`icb0` = node arrays;
`ib1`/`ibm1` = `(node, head)` with `head` a float / `"AVHEAD"` / `None`;
`icb1`/`icbm1` = `(node, conc)` with `conc` a length-MCOMP sequence /
`"AVCONC"` / `None`. Authoring validates: transport data requires BCT, conc
length must equal MCOMP, nodes must be 0-based.

Resolution (Stage 4.1):

- Decided **semantic** (was deferred at Card 4); Fortran grammar audited from
  `GWF2TIB1RP` and recorded above before implementing.
- Raw-body preservation kept as the `load` default + the `parse=True` fallback
  for unknown/complex files.

Tests (in `autotest/test_usg_transport.py`):

- Existing raw-body round-trip (`test_mfusgtib_roundtrip`) — unchanged.
- From-scratch non-transport authoring (`..._authoring_nontransport_from_scratch`).
- From-scratch transport authoring with conc blocks (`..._authoring_transport_from_scratch`).
- Semantic non-transport load (`..._semantic_load_nontransport`).
- Semantic write + reload, byte-stable rewrite (`..._semantic_write_reload_roundtrip`).
- Multi-node-per-line `U1DINT` list (covered by the authoring/round-trip tests).
- Raw fallback on unsupported `U1DINT` (`..._parse_falls_back_on_unsupported`).
- Explicit authoring rejections (`..._authoring_rejects_invalid`).

Acceptance:

- From-scratch authoring works without first loading a `.tib`; status is honest
  (raw round-trip remains the safe default and documented fallback).

Stage 4.1 review follow-up (resolved): (1) the three input modes
(`stress_period_data`/`blocks`/`raw_body`) are now mutually exclusive — supplying
more than one raises `ValueError` instead of `write_file` silently preferring
`raw_body`; (2) `parse=True` now raises on a premature EOF (file with fewer
headers than `nper`) and `load` falls back to the raw body rather than expanding
the file with synthetic no-op stress periods. Two tests added
(`test_mfusgtib_rejects_mixed_input_modes`,
`test_mfusgtib_parse_rejects_truncated_file`). Status unchanged: Full (authoring).

Polish: mode-exclusivity now counts a mode as supplied by **explicit presence**
(`is not None`) — the same test `write_file` uses to choose a branch — so an
explicitly-passed but empty mode (`raw_body=""`, `stress_period_data={}`,
`blocks={}`) also counts and can no longer slip past validation into a silent
`write_file` branch. A single explicitly-empty mode stays valid; mixing still
raises. Regression cases added to `test_mfusgtib_rejects_mixed_input_modes`.

---

## Priority 3 - Review Existing "Full" Or Near-Full Packages

These packages look strong, but should still be reviewed before upstream work.
The task is not to add features first; it is to verify that the package really
supports authoring from scratch.

### CHD / RIV / GHB / DRN

Classes: `MfUsgChd`, `MfUsgRiv`, `MfUsgGhb`, `MfUsgDrn`

Tasks:

- Confirm default dtypes match Fortran precision and field order.
- Confirm AUX registration and output order.
- Confirm stress-period reuse semantics.
- Confirm no legacy 1-based internal values remain.

Required tests:

- Authoring with 0-based nodes.
- File inspection for 1-based nodes.
- Reload to 0-based nodes.
- AUX concentration fields.
- Reuse `-1`.

Acceptance:

- Keep status `Full semantic` only if all tests remain green.

### WEL

Class: `MfUsgWel`

Fortran: `gwf2wel7u1.f`

Tasks:

- Keep conceptual boundary clear: WEL applies rates; CLN defines connectivity.
- Verify GWF rates and CLN-node rates have separate stress-period counts.
- Verify `ITMP NP ITMPCLN` behavior.
- Verify AUX concentration fields.

Required tests:

- Programmatic GWF-only WEL.
- Programmatic CLN-rate WEL.
- Mixed GWF and CLN rates.
- Reuse behavior.
- Reload nodes are 0-based.

Acceptance:

- Users can author pumping/injection rates on CLN nodes without defining CLN
  connectivity in WEL.

### RCH / EVT / ETS

Classes: `MfUsgRch`, `MfUsgEvt`, `MfUsgEts`

Tasks:

- Confirm transport concentration arrays and zones are authorable.
- Confirm fixed spacing and item order against Fortran.
- Confirm parameter handling is explicit.

Required tests:

- Programmatic recharge/evaporation with concentration arrays.
- Round-trip files with transport flags.
- Parameter behavior tests for ETS.

Acceptance:

- Boundary transport fields can be authored without prior load.

### CLN

Class: `MfUsgCln`

Fortran: `cln2basu1.f`, `cln2props1.f`

Status: patch candidate or implemented depending on branch review.

Tasks:

- Critically review the API for CLN segment connectivity versus IA/JA input.
- Verify `PROCESSCCF` output-unit handling.
- Verify `ISHAPE` dtype does not break old circular-only format.
- Verify `GENERAL_SEC` data representation is stable and documented.
- Verify structured and unstructured CLN-GWF connection records.

Required tests:

- Programmatic old circular CLN.
- Programmatic rectangular CLN.
- Programmatic general-section CLN.
- Programmatic IA/JA CLN network.
- Round-trip with `PROCESSCCF`.
- Round-trip with structured CLN-GWF connections.
- Round-trip with unstructured CLN-GWF connections.

Acceptance:

- CLN can define connectivity/geometries from scratch. WEL remains rate-only.

### DPF

Class: `MfUsgDpf`

Fortran: `gwf2dpf1u1.f`

Status: patch candidate or implemented depending on branch review.

Tasks:

- Verify `model.idpf` behavior with package add/remove/load.
- Verify layer conditionality for `SC2IM`.
- Verify `IUZONTABIM` in TABRICH mode.
- Verify immobile Richards arrays when `LAYCON=5`.
- Verify BCF and LPF interactions separately.

Required tests:

- Programmatic DPF all-confined.
- Programmatic DPF mixed confined/convertible.
- Programmatic DPF `LAYCON=5` with explicit Richards arrays.
- Programmatic DPF TABRICH with `IUZONTABIM`.
- Round-trip each of the above minimal cases.

Acceptance:

- DPF files authored from scratch do not contain shifted conditional arrays.

### BCT / DDF / TVM

Classes: `MfUsgBct`, `MfUsgDdf`, `MfUsgTvm`

Tasks:

- Maintain existing semantic coverage.
- Add authoring tests for any field that was only validated through real-model
  round-trip.
- Confirm multi-species behavior for BCT.
- Confirm TVM nper+1 boundary behavior.

Required tests:

- Minimal BCT authoring with one species.
- BCT multi-species authoring.
- BCT IDISP=1 and IDISP=2 authoring.
- DDF nonlinear table authoring.
- TVM with and without BCT.

Acceptance:

- Keep `Full semantic` only where authoring tests exist.

---

## Priority 4 - Base FloPy Packages Used In USG-T Context

These packages currently rely on base FloPy classes or text wrappers. The main
question is whether USG-T users need semantic authoring or only compatibility.

### LAK

Class: `MfUsgLak`

Fortran: `gwf2lak7u1.f`

Status: **Stage 4 LAK Fullness Card D executed; decision stays `✅` not Full
(hardened).** From-scratch authoring is now real and tested for the main
branches; six authoring bugs were fixed; bounded gaps keep it `✅`, not `Full`.

Audit (`gwf2lak7u1.f`, ≈4630 lines — the largest MODFLOW package): the option
line is parsed at `GWF2LAK7U1AR` (L83 `TABLEINPUT`, L89 `TRANSPORTBOUNDARY`);
the per-SP geometry/connectivity/gage data is read across `GWF2LAK7U1RP/RPS/RPU`.
Dataset 9b is the transport concentration block: with `TRANSPORTBOUNDARY`
(`ILKTRNSPT=0`) it is **one `CLAKE(1:NSOL)` line per lake** (`gwf2lak7u1.f:1060`);
classic transport (`ILKTRNSPT=1`) is per (lake, component) `CPPT, CRNF [, CAUG]`.

Card D (Stage 4) outcome: fixed six authoring bugs — `conc_data` mis-assignment
(`{0: sill_data}` → `{0: conc_data}`); `write_file` crash when `flux_data` is
None (now required); `conc_data` accessed when `mcomp==0`; `TRANSPORTBOUNDARY`
9b written per-component instead of one line per lake; load stored conc as
strings (now float); `transportboundary` flag not synced to the header keyword.
Added validation (flux_data required; TRANSPORTBOUNDARY needs transport; clake
nlakes×mcomp; transport needs conc_data). Five synthetic from-scratch/round-trip
tests (no-transport, classic transport, TRANSPORTBOUNDARY MCOMP>1, TABLEINPUT,
negatives) plus the `Ex8_Lake` real-model round-trip/run. See
`USGT_STAGE4_LAK_FULLNESS.md`.

Review follow-up (executed): hardened three from-scratch authoring inputs that
`__init__` accepted but `write_file()` then crashed on with a raw
`IndexError`/`KeyError`/`TypeError` — TABLEINPUT now requires exactly one
`tab_file` (and `tab_unit`, if supplied) per lake; `conc_data` is validated per
written period (every `(lake, component)` present; classic = 2 values when
`WTHDRW>=0` / 3 when `WTHDRW<0`; TRANSPORTBOUNDARY = a single value); `flux_data`
requires one dataset-9a entry per lake — all now clear `ValueError`s. Five tests
added (`-k mfusglak` 10 passed; `-k "mfusglak or Ex8"` 11; focused 165; exe 4;
combined 169 ARM). Decision unchanged.

Why still not `Full`: sill/connectivity (datasets 7/8) and multi-lake sublake
systems round-trip (and run via Ex8) but aren't authored from scratch in the
tests; `TABLEINPUT` bathymetry table *contents* are external (FloPy only
registers the per-lake tab unit / external file); GAGE coupling is a separate
package.

Tasks:

- Verify `TABLEINPUT` and `TRANSPORTBOUNDARY` authoring. (Done — Card D.)
- Verify lake transport concentration behavior. (Done — classic + boundary 9b
  authoring tests + Ex8 round-trip.)
- Decide whether current support is full enough for project use. (Yes; kept `✅`
  not Full with the gaps above listed precisely.)

### SFR / STR / GAGE / FHB / SUB / SWT

Classes: base FloPy classes

Fortran: `gwf2sfr7u1.f`, `gwf2str7u1.f`, `gwf2gag7u1.f`,
`gwf2fhb7u1.f`, `gwf2sub7u1.f`

Status: **Decision (2026-05-30)** — **out of scope for USG-T authoring;
compatibility-only.** These use the base MODFLOW-2005 FloPy classes. Their
USG-T unstructured/node-based record variants are not independently validated,
so they are marked `⚠️ Partial` in the roadmap (not silently "supported").
CLN is the preferred surface-water / conduit coupling in this project, so SFR
and STR are not pursued. A USG-T model that relies on unstructured records in
these packages should be validated by the user before use.

Tasks:

- Decide whether these are in scope for USG-T authoring. (Out of scope.)
- If out of scope, document explicitly as compatibility-only. (Done — roadmap
  rows say "base class used … not validated".)

Required tests if in scope:

- Minimal programmatic USG-T variant.
- Round-trip with unstructured nodes.
- NAM load registry behavior.

Acceptance:

- No package should appear silently supported if USG-T-specific records are not
  understood.

### GSF

Class: `MfUsgGsf`

Status: ✅ **Full (authoring)** — **done in Stage 4.2** (supersedes the earlier
"text round-trip is sufficient" decision). `MfUsgGsf` gained a semantic
`vertices` + `node_data` constructor, a `parse=True` loader, a semantic writer,
and `from_grid(model, grid, zverts)`, while keeping the raw `lines` round-trip
as the `load` default and parse fallback.

GSF is a grid-specification file produced by gridgen-style tools and consumed by
post-processors (FloPy's `UnstructuredGrid.from_gridspec`); it is **not** read by
the MODFLOW-USG / USG-T solver (confirmed: no `.gsf`/GRIDSPEC reader in the USG-T
source), so there is no Fortran reader to audit and no executable smoke test —
correctness is validated through `from_gridspec` and semantic write/reload.

Semantic model (0-based internal, 1-based file): header (`UNSTRUCTURED [GWF]`),
`nnodes`, `nlay`, `vertices` as `(nverts, 3)` `(x, y, z)`, and per-node records
`{node, xc, yc, zc, layer, vertices}`. `from_gridspec` reads only `NNODES` from
line 2 and derives the layer count from the per-node `LAY` column, so the
trailing line-2 gridgen integers are preserved verbatim on load and default to
`(1, 1)` when authoring. `from_grid` requires per-vertex `zverts` because
`UnstructuredGrid` does not retain vertex z (it is collapsed into per-cell
top/botm) — it fails explicitly rather than inventing elevations.

Tests (in `autotest/test_usg_transport.py`):

- existing raw round-trip + `to_grid()` smoke (`test_mfusggsf_load_stores_lines`,
  `test_mfusggsf_text_roundtrip`, `test_mfusggsf_to_grid`) — unchanged;
- semantic load (`test_mfusggsf_semantic_load`);
- from-scratch authoring + `to_grid` (`test_mfusggsf_authoring_from_scratch`);
- write/reload byte-stable, 0-based preserved (`test_mfusggsf_write_reload_roundtrip`);
- invalid vertex refs + ambiguous/mixed modes
  (`test_mfusggsf_rejects_invalid_and_mixed_modes`);
- `from_grid` round-trip incl. the missing-z failure (`test_mfusggsf_from_grid`).

Acceptance:

- From-scratch authoring works without first loading a `.gsf`; `load` defaults
  to the safe raw round-trip; status is honest (Full for authoring).

Stage 4.2 review follow-up (resolved): (1) fixed an operator-precedence bug in
`UnstructuredGrid.from_gridspec` so `UNSTRUCTURED GWF` headers parse (accepts
exactly `UNSTRUCTURED` / `UNSTRUCTURED GWF`); (2) `parse=True` rejects trailing
non-comment content (falls back to raw round-trip); (3) `from_grid` now supports
the USG-T top/bottom doubled-vertex convention (`top_zverts`+`bot_zverts`, the
teaching-notebook pattern that `from_gridspec(..., split_vertices=True)`
reconstructs), keeps single-surface `zverts` as explicit legacy, and auto-drops
per-cell closing-duplicate vertices; (4) added
`from_disv_gridprops(model, disv_gridprops, top, botm, skip_degenerate=...)`
(closing-vertex removal, degenerate-cell skip with node renumbering). Hardening:
header / unique-node-id / `nlay >= max(layer)+1` validation. Five tests added
(`test_mfusggsf_unstructured_gwf_header`,
`test_mfusggsf_parse_rejects_trailing_content`,
`test_mfusggsf_from_grid_top_bottom`, `test_mfusggsf_from_disv_gridprops`,
`test_mfusggsf_hardening_rejects`). Status unchanged: Full (authoring).

Vertex-mode follow-up (secondary reference `gridgen2gsf.f90`; primary spec stays
gwutil_a 2.17): added `vertex_mode` to `from_grid`/`from_disv_gridprops`
emulating the two GRIDGEN2GSF layouts — `"shared"`/`"parsimonious"` (default,
neighbours reuse vertex ids) and `"cell"`/`"nonparsimonious"` (unique top/bottom
vertices per cell, 8 per quad, never shared). Both keep top vertices then bottom
vertices so `from_gridspec(split_vertices=True)` recovers top/botm. Documented
the line-2 flags as `IZ IC` (spec 2.17; both must be 1). Test
`test_mfusggsf_vertex_modes` added. Status unchanged: Full (authoring).

Final follow-up (resolved): tightened the semantic contract to spec 2.17.
(1) `_parse_semantic` reuses `_normalize_header` — only `UNSTRUCTURED` /
`UNSTRUCTURED GWF` parse (`UNSTRUCTURED EXTRA GWF` → raw). (2) `IZ`/`IC`
validated: authoring defaults `(1,1)`, accepts omitted as assumed `1 1`, rejects
other values or bad lengths; on load `nnode nlay` and `nnode nlay 1 1` are
semantic while `0 1` / `1 0` / odd token counts → raw. (3) node ids must be
`0..nnodes-1` contiguous/ordered (authoring raises on gap/dup/reorder; a file
with out-of-order `inode` → raw). (4) `vertex_mode="cell"` documented as a
non-shared generalization, not byte-equivalent to the GRIDGEN2GSF quadtree
winding (only the top/bottom-half split is guaranteed). Three tests added
(`test_mfusggsf_parse_header_strict`, `test_mfusggsf_iz_ic_flags`,
`test_mfusggsf_inode_validation`). Status unchanged: Full (authoring).

Stage 4.3 — `gridgen2gsf` utility (separate, does not modify `MfUsgGsf`): a
non-interactive `gridgen_to_gsf(model, source, top, botm, vertex_mode, ...)` in
`flopy/mfusg/gridgen2gsf.py` (exported from `flopy.mfusg`) — the Python
equivalent of the `GRIDGEN2GSF` program. It builds an `MfUsgGsf` from a
`disv_gridprops` dict, a flopy `Gridgen` object (`get_gridprops_disv()`), or an
`UnstructuredGrid`, delegating all authoring to `MfUsgGsf` (no duplicated
geometry logic). Modes `shared`/`parsimonious` and `cell`/`nonparsimonious`;
`top`/`botm` are the GSF surfaces (default unit slab). Tests:
`test_gridgen_to_gsf_disv_modes`, `test_gridgen_to_gsf_source_types`,
`test_gridgen_to_gsf_validation`. See `USGT_STAGE4_03_GRIDGEN2GSF.md`.

Follow-up (resolved): shared/parsimonious now compacts unused vertices for
DISV/Gridgen sources (drops vertices not used by any surviving cell, remaps
`cell2d` ids and per-vertex `top`/`botm`; `skip_degenerate` drops exclusive
vertices too) — the GRIDGEN2GSF vertex-parsimonious behaviour; cell mode needs
no compaction. Docstring reworded to "inspired by" (not a full GRIDGEN2GSF port:
no interactive/definition/quadtree parsing, no refinement/threshold/rotation/
offset/quadtree logic). Focal test filter `-k gridgen_to_gsf`. Two tests added
(`test_gridgen_to_gsf_parsimonious_compacts`,
`test_gridgen_to_gsf_skip_degenerate_compacts`); `MfUsgGsf` unchanged.

### OC — Fullness Card A (Stage 4)

Audited `MfUsgOc` against USG-T 2.7 OC (`glo2basu1.f` SGWF2BAS7I/J/N). Implemented
the `BOOTSTRAPPING` header (authoring + round-trip; written on the **first** OC
line, the only place USG-T parses it — `SGWF2BAS7J` rejects a standalone
`BOOTSTRAPPING` record), and fixed `load` to preserve **layer-qualified**
`PRINT`/`SAVE` actions and **`DDREFERENCE`** on round-trip. Added explicit tests
for those plus per-SP `BOOTSTRAP`/`NOBOOTSTRAP`/`BOOTSTRAPSCALE`/
`NOBOOTSTRAPSCALE`, `SAVE`/`PRINT CONC`/`BUDGET`, and `SAVE IBOUND`.

Decision: **kept `✅ (intentionally not Full)`** with explicit gaps — `SAVE
IBOUND` is commented out in USG-T 2.7 (solver rejects; FloPy still preserves the
keyword); `FASTFORWARD`/`FASTFORWARDC` separate-line placement and
`BOOTSTRAPPING` execution are not exe-verified; numeric-format OC rewrites as
words. See `USGT_STAGE4_OC_FULLNESS.md`. (EVT / MDT / LAK untouched.)

Polish follow-up: `write_file` now emits a DDREFERENCE-only period line (it was
dropped before, leaking the flag to the next period); `check()` recognises the
USG-T OC actions (BOOTSTRAP toggles + DDREFERENCE as single words; ATS/solver
params `DELTAT`/`TMINAT`/`TMAXAT`/`TADJAT`/`TCUTAT`/`HCLOSE`/`BTOL`/`MXITER`
require a value, else flagged) instead of false "ignored" warnings, and emits a
specific warning that `SAVE IBOUND` is preserved by FloPy but rejected by
USG-T 2.7. Four tests added; `-k mfusgoc` 11 passed.

### EVT — Fullness Card B (Stage 4)

Audited `MfUsgEvt` against `gwf2evt8u1.f`. Implemented/fixed: from-scratch
authoring + round-trip for all `NEVTOP` modes (structured + unstructured
`NEVTOP=2` with `MXNDEVT`; `IEVT` 0-based internal / 1-based file with layer-range
validation), transport `IETFACTOR` 0/<0/>0 with per-`MCOMP` `ETFACTOR`, and per-SP
reuse (`-1`). Fixed a real round-trip bug (`load` dropped `IETFACTOR`) and made
the writer emit the 3-integer dataset-1 line whenever transport is active. Added
a USG-T 2.7 EVT executable smoke (`test_usgt_exe_evt_from_scratch`). Six
synthetic tests + one exe test.

Decision: **kept `✅ (intentionally not Full)`** — explicit gaps: ETS zonal
time-series (`ETS MXZNEVT`/`IZNEVT`) raises `NotImplementedError`; `NPEVT`
parameters are *Expanded valid write* (loaded as arrays, `NP=0` on write,
authoring-with-params unsupported) — same treatment as ETS. See
`USGT_STAGE4_EVT_FULLNESS.md`. (OC / MDT / LAK untouched beyond docs.)

Review follow-up: fixed a scalar-`ETFACTOR` crash (`__init__` normalizes
`self.etfactor` to a 1-D array, so `etfactor=2.5` with MCOMP=1 writes correctly)
and added unstructured `NEVTOP=2` `IEVT` node-range validation (0-based in
`[0, NODES-1]`, file 1-based; `< 0` or `>= NODES` raises). Two tests added;
`-k mfusgevt` 9 passed. Status unchanged.

### MDT — Fullness Card C (Stage 4)

Audited `MfUsgMdt` against `gwt2mdtu1.for`. Fixed three real bugs: `load`
dropped `FRAHK`/`FRADARCY` (it upper-cased the line but searched lowercase
keywords — e.g. the Ex7 Multispecies `FRAHK` was lost on round-trip);
`write_file(f=handle)` crashed (`f_obj` only assigned when `f is None`, and it
closed external handles); and a stray debug `print`. Added from-scratch
authoring + round-trip for all main branches (header options only when
`IDPF==0`, base arrays with `VOLFRACMD` skipped when `IDPF!=0`, per-species
`KDMD`/`DECAYMD`/`YIELDMD`/`DIFFMD`, `AIOLD1MD`/`AIOLD2MD` under `TSHIFTMD>0`,
multi-species) plus validation (FRAHK⊕FRADARCY, IDPF-options, MULTIFILE needs
`imdtcf>0`, per-`MCOMP` list lengths). Eight synthetic tests; the three Ex7
real-model round-trip/run tests kept.

Decision: **kept `✅ (intentionally not Full)`** — gaps: the per-species loop
uses `MCOMP` (chained-decay `NTCOMP>MCOMP` not independently verified) and the
AI1/AI2 output binaries (`MULTIFILE_MD`/`SEPARATE_AI2`) are authored but not
read. See `USGT_STAGE4_MDT_FULLNESS.md`. (OC / EVT / LAK untouched beyond docs.)

Review follow-up: fixed a TSHIFTMD misalignment — USG-T reads AIOLD only when
`TSHIFTMD > 1e-10`, but the writer used `> 0` for AIOLD and a `{:9.2f}` format
that rounded a small valid value (e.g. `1e-6`) to `0.00`. A module constant
`MDT_TSHIFT_THRESHOLD = 1e-10` is now used everywhere (constructor, writer,
IDPF-options check, load) and TSHIFTMD is written with a general format. One test
added; `-k mfusgmdt` 9 passed. Status unchanged.

---

## Priority 5 - Post-Processing And Validation

### MfusgTransportListBudget

Class: `MfusgTransportListBudget`

Status: **Done** — covered by `test_mfusg_transport_list_budget_old_format`,
`..._new_format`, and `..._species_isolation` (species-specific budget values
differ and stay isolated; old `VOLUMETRIC BUDGET` and new `MASS BUDGET`
keywords both parsed).

Tasks:

- Confirm old-format and new-format transport budget parsing. (Tested.)
- Confirm multiple species are isolated. (Tested.)
- Confirm flow and transport budget blocks are not mixed. (Tested.)

Required tests:

- Old USG-T format.
- New USG-T format.
- Multi-species listing.
- Species-specific budget values differ and remain isolated.

Acceptance:

- Users can inspect transport budgets per species without manual listing-file
  parsing.

### Slow real-model validation

Status: **Decision (2026-05-30; executable tier added Stage 3 Card 9).** Three
tiers are in place:

1. **Default CI (fast):** `test_usg_transport.py` has 101 tests — 85
   synthetic from-scratch authoring / round-trip tests that always run, plus 16
   `Ex1..Ex9` real-model run tests gated by `@requires_exe(USGT_EXE)`. With a
   USG-T executable that is **101 passed**; without one it is **85 passed, 16
   skipped**. Catches layout/API regressions in well under a minute.
2. **Executable end-to-end (opt-in via `USGT_EXE`):** the `Ex1..Ex9` tests and
   the new `autotest/test_usg_transport_exe.py` from-scratch tests share one
   contract — `@requires_exe(USGT_EXE)`, `USGT_EXE` defaulting to `mfusg_gsi`
   (a name on `PATH` or an absolute path; must be **USG-Transport 2.7**). When
   the executable resolves they **run** the models and check outputs —
   from-scratch flow (analytical heads + closed `.list` budget) and transport
   (`.con` produced + closed species-isolated mass budget), plus the nine real
   models (incl. Ex7 multi-species). They skip cleanly when the exe is absent.
3. **Manual run-validation (slow, out of CI):** executing USG-T on the full
   real-world Vistas models and comparing `LST`/`HDS`/`CON`/`CBB` bit-for-bit
   (Model A: BCT+CLN+TIB, ~112k nodes; Model B: BCT IDISP=2 + DDF, ~19k nodes,
   4322 time steps — all `max|diff| = 0`). Those models are large/proprietary
   and not bundled; results are recorded in `USGT_improvements.md`.

Tasks:

- Keep real-world models out of default CI. (Done — only load+write Ex* tests.)
- Create a slow/optional executable tier. (Done — `USGT_EXE`-gated suite.)
- Validate run outputs against reference outputs. (Done — see improvements.md.)

---

## Suggested Agent Work Packages

### Agent A - BCF/LPF TABRICH

Scope:

- `MfUsgBcf`, `MfUsgLpf`, TABRICH arrays only.

Deliverables:

- Fortran-derived TABRICH spec.
- Constructor/load/write support for `IUZONTAB` and `RETCRVS`.
- BCF and LPF authoring tests.
- BCF and LPF round-trip tests.

Do not:

- Refactor unrelated LPF/BCF flow properties.

### Agent B - DRT

Scope:

- USG-T DRT8 extensions only.

Deliverables:

- Fortran-derived DRT spec.
- `MfUsgDrt` design recommendation: subclass or new class.
- Semantic implementation or explicit staged plan.
- Tests for AUX concentration, `IQCHANGEC`, and `MXSPREADNDS`.

Do not:

- Change WEL/CLN behavior.

### Agent C - SGB

Scope:

- New SGB package.

Deliverables:

- Fortran-derived SGB spec.
- `MfUsgSgb` class.
- NAM registry mapping.
- Minimal authoring and round-trip tests.

Do not:

- Approximate SGB with GHB/CHD in code.

### Agent D - QRT

Scope:

- New QRT package.

Deliverables:

- Fortran-derived QRT spec.
- `MfUsgQrt` class.
- NAM registry mapping.
- Minimal authoring and round-trip tests.
- Shared helper proposal with DRT if useful.

Do not:

- Hide unsupported return-flow modes behind partial writes.

### Agent E - Critical Review Of Patch Candidates

Scope:

- Review current patch candidates for NAM, CHD/RIV/GHB/DRN/WEL, HFB, ETS, TIB,
  BAS, CLN, and DPF.

Deliverables:

- Findings ordered by severity.
- Fortran citations.
- Test gaps.
- Recommendation: accept, revise, or revert/rework.

Do not:

- Implement fixes unless explicitly assigned after review.
