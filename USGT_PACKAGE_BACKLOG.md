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

Status: patch candidate exists.

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

Status: patch candidates for CHD/RIV/GHB/DRN/WEL exist.

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

Status: partial. This is the clearest remaining package gap for programmatic
authoring of Richards/TABRICH models.

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

Status: partial. Base MODFLOW DRT is not enough for USG-T transport authoring.

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

Status: not implemented.

Problem:

- SGB appears in USG-T 2.7 but has no package class or load registry support.
- Models using SGB may be silently incomplete after `MfUsg.load()`.

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

Status: not implemented.

Problem:

- QRT is absent from FloPy USG-T registry.
- It overlaps conceptually with DRT but has its own package I/O and transport
  return-flow behavior.

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

Status: expanded valid write for parameterized input; semantic for
non-parametric ETS.

Problem:

- Parameterized files can be loaded and expanded.
- Original parameter syntax is not preserved.
- Programmatic `NPETS>0` should not write misleading incomplete syntax.

Tasks:

- Decide product direction: preserve parameter syntax or formally keep expanded
  valid write as v1.
- If preserving parameters, design parameter object storage and active-parameter
  records.
- If not preserving parameters, ensure explicit failure paths stay covered.

Required tests:

- `NETSEG=1` authoring.
- `NETSEG>1` authoring.
- `NETSOP=2` authoring.
- `IESFACTOR` authoring.
- Parameterized load writes expanded `NPETS=0`.
- Programmatic `npets>0` fails explicitly.

Acceptance:

- Users can author non-parametric ETS from scratch.
- Parameter behavior is explicit and documented.

### HFB - parameterized barriers

FloPy class: `MfUsgHfb`

Fortran: `gwf2hfb7u1.f`

Status: semantic for non-parametric static/transient HFB; partial for
parameterized HFB.

Problem:

- `NPHFB>0` is not supported.
- Non-parametric Fortran semantics around `IHFBRD` must remain protected.

Tasks:

- Audit parameter read/active behavior in Fortran.
- Decide expanded valid write versus parameter preservation.
- Add parameter storage or explicit unsupported failure.
- Preserve current non-parametric behavior.

Required tests:

- Static non-parametric authoring.
- Transient non-parametric authoring.
- `IHFBRD=-1`, `0`, and `>0`.
- Structured and DISU formats where applicable.
- Parameterized input behavior: expanded or explicit fail.

Acceptance:

- Parameterized HFB behavior is no longer ambiguous.

### BAS6 - niche options

FloPy class: `MfUsgBas`

Fortran: `glo2basu1.f`

Status: mostly covered for common USG-T options; niche options remain.

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

Status: partial.

Problem:

- `DLIM` condition was corrected, but air-water interface adsorption for the
  immobile domain needs audit.

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

Status: raw/text round-trip.

Problem:

- Raw preservation is useful for existing files but not enough for rich
  programmatic authoring.
- TIB contains `U1DINT` style lists that are easy to mis-parse.

Tasks:

- Decide whether TIB should remain raw/text v1 or gain semantic constructors.
- If semantic, write a Fortran-derived grammar/spec first.
- Keep raw-body preservation available for unknown or complex files.

Required tests:

- Existing raw-body round-trip.
- Programmatic minimal TIB if semantic mode is added.
- Multi-node-per-line `U1DINT` lists.

Acceptance:

- Status is honest: raw/text if no semantic authoring API exists.

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

Tasks:

- Verify `TABLEINPUT` and `TRANSPORTBOUNDARY` authoring.
- Verify lake transport concentration behavior.
- Decide whether current support is full enough for project use.

Required tests:

- Programmatic LAK with transport boundary.
- Round-trip LAK with table input.

### SFR / STR / GAGE / FHB / SUB / SWT

Classes: base FloPy classes

Fortran: `gwf2sfr7u1.f`, `gwf2str7u1.f`, `gwf2gag7u1.f`,
`gwf2fhb7u1.f`, `gwf2sub7u1.f`

Tasks:

- Decide whether these are in scope for USG-T authoring.
- If in scope, audit USG-T unstructured/node-based differences.
- If out of scope, document explicitly as compatibility-only or unsupported.

Required tests if in scope:

- Minimal programmatic USG-T variant.
- Round-trip with unstructured nodes.
- NAM load registry behavior.

Acceptance:

- No package should appear silently supported if USG-T-specific records are not
  understood.

### GSF

Class: `MfUsgGsf`

Status: text round-trip plus `to_grid()` helper.

Tasks:

- Decide whether text round-trip is sufficient.
- If semantic editing is needed, define a separate API from raw preservation.

Required tests:

- Raw round-trip.
- `to_grid()` smoke test.

Acceptance:

- Status remains honest: text wrapper unless semantic editing is added.

---

## Priority 5 - Post-Processing And Validation

### MfusgTransportListBudget

Class: `MfusgTransportListBudget`

Tasks:

- Confirm old-format and new-format transport budget parsing.
- Confirm multiple species are isolated.
- Confirm flow and transport budget blocks are not mixed.

Required tests:

- Old USG-T format.
- New USG-T format.
- Multi-species listing.
- Species-specific budget values differ and remain isolated.

Acceptance:

- Users can inspect transport budgets per species without manual listing-file
  parsing.

### Slow real-model validation

Tasks:

- Keep real-world models out of default CI.
- Create a slow marker or manual script.
- Validate run outputs against reference outputs:
  LST, HDS, CON, CBB, and package-specific outputs where relevant.

Acceptance:

- Synthetic tests catch layout/API regressions.
- Real models catch integration regressions.

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

