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
  (authoring)** in Stage 4.2; LAK validated via Ex8 (authoring deferred);
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
`CHANGEC`/`IDCHNGTYP`, AUX, reuse; `NPDRT>0` fails explicitly; structured
delegates to base. Authoring + round-trip + reuse tests added.

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

- **DONE** (2026-05-30): `MfUsgSgb` added (`flopy/mfusg/mfusgsgb.py`) and
  registered as `"sgb"`. Node-based `(node, gradient)`, AUX, `ITMP/-1` reuse,
  0-based internal / 1-based file, `NPSGB>0` explicit failure. Authoring +
  round-trip + NAM-registry + parameter-failure tests added.

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

- **DONE** (2026-05-30): `MfUsgQrt` added (`flopy/mfusg/mfusgqrt.py`) and
  registered as `"qrt"`. Per-sink `(node, q, rfprop)` + recipient-node lists,
  `CHANGEC`/`IQCHNGTYP`, AUX, reuse; `NPQRT>0` and `TRANSIENTQ` explicit
  failures; shares `_usgt_returnflow.py` with DRT. Authoring (minimal,
  return-flow concentration, multi-recipient) + round-trip + registry tests.

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
non-parametric ETS. **Authoring tests added** (2026-05-30): NETSEG=1,
NETSEG>1, NETSOP=2, IESFACTOR, and explicit `npets>0` write failure.

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
parameterized HFB. **Tests cover** (2026-05-30): static, structured static,
transient `IHFBRD=>0/0/-1`, and explicit `NPHFB>0` load failure.

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

Status: **Decision confirmed (Stage 3 Card 5): keep `✅` not Full.** Current
support is sufficient for project use; from-scratch authoring deferred.

Audit (`gwf2lak7u1.f`, ≈4630 lines — the largest MODFLOW package): the option
line is parsed at `GWF2LAK7U1AR` (L83 `TABLEINPUT`, L89 `TRANSPORTBOUNDARY`);
the per-SP geometry/connectivity/gage data is read across `GWF2LAK7U1RP/RPS/RPU`.
`MfUsgLak` already parses and preserves the `TABLEINPUT` option and the
`transportboundary` flag, and load/write round-trips the `Ex8_Lake` real model.

Why deferred: a from-scratch authoring test must build complete lake
connectivity, bathymetry/stage tables, and (with transport) per-lake boundary
concentrations — a large surface for a package GUIs normally generate. Per the
Stage 3 priority guidance ("a clean, well-tested core beats full support for
rare packages"), this is **deferred**, not blocked. If a target model needs
from-scratch LAK authoring, open a dedicated card to add a minimal synthetic
lake plus a round-trip test and only then promote toward `Full`.

Tasks:

- Verify `TABLEINPUT` and `TRANSPORTBOUNDARY` authoring. (Deferred — see above.)
- Verify lake transport concentration behavior. (Covered by Ex8 round-trip.)
- Decide whether current support is full enough for project use. (Yes.)

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
