# Stage 4.6 — Final USG-T 2.7 Gap Audit (no code changes)

Date: 2026-06-02 · Base: `develop` @ `e0ea8073` (after Stage 4.5B + polish + cleanup)

Purpose: before implementing more features, take stock of where FloPy stands on
**programmatic authoring of USG-T 2.7 packages from scratch** (the primary Stage
4 objective), classify every package, enumerate the concrete remaining gaps,
prioritize them by practical impact, and propose the next implementable cards.

This is a docs-only audit. No package code, tests, or binaries were changed.

> **Numbering note.** The master plan's Stage 4.6/4.7/4.8 are already written as
> `USGT_STAGE4_06_DPT_AW_ADSORBIM.md` / `_07_COMPAT_PACKAGES.md` /
> `_08_UPSTREAM_INFRA.md`. To avoid a filename clash this audit is filed as
> `_09_`, but it covers and **re-sequences** the remaining work (those three
> planned stages are folded into the prioritized backlog below). The proposed
> implementation cards are labelled **Stage 4.6A, 4.6B, …** per the review
> request; they supersede the old flat 4.6→4.8 ordering.

Method: cross-referenced `USGT_roadmap.md` (the authoritative status table),
`USGT_PACKAGE_BACKLOG.md`, `USGT_improvements.md`, all `USGT_STAGE4_*.md`, the
live `flopy/mfusg/*.py` classes and the `MfUsg.mfnam_packages` registry, and the
USG-T 2.7 Fortran (`USGT_V_2-7-0_Source_Code/`, incl. the `mfusg.f` `CUNIT`
array and `utl7u1.f`).

---

## 1. Headline findings

1. **Authoring from scratch is broad and solid.** ~20 packages are `Full`
   (Fortran-spec + from-scratch authoring tests + round-trip + explicit
   unsupported-mode failures). The flow core, transport core (BCT/PCB/DDF/DPF),
   CLN, TVM, and the boundary conditions all author from Python.
2. **No silent-wrong-write P0 remains.** The project's consistent
   *fail-explicitly* discipline (60+ `NotImplementedError` / writer `ValueError`
   sites across `flopy/mfusg`) means unsupported modes refuse to write rather
   than emit a partial/misread file. The one *latent* P0 this audit flagged — the
   base-class compatibility packages — was **closed by Stage 4.6A**: STR/SUB/SWT
   are now guarded (`mfusgcompat.py`) so unstructured-`MfUsg` use fails
   explicitly. SFR was *excluded* from the guard because it is empirically
   DISU-validated (`freyberg_usg` loads/writes/runs it via base `ModflowSfr2`);
   FHB/GAGE remain compatibility (Ex8). See §4 / §3.2–3.3.
3. **The dead `CUNIT` slots are truly dead.** `EVS`, `RTS`, `RES`, `IBS` appear
   in `mfusg.f`'s `CUNIT` array but have **no reader** in USG-T 2.7 (`RES`/`IBS`
   dispatch is commented out with `CSP`; `EVS`/`RTS` have no `CALL` sites). So
   FloPy not registering them is correct — there is no load-skip data-loss gap
   there. All *functional* USG-T 2.7 packages are registered in
   `mfnam_packages`.
4. **The largest authoring gap is from-scratch MODFLOW *parameters*.** ETS, HFB,
   DRT, SGB, QRT all *preserve* loaded `NP>0` definitions but raise
   `NotImplementedError` when asked to author parameter definitions from
   scratch. This is the main remaining "author a common feature from zero" gap.

---

## 2. Per-package classification

Buckets: **FA** = Full authoring from scratch · **FS-exec** = Full semantic
load/write but with execution/coverage gaps (intentionally not `Full`) ·
**ParamPreserve** = parameter-preserving only (no from-scratch param authoring)
· **Expanded** = expanded valid write only (controls not preserved) · **Raw** =
raw/text round-trip · **Compat** = compatibility-only (base class).

| Package | Bucket(s) | Notes |
|---|---|---|
| BAS6, DIS, DISU, BCF6, LPF, SMS | FA | Core; verified, real-model validated. |
| CHD, WEL, DRN, RIV, GHB, RCH | FA | Node-based BCs, AUX transport conc. |
| BCT, PCB, DDF, DPF, TVM | FA | Transport core + density + TVM. |
| CLN | FA | Circular/rect/`GENERAL_SEC`/`PROCESSCCF`/`ISHAPE`. |
| TIB | FA (`parse=True`) | `load` default is byte-exact **Raw**; semantic on `parse=True`; raw fallback on `U1DINT` EXTERNAL/OPEN-CLOSE. |
| GSF | FA (`parse=True`) | Not solver input. `load` default **Raw**; semantic + `from_grid`/`to_grid`. |
| OC | FS-exec | Broad authoring; gaps: `SAVE IBOUND` solver-rejected, FASTFORWARD/BOOTSTRAPPING exec unverified, numeric→words rewrite. |
| EVT | FS-exec | All `NEVTOP`, transport `IETFACTOR`, reuse; gaps: ETS-zonal time series `NotImplementedError`, `NPEVT` Expanded-only. |
| MDT | FS-exec | Audited authoring; gaps: chained-decay `NTCOMP>MCOMP` unverified, AI1/AI2 binaries not read. |
| LAK | FS-exec | No-transport/classic/`TRANSPORTBOUNDARY`/`TABLEINPUT` authoring; gaps: sill/connectivity (ds 7/8) + multi-lake not authored from scratch; bathymetry tables external. |
| GNC | FS-exec (✅) | Ghost-node helper; verified. |
| ETS | FA (incl. **from-scratch ETSR array-param authoring**, Stage 4.6C-D) + ParamPreserve (ETSR array) | Only ETSR parameterized (Fortran limit); parametric execution not exe-smoke-tested; opt-in expand fallback. |
| HFB | FA (incl. **from-scratch `NPHFB>0` authoring**, Stage 4.6C-A) + ParamPreserve (list) | Gaps: `TRANSIENT_HFB`+`NPHFB>0`, `INSTANCES`. |
| DRT | FA (incl. **from-scratch `NPDRT>0` authoring**, Stage 4.6B) + ParamPreserve + Expanded (list/recip controls) | Type-consistent; from-scratch parameter authoring supported; activated SPREAD = structural round-trip only (Fortran copies `DRTF` not `NodDRT`); `INSTANCES` unsupported. Unstructured only. |
| SGB | FA + **Param-definition preserve + from-scratch authoring** (Stage 4.6C-B; no activations) + **Expanded** | Definitions author/round-trip; active SGB params abort the Fortran (`PARTYP='SGB'` vs `'G'`) so `NP>0`/`active_params`/`INSTANCES` → `NotImplementedError`. |
| QRT | FA (incl. **from-scratch `NPQRT>0` authoring**, Stage 4.6C-C; structural) + ParamPreserve + Expanded + `TRANSIENTQ` | Active param value scales `NumRT` not `Q` (Fortran bug); `NodQRT` not copied → structural, not exec-guaranteed. `INSTANCES`/`TRANSIENTQ`+`NPQRT>0` → `NotImplementedError`. |
| DPT | FS-exec **except** `A-W_ADSORBIM` (explicit fail) | Immobile air-water adsorption unsupported (rare PFAS sub-mode). |
| SFR | **Compat (DISU-validated)** | Base `ModflowSfr2`, kept enabled: `freyberg_usg` (DISU+SFR) loads → writes → runs via the base class in `test_usg.py`. *Not* a Full USG-T transport claim; not guarded. |
| FHB, GAGE | **Compat** | Base classes; round-trip via base in the `Ex8` real model. Untouched by Stage 4.6A. |
| STR, SUB, SWT | **Compat (guarded, Stage 4.6A)** | Base layout unvalidated for USG-T unstructured; `mfusgcompat.py` wrappers raise `NotImplementedError` on load/authoring for an unstructured `MfUsg` model. Used by no target model. |
| HOB/DROB/RVOB/GBOB/CHOB/hyd, MNW1/2, UZF, GMG/PCG/SIP/DE4, LVDA/KDEP, SYF, lmt6, gwm, PATH/PTH | Intentionally absent | Observations/HUF/coupling/solvers; documented out of scope. |
| EVS, RTS, RES, IBS | n/a (dead) | `CUNIT` placeholders with no functional reader in USG-T 2.7. |

---

## 3. Gap inventory by category

### 3.1 Authoring-from-scratch gaps (primary objective)

- **A1 — From-scratch list/array *parameter* authoring — DONE.** The full
  family is complete: **DRT** (4.6B, pilot), **HFB** (4.6C-A, non-transient),
  **SGB definitions** (4.6C-B — active SGB params stay unsupported by Fortran),
  **QRT** (4.6C-C, structural — active params not exec-guaranteed by Fortran),
  and **ETS** (4.6C-D — the *array*-parameter form for ETSR, via a `ModflowParBc`
  builder, a different path from the list-parameter helpers). The biggest
  authoring gap from the audit is now closed.
- **A2 — LAK multi-lake + sill/connectivity (ds 7/8)** not authored from scratch
  (single-lake + basic transport authored; multi-lake systems only round-trip).
- **A3 — DPT `A-W_ADSORBIM`** immobile air-water adsorption unsupported (cascade
  of conditional arrays; rare PFAS use). Planned as old Stage 4.6.
- **A4 — EVT ETS-zonal time series** (`ETS MXZNEVT`/`IZNEVT`) `NotImplementedError`;
  **`NPEVT` parameters** Expanded-only (no param authoring).
- **A5 — HFB `TRANSIENT_HFB`+`NPHFB>0`** and **`INSTANCES`** (all param packages)
  unsupported.

### 3.2 Parser / load gaps

- **L1 — Silent load-skip of an unregistered NAM filetype** when
  `model.verbose` is False (`mfusg.py` prints "package load...skipped" only under
  verbose). Impact is *theoretical*: every functional USG-T 2.7 package is
  registered; only intentionally-absent obs/HUF/coupling packages would skip.
- **L2 — Compat-package load via base classes — addressed (Stage 4.6A).**
  - **SFR:** kept on base `ModflowSfr2`; the USG-T unstructured layout is
    empirically supported (`freyberg_usg` DISU+SFR loads via base in
    `test_usg.py`). Compatibility/base support; **no** Full USG-T transport claim.
  - **STR/SUB/SWT:** their USG-T unstructured layout is unvalidated, so they are
    now **guarded** — `mfusgcompat.py` raises `NotImplementedError` on load for an
    unstructured `MfUsg` model instead of mis-reading via the structured base.
  - **FHB/GAGE:** compatibility retained (proven via the `Ex8` real model).

### 3.3 Writer gaps

- **W1 — Compat-package write via base classes — resolved (Stage 4.6A).**
  STR/SUB/SWT are guarded so a write/authoring on an unstructured `MfUsg` model
  raises before any file is produced (no silent structured-format output). SFR
  stays enabled (DISU-validated via `freyberg_usg`); FHB/GAGE unchanged (Ex8).
  This closed the audit's one latent P0.
- **W2 — OC `SAVE IBOUND`** preserved on write though USG-T 2.7's solver rejects
  it (`check()` warns); numeric-format OC is rewritten as words.

### 3.4 Execution gaps (USGT_EXE)

- **E1 — QRT `NP>0` activation** not execution-guaranteed (Fortran bug: value
  scales `NumRT` not `Q`; `NodQRT` not copied). FloPy preserves structurally and
  does not apply the value to `Q`. *Not a FloPy defect.*
- **E2 — OC FASTFORWARD/FASTFORWARDC + BOOTSTRAPPING** placement/execution not
  exe-verified.
- **E3 — MDT chained-decay `NTCOMP>MCOMP`** not independently exe-verified.
- (E-closed — SGB active params abort the run; FloPy already refuses to write.)

### 3.5 Docs / tests / cleanup gaps

- **D1 — Stale `# todo` markers** in `mfusgdisu.py` ("todo: check",
  "model.structured = False # todo: why?") and `mfusgoc.py` (layer-list todo) —
  DISU/OC are validated, so these are cleanup-only.
- **D2 — Stage-doc numbering** drift (this audit filed as `_09_`; the master
  plan's 4.6→4.8 are superseded here).
- **D3 — Compat-package status is documented but not enforced** (see W1/L2).

---

## 4. Prioritized backlog

Priorities use the review's definitions:
**P0** = may write an invalid or semantically different file with no warning ·
**P1** = blocks from-scratch authoring of a common feature ·
**P2** = incomplete preservation or missing tests ·
**P3** = docs / lint / cleanup.

### P0 — silent invalid/wrong write

- **W1 — compat packages — RESOLVED (Stage 4.6A, executed).** A revised audit
  found the original premise was too broad: **SFR is empirically supported on
  DISU** — `examples/data/freyberg_usg` (DISU+SFR) loads → writes → runs via base
  `ModflowSfr2` in `autotest/test_usg.py`. So SFR is **kept on the base class**
  (compatibility, DISU-validated; *not* a Full USG-T transport claim). The real
  latent risk — **STR/SUB/SWT** (no DISU validation, used by no model) — is now
  closed: `flopy/mfusg/mfusgcompat.py` registers guarded wrappers
  (`MfUsgStr`/`MfUsgSub`/`MfUsgSwt`) so unstructured-`MfUsg` load **or** authoring
  raises `NotImplementedError` before any read/write. FHB/GAGE untouched
  (Ex8). No remaining P0.

  *No other P0s.* The explicit-failure discipline (unsupported modes raise before
  opening the file; writers validate and refuse partial output) has held through
  Stages 3–4.5.

### P1 — blocks from-scratch authoring of common features

- **A1 — from-scratch parameter authoring.** ✅ **DONE** — DRT (4.6B), HFB
  (4.6C-A), SGB defs (4.6C-B), QRT (4.6C-C, structural), and ETS array (4.6C-D)
  all author `NP*>0` from scratch.
- **A2 — LAK multi-lake + sill/connectivity (ds 7/8) authoring** — now the top
  P1 remaining.

### P2 — incomplete preservation / missing tests

- **A3 — DPT `A-W_ADSORBIM`** (implement or formally keep deferred — old 4.6).
- **A4 — EVT ETS-zonal time series + `NPEVT` param authoring.**
- **A5 — HFB `TRANSIENT_HFB`+`NPHFB>0`; `INSTANCES` (all param packages).**
- ~~**L2 — compat-package load validation.**~~ Addressed by Stage 4.6A
  (STR/SUB/SWT guarded; SFR DISU-validated via `freyberg_usg`; FHB/GAGE via Ex8).
- **W2 — OC `SAVE IBOUND` / numeric-format** behavior tightening.
- **E2/E3 — OC BOOTSTRAPPING & MDT chained-decay exe verification.**

### P3 — docs / cleanup

- **D1 — remove/resolve stale DISU/OC `todo` markers.**
- **L1 — make load emit a one-line warning for an unregistered NAM filetype**
  regardless of `verbose` (defensive; low impact).
- **D2 — reconcile stage-doc numbering** (this file + master plan).

---

## 5. Proposed implementation sequence (Stage 4.6A …)

Ordered P0-latent first, then by authoring impact. One reviewable card each.

- **Stage 4.6A — Compatibility-package USG-T guard — EXECUTED (revised scope).**
  Guards **STR/SUB/SWT only** (`flopy/mfusg/mfusgcompat.py` wrappers registered
  in `MfUsg.mfnam_packages`): unstructured-`MfUsg` load/authoring →
  `NotImplementedError` before any read/write. **SFR excluded** — DISU+SFR is
  base-validated by `freyberg_usg` (`test_usg.py`), so guarding it would break
  passing tests and contradict reality; SFR stays compatibility/base support.
  FHB/GAGE untouched. Structured `MfUsg` and plain `flopy.modflow.Modflow`
  unaffected. 3 negative/regression tests added (`-k compat`); freyberg SFR
  tests still pass; transport focused **241**, combined **245** (ARM).

- **Stage 4.6B — From-scratch list-parameter authoring: DRT (pilot) — EXECUTED.**
  `MfUsgDrt` builds a parameterized package from Python (`parameters=`,
  `active_params=`; ergonomic input normalized + validated before open) and
  writes a valid `NPDRT>0` file that reloads with the same semantics. Only
  `mfusgdrt.py` touched in code. `INSTANCES` still explicit-fail; HFB/SGB/QRT/ETS
  unchanged. 7 new tests + 1 repurposed; `-k mfusgdrt` **34**, focused **248**,
  combined **252** (ARM). See `USGT_STAGE4_10_DRT_PARAMETER_AUTHORING.md`.

- **Stage 4.6C-A — HFB from-scratch param authoring — EXECUTED.** `MfUsgHfb`
  builds a parameterized package from Python (`parameters=`/`acthfb_names`; auto
  counts; structured + unstructured) and writes a valid `NPHFB>0` file that
  reloads with the same semantics. Factored the shared name/`parval` validation
  into `_usgt_parameters` (DRT delegates). `TRANSIENT_HFB`+params / `INSTANCES`
  still fail. 10 new + 1 repurposed test; `-k mfusghfb` **30**, focused **263**,
  combined **267** (ARM). See `USGT_STAGE4_11_HFB_PARAMETER_AUTHORING.md`.
- **Stage 4.6C-B — SGB from-scratch *definition* authoring — EXECUTED.**
  `MfUsgSgb` authors `NPSGB>0` parameter definitions from Python
  (`parameters=`; `PARAMETER NPSGB MXS` + `UPARLSTRP` defs; counts auto;
  AUX via dtype; `NP=0` every period). Active SGB params stay unsupported
  (Fortran `PARTYP='SGB'` vs `'G'`). Only `mfusgsgb.py` + the shared
  `_usgt_parameters` docstring touched. 9 new + 1 repurposed test; `-k mfusgsgb`
  **24**, focused **272**, combined **276** (ARM). See
  `USGT_STAGE4_12_SGB_PARAMETER_AUTHORING.md`.
- **Stage 4.6C-C — QRT from-scratch list-parameter authoring (structural) —
  EXECUTED.** `MfUsgQrt` builds a parameterized package from Python
  (`parameters=`/`active_params=`; auto counts; recipient_nodes) and writes a
  valid `NPQRT>0` file that reloads with the same semantics. **Structural only**
  (Fortran scales `NumRT` not `Q`; `NodQRT` not copied). `INSTANCES` /
  `TRANSIENTQ`+`NPQRT>0` still fail. Only `mfusgqrt.py` touched. 10 new + 1
  repurposed test; `-k mfusgqrt` **47**, focused **282**, combined **286** (ARM).
  See `USGT_STAGE4_13_QRT_PARAMETER_AUTHORING.md`. **This completes the
  list-parameter from-scratch family (DRT/HFB/SGB/QRT).**
- **Stage 4.6C-D — ETS from-scratch *array*-parameter authoring — EXECUTED.**
  `MfUsgEts` authors `NPETS>0` ETSR array parameters from Python
  (`parameters={name: {parval, clusters | instances}}` + `evtr_parm`); a shared
  `build_array_parameter_bc_parms` turns the ergonomic dict into a `ModflowParBc`
  (a different path from the list-parameter helpers). Only ETSR is parameterized;
  first period must activate. Touched `mfusgets.py` + `_usgt_parameters.py`. 6 new
  + 1 repurposed test; `-k mfusgets` **16**, focused **288**, combined **292**
  (ARM). See `USGT_STAGE4_14_ETS_PARAMETER_AUTHORING.md`. **This completes the
  parameter from-scratch authoring family (DRT/HFB/SGB/QRT list + ETS array).**

- **Stage 4.6D — LAK multi-lake + sill/connectivity (ds 7/8) authoring.**
  *Acceptance:* author a 2-lake system with sill connectivity from scratch →
  write → reload; `Ex8` still round-trips/runs.

- **Stage 4.6E — DPT `A-W_ADSORBIM` decision (old Stage 4.6).**
  Either implement the immobile air-water adsorption cascade with from-scratch
  authoring, or ratify the explicit-deferral with a tighter Fortran-derived spec.

- **Stage 4.6F — EVT ETS-zonal + `NPEVT` authoring (old EVT gaps).**

- **Stage 4.6G — P3 cleanup** (DISU/OC todos; load-skip warning; doc numbering).

- **(Separate track) Stage 4.8 — Upstream infrastructure** stays as planned in
  `USGT_STAGE4_08_UPSTREAM_INFRA.md` (executable source/release/CI); it is not a
  FloPy authoring gap and is sequenced after package completeness.

---

## 6. Recommended next card

The **parameter from-scratch authoring family is complete** (A1 closed):
**Stage 4.6A** (STR/SUB/SWT guard; SFR kept), **4.6B** (DRT), **4.6C-A** (HFB),
**4.6C-B** (SGB definitions), **4.6C-C** (QRT, structural), and **4.6C-D** (ETS
array) are all **done**. The recommended next card is **Stage 4.6D — LAK
multi-lake + sill/connectivity (ds 7/8) authoring** (A2, the top remaining P1);
then the P2 items (DPT `A-W_ADSORBIM`, EVT ETS-zonal/`NPEVT`, OC/MDT exe
verification) per §4.
