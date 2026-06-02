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
2. **No silent-wrong-write P0 remains for in-scope packages.** The project's
   consistent *fail-explicitly* discipline (60+ `NotImplementedError` / writer
   `ValueError` sites across `flopy/mfusg`) means unsupported modes refuse to
   write rather than emit a partial/misread file. The single **latent** P0 is the
   six base-class compatibility packages (§4, W1): they are *documented* unused,
   not *guarded*.
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
| ETS | FA (non-param) + **ParamPreserve** (ETSR array) | From-scratch param authoring → `NotImplementedError`; opt-in expand fallback. |
| HFB | FA (non-param) + **ParamPreserve** (list) | Gaps: from-scratch param authoring, `TRANSIENT_HFB`+`NPHFB>0`, `INSTANCES`. |
| DRT | FA + **ParamPreserve** (NPDRT incl. activations + recipients) + **Expanded** (list/recip controls) | Type-consistent; activated SPREAD = structural round-trip only (Fortran copies `DRTF` not `NodDRT`). Unstructured only. |
| SGB | FA + **Param-definition-preserve** (no activations) + **Expanded** | Active SGB params abort the Fortran (`PARTYP='SGB'` vs `'G'`); `NP>0`/`INSTANCES` → `NotImplementedError`. |
| QRT | FA + **ParamPreserve** (structural, not exec-guaranteed) + **Expanded** + `TRANSIENTQ` | Active param value scales `NumRT` not `Q` (Fortran bug); `NodQRT` not copied. From-scratch param authoring → `NotImplementedError`. |
| DPT | FS-exec **except** `A-W_ADSORBIM` (explicit fail) | Immobile air-water adsorption unsupported (rare PFAS sub-mode). |
| SFR, STR, GAGE, FHB, SUB, SWT | **Compat** | Base MODFLOW-2005 classes; USG-T unstructured format not validated. FHB/GAGE round-trip via base in `Ex8`; SFR/STR/SUB/SWT used by no target model. |
| HOB/DROB/RVOB/GBOB/CHOB/hyd, MNW1/2, UZF, GMG/PCG/SIP/DE4, LVDA/KDEP, SYF, lmt6, gwm, PATH/PTH | Intentionally absent | Observations/HUF/coupling/solvers; documented out of scope. |
| EVS, RTS, RES, IBS | n/a (dead) | `CUNIT` placeholders with no functional reader in USG-T 2.7. |

---

## 3. Gap inventory by category

### 3.1 Authoring-from-scratch gaps (primary objective)

- **A1 — From-scratch list/array *parameter* authoring** (ETS, HFB, DRT, SGB,
  QRT). All five preserve loaded `NP>0` definitions but cannot author them from
  Python; each raises `NotImplementedError`. This is the single biggest
  remaining "author a common feature from zero" gap.
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
- **L2 — Compat packages load via base classes** (SFR/STR/SUB/SWT); the USG-T
  unstructured item layout is unvalidated, so a load could misread *if used*.
  FHB/GAGE are proven via the `Ex8` real model.

### 3.3 Writer gaps

- **W1 — Compat packages write via base classes** (SFR/STR/SUB/SWT). If a USG-T
  model used them, the base writer could emit a structured-format / semantically
  different file **without warning** — the one latent P0. Mitigated only by "no
  target/example model uses them"; not enforced by a guard.
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

- **(Latent) W1 — compat packages (SFR/STR/SUB/SWT)** via base MODFLOW-2005
  classes could silently write a non-USG-T file. **In practice P2** because no
  target/example model uses them, but it is the only place a USG-T model *could*
  get a silently-wrong file. Closing it is cheap (a guard) and removes the last
  theoretical P0 — recommended as Stage 4.6A.

  *No other P0s.* The explicit-failure discipline (unsupported modes raise before
  opening the file; writers validate and refuse partial output) has held through
  Stages 3–4.5.

### P1 — blocks from-scratch authoring of common features

- **A1 — from-scratch list-parameter authoring** for the type-consistent
  packages (DRT first, then HFB; QRT/SGB structural). The biggest authoring gap.
- **A2 — LAK multi-lake + sill/connectivity (ds 7/8) authoring.**

### P2 — incomplete preservation / missing tests

- **A3 — DPT `A-W_ADSORBIM`** (implement or formally keep deferred — old 4.6).
- **A4 — EVT ETS-zonal time series + `NPEVT` param authoring.**
- **A5 — HFB `TRANSIENT_HFB`+`NPHFB>0`; `INSTANCES` (all param packages).**
- **L2 — compat-package load validation** (or formal raw/text decision).
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

- **Stage 4.6A — Compatibility-package USG-T guard (closes the latent P0).**
  Make `MfUsgSfr`/`MfUsgStr`/`MfUsgSub`/`MfUsgSwt` *use* on an **unstructured**
  USG-T model fail or warn explicitly instead of silently using the structured
  base format. Keep FHB/GAGE working (proven via `Ex8`). Tiny, high-certainty.
  *Acceptance:* loading/authoring one of these on a DISU model raises a clear
  `NotImplementedError` (or documented guard warning); `Ex8` still round-trips;
  one negative test; roadmap "compatibility-only" rows updated to "guarded".

- **Stage 4.6B — From-scratch list-parameter authoring: DRT (pilot).**
  Extend `_usgt_parameters.py` + `MfUsgDrt` so a parameterized package can be
  built in Python (`parameters=`, `active_params=`) and written to a valid
  `NPDRT>0` file. DRT first because it is type-consistent and execution-valid.
  *Acceptance:* author `NPDRT>0` from scratch → write → reload equality; exe-run
  under `USGT_EXE`; `INSTANCES` still explicit-fail; ETS/HFB/SGB/QRT unchanged.

- **Stage 4.6C — Extend from-scratch param authoring to HFB, then SGB(defs)/QRT.**
  Reuse the 4.6B shared path. SGB stays definition-only (active params abort the
  Fortran); QRT stays structural (documented Fortran value-scaling bug).
  *Acceptance:* per-package from-scratch authoring tests; honest labels kept.

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

**Stage 4.6A — Compatibility-package USG-T guard.** It closes the only latent P0
(the sole place a USG-T model could be written silently wrong), is small and
low-risk, needs no Fortran reimplementation, and makes the existing
"compatibility-only" documentation *enforceable*. It should land before the
larger P1 authoring work (4.6B+), consistent with "P0 before P1".

If the reviewer prefers to lead with authoring impact instead, **Stage 4.6B**
(DRT from-scratch parameters) is the highest-value P1 and a clean pilot for the
shared parameter-authoring path that 4.6C reuses.
