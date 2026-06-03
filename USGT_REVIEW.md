# USG-T 2.7 FloPy — Review Handoff

Reviewer-facing summary of the USG-Transport 2.7 support work on the `develop`
branch of this fork. Pair this with `USGT_roadmap.md` (coverage table),
`USGT_improvements.md` (change log), and `USGT_PACKAGE_BACKLOG.md` (closeout
status).

**Bottom line:** all backlog tiers (P0–P5) are worked through. The four
highest-value gaps are implemented with parser + writer + from-scratch
authoring tests + round-trip tests, each verified against the USG-T 2.7
Fortran source. **80 tests pass** in `autotest/test_usg_transport.py`
(synthetic authoring/round-trip + `Ex1..Ex9` real-model load+write), ~40 s.

Fortran source of truth (used for every I/O claim below):
`USGT_V_2-7-0_Source_Code/` (USG-Transport 2.7.0, 2026.03.21).

---

## 1. What to review (commit range)

```
git log --oneline b0de152f..HEAD
```

| Commit | Scope |
|---|---|
| `b0de152f` | **Checkpoint of prior work** (P0/P2/P3 authoring already on the branch). Snapshot only — not part of this slice's new implementation. |
| `ad7de360` | **P1**: new `MfUsgSgb`, `MfUsgQrt`, `MfUsgDrt`; BCF/LPF `TABRICH` 1c/1d. |
| `036d714d` | **P2**: BAS `RICHARDS_HP`/`IHM`, DPT `A-W_ADSORBIM` explicit fail, ETS/HFB tests. |
| `6429dcc6` | **P3/P4/P5**: BCT/DDF from-scratch authoring tests + scope decisions. |
| `4de2e0bc` | Docs closeout summary + accept P0 guardrails. |

The new implementation to review is `b0de152f..HEAD`:
**15 files, +2435 / −117 lines, 29 new tests.**

```
git diff --stat b0de152f..HEAD
```

---

## 2. How to verify (run it)

```bash
# from the repo root, in the dev env (numpy/pandas/pytest/modflow_devtools)
pip install -e .

# whole USG-T suite (synthetic + Ex1..Ex9 real-model load+write): expect 80 passed
python -m pytest autotest/test_usg_transport.py -q

# fast synthetic-only subset (no real-model loads): expect 64 passed
python -m pytest autotest/test_usg_transport.py -k "not test_usg_load_Ex" -q

# per-package focus
python -m pytest autotest/test_usg_transport.py -k "sgb"     -q   # 7
python -m pytest autotest/test_usg_transport.py -k "qrt"     -q   # 5
python -m pytest autotest/test_usg_transport.py -k "drt"     -q   # 6
python -m pytest autotest/test_usg_transport.py -k "tabrich" -q   # 5 (incl. DPF)

# imports + NAM registry wired
python -c "from flopy.mfusg import MfUsg, MfUsgSgb, MfUsgQrt, MfUsgDrt; \
m=MfUsg(structured=False); \
print(all(m.mfnam_packages[k] is c for k,c in \
[('sgb',MfUsgSgb),('qrt',MfUsgQrt),('drt',MfUsgDrt)]))"   # True
```

---

## 3. New / changed code

**New modules**

| File | Purpose |
|---|---|
| `flopy/mfusg/mfusgsgb.py` | `MfUsgSgb` — Specified Gradient Boundary |
| `flopy/mfusg/mfusgqrt.py` | `MfUsgQrt` — Sink with Return Flow |
| `flopy/mfusg/mfusgdrt.py` | `MfUsgDrt` — Drain Return (DRT8), replaces base `ModflowDrt` |
| `flopy/mfusg/_usgt_returnflow.py` | shared `U1DINT` recipient-node read/write (DRT + QRT) |
| `flopy/mfusg/_tabrich.py` | shared TABRICH 1c/1d read/write (BCF + LPF) |

**Changed modules**

| File | Change |
|---|---|
| `flopy/mfusg/mfusg.py` | registry: add `"sgb"`, `"qrt"`; repoint `"drt"` → `MfUsgDrt` |
| `flopy/mfusg/__init__.py` | export `MfUsgSgb`, `MfUsgQrt`, `MfUsgDrt` |
| `flopy/mfusg/mfusgbcf.py` | TABRICH `iuzontab`/`retcrvs` constructor args + write/load |
| `flopy/mfusg/mfusglpf.py` | TABRICH 1c/1d; skip per-layer Richards arrays under TABRICH; fix parse bug (`nutabrows` token index + `int` cast) |
| `flopy/mfusg/mfusgbas.py` | `RICHARDS_HP`, `IHM [IUIHM]`; keep `_` in option-line cleaner |
| `flopy/mfusg/mfusgdpt.py` | explicit fail on `A-W_ADSORBIM` |

---

## 4. Fortran cross-reference (verify I/O against these)

| Package | Fortran file / subroutines | Record / dataset layout (unstructured, free format) |
|---|---|---|
| SGB | `glo2sgbu1.f` — `GLO2SGBU1AR/RP`, `ULSTRDU` | Header `MXACTS ISGBCB [AUX..][NOPRINT]`; per SP `ITMP NP`, rows `NODE GRADIENT [aux]`. `NSGBVL=5+NAUX`. |
| QRT | `gwf2QRT8u.f` — `GWF2QRT8U1AR/RP`, `SGWF2QRT8LR` | Header `MXAQRT MXRTCELLS IQRTCB NPQRT MXL`; line `ND Q NumRT [Rfprop] [IQCHNGTYP] [aux]`; recipient nodes (`NodQRT`) via `U1DINT` **after all sink lines**. `NQRTVL=5+NAUX+IQRTFL`. |
| DRT | `gwf2drt8u.f` — `GWF2DRT8U1AR/RP`, `SGWF2DRT8LR` | Header `MXADRT IDRTCB NPDRT MXL`; line `ND EL COND NR [PROP] [IDCHNGTYP] [aux]`. `NR>0` single inline recipient; `NR<0` ⇒ `-NumRT` spreading nodes via `U1DINT` **immediately after that line** (interleaved). `NDRTVL=5+NAUX+2+IDRTFL`. |
| BCF/LPF TABRICH | `gwf2bcf-lpf-u1.f` (BCF block ~L163–180, LPF block ~L3120) | After the option line, before LAYCON/LAYTYP: 1c `IUZONTAB` = `U1DINT(NODES)`; 1d `RETCRVS(3,NUTABROWS,NUZONES)` read zone-outer/row-inner, 3 free values per line = **capillary head / saturation / relative permeability**. Python shape `(nuzones, nutabrows, 3)`. |
| BAS | `glo2basu1.f` L161 (`RICHARDS_HP` ⇒ IUNSat=1, IPRES=1), L184–198 (`IHM` then `IUIHM`) | option keywords on the BAS option line |
| DPT | `gwt2dptu1.f` L135 (`A-W_ADSORBIM`), `dpt2aw_adsorb.f` | option triggers extra arrays — array-only branches implemented (Stage 4.6E); scalar/tabular branches explicit-fail |

Indexing contract (verify in the writers/loaders): internal `node` and
recipient values are **0-based**; files are read/written **1-based**.
`ITMP < 0` reuses the previous stress period.

---

## 5. Design decisions & explicit-failure scope cuts

Per the USG-T design rules, unsupported sub-modes **fail explicitly** rather
than write incomplete / mis-parse.

> **Updated (Stages 4.4–4.6):** several of the *original* scope-cuts listed below
> have since been closed — `NPSGB`/`NPQRT`/`NPDRT`/`NPHFB`/`NPETS`/`NPEVT` named
> parameters are now preserved **and** authored from scratch (Stages 4.4 /
> 4.6B–4.6F-A); QRT inline `TRANSIENTQ` (4.5A) and DRT/QRT recipient
> `EXTERNAL`/`OPEN-CLOSE` (4.5B, free-format) are supported; DPT `A-W_ADSORBIM`
> array-only branches are implemented (4.6E). The list below is the original
> scope; the **authoritative current state** is `USGT_roadmap.md` /
> `USGT_STAGE4_09_FINAL_GAP_AUDIT.md`. The explicit failures that *remain* are the
> per-package caveats recorded there: parameter `INSTANCES` (list packages), SGB
> active params (Fortran `SGB`/`G` conflict), `TRANSIENT_HFB`+`NPHFB>0`,
> `TRANSIENTQ`+`NPQRT>0` / external-unit `TRANSIENTQ`, DPT `A-W_ADSORBIM`
> scalar/tabular branches (`IAREA_FNIM ∈ {2,3,5}`, `IKAWI_FNIM ∈ {3,4}`), EVT
> ETS-zonal time-series (`ETS MXZNEVT`), and STR/SUB/SWT on unstructured `MfUsg`.

Original scope-cuts (the design intent; see the note above for current state):

- Named parameters: `NPSGB>0`, `NPQRT>0`, `NPDRT>0` → `NotImplementedError`.
- QRT `TRANSIENTQ` (transient-flow time series) → `NotImplementedError`.
- DRT/QRT recipient lists with `EXTERNAL`/`OPEN/CLOSE` → `NotImplementedError`.
- TABRICH write without both `iuzontab` and `retcrvs` → `ValueError`
  (construction is allowed so BCF/LPF can still serve as model setup).
- HFB `NPHFB>0` → `NotImplementedError` (load + write).
- DPT `A-W_ADSORBIM` → `NotImplementedError` (load).

Intentionally **Partial** / not `Full`, documented in the roadmap. Current
reasons (post-4.6): **ETS** — only ETSR is parameterizable (Fortran limit) and
parametric execution is not exe-smoke-tested (preservation + from-scratch
authoring done); **HFB** — `TRANSIENT_HFB`+`NPHFB>0` and parameter `INSTANCES`;
**DPT** — `A-W_ADSORBIM` scalar/tabular branches (array-only done); **LAK** —
TABLEINPUT table contents external, GAGE separate, multi-lake-with-sill
from-scratch execution not exe-smoke-tested (authoring done); **EVT** — ETS-zonal
time-series (ATS-coupled). **SFR/FHB/GAGE** are compatibility base classes;
**STR/SUB/SWT** are compatibility-only **guarded** on unstructured `MfUsg`
(Stage 4.6A).

`MfUsgDrt`: new class subclassing `ModflowDrt`; **structured (DIS) grids
delegate to the base class** (USG-T extensions are unstructured-only).

---

## 6. New tests (29)

```
git diff b0de152f..HEAD -- autotest/test_usg_transport.py | grep '^+def test'
```

- **SGB (5):** roundtrip, aux roundtrip, programmatic authoring (0-based→1-based),
  NAM registry, parameters-fail.
- **QRT (5):** minimal authoring, return-flow concentration (CHANGEC+AUX),
  multi-recipient + pure-sink, NAM registry, unsupported-modes-fail (NPQRT/TRANSIENTQ).
- **DRT (6):** inline single recipient, CHANGEC concentration, SPREAD multi-node,
  stress-period reuse, NAM registry, parameters-fail.
- **TABRICH (4):** BCF authoring+roundtrip, LPF authoring+roundtrip,
  incomplete-write-fails, RETCRVS shape validation.
- **P2 (5):** BAS RICHARDS_HP+IHM, ETS NETSOP=2, ETS IESFACTOR,
  HFB NPHFB>0 fail, DPT A-W_ADSORBIM fail.
- **P3 (4):** BCT minimal/IDISP=2/multi-species authoring, DDF NONLINEAR table.

Each new package has the backlog-required pair: a **from-scratch authoring**
test (build from Python, write, inspect the file for 1-based ids) **and** a
**round-trip** test (`load → write → load`, nodes back to 0-based), plus a
NAM-registry test and explicit-failure tests.

---

## 7. Suggested review focus (highest risk first)

1. **Indexing** in each new writer/loader: `node + 1` on write, `node - 1` on
   load; `ITMP < 0` reuse. One off-by-one silently moves the target cell.
2. **QRT vs DRT recipient ordering**: QRT reads all recipient `U1DINT` blocks
   *after* the sink list; DRT reads each spreading block *immediately after*
   its drain line. Confirm against `SGWF2QRT8LR` / `SGWF2DRT8LR`.
3. **TABRICH placement**: 1c/1d are read *before* LAYCON/LAYTYP, and the LPF
   per-layer Richards arrays must be *skipped* under TABRICH.
4. **DRT registry repoint**: `"drt"` now maps to `MfUsgDrt`; confirm structured
   models still load via the base class (no `Ex*` model uses DRT, so the
   structured path is delegation-only).
5. **BAS option-line cleaner** keeping `_`: confirm it doesn't change parsing of
   any other option (`RICHARDS_HP` is the only `_` keyword).

---

## 8. Known non-issues

- Untracked and **deliberately excluded** from all commits: `Mf6Tid_demo.ipynb`,
  `SPEC_MF6_TID.md`, `autotest/test_mf6_tid.py`, `flopy/mf6/tid.py` (a separate
  MF6-TID effort), and a generated `*.CBB` binary.
- `ruff` reports 2 pre-existing sort warnings in `flopy/mfusg/__init__.py`
  (`I001`/`RUF022`) and 4 pre-existing `E501` in the test file — **all predate
  this work** (the new modules and new test lines are ruff-clean). Left
  untouched per surgical-change policy; can be auto-fixed separately.
