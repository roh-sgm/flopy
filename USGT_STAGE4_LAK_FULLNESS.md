# LAK Fullness Card D — executed (USG-T 2.7)

Scope: audit `MfUsgLak` against USG-T 2.7 Lake Package, fix authoring bugs, and
decide honestly whether LAK can be promoted to `Full`. This closes the
`USGT_STAGE4_03_OC_EVT_MDT_LAK.md` theme (OC/EVT/MDT already done); OC/EVT/MDT are
not touched here beyond docs.

## Fortran audit (gwf2lak7u1.f)

- **Dataset 1a/1b**: keyword options `TABLEINPUT` and `TRANSPORTBOUNDARY`, then
  `NLAKES ILKCB`. `ILKTRNSPT` is set from `IUNITGWT` (transport); `TRANSPORTBOUNDARY`
  is rejected when transport is inactive (it sets `ILKTRNSPT=0`, "lake as a
  boundary").
- **Dataset 2/3**: `THETA [NSSITR SSCNCR] [SURFDEP]` (extra values when THETA<0 or
  steady); per lake `STAGES [SSMN SSMX]`, then `CLAKE(1:NSOL)` when transport,
  then the tab unit when `TABLEINPUT`.
- **Per stress period**: `ITMP ITMP1 LWRT`; `LKARR`/`BDLKNC` when `ITMP>0`; sill
  data (datasets 7/8); dataset 9a `PRCPLK EVAPLK RNF WTHDRW [SSMN SSMX]` per lake;
  dataset 9b concentrations per lake (when transport):
  - **`TRANSPORTBOUNDARY` (ILKTRNSPT=0):** one line per lake with `CLAKE(1:NSOL)`
    (all components on one line) — `gwf2lak7u1.f:1060-1064`;
  - **classic transport (ILKTRNSPT=1):** per (lake, component) `CPPT, CRNF
    [, CAUG]` (CAUG only when `WTHDRW<0`).

## Bugs fixed (`flopy/mfusg/mfusglak.py`)

1. **conc_data mis-assignment** — a non-dict `conc_data` was wrapped as
   `{0: sill_data}` instead of `{0: conc_data}`.
2. **`write_file` crash when `flux_data` is None** — `list(self.flux_data.keys())`
   raised; LAK now requires `flux_data` (dataset 9) with a clear error, and the
   key access is guarded.
3. **`conc_data` accessed without transport** — `ds9b = self.conc_data[kper]` ran
   even for `mcomp==0` (conc_data None); now read only when `mcomp>0`.
4. **TRANSPORTBOUNDARY dataset-9b shape** — the writer emitted one line per
   component; the Fortran reads **one line per lake** with all NSOL values. The
   writer now matches (and classic transport still writes per (lake, component)).
5. **TRANSPORTBOUNDARY load type** — concentrations were stored as strings; now
   `float`.
6. **`transportboundary` flag vs options** — `transportboundary=True` did not add
   the `TRANSPORTBOUNDARY` keyword to the header, so the file used the boundary
   9b layout without the keyword and did not reload. The flag and `options` are
   now kept in sync (both directions).

## Validation added

`__init__` now validates: `flux_data` is required; `TRANSPORTBOUNDARY` requires
active transport (`mcomp>0`); `clake` must be `nlakes x mcomp`; and active
transport requires `conc_data` (dataset 9b).

## Tests

Synthetic (`autotest/test_usg_transport.py`, `-k mfusglak`):

- `test_mfusglak_minimal_authoring_roundtrip` (no transport),
- `test_mfusglak_classic_transport_roundtrip` (CPPT/CRNF per lake-component),
- `test_mfusglak_transportboundary_roundtrip` (one CLAKE line per lake, MCOMP=2,
  header keyword + float round-trip),
- `test_mfusglak_tableinput_authoring` (TABLEINPUT keyword + per-lake tab unit +
  external registration),
- `test_mfusglak_rejects_invalid` (missing flux_data, TRANSPORTBOUNDARY without
  transport, bad clake length, transport without conc_data).

Real model: the `Ex8_Lake` test (load → write → run under the USG-T 2.7
executable) is kept. No from-scratch LAK executable smoke was added — a minimal
convergent lake model is not cheap/stable to build, and Ex8 already exercises LAK
execution on a real model.

The review follow-up adds five more (see below). `-k mfusglak/Ex8` **11 passed**;
focused **165 passed**; exe **4 passed**; combined **169 passed** under the USG-T
2.7 ARM binary.

## Decision: keep `✅ (intentionally not Full)`, hardened

From-scratch authoring is now real and tested for the main branches (no-transport,
classic transport, `TRANSPORTBOUNDARY` with MCOMP>1, `TABLEINPUT`), six authoring
bugs were fixed, and validation was added. But bounded gaps keep the honest label
at `✅ (intentionally not Full)`:

- **Sill / connectivity (datasets 7/8) and multi-lake sublake systems**: written
  and round-tripped (and exercised by Ex8) but not authored from scratch in the
  tests; complex connectivity authoring is not exhaustively covered.
- **`TABLEINPUT` bathymetry table content**: FloPy registers the per-lake tab
  unit / external file, but the table file *contents* are external and not
  authored here.
- **GAGE coupling** is a separate package.

## Review follow-up

Code review found three from-scratch authoring inputs that were accepted by
`__init__` but then crashed inside `write_file()` with a raw `IndexError` /
`KeyError` / `TypeError` instead of failing with a clear error. All three now
validate in `__init__`:

1. **P1 — TABLEINPUT with too few `tab_files`.** The constructor only built a
   dead `msg` (with a `# TODO`) when `len(tab_files) < nlakes`; `write_file()`
   then crashed with `IndexError` on `self.iunit_tab[n]`. Now `TABLEINPUT`
   requires **exactly one `tab_file` per lake**, and if `tab_units` is supplied
   it must also be **one per lake** — otherwise a clear `ValueError`.
2. **P1 — incomplete / malformed `conc_data`.** A missing `(lake, component)`
   entry (e.g. `mcomp=2` with only `(0,0)`) was accepted and crashed
   `write_file()` with `KeyError`. `__init__` now validates dataset 9b for every
   stress period where dataset 9 (`flux_data`) is written:
   - an entry must exist for **every `(lake, component)`**;
   - **classic transport:** each entry is a sequence of **2 values** (`CPPT,
     CRNF`) when `WTHDRW >= 0` and **3 values** (`CPPT, CRNF, CAUG`) when
     `WTHDRW < 0` (matches the Fortran, which reads CAUG only for augmentation);
   - **TRANSPORTBOUNDARY:** each entry is a **single** concentration (the writer
     still emits one `CLAKE(1:NSOL)` line per lake);
   - any other shape raises a clear `ValueError`, never `KeyError`/`TypeError`.
3. **P2 — partial `flux_data` validation.** Dataset 9a now requires **one entry
   per lake** (a missing lake raises `ValueError` instead of `KeyError`), and the
   per-period length check (≥4 values transient/first period, ≥6 for steady
   periods after the first) raises `ValueError` rather than a bare `Exception`.

Five tests were added (all green): TABLEINPUT with one `tab_file` for two lakes
→ `ValueError`; TABLEINPUT with a wrong-length `tab_units` → `ValueError`; classic
`mcomp=2` with incomplete `conc_data` → `ValueError`; classic `WTHDRW<0` with two
values → `ValueError` and with three values → writes/reloads; TRANSPORTBOUNDARY
`mcomp=2` with incomplete `conc_data` → `ValueError`; and `flux_data` missing a
required lake entry → `ValueError`. The five original positive/negative tests and
the `Ex8_Lake` round-trip are unchanged. Decision is unchanged: LAK stays
`✅ (intentionally not Full)`.

## Stage 4.6D — multi-lake + sill/connectivity (datasets 7/8) authoring

Card D left this gap: "Sill/connectivity (datasets 7/8) and multi-lake sublake
systems: written and round-tripped (and exercised by Ex8) but not authored from
scratch in the tests." **Stage 4.6D closes the authoring half of that gap.**

The writer already emitted datasets 7/8 (only when `ITMP>0`, matching
`gwf2lak7u1.f` `GWF2LAK7U1RPU`: `ITMP<=0` skips datasets 5/6/7/8). What was
missing was validation and from-scratch tests. Added (`flopy/mfusg/mfusglak.py`):

- `_validate_sill_data(sill_data, nper)` — validates datasets 7/8 in `__init__`
  before any file is opened: `kper` in range; the period has `ITMP>0` (else the
  sill would be silently dropped); each system is `(ds8a, sillvt)` with
  `ds8a=[IC, lake1..lakeIC]`, `IC>=2`, `IC==len(lakes)`, lake numbers **1-based**
  in `[1, NLAKES]` and unique, and `len(sillvt)==IC-1`.
- Rewrote the `sill_data` docstring (exact `{kper: [([IC, lake1..lakeIC],
  [sill1..sill_{IC-1}]), ...]}` layout; 1-based lakes, center first; `ITMP>0`
  rule). The old docstring mislabeled `IC` as "the number of sublakes" — it is
  the total number of lakes in the system.

**Convention:** lake numbers in `sill_data` stay **1-based** (as `load` already
stores them and as `LKARR`/lake-IDs use everywhere); not changed (would break the
Ex8 round-trip). Documented in the docstring + validator.

Three tests added: `test_mfusglak_multilake_sill_authoring_roundtrip`
(two lakes, no transport, datasets 7/8 written 1-based + reload preserves
`IC`/`ISUB`/`SILLVT`); `test_mfusglak_multilake_sill_transport_roundtrip` (two
lakes + sill + classic `mcomp=1`); `test_mfusglak_sill_rejects_invalid` (7
negatives, no partial file). `Ex8_Lake` stays green. `-k "mfusglak or Ex8"`
**14**, focused **291**, exe **4**, combined **295** (ARM). See
`USGT_STAGE4_15_LAK_MULTILAKE_CONNECTIVITY.md`.

**Decision: LAK stays `✅ (intentionally not Full)`** — the sill/multi-lake
*authoring* gap is now closed, but the honest gaps that remain are
`TABLEINPUT` bathymetry-table **contents** (external), **GAGE** (separate
package), and multi-lake-with-sill from-scratch **execution** (not exe-smoke-
tested; manual tier — `Ex8_Lake` covers LAK execution on a real model).

## Validation

```bash
python -m pytest autotest/test_usg_transport.py -k "mfusglak or Ex8" -q
python -m pytest autotest/test_usg_transport.py -q
python -m pytest autotest/test_usg_transport_exe.py -q
USGT_EXE=/Users/roh.sgm/Documents/GitHub_Projects/GW-Software-Compiled-via-Claude/gfortran/usgt_2.7/usgt_270_arm python -m pytest autotest/test_usg_transport.py autotest/test_usg_transport_exe.py -q
git diff --check
git status --short
```
