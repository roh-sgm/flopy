# Stage 4.6D — LAK multi-lake + sill/connectivity (datasets 7/8) authoring

Date: 2026-06-02 · Base: `develop` @ `be2102cc`

Goal: let a user build a multi-lake `MfUsgLak` **from scratch** with the
connected-lake (sublake / sill) systems — datasets 7/8 — that Card D
(`USGT_STAGE4_LAK_FULLNESS.md`) left as "written and round-tripped (via Ex8) but
not authored from scratch / validated". The writer already emitted datasets 7/8;
this stage adds **validation** (so authoring fails with an actionable
`ValueError` before any file is opened) and **from-scratch authoring tests**
(multi-lake, with and without transport), without regressing the Ex8 real-model
round-trip/run.

## Fortran audit (`gwf2lak7u1.f`)

Three read/prepare subroutines: `GWF2LAK7U1RP` (dispatch), `GWF2LAK7U1RPS`
(structured), and **`GWF2LAK7U1RPU`** (unstructured / DISU — the USG-T path).
Datasets 7/8 are read identically in RPS (`:986-1024`) and RPU (`:1599-1638`).

### Datasets 7/8 — "linkage parameters for coalescing lakes"

The Fortran comment calls these *connected lake systems* (sublakes that coalesce
with a center lake once the stage rises above a sill):

- **Dataset 7** — `NSLMS` (`:1603/1605`): number of connected-lake systems. If
  `NSLMS <= 0`, datasets 8 are skipped (`GO TO 760`).
- **Dataset 8a** — per system `IS = 1..NSLMS`: `IC, (ISUB(IS,I), I=1,IC)`
  (`:1613/1615`). `IC` is the **number of lakes** in the system; `ISUB(IS,1)` is
  the **center lake** and `ISUB(IS,2..IC)` are its sublakes. `IF(IC.LE.0) GO TO
  750` — `IC<=0` is the **end-of-list sentinel**, so a real system has `IC>=2`.
- **Dataset 8b** — `(SILLVT(IS,I), I=1,IC-1)` (`:1622/1624`): `IC-1` sill
  elevations, one per sublake in `ISUB` order. Used at `:2887-2928`: when both
  the center and a sublake stage exceed `SILLVT`, the lakes coalesce.

Lake numbers in `ISUB` are **1-based** (as in the file) — the same convention as
`LKARR` (the lake-ID grid array, where 0 = no lake).

### Read condition (critical for the writer)

Per stress period the package reads dataset 4 `ITMP ITMP1 LWRT` (`:1286/1288`),
then:

- `IF(ITMP.GE.0) GO TO 50` then `IF(ITMP.EQ.0) GOTO 900` (`:1293`, `:1302`):
  **`ITMP <= 0` skips datasets 5/6/7/8** (and, for `ITMP==0`, dataset 9 too —
  full reuse of the previous period).
- `ITMP > 0` reads dataset 5 (`LKARR` via `U1DINT`, `:1320`), dataset 6
  (`BDLKNC` via `U1DREL`, `:1350`), **then datasets 7/8 (sill)**, then dataset 9.

So **datasets 7/8 are read only when `ITMP > 0`** — the same condition under which
`LKARR`/`BDLKNC` are (re)specified. FloPy's writer already nests the sill block
inside `if itmp > 0` (`mfusglak.py`), so the new validation enforces the matching
rule on input: `sill_data` for a period with `ITMP<=0` is rejected (it would
otherwise be silently dropped).

### Single-lake vs multi-lake / TABLEINPUT / transport

- `NLAKES` (dataset 1b) drives the per-lake loops: dataset 3 (`STAGES [SSMN SSMX]
  [CLAKE] [LAKTAB]`) and dataset 9a (`PRCPLK EVAPLK RNF WTHDRW [SSMN SSMX]`) are
  read once per lake. Multi-lake is just `NLAKES > 1`; datasets 7/8 only make
  sense with `NLAKES > 1` (they connect lakes).
- `TABLEINPUT` (dataset 1a keyword, `IRDTAB>0`): dataset 3 also reads `LAKTAB(LM)`
  (a per-lake bathymetry-table unit). The table *contents* are an external file —
  out of scope (FloPy registers the unit/file only).
- `TRANSPORTBOUNDARY` (`ILKTRNSPT=0`): dataset 9b is one `CLAKE(1:NSOL)` line per
  lake. **Classic transport** (`ILKTRNSPT=1`): dataset 9b is `CPPT, CRNF [, CAUG]`
  per (lake, component), `CAUG` only when `WTHDRW < 0`. (Both unchanged from Card
  D; this stage exercises classic transport together with datasets 7/8.)

## What changed (code, `flopy/mfusg/mfusglak.py`)

No format change — datasets 7/8 were already written/read. Added:

- **`_validate_sill_data(sill_data, nper)`** (called from `__init__`, before any
  file is opened). For each `kper → [(ds8a, sillvt), ...]` system it checks:
  `kper` in `[0, nper)`; the period has `ITMP>0` (`lakarr` (re)specified) else the
  sill would be dropped; each system is a `(ds8a, sillvt)` pair; `ds8a =
  [IC, lake1..lakeIC]` with `IC>=2` and `IC == len(lakes)`; every lake number in
  `[1, NLAKES]` (1-based) and unique within the system; `len(sillvt) == IC-1`.
- Replaced the old silent normalization (a bare list → `{0: ...}` only, no
  validation) with normalization **+** `_validate_sill_data`.
- Rewrote the `sill_data` docstring to state the exact layout
  (`{kper: [([IC, lake1..lakeIC], [sill1..sill_{IC-1}]), ...]}`), the **1-based**
  lake-number convention, that the center lake is listed first, and the `ITMP>0`
  rule. (The old docstring mislabeled `IC` as "the number of sublakes"; it is the
  total number of lakes in the system.)

### Convention decision

`sill_data` keeps lake numbers **1-based** (as in the file and as `load` already
stores them), *not* 0-based. This matches `LKARR`/lake-IDs everywhere in LAK and
preserves the Ex8 round-trip; changing it would be a silent breaking change. The
docstring and validator state this explicitly.

## Validations (all before the file is opened — no partial file)

| Condition | Result |
|---|---|
| `sill_data` key (`kper`) not in `[0, nper)` | `ValueError` |
| `sill_data` for a period with `ITMP<=0` (no `lakarr`) | `ValueError` |
| a system is not a `(ds8a, sillvt)` pair | `ValueError` |
| `IC < 2` (no sublake; `IC<=0` is the Fortran sentinel) | `ValueError` |
| `IC != len(listed lake numbers)` | `ValueError` |
| a lake number out of `[1, NLAKES]` (1-based) | `ValueError` |
| a lake number repeated within a system | `ValueError` |
| `len(sillvt) != IC - 1` | `ValueError` |

(The Card-D `flux_data`/`conc_data`/`TABLEINPUT`/`TRANSPORTBOUNDARY` validations
are unchanged.)

## Tests (`autotest/test_usg_transport.py`)

New (`-k mfusglak`):

- `test_mfusglak_multilake_sill_authoring_roundtrip` — two lakes, no transport,
  one connected system `([2, 1, 2], [95.0])`; asserts the written `NSLMS=1` /
  `8a = 2 1 2` (1-based) and that load → reload preserves `IC`, `ISUB`, `SILLVT`.
- `test_mfusglak_multilake_sill_transport_roundtrip` — two lakes + sill with
  classic transport (`mcomp=1`); datasets 7/8 and dataset 9b (`CPPT, CRNF` per
  lake) both round-trip.
- `test_mfusglak_sill_rejects_invalid` — 7 negatives (bad lake id, IC/list
  mismatch, `IC<2`, wrong sill count, duplicate lake, `kper` out of range, sill
  for an `ITMP<=0` period) all raise `ValueError` and leave no partial file.

Unchanged: the five Card-D positive/negative LAK tests, the Card-D review
follow-up tests, and the **`Ex8_Lake` real-model round-trip/run** (kept green
under the USG-T 2.7 ARM binary).

### No multi-lake from-scratch executable smoke

A convergent multi-lake-with-sill model built from scratch is not cheap or stable
to set up (the sill logic only engages once both stages exceed `SILLVT`, so an
unbalanced synthetic model either fails to converge or never exercises
coalescence). `Ex8_Lake` already exercises LAK execution on a real model, and the
writer is audited line-by-line against `gwf2lak7u1.f` and round-trips. Multi-lake
from-scratch *execution* therefore stays a documented **manual / exe tier**, same
as the Card-D decision.

## Results

| Suite | Result |
|---|---|
| `-k "mfusglak or Ex8"` | **14 passed** |
| focused (`test_usg_transport.py`) | **291 passed** |
| exe (`test_usg_transport_exe.py`) | **4 passed** |
| combined (ARM, incl. `Ex8_Lake` run) | **295 passed** |
| `git diff --check` | clean |

## Status after Stage 4.6D

LAK stays **`✅ (intentionally not Full)`**, but the **sill/connectivity +
multi-lake from-scratch authoring** gap from Card D is **closed** (real authoring,
validated, with tests). The honest gaps that keep it not-`Full`:

- **`TABLEINPUT` bathymetry-table contents** are external (FloPy registers the
  per-lake unit/file only).
- **GAGE** lake gauging is a separate package.
- **multi-lake-with-sill from-scratch *execution*** is not exe-smoke-tested here
  (manual tier; `Ex8_Lake` covers LAK execution on a real model).

## Validation

```bash
python -m pytest autotest/test_usg_transport.py -k "mfusglak or Ex8" -q   # 14
python -m pytest autotest/test_usg_transport.py -q                        # 291
python -m pytest autotest/test_usg_transport_exe.py -q                    # 4
USGT_EXE=.../usgt_270_arm python -m pytest autotest/test_usg_transport.py autotest/test_usg_transport_exe.py -q   # 295
git diff --check
```
