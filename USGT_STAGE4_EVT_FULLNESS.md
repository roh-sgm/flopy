# EVT Fullness Card B — executed (USG-T 2.7)

Scope: audit `MfUsgEvt` against USG-T 2.7 Evapotranspiration and decide honestly
whether EVT can be promoted to `Full`. This is the EVT slice of the
`USGT_STAGE4_03_OC_EVT_MDT_LAK.md` pass; OC / MDT / LAK are **not** touched here
(beyond docs).

## Fortran audit (gwf2evt8u1.f)

`GWF2EVT8U1AR` (dataset 1):

- without transport: `NEVTOP IEVTCB`;
- with transport (`INBCT>0`): `NEVTOP IEVTCB IETFACTOR` — always 3 integers;
- optional `ETS MXZNEVT` keyword on the same line (requires ATS, else STOP);
- `NEVTOP` must be 1..3 (else abort);
- unstructured + `NEVTOP=2`: a `MXNDEVT` line follows;
- transport `ETFACTOR(MCOMP)`: `IETFACTOR==0` → 0.0, `IETFACTOR<0` → 1.0,
  `IETFACTOR>0` → read MCOMP values on a separate line;
- `NPEVT>0`: named parameters (`UPARARRAL`/`UPARARRRP`).

`GWF2EVT8U1RP` (per stress period): `INSURF INEVTR INEXDP [INIEVT if NEVTOP=2]`
(plus an optional `INEVTZONES` keyword). A negative flag reuses the previous
array. `SURF`/`EVTR`/`EXDP` via `U2DREL`; `IEVT` via `U2DINT` — structured: a
**1-based layer index** validated in 1..NLAY; unstructured: a **1-based node
index** validated ≤ NODES.

## What was fixed / implemented (`flopy/mfusg/mfusgevt.py`)

1. **IETFACTOR round-trip bug** — `load` read `IETFACTOR` but never passed it to
   the constructor, so it reset to 0 on reload (and the `ETFACTOR` array would be
   dropped on the next write). Now preserved.
2. **`write_file` IETFACTOR** — written whenever transport is active (was only
   when `!= 0`), matching the Fortran's 3-integer dataset-1 read for `INBCT>0`.
   The `ETFACTOR` array is now written when `ietfactor > 0` (was `== 1`).
3. **Validation** — `NEVTOP` must be 1..3 (`ValueError`); `ETFACTOR` must have
   `MCOMP` values when `ietfactor>0` (`ValueError`); `NEVTOP=2` structured `IEVT`
   (0-based layer) must be in `[0, NLAY-1]` on write (`ValueError`).
4. **ETS zonal time-series** — `ETS MXZNEVT`/per-SP `IZNEVT` is **not** authored
   or parsed, so `mxetzones>0` now raises `NotImplementedError` rather than
   writing/loading an incomplete file.

0-based internal / 1-based file indexing for `NEVTOP=2 IEVT` (structured layer
and unstructured node) is preserved, and per-SP reuse (`-1`) works through
`Transient2d`.

## Tests

Synthetic (`autotest/test_usg_transport.py`, `-k mfusgevt`):

- `test_mfusgevt_nevtop1_and_3_authoring_roundtrip`
- `test_mfusgevt_nevtop2_structured_roundtrip` (0-based ↔ 1-based)
- `test_mfusgevt_nevtop2_unstructured_mxndevt`
- `test_mfusgevt_multi_period_reuse`
- `test_mfusgevt_transport_ietfactor_roundtrip` (0 / <0 / >0, MCOMP>1)
- `test_mfusgevt_rejects_invalid` (NEVTOP, ETS, ETFACTOR length, IEVT range)
- `test_mfusgevt_transport_etfactor_authoring` (existing).

Executable (`autotest/test_usg_transport_exe.py`, USG-T 2.7):

- `test_usgt_exe_evt_from_scratch` — a from-scratch `NEVTOP=1` EVT on the steady
  CHD flow model runs to normal termination, reports an `ET` budget term, and
  closes. (A cheap acceptance smoke; EVT-perturbed heads are not asserted
  analytically.)

`-k mfusgevt` **7 passed**; focused **144 passed**; exe **4 passed**; combined
**148 passed** under the USG-T 2.7 ARM binary.

## Decision: keep `✅ (intentionally not Full)`, hardened

Core EVT authoring is now broad and round-trip tested (all three `NEVTOP`
options, structured and unstructured `NEVTOP=2`, transport `IETFACTOR` 0/<0/>0
with per-`MCOMP` `ETFACTOR`, per-SP reuse), a real `IETFACTOR` round-trip bug was
fixed, and the executable accepts a from-scratch EVT. But two bounded gaps keep
EVT honest at `✅ (intentionally not Full)`:

- **ETS zonal time-series** (`ETS MXZNEVT` / per-SP `IZNEVT`): unsupported —
  raises `NotImplementedError`. (Requires ATS and per-SP zone arrays that are not
  authored or parsed.)
- **NPEVT named parameters**: a parameterized EVT file loads as expanded arrays
  (`parameter_bcfill`) and is rewritten non-parametrically (`NP=0`) — i.e.
  *Expanded valid write*, the same treatment as the ETS package. Authoring with
  parameters is not supported.

These mirror how the sibling ETS package is classified, so promoting EVT to
`Full` would overclaim.

## Review follow-up (resolved)

1. **Scalar ETFACTOR no longer crashes (P1).** `MfUsgEvt(..., ietfactor=1,
   etfactor=2.5)` previously passed the length check (`np.atleast_1d(2.5).size
   == 1`) but crashed in `write_file` at `self.etfactor[icomp]`. `__init__` now
   normalizes `self.etfactor` to a 1-D array, so a scalar (MCOMP=1) and an array
   (MCOMP>1) both author/load/reload correctly; an incompatible length still
   raises `ValueError`.
2. **Unstructured IEVT node range validated (P2).** For unstructured `NEVTOP=2`,
   `IEVT` is a 0-based node index written 1-based. `write_file` now validates it
   against the total node count (via `node_count`, i.e. DISU `nodes` or the
   DIS fallback): a value `< 0` or `>= NODES` raises a clear `ValueError`.
   Structured `NEVTOP=2` continues to validate the layer index in `[0, NLAY-1]`.

Tests added: `test_mfusgevt_etfactor_scalar_and_array_roundtrip`,
`test_mfusgevt_nevtop2_unstructured_ievt_out_of_range`. `-k mfusgevt`
**9 passed**; focused **146 passed**; exe **4 passed**; combined **150 passed**
under the USG-T 2.7 ARM binary. EVT status unchanged: `✅ (intentionally not
Full)`.

## Validation

```bash
python -m pytest autotest/test_usg_transport.py -k mfusgevt -q
python -m pytest autotest/test_usg_transport.py -q
python -m pytest autotest/test_usg_transport_exe.py -q
USGT_EXE=/Users/roh.sgm/Documents/GitHub_Projects/GW-Software-Compiled-via-Claude/gfortran/usgt_2.7/usgt_270_arm python -m pytest autotest/test_usg_transport.py autotest/test_usg_transport_exe.py -q
git diff --check
git status --short
```
