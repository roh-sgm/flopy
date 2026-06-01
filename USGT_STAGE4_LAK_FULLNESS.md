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

`-k mfusglak/Ex8` **6 passed**; focused **160 passed**; exe **4 passed**;
combined **164 passed** under the USG-T 2.7 ARM binary.

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

## Validation

```bash
python -m pytest autotest/test_usg_transport.py -k "mfusglak or Ex8" -q
python -m pytest autotest/test_usg_transport.py -q
python -m pytest autotest/test_usg_transport_exe.py -q
USGT_EXE=/Users/roh.sgm/Documents/GitHub_Projects/GW-Software-Compiled-via-Claude/gfortran/usgt_2.7/usgt_270_arm python -m pytest autotest/test_usg_transport.py autotest/test_usg_transport_exe.py -q
git diff --check
git status --short
```
