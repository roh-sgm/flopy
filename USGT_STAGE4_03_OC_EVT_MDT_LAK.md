# Stage 4.3 - OC / EVT / MDT / LAK Fullness Pass

Goal: promote the remaining project-relevant `OK but not Full` packages only
where Fortran audit, semantic authoring, round-trip tests, and unsupported-mode
behavior justify it.

## Package Order

1. `OC`
2. `EVT`
3. `MDT`
4. `LAK`

This order keeps the smaller output/control packages ahead of large LAK.

## Current State

### OC

Status: implemented and tested, not Full.

Known gap:

- full output-block round-trip and authoring coverage is not exhaustive,
- `BOOTSTRAPPING` is parsed on load but not fully protected by from-scratch
  authoring tests.

### EVT

Status: implemented and tested, not Full.

Known gap:

- transport `IETFACTOR`/`ETFACTOR` authoring exists,
- full per-`NEVTOP` array authoring combinations are not exhaustively tested.

### MDT

Status: implemented and round-trip validated via Ex7, not Full.

Known gap:

- field order is not independently Fortran-audited,
- no from-scratch authoring test.

### LAK

Status: implemented enough for Ex8 load/write, not Full.

Known gap:

- from-scratch authoring of `TABLEINPUT`,
- `TRANSPORTBOUNDARY`,
- complete lake connectivity,
- bathymetry/stage tables,
- transport boundary concentrations.

## Work Plan

Treat each package as its own mini-card. Do not try to promote all four in one
large commit unless the changes are documentation-only.

### Card A - OC Fullness — EXECUTED

Required:

- audit `mfusgoc.py` against USG-T OC behavior,
- add from-scratch `BOOTSTRAPPING` test if supported,
- add output-block coverage for common save/print concentration/head/budget
  combinations,
- decide whether compact/non-compact or rare modes are in scope.

**Done.** Full audit, implementation, tests, and the honest Full/not-Full
decision are recorded in **`USGT_STAGE4_OC_FULLNESS.md`**: implemented the
`BOOTSTRAPPING` header (first-line placement + round-trip) and fixed
layer-qualified / `DDREFERENCE` round-trip; added six `-k mfusgoc` tests; kept
OC `✅ (intentionally not Full)` with explicit gaps (`SAVE IBOUND`
solver-rejected, FASTFORWARD/BOOTSTRAPPING execution not exe-verified).
**(Qualified by Stage 4.7B: `ATSA` is now exe-verified; FASTFORWARD/FASTFORWARDC/
BOOTSTRAPPING stay manual tier — external prior-run head/conc file.)** A polish
follow-up then fixed the DDREFERENCE-only period write, taught `check()` the
USG-T OC actions (single-word toggles vs keyword-value params that require a
value), and added a specific `check()` warning for `SAVE IBOUND` (4 more tests;
`-k mfusgoc` 11 passed). MDT / LAK (Cards C–D) remain not started.

### Card B - EVT Fullness — EXECUTED

Required:

- audit `gwf2evt8u1.f`,
- add `NEVTOP` combination tests,
- add authoring tests with and without transport,
- verify 0-based/1-based behavior where node/layer indices apply,
- protect `IETFACTOR`/`ETFACTOR` reload.

**Done.** Audit, fixes, tests, and the Full/not-Full decision are recorded in
**`USGT_STAGE4_EVT_FULLNESS.md`**: from-scratch authoring + round-trip for all
NEVTOP modes (structured + unstructured `NEVTOP=2` with `MXNDEVT`; `IEVT`
0-based/1-based + layer-range validation) and transport `IETFACTOR` 0/<0/>0 with
per-`MCOMP` `ETFACTOR`; fixed the `IETFACTOR` round-trip bug (it was dropped on
load) and the dataset-1 3-integer write under transport; ETS zonal raises
`NotImplementedError`; added a USG-T 2.7 EVT executable smoke. A review
follow-up then fixed a scalar-ETFACTOR crash and added unstructured `IEVT`
node-range validation. Eight synthetic tests + one exe test; `-k mfusgevt` 9
passed. Kept EVT `✅ (intentionally not Full)` (gaps: ETS zonal, NPEVT parameter
preservation). LAK (Card D) remains not started.

### Card C - MDT Fullness — EXECUTED

Required:

- audit `gwt2mdtu1.for`,
- document exact dataset order and conditional reads,
- add from-scratch MDT authoring test,
- reload written MDT and compare all semantic fields,
- include multi-species if Fortran supports/uses it.

**Done.** Audit, bug fixes, tests, and the Full/not-Full decision are recorded
in **`USGT_STAGE4_MDT_FULLNESS.md`**: Fortran-audited the dataset order +
conditional reads; from-scratch authoring + round-trip for all main branches
(options, base arrays, per-species, AIOLD under TSHIFTMD, IDPF on/off,
multi-species); fixed three real bugs (`load` dropped `FRAHK`/`FRADARCY`,
`write_file(f=)` crash, debug print) and added option/length validation; kept
the Ex7 real-model round-trip/run. A review follow-up then fixed the TSHIFTMD
threshold (use 1e-10 consistently + a general number format so small valid
values are not rounded to 0.0). Nine synthetic tests; `-k mfusgmdt` 9 passed.
Kept MDT `✅ (intentionally not Full)` (gaps: species count uses MCOMP, AI1/AI2
output binaries not read). LAK (Card D) remains not started.

### Card D - LAK Fullness — EXECUTED

Required:

- audit `gwf2lak7u1.f`,
- define a minimal valid synthetic lake model,
- add from-scratch authoring for `TABLEINPUT` and `TRANSPORTBOUNDARY`,
- test load/write/reload,
- only add executable validation if the synthetic LAK model is stable and cheap.

**Done.** Audit, bug fixes, tests, and the Full/not-Full decision are recorded in
**`USGT_STAGE4_LAK_FULLNESS.md`**: Fortran-audited the datasets; fixed six
authoring bugs (conc_data mis-assignment; `flux_data=None` crash; conc_data read
without transport; TRANSPORTBOUNDARY 9b written per-component instead of
one-line-per-lake; load stored conc as strings; `transportboundary` flag not
synced to the header keyword); added validation; from-scratch authoring +
round-trip for no-transport, classic transport, `TRANSPORTBOUNDARY` (MCOMP>1),
and `TABLEINPUT`. Kept the `Ex8_Lake` real-model round-trip/run (no synthetic LAK
exe smoke — Ex8 covers execution). Kept LAK `✅ (intentionally not Full)` (gaps:
sill/connectivity + multi-lake systems not authored from scratch; TABLEINPUT
table contents external; GAGE separate). **(Superseded by Stage 4.6D:
sill/connectivity + multi-lake are now authored from scratch and validated; the
remaining gaps are TABLEINPUT contents, GAGE, and from-scratch multi-lake
execution — see `USGT_STAGE4_15_LAK_MULTILAKE_CONNECTIVITY.md`.)**

**Review follow-up (executed).** Hardened three from-scratch authoring inputs
that `__init__` accepted but `write_file()` then crashed on with a raw
`IndexError`/`KeyError`/`TypeError`: TABLEINPUT now requires exactly one
`tab_file` (and `tab_unit`, if supplied) per lake; `conc_data` is validated per
written period (every `(lake, component)` present; classic entries 2 values when
`WTHDRW>=0` / 3 when `WTHDRW<0`; TRANSPORTBOUNDARY entries a single value);
`flux_data` requires one dataset-9a entry per lake — all now clear `ValueError`s.
Five tests added. `-k mfusglak` 10 passed (5 original + 5 follow-up); `-k
"mfusglak or Ex8"` 11; focused 165; exe 4; combined 169 under the USG-T 2.7 ARM
binary. Decision unchanged. See `USGT_STAGE4_LAK_FULLNESS.md` → "Review
follow-up".

**Theme complete:** OC, EVT, MDT, and LAK (Cards A–D) are all executed.

## Documentation Updates

After each card, update:

- `USGT_roadmap.md`,
- `USGT_PACKAGE_BACKLOG.md`,
- `USGT_improvements.md`.

Only promote a package to `Full` after its own card is complete.

## Validation

```bash
python -m pytest autotest/test_usg_transport.py -q
python -m pytest autotest/test_usg_transport_exe.py -q
USGT_EXE=/Users/roh.sgm/Documents/GitHub_Projects/GW-Software-Compiled-via-Claude/gfortran/usgt_2.7/usgt_270_arm python -m pytest autotest/test_usg_transport.py autotest/test_usg_transport_exe.py -q
git diff --check
git status --short
```

## First Agent Prompt For This Theme

```text
Goal: complete the OC package fullness pass for USG-T 2.7 without touching EVT/MDT/LAK yet.

Read first:
- USGT_STAGE4_MASTER_PLAN.md
- USGT_STAGE4_03_OC_EVT_MDT_LAK.md
- USGT_roadmap.md
- flopy/mfusg/mfusgoc.py
- autotest/test_usg_transport.py
- USG-T 2.7 Fortran source files relevant to OC/output control.

Task:
Audit OC behavior, then add focused semantic authoring/round-trip tests for the missing OC modes, especially BOOTSTRAPPING and output-block combinations. Promote OC only if the tests justify it; otherwise keep the status honest and document exactly what remains.

Validation:
python -m pytest autotest/test_usg_transport.py -q
python -m pytest autotest/test_usg_transport_exe.py -q
USGT_EXE=/Users/roh.sgm/Documents/GitHub_Projects/GW-Software-Compiled-via-Claude/gfortran/usgt_2.7/usgt_270_arm python -m pytest autotest/test_usg_transport.py autotest/test_usg_transport_exe.py -q
git diff --check
git status --short

Constraints:
- Push only to origin, never upstream.
- Do not touch MF6-TID.
- Do not bundle OC, EVT, MDT, and LAK changes in one commit unless docs only.
```
