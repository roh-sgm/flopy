# USGT Agent Brief

## Summary

Goal: make this FloPy fork's USG-Transport support more upstream-ready by
correcting design drift, closing verified gaps against USG-T 2.7 Fortran, and
expanding autotests for existing and missing features.

Primary source of truth for I/O and behavior:

`/Users/roh.sgm/Documents/SGM Co/GDM/01_Software/GSI/USG-Transport 2.7.0 (2026.03.21)/USGT_V_2-7-0_Source_Code`

The USG-T 2.7 I/O manual is secondary. If the manual and Fortran disagree, use
the Fortran.

## Design Criteria

- Keep FloPy conventions: internal indices are 0-based (`node`, `k/i/j`,
  `node1/node2`); USG-T files are read and written 1-based.
- Do not change global FloPy behavior without clear regression tests; be careful
  around NAM files, external files, and `MfList`.
- Treat programmatic authoring as a first-class requirement. A semantic package
  must be constructible from Python/numpy inputs and able to write a valid USG-T
  file without first loading an existing input file.
- Classify each package honestly as `Full semantic`, `Expanded valid write`, or
  `Raw/text round-trip`. Do not mark a package full without parser and tests.
- For MODFLOW parameters, v1 may load parameters and write expanded
  non-parametric arrays. If parameter preservation is not supported, fail
  explicitly when asked to preserve/write parametric syntax.
- Keep round-trip tests separate from programmatic-creation tests. New packages
  need both.
- Update `USGT_improvements.md` and `USGT_roadmap.md` whenever scope or status
  changes.

## Implemented In This Slice

- NAM: preserve subdirectories for external input `DATA` files; continue rebasing
  output `DATA`/`DATA(BINARY)` entries to basenames.
- CHD/RIV/GHB/DRN: use internal 0-based `node` values and write USG-T files
  1-based.
- HFB: match `gwf2hfb7u1.f` for non-parametric static and `TRANSIENT_HFB`
  layouts. `IHFBRD` is treated as a read/reuse flag, not a row count. `NPHFB>0`
  fails explicitly.
- ETS: programmatic `NPETS>0` write now fails explicitly. Parameterized loads are
  expanded to concrete arrays and written as valid non-parametric `NPETS=0`.
- TIB: load/write now preserves the raw body instead of splitting with regexes,
  which is safer for `U1DINT` node-list continuation lines.
- BCT: `ICBUND` loads as integer data, preventing invalid float constants on
  round-trip writes.
- BAS/WEL: added authoring-focused coverage for BAS options and WEL rates on
  GWF and CLN nodes, including `ITMP NP ITMPCLN` headers and AUX
  concentrations. CLN connectivity remains owned by the CLN package.
- CLN: added semantic support and authoring tests for `PROCESSCCF`/`ICLNGWCB`,
  `GENERAL_SEC`, and the 9-field `ISHAPE` node-property records used by
  rectangular/general conduit shapes.
- DPF: added authoring and load/write coverage for `FRAHK`, `IUZONTABIM`,
  conditional `SC2IM` by layer, immobile Richards arrays, and the programmatic
  `model.idpf` flag.
- Autotests added/updated for NAM external paths, 0-based boundary nodes,
  non-parametric HFB static/transient/structured layouts, TIB raw preservation,
  TVM semantic expectations, boundary package programmatic creation, CLN/DPF
  authoring, and ETS explicit parameter-write failure.

## Remaining Implementation Plan

1. Complete the Fortran audit for these files:
   `gwf2drt8u.f`, `glo2sgbu1.f`, and `gwf2QRT8u.f`.
2. Implement TABRICH BCF/LPF when `IUZONTAB` and `RETCRVS` are present.
3. Implement DRT USG-T transport extensions.
4. Add minimal semantic `MfUsgSgb` and `MfUsgQrt` packages.
5. Add slow/optional end-to-end validation against real USG-T models after the
   focused synthetic suite is green.

## Autotest Plan

- Install `modflow_devtools` if it is missing, then run targeted `pytest` and the
  broader USG-T suite.
- Add small synthetic fixtures per package; do not depend on large real models
  for default CI.
- Every semantic package needs at least one from-scratch authoring test that
  instantiates the package directly, writes it, and verifies 0-based internal /
  1-based file behavior without using `load()` as setup.
- Required focused tests:
  - CHD/RIV/GHB/DRN: API 0-based, file 1-based, AUX, reuse `-1`.
  - HFB: static, transient, `IHFBRD=-1/0/>0`, structured and DISU where useful.
  - ETS: `NETSEG=1`, `NETSEG>1`, `NETSOP=2`, `IESFACTOR`, parameterized expanded
    write behavior.
  - TIB: `U1DINT` lists with multiple nodes per continuation line.
  - NAM/external paths: input subdirectories preserved; outputs rebased.
  - BAS, WEL, CLN, and DPF regression coverage; DRT, SGB, and QRT according to
    remaining gaps.
  - `MfusgTransportListBudget`: old/new format, multi-species, species
    isolation.

## Assumptions

- Priority is upstream-readiness before broad new feature surface.
- The official Fortran source above replaces `src_patched` as the reference.
- Existing large-model regression runs remain valuable, but should be treated as
  slow/optional until the small synthetic suite is robust.
