# USG-Transport FloPy — Package Coverage Reference

This document is the authoritative record of which USG-Transport 2.7 packages
and features are supported in this fork, which are partial, and which are absent.

All entries marked **verified** were cross-referenced against the actual
USG-Transport 2.7.0 Fortran source (`USGT_V_2-7-0_Source_Code/`) in May 2026,
in addition to the I/O manual. Entries without that mark rely on the I/O manual
and/or round-trip testing only.

---

## Status legend

| Symbol | Meaning |
|--------|---------|
| ✅ Full | Implemented, verified against Fortran source, and round-trip tested |
| ✅ | Implemented and confirmed working, not exhaustively verified against source |
| ⚠️ Partial | Implemented but with known gaps (see Known Gaps section) |
| ❌ | Not implemented |

For USG-T packages, `Full` requires both semantic load/write behavior and
from-scratch programmatic authoring tests. A package that only preserves existing
text is classified as raw/text round-trip, even if it can write a runnable file.

---

## Package coverage table

All packages are listed in USG-T CUNIT array order from `mfusg.f`.

| Package | CUNIT | FloPy class | Status | Fortran file | Notes |
|---------|-------|-------------|--------|--------------|-------|
| BAS6 | `BAS6` | `MfUsgBas` | ✅ Full | `gwf2basu1.f` | All USG-T options load/write tested: UNSTRUCTURED, PRINTTIME, SHOWPROGRESS, DPIN/DPOUT/DPIO, SY-ALL, CONVERGE, plus `RICHARDS_HP` (implies Richards mode) and `IHM [IUIHM]` (Gap §1 resolved). **Verified** |
| DIS | `DIS` | `MfUsgDis` | ✅ | `mfusg.f` | Structured grid discretization |
| DISU | `DISU` | `MfUsgDisU` | ✅ | `mfusg.f` | Unstructured. `free_format_npl=10` default prevents buffer overflow for large grids. **Verified** |
| BCF6 | `BCF6` | `MfUsgBcf` | ✅ Full | `gwf2bcf-lpf-u1.f` | TABRICH items 1c (`IUZONTAB`) and 1d (`RETCRVS`, shape `(nuzones, nutabrows, 3)` = caphead/saturation/relperm) now authored/loaded/written; incomplete TABRICH write fails explicitly. See Gap §2 (resolved). **Verified** |
| LPF | `LPF` | `MfUsgLpf` | ✅ Full | `gwf2bcf-lpf-u1.f` | TABRICH 1c/1d implemented (RETCRVS replaces per-layer alpha/beta/sr/brook, which are skipped under TABRICH); parse bug fixed (nutabrows token, int cast). Richards (LAYTYP=5) arrays load/write after `Util2d.__eq__` fix. See Gap §2 (resolved). **Verified** |
| SMS | `SMS` | `MfUsgSms` | ✅ Full | `glo2sms-u1.f` | All solver options. **Verified** |
| OC | `OC` | `MfUsgOc` | ✅ | — | ATS adaptive time-stepping and BOOTSTRAPPING supported |
| CHD | `CHD` | `MfUsgChd` | ✅ Full | `gwf2chd7u1.f` | Node-based unstructured. Internal nodes 0-based, file I/O 1-based. AUX transport concentrations. Programmatic authoring tested. float64 precision. **Verified** |
| WEL | `WEL` | `MfUsgWel` | ✅ Full | `gwf2wel7u1.f` | Rates on GWF and CLN nodes, AUX transport concentrations, `ITMPCLN`, and programmatic authoring tested. CLN connectivity is handled by `MfUsgCln`. **Verified** |
| DRN | `DRN` | `MfUsgDrn` | ✅ Full | `gwf2drn7u1.f` | Node-based unstructured. Internal nodes 0-based, file I/O 1-based. AUX transport concentrations. Programmatic authoring tested. **Verified** |
| RIV | `RIV` | `MfUsgRiv` | ✅ Full | `gwf2riv7u1.f` | Node-based, internal nodes 0-based, file I/O 1-based, AUX concentrations, optional `irch` column, programmatic authoring tested, float64 precision. **Verified** |
| GHB | `GHB` | `MfUsgGhb` | ✅ Full | `gwf2ghb7u1.f` | Node-based, internal nodes 0-based, file I/O 1-based, AUX transport concentrations. Programmatic authoring tested. **Verified** |
| RCH | `RCH` | `MfUsgRch` | ✅ Full | `gwf2rch8u1.f` | INRCHZONES crash bug fixed (2026-05-20). INCONC/INIRCH spacing fixed. **Verified** |
| EVT | `EVT` | `MfUsgEvt` | ✅ | `gwf2evt8u1.f` | **Verified** |
| ETS | `ETS` | `MfUsgEts` | ⚠️ Partial | `gwf2ets8u1.f` | Authoring tested for NETSEG=1, NETSEG>1 (PXDP/PETM), NETSOP=2 (IEVT), and IESFACTOR transport flag. Parameterized files load as expanded non-parametric arrays (`NPETS=0` on write); programmatic `npets>0` fails explicitly. Parameter-syntax preservation intentionally not supported (Expanded valid write). **Verified** |
| HFB | `HFB6` | `MfUsgHfb` | ⚠️ Partial | `gwf2hfb7u1.f` | Node-based. Non-parametric static, structured static, and `TRANSIENT_HFB` IHFBRD = >0/0/-1 semantics implemented and tested. `NPHFB>0` (named parameters) fails explicitly on load/write. **Verified** |
| GNC | `GNC` | `MfUsgGnc` | ✅ | `disu2gncn1.f` | **Verified** |
| LAK | `LAK` | `MfUsgLak` | ✅ | `gwf2lak7u1.f` | TABLEINPUT and TRANSPORTBOUNDARY options; lake transport coupling. Load+write validated via the `Ex8_Lake` round-trip; from-scratch authoring of those options deferred (not `Full`). **Verified (header + Ex8 round-trip)** |
| CLN | `CLN` | `MfUsgCln` | ✅ | `cln2basu1.f`, `cln2props1.f` | `PROCESSCCF`/`ICLNGWCB`, `ISHAPE` node records, and `GENERAL_SEC` tabular shape authoring/load/write tested. **Verified** |
| DDF | `DDF` | `MfUsgDdf` | ✅ Full | `density.f` | RHOFRESH/RHOSTD/CSTD/ITHICKAV/IMPHDD/ISHARP + NONLINEAR table. From-scratch NONLINEAR-table authoring test added (2026-05-30). **Verified** |
| BCT | `BCT` | `MfUsgBct` | ✅ Full | `glo2btnu1.f` | IDISP=1 and IDISP=2 (DLX/DLY/DLZ/DTXY/DTYZ/DTXZ), all 19 item 1a fields, A-W_ADSORB, ICHAIN, ISPRCT, ISOLUBILITY, IMULTI, IMASSWR options. From-scratch authoring tests (1-species, IDISP=2, multi-species) added 2026-05-30. **Verified** |
| PCB | `PCB` | `MfUsgPcb` | ✅ | — | **Verified (field order)** |
| MDT | `MDT` | `MfUsgMdt` | ✅ | — | Not independently verified against Fortran source |
| DPF | `DPF` | `MfUsgDpf` | ✅ | `gwf2dpf1u1.f` | `FRAHK`, `IUZONTABIM`, conditional `SC2IM`, immobile Richards arrays, and programmatic `model.idpf` covered by focused tests. f_obj bug fixed 2026-05-20. **Verified** |
| DPT | `DPT` | `MfUsgDpt` | ⚠️ Partial | `gwt2dptu1.f` | DLIM condition checks both IDPF and IDISPIM. Immobile-domain air-water adsorption (`A-W_ADSORBIM`) now fails explicitly on load instead of silently shifting reads (Gap §6). **Verified** |
| TIB | `TIB` | `MfUsgTib` | ✅ Raw/text round-trip | `glo2basu1.f` | Raw-body text round-tripper; avoids parsing `U1DINT` node-list continuation lines. Sufficient for load/write of existing files. **No semantic constructor** — not authoring-ready by design (see Stage 3 audit / Card 4). |
| TVM | `TVM` | `MfUsgTvm` | ✅ Full | `tvmu2.f` | Full semantic implementation (2026-05-20). HK/VKA/SS/SY/DDFTR/POR; transport-aware field counts; nper+1 boundaries; 20 autotests pass. **Verified** |
| GSF | `GSF` | `MfUsgGsf` | ✅ Raw/text round-trip | — | Text round-trip (stores raw lines). `to_grid()` delegates to `UnstructuredGrid.from_gridspec()`. No semantic editing API by design. |
| SFR | `SFR` | base `ModflowSfr2` | ⚠️ Partial | `gwf2sfr7u1.f` | Base class used. Unstructured node-indexed items not independently validated in USG-T context |
| STR | `STR` | base `ModflowStr` | ⚠️ Partial | `gwf2str7u1.f` | Base class used. Unstructured format not validated |
| GAG | `GAGE` | base `ModflowGage` | ⚠️ Partial | `gwf2gag7u1.f` | Base class used. Depends on SFR/LAK; not independently validated |
| FHB | `FHB` | base `ModflowFhb` | ⚠️ Partial | `gwf2fhb7u1.f` | Base class used. Unstructured format not validated |
| DRT | `DRT` | `MfUsgDrt` | ✅ Full (authoring) / Expanded valid write (list controls) | `gwf2drt8u.f` | Node-based DRT8. EL+COND, RETURNFLOW (inline single recipient `NR>0` or `SPREAD` multi-node `NR<0` via U1DINT), `CHANGEC`/`IDCHNGTYP` transport, AUX, `ITMP/-1` reuse; 0-based internal / 1-based file. Reads `SFAC`/`OPEN-CLOSE`/`EXTERNAL` list controls (expanded inline on write). `recipient_nodes` validated against the stress list. `NPDRT>0` fails explicitly. **Unstructured only** — structured construction raises (use `ModflowDrt`); structured load delegates to base. **Verified (Phase 2)** |
| SUB | `SUB` | base `ModflowSub` | ⚠️ Partial | `gwf2sub7u1.f` | Base class used. Unstructured items not validated |
| SWT | `SWT` | base `ModflowSwt` | ⚠️ Partial | — | Base class used. Unstructured format not validated |
| SGB | `SGB` | `MfUsgSgb` | ✅ Full (authoring) / Expanded valid write (list controls) | `glo2sgbu1.f` | Specified Gradient Boundary. Node-based `(node, gradient)` list, AUX transport concentrations, `ITMP/-1` reuse; 0-based internal / 1-based file. Reads `SFAC` (inert on gradient per Fortran ISCLOC), `OPEN-CLOSE`, `EXTERNAL` list controls (expanded inline on write). `NPSGB>0` fails explicitly. Registered in `MfUsg.load()`. **Verified (Phase 2)** |
| QRT | `QRT` | `MfUsgQrt` | ✅ Full (authoring) / Expanded valid write (list controls) | `gwf2QRT8u.f` | Sink with Return Flow. Node-based `(node, q, rfprop)` + variable-length recipient-node lists (`NodQRT` via U1DINT), `CHANGEC`/`IQCHNGTYP` transport, AUX, `ITMP/-1` reuse. Reads `SFAC` (scales Q)/`OPEN-CLOSE`/`EXTERNAL` (expanded inline on write); `recipient_nodes` validated against the stress list. `AUTOFLOWREDUCE` preserved. `NPQRT>0` and `TRANSIENTQ` fail explicitly. Registered in `MfUsg.load()`. **Verified (Phase 2)** |

---

## Stage 3 status audit (2026-05-31)

Honest classification of every package that is not already a tested
`✅ Full`. "Decision" is the Stage 3 disposition; "Target" is the intended end
state. Cards refer to `USGT_STAGE3_COMPLETION_PLAN.md`.

| Package | Current | Decision | Target | Card |
|---------|---------|----------|--------|------|
| DIS | ✅ | finish now | Full semantic (foundational; exercised by all tests) | 8 |
| DISU | ✅ | finish now | Full semantic (protect large-grid `free_format_npl` formatting) | 8 |
| OC | ✅ | finish now | Full semantic (ATS / BOOTSTRAPPING authoring tests) | 8 |
| EVT | ✅ | finish now | Full semantic (transport-concentration array authoring) | 8 |
| CLN | ✅ | finish now | Full semantic (circular/rect/general-sec/PROCESSCCF tested) | 8 |
| DPF | ✅ | finish now | Full semantic (FRAHK/IUZONTABIM/SC2IM/immobile Richards tested) | 8 |
| PCB | ✅ | finish now | Full semantic (verify dataset order vs Fortran) | 8 |
| MDT | ✅ (unverified) | verify or demote | Full semantic *or* honest demotion | 8 |
| GNC | ✅ | keep | ✅ (GNC unstructured helper; verified) | — |
| LAK | ✅ | explicitly defer | `✅ not Full` — Ex8-validated; from-scratch authoring deferred | 5 |
| TIB | ✅ Raw/text | explicitly defer | Raw/text round-trip (no semantic constructor by design) | 4 |
| GSF | ✅ Raw/text | explicitly defer | Raw/text round-trip (+ `to_grid()`) | — |
| ETS | ⚠️ Partial | keep | Expanded valid write (param preservation deferred) | 2 |
| HFB | ⚠️ Partial | keep | Partial (non-param Full; `NPHFB>0` deferred) | 2 |
| DPT | ⚠️ Partial | keep | Partial (`A-W_ADSORBIM` explicit-unsupported) | 3 |
| SFR/STR/GAGE/FHB/SUB/SWT | ⚠️ Partial | compatibility-only | Compatibility-only (base FloPy classes; CLN preferred) | 6 |

Rule applied: no package is labeled `Full semantic` without a Fortran-derived
spec, from-scratch authoring tests, round-trip tests, and explicit failure for
unsupported modes. Plain `✅` rows above are honest "loads/writes, not yet
promoted" markers until their Card lands.

---

## Not applicable — solver and internal files

These Fortran files are solver internals, Fortran MODULE definitions, or utility
routines. They do not correspond to user-facing input packages and require no
FloPy class.

| Fortran file | Purpose |
|---|---|
| `glo2sms-u1.f`, `pcgu7.f`, `xmd.f`, `xmdlib_2.f`, `sparse.f` | Linear solvers (SMS, PCGU, XMD) |
| `mfusg.f` | Main driver, CUNIT registry, time-stepping loop |
| `utl7u1.f` | Array/token utility routines (U1DREL, URWORD, etc.) |
| `gmodules.f` | Global Fortran MODULE declarations |
| `gwf2bcf-lpf-u1.f` | BCF/LPF shared internals (TABRICH curves, LAYCON helpers) |
| `lak_gag_sfr_modules.f` | Fortran MODULE data structures for LAK/GAG/SFR |
| `cln2props1.f` | CLN property computation helpers |
| `disu2gncb1.f`, `disu2gncn1.f` | GNC unstructured helpers |
| `gwt2bndsu1.f` | Transport boundary budget kernel (not a package; called internally for WEL, GHB, DRN, DRT, QRT, EVT, ETS, RCH, RIV, LAK, PCB) |
| `glo2btnu1.f` | BCT core solver |
| `gwt2dptu1.f` | DPT solver |
| `density.f` | DDF solver |
| `tvmu1.f`, `tvmu2.f` | TVM interpolation engine |
| `dpt2aw_adsorb.f`, `gwt2aw_adsorb.f` | Air-water interface adsorption (exposed via BCT `A-W_ADSORB` option — implemented) |
| `parutl7.f` | PARAMETER keyword parser shared across packages |

---

## Known gaps within implemented packages

### Gap §1 — BAS6: niche options — RESOLVED

All previously-at-risk options are covered by programmatic write + load tests:
`UNSTRUCTURED`, `PRINTFV`, `CONVERGE`, `FREE`, `PRINTTIME`, `SHOWPROGRESS`,
`RICHARDS`, `DPIN`, `DPOUT`, `DPIO`, `SY-ALL`, and `STOPERROR`.

The two remaining niche options are now implemented (2026-05-30):

- `RICHARDS_HP`: Richards mode with pressure-head initial values. Implemented
  as `richards_hp=True`; it sets the effective Richards mode so BCF/LPF
  `LAYTYP=5` dependencies still hold. The BAS option-line cleaner now keeps
  `_` so `RICHARDS_HP` is not split into `RICHARDS HP`.
- `IHM [IUIHM]`: integrated-hydrologic-model coupling flag plus debug unit.
  Implemented as `ihm=True, iuihm=<unit>`.

The LPF Richards bug (LAYTYP=5 condition always False due to `Util2d.__eq__`) was
fixed (2026-05-20) — see `USGT_improvements.md`.

### Gap §2 — BCF/LPF: TABRICH items 1c and 1d — RESOLVED

**Resolved.** `MfUsgBcf` and `MfUsgLpf` now author/load/write items 1c
(`IUZONTAB`) and 1d (`RETCRVS`) via the shared `flopy/mfusg/_tabrich.py`
helper. `RETCRVS` is an ndarray of shape `(nuzones, nutabrows, 3)`
(capillary head / saturation / relative permeability), read zone-outer /
row-middle as the Fortran does. For LPF the per-layer Richards arrays
(alpha/beta/sr/brook) are skipped under TABRICH, matching `ITABRICH/=0`.
Writing `tabrich=True` without both arrays fails explicitly instead of
emitting an incomplete file.

Original gap (for reference): when `tabrich=True`, the Fortran reads two
additional datasets after the standard array block:

- **Item 1c**: `IUZONTAB` — integer zone map array (one zone index per node)
- **Item 1d**: `RETCRVS(NUZONES, NUTABROWS, 3)` — tabular moisture-retention curves
  (one row per zone, NUTABROWS points, columns: water content / pressure head / conductivity)

FloPy reads and stores `TABRICH`, `NUZONES`, and `NUTABROWS`, but does not load
or write items 1c/1d. Writing `tabrich=True` produces a syntactically incomplete
BCF/LPF file that USG-T cannot run.

**Source**: `gwf2bcf-lpf-u1.f` lines that read `IUZONTAB` and call
`U2DREL(RETCRVS, ...)` in the TABRICH block.

### Gap §3 — WEL: ITMPCLN resolved

USG-T WEL writes a 3-token per-SP header `ITMP NP ITMPCLN` when CLN is active.
`MfUsgWel` now has focused tests proving from-scratch authoring of stress-period
rates on GWF and CLN nodes, correct `ITMPCLN` headers, AUX concentration output,
and 0-based internal node storage after reload. Connectivity, geometry, and
CLN-GWF links remain the responsibility of `MfUsgCln`.

### Gap §4 — CLN: PROCESSCCF, ISHAPE, and GENERAL_SEC resolved

CLN now supports the USG-T 2.7 features verified against `cln2basu1.f` and
`cln2props1.f`:

- `PROCESSCCF` / `ICLNGWCB` on the `OPTIONS` line.
- 9-field node-property records with `ISHAPE` when `NRECTYP>0` or
  `NGENSHPTYP>0`.
- `GENERAL_SEC NGENSHPTYP NGENTABROWS` plus per-type
  `(FDEPTH, FAREA, FWETPERI, FTOPWID)` tables.

Focused tests cover from-scratch authoring, file output, and reload.

### Gap §5 — DPF: SC2IM, Richards arrays, and IUZONTABIM resolved

DPF now matches the key conditional reads in `gwf2dpf1u1.f`:

- `SC2IM` is written/loaded only for layers where `LAYCON(K) != 0`.
- For `LAYCON(K)=5` with `TABRICH` off, `alphaIM`, `betaIM`, `srIM`, `brookIM`,
  and optional `bPIM` are written/loaded.
- For `ITABRICH != 0`, the immobile-domain `IUZONTABIM` node map is
  written/loaded before `IBOUNDIM`.
- Programmatic construction sets `model.idpf = 1`, so dependent packages see
  DPF as active.

Focused tests cover from-scratch authoring, file output, and reload.

### Gap §6 — DPT: air-water interface adsorption immobile domain (Stage 3 Card 3)

**Decision: explicitly unsupported (deferred); rare PFAS-type sub-mode.**

Fortran-derived spec (`gwt2dptu1.f` option block + `dpt2aw_adsorb.f`
`AW_ADSORBIM1AL/RP1/RP2`): the DPT option line carries
`A-W_ADSORBIM IAREA_FNIM IKAWI_FNIM`, then a cascade of conditional reads:

- If `IAREA_FNIM==5` or `IKAWI_FNIM==4` (tabular): `NAZONESIM NATABROWSIM`, a
  per-node zone map `IAWIZONMAPIM` (`U1DINT`, NODES), and a tabular area array
  `AWI_AREA_TABIM(2, NATABROWSIM, NAZONESIM)`.
- `IAREA_FNIM==1/3`: `AWAMAXIM(NODES)`; `==4`: `AWAREA_X2IM/X1IM/X0IM(NODES)`;
  else a `ROG_SIGMAIM` constant.
- Langmuir isotherm arrays `ALANGAWIM(NODES,MCOMP)`, `BLANGAWIM(NODES,MCOMP)`,
  plus per-stress-period reads (`AW_ADSORBIM1RP1/RP2`).

This many conditional arrays make it a substantial, rarely-used sub-mode, so it
is **not** modeled in v1. To prevent silently shifting every subsequent read,
`MfUsgDpt.load` **fails explicitly** with `NotImplementedError` the moment the
`A-W_ADSORBIM` token is seen on the option line — before any of the extra reads
above — for both the bare keyword and the `IAREA_FNIM IKAWI_FNIM` form (tested).
The DLIM conditional read correctly requires both `IDPF/=0` and `IDISPIM/=0`,
and is independent of this option.

### Gap §7 — DRT: USG-T 2.7 transport extensions — RESOLVED

**Resolved.** `MfUsgDrt` (subclass of `ModflowDrt`) implements the node-based
DRT8 unstructured format with the USG-T extensions:

- `CHANGEC` / `IDCHNGTYP` per-cell return-flow concentration-change type.
- `RETURNFLOW` recipient nodes: a single inline recipient (`NR>0`) or a
  `SPREAD` multi-node spreading ground (`NR<0`, recipients via a `U1DINT`
  block written immediately after the drain line). `MXSPREADNDS` is
  recomputed from the data on write.
- AUX concentration variables and `ITMP/-1` stress-period reuse.

`NPDRT>0` (named parameters) fails explicitly. Structured (DIS) models
delegate to the base `ModflowDrt`. Return-flow node lists with
`EXTERNAL`/`OPEN/CLOSE` control records are not supported (explicit failure).

### Gap §8 — ETS: parameter syntax is expanded, not preserved

USG-T ETS supports named parameters (`NPETS > 0`) for the ETSR array. FloPy now
loads those definitions through the shared MODFLOW parameter reader and expands
them into concrete ETSR arrays. The returned `MfUsgEts` object writes a valid
non-parametric ETS file with `NPETS=0`.

This is intentionally classified as partial: preserving the original parameter
definitions and active-parameter records is not implemented. Programmatic
`MfUsgEts(..., npets > 0).write_file()` raises `NotImplementedError` rather than
writing a misleading or incomplete parametric package.

### Gap §9 — HFB: parameterized barriers are not supported

For non-parametric HFB, the implementation matches `gwf2hfb7u1.f`: static HFB
reads/writes `NHFBNP` rows once, and transient HFB reads `IHFBRD` for every stress
period, using `IHFBRD <= 0` as reuse/no-read and `IHFBRD > 0` as the flag to read
exactly `NHFBNP` rows.

Parameterized HFB (`NPHFB > 0`) is not implemented yet. Load and write fail
explicitly for that case until parameter expansion/preservation is added.

---

## Previously-missing packages — now implemented

### SGB — Specified Gradient Boundary (`glo2sgbu1.f`) — IMPLEMENTED

`MfUsgSgb` applies a specified hydraulic gradient at boundary nodes (distinct
from GHB, which specifies a head-to-head conductance). Node-based
`(node, gradient)` list with AUX transport concentrations and `ITMP/-1`
reuse; internal nodes 0-based, file 1-based. Registered as `"sgb"` in
`MfUsg.load()`, so SGB models are no longer silently skipped. `NPSGB>0`
(named parameters) fails explicitly.

### QRT — Sink with Return Flow (`gwf2QRT8u.f`) — IMPLEMENTED

`MfUsgQrt` extracts water at sink nodes and returns a proportion to one or
more recipient nodes (analogous to DRT). Per-sink `(node, q, rfprop)` with
variable-length recipient-node lists (`NodQRT`, read/written via `U1DINT`),
`CHANGEC`/`IQCHNGTYP` return-flow concentration-change type, AUX, and
`ITMP/-1` reuse. `AUTOFLOWREDUCE` is preserved as an option. `NPQRT>0` and
the `TRANSIENTQ` transient-flow time-series option fail explicitly.
Registered as `"qrt"` in `MfUsg.load()`.

---

## Packages intentionally not implemented

| Package | Reason |
|---------|--------|
| STR, SFR | CLN is the preferred surface-water coupling for USG-T models in this project |
| MNW1, MNW2 | Replaced by CLN in USG-T applications |
| OBS, HOB, DROB, RVOB, GBOB, CHOB | Observation packages — low priority |
| UZF | Unsaturated zone flow. Potentially relevant but out of scope |
| IHM / SYF | Runtime coupling to external HSPF. No standalone input package needed |
| GMG, PCG, SIP, DE4 | Solvers not used with USG-T; SMS and XMD are the relevant solvers |

---

## Validation record

| Model | Grid | Version | Packages | Result |
|-------|------|---------|----------|--------|
| MD26 BC2 (real-world, ~112k nodes, 493 SPs) | DISU | USG-T 1.8 | BAS6, DISU, LPF, WEL, GHB, DRN, CHD, RIV, RCH, EVT, CLN, BCT, TIB, OC, SMS | Bit-for-bit identical on all budget columns (flow + transport) |
| Sóncor (real-world, ~19k nodes, 1382 SPs, 4322 TSs) | DISU | USG-T 2.7 ARM | BAS6, DISU, LPF, WEL, GHB, DRN, CHD, RIV, RCH, EVT, CLN, BCT (IDISP=2), DDF, OC, SMS | Bit-for-bit identical on LST, HDS, CON, CBB (7 record types) |
