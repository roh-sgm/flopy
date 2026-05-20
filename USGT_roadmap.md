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

---

## Package coverage table

All packages are listed in USG-T CUNIT array order from `mfusg.f`.

| Package | CUNIT | FloPy class | Status | Fortran file | Notes |
|---------|-------|-------------|--------|--------------|-------|
| BAS6 | `BAS6` | `MfUsgBas` | ⚠️ Partial | `gwf2basu1.f` | All standard options implemented. See Gap §1 for write-only omissions. **Verified** |
| DIS | `DIS` | `MfUsgDis` | ✅ | `mfusg.f` | Structured grid discretization |
| DISU | `DISU` | `MfUsgDisU` | ✅ | `mfusg.f` | Unstructured. `free_format_npl=10` default prevents buffer overflow for large grids. **Verified** |
| BCF6 | `BCF6` | `MfUsgBcf` | ⚠️ Partial | `gwf2bcf-lpf-u1.f` | TABRICH items 1c/1d (zone map + retention curves) not loaded/written. See Gap §2. **Verified** |
| LPF | `LPF` | `MfUsgLpf` | ⚠️ Partial | `gwf2bcf-lpf-u1.f` | Same TABRICH gap as BCF. Richards (LAYTYP=5) arrays correctly load/write after `Util2d.__eq__` fix (2026-05-20). **Verified** |
| SMS | `SMS` | `MfUsgSms` | ✅ Full | `glo2sms-u1.f` | All solver options. **Verified** |
| OC | `OC` | `MfUsgOc` | ✅ | — | ATS adaptive time-stepping and BOOTSTRAPPING supported |
| CHD | `CHD` | `MfUsgChd` | ✅ Full | `gwf2chd7u1.f` | Node-based unstructured. AUX transport concentrations. float64 precision. **Verified** |
| WEL | `WEL` | `MfUsgWel` | ⚠️ Partial | `gwf2wel7u1.f` | See Gap §3 (ITMPCLN) |
| DRN | `DRN` | `MfUsgDrn` | ✅ Full | `gwf2drn7u1.f` | Node-based unstructured. AUX transport concentrations. **Verified** |
| RIV | `RIV` | `MfUsgRiv` | ✅ Full | `gwf2riv7u1.f` | Node-based, AUX concentrations, optional `irch` column, float64 precision. **Verified** |
| GHB | `GHB` | `MfUsgGhb` | ✅ Full | `gwf2ghb7u1.f` | Node-based, AUX transport concentrations. **Verified** |
| RCH | `RCH` | `MfUsgRch` | ✅ Full | `gwf2rch8u1.f` | INRCHZONES crash bug fixed (2026-05-20). INCONC/INIRCH spacing fixed. **Verified** |
| EVT | `EVT` | `MfUsgEvt` | ✅ | `gwf2evt8u1.f` | **Verified** |
| ETS | `ETS` | `MfUsgEts` | ✅ Full | `gwf2ets8u1.f` | NETSEG>1 segment arrays, IESFACTOR transport flag. **Verified** |
| HFB | `HFB6` | `MfUsgHfb` | ✅ Full | `gwf2hfb7u1.f` | Node-based. `TRANSIENT_HFB` keyword and per-SP IHFBRD. **Verified** |
| GNC | `GNC` | `MfUsgGnc` | ✅ | `disu2gncn1.f` | **Verified** |
| LAK | `LAK` | `MfUsgLak` | ✅ | `gwf2lak7u1.f` | TABLEINPUT and TRANSPORTBOUNDARY options; lake transport coupling. **Verified (header)** |
| CLN | `CLN` | `MfUsgCln` | ⚠️ Partial | `cln2basu1.f` | See Gap §4 (ISHAPE, PROCESSCCF, NGENSHPTYP). **Verified** |
| DDF | `DDF` | `MfUsgDdf` | ✅ Full | `density.f` | RHOFRESH/RHOSTD/CSTD/ITHICKAV/IMPHDD/ISHARP + NONLINEAR table. ISHARP added 2026-05-20. **Verified** |
| BCT | `BCT` | `MfUsgBct` | ✅ Full | `glo2btnu1.f` | IDISP=1 and IDISP=2 (DLX/DLY/DLZ/DTXY/DTYZ/DTXZ), all 19 item 1a fields, A-W_ADSORB, ICHAIN, ISPRCT, ISOLUBILITY, IMULTI, IMASSWR options. **Verified** |
| PCB | `PCB` | `MfUsgPcb` | ✅ | — | **Verified (field order)** |
| MDT | `MDT` | `MfUsgMdt` | ✅ | — | Not independently verified against Fortran source |
| DPF | `DPF` | `MfUsgDpf` | ⚠️ Partial | `gwf2dpf1u1.f` | See Gap §5 (SC2IM conditional, Richards arrays). f_obj bug fixed 2026-05-20. **Verified** |
| DPT | `DPT` | `MfUsgDpt` | ⚠️ Partial | `gwt2dptu1.f` | DLIM condition fixed 2026-05-20 (now checks both IDPF and IDISPIM). See Gap §6 for remaining items. **Verified** |
| TIB | `TIB` | `MfUsgTib` | ✅ | — | Text round-tripper; sufficient for all load/write of existing files |
| TVM | `TVM` | `MfUsgTvm` | ✅ Full | `tvmu2.f` | Full semantic implementation (2026-05-20). HK/VKA/SS/SY/DDFTR/POR; transport-aware field counts; nper+1 boundaries; 20 autotests pass. **Verified** |
| GSF | `GSF` | `MfUsgGsf` | ✅ | — | Text round-trip. `to_grid()` delegates to `UnstructuredGrid.from_gridspec()` |
| SFR | `SFR` | base `ModflowSfr2` | ⚠️ Partial | `gwf2sfr7u1.f` | Base class used. Unstructured node-indexed items not independently validated in USG-T context |
| STR | `STR` | base `ModflowStr` | ⚠️ Partial | `gwf2str7u1.f` | Base class used. Unstructured format not validated |
| GAG | `GAGE` | base `ModflowGage` | ⚠️ Partial | `gwf2gag7u1.f` | Base class used. Depends on SFR/LAK; not independently validated |
| FHB | `FHB` | base `ModflowFhb` | ⚠️ Partial | `gwf2fhb7u1.f` | Base class used. Unstructured format not validated |
| DRT | `DRT` | base `ModflowDrt` | ⚠️ Partial | `gwf2drt8u.f` | Base class used. DRT8 USG-T extensions (transport AUX, IQCHANGEC, MXSPREADNDS) not in base class. **Verified (gap documented)** |
| SUB | `SUB` | base `ModflowSub` | ⚠️ Partial | `gwf2sub7u1.f` | Base class used. Unstructured items not validated |
| SWT | `SWT` | base `ModflowSwt` | ⚠️ Partial | — | Base class used. Unstructured format not validated |
| **SGB** | — | ❌ **Missing** | ❌ | `glo2sgbu1.f` | Specified Gradient Boundary. Not in FloPy registry. If present in a NAM file, `MfUsg.load()` will silently skip it. **Verified (absent)** |
| **QRT** | — | ❌ **Missing** | ❌ | `gwf2QRT8u.f` | Sink with Return Flow. Not in FloPy registry. Transport return-flow concentrations not supported. **Verified (absent)** |

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

### Gap §1 — BAS6: options parsed on load but not written

The following BAS6 options are correctly detected and stored on `load`, but
`write_file` never emits them. Models loaded with any of these flags will lose
them silently on round-trip:

`DPIN`, `DPOUT`, `DPIO` (double-precision binary output),
`PRINTTIME`, `SHOWPROGRESS`, `SY-ALL`

`RICHARDS` is separately handled: it IS written by `MfUsgBas` (confirmed). The
LPF Richards bug (LAYTYP=5 condition always False due to `Util2d.__eq__`) was
fixed (2026-05-20) — see `USGT_improvements.md`.

### Gap §2 — BCF/LPF: TABRICH items 1c and 1d not implemented

When `tabrich=True`, the Fortran reads two additional datasets after the standard
array block:

- **Item 1c**: `IUZONTAB` — integer zone map array (one zone index per node)
- **Item 1d**: `RETCRVS(NUZONES, NUTABROWS, 3)` — tabular moisture-retention curves
  (one row per zone, NUTABROWS points, columns: water content / pressure head / conductivity)

FloPy reads and stores `TABRICH`, `NUZONES`, and `NUTABROWS`, but does not load
or write items 1c/1d. Writing `tabrich=True` produces a syntactically incomplete
BCF/LPF file that USG-T cannot run.

**Source**: `gwf2bcf-lpf-u1.f` lines that read `IUZONTAB` and call
`U2DREL(RETCRVS, ...)` in the TABRICH block.

### Gap §3 — WEL: ITMPCLN missing from per-SP header

USG-T WEL writes a 3-token per-SP header `MXACT ITMPCWL ITMPCLN`. FloPy emits
only `MXACT ITMPCWL`. USG-T defaults `ITMPCLN=0` when absent, so this is
harmless if there are no CLN injection wells, but the written file does not
exactly match Vistas output.

### Gap §4 — CLN: three missing features (verified against `cln2basu1.f`)

1. **ISHAPE in 9-field node records** (`NRECTYP>0` or `NGENSHPTYP>0`): When
   conduit types use the new format, the Fortran reads a 9-field record
   `(IFNO, ISHAPE, IFTYPE, IFDIR, FLENG, FELEV, FANGLE, IFLIN, ICCWADI)`.
   The FloPy `cln_dtypes.py` dtype has 8 fields — ISHAPE is absent. Any CLN file
   with `nrectyp>0` will fail to round-trip: FloPy writes 8 fields, Fortran reads
   9, causing a field shift for all subsequent node records.

2. **PROCESSCCF option**: When present on the OPTIONS2 line, this flag activates a
   separate CLN-GWF exchange flow budget file (`ICLNGWCB` unit number). FloPy CLN
   has no parameter for it and neither reads nor writes it.

3. **NGENSHPTYP**: The general conduit shape type (tabular depth/area/wetted-perimeter
   input) is not implemented. The Fortran reads `NGENSHPTYP` from OPTIONS2 and then
   reads a table `(FDEPTH, FAREA, FWETPERI, FTOPWID)` per shape type. FloPy CLN
   has no parameter or parsing for this.

### Gap §5 — DPF: SC2IM conditionality and Richards arrays (verified against `gwf2dpf1u1.f`)

1. **SC2IM written unconditionally**: Fortran reads SC2IM only when `LAYCON(K)≠0`
   (line 174). FloPy writes it for every layer. For a DPF model with confined layers
   (`LAYCON=0`), the written file contains an extra array that Fortran does not
   consume, causing a read offset and corrupting all arrays that follow. This bug
   does not affect models where all layers are convertible (LAYCON=1,3,4), which
   is the common case.

2. **Richards arrays absent**: When `LAYCON(K)=5`, Fortran reads `alphaIM`,
   `betaIM`, `srIM`, `brookIM`, and optionally `bPIM` (lines 180–187). These are
   declared in the Fortran module but commented out in FloPy. A DPF model with
   Richards-equation layers cannot produce a valid DPF file with this implementation.

3. **IUZONTABIM**: Tabular moisture retention for the immobile domain is also
   commented out. Affects only `ITABRICH≠0` DPF models.

### Gap §6 — DPT: air-water interface adsorption immobile domain

`dpt2aw_adsorb.f` implements air-water interface adsorption for the immobile
domain. The physics is activated via a BCT option (`A-W_ADSORB`), but DPT may
need complementary parameters that have not been validated. Low priority —
`A-W_ADSORB` models are very uncommon.

### Gap §7 — DRT: USG-T 2.7 transport extensions (verified against `gwf2drt8u.f`)

The base `ModflowDrt` class is used. It does not support:

- Transport concentration auxiliary variables for inflowing return-flow water
- `IQCHANGEC` / `IQCHNGTYP` flow-reduction behavior flags
- `MXSPREADNDS` (spreading return flows to multiple nodes)

These are DRT8 USG-T additions not present in MODFLOW-2005 DRT. A model using
DRT with transport will load correctly (the static barriers are format-compatible)
but cannot specify return-flow concentrations.

---

## Missing packages

### SGB — Specified Gradient Boundary (`glo2sgbu1.f`)

SGB is listed in the USG-T 2.7 CUNIT array. It applies a specified hydraulic
gradient at the model boundary (distinct from GHB, which specifies head-to-head
conductance). It has its own I/O format and Fortran subroutines. There is no
FloPy class for it. If a USG-T model uses SGB, `MfUsg.load()` will silently
skip the package and the written model will be incomplete.

**Workaround**: SGB can often be approximated by GHB or CHD.

### QRT — Sink with Return Flow (`gwf2QRT8u.f`)

QRT8 is a flow package (analogous to DRT) that allows extracted water to be
returned to other nodes, with optional transport concentration for the returned
water, `IQCHANGEC`/`IQCHNGTYP` flow-reduction control, and `MXSPREADNDS`.
There is no FloPy class. If a model uses QRT, `MfUsg.load()` silently skips it.

**Workaround**: None available through FloPy. QRT-using models require manual
file management.

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
