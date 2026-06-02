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
| DIS | `DIS` | `MfUsgDis` | ✅ Full | `mfusg.f` | Structured grid discretization. Authored in essentially every synthetic test and loaded in the structured Ex* models. **Verified** |
| DISU | `DISU` | `MfUsgDisU` | ✅ Full | `mfusg.f` | Unstructured grid. Authored/loaded across the suite; `free_format_npl=10` default prevents large-grid buffer overflow, with dedicated round-trip formatting tests in `autotest/test_usg.py` (`test_free_format_npl`). **Verified** |
| BCF6 | `BCF6` | `MfUsgBcf` | ✅ Full | `gwf2bcf-lpf-u1.f` | TABRICH items 1c (`IUZONTAB`) and 1d (`RETCRVS`, shape `(nuzones, nutabrows, 3)` = caphead/saturation/relperm) now authored/loaded/written; incomplete TABRICH write fails explicitly. See Gap §2 (resolved). **Verified** |
| LPF | `LPF` | `MfUsgLpf` | ✅ Full | `gwf2bcf-lpf-u1.f` | TABRICH 1c/1d implemented (RETCRVS replaces per-layer alpha/beta/sr/brook, which are skipped under TABRICH); parse bug fixed (nutabrows token, int cast). Richards (LAYTYP=5) arrays load/write after `Util2d.__eq__` fix. See Gap §2 (resolved). **Verified** |
| SMS | `SMS` | `MfUsgSms` | ✅ Full | `glo2sms-u1.f` | All solver options. **Verified** |
| OC | `OC` | `MfUsgOc` | ✅ | `glo2basu1.f` (SGWF2BAS7I/J/N) | OC Fullness Card A (Stage 4): from-scratch authoring + round-trip tested for `ATSA`/`NPTIMES`/`NPSTPS`, the `BOOTSTRAPPING` header (now written on the first OC line where USG-T parses it, and preserved on load), per-SP `BOOTSTRAP`/`NOBOOTSTRAP`/`BOOTSTRAPSCALE`/`NOBOOTSTRAPSCALE`, layer-qualified `PRINT`/`SAVE HEAD`/`DRAWDOWN`/`CONC`, `DDREFERENCE`, `SAVE`/`PRINT CONC`/`BUDGET`, `COMPACT BUDGET AUX`. Kept `✅` not `Full` — explicit gaps: `SAVE IBOUND` is commented out in USG-T 2.7 (solver rejects it; FloPy preserves the keyword and `check()` warns about it); `FASTFORWARD`/`FASTFORWARDC` separate-line placement not exe-verified; `BOOTSTRAPPING` execution not exe-verified; numeric-format OC rewrites as words. See `USGT_STAGE4_OC_FULLNESS.md`. |
| CHD | `CHD` | `MfUsgChd` | ✅ Full | `gwf2chd7u1.f` | Node-based unstructured. Internal nodes 0-based, file I/O 1-based. AUX transport concentrations. Programmatic authoring tested. float64 precision. **Verified** |
| WEL | `WEL` | `MfUsgWel` | ✅ Full | `gwf2wel7u1.f` | Rates on GWF and CLN nodes, AUX transport concentrations, `ITMPCLN`, and programmatic authoring tested. CLN connectivity is handled by `MfUsgCln`. **Verified** |
| DRN | `DRN` | `MfUsgDrn` | ✅ Full | `gwf2drn7u1.f` | Node-based unstructured. Internal nodes 0-based, file I/O 1-based. AUX transport concentrations. Programmatic authoring tested. **Verified** |
| RIV | `RIV` | `MfUsgRiv` | ✅ Full | `gwf2riv7u1.f` | Node-based, internal nodes 0-based, file I/O 1-based, AUX concentrations, optional `irch` column, programmatic authoring tested, float64 precision. **Verified** |
| GHB | `GHB` | `MfUsgGhb` | ✅ Full | `gwf2ghb7u1.f` | Node-based, internal nodes 0-based, file I/O 1-based, AUX transport concentrations. Programmatic authoring tested. **Verified** |
| RCH | `RCH` | `MfUsgRch` | ✅ Full | `gwf2rch8u1.f` | INRCHZONES crash bug fixed (2026-05-20). INCONC/INIRCH spacing fixed. **Verified** |
| EVT | `EVT` | `MfUsgEvt` | ✅ | `gwf2evt8u1.f` | EVT Fullness Card B (Stage 4): from-scratch authoring + round-trip tested for `NEVTOP`=1/2/3 (structured + unstructured `NEVTOP=2` with `MXNDEVT`; `IEVT` 0-based internal / 1-based file, layer-range (structured) and node-range (unstructured) validated), transport `IETFACTOR` 0/<0/>0 with per-`MCOMP` `ETFACTOR` (scalar or array), and per-SP reuse (`-1`). Fixed a real round-trip bug (`load` dropped `IETFACTOR`) and made the writer emit 3 dataset-1 ints whenever transport is active. Plus a USG-T 2.7 executable smoke (`test_usgt_exe_evt_from_scratch`). Kept `✅` not `Full` — explicit gaps: ETS zonal time-series (`ETS MXZNEVT`/`IZNEVT`) raises `NotImplementedError`; `NPEVT` parameters are expanded-valid-write (loaded as arrays, `NP=0` on write), authoring-with-params unsupported. See `USGT_STAGE4_EVT_FULLNESS.md`. **Verified** |
| ETS | `ETS` | `MfUsgEts` | ⚠️ Partial | `gwf2ets8u1.f` | Authoring tested for NETSEG=1, NETSEG>1 (PXDP/PETM), NETSOP=2 (IEVT), and IESFACTOR transport flag. **Stage 4.4A: ETSR array parameters are now preserved** (load → write → reload keeps `NPETS>0` in item 2a, the definition blocks, and per-period activation records, incl. `INSTANCES`; ETSS/ETSX/IETS/PXDP/PETM stay plain arrays — the mix the Fortran allows). USG-T reads `NPETS` from item 2a (`UPARARRAL` called with `IN=-1`), so no `PARAMETER` line is written. Opt-in expanded fallback (`expand_parameters=True`) still writes `NPETS=0`; from-scratch parameter *authoring* (`npets>0` without loaded defs) fails explicitly. Shared write helper `flopy/mfusg/_usgt_parameters.py` reuses `ModflowParBc` (no second parser). Not `Full` (no from-scratch param authoring). See `USGT_STAGE4_04_PARAMETERS_ETS.md`. **Verified** |
| HFB | `HFB6` | `MfUsgHfb` | ⚠️ Partial | `gwf2hfb7u1.f` | Node-based. Non-parametric static, structured static, and `TRANSIENT_HFB` IHFBRD = >0/0/-1 semantics implemented and tested. **Stage 4.4B: HFB list parameters (`NPHFB>0`) are now preserved** (load → write → reload keeps the definition blocks with `NLST` barrier rows, `NACTHFB`, and the active names; barriers 0-based internal / 1-based file; parameter-defined and non-parametric barriers can mix). Uses the shared list-parameter helpers in `flopy/mfusg/_usgt_parameters.py` (`UPARLSTRP`/`UPARLSTSUB`). Barrier lists (parameter `NLST` rows and non-parametric `NHFBNP` rows) also support leading `SFAC`/`OPEN/CLOSE`/`EXTERNAL` list controls via the shared `_usgt_list.begin_list_block` (list-control follow-up; `SFAC` scales `HYDCHR`, `EXTERNAL` via `ext_unit_dict` else `NotImplementedError`). Not `Full`: from-scratch parameter authoring, `TRANSIENT_HFB`+`NPHFB>0`, and parameter `INSTANCES` raise `NotImplementedError` (the last matches the Fortran). See `USGT_STAGE4_04_PARAMETERS_HFB.md`. **Verified** |
| GNC | `GNC` | `MfUsgGnc` | ✅ | `disu2gncn1.f` | **Verified** |
| LAK | `LAK` | `MfUsgLak` | ✅ | `gwf2lak7u1.f` | LAK Fullness Card D (Stage 4): from-scratch authoring + round-trip tested for no-transport, classic transport (`CPPT`/`CRNF`), `TRANSPORTBOUNDARY` (one `CLAKE(1:NSOL)` line per lake, MCOMP>1), and `TABLEINPUT` (per-lake tab unit + external registration). Fixed six authoring bugs (conc_data mis-assignment; `flux_data=None` crash; conc_data accessed without transport; TRANSPORTBOUNDARY 9b written per-component instead of one-line-per-lake; load stored conc as strings; `transportboundary` flag not synced to the header keyword) and added validation (flux_data required; TRANSPORTBOUNDARY needs transport; clake nlakes×mcomp; transport needs conc_data). A review follow-up further validates authoring inputs in `__init__` so malformed input fails with a clear `ValueError` instead of crashing in `write_file()`: TABLEINPUT requires exactly one `tab_file`/`tab_unit` per lake; `conc_data` is checked per written period (every `(lake, component)` present; classic 2 values for `WTHDRW>=0` / 3 for `WTHDRW<0`; TRANSPORTBOUNDARY a single value); `flux_data` requires one dataset-9a entry per lake. `Ex8_Lake` real-model round-trip/run kept. Kept `✅` not `Full` — gaps: sill/connectivity (ds 7/8) and multi-lake systems round-trip but aren't authored from scratch in tests; `TABLEINPUT` bathymetry table contents are external; GAGE coupling separate. See `USGT_STAGE4_LAK_FULLNESS.md`. |
| CLN | `CLN` | `MfUsgCln` | ✅ Full | `cln2basu1.f`, `cln2props1.f` | From-scratch authoring + round-trip tested for `PROCESSCCF`/`ICLNGWCB`, the 9-field `ISHAPE` node records, and `GENERAL_SEC` tabular shapes; also exercised by the Ex3 CLN conduit real models. Promoted to Full (Stage 3 Card 8). **Verified** |
| DDF | `DDF` | `MfUsgDdf` | ✅ Full | `density.f` | RHOFRESH/RHOSTD/CSTD/ITHICKAV/IMPHDD/ISHARP + NONLINEAR table. From-scratch NONLINEAR-table authoring test added (2026-05-30). **Verified** |
| BCT | `BCT` | `MfUsgBct` | ✅ Full | `glo2btnu1.f` | IDISP=1 and IDISP=2 (DLX/DLY/DLZ/DTXY/DTYZ/DTXZ), all 19 item 1a fields, A-W_ADSORB, ICHAIN, ISPRCT, ISOLUBILITY, IMULTI, IMASSWR options. From-scratch authoring tests (1-species, IDISP=2, multi-species) added 2026-05-30. **Verified** |
| PCB | `PCB` | `MfUsgPcb` | ✅ Full | `gwt2bndsu1.f` | Prescribed Concentration Boundary. Node-based `(node, iSpec, conc)`; 0-based internal / 1-based file. Field order verified vs Fortran; from-scratch authoring + round-trip tested (Stage 3 Card 8); exercised via the Ex transport models. **Verified** |
| MDT | `MDT` | `MfUsgMdt` | ✅ | `gwt2mdtu1.for` | Matrix Diffusion Transport. MDT Fullness Card C (Stage 4): field order Fortran-audited and from-scratch authoring + round-trip tested for the header options (`FRAHK`/`FRADARCY`/`TSHIFTMD`/`SEPARATE_AI2`/`MULTIFILE_MD`, only when `IDPF==0`), base arrays (`MDFLAG`/`VOLFRACMD`/`PORMD`/`RHOBMD`/`DIFFLENMD`/`TORTMD`; `VOLFRACMD` skipped when `IDPF!=0`), per-species `KDMD`/`DECAYMD`/`YIELDMD`/`DIFFMD`, and `AIOLD1MD`/`AIOLD2MD` under `TSHIFTMD > 1e-10` (the Fortran threshold, applied consistently; TSHIFTMD written in general format so small valid values are not rounded to 0). Fixed real bugs (`load` dropped `FRAHK`/`FRADARCY` via a case mismatch; `write_file(f=handle)` crashed; stray debug print) and added validation (FRAHK⊕FRADARCY, IDPF-options, MULTIFILE needs `imdtcf>0`, per-`MCOMP` list lengths). Still **round-trip-run via the three Ex7 real models**. Kept `✅` not `Full` — gaps: species loop uses `MCOMP` (chained-decay `NTCOMP>MCOMP` not independently verified); AI1/AI2 output binaries authored but not read. See `USGT_STAGE4_MDT_FULLNESS.md`. |
| DPF | `DPF` | `MfUsgDpf` | ✅ Full | `gwf2dpf1u1.f` | From-scratch authoring + round-trip tested: `FRAHK`, `IUZONTABIM`, conditional `SC2IM` by layer, immobile Richards arrays, and the programmatic `model.idpf` flag. f_obj bug fixed 2026-05-20. Promoted to Full (Stage 3 Card 8). **Verified** |
| DPT | `DPT` | `MfUsgDpt` | ⚠️ Partial | `gwt2dptu1.f` | DLIM condition checks both IDPF and IDISPIM. Immobile-domain air-water adsorption (`A-W_ADSORBIM`) now fails explicitly on load instead of silently shifting reads (Gap §6). **Verified** |
| TIB | `TIB` | `MfUsgTib` | ✅ Full (authoring) | `glo2basu1.f` | Transient Ibound. Semantic `stress_period_data` authoring + `parse=True` load covering the full `GWF2TIB1RP` grammar: flow blocks (`NIB0` inactivate / `NIB1` activate / `NIBM1` prescribed-head, with `HEAD`/`AVHEAD`/bare-reuse) and transport blocks (`NICB0`/`NICB1`/`NICBM1`, with multi-component `CONC`/`AVCONC`/bare); `NIB0`/`NICB0` node lists via `U1DINT`; 3- vs 6-int header keyed to BCT presence; 0-based internal / 1-based file. From-scratch authoring, semantic write/reload, multi-node `U1DINT`, transport-conc, and explicit-rejection tests. `MfUsgTib.load` defaults to byte-exact raw round-trip (what `MfUsg.load` uses); `load(parse=True)` returns semantic data and falls back to raw on `U1DINT` `EXTERNAL`/`OPEN-CLOSE` rather than writing partial data. **Verified (Stage 4.1)** |
| TVM | `TVM` | `MfUsgTvm` | ✅ Full | `tvmu2.f` | Full semantic implementation (2026-05-20). HK/VKA/SS/SY/DDFTR/POR; transport-aware field counts; nper+1 boundaries; 20 autotests pass. **Verified** |
| GSF | `GSF` | `MfUsgGsf` | ✅ Full (authoring) | — | Grid Specification File (not solver input — consumed by post-processors / `UnstructuredGrid.from_gridspec`, so no Fortran reader). Stage 4.2: semantic `vertices` + `node_data` authoring and `parse=True` load — header (`UNSTRUCTURED`/`UNSTRUCTURED GWF`, validated), `nnodes`/`nlay`, vertex `(x,y,z)`, per-node `(node, xc, yc, zc, layer, vertices)`; 0-based internal nodes/vertices/layers, 1-based file ids; trailing line-2 gridgen flags preserved (default `(1,1)`). `from_grid(model, grid, top_zverts, bot_zverts)` writes the USG-T top/bottom doubled-vertex convention (single-surface `zverts` kept as legacy; per-vertex z required — the grid does not retain it); `from_disv_gridprops(model, disv_gridprops, top, botm, skip_degenerate=)` maps a MF6 DISV 2D template to a single-layer GSF (closing-vertex removal, degenerate-cell skip). Both accept `vertex_mode="shared"`/`"parsimonious"` (default; neighbours reuse vertex ids) or `"cell"`/`"nonparsimonious"` (unique top/bottom vertices per cell, 8 per quad) — the two GRIDGEN2GSF layouts (`gridgen2gsf.f90`, secondary ref; primary spec gwutil_a 2.17, line 2 = `NNODES NLAY IZ IC`). `to_grid()` delegates to `from_gridspec` (which now correctly accepts `UNSTRUCTURED GWF`). `load` defaults to byte-faithful raw round-trip; `parse=True` is strict (header exactly `UNSTRUCTURED`/`UNSTRUCTURED GWF`; `IZ IC` omitted or `1 1`; `inode` = `1..nnode` in order; trailing content) → raw fallback otherwise. Constructor validates header / `IZ IC`=`(1,1)` / contiguous-ordered node ids `0..nnodes-1` / `nlay >= max(layer)+1`. From-scratch authoring, write/reload, 0-based↔1-based, `to_grid`, invalid-ref, mode-exclusivity, top/bottom, DISV, vertex-mode, strict-header, IZ-IC, and inode tests. Companion utility `gridgen_to_gsf` (`flopy/mfusg/gridgen2gsf.py`, Stage 4.3) builds a GSF from a flopy `Gridgen` / `disv_gridprops` / `UnstructuredGrid` source (a non-interactive helper inspired by `GRIDGEN2GSF`; parsimonious mode compacts unused vertices; delegates to `MfUsgGsf`). **Verified (Stage 4.2 + follow-ups; Stage 4.3 utility)** |
| SFR | `SFR` | base `ModflowSfr2` | Compatibility (DISU-validated) | `gwf2sfr7u1.f` | Base MODFLOW-2005 class, **kept enabled**: `examples/data/freyberg_usg` is a DISU model that uses SFR and is loaded → written → run via the base class in `autotest/test_usg.py` (`test_freyburg_usg`/`_external`), so unstructured SFR is empirically supported. **Not** promoted to Full USG-T transport (no transport-coupled authoring/tests); **not guarded** (Stage 4.6A deliberately leaves SFR on the base class). CLN remains the preferred coupling. |
| STR | `STR` | `MfUsgStr` (guards base `ModflowStr`) | Compatibility-only guarded | `gwf2str7u1.f` | **Stage 4.6A:** USG-T DISU layout unvalidated; using STR on an unstructured `MfUsg` model (load **or** authoring) raises `NotImplementedError` before any read/write rather than silently producing a non-USG-T file. Structured `MfUsg` and plain `flopy.modflow.Modflow` are unaffected (delegate to / use base `ModflowStr`). Not used by any target/example model. |
| GAG | `GAGE` | base `ModflowGage` | Compatibility-only | `gwf2gag7u1.f` | Base class. **Loads via base in the Ex8 real model** (round-trip exercised there); USG-T-specific records not independently validated. Compatibility-only (Card 6). |
| FHB | `FHB` | base `ModflowFhb` | Compatibility-only | `gwf2fhb7u1.f` | Base class. **Loads via base in the Ex8 real model** (`m.fhb` asserted; round-trip exercised there); USG-T unstructured format not independently validated. Compatibility-only (Card 6). |
| DRT | `DRT` | `MfUsgDrt` | ✅ Full (authoring) / Parameter-preserving (NPDRT, incl. activations + recipients) / Expanded valid write (list controls) | `gwf2drt8u.f` | Node-based DRT8. EL+COND, RETURNFLOW (inline single recipient `NR>0` or `SPREAD` multi-node `NR<0` via U1DINT), `CHANGEC`/`IDCHNGTYP` transport, AUX, `ITMP/-1` reuse; 0-based internal / 1-based file. Reads `SFAC`/`OPEN-CLOSE`/`EXTERNAL` list controls (expanded inline on write). `recipient_nodes` validated against the stress list. **Stage 4.5B: the spreading-recipient `U1DINT` block (`NR<0`) also reads `EXTERNAL <unit>` / `OPEN-CLOSE <fname>` (resolved via `ext_unit_dict`/`model_ws`, quote-aware, **free-format `(FREE)` only** — a fixed `FMTIN` raises `NotImplementedError`; expanded inline on write; unresolved `EXTERNAL` raises `NotImplementedError`).** **Stage 4.4D: `NPDRT>0` list parameters are now preserved** (load → write → reload keeps `NPDRT`/`MXL` in item 1, definitions with `NLST` drain rows + their RETURNFLOW recipients/spreading blocks, and per-SP `ITMP NP` + active names). DRT is internally consistent — definition and activation both use `PARTYP='DRT'` (`gwf2drt8u.f:135`/`:1108`), so active parameters are valid (unlike SGB); parameter value scales `COND` (`IPVL=5`). Not `Full` on the parametric axis: from-scratch authoring and `INSTANCES` (Fortran-supported) raise `NotImplementedError`; `MXL`/consistency validated (`ValueError`, no partial file). **Unstructured only** — structured construction raises (use `ModflowDrt`); structured load delegates to base. See `USGT_STAGE4_04_PARAMETERS_DRT.md`. **Verified** |
| SUB | `SUB` | `MfUsgSub` (guards base `ModflowSub`) | Compatibility-only guarded | `gwf2sub7u1.f` | **Stage 4.6A:** USG-T items unvalidated; using SUB on an unstructured `MfUsg` model (load or authoring) raises `NotImplementedError` before any read/write. Structured `MfUsg` and plain `flopy.modflow.Modflow` unaffected. Not used by any target/example model. |
| SWT | `SWT` | `MfUsgSwt` (guards base `ModflowSwt`) | Compatibility-only guarded | — | **Stage 4.6A:** USG-T format unvalidated; using SWT on an unstructured `MfUsg` model (load or authoring) raises `NotImplementedError` before any read/write. Structured `MfUsg` and plain `flopy.modflow.Modflow` unaffected. Not used by any target/example model. |
| SGB | `SGB` | `MfUsgSgb` | ✅ Full (authoring) / Parameter-definition-preserving (NPSGB, no activations) / Expanded valid write (list controls) | `glo2sgbu1.f` | Specified Gradient Boundary. Node-based `(node, gradient)` list, AUX transport concentrations, `ITMP/-1` reuse; 0-based internal / 1-based file. Reads `SFAC` (inert on gradient per Fortran ISCLOC), `OPEN-CLOSE`, `EXTERNAL` list controls (expanded inline on write). **Stage 4.4C (+ review follow-up): `NPSGB>0` parameter *definitions* are preserved** (load → write → reload keeps the `PARAMETER NPSGB MXS` record + definitions with `NLST` `NODE GRADIENT [aux]` rows; shared helpers in `_usgt_parameters.py`). **Active SGB parameters are unsupported**: USG-T 2.7 defines them as `PARTYP='SGB'` (`UPARLSTRP`, glo2sgbu1.f:97) but activates them as `PTYP='G'` (`UPARLSTSUB`, glo2sgbu1.f:185) — a type conflict (`parutl7.f:684/800`) that aborts the run — so per-SP `NP>0` and `INSTANCES` raise `NotImplementedError`. Writer validates `MXS`; inconsistent state raises `ValueError` (no partial file). Registered in `MfUsg.load()`. See `USGT_STAGE4_04_PARAMETERS_SGB.md`. **Verified** |
| QRT | `QRT` | `MfUsgQrt` | ✅ Full (authoring) / Parameter structural-preserving (NPQRT, not execution-guaranteed) / Expanded valid write (list controls) | `gwf2QRT8u.f` | Sink with Return Flow. Node-based `(node, q, rfprop)` + variable-length recipient-node lists (`NodQRT` via U1DINT), `CHANGEC`/`IQCHNGTYP` transport, AUX, `ITMP/-1` reuse. Reads `SFAC` (scales Q)/`OPEN-CLOSE`/`EXTERNAL` (expanded inline on write); `recipient_nodes` validated against the stress list. **Stage 4.5B: the recipient `U1DINT` block also reads `EXTERNAL <unit>` / `OPEN-CLOSE <fname>` (resolved via `ext_unit_dict`/`model_ws`, quote-aware, **free-format `(FREE)` only** — a fixed `FMTIN` raises `NotImplementedError`; expanded inline on write; unresolved `EXTERNAL` raises `NotImplementedError`).** `AUTOFLOWREDUCE` preserved. **Stage 4.4E: `NPQRT>0` list parameters are structurally preserved** (load → write → reload of definitions + recipient blocks + per-SP `ITMP NP` activations; `MXAQRT` covers the active total; `MXRTCELLS` reflects definition recipients). QRT is type-consistent (`PARTYP='QRT'` both sides, `gwf2QRT8u.f:183`/`:1118`), but **not execution-guaranteed**: the parameter value scales `QRTF(5)=NumRT` not `Q` (`IPVL=5`; a Fortran bug — `SFAC` scales `Q` at `ISCLOC=4`), and `NodQRT` is not copied on activation. From-scratch parameter authoring and `INSTANCES` raise `NotImplementedError`; `MXL`/consistency validated. **Stage 4.5A: the inline `TRANSIENTQ` transient extraction-flow time series is supported** (load → write → reload + from-scratch authoring via `transientq_*`; overrides `Q`, not recipients); `TRANSIENTQ`+`NPQRT>0` and external-unit `TRANSIENTQ` data fail explicitly. Registered in `MfUsg.load()`. See `USGT_STAGE4_04_PARAMETERS_QRT.md` / `USGT_STAGE4_05_QRT_TRANSIENTQ.md`. **Verified** |

---

## Stage 3 status audit (2026-05-31) — final state

Final, post-Stage-3 classification of every package that is not a straightforward
tested `✅ Full`. This is the resolved disposition (all Stage 3 cards have
landed), consistent with the coverage table above.

| Package | Final state | Basis |
|---------|-------------|-------|
| DIS | ✅ Full | Foundational structured grid; authored across the suite, loaded in structured Ex* models. |
| DISU | ✅ Full | Unstructured grid; authored/loaded across the suite; large-grid `free_format_npl` formatting round-trip tested (`test_usg.py`). |
| CLN | ✅ Full | Circular/rect/`GENERAL_SEC`/`PROCESSCCF`/`ISHAPE` authoring + round-trip tests; Ex3 conduit models. |
| DPF | ✅ Full | `FRAHK`/`IUZONTABIM`/conditional `SC2IM`/immobile-Richards authoring + round-trip tests. |
| PCB | ✅ Full | `(node, iSpec, conc)` authoring + round-trip; dataset order verified vs Fortran; exercised via Ex transport models. |
| OC | ✅ (intentionally not Full) | OC Fullness Card A (Stage 4): broad authoring/round-trip now tested (BOOTSTRAPPING header on line 1, per-SP bootstrap toggles, layer-qualified actions, DDREFERENCE, CONC/BUDGET, ATSA). Explicit gaps keep it not-Full: `SAVE IBOUND` solver-rejected by USG-T 2.7; FASTFORWARD/FASTFORWARDC + BOOTSTRAPPING execution not exe-verified; numeric→words rewrite. See `USGT_STAGE4_OC_FULLNESS.md`. |
| EVT | ✅ (intentionally not Full) | EVT Fullness Card B (Stage 4): all NEVTOP modes (struct + unstruct), transport IETFACTOR 0/<0/>0, per-SP reuse authored + round-trip tested; IETFACTOR round-trip bug fixed; exe smoke added. Explicit gaps keep it not-Full: ETS zonal time-series raises `NotImplementedError`; NPEVT params are expanded-valid-write. See `USGT_STAGE4_EVT_FULLNESS.md`. |
| MDT | ✅ (intentionally not Full) | MDT Fullness Card C (Stage 4): field order Fortran-audited; from-scratch authoring + round-trip tested for all main branches; fixed FRAHK/FRADARCY load, write_file(f=) crash, debug print; added option/length validation. Ex7 real-model round-trip/run kept. Gaps keeping it not-Full: species count uses MCOMP (chained-decay NTCOMP>MCOMP unverified); AI1/AI2 output binaries not read. See `USGT_STAGE4_MDT_FULLNESS.md`. |
| LAK | ✅ (intentionally not Full) | LAK Fullness Card D (Stage 4): from-scratch authoring + round-trip tested (no-transport, classic transport, TRANSPORTBOUNDARY MCOMP>1, TABLEINPUT); six authoring bugs fixed; validation added; Ex8 real-model round-trip/run kept. Gaps keeping it not-Full: sill/connectivity + multi-lake systems not authored from scratch; TABLEINPUT table contents external; GAGE separate. See `USGT_STAGE4_LAK_FULLNESS.md`. |
| GNC | ✅ | GNC unstructured ghost-node helper; verified. |
| TIB | ✅ Full (authoring) | Promoted in Stage 4.1: semantic `stress_period_data` authoring + `parse=True` load for flow (`NIB0/NIB1/NIBM1`) and transport (`NICB0/NICB1/NICBM1`) blocks; `load` default stays byte-exact raw round-trip, with raw fallback on unsupported `U1DINT` (`EXTERNAL`/`OPEN-CLOSE`). |
| GSF | ✅ Full (authoring) | Promoted in Stage 4.2: semantic `vertices`+`node_data` authoring + `parse=True` load (0-based internal / 1-based file), `from_grid` (requires per-vertex z), `to_grid()` via `from_gridspec`; `load` default stays raw round-trip. GSF is not solver input. |
| ETS | Parameter-preserving (ETSR) | Stage 4.4A: ETSR array parameters load → write → reload with syntax intact (`NPETS>0`, defs, per-period activation records, `INSTANCES`). Opt-in expanded fallback (`expand_parameters=True` → `NPETS=0`); from-scratch parameter authoring fails explicitly. See `USGT_STAGE4_04_PARAMETERS_ETS.md`. |
| HFB | Parameter-preserving (list) | Stage 4.4B: `NPHFB>0` list parameters load → write → reload with syntax intact (defs + `NLST` barrier rows + `NACTHFB` + active names). From-scratch authoring, `TRANSIENT_HFB`+`NPHFB>0`, and `INSTANCES` fail explicitly. See `USGT_STAGE4_04_PARAMETERS_HFB.md`. |
| DPT | Partial | `A-W_ADSORBIM` immobile air-water adsorption explicitly unsupported (fails before any shifted read); `DLIM` conditional correct. |
| SFR / STR / GAGE / FHB / SUB / SWT | Compatibility (SFR/FHB/GAGE) · guarded (STR/SUB/SWT) | Base MODFLOW-2005 classes; CLN is the project coupling. **SFR** kept enabled — DISU+SFR loads/writes/runs via base in `freyberg_usg` (`test_usg.py`). **FHB/GAGE** load via base in Ex8. **STR/SUB/SWT** guarded (Stage 4.6A): unstructured-`MfUsg` use raises `NotImplementedError`. |

Rule applied throughout: a package is labeled `✅ Full` only with a
Fortran-derived spec, from-scratch authoring tests, round-trip tests, and
explicit failure for unsupported modes. The `✅ (intentionally not Full)` rows
load/write correctly and are tested as noted, but are deliberately not promoted
(documented reason in the coverage table). No row is a pending TODO.

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

`NPDRT>0` (named parameters) are **preserved** (Stage 4.4D): definitions +
activations + per-row recipients round-trip; `MXADRT` is sized to the active
total per period; from-scratch authoring and `INSTANCES` raise
`NotImplementedError`. Activated SPREAD (`NR<0`) recipients are
structural-round-trip only, not execution-guaranteed (USG-T copies `DRTF` but
not `NodDRT` on activation). Structured (DIS) models delegate to the base
`ModflowDrt`. Return-flow node lists with `EXTERNAL`/`OPEN/CLOSE` control
records are not supported (explicit failure). See
`USGT_STAGE4_04_PARAMETERS_DRT.md`.

### Gap §8 — ETS: parameter preservation (RESOLVED for ETSR, Stage 4.4A)

USG-T ETS supports named parameters (`NPETS > 0`) for the ETSR array. FloPy now
**preserves** them: `MfUsgEts.load` keeps the parsed parameter definitions
(`self.parameters`) and the per-period active-parameter records
(`self.evtr_parm`), and `write_file` re-emits `NPETS>0` in item 2a, the
definition blocks, and the activation records (with `INSTANCES`), while
ETSS/ETSX/IETS/PXDP/PETM stay plain arrays. USG-T reads `NPETS` from item 2a
(`UPARARRAL` is called with `IN=-1`, so there is no `PARAMETER` line). The
shared write helper is `flopy/mfusg/_usgt_parameters.py`, reusing
`ModflowParBc` as the parser.

Remaining (keeps ETS not `Full`): from-scratch parameter *authoring*
(`MfUsgEts(..., npets>0)` with no loaded definitions) is unsupported and raises
`NotImplementedError`. The legacy expand path is still available opt-in via
`MfUsgEts.load(..., expand_parameters=True)` (writes `NPETS=0`). See
`USGT_STAGE4_04_PARAMETERS_ETS.md`.

### Gap §9 — HFB: parameterized barriers (RESOLVED via preservation, Stage 4.4B)

For non-parametric HFB, the implementation matches `gwf2hfb7u1.f`: static HFB
reads/writes `NHFBNP` rows once, and transient HFB reads `IHFBRD` for every stress
period, using `IHFBRD <= 0` as reuse/no-read and `IHFBRD > 0` as the flag to read
exactly `NHFBNP` rows.

Parameterized HFB (`NPHFB > 0`) is now **preserved**: `MfUsgHfb.load` keeps the
list-parameter definitions (`self.parameters`), the non-parametric barriers
(`self.hfb_data`), and the active-parameter names (`self.acthfb_names`); and
`write_file` re-emits the definition blocks, the non-parametric barriers, and the
`NACTHFB` + active-name records. Barrier rows are 0-based internally / 1-based on
file. Uses the shared list-parameter helpers (`UPARLSTRP`/`UPARLSTSUB`).

Remaining (keeps HFB not `Full`): from-scratch parameter authoring,
`TRANSIENT_HFB` combined with `NPHFB>0` (the Fortran would redefine parameters
each stress period under `ITERP=1`), and parameter `INSTANCES` (the Fortran
aborts) all raise `NotImplementedError`. See `USGT_STAGE4_04_PARAMETERS_HFB.md`.

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
`ITMP/-1` reuse. `AUTOFLOWREDUCE` is preserved as an option. `NPQRT>0` list
parameters are structurally preserved (Stage 4.4E). **Stage 4.5A: the inline
`TRANSIENTQ` transient extraction-flow time series is now supported** — loaded,
written, and authorable from scratch via the `transientq_*` attributes (it
overrides `QRTF(4)=Q` each step; recipients are untouched). `TRANSIENTQ`
combined with `NPQRT>0` and external-unit `TRANSIENTQ` data fail explicitly, as
do from-scratch parameter authoring and `INSTANCES`. Registered as `"qrt"` in
`MfUsg.load()`. See `USGT_STAGE4_05_QRT_TRANSIENTQ.md`.

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
