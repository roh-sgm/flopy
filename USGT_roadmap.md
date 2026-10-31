# USG-Transport FloPy — Package Coverage Reference

This document is the authoritative record of which USG-Transport 2.7 packages
and features are supported in this fork, which are partial, and which are absent.

All entries marked **verified** were cross-referenced against the actual
USG-Transport 2.7.0 Fortran source (`USGT_V_2-7-0_Source_Code/`) in May 2026,
in addition to the I/O manual. Entries without that mark rely on the I/O manual
and/or round-trip testing only.

Stage 5 entries were audited against the USG-Transport 2.8.1 source.
Executables are named only where an entry states them; no executable run is
implied for the others.

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
| DIS | `DIS` | `MfUsgDis` | ✅ Full | `mfusg.f` | Structured grid discretization. Authored in essentially every synthetic test and loaded in the structured Ex* models. **Verified** G06: stress periods `SS`/`TR`/`TRTOSS stofrac` via `steady` + `trtoss` (Gap §11). |
| DISU | `DISU` | `MfUsgDisU` | ✅ Full | `mfusg.f` | Unstructured grid. Authored/loaded across the suite; `free_format_npl=10` default prevents large-grid buffer overflow, with dedicated round-trip formatting tests in `autotest/test_usg.py` (`test_free_format_npl`). **Verified** G06: stress periods `SS`/`TR`/`TRTOSS stofrac` via `steady` + `trtoss` (Gap §11). |
| BCF6 | `BCF6` | `MfUsgBcf` | ✅ | `gwf2bcf-lpf-u1.f` | TABRICH items 1c (`IUZONTAB`) and 1d (`RETCRVS`, shape `(nuzones, nutabrows, 3)` = caphead/saturation/relperm) now authored/loaded/written; incomplete TABRICH write fails explicitly. See Gap §2 (resolved). Stage 5.02-B2 (checked against USG-T 2.8.1 `READBCF`/`SGWF2BCFU1S`/`SGWF2BCFU1G`): LAYCON codes `10*INTERCELLT+LAYCON` written as one integer (free) or `40I2` (fixed, 40 per line); fixed-format item 1 written and loaded with `I10,F10.3,I10,F10.3,4I10` (IKCFLAG only on DISU; F10.3 reals written with a decimal point in the mantissa, e.g. `-1.E+30`, because an exponent alone keeps the three implied decimals, and read with the same rule; checked against an F10.3 read compiled with gfortran and the USG-T listing echo), values checked to fit before the file is opened; fixed-format files with header options are rejected; TRPY only with `IKCFLAG=0`, DISU ANGLEX on `IHANISO`, connection array for any `IKCFLAG!=0`, IKCFLAG treated as 0 on structured grids; SF2, TRAN/HY and WETDRY follow the effective LAYCON (BAS `RICHARDS`/`RICHARDS_HP` promotion); new optional `alpha`/`beta`/`sr`/`brook`/`bp` arrays (appended to the signature); KV for every layer with `IKVFLAG!=0`; `write_file(f=handle)` leaves the handle open. Tests: writer order, hand-written loads (free/fixed, `OPEN/CLOSE` + shared `EXTERNAL`), round-trips, executable comparisons with sensitivity controls, and Ex3/Ex4/Ex8 rewrite = original heads. Not Full: DISU `IDSYMRD=1`, option detection when loading fixed-format files, free-format header precision, and the engine abort seen with transient BCF models whose only convertible layers are LAYCON 5 (cause not isolated). **Verified** |
| LPF | `LPF` | `MfUsgLpf` | ✅ | `gwf2bcf-lpf-u1.f` | TABRICH 1c/1d implemented (RETCRVS replaces per-layer alpha/beta/sr/brook, which are skipped under TABRICH); parse bug fixed (nutabrows token, int cast). Stage 5.02-B1 (checked against USG-T 2.8.1 `READLPF`/`SGWF2LPFU1S`/`SGWF2LPFU1G`): Richards arrays follow the effective LAYCON (=5, with BAS `RICHARDS`/`RICHARDS_HP` promoting convertible layers) instead of every layer; `WETDRY` is loaded whenever `LAYWET!=0`; `VKCB`/`WETDRY` are written in steady-state models; `SY` uses `LAYTYP` after the `THICKSTRT` reset; the DISU connection array is written/loaded for any `IKCFLAG!=0`; item 7 is present when any `LAYWET!=0`; `write_file(f=handle)` works and leaves the handle open. Tests: writer order, hand-written loads (incl. `OPEN/CLOSE` + shared `EXTERNAL`), round-trips, and executable comparisons (FloPy-written = hand-written, with sensitivity controls). Not Full: DISU `IDSYMRD=1` connection-array lengths and LPF named parameters (`NPLPF`) are not audited. **Verified** |
| SMS | `SMS` | `MfUsgSms` | ⚠️ Partial | `glo2sms-u1.f` | Records 1a/1b, item 2, XMD and PCGU. G04: SHIFT and the other record 1b options kept in any case, NUMTRACK=0 records, URWORD REAL syntax and REAL-preserving output (Gap §12). Not supported: IMS (LINMETH=5), SOLVEROUTPUTFILE (dropped with a warning), STRESSPERIODSETTINGS (rejected). |
| OC | `OC` | `MfUsgOc` | ✅ | `glo2basu1.f` (SGWF2BAS7I/J/N) | OC Fullness Card A (Stage 4): from-scratch authoring + round-trip tested for `ATSA`/`NPTIMES`/`NPSTPS`, the `BOOTSTRAPPING` header (now written on the first OC line where USG-T parses it, and preserved on load), per-SP `BOOTSTRAP`/`NOBOOTSTRAP`/`BOOTSTRAPSCALE`/`NOBOOTSTRAPSCALE`, layer-qualified `PRINT`/`SAVE HEAD`/`DRAWDOWN`/`CONC`, `DDREFERENCE`, `SAVE`/`PRINT CONC`/`BUDGET`, `COMPACT BUDGET AUX`. **Stage 4.7B: `ATSA` adaptive time-stepping is now exe-verified** (a transient model with `MfUsgOc(atsa=1)` runs under USG-T 2.7 and reaches the analytical gradient; the Fortran overrides `NSTP` and uses built-in ATS constants). Kept `✅` not `Full` — explicit gaps: `SAVE IBOUND` is commented out in USG-T 2.7 (solver rejects it; FloPy preserves the keyword and `check()` warns about it); `FASTFORWARD`/`FASTFORWARDC` and `BOOTSTRAPPING` execution stay **manual tier** (they read GWF/CLN/DDF heads or concentrations from an external prior-run file); numeric-format OC rewrites as words. G05: TIMOT read and written as USG-T reads it, ATS parameters kept, numeric ATS records converted or rejected (Gap §13). See `USGT_STAGE4_OC_FULLNESS.md` / `USGT_STAGE4_18_EXE_VERIFICATION.md`. |
| CHD | `CHD` | `MfUsgChd` | ✅ Full | `gwf2chd7u1.f` | Node-based unstructured. Internal nodes 0-based, file I/O 1-based. AUX transport concentrations. Programmatic authoring tested. float64 precision. **Verified** |
| WEL | `WEL` | `MfUsgWel` | ✅ Full | `gwf2wel7u1.f` | Rates on GWF and CLN nodes, AUX transport concentrations, `ITMPCLN`, and programmatic authoring tested. CLN connectivity is handled by `MfUsgCln`. **Verified** |
| DRN | `DRN` | `MfUsgDrn` | ✅ Full | `gwf2drn7u1.f` | Node-based unstructured. Internal nodes 0-based, file I/O 1-based. AUX transport concentrations. Programmatic authoring tested. **Verified** |
| RIV | `RIV` | `MfUsgRiv` | ✅ | `gwf2riv7u1.f` | Node-based, internal nodes 0-based, file I/O 1-based, AUX concentrations. G02 (USG-T 2.8.1 `GWF2RIV7U1AR`): the DISU header is `MXACTR IRIVCB [options]` with IRIVCB = `ipakcb` (it was written from `irdflag`, 50 by default, so the user's CBC unit was lost); load keeps IRIVCB and the CBC output entry. `irdflag` stays in the signature but is ignored (warns if it differs from `ipakcb`). No undeclared `irch` column: USG-T reads node, stage, cond, rbot and the declared AUX values in declaration order (ULSTRDU); trailing values are ignored on load with a warning, and an `irch` field that is not declared as AUX is rejected before writing. A variable declared as `AUX IRCH`/`AUXILIARY IRCH` (any case, whole-token match) is an ordinary AUX variable and is kept on load, write and reload; FloPy does not declare it automatically. Executable checks (DISU): CBC on the requested unit (55, 61) = hand-written file, Q = C(stage−h); IRIVCB 0 and −1 controls; hand-written `AUX IRCH` file = FloPy load+rewrite (heads, RIV flow, CBC AUX value). Not Full: DISU writer ignores fixed format (IFREFM=0); PARAMETER/SFAC/OPEN-CLOSE/EXTERNAL lists and BINARY input not supported (G15/F5k/F5i). |
| GHB | `GHB` | `MfUsgGhb` | ✅ Full | `gwf2ghb7u1.f` | Node-based, internal nodes 0-based, file I/O 1-based, AUX transport concentrations. Programmatic authoring tested. **Verified** |
| RCH | `RCH` | `MfUsgRch` | ⚠️ Partial | `gwf2rch8u1.f` | INRCHZONES crash bug fixed (2026-05-20). INCONC/INIRCH spacing fixed. G01 (partial; source USG-T 2.8.1, run with 2.8.1 and 2.7.0): the `selev`/`iznrch` arrays are kept and the RTS zone count is written as an integer (F39/F40); SEEPELEV and one-zone RTS models match hand-written controls. Open: IRCHCONC per component (F41), INRCHZONES/INSELEV flags fixed at 1, so multi-zone RTS rewrites are wrong (F42), IZNRCH without RTS (F43), load with INSELEV -1/0, DISU NRCHOP=2, INCONC given once. Partial is a correction of the documented status, not a loss of function. |
| RTS | `RTS` | — (name-file entry kept) | ⚠️ Partial | `glo2basu1.f` (SGWF2BAS8OPEN), `gwf2rch8u1.f` | Recharge time-series file read by RCH with the `RTS` option (IURTS = IUNIT(9)). G03: the name-file entry is kept as a reference to an opaque file; see Gap §10. The file contents are not read. Since G01 (F39) RCH with the `RTS` option is written; a one-zone RTS model matches its control with USG-T 2.8.1 and 2.7.0, while multi-zone rewrites still write `INRCHZONES 1` (F42, see RCH). Not Full. |
| EVT | `EVT` | `MfUsgEvt` | ✅ | `gwf2evt8u1.f` | EVT Fullness Card B (Stage 4): from-scratch authoring + round-trip tested for `NEVTOP`=1/2/3 (structured + unstructured `NEVTOP=2` with `MXNDEVT`; `IEVT` 0-based internal / 1-based file, layer-range (structured) and node-range (unstructured) validated), transport `IETFACTOR` 0/<0/>0 with per-`MCOMP` `ETFACTOR` (scalar or array), and per-SP reuse (`-1`). Fixed a real round-trip bug (`load` dropped `IETFACTOR`) and made the writer emit 3 dataset-1 ints whenever transport is active. Plus a USG-T 2.7 executable smoke (`test_usgt_exe_evt_from_scratch`). **Stage 4.6F-A: `NPEVT` EVTR array parameters are now preserved (load → write → reload) and authored from scratch** — optional `PARAMETER NPEVT` line + `UPARARRRP` defs (PTYP='EVT', `INSTANCES`) + per-SP `UPARARRSUB2` activations; `parameters={name: {parval, clusters \| instances}}`/`evtr_parm`, auto `NPEVT`, first period must activate, validated before opening. `expand_parameters=True` keeps the legacy expanded `NPEVT=0` write. Reuses the shared ETS array-parameter helpers. **Stage 4.6F-B: ETS zonal time-series (`ETS MXZNEVT`/`INEVTZONES`) is deferred** — the Fortran audit shows it is an ATS-coupled dynamic execution mode (requires `IATS>0`, else `STOP`; an external `IUETS` time-series read progressively during the run that supersedes the EVTR array), not static I/O; it now fails explicitly **per branch** (authoring + load of `ETS MXZNEVT` and of a per-SP `INEVTZONES`), no partial parse, with the full Fortran spec in `USGT_STAGE4_16_EVT_ETS_ZONAL.md`. Kept `✅` not `Full` for that one honest gap. See `USGT_STAGE4_EVT_FULLNESS.md`. **Verified** |
| EVS | `EVS` | — (name-file entry kept) | ⚠️ Partial | `glo2basu1.f` (SGWF2BAS8OPEN), `gwf2evt8u1.f` | ET time-series file read by EVT with the `ETS` option (IUETS = IUNIT(6)). G03: the name-file entry is kept by `MfUsg.load`/`write_input` as a reference to an opaque file (contents not read); see Gap §10. EVT with `ETS` still fails to load (`NotImplementedError`), so a model using EVS is rewritten only with EVT skipped by `load_only`. Not Full. |
| ETS | `ETS` | `MfUsgEts` | ⚠️ Partial | `gwf2ets8u1.f` | Authoring tested for NETSEG=1, NETSEG>1 (PXDP/PETM), NETSOP=2 (IEVT), and IESFACTOR transport flag. **Stage 4.4A: ETSR array parameters are now preserved** (load → write → reload keeps `NPETS>0` in item 2a, the definition blocks, and per-period activation records, incl. `INSTANCES`; ETSS/ETSX/IETS/PXDP/PETM stay plain arrays — the mix the Fortran allows). USG-T reads `NPETS` from item 2a (`UPARARRAL` called with `IN=-1`), so no `PARAMETER` line is written. Opt-in expanded fallback (`expand_parameters=True`) still writes `NPETS=0`. **Stage 4.6C-D: from-scratch ETSR array-parameter authoring is now supported** — pass `parameters={name: {parval, clusters | instances}}` + `evtr_parm={kper: [(name, instance_or_None), ...]}` (ergonomic dict built into a `ModflowParBc` via `build_array_parameter_bc_parms`; `NPETS`/`nclu` auto; partyp defaults/validated `ets`; first period must activate); normalize + validate before opening (`ValueError`, no partial file). Shared write helper `flopy/mfusg/_usgt_parameters.py` reuses `ModflowParBc` (no second parser). Still not `Full` (only ETSR parameterized; parametric execution not USG-T smoke-tested). See `USGT_STAGE4_04_PARAMETERS_ETS.md` / `USGT_STAGE4_14_ETS_PARAMETER_AUTHORING.md`. **Verified** |
| HFB | `HFB6` | `MfUsgHfb` | ⚠️ Partial | `gwf2hfb7u1.f` | Node-based. Non-parametric static, structured static, and `TRANSIENT_HFB` IHFBRD = >0/0/-1 semantics implemented and tested. **Stage 4.4B: HFB list parameters (`NPHFB>0`) are now preserved** (load → write → reload keeps the definition blocks with `NLST` barrier rows, `NACTHFB`, and the active names; barriers 0-based internal / 1-based file; parameter-defined and non-parametric barriers can mix). Uses the shared list-parameter helpers in `flopy/mfusg/_usgt_parameters.py` (`UPARLSTRP`/`UPARLSTSUB`). Barrier lists (parameter `NLST` rows and non-parametric `NHFBNP` rows) also support leading `SFAC`/`OPEN/CLOSE`/`EXTERNAL` list controls via the shared `_usgt_list.begin_list_block` (list-control follow-up; `SFAC` scales `HYDCHR`, `EXTERNAL` via `ext_unit_dict` else `NotImplementedError`). **Stage 4.6C-A: from-scratch `NPHFB>0` parameter authoring is now supported** — pass `parameters={name: {parval, data, ...}}` + `acthfb_names` (ergonomic input: `data` array-like, `nphfb`/`mxfb`/`nlst`/`nacthfb` auto; partyp defaults/validated `HFB`; structured + unstructured); normalize + validate before opening (`ValueError`, no partial file). Still `⚠️ Partial` / not `Full`: `TRANSIENT_HFB`+`NPHFB>0` and parameter `INSTANCES` raise `NotImplementedError` (the last matches the Fortran). See `USGT_STAGE4_04_PARAMETERS_HFB.md` / `USGT_STAGE4_11_HFB_PARAMETER_AUTHORING.md`. **Verified** |
| GNC | `GNC` | `MfUsgGnc` | ✅ | `disu2gncn1.f` | **Verified** |
| LAK | `LAK` | `MfUsgLak` | ✅ | `gwf2lak7u1.f` | LAK Fullness Card D (Stage 4): from-scratch authoring + round-trip tested for no-transport, classic transport (`CPPT`/`CRNF`), `TRANSPORTBOUNDARY` (one `CLAKE(1:NSOL)` line per lake, MCOMP>1), and `TABLEINPUT` (per-lake tab unit + external registration). Fixed six authoring bugs (conc_data mis-assignment; `flux_data=None` crash; conc_data accessed without transport; TRANSPORTBOUNDARY 9b written per-component instead of one-line-per-lake; load stored conc as strings; `transportboundary` flag not synced to the header keyword) and added validation (flux_data required; TRANSPORTBOUNDARY needs transport; clake nlakes×mcomp; transport needs conc_data). A review follow-up further validates authoring inputs in `__init__` so malformed input fails with a clear `ValueError` instead of crashing in `write_file()`: TABLEINPUT requires exactly one `tab_file`/`tab_unit` per lake; `conc_data` is checked per written period (every `(lake, component)` present; classic 2 values for `WTHDRW>=0` / 3 for `WTHDRW<0`; TRANSPORTBOUNDARY a single value); `flux_data` requires one dataset-9a entry per lake. `Ex8_Lake` real-model round-trip/run kept. **Stage 4.6D: multi-lake + sill/connectivity (datasets 7/8) are now authored from scratch and validated** — `sill_data={kper: [([IC, lake1..lakeIC], [sill1..sill_{IC-1}]), ...]}` (lake numbers 1-based, center first; `IC>=2`); `_validate_sill_data` checks `kper` range, `ITMP>0` (the Fortran only reads ds 7/8 when `LKARR`/`BDLKNC` are (re)specified), IC/list/sillvt consistency, and lake-id range/uniqueness before opening the file. Two from-scratch round-trip tests (no-transport and classic `mcomp=1`) + 7 negatives. Kept `✅` not `Full` — remaining gaps: `TABLEINPUT` bathymetry table **contents** are external; GAGE coupling separate; multi-lake-with-sill from-scratch *execution* is not exe-smoke-tested (manual tier; `Ex8_Lake` covers LAK execution). See `USGT_STAGE4_LAK_FULLNESS.md` / `USGT_STAGE4_15_LAK_MULTILAKE_CONNECTIVITY.md`. |
| CLN | `CLN` | `MfUsgCln` | ✅ Full | `cln2basu1.f`, `cln2props1.f` | From-scratch authoring + round-trip tested for `PROCESSCCF`/`ICLNGWCB`, the 9-field `ISHAPE` node records, and `GENERAL_SEC` tabular shapes; also exercised by the Ex3 CLN conduit real models. Promoted to Full (Stage 3 Card 8). **Verified** |
| DDF | `DDF` | `MfUsgDdf` | ✅ Full | `density.f` | RHOFRESH/RHOSTD/CSTD/ITHICKAV/IMPHDD/ISHARP + NONLINEAR table. From-scratch NONLINEAR-table authoring test added (2026-05-30). **Verified** |
| BCT | `BCT` | `MfUsgBct` | ✅ Full | `glo2btnu1.f` | IDISP=1 and IDISP=2 (DLX/DLY/DLZ/DTXY/DTYZ/DTXZ), all 19 item 1a fields, A-W_ADSORB, ICHAIN, ISPRCT, ISOLUBILITY, IMULTI, IMASSWR options. From-scratch authoring tests (1-species, IDISP=2, multi-species) added 2026-05-30. **Verified** |
| PCB | `PCB` | `MfUsgPcb` | ✅ Full | `gwt2bndsu1.f` | Prescribed Concentration Boundary. Node-based `(node, iSpec, conc)`; 0-based internal / 1-based file. Field order verified vs Fortran; from-scratch authoring + round-trip tested (Stage 3 Card 8); exercised via the Ex transport models. **Verified** |
| MDT | `MDT` | `MfUsgMdt` | ✅ | `gwt2mdtu1.for` | Matrix Diffusion Transport. MDT Fullness Card C (Stage 4): field order Fortran-audited and from-scratch authoring + round-trip tested for the header options (`FRAHK`/`FRADARCY`/`TSHIFTMD`/`SEPARATE_AI2`/`MULTIFILE_MD`, only when `IDPF==0`), base arrays (`MDFLAG`/`VOLFRACMD`/`PORMD`/`RHOBMD`/`DIFFLENMD`/`TORTMD`; `VOLFRACMD` skipped when `IDPF!=0`), per-species `KDMD`/`DECAYMD`/`YIELDMD`/`DIFFMD`, and `AIOLD1MD`/`AIOLD2MD` under `TSHIFTMD > 1e-10` (the Fortran threshold, applied consistently; TSHIFTMD written in general format so small valid values are not rounded to 0). Fixed real bugs (`load` dropped `FRAHK`/`FRADARCY` via a case mismatch; `write_file(f=handle)` crashed; stray debug print) and added validation (FRAHK⊕FRADARCY, IDPF-options, MULTIFILE needs `imdtcf>0`, per-`MCOMP` list lengths). Still **round-trip-run via the three Ex7 real models**. **Stage 4.7B: MDT execution is covered by those Ex7 models**; the remaining gaps are **manual tier** — chained-decay `NTCOMP>MCOMP` (the array reads loop over `NTCOMP`, but exercising it needs a tailored BCT chained-decay setup) and the `SEPARATE_AI2`/`MULTIFILE_MD` per-component binary AI-output files. Kept `✅` not `Full`. See `USGT_STAGE4_MDT_FULLNESS.md` / `USGT_STAGE4_18_EXE_VERIFICATION.md`. |
| DPF | `DPF` | `MfUsgDpf` | ✅ | `gwf2dpf1u1.f` | Stage 5.01 (checked against USG-T 2.8.1 `GWF2DPFU1AR`, lines 111–191, identical in 2.7.0): `PHIF`…`BPIM` are written and loaded layer by layer (`DO 200`), after the complete `IUZONTABIM`, `IBOUNDIM` and `HNEWIM` blocks. Earlier multi-layer files were written property by property, which USG-T mis-reads (and FloPy itself mis-read for ≥2 Richards layers); they are not auto-detected and must be regenerated from the original data. Conditional reads follow the effective LAYCON (LPF `LAYTYP` mapping incl. `THICKSTRT`; BAS `RICHARDS`/`RICHARDS_HP` turn every LAYCON≠0 into 5); `TABRICH`/`BUBBLEPT` are taken from BCF or LPF; the default unit is passed to the package so `write_input()` works. Tests: layer order, manual Fortran-order load, two Richards layers ±`BUBBLEPT`, mixed layers ±`TABRICH`, DISU `nodelay=[3,2]`, `OPEN/CLOSE` + shared `EXTERNAL`, name-file unit; executable check (confined two-layer model): FloPy-written DPF = hand-written Fortran-order DPF (atol 1e-6). Not Full: Richards (LAYCON=5) DPF has comparative executable evidence only (FloPy-written and hand-written DPF gave identical heads with an all-LAYCON-5 LPF model, with and without `BUBBLEPT`). The 9–13 % flow-budget discrepancy of those runs was traced, for that small case with USG-T 2.8.1, to the immobile storage of the CHD cells, which the CONSTANT HEAD term does not include (see Gap §5). D7 stage 1 (USG-T 2.8.1 `U1DREL8`): with BAS `DPIN`/`DPIO` `HNEWIM` is float64 (17-digit text, REAL `CNSTNT`/`CONSTANT`, binary loaded only through a negative unit with the 52-byte header and no multiplier). D7 stage 2: binary `HNEWIM` output through `EXTERNAL -n` and a `DATA(BINARY)` file of the package (REAL*8 with `DPIN`/`DPIO`, REAL otherwise), without replacing files or name-file entries of other components (see Gap §5). **Verified** |
| DPT | `DPT` | `MfUsgDpt` | ⚠️ Partial | `gwt2dptu1.f` | DLIM condition checks both IDPF and IDISPIM. **Stage 4.6E: immobile-domain air-water adsorption (`A-W_ADSORBIM`) is now authored/loaded/written** for the array-only function-index branches — `IAREA_FNIM ∈ {1 (AMAX), 4 (X2/X1/X0)}` × `IKAWI_FNIM ∈ {1, 2}` (Langmuir A/B per species). Reads the option-line indices, RP1 area arrays (after heat), and RP2 ALANG/BLANG per species (before ADSORBIM), matching `dpt2aw_adsorb.f`. The scalar/tabular branches `IAREA_FNIM ∈ {2,3,5}` and `IKAWI_FNIM ∈ {3,4}` (zone map + `ROG_SIGMA`/`SIGMA_RT` + tables) raise a specific `NotImplementedError` (authoring + load, before any array read). See `USGT_STAGE4_06_DPT_AW_ADSORBIM.md`. **Verified** |
| TIB | `TIB` | `MfUsgTib` | ✅ Full (authoring) | `glo2basu1.f` | Transient Ibound. Semantic `stress_period_data` authoring + `parse=True` load covering the full `GWF2TIB1RP` grammar: flow blocks (`NIB0` inactivate / `NIB1` activate / `NIBM1` prescribed-head, with `HEAD`/`AVHEAD`/bare-reuse) and transport blocks (`NICB0`/`NICB1`/`NICBM1`, with multi-component `CONC`/`AVCONC`/bare); `NIB0`/`NICB0` node lists via `U1DINT`; 3- vs 6-int header keyed to BCT presence; 0-based internal / 1-based file. From-scratch authoring, semantic write/reload, multi-node `U1DINT`, transport-conc, and explicit-rejection tests. `MfUsgTib.load` defaults to byte-exact raw round-trip (what `MfUsg.load` uses); `load(parse=True)` returns semantic data and falls back to raw on `U1DINT` `EXTERNAL`/`OPEN-CLOSE` rather than writing partial data. **Verified (Stage 4.1)** |
| TVM | `TVM` | `MfUsgTvm` | ✅ Full | `tvmu2.f` | Full semantic implementation (2026-05-20). HK/VKA/SS/SY/DDFTR/POR; transport-aware field counts; nper+1 boundaries; 20 autotests pass. **Verified** |
| GSF | `GSF` | `MfUsgGsf` | ✅ Full (authoring) | — | Grid Specification File (not solver input — consumed by post-processors / `UnstructuredGrid.from_gridspec`, so no Fortran reader). Stage 4.2: semantic `vertices` + `node_data` authoring and `parse=True` load — header (`UNSTRUCTURED`/`UNSTRUCTURED GWF`, validated), `nnodes`/`nlay`, vertex `(x,y,z)`, per-node `(node, xc, yc, zc, layer, vertices)`; 0-based internal nodes/vertices/layers, 1-based file ids; trailing line-2 gridgen flags preserved (default `(1,1)`). `from_grid(model, grid, top_zverts, bot_zverts)` writes the USG-T top/bottom doubled-vertex convention (single-surface `zverts` kept as legacy; per-vertex z required — the grid does not retain it); `from_disv_gridprops(model, disv_gridprops, top, botm, skip_degenerate=)` maps a MF6 DISV 2D template to a single-layer GSF (closing-vertex removal, degenerate-cell skip). Both accept `vertex_mode="shared"`/`"parsimonious"` (default; neighbours reuse vertex ids) or `"cell"`/`"nonparsimonious"` (unique top/bottom vertices per cell, 8 per quad) — the two GRIDGEN2GSF layouts (`gridgen2gsf.f90`, secondary ref; primary spec gwutil_a 2.17, line 2 = `NNODES NLAY IZ IC`). `to_grid()` delegates to `from_gridspec` (which now correctly accepts `UNSTRUCTURED GWF`). `load` defaults to byte-faithful raw round-trip; `parse=True` is strict (header exactly `UNSTRUCTURED`/`UNSTRUCTURED GWF`; `IZ IC` omitted or `1 1`; `inode` = `1..nnode` in order; trailing content) → raw fallback otherwise. Constructor validates header / `IZ IC`=`(1,1)` / contiguous-ordered node ids `0..nnodes-1` / `nlay >= max(layer)+1`. From-scratch authoring, write/reload, 0-based↔1-based, `to_grid`, invalid-ref, mode-exclusivity, top/bottom, DISV, vertex-mode, strict-header, IZ-IC, and inode tests. Companion utility `gridgen_to_gsf` (`flopy/mfusg/gridgen2gsf.py`, Stage 4.3) builds a GSF from a flopy `Gridgen` / `disv_gridprops` / `UnstructuredGrid` source (a non-interactive helper inspired by `GRIDGEN2GSF`; parsimonious mode compacts unused vertices; delegates to `MfUsgGsf`). **Verified (Stage 4.2 + follow-ups; Stage 4.3 utility)** |
| SFR | `SFR` | base `ModflowSfr2` | Compatibility (DISU-validated) | `gwf2sfr7u1.f` | Base MODFLOW-2005 class, **kept enabled**: `examples/data/freyberg_usg` is a DISU model that uses SFR and is loaded → written → run via the base class in `autotest/test_usg.py` (`test_freyburg_usg`/`_external`), so unstructured SFR is empirically supported. **Not** promoted to Full USG-T transport (no transport-coupled authoring/tests); **not guarded** (Stage 4.6A deliberately leaves SFR on the base class). CLN remains the preferred coupling. |
| STR | `STR` | `MfUsgStr` (guards base `ModflowStr`) | Compatibility-only guarded | `gwf2str7u1.f` | **Stage 4.6A:** USG-T DISU layout unvalidated; using STR on an unstructured `MfUsg` model (load **or** authoring) raises `NotImplementedError` before any read/write rather than silently producing a non-USG-T file. Structured `MfUsg` and plain `flopy.modflow.Modflow` are unaffected (delegate to / use base `ModflowStr`). Not used by any target/example model. |
| GAG | `GAGE` | base `ModflowGage` | Compatibility-only | `gwf2gag7u1.f` | Base class. **Loads via base in the Ex8 real model** (round-trip exercised there); USG-T-specific records not independently validated. Compatibility-only (Card 6). |
| FHB | `FHB` | base `ModflowFhb` | Compatibility-only | `gwf2fhb7u1.f` | Base class. **Loads via base in the Ex8 real model** (`m.fhb` asserted; round-trip exercised there); USG-T unstructured format not independently validated. Compatibility-only (Card 6). |
| DRT | `DRT` | `MfUsgDrt` | ✅ Full (authoring) / Parameter-preserving (NPDRT, incl. activations + recipients) / Expanded valid write (list controls) | `gwf2drt8u.f` | Node-based DRT8. EL+COND, RETURNFLOW (inline single recipient `NR>0` or `SPREAD` multi-node `NR<0` via U1DINT), `CHANGEC`/`IDCHNGTYP` transport, AUX, `ITMP/-1` reuse; 0-based internal / 1-based file. Reads `SFAC`/`OPEN-CLOSE`/`EXTERNAL` list controls (expanded inline on write). `recipient_nodes` validated against the stress list. **Stage 4.5B: the spreading-recipient `U1DINT` block (`NR<0`) also reads `EXTERNAL <unit>` / `OPEN-CLOSE <fname>` (resolved via `ext_unit_dict`/`model_ws`, quote-aware, **free-format `(FREE)` only** — a fixed `FMTIN` raises `NotImplementedError`; expanded inline on write; unresolved `EXTERNAL` raises `NotImplementedError`).** **Stage 4.4D: `NPDRT>0` list parameters are now preserved** (load → write → reload keeps `NPDRT`/`MXL` in item 1, definitions with `NLST` drain rows + their RETURNFLOW recipients/spreading blocks, and per-SP `ITMP NP` + active names). DRT is internally consistent — definition and activation both use `PARTYP='DRT'` (`gwf2drt8u.f:135`/`:1108`), so active parameters are valid (unlike SGB); parameter value scales `COND` (`IPVL=5`). **Stage 4.6B: from-scratch `NPDRT>0` parameter authoring is now supported** — pass `parameters={name: {parval, data, ...}}` + `active_params={kper: [...]}` (ergonomic input: `data` array-like, `nlst`/`MXL`/`recipient_nodes` auto-computed; partyp defaults/validated `DRT`); normalize + validate before opening (ValueError, no partial file). `INSTANCES` still raises `NotImplementedError`; an activated `SPREAD` (`NR<0`) parameter round-trips structurally but is not execution-guaranteed (`NodDRT` not copied on activation). See `USGT_STAGE4_10_DRT_PARAMETER_AUTHORING.md`. **Unstructured only** — structured construction raises (use `ModflowDrt`); structured load delegates to base. See `USGT_STAGE4_04_PARAMETERS_DRT.md`. **Verified** |
| SUB | `SUB` | `MfUsgSub` (guards base `ModflowSub`) | Compatibility-only guarded | `gwf2sub7u1.f` | **Stage 4.6A:** USG-T items unvalidated; using SUB on an unstructured `MfUsg` model (load or authoring) raises `NotImplementedError` before any read/write. Structured `MfUsg` and plain `flopy.modflow.Modflow` unaffected. Not used by any target/example model. |
| SWT | `SWT` | `MfUsgSwt` (guards base `ModflowSwt`) | Compatibility-only guarded | — | **Stage 4.6A:** USG-T format unvalidated; using SWT on an unstructured `MfUsg` model (load or authoring) raises `NotImplementedError` before any read/write. Structured `MfUsg` and plain `flopy.modflow.Modflow` unaffected. Not used by any target/example model. |
| SGB | `SGB` | `MfUsgSgb` | ✅ Full (authoring) / Parameter-definition preserving + from-scratch authoring (NPSGB, no activations) / Expanded valid write (list controls) | `glo2sgbu1.f` | Specified Gradient Boundary. Node-based `(node, gradient)` list, AUX transport concentrations, `ITMP/-1` reuse; 0-based internal / 1-based file. Reads `SFAC` (inert on gradient per Fortran ISCLOC), `OPEN-CLOSE`, `EXTERNAL` list controls (expanded inline on write). **Stage 4.4C: `NPSGB>0` parameter *definitions* are preserved** (`PARAMETER NPSGB MXS` + definitions with `NLST` `NODE GRADIENT [aux]` rows; shared helpers in `_usgt_parameters.py`). **Stage 4.6C-B: from-scratch *definition* authoring** — pass `parameters={name: {parval, data, ...}}` (ergonomic input: `data` array-like, `MXS`/`nlst` auto, partyp defaults/validated `SGB`, AUX via dtype); normalize + validate before opening (`ValueError`, no partial file); `NP=0` written every period. **Active SGB parameters are unsupported**: USG-T 2.7 defines them as `PARTYP='SGB'` (`UPARLSTRP`, glo2sgbu1.f:97) but activates them as `PTYP='G'` (`UPARLSTSUB`, glo2sgbu1.f:185) — a type conflict (`parutl7.f:684/800`) that aborts the run — so per-SP `NP>0`, `active_params`, and `INSTANCES` raise `NotImplementedError`. Registered in `MfUsg.load()`. See `USGT_STAGE4_04_PARAMETERS_SGB.md` / `USGT_STAGE4_12_SGB_PARAMETER_AUTHORING.md`. **Verified** |
| QRT | `QRT` | `MfUsgQrt` | ✅ Full (authoring) / Parameter structural-preserving (NPQRT, not execution-guaranteed) / Expanded valid write (list controls) | `gwf2QRT8u.f` | Sink with Return Flow. Node-based `(node, q, rfprop)` + variable-length recipient-node lists (`NodQRT` via U1DINT), `CHANGEC`/`IQCHNGTYP` transport, AUX, `ITMP/-1` reuse. Reads `SFAC` (scales Q)/`OPEN-CLOSE`/`EXTERNAL` (expanded inline on write); `recipient_nodes` validated against the stress list. **Stage 4.5B: the recipient `U1DINT` block also reads `EXTERNAL <unit>` / `OPEN-CLOSE <fname>` (resolved via `ext_unit_dict`/`model_ws`, quote-aware, **free-format `(FREE)` only** — a fixed `FMTIN` raises `NotImplementedError`; expanded inline on write; unresolved `EXTERNAL` raises `NotImplementedError`).** `AUTOFLOWREDUCE` preserved. **Stage 4.4E: `NPQRT>0` list parameters are structurally preserved** (load → write → reload of definitions + recipient blocks + per-SP `ITMP NP` activations; `MXAQRT` covers the active total; `MXRTCELLS` reflects definition recipients). QRT is type-consistent (`PARTYP='QRT'` both sides, `gwf2QRT8u.f:183`/`:1118`), but **not execution-guaranteed**: the parameter value scales `QRTF(5)=NumRT` not `Q` (`IPVL=5`; a Fortran bug — `SFAC` scales `Q` at `ISCLOC=4`), and `NodQRT` is not copied on activation. **Stage 4.6C-C: from-scratch `NPQRT>0` parameter authoring is now supported** — pass `parameters={name: {parval, data, ...}}` + `active_params={kper: [...]}` (ergonomic input: `data` array-like, `MXL`/`nlst`/`recipient_nodes` auto, partyp defaults/validated `QRT`); normalize + validate before opening (`ValueError`, no partial file). **Structural only** — an activated parameter is not execution-guaranteed (caveat above). `INSTANCES` raises `NotImplementedError`. **Stage 4.5A: the inline `TRANSIENTQ` transient extraction-flow time series is supported** (load → write → reload + from-scratch authoring via `transientq_*`; overrides `Q`, not recipients); `TRANSIENTQ`+`NPQRT>0` and external-unit `TRANSIENTQ` data fail explicitly. Registered in `MfUsg.load()`. See `USGT_STAGE4_04_PARAMETERS_QRT.md` / `USGT_STAGE4_05_QRT_TRANSIENTQ.md`. **Verified** |

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
| DPF | ✅ (intentionally not Full) | Stage 5.01: material arrays layer by layer per `GWF2DPFU1AR`, effective LAYCON, `TABRICH`/`BUBBLEPT` from BCF or LPF, default unit; Python regressions + confined executable check. D7 stages 1–2: double-precision `HNEWIM` with `DPIN`/`DPIO` and binary `HNEWIM` output (technically closed in Gate C, 2026-10-07; run with USG-T 2.8.1 only, 2.7 not run for D7; no full atomicity). Gaps keeping it not-Full: binary rules of the other DPF arrays; the CHD/immobile budget behaviour is characterized only for a small LPF/Richards case. |
| PCB | ✅ Full | `(node, iSpec, conc)` authoring + round-trip; dataset order verified vs Fortran; exercised via Ex transport models. |
| OC | ✅ (intentionally not Full) | OC Fullness Card A (Stage 4): broad authoring/round-trip now tested (BOOTSTRAPPING header on line 1, per-SP bootstrap toggles, layer-qualified actions, DDREFERENCE, CONC/BUDGET, ATSA). **Stage 4.7B: `ATSA` exe-verified.** Explicit gaps keep it not-Full: `SAVE IBOUND` solver-rejected by USG-T 2.7; FASTFORWARD/FASTFORWARDC + BOOTSTRAPPING execution = manual tier (external prior-run file); numeric→words rewrite (G05: ATS cases converted or rejected, Gap §13). See `USGT_STAGE4_OC_FULLNESS.md` / `USGT_STAGE4_18_EXE_VERIFICATION.md`. |
| EVT | ✅ (intentionally not Full) | EVT Fullness Card B (Stage 4): all NEVTOP modes (struct + unstruct), transport IETFACTOR 0/<0/>0, per-SP reuse authored + round-trip tested; IETFACTOR round-trip bug fixed; exe smoke added. **Stage 4.6F-A: NPEVT EVTR array parameters preserved + authored from scratch** (`parameters=`/`evtr_parm=`, PARAMETER line, INSTANCES; `expand_parameters=True` keeps legacy NPEVT=0). **Stage 4.6F-B: ETS zonal time-series deferred** — an ATS-coupled dynamic mode (external time-series superseding EVTR), fails explicitly per branch (authoring + load) with the full spec recorded. That one honest gap keeps EVT not-Full. See `USGT_STAGE4_EVT_FULLNESS.md` / `USGT_STAGE4_16_EVT_ETS_ZONAL.md`. |
| MDT | ✅ (intentionally not Full) | MDT Fullness Card C (Stage 4): field order Fortran-audited; from-scratch authoring + round-trip tested for all main branches; fixed FRAHK/FRADARCY load, write_file(f=) crash, debug print; added option/length validation. Ex7 real-model round-trip/run kept. **Stage 4.7B: execution covered by the Ex7 models**; gaps kept manual tier (not-Full): chained-decay NTCOMP>MCOMP (tailored BCT setup) and AI1/AI2 (`SEPARATE_AI2`/`MULTIFILE_MD`) per-component binary outputs. See `USGT_STAGE4_MDT_FULLNESS.md` / `USGT_STAGE4_18_EXE_VERIFICATION.md`. |
| LAK | ✅ (intentionally not Full) | LAK Fullness Card D (Stage 4): from-scratch authoring + round-trip tested (no-transport, classic transport, TRANSPORTBOUNDARY MCOMP>1, TABLEINPUT); six authoring bugs fixed; validation added; Ex8 real-model round-trip/run kept. **Stage 4.6D: multi-lake + sill/connectivity (datasets 7/8) authored from scratch + validated** (`sill_data` 1-based lakes, `IC>=2`, `ITMP>0` rule; 2 round-trip tests + 7 negatives). Gaps keeping it not-Full: `TABLEINPUT` table contents external; GAGE separate; multi-lake-with-sill from-scratch *execution* not exe-smoke-tested (manual tier). See `USGT_STAGE4_LAK_FULLNESS.md` / `USGT_STAGE4_15_LAK_MULTILAKE_CONNECTIVITY.md`. |
| GNC | ✅ | GNC unstructured ghost-node helper; verified. |
| TIB | ✅ Full (authoring) | Promoted in Stage 4.1: semantic `stress_period_data` authoring + `parse=True` load for flow (`NIB0/NIB1/NIBM1`) and transport (`NICB0/NICB1/NICBM1`) blocks; `load` default stays byte-exact raw round-trip, with raw fallback on unsupported `U1DINT` (`EXTERNAL`/`OPEN-CLOSE`). |
| GSF | ✅ Full (authoring) | Promoted in Stage 4.2: semantic `vertices`+`node_data` authoring + `parse=True` load (0-based internal / 1-based file), `from_grid` (requires per-vertex z), `to_grid()` via `from_gridspec`; `load` default stays raw round-trip. GSF is not solver input. |
| ETS | Parameter-preserving + from-scratch authoring (ETSR array) | Stage 4.4A: ETSR array parameters load → write → reload with syntax intact (`NPETS>0`, defs, per-period activation records, `INSTANCES`). **Stage 4.6C-D: from-scratch ETSR array-parameter authoring** (`parameters={name: {parval, clusters \| instances}}`/`evtr_parm`, auto `NPETS`/`nclu`, partyp validated `ets`, first period must activate; validated before opening, no partial file). Opt-in expanded fallback (`expand_parameters=True` → `NPETS=0`). `npets>0`/`evtr_parm` with no defs raises `ValueError`. See `USGT_STAGE4_04_PARAMETERS_ETS.md` / `USGT_STAGE4_14_ETS_PARAMETER_AUTHORING.md`. |
| HFB | Parameter-preserving + from-scratch authoring (list) | Stage 4.4B: `NPHFB>0` list parameters load → write → reload with syntax intact (defs + `NLST` barrier rows + `NACTHFB` + active names). **Stage 4.6C-A: from-scratch `NPHFB>0` authoring** (`parameters=`/`acthfb_names`, auto counts, normalized + validated; structured + unstructured). `TRANSIENT_HFB`+`NPHFB>0` and `INSTANCES` still fail explicitly. See `USGT_STAGE4_04_PARAMETERS_HFB.md` / `USGT_STAGE4_11_HFB_PARAMETER_AUTHORING.md`. |
| DPT | Partial | **Stage 4.6E:** `A-W_ADSORBIM` immobile air-water adsorption authored/loaded/written for the array-only branches (`IAREA_FNIM ∈ {1,4}` × `IKAWI_FNIM ∈ {1,2}`); scalar/tabular branches (`IAREA_FNIM ∈ {2,3,5}`, `IKAWI_FNIM ∈ {3,4}`) raise a specific `NotImplementedError` before any shifted read. `DLIM` conditional correct. |
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

**Stage 5.02-B1 (LPF, USG-T 2.8.1 `gwf2bcf-lpf-u1.f`).** Per layer the reader
takes HK, HANI (`CHANI<=0`), VKA (DISU: only with `IKCFLAG=0`); SS when
transient; SY when transient and `LAYTYP!=0` after the `THICKSTRT` reset
(3198–3209, 3653/3918); VKCB when `LAYCBD!=0` (DISU: and `IKCFLAG=0`,
3676/3934); WETDRY when `LAYWET!=0` (3698/3950); the Richards arrays when the
effective LAYCON is 5 and `TABRICH` is off (3708/3955), where `RICHARDS` or
`RICHARDS_HP` turns every LAYCON!=0 into 5 (3334–3349); item 7 when any
`LAYWET!=0` (3286–3299); the connection array when `IKCFLAG!=0` (3974).
`MfUsgLpf` now follows these conditions. Not covered: DISU `IDSYMRD=1`
connection-array lengths (and their interaction with CLN/GNC), and LPF named
parameters.

**Stage 5.02-B2 (BCF, USG-T 2.8.1 `gwf2bcf-lpf-u1.f`).** Item 1 is read in
free format or with `FORMAT(I10,F10.3,I10,F10.3,4I10)` when the BAS file has no
`FREE` option (91–113; IKCFLAG only on DISU; header options are not parsed from
a fixed-format line). In an F10.3 field without a decimal point in the mantissa
the last three digits are fractional even when an exponent follows (`-1E+30`
is read as `-1E+27`, `123E+2` as `12.3`); FloPy writes these reals with a
decimal point (`-1.E+30`) and fails before writing if a value cannot keep its
single-precision value within ten characters. LAYCON codes `10*INTERCELLT+LAYCON` are read with
`READ(IN,*)` or `FORMAT(40I2)` (231–237) and `RICHARDS`/`RICHARDS_HP` turns
every LAYCON!=0 into 5 (335–350). TRPY is read with `IKCFLAG=0` (353–361) and
DISU ANGLEX when some `|TRPY-1|>1e-6` (747–749). Per layer: SF1 when transient;
with `IKCFLAG=0`, TRAN for LAYCON 0/2 and HY otherwise, then VCONT
(`IKVFLAG=0`, all but the last layer) or KV (`IKVFLAG!=0`, every layer of a
multi-layer model); SF2 when transient and LAYCON is 2–5; WETDRY when
`IWDFLG!=0` and LAYCON is 1 or 3; the Richards arrays when LAYCON is 5 and
`TABRICH` is off (S 558–687, G 753–816); DISU connection array when
`IKCFLAG!=0` (819). `MfUsgBcf` now follows these conditions. Real
fixed-format models (Ex3, Ex4) rewritten by FloPy now run with the same heads
as the original files. Not covered: DISU `IDSYMRD=1`, option detection when
loading fixed-format files, the free-format header precision of HDRY/WETFCT,
and an engine abort with transient models whose only convertible BCF layers
are LAYCON 5 (seen with hand-written input; cause not isolated).

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

**Stage 5.01 correction (USG-T 2.8.1 `GWF2DPFU1AR`, `gwf2dpf1u1.f` 111–191).**
The material arrays are read inside `DO 200 K=1,NLAY`: for each layer `PHIF`,
`DDFTR`, `SC1IM`, then `SC2IM` if `LAYCON(K)!=0`, then `alphaIM`, `betaIM`,
`srIM`, `brookIM` if `LAYCON(K)=5` and `ITABRICH=0`, plus `bPIM` if `IBPN>0`.
`IUZONTABIM` (one `NODES` record), all `IBOUNDIM` layers and all `HNEWIM`
layers come first, each complete. `MfUsgDpf` previously wrote and loaded the
material arrays property by property, so any multi-layer file was mis-read by
USG-T (confirmed with the executable) and two Richards layers did not even
survive a FloPy round-trip. Old multi-layer DPF files from this fork must be
regenerated; they are not detected or reinterpreted.

The conditions use the effective LAYCON of the flow package: BCF `LAYCON`
(already separated from the interblock code); LPF `LAYTYP` 0→0, 1–3→3, 4→4,
5→5, negative→0 with `THICKSTRT` and 3 without; then BAS `RICHARDS` or
`RICHARDS_HP` turns every LAYCON≠0 into 5 (`gwf2bcf-lpf-u1.f` 343–348,
3166–3179, 3345–3348). `TABRICH` and `BUBBLEPT` are taken from BCF or LPF.

DPF with Richards layers has comparative executable evidence (FloPy-written
and hand-written DPF files produced identical heads in an all-LAYCON-5 LPF
model, with and without `BUBBLEPT`).

**Budget with CHD cells (USG-T 2.8.1).** The 9–13 % flow-budget discrepancy of
those runs comes, in the case studied, from how the budget is printed:
`GWF2DPFU1FM` also applies the mobile-immobile exchange at CHD cells
(`IBOUNDIM!=0` and `IBOUND!=0`, `gwf2dpf1u1.f` 381/393) and `GWF2DPFU1BDS`
reports their IMMOBILE STORAGE, while the CONSTANT HEAD term
(`GWF2BCFU1BDCH`/`GWF2BCFU1BDADJ`, `gwf2bcf-lpf-u1.f`) only sums the flows
between cells. The case is small: LPF `LAYTYP=[5,5]` with BAS `RICHARDS`, a
1x3x2 grid, two CHD cells, recharge, DPF with LAYCON 5 and five 1-day steps.
Let `D = IN - OUT` be the printed rate residual of a step and `q` the signed
IMMOBILE STORAGE of the IDPFCB record summed over the CHD cells (negative
when water goes into immobile storage). In every step of that case `D`
coincides with `q` within 2e-4 L3/T (four printed decimals, single-precision
cell values), so the reconstructed residual `D - q` is of that order; the
printed percentage is `100 D / ((IN + OUT) / 2)`. A cumulative reconstruction
needs the rates integrated over the step lengths. This is a diagnostic
reconstruction: the engine budget is left as printed and FloPy applies no
correction. Whether the omission is intended has not been determined, the
result is not generalized to other DPF models or boundary types, and it has
only been run with USG-T 2.8.1.

Setting `IBOUNDIM=0` under the CHD cells closes the printed budget in that
case and the mobile heads coincide, but it removes the immobile domain at
those cells, i.e. it is a different model; it is only an option where that is
physically acceptable, and it is not shown to be equivalent in general
(particularly for transport). The executable tests
`test_usgt_exe_dpf_chd_*` and `test_usgt_exe_dpf_no_immobile_under_chd_*`
characterize this behaviour for the case above.

**D7 stage 1 — `HNEWIM` with BAS `DPIN`/`DPIO` (USG-T 2.8.1).** With
`IDPIN=1` `GWF2DPFU1AR` reads `HNEWIM` with `U1DREL8` (`gwf2dpf1u1.f` 131–157,
`utl7u1.f` 786–922), otherwise with `U1DREL`. `MfUsgDpf` now stores `hnewim`
as float64 when the model's BAS6 has `DPIN` or `DPIO` at construction (or,
for a direct `MfUsgDpf.load` before BAS6, when the BAS6 options line in
`ext_unit_dict` has them; the BAS6 handle position is restored). Rules:

- Text is written with 17 significant digits (`E24.16`). A uniform value is
  written as `CONSTANT` only if it is exact as REAL, because `CNSTNT` is REAL in
  `U1DREL8`; otherwise as `INTERNAL`. An explicit `how="constant"` with such a
  value is rejected. Multipliers are REAL (`A*CNSTNT`, `CNSTNT=0` means 1),
  also after assigning `cnstnt` to a layer or to the whole array; a loaded
  `CONSTANT` keeps its REAL value.
- Binary data are loaded only through a negative unit (`EXTERNAL -n` or a
  fixed-format negative `LOCAT`; `FMTIN`/`IPRN` optional): 52-byte header
  with REAL*8 `PERTIM`/`TOTIM`, then REAL*8 values per layer, without the
  multiplier (`U1DREL8` returns before applying it). Later data in a shared
  file are left unread; truncated data raise an error. A header written in
  the other precision is not always detectable. `OPEN/CLOSE … (BINARY)` and
  `EXTERNAL n … (BINARY)` with a positive unit are rejected for `HNEWIM`
  (USG-T 2.8.1 tries a formatted read with format `(BINARY)`).
- Before the DPF file is opened (or a caller's handle written), `write_file`
  rejects: `DPIN`/`DPIO` active with single-precision `hnewim` (DPF created
  before BAS6, or the option set later; create `MfUsgDpf` after `MfUsgBas`),
  and the invalid `CONSTANT` above. Clearing `DPIN`/`DPIO` later only warns (USG-T then reads REAL).
  Single-precision `Util2d`/`Util3d` layers are rejected under `DPIN`, also
  when they are assigned after construction (`dpf.hnewim[k] = Util2d(...)`,
  a list of `Util2d`, or `hnewim.util_2ds[k]`): `write_file` checks every
  layer, reports all single-precision layers and writes nothing. A float32
  ndarray passed as values is promoted exactly, which does not restore
  digits already rounded. Without `DPIN`/`DPIO` a float64 layer only warns.
- Executable tests (`DDFTR=0`, REAL*8 `IMMOBILE HEAD` output) show the values
  read for FloPy text, a hand-written binary record, fixed-format records and
  an irregular DISU grid. A custom `fmtin` with fewer digits is honored and
  rounds the values written.

**D7 stage 2 — binary `HNEWIM` output.** Layers with `format.binary = True`
are written, in layer order, to one file of the package (`<model>.hnewim.bin`,
under `external_path` if set) registered as `DATA(BINARY) n`, and read with
`EXTERNAL -n 1.0 (BINARY) -1` (fixed format: negative `LOCAT`). Per layer: 52
bytes of header with REAL*8 times and REAL*8 values under `DPIN`/`DPIO`, 44
bytes and REAL values otherwise; `NSTRT`/`NNDLAY` are the global node range
of the layer and `ILAY` the layer. Text and binary layers can be mixed;
irregular DISU grids are supported. The values written are the effective
ones (`CNSTNT` applied once) because the 1-D binary branch does not apply
the multiplier; the record multiplier is 1.0.

- Precision: REAL*8 records need float64 layers and REAL records float32
  layers; a layer of the other precision is rejected (FloPy does not round
  float64 values to write REAL).
- Destinations: before anything is written, the package compares the
  normalized destination (absolute, symbolic links resolved, case folded on
  case-insensitive file systems) with the name file, package files, output
  files and external files of the model, and with existing files. A taken
  name, or a file of another component (registered or not), is never
  overwritten and an entry of another component is never replaced: the next
  free `<model>.hnewim_<i>.bin` is used. A new file is created exclusively.
  Units avoid package, output and external units.
- Ownership: the instance remembers the identity (device and inode) of the
  file it created and the name-file entry it created, not only the path. A
  later write rewrites a destination only if it is still that regular,
  singly linked file and no destination of another component refers to it
  (the final open does not follow links and checks the identity again
  before truncating); a path replaced by a symbolic link or another file, a
  directory link that now resolves elsewhere, or a registered hard link
  makes the package use a free name instead. A normal rewrite keeps its
  file and unit; after a
  change of `external_path`, model name or workspace the new destination is
  used, the own stale entry is dropped and the previous file is kept. If
  another component takes the unit, its entry is kept and a new unit is
  used. Disabling binary output drops only the own entry; files are never
  deleted. A loaded package does not adopt the original file.
- `write_file()` registers the file in the model but does not rewrite the
  name file (`write_input()` does). Checks (precision per layer, number of
  values, text-layer records, folder, destination and unit) run before any
  file is written; there is no full atomicity: an I/O error after the binary
  file is written can leave it on disk and registered, and changes made
  concurrently by other processes are not fully excluded.

Still open: the other DPF arrays (`U1DREL`) keep the generic binary handling; BAS `STRT` under `DPIN` (`U2DREL8`/`U1DREL8`) is still
single precision; global float32 output formats are unchanged; USG-T 2.7 was
not executed for D7.

### Gap §6 — DPT: air-water interface adsorption immobile domain (Stage 3 Card 3)

**Decision (Stage 4.6E): array-only branches implemented; scalar/tabular
branches kept as explicit `NotImplementedError`.**

Fortran-derived spec (`gwt2dptu1.f` option block + `dpt2aw_adsorb.f`
`AW_ADSORBIM1AL/RP1/RP2`): the DPT option line carries
`A-W_ADSORBIM IAREA_FNIM IKAWI_FNIM`, then a cascade of conditional reads:

- If `IAREA_FNIM==5` or `IKAWI_FNIM==4` (tabular): `NAZONESIM NATABROWSIM`, a
  per-node zone map `IAWIZONMAPIM` (`U1DINT`, NODES), and tabular area / K_AWI
  arrays. **Deferred** — `NotImplementedError`.
- `IAREA_FNIM==1`: `AWAMAXIM(NODES)`; `==4`: `AWAREA_X2IM/X1IM/X0IM(NODES)`
  (**implemented**, RP1, after heat); `==2` grain diameter / `==3` `ROG_SIGMA`
  constant (**deferred**).
- Langmuir isotherm arrays `ALANGAWIM(NODES,MCOMP)`, `BLANGAWIM(NODES,MCOMP)`
  (RP2, per species before `ADSORBIM`) for `IKAWI_FNIM ∈ {1,2}` (**implemented**;
  `==1`/`==2` share the file layout). `==3` (Brusseau, `SIGMA_RT` scalar) and
  `==4` (tabular K_AWI) are **deferred**.

The implemented subset (`IAREA_FNIM ∈ {1,4}` × `IKAWI_FNIM ∈ {1,2}`) reads only
plain arrays — no zone map, scalar, or table — so authoring/load/write is exact.
The deferred branches raise a specific `NotImplementedError` naming the index
(authoring **and** load, the latter before any array is read, so subsequent
reads never shift). The DLIM conditional read correctly requires both `IDPF/=0`
and `IDISPIM/=0`, and is independent of this option. See
`USGT_STAGE4_06_DPT_AW_ADSORBIM.md`.

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
total per period. From-scratch `NPDRT>0` authoring is supported (Stage 4.6B;
`parameters=`/`active_params=`, auto counts; no-defs → `ValueError`); parameter
`INSTANCES` still raises `NotImplementedError`. Activated SPREAD (`NR<0`)
recipients are
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

From-scratch ETSR parameter *authoring* is now supported (Stage 4.6C-D): pass
`parameters={name: {parval, clusters | instances}}` + `evtr_parm`; see
`USGT_STAGE4_14_ETS_PARAMETER_AUTHORING.md`. Remaining (keeps ETS not `Full`):
only the ETSR array can be parameterized (the Fortran limit), parametric
*execution* is not USG-T smoke-tested, and `npets>0`/`evtr_parm` with no
definitions now raises `ValueError`. The legacy expand path is still available
opt-in via `MfUsgEts.load(..., expand_parameters=True)` (writes `NPETS=0`). See
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

From-scratch `NPHFB>0` parameter authoring is supported (Stage 4.6C-A;
`parameters=`/`acthfb_names`, auto counts; no-defs → `ValueError`). Remaining
(keeps HFB not `Full`): `TRANSIENT_HFB` combined with `NPHFB>0` (the Fortran
would redefine parameters each stress period under `ITERP=1`) and parameter
`INSTANCES` (the Fortran aborts) still raise `NotImplementedError`. See
`USGT_STAGE4_04_PARAMETERS_HFB.md`.

### Gap §10 — EVS/RTS name-file entries: kept as opaque files (G03)

USG-T opens `EVS` (IUNIT 6) and `RTS` (IUNIT 9) as formatted, read-only input
(`SGWF2BAS8OPEN`). EVT with `ETS` and RCH with `RTS` read them as time-series
rows (`TSTART TEND FACTOR value(1..nzones)`; `gwf2evt8u1.f`, `gwf2rch8u1.f`,
`ATS1AJ` in `glo2basu1.f`), which needs ATS. FloPy has no class for these
types and used to drop the entries from the name file on rewrite; USG-T then
stops when RCH/EVT read the series (unit 0).

`MfUsg.load` now keeps the `EVS`/`RTS` entries with any `load_only`
(`preserved_nam_entries`, read-only copies) and `write_input`/`write_name_file`
write them back once:

- Parsing follows URWORD: a name in single quotes can hold spaces and commas
  and is written back quoted. An unclosed quote, an empty name, a missing
  name or a line longer than 299 characters is rejected on load.
- Options `FLUSH` and `NO-ARCHIVE` are kept; other options are dropped with a
  warning (`SGWF2BAS8OPEN` does not use them for an input file).
- A unit <= 0 or a missing file warns on load and is rejected on write, in
  every mode, before any file is written.
- An EVS/RTS unit shared with another entry of the original name file
  (another EVS/RTS, a package, DATA or an output) is rejected on load, with
  the unit, types and line numbers, before the unit-keyed name-file reader
  would keep only one of them. Shared units between other types keep the
  previous behaviour.
- Same workspace: the path is written as it was; the files are not touched.
- New workspace: a relative path inside the original workspace is copied to
  the same relative path, with exclusive creation; an identical destination
  is left as is (bytes and mtime); a different one raises before any package
  or the name file is written.
- An absolute path, or a relative path outside the original workspace, is
  written as an absolute reference and not copied (warning): the model
  still depends on that file.
- Units and files are checked against what the name file lists (LIST, the
  name file, package entries, entries skipped by `load_only`, DATA and
  outputs as written, i.e. outputs by basename) and between EVS and RTS,
  with resolved paths (case-insensitive on macOS/Windows), from
  `write_input` and from a direct `write_name_file`, before writing.

Tested with hand-written fixtures (3-node DISU, ATS, RCH `RTS` 2 zones, EVT
`ETS` 1 zone): the control matches the analytical recharge/ET rates, step
ends and cumulative volumes with USG-T 2.8.1 and 2.7.0; removing the series
changes the results; dropping the entries stops the run; FloPy rewrites
(same and new workspace, RTS only, EVS only, both, quoted names) give the
same budgets and heads as the control.

Limits (open):

- The contents are not interpreted: no RTS/ETS grammar in RCH/EVT.
- Full load of EVT with `ETS` raises `NotImplementedError`. Writing RCH with
  `RTS` used to fail (G01/F39) and is fixed in G01; one-zone RTS rewrites
  were then checked with USG-T 2.8.1 and 2.7.0, and multi-zone rewrites still
  write `INRCHZONES 1` (F42). The G03 rewrites were tested with RCH/EVT
  skipped by `load_only`.
- Entries skipped by `load_only` (e.g. RCH/EVT) are written by basename and
  are not copied to a new workspace (existing behaviour, separate item); the
  tests copy them when running in a new workspace.
- Not a transaction: a later error in another package or an I/O error can
  leave files already written or copied.
- Other types the name-file reader accepts are unchanged; `UZF` (active in
  USG-T, no FloPy class) is still dropped on rewrite (open finding).

### Gap §11 — DIS/DISU: TRTOSS stress periods (G06)

USG-T reads each stress-period record as `PERLEN NSTP TSMULT SS|TR|TRTOSS
[stofrac]` (`SGLO2BAS8ARDIS`, `glo2basu1.f`). `TRTOSS` runs the period as
transient (ISSFLG=0) and, with adaptive time stepping, ends it one step after
the storage rate falls below `stofrac` times the total budget rate
(`GWF2BAS7OT`, `ATS1AD`). FloPy used to load any flag other than `TR` as
`SS`, so `TRTOSS` and its fraction were lost on rewrite.

`MfUsgDis` and `MfUsgDisU` keep `steady` (bool per period) and add `trtoss`
(last argument; storage fraction per period, None/NaN = not TRTOSS, scalar or
`nper` values):

- `steady=True` -> `SS`; `steady=False` -> `TR`; `steady=False` + a
  fraction -> `TRTOSS`. A fraction on a `steady=True` period (also the
  default `steady=True`) is rejected.
- Load accepts `SS`/`TR`/`TRTOSS` in any case and the number syntax of URWORD
  (E/D/Q or no exponent letter, leading dot, sign, comma separator, word of
  at most 30 characters). `TRTOSS` without a fraction loads as 0.0 with a
  warning (USG-T reads 0.0) and is written as 0.0. Other flags, NaN/Inf
  words and values beyond the REAL range are rejected.
- The fraction is written as the shortest text that USG-T reads as the same
  single-precision value (`numpy.float32` of the value; a loaded value keeps
  the REAL value read from the file). Checked against a gfortran conversion
  with the URWORD rules (7014 words, incl. words next to a rounding tie, and
  5006 floats) and with the engine at a termination threshold (two adjacent
  REAL fractions ending at different steps keep their side on rewrite).
- `SS` and `TRTOSS` periods without any `TR` period are rejected on
  construction, load (at DIS/DISU, before LPF/BCF are read) and write:
  USG-T 2.8.1 then sets ITRSS=0 (storage arrays not read, `SC1(1)` used in
  the TRTOSS period), while USG-T 2.7.0 counts TRTOSS as TR and reads them.
- `steady`/`trtoss` are checked again before the file is opened.
- `period_type` and `itrss` are read-only views; `itrss` follows the USG-T
  2.8.1 rule (TRTOSS counts as neither SS nor TR), not 2.7.0.
- The former `Trtoss` argument of `MfUsgDis` has no effect and warns when it
  is set; use `trtoss`.

Tested (hand-written fixtures, DIS/DISU, LPF/BCF): termination one step
after the criterion with the USG-T message, `1e-6` ending later than `1e-3`,
TRTOSS without ATS = TR, a short TRTOSS that does not reach the criterion =
TR and different from SS, FloPy load + rewrite = hand-written control
(steps, times, budgets, heads) with USG-T 2.8.1 and 2.7.0; 2-layer LPF and
BCF with every period TRTOSS keep their storage arrays in place.

Limits (open):

- Early termination depends on ATS in OC; it is documented, not checked.
- FloPy times (`get_totim`) still use PERLEN; the actual end of a TRTOSS
  period is only known from the run.
- ITRSS differs between USG-T 2.7.0 and 2.8.1 when every period is TRTOSS
  (TRANSIENT vs COMBINED); storage is read in both.
- LAK with every period TRTOSS: USG-T 2.8.1 reads NSSITR/SSCNCR (ITRSS=-1)
  while `MfUsgLak` writes them only when a period is steady; not addressed
  (separate gate). DPF, CLN storage and transport with TRTOSS are not
  validated.

### Gap §12 — SMS: options, NUMTRACK records and REAL values (G04)

`SMS7U1AR` (`glo2sms-u1.f`, USG-T 2.8.1) reads the option words of record 1b
with URWORD (any case), the item-2 values BTOL, BREDUC and RES_LIM only when
NUMTRACK > 0, and every REAL through a single-precision REAL. `MfUsgSms`
used to turn `SHIFT` into `DAMPBOT`, drop lower-case options when there was
no preset, fail on NUMTRACK=0 records and on D exponents, and write REALs
with three significant digits.

`MfUsgSms` (constructor unchanged):

- keeps SOLVEACTIVE, DAMPBOT, SHIFT, TRUNCATEDNEWTON and
  TRUNCATEDNEWTONCUTOFF value in their order, in any case, with or without
  a SIMPLE/MODERATE/COMPLEX preset; other words of record 1b are ignored
  like in USG-T (labels);
- reads and writes BTOL, BREDUC and RES_LIM only when NUMTRACK > 0;
- reads REALs with E, D or Q exponents; writes 9.2e when that keeps the REAL
  value (existing files unchanged), otherwise the shortest text USG-T reads
  as the same REAL (checked with a gfortran conversion following URWORD on
  words next to rounding ties, incl. TRUNCATEDNEWTONCUTOFF);
- rejects missing, non-finite or out-of-range REALs before the file is
  opened, and names the SMS field;
- `STRESSPERIODSETTINGS` raises NotImplementedError on load and write (its
  PERIOD blocks are not kept); `SOLVEROUTPUTFILE n` warns on load and is not
  kept (the separate solver output is not written).
- load and write check the option words the same way: split on spaces,
  tabs and commas like URWORD, every TRUNCATEDNEWTONCUTOFF value is checked
  (repeated cutoffs are kept in order; a cutoff with no word after it is
  read as 0.0, like USG-T), so `options2` changed through the API is
  checked before the file is opened.

Tested with hand-written models (DIS and DISU, LPF and a confined BCF,
Picard and Newton, XMD and PCGU, presets incl. SIMPLE + PCGU, NUMTRACK 0 and
20): a FloPy rewrite gives the same engine echo, option messages, outer
iterations, times and heads as the hand-written control, with USG-T 2.8.1
and 2.7.0. The SHIFT cases with PCGU reproduce an abort or the same dry cells
(behaviour kept, not a usable model).

The REAL helpers are shared with the DIS stress-period records
(`flopy/mfusg/mfusgdis.py`, Gap §11).

---

### Gap §13 — OC: ATS, NPTIMES/TIMOT and numeric ATS records (G05)

USG-T 2.8.1 (`glo2basu1.f`) turns ATS on with `ATSA`, `NPTIMES` or `NPSTPS`
on the first OC record (URWORD stops scanning options once a word ends at
column 198), reads the NPTIMES output times TIMOT right after that record
as DOUBLE values with a list-directed READ (any number of lines, blanks,
tabs or commas) and stops when they do not increase by more than
TMINAT=1.0E-10. With ATS, output control is read once per stress period:
the time step of a PERIOD record is skipped and DDREFERENCE is read only
right after the period number. Numeric OC reads, per period, the ATS1RP
record DELTAT TMINAT TMAXAT TADJAT TCUTAT NPTIMES NPSTPS IBOOT IBOOTSCALE;
its NPTIMES and NPSTPS replace the first-line values, values below the
ATS1RP limits get defaults, DELTAT is limited by TIMOT and set to PERLEN in
steady-state periods. The word form (SGWF2BAS7N) does none of this.
`MfUsgOc` used to write numeric ATS records as literal `{float(...)}`
text, read TIMOT as integers (decimals lost, then a write error), write
TIMOT twice (pushing NPSTPS past column 200, where USG-T ignores it) and
round ATS parameters to four digits.

`MfUsgOc` (constructor unchanged):

- reads TIMOT like USG-T (decimals, E/D/Q exponents, several lines, commas)
  and writes it once after the first line with its DOUBLE value; repeat
  counts, null values, `/`, NaN/Inf, values beyond NPTIMES and OC records
  where times are expected are rejected; options USG-T does not scan on the
  first record are ignored on load (files written by older FloPy);
- turns ATS on with nptimes or npsteps (load and write); a words OC file
  whose first word is NPTIMES, NPSTPS or ATS is rejected (USG-T stops);
- keeps the text of DELTAT, TMINAT, TMAXAT, TADJAT and TCUTAT once it is a
  REAL USG-T reads; name and value are separated like URWORD (blanks, tabs
  or commas, any case) on load and in the check before writing (HCLOSE and
  BTOL are still written with four digits);
- with ATS writes `PERIOD k [DDREFERENCE]` without STEP, and rejects on load
  and write a stress period given twice or out of order and DDREFERENCE
  after STEP; without ATS the records are unchanged;
- checks nptimes, npsteps, TIMOT (count, one-dimensional, finite, gap) and
  the ATS blocks before the file is opened;
- converts numeric ATS records of free-format flow models to words with the
  values USG-T uses (ATS1RP defaults, DELTAT at or above PERLEN in
  steady-state periods, cut to PERLEN by USG-T); an NPSTPS common to all
  records becomes npsteps, 0 meaning no output cadence (a numeric file with
  NPSTPS 0 makes USG-T evaluate MOD(KSTP,0), whose result depends on the
  platform: the conversion is compared with a words control, not with that
  file). Rejected (NotImplementedError): per-period NPTIMES different from
  the first line (also records with five values), NPSTPS that changes,
  IBOOT/IBOOTSCALE other than 0, INCODE < 0, ICBCFL = 2, IHDDFL = 0 with
  head or drawdown flags, DELTAT limited by TIMOT or TIMOT read beyond
  NPTIMES+1, fixed format, transport, CLN, DPF, FASTFORWARD(C) and
  BOOTSTRAPPING.

Tested with hand-written models (DIS and DISU, words and numeric OC, TIMOT
integer/decimal/exponent/multi-line/60 values, NPSTPS, ATS parameters with
many digits, ATS1RP limits, steady-state PERLEN 5.0 and 0.1, three periods,
per-layer flags, controls without ATS): the FloPy rewrite of the OC file
gives the same heads/budget output times, time steps, budgets, convergence
and heads as the control, with USG-T 2.8.1 and 2.7.0. Not covered: numeric
OC with transport, CLN, DPF, FASTFORWARD or BOOTSTRAPPING; IHDDFL = 0 in
numeric OC without ATS (unchanged); fixed-format numeric ATS records.

---

### Gap §14 — Voronoi grids to DISU (VORONOI)

`flopy.utils.voronoi.get_disu5_gridprops` and
`VoronoiGrid.get_disu5_gridprops` build DISU arrays (IVSD=0, IDSYMRD=0) for
one region without holes, extruded into contiguous layers with the same
plan view, for LPF flow. Source audited: USG-T 2.8.1. Run with USG-T 2.8.1
and 2.7.0 (linear flow, vertical series, two layers against an independent
solution) and on the Voronoi USG example (7,062 nodes); not run with 2.6.2.

Limits (open):

- Holes and separate regions are rejected; a `VoronoiGrid` built with holes
  does not cover its domain (`tri2vor`, separate item).
- Pinch-outs, different plan views per layer, LAYCBD, IVSD other than 0,
  IDSYMRD=1, CLN, GNC, anisotropy angles (CHANI other than 1) and DISU6 are
  not supported. Transport (BCT dispersion with vertical CL12) is not
  validated.
- Cell centers must be the Voronoi generators; gaps along the outer boundary
  are not detected.
- `MfUsgDisU` converts `ja` in place; `top`/`bot` are stored as float32; the
  `Gridgen.get_fahl` docstring is not corrected.

---

## Previously-missing packages — now implemented

### SGB — Specified Gradient Boundary (`glo2sgbu1.f`) — IMPLEMENTED

`MfUsgSgb` applies a specified hydraulic gradient at boundary nodes (distinct
from GHB, which specifies a head-to-head conductance). Node-based
`(node, gradient)` list with AUX transport concentrations and `ITMP/-1`
reuse; internal nodes 0-based, file 1-based. Registered as `"sgb"` in
`MfUsg.load()`, so SGB models are no longer silently skipped. `NPSGB>0`
parameter *definitions* are preserved and authorable from scratch (Stage
4.6C-B); active SGB parameters stay unsupported (the USG-T 2.7 `PARTYP='SGB'`
vs `'G'` Fortran type conflict aborts the run), so per-SP `NP>0`/`INSTANCES`
fail explicitly.

### QRT — Sink with Return Flow (`gwf2QRT8u.f`) — IMPLEMENTED

`MfUsgQrt` extracts water at sink nodes and returns a proportion to one or
more recipient nodes (analogous to DRT). Per-sink `(node, q, rfprop)` with
variable-length recipient-node lists (`NodQRT`, read/written via `U1DINT`),
`CHANGEC`/`IQCHNGTYP` return-flow concentration-change type, AUX, and
`ITMP/-1` reuse. `AUTOFLOWREDUCE` is preserved as an option. `NPQRT>0` list
parameters are structurally preserved (Stage 4.4E). **Stage 4.5A: the inline
`TRANSIENTQ` transient extraction-flow time series is now supported** — loaded,
written, and authorable from scratch via the `transientq_*` attributes (it
overrides `QRTF(4)=Q` each step; recipients are untouched). **Stage 4.6C-C:
from-scratch `NPQRT>0` parameter authoring is supported** (structural-only — an
activated parameter is not execution-guaranteed, see the QRT row above).
`TRANSIENTQ` combined with `NPQRT>0`, external-unit `TRANSIENTQ` data, and
parameter `INSTANCES` still fail explicitly. Registered as `"qrt"` in
`MfUsg.load()`. See `USGT_STAGE4_05_QRT_TRANSIENTQ.md`.

---

## Packages intentionally not implemented

| Package | Reason |
|---------|--------|
| STR, SFR | CLN is the preferred surface-water coupling for USG-T models in this project |
| MNW1, MNW2 | Replaced by CLN in USG-T applications |
| OBS, HOB, DROB, RVOB, GBOB, CHOB | Observation packages — low priority |
| UZF | Unsaturated zone flow. Potentially relevant but out of scope. Note: USG-T 2.8.1 reads it (IUNIT 55) and a `UZF` name-file entry is dropped on rewrite (open finding, not addressed) |
| IHM / SYF | Runtime coupling to external HSPF. No standalone input package needed |
| GMG, PCG, SIP, DE4 | Solvers not used with USG-T; SMS and XMD are the relevant solvers |

---

## Separate findings (not addressed)

- `MfUsgWel` on DISU writes `0` (no wells) for stress periods not given in
  `stress_period_data` instead of `-1` (reuse), unlike the structured WEL
  (found while building the G06 fixtures).
- `perlen` is stored as float32 in `MfUsgDis`/`MfUsgDisU`, while USG-T reads
  PERLEN in double precision.
- Engine: with SMS `SHIFT` the shift is also applied to constant-head cells
  (`INVSHIFTVECTOR` skips only IBOUND=0); in small PCGU models with CHD this
  ends in "CONSTANT-HEAD CELL WENT DRY" or in dry cells (USG-T 2.8.1 and
  2.7.0). FloPy keeps SHIFT as written.

---

## Future note (outside G03 and the current estimates)

- mod-PATH3DU support, after USG-T 2.8.1 conformance.

---

## Validation record

| Model | Grid | Version | Packages | Result |
|-------|------|---------|----------|--------|
| MD26 BC2 (real-world, ~112k nodes, 493 SPs) | DISU | USG-T 1.8 | BAS6, DISU, LPF, WEL, GHB, DRN, CHD, RIV, RCH, EVT, CLN, BCT, TIB, OC, SMS | Bit-for-bit identical on all budget columns (flow + transport) |
| Sóncor (real-world, ~19k nodes, 1382 SPs, 4322 TSs) | DISU | USG-T 2.7 ARM | BAS6, DISU, LPF, WEL, GHB, DRN, CHD, RIV, RCH, EVT, CLN, BCT (IDISP=2), DDF, OC, SMS | Bit-for-bit identical on LST, HDS, CON, CBB (7 record types) |
