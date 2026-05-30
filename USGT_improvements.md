# USG-Transport improvements on this fork

This fork of [modflowpy/flopy](https://github.com/modflowpy/flopy) carries
USG-Transport additions on `develop`, aimed at improving how
`MfUsg.load(...).write_input()` round-trips existing USG-T models and
extending post-processing utilities for transport output. Nothing else is
changed.

Tested against **USG-T 1.8** and **USG-T 2.7** binaries on both structured and
unstructured grids.

These are incremental improvements on top of the USG-T work already in
upstream flopy. They live here while testing continues.

## What's added

### Round-trip fixes (original five)

| Branch | What it does |
|---|---|
| `usgt/rch-transport-fix` | In `MfUsgRch.write_file`, `"# Stress period {kper + 1}"` becomes an f-string so the comment substitutes per SP, and `INRECH` is followed by `INIRCH` only when `NRCHOP == 2` (matching the read path). Both details were harmless at runtime but left the written RCH file not round-trippable. |
| `usgt/cln-load-none-unit` | In `MfUsgCln.__init__`, treat `None` entries in `unitnumber` as 0 instead of calling `int(None)`. This comes up when a CLN-declared output unit is not declared in the NAM's `ext_unit_dict`. |
| `usgt/tib-package` | Adds a new `MfUsgTib` class for the Transient Ibound package. Raw-body text-preserving round-trip — enough to load + write an existing TIB without mis-parsing `U1DINT` continuation lines. |
| `usgt/bas-preserve-unstructured` | `MfUsgBas.write_file` re-emits the `UNSTRUCTURED` keyword when `parent.structured is False`. The load side already reads the token; adding it to write closes the round-trip. |
| `usgt/nam-rebase-output-paths` | `BaseModel._reset_external` stores the basename of output files on `change_model_ws`, and `Modflow.write_name_file` preserves subdirectories for external input files. Previously the NAM writer also reduced inputs to basenames, breaking valid `DATA` paths in subfolders. |

### New packages and utilities

| What | File | What it does |
|---|---|---|
| `MfUsgChd` | `flopy/mfusg/mfusgchd.py` | CHD package for unstructured USG-T grids. Node-based (replaces k/i/j), supports AUX concentration variables. Internal `node` values are 0-based and file I/O is 1-based. Full `load` and `write_file` for the USG-T format (`NACT    Stress Period N` headers, `-1` reuse). |
| `MfUsgRiv` | `flopy/mfusg/mfusgriv.py` | RIV package for unstructured USG-T grids. Internal `node` values are 0-based and file I/O is 1-based. Supports AUX concentration and a trailing reach-ID column (`irch`) that is written positionally without being declared as AUX. `irch` is auto-detected from the first data row when loading. |
| `MfUsgEts` | `flopy/mfusg/mfusgets.py` | Segmented Evapotranspiration (ETS) package for USG-T 2.7. Supports `NETSEG > 1` with per-SP `PXDP`/`PETM` segment arrays and the `IESFACTOR` transport flag. Parameterized ETS files load by expanding parameters to concrete arrays and then write as valid non-parametric `NPETS=0`; preserving parameter syntax is not yet supported. |
| `MfUsgGhb` | `flopy/mfusg/mfusgghb.py` | GHB package for unstructured USG-T grids. Internal `node` values are 0-based and file I/O is 1-based. Supports AUX concentration variables, stores `ipakcb` (CBC unit). Full `load` and `write_file`. Load registry now maps `"ghb"` to `MfUsgGhb`. |
| `MfUsgDrn` | `flopy/mfusg/mfusgdrn.py` | DRN package for unstructured USG-T grids. Same pattern as GHB: internal `node` values are 0-based, file I/O is 1-based, AUX support, `ipakcb`. Full `load` and `write_file`. Load registry maps `"drn"` to `MfUsgDrn`. |
| `MfUsgTvm` | `flopy/mfusg/mfusgtvm.py` | TVM2 (Time-Variant Materials) package. Semantic implementation: global interpolation controls plus nper+1 stress-period boundary records, with 0-based internal nodes and 1-based file I/O. Missing boundaries emit all-zero headers on write. |
| `MfUsgGsf` | `flopy/mfusg/mfusgsf.py` | Grid Specification File wrapper. Text round-trip (stores raw lines). `to_grid()` delegates to `UnstructuredGrid.from_gridspec()` for full geometric parsing. Load registry maps `"gsf"` to `MfUsgGsf`. |
| `MfusgTransportListBudget` | `flopy/utils/mflistfile.py` | Reads transport species budget from a USG-T listing file for a single species. Handles both **old** USG-T format (transport blocks use `VOLUMETRIC BUDGET`, same keyword as flow) and **new** format (transport blocks use `MASS BUDGET`). Instantiate once per species: `MfusgTransportListBudget("model.lst", species=2)`. Returns the same recarrays / DataFrames as `MfusgListBudget`. |

Authoring note: semantic packages are expected to support direct construction
from Python/numpy inputs, not only `load()` + `write_file()` round-trips. The
USG-T boundary package tests now include from-scratch CHD/RIV/GHB/DRN creation
with 0-based internal nodes and 1-based file output.

### Authoring-focused fixes (2026-05-29)

| File | Fix |
|---|---|
| `flopy/mfusg/mfusgbas.py` | `MfUsgBas.load` now preserves the `UNSTRUCTURED` option on the constructed package. Programmatic tests cover write/load of `PRINTFV`, `CONVERGE`, `UNSTRUCTURED`, `FREE`, `PRINTTIME`, `SHOWPROGRESS`, `RICHARDS`, `DPIN`, `DPOUT`, `DPIO`, `SY-ALL`, and `STOPERROR`. |
| `flopy/mfusg/mfusgwel.py` | `options=None` is normalized before AUX auto-registration, so programmatic WEL authoring with AUX fields no longer crashes. Tests cover WEL rates assigned to GWF and CLN nodes, `ITMP NP ITMPCLN` headers, AUX concentrations, and 0-based internal node storage after reload. CLN connectivity stays in the CLN package. |
| `flopy/mfusg/mfusgcln.py`, `flopy/mfusg/cln_dtypes.py` | CLN now supports direct authoring and reload of `PROCESSCCF`/`ICLNGWCB`, `GENERAL_SEC`, and the 9-field `ISHAPE` node-property format used by rectangular/general conduit shapes. |
| `flopy/mfusg/mfusgdpf.py` | DPF now sets `model.idpf=1` when constructed programmatically, writes/loads optional `FRAHK`, supports `IUZONTABIM` for TABRICH models, writes `SC2IM` only for convertible layers, and supports immobile Richards arrays (`alphaIM`, `betaIM`, `srIM`, `brookIM`, optional `bPIM`). |

### Density-coupled round-trip fixes (2026-05-19)

Four bugs found during round-trip testing against a real-world Vistas-generated
density-coupled USG-T 2.7 model (BCT IDISP=2, DDF, 1382 stress periods):

| File | Fix |
|---|---|
| `flopy/mfusg/mfusgddf.py` | `ithickav` default: `default_val=1` → `default_val=0`. USG-T treats the absent ITHICKAV field as 0 (arithmetic averaging). Loading it as 1 changed transmissivity in the density layer and caused transport divergence. |
| `flopy/mfusg/mfusgchd.py` | `shead`/`ehead` dtype: `np.float32` → `np.float64`. Float32 round-trip loss (~0.04 mm per CHD node) accumulated over 469 density-coupled SPs into concentration differences that marginally failed outer-loop convergence (ΔC > CICLOSE=1e-8). |
| `flopy/mfusg/mfusgriv.py` | `stage` dtype: `np.float32` → `np.float64`, same reason. |
| `flopy/mfusg/mfusgbas.py` | `CONVERGE` option: (1) `converge=converge` was missing from the `cls(...)` call in `load`, so `self.converge` was always `False`; (2) `write_file` never emitted `CONVERGE` even when `self.converge=True`. The `CONVERGE` keyword in BAS6 tells USG-T to use the coupled flow–transport convergence criterion in the outer nonlinear loop — without it the model failed to accept time step 2 of the first pumping stress period after 250 iterations. |

### Minor fixes

| File | Fix |
|---|---|
| `flopy/mfusg/mfusgbct.py` | Removed a stray `print()` left from development. Fixed file-handle management in `write_file`: the file is now closed only when opened internally (`close_on_exit` flag), so callers that pass an open handle are not surprised. `ICBUND` now loads as `np.int32`, so round-trip writes integer constants instead of `1.000000E+00` values that USG-T rejects. |

### Fortran source audit, new implementations, and bug fixes (2026-05-20)

All packages cross-referenced against the USG-Transport 2.7.0 Fortran source code
(`USGT_V_2-7-0_Source_Code/`). See `USGT_roadmap.md` for the full coverage table.

#### New: `MfUsgHfb` — HFB with TRANSIENT_HFB support

| File | What it does |
|---|---|
| `flopy/mfusg/mfusghfb.py` | HFB6 (Hydraulic Flow Barrier) implementation for unstructured USG-T grids. Extends `ModflowHfb` with: (a) node-based `(node1, node2, hydchr)` unstructured format, (b) `TRANSIENT_HFB` keyword and the Fortran `IHFBRD` reuse/read flag semantics. Non-parametric static and transient HFB load/write are semantic; parameterized HFB (`NPHFB > 0`) now fails explicitly until preservation/expansion is implemented. |

#### New: `MfUsgTvm` — TVM2 full semantic implementation

`mfusgtvm.py` was rewritten from a verbatim text round-tripper into a complete
semantic implementation verified against `tvmu2.f`:

- Parameters: `itvmprint`, `tvmlogbasehk/vka/ss/sy`, `tvmddftr`, `tvmlogbasepor`
- `stress_period_data`: `dict[int, dict[str, np.recarray]]` — boundary index 0..nper,
  properties `hk/vka/ss/sy/ddftr/por`, nodes 0-based
- Transport-aware: with BCT → 7 global + 6 SP fields; without BCT → 6 + 5 fields
- nper+1 boundary blocks confirmed against `TVMU2AR` + `TVMU2RP` call structure
- Fortran fixed-format `(I10,F10.0)` per record confirmed
- 20 autotests in `autotest_local/test_mfusg_tvm.py` (all pass)

#### LPF Richards fix

| File | Fix |
|---|---|
| `flopy/mfusg/mfusglpf.py` | `Util2d.__eq__` returns `False` for any non-`Util2d` argument, so `if self.laytyp == 5:` was always `False`. This silently disabled Richards-equation array initialization (ALPHA/BETA/SR/BROOK) for any model with LAYTYP=5, even though the load path read them correctly. Fixed: replaced with `np.any(self.laytyp.array == 5)`. |

#### Bugs found in Fortran source audit — all fixed

| File | Bug | Severity | Fortran reference |
|---|---|---|---|
| `mfusgrch.py:396` | `t.index("INIZNRCH")` — `INIZNRCH` is an internal Fortran variable, not the file keyword. Raises `ValueError` on any model using RTS recharge zones. | **Critical** | `gwf2rch8u1.f` keyword is `INRCHZONES` |
| `mfusgdpt.py:355,507` | DLIM written/loaded when only `idpf` is True. Fortran requires `IDPF≠0 AND IDISPIM≠0`. Models with `idpf=1, idispim=0` wrote an extra array, shifting all subsequent reads. | High | `gwt2dptu1.f` line 281: `IF(IDPF.NE.0.AND.IDISPIM.NE.0)THEN` |
| `mfusgddf.py` | ISHARP (sharp-interface model flag) missing entirely from `__init__`, `write_file`, and `load`. ISHARP is the 6th numeric field after IMPHDD. Loading any file with ISHARP≠0 silently lost the flag; round-trip discarded it. | High | `density.f` line 64: `CALL URWORD(...,ISHARP,...)` |
| `mfusgdpf.py:191` | `f_obj` unbound `NameError` when `write_file(f=<open handle>)` called directly. The `else: f_obj = f` branch was missing. | Low | — |

## Install

```bash
conda create -n usgt -y python=3.12 numpy pandas matplotlib jupyter shapely
conda activate usgt
git clone --branch develop https://github.com/roh-sgm/flopy.git
cd flopy
pip install -e .
```

## Quick check

```python
import flopy
print(flopy.__version__)                              # 3.11.0.dev0 (or newer)
print('MfUsgTib' in dir(flopy.mfusg))                # True
print('MfUsgChd' in dir(flopy.mfusg))                # True
print('MfUsgRiv' in dir(flopy.mfusg))                # True
print('MfUsgGhb' in dir(flopy.mfusg))                # True
print('MfUsgDrn' in dir(flopy.mfusg))                # True
print('MfUsgTvm' in dir(flopy.mfusg))                # True
print('MfUsgGsf' in dir(flopy.mfusg))                # True
print('MfUsgEts' in dir(flopy.mfusg))                # True
print('MfusgTransportListBudget' in dir(flopy.utils)) # True

# Multi-species transport budget reader
from flopy.utils import MfusgTransportListBudget
lst_s1 = MfusgTransportListBudget("model.lst", species=1)
lst_s2 = MfusgTransportListBudget("model.lst", species=2)
inc, cum = lst_s1.get_budget()
df_inc, df_cum = lst_s2.get_dataframes(start_datetime="2000-01-01")
```

All other flopy APIs behave exactly as upstream; the patches only change the
text that `MfUsg.load(...).write_input()` produces for USG-T models.

## Validation

Exercised end-to-end against two real-world unstructured USG-T models built
with Groundwater Vistas, using their Vistas-generated input packages as the
reference:

**Model A — BCT + CLN + TIB transport** (≈ 112 k nodes, 493 stress periods, USG-T 1.8):
FloPy-rewritten input produces a listing file matching the reference bit-for-bit on every
budget column — flow and transport (`max|diff| = 0`, 20 flow + 16 mass components).
Wall-clock runtime equivalent to the reference.

**Model B — BCT IDISP=2 + DDF density-coupled transport** (≈ 19 k nodes, 1382 stress periods,
4322 time steps, USG-T 2.7 ARM):

| Output | Time steps | Max difference |
|--------|-----------|----------------|
| LST budget | 4322 / 4322 ✓ | 0.000e+00 |
| HDS heads | 4322 / 4322 ✓ | 0.000e+00 m |
| CON concentrations | 4322 / 4322 ✓ | 0.000e+00 |
| CBB (7 record types) | 4322 / 4322 ✓ | 0.000e+00 |

**Both models: bit-for-bit identical to the Vistas reference on all outputs.**

Tested binaries:
- USG-T 2.7.0 (ARM) — Model B full validation
- USG-T 1.8 (ARM and x86) — Model A full validation

## Upstream

Held on the fork while testing continues. No PRs to `modflowpy/flopy` have
been opened yet.
