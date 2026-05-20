# USG-Transport FloPy — Package Roadmap

Tracking which USG-T packages are supported in this fork, what is missing, and
implementation dates. Packages already in upstream flopy/mfusg are not listed
unless we changed them.

---

## Status legend

| Symbol | Meaning |
|--------|---------|
| ✅ | Implemented and tested in this fork |
| 🔲 | Not yet started |
| ⏭ | Intentionally deferred (see notes) |

---

## Packages already in upstream flopy/mfusg (unchanged here)

These were already present before any work in this fork.

| Package | Class | Notes |
|---------|-------|-------|
| BAS6 | `MfUsgBas` | `UNSTRUCTURED` write fix (2026-05); `CONVERGE` load+write fix (2026-05-19) |
| DIS | `MfUsgDis` | Structured grid |
| DISU | `MfUsgDisU` | Unstructured grid |
| BCF6 | `MfUsgBcf` | |
| LPF | `MfUsgLpf` | |
| WEL | `MfUsgWel` | |
| RCH | `MfUsgRch` | CONC transport fixes (2026-05): INCONC keyword spacing, trailing `1` removed; see RCH note |
| EVT | `MfUsgEvt` | |
| ETS | `MfUsgEts` | ✅ Full implementation including NETSEG>1 (2026-05) |
| GNC | `MfUsgGnc` | |
| LAK | `MfUsgLak` | |
| OC | `MfUsgOc` | |
| CLN | `MfUsgCln` | None-unitnumber fix (2026-05) |
| BCT | `MfUsgBct` | ✅ IDISP=1 DISU verified (2026-05); IDISP=2 DISU+DLX-DTXZ confirmed present and correct (2026-05-18) |
| PCB | `MfUsgPcb` | |
| TIB | `MfUsgTib` | ✅ New (2026-05) |
| MDT | `MfUsgMdt` | |
| DPF | `MfUsgDpf` | |
| DPT | `MfUsgDpt` | |
| SMS | `MfUsgSms` | |
| DDF | `MfUsgDdf` | ✅ NONLINEAR option + index bugfix (2026-05); `ithickav` default=0 fix (2026-05-19) |

---

## New packages added in this fork

| Package | Class | Status | Implemented |
|---------|-------|--------|-------------|
| CHD | `MfUsgChd` | ✅ | 2026-05-15; `shead`/`ehead` float64 precision fix (2026-05-19) |
| RIV | `MfUsgRiv` | ✅ | 2026-05-15; `stage` float64 precision fix (2026-05-19) |
| GHB | `MfUsgGhb` | ✅ | 2026-05-16 |
| DRN | `MfUsgDrn` | ✅ | 2026-05-16 |
| TVM2 | `MfUsgTvm` | ✅ | 2026-05-16 |
| GSF | `MfUsgGsf` | ✅ | 2026-05-16 |

---

## New utilities added in this fork

| Utility | Class | Status | Implemented |
|---------|-------|--------|-------------|
| Multi-species transport budget | `MfusgTransportListBudget` | ✅ | 2026-05-15 |

---

## Deferred / out of scope

| Package | Notes |
|---------|-------|
| STR | User prefers CLN for surface water |
| SFR2 | User prefers CLN for surface water |
| OBS | Low priority |
| HOB | Low priority |

---

## Notes

- All new USG-T packages support **unstructured grids** (node-based) and **AUX
  concentration variables** for transport.
- GHB and DRN follow the same read/write pattern as RIV and CHD: `MfList`-based
  stress period data, `ITMP NP` per-SP header, AUX names preserved exactly as
  written in the file.
- TVM2 uses a text-block round-trip (same approach as TIB): the global header
  and per-SP content are stored verbatim so the file survives load/write without
  semantic parsing of the material arrays.
- GSF is a utility wrapper: text round-trip for the Grid Specification File, with
  a `to_grid()` method that delegates to `UnstructuredGrid.from_gridspec()`.
- All new packages are registered in `MfUsg.mfnam_packages` and exported from
  `flopy.mfusg`.

### RCH transport note (2026-05-18)

Two bugs fixed in `mfusgRCH.py` per-SP write for transport models:

1. **`INCONC` keyword spacing** (`mfusgrch.py` line 225): the `# Stress period N`
   comment was appended directly to `INCONC` with no separator, producing
   `INCONC# Stress period N`. USG-T parses `INCONC#` as one unrecognized token
   and silently ignores the flag, then crashes looking for the next SP header where
   the RCHCONC array sits. Fix: changed `f_rch.write(f"# Stress period ...")` to
   `f_rch.write(f" # Stress period ...")` (space before `#`).

2. **`INCONC 1` → `INCONC`** (line 224): the manual says `INCONC` is a bare
   keyword with no value. The trailing `1` was removed.

### NAM file write fixes (2026-05-18)

Three bugs fixed in the round-trip NAM writer (`mf.py` + `mfusg.py`):

1. **`REPLACE` keyword removed**: USG-T does not use `REPLACE` in NAM. FloPy was
   writing `DATA(BINARY) 50 MDV.cbb REPLACE` for output unit files. Removed from
   the `output_units` write loop.

2. **AFR filename**: WEL's `add_output_file(iunitafr, fname=None, extension="afr")`
   generated `MDV.afr` instead of the actual `MDV_FlowReduction.dat`. Fixed in
   `_prepare_external_files` by overriding `output_fnames[idx]` with the real
   filename from ext_unit_dict when the unit is already registered.

3. **SMS/OC in NAM**: Packages excluded from `load_only` were not written to NAM.
   Fixed by populating `model._skipped_nam_entries` for skipped packages so they
   appear verbatim in the NAM.

### DISU free_format_npl (2026-05-17)

For DISU models, FloPy now defaults `free_format_npl = 10` in `_load_packages`.
Without this, arrays with 1942+ nodes land on a single line (~29 000 chars),
which exceeds the internal buffer in `usgt_270_arm` (EOF errors at SP 453).

### BCT IDISP=2 DISU (2026-05-18)

Confirmed: `DLX`/`DLY`/`DLZ`/`DTXY`/`DTYZ`/`DTXZ` arrays ARE present in the
reference Soncor BCT (lines 19,434–27,293 of MDV.bct, DISU+IDISP=2). The I/O
manual is correct. An earlier attempt to guard these with `model.structured` was
reverted. FloPy now loads and writes them unconditionally for IDISP=2, which is
the correct behaviour.

### Density-coupled round-trip fixes (2026-05-19)

Four bugs surfaced during the soncor_ddf_test (BCT IDISP=2, DDF, 1382 SPs, density-coupled):

1. **DDF `ithickav` default** (`mfusgddf.py`): `type_from_iterable(..., default_val=1)` →
   `default_val=0`. USG-T treats the absent ITHICKAV field as 0 (arithmetic averaging);
   FloPy was loading it as 1 (thickness-weighted), changing transmissivity in the
   density layer → transport divergence at SP 469.

2. **CHD `shead`/`ehead` precision** (`mfusgchd.py`): dtype changed from `np.float32`
   to `np.float64`. Float32 round-trip of CHD heads (e.g. 2298.46 → 2298.459961)
   introduced ~0.04 mm errors per node. Over 469 density-coupled SPs this accumulated
   into local concentration differences that caused the outer nonlinear loop to
   marginally fail transport convergence (ΔC > CICLOSE=1e-8) at SP 470 TS 2.

3. **RIV `stage` precision** (`mfusgriv.py`): dtype changed from `np.float32` to
   `np.float64` for the same reason (e.g. stage 2300.002 → 2300.001953 in float32).

4. **BAS6 `CONVERGE` option** (`mfusgbas.py`): two bugs in one:
   - `load`: `converge=converge` was missing from the `cls(...)` constructor call, so
     `self.converge` was always `False` after loading a file that declared `CONVERGE`.
   - `write_file`: the `CONVERGE` keyword was never emitted even when `self.converge=True`.
   The `CONVERGE` option in BAS6 tells USG-T to use the coupled flow–transport
   convergence criterion in the outer nonlinear loop. Without it, SP 470 TS 2 failed
   after 250 iterations (reference accepted it) → time-step cascade → run failure.

---

## Test results (2026-05-18)

### md_bct_test — FloPy BCT (IDISP=1, DISU) vs data_2 reference (x86 mfusg_gsi_1_8)

Run complete. Compared via `usgmod2smp` (SMP files from MainDam.hds / MainDam.con):

| Output | n obs | max\|Δ\| | mean\|Δ\| | RMSE |
|--------|-------|---------|---------|------|
| Heads [m] | 64,090 | 0.0000e+00 | 0.0000e+00 | 0.0000e+00 |
| Conc [mg/L] | 54,723 | 0.0000e+00 | 0.0000e+00 | 0.0000e+00 |

**Result: BCT IDISP=1 DISU round-trip is bit-for-bit identical to the reference run.**

CLN heads/conc also compared via `MainDam_cln.hds` / `MainDam_cln.con`: **all zeros** (not shown above).

### soncor_ddf_test — FloPy BCT+DDF (IDISP=2, DISU) vs reference (ARM usgt_270_arm)

**Status: PASSED (2026-05-19) — 1382 stress periods, 4322 time steps, zero convergence failures.**

Tested against an unstructured USG-T 2.7 model with density-coupled flow and transport
(BCT IDISP=2, DDF with extreme density contrast ρ_fresh=1050 / ρ_std=1230 kg/m³).
The FloPy-rewritten input was run through `usgt_270_arm` and compared to the Vistas-generated
reference run.

| Output | Time steps | Max difference |
|--------|-----------|----------------|
| LST budget (flow) | 4322 / 4322 ✓ | 0.000e+00 |
| HDS heads | 4322 / 4322 ✓ | 0.000e+00 m |
| CON concentrations | 4322 / 4322 ✓ | 0.000e+00 |
| CBB (7 record types) | 4322 / 4322 ✓ | 0.000e+00 |

**Result: bit-for-bit identical to the Vistas reference on all outputs.**

Four bugs were found and fixed during this test (see fixes below).

**Binary comparison notes for DISU outputs:**

- Use `flopy.utils.HeadUFile` (not `HeadFile`) for DISU `.hds` files.
- Use `HeadUFile(..., text='CONC')` for DISU `.con` (concentration) files — the default
  `text='headu'` filter skips all records.
- `HeadFile` fails silently for DISU because its binary header stores `nstrt`/`nend`
  (1-based node range), not `ncol`/`nrow`; `HeadUFile` handles this correctly.

---

## Known gaps and missing features (as of 2026-05-19)

### Packages not implemented

| Package | Notes |
|---------|-------|
| HFB | Horizontal Flow Barrier. `modflow/mfhfb.py` exists upstream but no `mfusghfb.py`. Needs node-based adaptation for DISU. |
| FHB | Flow and Head Boundary. Not in `flopy/mfusg/`. Rarely used; can often be replaced by WEL+CHD. |
| MNW2 | Multi-Node Well 2. Not in `flopy/mfusg/`. |
| SFT | Streamflow Transport (companion to STR/SFR2). Not implemented; CLN is the preferred substitute. |
| LKT | Lake Transport (companion to LAK). `mfusglak.py` exists but the transport BC side is untested. |

### Existing packages with known gaps

| Package | Gap |
|---------|-----|
| BAS6 | `DPIN`/`DPOUT`/`DPIO` (double-precision I/O), `RICHARDS`, `PRINTTIME`, `SHOWPROGRESS`, `SY-ALL` are parsed on load but **never written**. Models using any of these flags will silently lose them on round-trip. Fix: add corresponding `if self.X: opts.append("X")` lines to `write_file`. |
| WEL | Per-SP header is written as `MXACT ITMPCWL` (2 tokens); reference files include a third token `ITMPCLN` (transport flag for CLN wells). USG-T defaults `ITMPCLN=0` when absent, so this is harmless when there are no CLN injection wells, but the output doesn't match Vistas exactly. |
| WEL | Injection-well AUX concentration columns: `_check_for_aux` exists but has not been round-trip tested against a USG-T file that includes `AUX` in the WEL header. |
| PCB, MDT, DPF, DPT | Code exists (upstream, 300–600 lines each) but these packages have not been tested in any round-trip run. |
| CLN | CLN transport boundary concentrations (injection through CLN pipes) not verified. |
