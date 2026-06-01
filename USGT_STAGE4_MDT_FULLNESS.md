# MDT Fullness Card C — executed (USG-T 2.7)

Scope: audit `MfUsgMdt` against USG-T 2.7 Matrix Diffusion Transport and decide
honestly whether MDT can be promoted to `Full`. This is the MDT slice of the
`USGT_STAGE4_03_OC_EVT_MDT_LAK.md` pass; OC / EVT / LAK are **not** touched here
(beyond docs).

## Fortran audit (gwt2mdtu1.for, `gwt2mdtu1ar`)

- **Header line 1**: `IMDTCB IMDTCF`. When `IDPF == 0` (single-domain flow) the
  reader then parses, on the same line, the options `FRAHK` (IFRAHK=1),
  `FRADARCY` (IFRAHK=2), `TSHIFTMD <ts>`, `SEPARATE_AI2 <unit>` (IMDTCF2), and
  `MULTIFILE_MD <crootname>`. With `IDPF != 0` the option loop is **skipped**.
  `MULTIFILE_MD` requires `IMDTCF > 0` (the Fortran STOPs otherwise).
- **Base arrays** (per layer): `MDFLAG` (U1DINT, types 1–7); `VOLFRACMD`
  (U1DREL) **only when `IDPF == 0`** — with `IDPF != 0` it comes from DPF/PHIF
  and is neither written nor read here; then `PORMD`, `RHOBMD`, `DIFFLENMD`,
  `TORTMD`.
- **Per-species arrays** (loop over components, per layer): `KDMD`, `DECAYMD`,
  `YIELDMD`, `DIFFMD`; and, only when `TSHIFTMD > 1e-10`, `AIOLD1MD` and
  `AIOLD2MD`.
- `FRAHK` and `FRADARCY` both set the single `IFRAHK` flag (1 vs 2), so they are
  mutually exclusive in practice.

## Bugs fixed (`flopy/mfusg/mfusgmdt.py`)

1. **`load` dropped `FRAHK`/`FRADARCY`.** It upper-cased the line but searched
   for the lowercase keywords, so the options were never reloaded (e.g. the Ex7
   Multispecies model, whose header is `0 0 FRAHK`, lost the flag on round-trip).
   Keyword matching is now case-insensitive.
2. **Debug `print(f"line={line}")`** removed from `load`.
3. **`write_file(f=handle)` crashed** — `f_obj` was only assigned when
   `f is None`, raising `NameError` for a caller-supplied handle (and it closed
   an external handle). It now uses the handle and only closes what it opened.
4. The rootname (`MULTIFILE_MD`) case is preserved on load (it was upper-cased).

## Validation / behavior added

- `__init__` validates: `FRAHK` + `FRADARCY` are mutually exclusive; the header
  options require `IDPF == 0`; `MULTIFILE_MD` requires `imdtcf > 0`; and the
  per-species lists (`kdmd`/`decaymd`/`yieldmd`/`diffmd`, plus `aiold1md`/
  `aiold2md` when `tshiftmd>0`) must have `MCOMP` entries (clear `ValueError`
  instead of an `IndexError`).
- `load` parses the header options only when `IDPF == 0`, matching the Fortran.

## Tests

Synthetic (`autotest/test_usg_transport.py`, `-k mfusgmdt`):

- `test_mfusgmdt_minimal_authoring_roundtrip` (mcomp=1, tshiftmd=0; semantic
  reload + rewrite stability),
- `test_mfusgmdt_tshift_aiold_roundtrip` (TSHIFTMD>0 with AIOLD1/AIOLD2),
- `test_mfusgmdt_multispecies_roundtrip` (mcomp=2, distinct per-component values),
- `test_mfusgmdt_frahk_fradarcy_roundtrip` (each option; regression for the load
  case bug),
- `test_mfusgmdt_output_options_roundtrip` (SEPARATE_AI2 + MULTIFILE_MD, rootname
  case preserved),
- `test_mfusgmdt_idpf_skips_volfracmd` (VOLFRACMD dropped under IDPF!=0),
- `test_mfusgmdt_external_handle_write` (regression for the f_obj bug),
- `test_mfusgmdt_rejects_invalid` (frahk+fradarcy, list length, MULTIFILE needs
  imdtcf, IDPF+options).

Real models: the three Ex7 Matrix-Diffusion tests (SandTank, Multispecies,
DiscreteFracture) remain as load → write → run round-trips under the USG-T 2.7
executable, so MDT execution is exercised on real DISU models. No separate
from-scratch MDT executable smoke was added: a minimal convergent BCT+MDT model
is not cheap/stable to build, and the Ex7 exe runs already cover MDT execution.

`-k mfusgmdt` **8 passed**; focused **154 passed**; exe **4 passed**; combined
**158 passed** under the USG-T 2.7 ARM binary.

## Decision: keep `✅ (intentionally not Full)`, strongly hardened

The two prior reasons for not-Full are now resolved: the field order is
Fortran-audited (above) and there is from-scratch authoring with tests for every
main branch (options, base arrays, per-species arrays, AIOLD under TSHIFTMD,
IDPF on/off, structured + DISU via Ex7). Three real bugs were fixed. Two bounded
gaps keep the honest label at `✅ (intentionally not Full)`:

- **Species count.** The per-species loop uses `model.mcomp` (mobile
  components). USG-T's `NTCOMP` can exceed `MCOMP` for chained-decay /
  parent–daughter setups; that case is not independently verified (the tested
  models, including Ex7 Multispecies, have `NTCOMP == MCOMP`).
- **Output binaries.** `MULTIFILE_MD`/`SEPARATE_AI2` AI1/AI2 binary outputs are
  authored (keywords/units round-trip) but not read or post-processed by FloPy.

## Review follow-up (resolved)

P1 — TSHIFTMD could produce misaligned MDT files. USG-T reads AIOLD1MD/AIOLD2MD
only when `TSHIFTMD > 1e-10`, but the writer used `tshiftmd > 0` for the arrays
and a fixed `{tshiftmd:9.2f}` format for the keyword. So `tshiftmd=1e-6` wrote
`TSHIFTMD 0.00` *and* the AIOLD arrays — USG-T then read TSHIFTMD as 0.0, did not
expect AIOLD, and the file desynced; `tshiftmd=1e-12` likewise wrote AIOLD below
the solver threshold. Fixed:

- a module constant `MDT_TSHIFT_THRESHOLD = 1.0e-10` is now used consistently in
  the constructor (AIOLD built/required), `write_file` (TSHIFTMD keyword + AIOLD
  written), the IDPF-options check, and `load` (AIOLD read);
- the TSHIFTMD keyword is written with a general format (`{tshiftmd:15.7g}`) so a
  small valid value like `1e-6` is not rounded to `0.00`.

Test added: `test_mfusgmdt_tshiftmd_threshold` (tshiftmd=2.0 writes TSHIFTMD +
AIOLD and round-trips; 1e-6 writes a non-zero readable value + AIOLD; 1e-12
writes neither; IDPF!=0 with 1e-12 is inactive/no error, with 1e-6 raises).
`-k mfusgmdt` **9 passed**; focused **155 passed**; exe **4 passed**; combined
**159 passed** under the USG-T 2.7 ARM binary. Status unchanged
(`✅ intentionally not Full`).

## Validation

```bash
python -m pytest autotest/test_usg_transport.py -k mfusgmdt -q
python -m pytest autotest/test_usg_transport.py -q
python -m pytest autotest/test_usg_transport_exe.py -q
USGT_EXE=/Users/roh.sgm/Documents/GitHub_Projects/GW-Software-Compiled-via-Claude/gfortran/usgt_2.7/usgt_270_arm python -m pytest autotest/test_usg_transport.py autotest/test_usg_transport_exe.py -q
git diff --check
git status --short
```
