# OC Fullness Card A — executed (USG-T 2.7)

Scope: audit `MfUsgOc` against USG-T 2.7 Output Control and decide honestly
whether OC can be promoted to `Full`. This is the OC slice of the
`USGT_STAGE4_03_OC_EVT_MDT_LAK.md` pass; EVT / MDT / LAK are **not** touched
here. (Named separately to avoid renumbering the Stage 4 docs after
`gridgen_to_gsf` took the Stage 4.3 number.)

## Fortran audit (glo2basu1.f)

Output Control is read in two phases:

* **`SGWF2BAS7I`** reads the **first OC line** and parses the header options
  token-by-token: `ATS`/`ATSA`, `NPTIMES n`, `NPSTPS n`, `FASTFORWARD …`,
  `FASTFORWARDC …`, `BOOTSTRAPPING iug [iuc if CLN] [iud if DPF]`. The
  numeric-vs-alphabetic test also keys off this first line.
* **`SGWF2BAS7J`** reads the subsequent **alphabetic format records**
  (`HEAD/DRAWDOWN/IBOUND/CONC PRINT/SAVE FORMAT/UNIT`, `COMPACT …`) until
  `PERIOD`. It accepts `ATSA`/`FASTFORWARD`/`FASTFORWARDC` as standalone lines
  (skipped) but **not** `BOOTSTRAPPING` — a standalone `BOOTSTRAPPING` line
  makes USG-T STOP with "ERROR READING OUTPUT CONTROL". So `BOOTSTRAPPING` must
  live on the first OC line.
* **`SGWF2BAS7N`** reads the per-stress-period actions: `PRINT/SAVE BUDGET`,
  `PRINT/SAVE HEAD|DRAWDOWN|CONC[ENTRATION]` (with optional layer lists),
  `DELTAT/TMINAT/TMAXAT/TADJAT/TCUTAT`, `HCLOSE/BTOL/MXITER`, and the bootstrap
  toggles `BOOTSTRAP/NOBOOTSTRAP/BOOTSTRAPSCALE/NOBOOTSTRAPSCALE`. `DDREFERENCE`
  rides on the `PERIOD … STEP …` line. **`SAVE IBOUND` is commented out** in
  this USG-T 2.7 build, so the solver does not accept it.

## What was implemented

`flopy/mfusg/mfusgoc.py`:

1. **`BOOTSTRAPPING` header** — `__init__` now accepts `bootstrapping` +
   `iugboot`/`iucboot`/`iudboot` (and infers `bootstrapping` from a non-zero
   `iugboot` on load); `write_file` emits `BOOTSTRAPPING …` **on the first OC
   line** (appended to the `ATSA` line, or as line 1 when there is no `ATSA`),
   matching `SGWF2BAS7I`. Previously `load` read the units but `__init__` /
   `write_file` dropped them, so the header was lost on round-trip — now fixed.
2. **Layer-qualified actions** — `load` keeps the whole action line, so
   `SAVE/PRINT HEAD|DRAWDOWN|CONC <layers>` survive (it previously truncated to
   the first two tokens and lost the layer list).
3. **`DDREFERENCE`** — `load` recaptures it from the `PERIOD` line as an action,
   so it round-trips (previously dropped).

Per-stress-period `BOOTSTRAP/NOBOOTSTRAP/BOOTSTRAPSCALE/NOBOOTSTRAPSCALE` and
`SAVE/PRINT CONC|BUDGET` already round-tripped through the verbatim action-list
mechanism; this card adds explicit tests.

## Tests (`autotest/test_usg_transport.py`, `-k mfusgoc`)

`test_mfusgoc_bootstrapping_header_authoring_roundtrip`,
`test_mfusgoc_bootstrap_stress_period_actions`,
`test_mfusgoc_output_block_combinations`,
`test_mfusgoc_save_ibound_roundtrips_but_usgt_rejects`,
`test_mfusgoc_layer_qualified_roundtrip`, `test_mfusgoc_ddreference_roundtrip`
(plus the existing `test_mfusgoc_atsa_authoring_roundtrip`). `-k mfusgoc`
**7 passed**; focused suite **134 passed**; exe **3 passed**; combined
**137 passed** under the USG-T 2.7 ARM binary.

## Decision: keep `✅ (intentionally not Full)`, hardened

OC authoring/round-trip is now broad (BOOTSTRAPPING header, per-SP bootstrap
toggles, layer-qualified actions, DDREFERENCE, CONC/BUDGET, ATSA/NPTIMES/NPSTPS,
COMPACT), but explicit gaps remain — so the honest status stays
`✅ (intentionally not Full)`:

- **`SAVE IBOUND`** — FloPy preserves the keyword (standard MODFLOW), but USG-T
  2.7's OC reader has it commented out and will reject it. Do not use with
  USG-T 2.7.
- **`FASTFORWARD` / `FASTFORWARDC`** — written on their own lines (pre-existing
  behaviour). Per the audit, USG-T parses these options from the first OC line;
  separate-line placement is **not verified** to activate the fast-forward, and
  was not changed in this card.
- **`BOOTSTRAPPING`** — placement is now Fortran-correct (first line) and the
  FloPy round-trip is tested, but execution was **not** validated against the
  USG-T 2.7 executable (it needs a transient run with a bootstrap heads file,
  out of scope for a cheap smoke test).
- A numeric-format OC file loads correctly but is rewritten in the words format
  (semantically equivalent, not byte-identical).

These are documented, bounded gaps rather than silent failures.

## Polish follow-up (resolved)

1. **DDREFERENCE-only period** — `write_file` now emits the
   `period <kper> step <kstp> ddreference` line even when that period has no
   other action (previously the period was dropped and the flag leaked into the
   next written period). Test: `test_mfusgoc_ddreference_only_roundtrip`.
2. **`check()` recognises USG-T OC actions** — no more false "ignored" warnings
   for `BOOTSTRAP`/`NOBOOTSTRAP`/`BOOTSTRAPSCALE`/`NOBOOTSTRAPSCALE`,
   `DDREFERENCE`, the ATS params (`DELTAT`/`TMINAT`/`TMAXAT`/`TADJAT`/`TCUTAT`),
   and the solver params (`HCLOSE`/`BTOL`/`MXITER`). Genuinely unknown actions
   are still flagged. Test: `test_mfusgoc_check_accepts_usgt_actions`.
3. **`SAVE IBOUND` check warning** — `check()` now emits a specific warning,
   "SAVE IBOUND is preserved by FloPy but rejected by USG-T 2.7", so the
   documented solver gap surfaces at check time (write/load still preserve the
   keyword). Test: `test_mfusgoc_check_warns_save_ibound`.

After the polish: `-k mfusgoc` **10 passed**; focused **137 passed**; exe
**3 passed**; combined **140 passed** under the USG-T 2.7 ARM binary.

## Validation

```bash
python -m pytest autotest/test_usg_transport.py -k mfusgoc -q
python -m pytest autotest/test_usg_transport.py -q
python -m pytest autotest/test_usg_transport_exe.py -q
USGT_EXE=/Users/roh.sgm/Documents/GitHub_Projects/GW-Software-Compiled-via-Claude/gfortran/usgt_2.7/usgt_270_arm python -m pytest autotest/test_usg_transport.py autotest/test_usg_transport_exe.py -q
git diff --check
git status --short
```
