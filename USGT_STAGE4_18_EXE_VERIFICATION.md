# Stage 4.7B — OC / MDT executable verification (USG-T 2.7)

Date: 2026-06-03 · Base: `develop` @ `a7ea458f`

Goal: close or honestly classify the remaining **executable-verification** gaps
for OC and MDT (the post-A5 items, not authoring gaps). Add a cheap exe smoke
where it is stable; document the rest as a **manual tier** with the exact Fortran
reason. No package features changed.

Executable used: `USGT_EXE=…/usgt_2.7/usgt_270_arm` (USG-T 2.7 ARM), via
`@requires_exe(USGT_EXE)` / `exe_name=USGT_EXE` (skips cleanly when absent).

## OC audit (`glo2basu1.f` `SGWF2BAS7I`)

The OC option line is parsed at `SGWF2BAS7I` (`:954`+). Findings:

- **`ATSA` / `ATS`** (`:1068-1071`, `IATS=1`): adaptive time-stepping. When set,
  the Fortran **overrides** `NSTP` (1 for steady, 1e6 for transient) and uses its
  own built-in `DELTAT=1.0`, `TMINAT=1e-10`, `TMAXAT=1e10`, `TADJAT=2.0`,
  `TCUTAT=5.0` (`:1155-1171`) — the OC file only needs the `ATSA` keyword (no
  extra ATS variables are read from it). **Cheap + stable → added an exe smoke.**
- **`FASTFORWARD` / `FASTFORWARDC`** (`:1080-1125`): read GWF (and CLN/DDF) heads
  — or, for `FASTFORWARDC`, concentrations — from a **separate external file**
  produced by a *prior* run, to fast-forward to a stress period. **Manual tier:**
  a faithful smoke needs a two-stage fixture (run model A, capture its head/conc
  file, run model B reading it); a synthetic one would be brittle.
- **`BOOTSTRAPPING`** (`:1127-1147`): reads BOOTSTRAP GWF (and CLN/DDF) heads from
  a separate external file and allocates bootstrap arrays for a **transient**
  run, to seed the first-iteration head estimate. **Manual tier:** same external
  prior-run state as FASTFORWARD; not reasonable to synthesize cheaply/stably.
- **`SAVE IBOUND`**: USG-T 2.7's solver **rejects** it (the dispatch is commented
  out), so it is not exe-testable. FloPy preserves the keyword on write and
  `check()` warns. Decision unchanged (honest; documented in the OC Fullness
  card).

## MDT audit (`gwt2mdtu1.for`)

- The MDT array reads loop over **`NTCOMP`** (`:251`, `KDMD`/`DECAYMD`/`YIELDMD`/
  `DIFFMD`, `AIOLD1/2`), while the FM/BD loop uses `MCOMP` (`:451`). The
  `NTCOMP>MCOMP` chained-decay case is therefore a **BCT chained-decay**
  configuration (immobile/extra components) flowing through MDT, not an MDT-only
  knob.
- **`SEPARATE_AI2` / `MULTIFILE_MD`** (`:113-159`): per-component **binary**
  air-interface output files (`IMDTCF`/`IMDTCF2`, opened per `ICOMP` over
  `NTCOMP`). Validating them means reading USG-T node-based binaries per
  component.
- **`TSHIFTMD`**: already fixed and round-trip tested (Stage 4 Card C).

**Manual tier for MDT.** A convergent matrix-diffusion model authored from
scratch is not cheap (it needs BCT + dual-domain config + many MDT arrays + a
non-trivial solve), and **the three real Ex7 matrix-diffusion models already
load → write → run under USG-T 2.7**, so MDT *execution* is already covered. The
specific open items — `NTCOMP>MCOMP` chained decay and the `AI1/AI2` binary
outputs — are not cheap/stable to smoke (chained decay needs a tailored BCT
setup; the AI outputs are per-component binaries). They stay a documented manual
tier rather than a brittle synthetic test.

## What changed

- `autotest/test_usg_transport_exe.py`: new `test_usgt_exe_oc_atsa_from_scratch`
  — a transient 1-D model with `MfUsgOc(atsa=1)` runs under USG-T 2.7 and reaches
  the analytical linear gradient `[8, 6.5, 5, 3.5, 2]` with fixed CHD boundaries.
  (The budget of an ATS-transient listing has many adaptive steps and is not
  asserted; the head field is the reliable cross-check.)
- Docs only otherwise: this stage doc + gap-audit / backlog / improvements
  updates + the OC roadmap row qualified (ATSA now exe-verified).

No package code changed; no MF6-TID, binaries, or `examples/data` touched.

## Tests

```bash
python -m pytest autotest/test_usg_transport_exe.py -q                  # 5 passed
USGT_EXE=.../usgt_270_arm python -m pytest autotest/test_usg_transport.py autotest/test_usg_transport_exe.py -q   # 309 passed
USGT_EXE=.../usgt_270_arm python -m pytest autotest/test_usg_transport.py autotest/test_usg_transport_exe.py -k "mfusgoc or mfusgmdt or usgt_exe" -q   # 25 passed
git diff --check
```

## Status — exe verification

| Item | Verdict |
|---|---|
| OC `ATSA` execution | ✅ exe smoke (`test_usgt_exe_oc_atsa_from_scratch`) |
| OC `FASTFORWARD`/`FASTFORWARDC` | **manual tier** — external prior-run head/conc file |
| OC `BOOTSTRAPPING` | **manual tier** — external prior-run head file + transient |
| OC `SAVE IBOUND` | not exe-testable — USG-T 2.7 solver rejects it (check() warns) |
| MDT execution (general) | ✅ covered by the three real Ex7 models (round-trip/run) |
| MDT `NTCOMP>MCOMP` chained decay | **manual tier** — needs a tailored BCT chained-decay setup |
| MDT `AI1/AI2` (`SEPARATE_AI2`/`MULTIFILE_MD`) | **manual tier** — per-component binary outputs |

OC and MDT remain `✅ (intentionally not Full)`: the open items are either
exe-verified now (ATSA), covered by the real Ex models (MDT execution), or
honestly documented as manual tier (external-state / binary / solver-rejected
cases). No `Full` over-claim.
