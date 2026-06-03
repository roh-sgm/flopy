# Stage 4.6F-B — EVT ETS zonal time-series (`ETS MXZNEVT` / `INEVTZONES`)

Date: 2026-06-03 · Base: `develop` @ `b7df0fe8`

Goal: close the last `MfUsgEvt` gap — the `ETS MXZNEVT` / `INEVTZONES`
zonal ET time-series. **Decision after the Fortran audit: deferred** — it is a
dynamic, ATS-coupled *execution* mode (an external time-series file that
supersedes the EVTR array), not the static array I/O FloPy models. This stage
records the full Fortran spec and replaces the single generic
`NotImplementedError` with **specific, per-branch explicit failures (no partial
parse)**.

## Fortran audit (`gwf2evt8u1.f`, ATS via `glo2basu1.f`)

### `GWF2EVT8U1AR` — item 1/2 (`:76-95`)

After `NEVTOP IEVTCB [IETFACTOR]`, the option scanner accepts an `ETS MXZNEVT`
keyword:

```
IF(LINE(ISTART:ISTOP).EQ.'ETS') THEN
  IF(IATS.EQ.0) THEN ; WRITE(...,'ET TIME-SERIES NEEDS ADAPTIVE TIME-STEPPING.STOPPING') ; STOP ; ENDIF
  ALLOCATE(MXZNEVT) ; CALL URWORD(...,MXZNEVT,...) ; IETSOPT = 1
END IF
```

- **`ETS` requires adaptive time-stepping** — the Fortran **`STOP`s** when
  `IATS==0`. `IATS` is set by the OC/BAS ATS sub-system (`glo2basu1.f`), which
  FloPy does not model as execution.
- `MXZNEVT` = max number of ET zones; `IETSOPT=1` and the banner states
  *"TIME-SERIES OF ET SUPERSEDES ARRAY INPUT VIA STRESS PERIODS"*.

### `GWF2EVT8U1RP(IN,IUETS,KPER)` — per stress period (`:179-461`)

- The header may begin with `INEVTZONES <INIZNEVT>` (`:217-223`), read before the
  positional `INSURF INEVTR INEXDP [INIEVT]`.
- `IETSOPT==0` → skip the zonal block. Otherwise (`:427-460`):
  - `INIZNEVT < 0`: reuse the previous `IZNEVT`.
  - `INIZNEVT >= 0`: **`MXZNEVT = INIZNEVT`** (the per-SP flag *redefines* the zone
    count), then read the `IZNEVT` zone-index array per cell (`U2DINT`,
    `NROW×NCOL` or `1×NIEVT`).
  - On `KPER==1`: `inets = IUETS`; read **one** record from the external ETS
    file `Tstart, Tend, Factor, Ets(1..MXZNEVT)`; set `TIMEVT = TENDEVT`.
  - **Update EVTR dynamically**: `EVTR(NN) = etsevt(IZNEVT(n)) * AREA(N) *
    Factrevt` — the rate is taken from the time-series by zone.

The ETS file is consumed **progressively during the run** (the next record is
read when simulation time reaches `TENDEVT`), driven by the ATS machinery — not
in one pass at load time.

## Why it is deferred (not safe as static I/O)

1. **Hard ATS dependency** — `STOP` when `IATS==0`. Writing an `ETS` EVT without
   an ATS-enabled OC would produce a Fortran `STOP`; FloPy does not model ATS
   execution.
2. **External progressive time-series file (`IUETS`)** — records
   `Tstart Tend Factor Ets(MXZNEVT)` read step-by-step as the clock advances. This
   is execution state, not a package array.
3. **`IETSOPT=1` supersedes the array EVTR** — EVTR is recomputed every step from
   the time-series; there is no static EVTR array to preserve or round-trip.
4. **`INEVTZONES` redefines `MXZNEVT` per SP** and (re)reads `IZNEVT`, an unusual,
   stateful per-period control.
5. **Conflicts with `NPEVT`** — both define EVTR (parameters vs time-series); the
   modes are mutually exclusive in practice.

Modeling this faithfully would mean authoring/parsing an external time-series
file *and* simulating the ATS-coupled progressive read — outside FloPy's static
EVT I/O. Implementing only the array side would be a partial parse that
mis-represents the file, which the project policy forbids.

## What changed (code, `flopy/mfusg/mfusgevt.py`)

No semantic ETS-zonal support added; the generic `NotImplementedError` is
replaced by **specific, per-branch** failures, each before any partial read:

- **`__init__`** (`mxetzones>0`): a detailed `NotImplementedError` citing the ATS
  requirement, the external `IUETS` time-series, `IETSOPT` superseding EVTR, and
  the per-SP `INEVTZONES`/`IZNEVT`.
- **`load`** — `ETS MXZNEVT` on item 2: raises **before** reading any
  stress-period data (closes the file first; no partial parse).
- **`load`** — a per-SP `INEVTZONES` flag in a stress-period header: raises
  explicitly (rather than the previously commented-out / ignored handling).

`NPEVT` (Stage 4.6F-A) and all plain EVT paths are unaffected.

## Tests (`autotest/test_usg_transport.py`, `-k mfusgevt`)

- `test_mfusgevt_ets_zonal_deferred` — authoring `mxetzones>0` →
  `NotImplementedError` (cites ATS); load of `ETS MXZNEVT` (item 2) → raises
  before any SP parse; load of a per-SP `INEVTZONES` header → raises; a plain EVT
  (no ETS) still authors and round-trips.

The Stage 4.4 Card-B EVT tests, the Stage 4.6F-A NPEVT tests (+ hardening), and
the `test_usgt_exe_evt_from_scratch` exe smoke stay green. `-k mfusgevt` **16**,
focused **303**, exe **4**, combined **307** (USG-T 2.7 ARM).

## Status after Stage 4.6F-B — EVT stays `✅ (intentionally not Full)`

With NPEVT closed (4.6F-A) and ETS-zonal deferred (here), EVT's static I/O is
complete and broadly tested. EVT is **not promoted to `Full`** for one honest
reason: the `ETS MXZNEVT` zonal time-series is an **ATS-coupled dynamic execution
mode** (external time-series file superseding EVTR) that FloPy does not model;
it now fails explicitly per branch with the full spec recorded here. (Promoting
to `Full` would overclaim coverage of a feature that is intentionally
unsupported.)

## Validation

```bash
python -m pytest autotest/test_usg_transport.py -k mfusgevt -q   # 16 passed
python -m pytest autotest/test_usg_transport.py -q               # 303 passed
python -m pytest autotest/test_usg_transport_exe.py -q           # 4 passed
USGT_EXE=.../usgt_270_arm python -m pytest autotest/test_usg_transport.py autotest/test_usg_transport_exe.py -q   # 307 passed
git diff --check
```
