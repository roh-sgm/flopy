# EVT Fullness Card B — executed (USG-T 2.7)

Scope: audit `MfUsgEvt` against USG-T 2.7 Evapotranspiration and decide honestly
whether EVT can be promoted to `Full`. This is the EVT slice of the
`USGT_STAGE4_03_OC_EVT_MDT_LAK.md` pass; OC / MDT / LAK are **not** touched here
(beyond docs).

## Fortran audit (gwf2evt8u1.f)

`GWF2EVT8U1AR` (dataset 1):

- without transport: `NEVTOP IEVTCB`;
- with transport (`INBCT>0`): `NEVTOP IEVTCB IETFACTOR` — always 3 integers;
- optional `ETS MXZNEVT` keyword on the same line (requires ATS, else STOP);
- `NEVTOP` must be 1..3 (else abort);
- unstructured + `NEVTOP=2`: a `MXNDEVT` line follows;
- transport `ETFACTOR(MCOMP)`: `IETFACTOR==0` → 0.0, `IETFACTOR<0` → 1.0,
  `IETFACTOR>0` → read MCOMP values on a separate line;
- `NPEVT>0`: named parameters (`UPARARRAL`/`UPARARRRP`).

`GWF2EVT8U1RP` (per stress period): `INSURF INEVTR INEXDP [INIEVT if NEVTOP=2]`
(plus an optional `INEVTZONES` keyword). A negative flag reuses the previous
array. `SURF`/`EVTR`/`EXDP` via `U2DREL`; `IEVT` via `U2DINT` — structured: a
**1-based layer index** validated in 1..NLAY; unstructured: a **1-based node
index** validated ≤ NODES.

## What was fixed / implemented (`flopy/mfusg/mfusgevt.py`)

1. **IETFACTOR round-trip bug** — `load` read `IETFACTOR` but never passed it to
   the constructor, so it reset to 0 on reload (and the `ETFACTOR` array would be
   dropped on the next write). Now preserved.
2. **`write_file` IETFACTOR** — written whenever transport is active (was only
   when `!= 0`), matching the Fortran's 3-integer dataset-1 read for `INBCT>0`.
   The `ETFACTOR` array is now written when `ietfactor > 0` (was `== 1`).
3. **Validation** — `NEVTOP` must be 1..3 (`ValueError`); `ETFACTOR` must have
   `MCOMP` values when `ietfactor>0` (`ValueError`); `NEVTOP=2` structured `IEVT`
   (0-based layer) must be in `[0, NLAY-1]` on write (`ValueError`).
4. **ETS zonal time-series** — `ETS MXZNEVT`/per-SP `IZNEVT` is **not** authored
   or parsed, so `mxetzones>0` now raises `NotImplementedError` rather than
   writing/loading an incomplete file.

0-based internal / 1-based file indexing for `NEVTOP=2 IEVT` (structured layer
and unstructured node) is preserved, and per-SP reuse (`-1`) works through
`Transient2d`.

## Tests

Synthetic (`autotest/test_usg_transport.py`, `-k mfusgevt`):

- `test_mfusgevt_nevtop1_and_3_authoring_roundtrip`
- `test_mfusgevt_nevtop2_structured_roundtrip` (0-based ↔ 1-based)
- `test_mfusgevt_nevtop2_unstructured_mxndevt`
- `test_mfusgevt_multi_period_reuse`
- `test_mfusgevt_transport_ietfactor_roundtrip` (0 / <0 / >0, MCOMP>1)
- `test_mfusgevt_rejects_invalid` (NEVTOP, ETS, ETFACTOR length, IEVT range)
- `test_mfusgevt_transport_etfactor_authoring` (existing).

Executable (`autotest/test_usg_transport_exe.py`, USG-T 2.7):

- `test_usgt_exe_evt_from_scratch` — a from-scratch `NEVTOP=1` EVT on the steady
  CHD flow model runs to normal termination, reports an `ET` budget term, and
  closes. (A cheap acceptance smoke; EVT-perturbed heads are not asserted
  analytically.)

`-k mfusgevt` **7 passed**; focused **144 passed**; exe **4 passed**; combined
**148 passed** under the USG-T 2.7 ARM binary.

## Decision: keep `✅ (intentionally not Full)`, hardened

Core EVT authoring is now broad and round-trip tested (all three `NEVTOP`
options, structured and unstructured `NEVTOP=2`, transport `IETFACTOR` 0/<0/>0
with per-`MCOMP` `ETFACTOR`, per-SP reuse), a real `IETFACTOR` round-trip bug was
fixed, and the executable accepts a from-scratch EVT. But two bounded gaps keep
EVT honest at `✅ (intentionally not Full)`:

- **ETS zonal time-series** (`ETS MXZNEVT` / per-SP `INEVTZONES`/`IZNEVT`):
  **deferred — Stage 4.6F-B** (executed). It is an ATS-coupled dynamic execution
  mode (an external time-series file superseding the EVTR array), not static I/O,
  so it is not modeled; it now fails explicitly per branch (authoring + load),
  with the full Fortran spec recorded in `USGT_STAGE4_16_EVT_ETS_ZONAL.md`.
- **NPEVT named parameters**: **resolved in Stage 4.6F-A** (below) — the EVTR
  array parameters now load → write → reload with their syntax intact and can be
  authored from scratch; `expand_parameters=True` keeps the legacy expanded
  `NPEVT=0` path.

The ETS-zonal gap (4.6F-B) keeps EVT honest at `✅ (intentionally not Full)`.

## Review follow-up (resolved)

1. **Scalar ETFACTOR no longer crashes (P1).** `MfUsgEvt(..., ietfactor=1,
   etfactor=2.5)` previously passed the length check (`np.atleast_1d(2.5).size
   == 1`) but crashed in `write_file` at `self.etfactor[icomp]`. `__init__` now
   normalizes `self.etfactor` to a 1-D array, so a scalar (MCOMP=1) and an array
   (MCOMP>1) both author/load/reload correctly; an incompatible length still
   raises `ValueError`.
2. **Unstructured IEVT node range validated (P2).** For unstructured `NEVTOP=2`,
   `IEVT` is a 0-based node index written 1-based. `write_file` now validates it
   against the total node count (via `node_count`, i.e. DISU `nodes` or the
   DIS fallback): a value `< 0` or `>= NODES` raises a clear `ValueError`.
   Structured `NEVTOP=2` continues to validate the layer index in `[0, NLAY-1]`.

Tests added: `test_mfusgevt_etfactor_scalar_and_array_roundtrip`,
`test_mfusgevt_nevtop2_unstructured_ievt_out_of_range`. `-k mfusgevt`
**9 passed**; focused **146 passed**; exe **4 passed**; combined **150 passed**
under the USG-T 2.7 ARM binary. EVT status unchanged: `✅ (intentionally not
Full)`.

## Stage 4.6F-A — NPEVT parameter preservation + from-scratch authoring

Closes the `NPEVT` gap above. EVT parameterizes only the **EVTR** (max ET-rate)
array, via the **same** MODFLOW array-parameter machinery as ETS.

**Fortran audit (`gwf2evt8u1.f` + `parutl7.f`):** item 1 is read with
`UPARARRAL(IN,IOUT,LINE,NPEVT)` — with `IN>0` it decodes an **optional**
`PARAMETER NPEVT` line that *precedes* item 2 (`NEVTOP IEVTCB [IETFACTOR]`); this
is the one structural difference from ETS, which carries `NPETS` in item 2a with
no `PARAMETER` line. When `NPEVT>0`, `NPEVT` definitions follow (`UPARARRRP`,
`PTYP='EVT'`, `ITVP=1` so `INSTANCES` are allowed). Per stress period, when
`INEVTR>=0` and `NPEVT>0`, `INEVTR` is the **count of active EVTR parameters**
(`UPARARRSUB2`, `'EVT'`); `INEVTR<0` reuses the previous period's EVTR. Only EVTR
is parameterized (`SURF`/`EXDP`/`IEVT` stay plain arrays). The `ETS MXZNEVT`
zonal sub-mode is independent of `NPEVT` (separate option, deferred to 4.6F-B).

**Implementation (`flopy/mfusg/mfusgevt.py`, mirrors `MfUsgEts`):**

- `load(..., expand_parameters=False)` (default) **preserves**: keeps `npevt`,
  the parsed `ModflowParBc` (`self.parameters`), and the per-SP activation
  records (`self.evtr_parm = {kper: [(name, instance_or_None), ...]}`).
  `expand_parameters=True` keeps the legacy expand-to-arrays / `NPEVT=0` path.
- `write_file` re-emits the `PARAMETER NPEVT` line, the definition blocks, and
  per-period `INEVTR = len(records)` (or `-1` to reuse) with the activation
  records in place of the EVTR array — `SURF`/`EXDP`/`IEVT` stay plain arrays.
- **From-scratch authoring:** pass `parameters={name: {parval, clusters |
  instances}}` + `evtr_parm={kper: [(name, instance_or_None), ...]}`. The
  ergonomic dict is built into a `ModflowParBc` via the shared
  `build_array_parameter_bc_parms` (`partyp='evt'`); `npevt` auto-computed; the
  first stress period must activate (USG-T cannot reuse an uninitialized EVTR).
  `_resolve_parameters` / `_validate_active_params` validate **before opening**
  the file (no partial file): `npevt>0`/`evtr_parm` with no defs → `ValueError`;
  per-SP duplicate/undefined/instance-mismatch → `ValueError`; duplicate
  definition names (case-insensitive) → `ValueError`.

Only `mfusgevt.py` was touched; the shared array-parameter helpers in
`_usgt_parameters.py` are reused unchanged, and `MfUsgEts` is not modified.

**Review follow-up (hardening):** malformed inputs now fail with a clear
`ValueError` before any file is opened (no `AttributeError`/tuple-unpack, no
partial file). `_check_evtr_parm` validates the `evtr_parm` structure (must be a
dict keyed by integer stress period; each value a list of `(name,
instance_or_None)` pairs — a string is not a valid record sequence; name a
non-empty string; instance a string or None), run at the top of
`_resolve_parameters` so even a non-dict `evtr_parm` with `parameters=None` is
caught before `.values()`. `parameters={}` now raises "parameters is empty"
directly instead of tripping the first-period check. On `load`, a `PARAMETER`
line without an integer count (`PARAMETER`, `PARAMETER abc`) raises a clear
`ValueError`; `PARAMETER 0` still loads as plain non-parametric.

**Tests:** `test_mfusgevt_npevt_authoring_and_roundtrip`,
`_npevt_instances_roundtrip`, `_npevt_expand_parameters`,
`_npevt_rejects_invalid`, `_npevt_validation_hardening`,
`_npevt_parameter0_loads_nonparametric`. The 9 Card-B EVT tests and the
`test_usgt_exe_evt_from_scratch` exe smoke stay green. `-k mfusgevt` **15**,
focused **302**, exe **4**, combined **306** (USG-T 2.7 ARM).

**Status:** EVT stays `✅ (intentionally not Full)` — the NPEVT gap is closed, but
the ETS-zonal time-series (4.6F-B) remains `NotImplementedError`.

## Stage 4.6F-B — ETS zonal time-series (deferred, explicit per-branch failure)

Fortran audit (`gwf2evt8u1.f`) confirms `ETS MXZNEVT` is a **dynamic, ATS-coupled
execution mode**, not static I/O: it **requires** adaptive time-stepping (the
Fortran `STOP`s when `IATS==0`), reads an **external time-series file** (`IUETS`)
record-by-record during the run (`Tstart Tend Factor Ets(1..MXZNEVT)`), and with
`IETSOPT=1` the time-series **supersedes** the EVTR array (EVTR is recomputed each
step as `etsevt(IZNEVT(n))*AREA*Factor`); a per-SP `INEVTZONES` flag redefines
`MXZNEVT` and (re)reads the `IZNEVT` zone-index array. **Decision: deferred** —
FloPy models static EVT I/O only.

The single generic `NotImplementedError` was replaced by specific, per-branch
explicit failures (no partial parse): `__init__` (`mxetzones>0`) cites the full
spec; `load` rejects `ETS MXZNEVT` on item 2 before reading any stress period;
`load` rejects a per-SP `INEVTZONES` header. Plain EVT and NPEVT are unaffected.
Test: `test_mfusgevt_ets_zonal_deferred`. `-k mfusgevt` **16**, focused **303**,
exe **4**, combined **307** (ARM). Full spec + rationale:
`USGT_STAGE4_16_EVT_ETS_ZONAL.md`.

**EVT stays `✅ (intentionally not Full)`**: static I/O is complete (NPEVT closed,
ETS-zonal deferred with the spec recorded), but the ATS-coupled ETS time-series
is intentionally unsupported, so promoting to `Full` would overclaim.

## Validation

```bash
python -m pytest autotest/test_usg_transport.py -k mfusgevt -q
python -m pytest autotest/test_usg_transport.py -q
python -m pytest autotest/test_usg_transport_exe.py -q
USGT_EXE=/Users/roh.sgm/Documents/GitHub_Projects/GW-Software-Compiled-via-Claude/gfortran/usgt_2.7/usgt_270_arm python -m pytest autotest/test_usg_transport.py autotest/test_usg_transport_exe.py -q
git diff --check
git status --short
```
