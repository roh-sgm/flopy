# Stage 4.5A — QRT `TRANSIENTQ` support

Goal: replace `MfUsgQrt`'s explicit `NotImplementedError` for the `TRANSIENTQ`
option with real semantic support (authoring from scratch + load → write →
reload), if the Fortran audit confirms a stable format. It does — for the
**inline** form — so this stage promotes `TRANSIENTQ` to a preserved,
authorable feature, while keeping the genuinely fragile combination
(`TRANSIENTQ` + `NPQRT>0`) an explicit failure.

## Fortran audit (`gwf2QRT8u.f`, USG-T 2.7)

All line numbers refer to `USGT_V_2-7-0_Source_Code/gwf2QRT8u.f`.

### Item 1 option `TRANSIENTQ NBDQTIM` (`GWF2QRT8U1AR`, lines 129–146)

```fortran
ELSEIF(LINE(ISTART:ISTOP).EQ.'TRANSIENTQ') THEN
  CALL URWORD(LINE,LLOC,ISTART,ISTOP,2,NBDQTIM,R,IOUT,INOC)   ! signed int
  IF(NBDQTIM.LT.0) THEN
    ISTEPQ = 1                 ! staircase (no interpolation)
    NBDQTIM = -NBDQTIM         ! magnitude = number of time points
    ALLOCATE(TIMQRT); TIMQRT=0
  ENDIF
  ALLOCATE (BDQTIM(NBDQTIM), BDQV(NBDQTIM,MXAQRT))
  ALLOCATE (IQRTN(MXAQRT))
ENDIF
```

- `NBDQTIM` is read with `URWORD` type 2 (integer), **signed**:
  - `NBDQTIM < 0` → `ISTEPQ = 1` (**staircase**: the value at the lower-bound
    time index is held, no interpolation); `NBDQTIM` is then set to its
    magnitude.
  - `NBDQTIM > 0` → `ISTEPQ = 0` (**interpolation**, the default).
- Arrays allocated: `BDQTIM(NBDQTIM)` (the times), `BDQV(NBDQTIM, MXAQRT)` (a
  flow value per time per active sink), `IQRTN(MXAQRT)` (a node tag per sink).
- **Ordering constraint:** the `TRANSIENTQ` branch is the only option branch
  with **no `GOTO 10`** back to the option loop. Once `TRANSIENTQ` is parsed the
  option scan ends, so any option written *after* `TRANSIENTQ` on item 1 is
  silently ignored. **`TRANSIENTQ` must be the last option on item 1.** FloPy
  writes it last.

### Data block (`GWF2QRT8U1AR`, lines 207–246) — read once, in AR

The block is read inside `GWF2QRT8U1AR`, **after** the `NPQRT` named-parameter
block (lines 175–206) and **before** any stress period (`GWF2QRT8U1RP`):

```fortran
IF(NBDQTIM.NE.0)THEN
  READ(IN,*) IQRTUN,CNSTM                         ! times control line
  READ(IQRTUN,*) (BDQTIM(L),L=1,NBDQTIM)          ! the NBDQTIM times
  DO L=1,NBDQTIM;  BDQTIM(L)=BDQTIM(L)*CNSTM;  ENDDO
  READ(IN,*) IQRTUN,CNSTM                         ! values control line
  DO N = 1,MXAQRT                                 ! EXACTLY MXAQRT rows
    READ(IQRTUN,*) IQRTN(N),(BDQV(L,N),L=1,NBDQTIM)
    DO L=1,NBDQTIM;  BDQV(L,N)=BDQV(L,N)*CNSTM;  ENDDO
  ENDDO
ENDIF
```

- Two control lines, each `IQRTUN CNSTM`: a Fortran unit and a multiplier.
  `BDQTIM` and `BDQV` are scaled by their respective `CNSTM` on read.
- **Inline vs external:** `READ(IQRTUN,*)` reads from unit `IQRTUN`. When
  `IQRTUN` equals the unit the NAM connected the QRT file to, the data follow
  **inline** in the QRT file; otherwise they live in a separate file on that
  unit. FloPy supports the **inline** form only (it writes `IQRTUN` = the QRT
  package's own unit number; on load it reads the data from the same stream).
  External-unit references are not resolved (documented limitation).
- **Exactly `MXAQRT` value-rows** are read (`DO N=1,MXAQRT`), regardless of how
  many sinks are active in any stress period. Each row is
  `IQRTN(N)  BDQV(1..NBDQTIM)`.

### What varies in time, and the row→sink mapping (`GWF2QRT8U1AD`, lines 1306–1398)

```fortran
DO 90 NF=1,MXQRT
  QA=BDQV(IB1,NF); QB=BDQV(IB2,NF); ...            ! interpolate over time
  QRTF(4,NF)= ...                                  ! <-- writes Q only
```

- **Only `QRTF(4,NF) = Q` is updated.** `QRTF(5)=NumRT` (recipient count) and
  `QRTF(6)=Rfprop` (return proportion) are **never touched** by the time series.
  → audit Q3 answered: **`TRANSIENTQ` varies `Q`, not recipients.**
- **The mapping is positional.** `BDQV(L,NF)` drives QRT row `NF` directly by
  position (`NF=1..MXQRT`). `IQRTN(N)` is read (line 226) and printed (line 242)
  and imported into `AD` (line 1315) but **never used** in the computation — it
  is an *informational node tag*, not a key. So value-row `N` applies to the
  `N`-th QRT sink (in `QRTF`/item order), and `IQRTN(N)` is preserved as-is.
- `NBDQTIM == 1` is a special case (lines 1324–1331): `QRTF(4,NF)=BDQV(1,NF)`
  (a constant override). `NBDQTIM >= 2` interpolates (or staircases) between
  bracketing times. Staircase (`ISTEPQ=1`) sets the interpolation factors to 0
  (lines 1352, 1369).

### Interactions

- **`ITMP` / per-SP rows (`GWF2QRT8U1RP`, lines 274–322):** unchanged by
  `TRANSIENTQ`. Each stress period still reads `ITMP` and its sink rows (node,
  `NumRT`, `Rfprop`, recipient `U1DINT` blocks). `AD` overwrites `QRTF(4)=Q`
  every time step, so the `Q` column in the per-SP rows is effectively ignored
  while `TRANSIENTQ` is active. FloPy still writes the per-SP blocks normally;
  the `TRANSIENTQ` block sits between item 1 and the first stress period.
- **`NPQRT > 0` (parameters):** the Fortran *allows* the combination, but it is
  **fragile**: `BDQV` is allocated `(NBDQTIM, MXAQRT)` (line 144) while `AD`
  loops `NF=1..MXQRT` where `MXQRT = MXAQRT + MXL` (line 153). With `MXL>0`
  (i.e. `NPQRT>0`) the loop reads `BDQV(L,NF)` past its `MXAQRT` columns —
  out-of-bounds for the parametric rows. Per the task's instruction #2, FloPy
  **rejects `TRANSIENTQ` + `NPQRT>0`** with `NotImplementedError` (on both load
  and write), rather than emit/parse a fragile file.
- **`RETURNFLOW`:** compatible. The time-varying `Q` is partitioned onto the
  recipient nodes each step (`AD` lines 1406–1421) using the unchanged
  `NumRT`/`Rfprop`. `TRANSIENTQ` + `RETURNFLOW` round-trips.
- **`AUTOFLOWREDUCE` / `IUNIT_AFR_QRT`:** independent header options; no
  conflict.
- **`SFAC`:** scales `Q` during the per-SP row read (`ISCLOC=4`), but that `Q`
  is overwritten by the time series in `AD`, so `SFAC` on `Q` is moot under
  `TRANSIENTQ`. Not a blocker.

## FloPy design (`MfUsgQrt`)

New attributes (all default to the "no TRANSIENTQ" state):

| attribute | meaning | indexing |
|---|---|---|
| `transientq_times` | `BDQTIM` (length `NBDQTIM`), **raw** (pre-`CNSTM`) | — |
| `transientq_values` | `BDQV` as shape `(MXAQRT, NBDQTIM)` = (sink, time), **raw** | — |
| `transientq_nodes` | `IQRTN` (length `MXAQRT`), informational node tag | **0-based internal / 1-based file** |
| `transientq_staircase` | `ISTEPQ` (the file's `NBDQTIM < 0`) | — |
| `transientq_times_mult` | `CNSTM` for the times line (default `1.0`) | — |
| `transientq_values_mult` | `CNSTM` for the values line (default `1.0`) | — |

`TRANSIENTQ` is **active** iff `transientq_times is not None`. It is *not*
expressed through the `options` string list — it carries structured data, so it
has dedicated attributes; on write FloPy derives the `TRANSIENTQ <±NBDQTIM>`
token from the arrays and appends it **last** on item 1.

- **Raw values + multipliers are preserved verbatim** (the block round-trips
  textually, including a non-unit `CNSTM`). For from-scratch authoring the
  multipliers default to `1.0`, so `transientq_times`/`transientq_values` are
  the effective values.
- **0-based internal / 1-based file** for `transientq_nodes` (`+1` on write,
  `-1` on load), matching every other QRT node field.
- The writer emits `IQRTUN` = `self.unit_number[0]` on both control lines, so
  the data are inline and the file is runnable.

### Explicit failures (no partial write)

On **write**, validated **before the file is opened**:

- `TRANSIENTQ` + `NPQRT>0` (loaded params or any active params) →
  `NotImplementedError` (the fragile Fortran combination above).
- `TRANSIENTQ` active with `MXAQRT == 0` (no active sinks) → `ValueError`.
- `len(transientq_nodes) != MXAQRT` or `transientq_values.shape != (MXAQRT,
  NBDQTIM)` → `ValueError` (the Fortran reads exactly `MXAQRT` rows of
  `NBDQTIM` values).
- A `transientq_nodes` entry that is not a non-negative 0-based integer →
  `ValueError` (the file is 1-based, written as `node+1`). The upper bound is
  not range-checked against a node count, consistent with `recipient_nodes`
  (USG node numbering is not tied to a synthetic model's structured shape).
- `NBDQTIM < 1` → `ValueError`.

On **load** (Stage 4.5A review follow-up closes two contract gaps):

- **External-unit `TRANSIENTQ` data** — a control-line `IQRTUN` that is not the
  QRT package's own (inline) unit → `NotImplementedError`, raised *before* any
  data is interpreted as inline (checked on **both** the times and values
  control lines). The inline unit is resolved from the NAM via
  `model.get_ext_dict_attr(..., pop_key=False)` when an `ext_unit_dict` is
  available, else `MfUsgQrt._defaultunit()` — matching what `write_file` emits
  (`self.unit_number[0]`).
- **`TRANSIENTQ` not last on item 1** — any trailing token after
  `TRANSIENTQ <±NBDQTIM>` → `ValueError` (the Fortran branch does not loop back,
  so USG-T would silently ignore it; FloPy rejects rather than reorder).
- `TRANSIENTQ` + `NPQRT>0` → `NotImplementedError`; an unexpected EOF inside the
  block → `ValueError`.

## Tests (`autotest/test_usg_transport.py`)

- `test_mfusgqrt_transientq_authoring` — from-scratch QRT with a minimal
  `TRANSIENTQ` (interpolation, `NBDQTIM=2`); inspects the written text for the
  `TRANSIENTQ 2` item-1 token, inline control lines, and 1-based node tags.
- `test_mfusgqrt_transientq_roundtrip` — load → write → reload preserves the
  block (times, values, nodes, staircase flag, multipliers).
- `test_mfusgqrt_transientq_staircase` — `NBDQTIM<0` round-trips with
  `transientq_staircase=True` and a `-N` token on item 1.
- `test_mfusgqrt_transientq_nodes_1based` — internal nodes are 0-based, the
  file is 1-based (text inspection + reload equality).
- `test_mfusgqrt_transientq_with_params_fails` — `TRANSIENTQ` + `NPQRT>0` raises
  `NotImplementedError` with no partial file.
- `test_mfusgqrt_transientq_dim_mismatch_fails` — wrong `transientq_values`
  shape raises `ValueError` with no partial file.

Review follow-up tests (contract guards):

- `test_mfusgqrt_transientq_external_times_unit_fails` /
  `test_mfusgqrt_transientq_external_values_unit_fails` — an external `IQRTUN`
  on the times / values control line raises `NotImplementedError`.
- `test_mfusgqrt_transientq_trailing_option_fails` — a token after
  `TRANSIENTQ <N>` (e.g. `AUX C01` or `NOPRINT`) raises `ValueError`.
- `test_mfusgqrt_transientq_negative_node_fails` — a negative `transientq_nodes`
  entry raises `ValueError` with no partial file.
- `test_mfusgqrt_transientq_multipliers_roundtrip` — non-unit `CNSTM`
  multipliers survive load → write → reload (stored times/values stay raw).

The existing QRT suite (parameters, return flow, CHANGEC, list controls) stays
green: `-k mfusgqrt` **31 passed**; focused **229**; exe **4**; combined
**233** under the USG-T 2.7 ARM binary.

## Validation

```bash
python -m pytest autotest/test_usg_transport.py -k mfusgqrt -q
python -m pytest autotest/test_usg_transport.py -q
python -m pytest autotest/test_usg_transport_exe.py -q
USGT_EXE="/Users/roh.sgm/Documents/GitHub_Projects/GW-Software-Compiled-via-Claude/gfortran/usgt_2.7/usgt_270_arm" python -m pytest autotest/test_usg_transport.py autotest/test_usg_transport_exe.py -q
git diff --check
git status --short
```

## Final QRT status

QRT is now **structural-preserving (Stage 4.4E) + `TRANSIENTQ` supported (Stage
4.5A)**: non-parametric sinks with return flow / CHANGEC / list controls are
full-fidelity; named parameters (`NPQRT>0`) are structurally preserved but not
execution-guaranteed (the parameter value scales `NumRT`, not `Q`, and `NodQRT`
is not copied on activation — both Fortran issues, see
`USGT_STAGE4_04_PARAMETERS_QRT.md`); and the **inline `TRANSIENTQ`** transient
extraction-flow time series now loads → writes → reloads and can be authored
from scratch. Unsupported (explicit failure): `TRANSIENTQ` + `NPQRT>0`,
external-unit `TRANSIENTQ` data, parameter `INSTANCES`, and from-scratch
parameter authoring.
