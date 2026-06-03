# Stage 4.6E — DPT A-W_ADSORBIM (immobile-domain air-water adsorption)

Date: 2026-06-03 · Base: `develop` @ `b6cd1505`

Goal: close the DPT immobile-domain air-water interface adsorption
(`A-W_ADSORBIM`) gap. Previously `MfUsgDpt.load` raised `NotImplementedError` at
the option line. This stage implements the **array-only** function-index
branches end-to-end (constructor → write → load → reload + validation) and keeps
every other branch as an explicit, Fortran-grounded `NotImplementedError`.

## Fortran audit (`gwt2dptu1.f`, `dpt2aw_adsorb.f`)

### Item 1a (option line, `gwt2dptu1.f:135-148`)

After the seven DPT integers (`IPAKCB IDPTCON IC_BNDIM_FLG IADSORBIM IDISPIM
IZODIM IFODIM`) and the text options, a `URWORD` loop scans for the keyword
`A-W_ADSORBIM`, then reads **two integers** on the same line: `IAREA_FNIM`
(air-water area function index) and `IKAWI_FNIM` (K_AWI partition function
index). Sets `IAW_ADSORBIM=1`.

### `AW_ADSORBIM1AL` (dimensions/constants, `dpt2aw_adsorb.f:3-91`)

Always allocates `AREA_AWIIM(NODES)`, `AK_AWIIM(NODES,MCOMP)`. Then:

- **Zone map** — *iff* `IAREA_FNIM==5 .OR. IKAWI_FNIM==4`: set `ITAB_AWIIM=1`,
  read `NAZONESIM NATABROWSIM` (one line), then `IAWIZONMAPIM(NODES)` via
  `U1DINT`.
- `IAREA_FNIM==2`: read scalar `ROG_SIGMAIM` (via `URWORD`). *(Note: the Fortran
  guard here is the mobile-domain `IAREA_FN.EQ.2`, a likely typo for
  `IAREA_FNIM` — one reason this branch is risky to model.)*
- `IKAWI_FNIM==3`: read scalar `SIGMA_RTIM` (via `URWORD`).
- `IAREA_FNIM==5`: read the area-vs-saturation table
  `AWI_AREA_TABIM(2,NATABROWSIM,NAZONESIM)` here (per zone × row, 2 values).

### `AW_ADSORBIM1RP1` (area arrays, `dpt2aw_adsorb.f:93-245`)

Read **after** the heat parameters, **before** the per-species loop
(`gwt2dptu1.f:343-345`). Per `IAREA_FNIM`:

| `IAREA_FNIM` | reads | meaning |
|---|---|---|
| 1 | `AWAMAXIM` (1 array) | `A = AMAX·(1−Sw)` |
| 4 | `AWAREA_X2IM`,`X1IM`,`X0IM` (3 arrays) | `A = X2·Sw² + X1·Sw + X0` |
| 2 | grain-diameter array → `AWAMAXIM = 3.9/d^1.2` | (1 array) |
| 3 | none (computes `AWAMAXIM = ROG_SIGMA·porosity`) | needs the `==2` scalar |
| 5 | none (table already read in AL) | tabular |

### `AW_ADSORBIM1RP2(ICOMP)` (Langmuir/K_AWI per species, `dpt2aw_adsorb.f:247-353`)

Called **inside the species loop, before `ADSORBIM`** (`gwt2dptu1.f:379-381`),
for each mobile species. If `IKAWI_FNIM != 4`: read `ALANGAWIM(:,ICOMP)` and
`BLANGAWIM(:,ICOMP)` (2 arrays). Then prepare per `IKAWI_FNIM`:

| `IKAWI_FNIM` | I/O | prep (solver-internal) |
|---|---|---|
| 1 | 2 arrays (A, B) | none |
| 2 | 2 arrays (Cmax, KL) | `A = A + B` |
| 3 | 2 arrays (A_aw, B_aw) | Brusseau (uses `SIGMA_RTIM`) |
| 4 | table `AWI_KAWI_TABIM(2,NATABROWSIM,NAZONESIM,ICOMP)` per zone × row | tabular |

So `IKAWI_FNIM` 1/2/3 have **identical file layout** (two arrays per species);
they differ only in a solver-internal transform. `==4` reads a table instead.

### Full per-model I/O order

option line (`A-W_ADSORBIM iarea ikawi`) → AL (zone map / scalars / area table,
only for the tabular/scalar branches) → DPT base arrays (PHIF, PRSITYIM, …,
DDTTR, heat) → **RP1 area arrays** → species loop: **RP2 (ALANG/BLANG)** →
ADSORBIM → … → CONCIM.

## Decision — implement the array-only subset

FloPy implements the branches that read **only plain arrays** (no zone map, no
scalar constants, no tabular functions), so authoring/load/write is exact and no
later read can be shifted:

- **`IAREA_FNIM ∈ {1, 4}`** — AMAX (1 array) / X2,X1,X0 (3 arrays).
- **`IKAWI_FNIM ∈ {1, 2}`** — Langmuir A/B (2 arrays per mobile species).

These cover the common non-tabular usage and require **no** AL reads.

### Deferred (explicit `NotImplementedError`, Fortran-grounded)

- `IAREA_FNIM == 2` (grain diameter; also has the `IAREA_FN` typo guard),
  `== 3` (`ROG_SIGMA·porosity`, needs the scalar), `== 5` (tabular area + zone
  map).
- `IKAWI_FNIM == 3` (Brusseau, needs the `SIGMA_RT` scalar), `== 4` (tabular
  K_AWI + zone map).

`IAREA_FNIM==5` and `IKAWI_FNIM==4` (the zone-map/tabular forms the prompt calls
out) raise a specific `NotImplementedError` naming the index, on **both**
authoring (`__init__`) and load — on load *before any array is read*, so the
rest of the file is never mis-parsed.

## Implementation (`flopy/mfusg/mfusgdpt.py`)

- `__init__` gains `aw_adsorbim`, `iarea_fnim`, `ikawi_fnim`, `awamaxim`,
  `awarea_x2im/x1im/x0im`, `alangawim`, `blangawim`. When `aw_adsorbim`: validate
  the branch (`_check_aw_adsorbim_supported`), require `mcomp>0` (it is a
  per-species Langmuir option), and build the `Util3d` area array(s) + one
  Langmuir A/B `Util3d` pair per species.
- `write_file`: append `A-W_ADSORBIM <iarea> <ikawi>` to the option line; write
  the RP1 area arrays after the heat block; write the RP2 ALANG/BLANG pair at the
  top of each species iteration, before `ADSORBIM` — matching the Fortran order.
- `load`: parse the two indices off the option line; reject unsupported branches
  before reading; read RP1/RP2 in the matching positions.

## Validations (actionable, before any write / before shifting reads)

| Condition | Result |
|---|---|
| `IAREA_FNIM ∉ {1,4}` (2/3/5) | `NotImplementedError` naming `IAREA_FNIM` |
| `IKAWI_FNIM ∉ {1,2}` (3/4) | `NotImplementedError` naming `IKAWI_FNIM` |
| `A-W_ADSORBIM` with `mcomp==0` | `ValueError` (per-species option needs transport) |
| `A-W_ADSORBIM` on load without the two indices | `ValueError` (malformed) |

Authoring failures occur in `__init__` (before `write_file`), so no partial
file is produced; load failures occur at the option line (before any array
read).

## Tests (`autotest/test_usg_transport.py`, `-k mfusgdpt`)

- `test_mfusgdpt_aw_adsorbim_authoring_roundtrip` — `IAREA_FNIM=1`,
  `IKAWI_FNIM=1`: writes `A-W_ADSORBIM 1 1`, reload preserves the indices and the
  `AWAMAXIM`/`ALANGAWIM`/`BLANGAWIM` arrays.
- `test_mfusgdpt_aw_adsorbim_iarea4_roundtrip` — `IAREA_FNIM=4`, `IKAWI_FNIM=2`,
  two species: X2/X1/X0 + per-species Langmuir A/B round-trip.
- `test_mfusgdpt_aw_adsorbim_rejects_unsupported` — `IAREA_FNIM ∈ {2,3,5}` and
  `IKAWI_FNIM ∈ {3,4}` (authoring → `NotImplementedError`); `mcomp=0` →
  `ValueError`; load of `5 4` (raises on `IAREA_FNIM=5`) and `1 4` (raises on
  `IKAWI_FNIM=4`) before any array read; a bare keyword → `ValueError`.

(The old `test_mfusgdpt_aw_adsorbim_fails_explicitly` is superseded by these.)

No `A-W_ADSORBIM` executable smoke was added: a convergent dual-porosity
air-water adsorption model is not cheap to build, and the writer is audited
line-by-line against the Fortran and round-trips in FloPy.

## Status after Stage 4.6E

DPT stays **⚠️ Partial** — but `A-W_ADSORBIM` is upgraded from *explicit failure*
to **array-only authoring + load + round-trip** (`IAREA_FNIM ∈ {1,4}`,
`IKAWI_FNIM ∈ {1,2}`). Remaining (kept `NotImplementedError`): the scalar/tabular
branches `IAREA_FNIM ∈ {2,3,5}` and `IKAWI_FNIM ∈ {3,4}` (these add a zone map,
`ROG_SIGMA`/`SIGMA_RT` scalars, and area / K_AWI tables), and from-scratch
parametric *execution* is not USG-T smoke-tested.

## Validation

```bash
python -m pytest autotest/test_usg_transport.py -k mfusgdpt -q   # 3 passed
python -m pytest autotest/test_usg_transport.py -q               # 294 passed
python -m pytest autotest/test_usg_transport_exe.py -q           # 4 passed
USGT_EXE=.../usgt_270_arm python -m pytest autotest/test_usg_transport.py autotest/test_usg_transport_exe.py -q   # 298 passed
git diff --check
```
