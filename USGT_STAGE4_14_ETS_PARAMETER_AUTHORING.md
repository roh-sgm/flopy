# Stage 4.6C-D — ETS from-scratch *array*-parameter authoring (`NPETS>0`)

Date: 2026-06-02 · Base: `develop` @ `f975fe32`

Goal: let a user build a parameterized `MfUsgEts` from Python (`NPETS>0`) for the
ETSR (max ET-rate) array, write a valid USG-T 2.7 file, and reload it with the
same semantics. ETS uses MODFLOW **array** parameters (`ModflowParBc` /
`UPARARRRP`), **not** the list-parameter family (DRT/HFB/SGB/QRT) — so this stage
deliberately does *not* reuse the list-parameter `_normalize_param` pattern; it
builds a `ModflowParBc` from an ergonomic dict instead. This is the last of the
parameter-authoring cards.

## Audit

The on-file format is identical whether the definitions came from `load` or the
user — only the *source* of `self.parameters` differs. The write side already
used the array-parameter writers (`write_array_parameter_defs` → `UPARARRRP`;
`write_active_array_parameters` → `UPARARRSUB2` per-period `PNAME [INSTANCE]`
records) and the preservation path round-trips (Stage 4.4A). So this is a
FloPy-side ergonomics + validation change in `mfusgets.py` (+ a shared builder).

Fortran (`gwf2ets8u1.f`): USG-T reads `NPETS` from **item 2a** and calls
`UPARARRAL` with `IN=-1` (`parutl7.f`), so — unlike MODFLOW-2005 ETS — there is
**no** separate `PARAMETER N` line. Definitions follow item 2a / the optional
`ESFACTOR` record. Per stress period, `INETSR` is the count of active ETSR
parameter records (`UPARARRSUB2`), or `< 0` to reuse the previous period's ETSR.
Only **ETSR** can be parameterized; ETSS/ETSX/IETS/PXDP/PETM are always plain
arrays, so a parameterized ETS naturally mixes a parameterized ETSR with
non-parameterized arrays. The **first** parametric period must activate (an
`INETSR=-1` first period would reuse an uninitialized ETSR).

`ModflowParBc.bc_parms` layout (from `loadarray`):
`bc_parms[name.lower()] = [{partyp, parval, nclu, timevarying}, {instnam.lower():
[[mltarr, zonarr, [izone, ...]], ...nclu clusters]}]`; a static parameter uses the
single instance key `"static"`.

## What changed (code)

- **`flopy/mfusg/_usgt_parameters.py`**: new `build_array_parameter_bc_parms(
  parameters, partyp, prefix)` (+ `_build_array_clusters`) — builds a validated
  `bc_parms` dict from an ergonomic mapping, reusing the shared
  `check_parameter_name` / `check_parval`. Validates `partyp`, names (single
  token ≤10, unique case-insensitively), instance names, per-instance `nclu`
  consistency, non-empty clusters, and positive-integer zones (`ALL` → empty
  zone list).
- **`flopy/mfusg/mfusgets.py`**: new `_resolve_parameters` (accepts an ergonomic
  `dict` → builds `ModflowParBc`, or a preserved `ModflowParBc`; auto-computes
  `NPETS`) and `_validate_active_params` (first period must activate; names
  defined case-insensitively; no per-period duplicates; time-varying params must
  name an existing instance; static params must not name a non-static instance;
  `kper` in range). `write_file` resolves+validates **before opening** the file,
  writes the auto `NPETS` in item 2a, and the old from-scratch
  `NotImplementedError` is replaced (npets/active without defs → `ValueError`).

## API

```python
# static parameter:
ets = MfUsgEts(
    mfusg_model, netsop=1,
    parameters={"etsrate": {"partyp": "ets", "parval": "5e-4",
                            "clusters": [("NONE", "ALL", [])]}},
    evtr_parm={0: [("etsrate", None)]},   # activate per stress period
)
# time-varying (INSTANCES):
parameters={"etsrate": {"parval": "5e-4",
                        "instances": {"spring": [("NONE", "ALL", [])],
                                      "fall":   [("NONE", "ALL", [])]}}}
evtr_parm={0: [("etsrate", "spring")]}
```

A cluster is `(MLTARR, ZONARR[, zones])`; `ZONARR="ALL"` means the whole grid
(empty zone list). `partyp` defaults to `"ets"`; `nclu` is computed from the
clusters; `NPETS = len(parameters)` when omitted.

## Validations (all before the file is opened — no partial file)

| Condition | Result |
|---|---|
| `evtr_parm` / `npets>0` with no `parameters` | `ValueError` |
| first stress period does not activate (`evtr_parm[0]` missing) | `ValueError` |
| `partyp` not `ets` | `ValueError` |
| name not a single token ≤10 chars; duplicate def name (case-insensitive) | `ValueError` |
| `parval` missing/blank or multi-token | `ValueError` |
| empty clusters; non-positive-integer zone | `ValueError` |
| `nclu` given ≠ clusters; instances with differing `nclu` | `ValueError` |
| active name undefined / duplicate per period; `kper` out of range | `ValueError` |
| time-varying param activated without / with unknown instance | `ValueError` |
| static param activated with a non-static instance | `ValueError` |
| explicit `npets` ≠ `len(parameters)` | `ValueError` |

## Tests (`autotest/test_usg_transport.py`)

`test_mfusgets_parameter_authoring_from_scratch_static`, `_instances`,
`_netseg2_mixed_plain_arrays`, `_netsop2`, `_auto_npets`, and
`_negatives` (13 validations); `test_mfusgets_parameterized_write_npets_without_
defs_fails` (repurposed from `_write_fails_explicitly`, now `ValueError`). The
Stage 4.4A preservation tests (`_load_preserves`, `_write_preserves_syntax`,
`_roundtrip`, `_instances_roundtrip`, `_load_expands_to_npets0`) stay green.

Results: `-k mfusgets` **16**, `-k "mfusgets or mfusgqrt or mfusgdrt or
mfusgsgb"` **126**, focused **288**, exe **4**, combined **292** (USG-T 2.7 ARM).

## Status after Stage 4.6C-D

ETS: **ETSR array-parameter preserving + from-scratch authoring**. Still not
`Full`: only ETSR is parameterized (the Fortran limit); from-scratch parametric
*execution* is not USG-T smoke-tested here; the opt-in `expand_parameters=True`
load path (writes `NPETS=0`) is unchanged. **This completes the parameter
from-scratch authoring family** (ETS array + DRT/HFB/SGB/QRT list).

## Validation

```bash
python -m pytest autotest/test_usg_transport.py -k mfusgets -q   # 16 passed
python -m pytest autotest/test_usg_transport.py -q               # 288 passed
python -m pytest autotest/test_usg_transport_exe.py -q           # 4 passed
USGT_EXE=.../usgt_270_arm python -m pytest autotest/test_usg_transport.py autotest/test_usg_transport_exe.py -q   # 292 passed
git diff --check
```
