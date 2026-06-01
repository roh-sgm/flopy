# Stage 4.4B — MODFLOW Parameter Preservation: HFB (executed)

Scope: implement load → write → reload **preservation** of HFB named parameters,
audited against `gwf2hfb7u1.f` and `parutl7.f`. Only HFB and the shared helper
are touched; SGB/DRT/QRT are out of scope (they keep their explicit
`NotImplementedError` on `NP*>0`).

## Fortran audit

### `gwf2hfb7u1.f`

- **Item 1 (AR, `GWF2HFB7U1AR`)**: `NPHFB MXFBP NHFBNP [NOPRINT] [TRANSIENT_HFB]`.
  `NPHFB` = number of barriers defined by parameters, `MXFBP` = max barriers from
  parameters, `NHFBNP` = barriers not defined by parameters. No parameters are
  read in AR.
- **Per stress period (RP, `GWF2HFB7U1RP`)**:
  - `IHFBRD` controls whether the block is (re)read. For `TRANSIENT_HFB`
    (`ITRHFB>0`) it is read from the file each period; otherwise it is set
    internally (`1` for KPER=1, `0` after). `IHFBRD<=0` → reuse previous period
    (read nothing). So a **non-transient** file has **no `IHFBRD` line**.
  - **Items 2-3** (`IF NPHFB>0`): `NPHFB` list-parameter definitions via
    `UPARLSTRP(...,'HFB ','HFB ',1,NUMINST)`. **`NUMINST>0` is rejected** — HFB
    does not support `INSTANCES` (`gwf2hfb7u1.f:135-138`). Each definition is a
    header `PARNAM PARTYP PARVAL NLST` followed by `NLST` barrier rows
    (`SGWF2HFB7RL`/`RLU`).
  - **Item 4** (`IF NHFBNP>0`): `NHFBNP` non-parametric barrier rows.
  - **Items 5-6** (`IF NPHFB>0`): `NACTHFB`, then `NACTHFB` activation lines, each
    a parameter name (`SGWF2HFB7SUB`). Activation scales the barrier `HYDCHR` by
    the parameter value `B(IP)` (`gwf2hfb7u1.f:691`).
  - A barrier row is `LAYER IROW1 ICOL1 IROW2 ICOL2 FACTOR` (structured) or
    `NODE1 NODE2 FACTOR` (unstructured).

### `parutl7.f` (list parameters)

HFB uses the **list**-parameter machinery (`UPARLSTRP` / `UPARLSTSUB`), not the
array-parameter machinery used by ETS. `UPARLSTRP` is called with `ITERP=1`, so
re-reading the same parameter definitions in a later stress period would trip
the "Duplicate parameter name" guard (`parutl7.f:604`). This is why
`TRANSIENT_HFB` + `NPHFB>0` is a problematic combination and is not supported
here (see Decision).

## Design — shared helper (list-parameter path)

`flopy/mfusg/_usgt_parameters.py` gains the **list-parameter** write/read path
alongside the existing array-parameter (ETS) path:

- `read_list_parameter_header(line)` → `(name, partyp, parval, nlst, numinst)`
  (`UPARLSTRP` header; `parval` kept as the original string).
- `write_list_parameter_header(f, name, partyp, parval, nlst)`.
- `read_active_list_parameters(f, count)` → list of names (`UPARLSTSUB`; HFB
  activations carry only the name).
- `write_active_list_parameters(f, names)`.

The `NLST` barrier rows are package-specific, so the row reader/writer stays in
HFB (`_read_hfb_rows` / `_write_hfb_rows`, which already convert 0-based↔1-based).
These helpers are reusable by SGB/DRT/QRT when their preservation lands.

## HFB implementation (`MfUsgHfb`)

- `load`: when `NPHFB>0` (non-transient), read the `NPHFB` definitions
  (header + `NLST` rows, **0-based internally**) into
  `self.parameters = {name: {"partyp", "parval", "nlst", "data"}}`, then the
  `NHFBNP` non-parametric rows (`self.hfb_data`), then `NACTHFB` and the active
  names (`self.acthfb_names`). `MXFBP`→`self.mxfb` and `NPHFB`→`self.nphfb` are
  preserved. `INSTANCES` (`NUMINST>0`) → `NotImplementedError`.
- `write_file`: when `NPHFB>0` with loaded definitions, write item 1 with
  `NPHFB`, the definition blocks, the non-parametric barriers, then `NACTHFB` and
  the active names (1-based barrier rows). The parameterized barriers and the
  non-parametric barriers are written together — the mix the Fortran allows.
- Non-parametric HFB (static and `TRANSIENT_HFB` with `NPHFB=0`,
  `IHFBRD = -1/0/>0`) is unchanged.
- Indexing: internal 0-based, file 1-based, for both parameter and
  non-parametric barrier rows.

## Decision: HFB stays ⚠️ Partial (not `Full`)

Upgraded from "`NPHFB>0` fails explicitly" to **parameter-preserving
(load → write → reload)** for HFB list parameters. Still **not `Full`** because:

- **From-scratch parameter authoring** (`NPHFB>0` with no loaded definitions)
  is unsupported → `write_file` raises `NotImplementedError`.
- **`TRANSIENT_HFB` + `NPHFB>0`** is unsupported (the Fortran would
  re-read/redefine the parameters each stress period under `ITERP=1`) → both
  `load` and `write_file` raise `NotImplementedError`.
- `INSTANCES` are unsupported (matches the Fortran, which aborts).

`SFAC` / `OPEN/CLOSE` / `EXTERNAL` list controls inside barrier lists are
supported as of the list-control follow-up (below).

## Tests (Stage 4.4B)

New (Stage 4.4B):

- `test_mfusghfb_parameterized_roundtrip` — unstructured; preserves defs +
  `NACTHFB` + activation names; 0-based internal / 1-based file.
- `test_mfusghfb_parameterized_structured_roundtrip` — structured; verifies
  `k/i/j` 0-based internal and 1-based on file.
- `test_mfusghfb_parameterized_with_nonparam` — parameter-defined barriers mixed
  with non-parametric ones (item 2-3 + item 4).
- `test_mfusghfb_parameter_authoring_from_scratch_fails` — `NPHFB>0` without
  defs → `NotImplementedError`.
- `test_mfusghfb_transient_with_parameters_fails` — `TRANSIENT_HFB` + `NPHFB>0`
  → `NotImplementedError`.

Unchanged (non-parametric regression): `test_mfusghfb_static_fortran_layout`,
`test_mfusghfb_structured_static_fortran_layout`,
`test_mfusghfb_transient_fortran_layout` (covers `IHFBRD = 1` read and `-1`/`0`
reuse).

The replaced `test_mfusghfb_parameterized_fails_explicitly` is superseded: HFB no
longer fails on parameterized **load**; it preserves. The explicit-failure
contract now lives in the two negative tests above (authoring / transient+params).

No from-scratch parametric **executable** smoke was added (a convergent
parameterized HFB model is not cheap); the writer is audited line-by-line against
the Fortran and round-trips in FloPy. Executable validation remains a documented
manual tier.

## List-control follow-up (executed)

Re-audit of `SGWF2HFB7RL` / `SGWF2HFB7RLU` showed that **every** HFB barrier list
(both the per-parameter `NLST` rows in item 2-3 and the non-parametric `NHFBNP`
rows in item 4) may begin with the standard MODFLOW list controls before the
data rows: `SFAC <factor>`, `OPEN/CLOSE <fname>`, `EXTERNAL <unit>`. The previous
loader assumed inline numeric rows and raised `ValueError` on a valid file that
started a block with `SFAC`.

`MfUsgHfb._read_hfb_rows` now consumes these via the shared
`flopy/mfusg/_usgt_list.py::begin_list_block` (the same helper SGB/QRT/DRT use —
**reused unchanged**):

- `SFAC` scales the barrier `HYDCHR`/`FACTOR` (`HFB(6,II) = FACTOR*SFAC`), applied
  on read to each block independently (a parameter block and the non-parametric
  block can each carry their own `SFAC`).
- `OPEN/CLOSE` reads the rows from a named file, resolved against
  `model.model_ws`; quoted filenames (with spaces) are supported.
- `EXTERNAL` reads the rows from a NAM-declared unit resolved via
  `ext_unit_dict`; an unresolvable unit raises `NotImplementedError` (not
  `ValueError`).

The writer is unchanged: it emits expanded inline rows with the scaled values
("Expanded valid write") — `SFAC`/`OPEN/CLOSE`/`EXTERNAL` are not reproduced.
Round-trip is semantically exact: for a parameter, `FACTOR*SFAC` is baked into
the written row and `PARVAL` is preserved separately, so the effective
`FACTOR*SFAC*PARVAL` is unchanged.

`write_file` also gained a guard: for a parameterized HFB, `nacthfb` must equal
`len(acthfb_names)`, otherwise a clear `ValueError` is raised (no partial write).

New tests (Stage 4.4B follow-up): `test_mfusghfb_sfac_scales_nonparam`,
`test_mfusghfb_sfac_parameterized_roundtrip`,
`test_mfusghfb_sfac_mixed_param_and_nonparam`, `test_mfusghfb_open_close`
(plain + quoted), `test_mfusghfb_external_via_ext_unit_dict`,
`test_mfusghfb_external_without_dict_fails`,
`test_mfusghfb_nacthfb_mismatch_fails`. The SGB/QRT/DRT `_usgt_list` regressions
stay green (the helper was reused, not modified).

## Validation

```bash
python -m pytest autotest/test_usg_transport.py -k mfusghfb -q
python -m pytest autotest/test_usg_transport.py -q
python -m pytest autotest/test_usg_transport_exe.py -q
USGT_EXE=/Users/roh.sgm/Documents/GitHub_Projects/GW-Software-Compiled-via-Claude/gfortran/usgt_2.7/usgt_270_arm python -m pytest autotest/test_usg_transport.py autotest/test_usg_transport_exe.py -q
git diff --check
git status --short
```

Results (after the list-control follow-up): `-k "mfusghfb or usgt_list"`
**18 passed**; `-k mfusghfb` **15 passed**; focused **180 passed**; exe
**4 passed**; combined **184 passed** under the USG-T 2.7 ARM binary.
