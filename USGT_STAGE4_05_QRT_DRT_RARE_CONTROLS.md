# Stage 4.5 - QRT/DRT Rare Controls

Goal: close the remaining rare-control gaps in QRT and DRT: QRT `TRANSIENTQ`
and recipient-node `U1DINT` `EXTERNAL`/`OPEN/CLOSE` controls.

## Current State

QRT and DRT are strong for common authoring:

- node-based records,
- return-flow recipients,
- AUX concentrations,
- `CHANGEC`,
- `ITMP/-1` reuse,
- main-list `SFAC`/`OPEN-CLOSE`/`EXTERNAL` load with expanded valid write.

Remaining gaps — **both now closed**:

- ~~`MfUsgQrt` rejects `TRANSIENTQ`.~~ **Done (Stage 4.5A)** — the inline
  `TRANSIENTQ` time series loads/writes/reloads and is authorable from scratch;
  see `USGT_STAGE4_05_QRT_TRANSIENTQ.md`.
- ~~`_usgt_returnflow` rejects recipient `EXTERNAL` / `OPEN/CLOSE`.~~ **Done
  (Stage 4.5B)** — see Outcome below.
- Parameter preservation is handled separately in
  `USGT_STAGE4_04_PARAMETERS.md`.

## Outcome — Stage 4.5B (recipient U1DINT EXTERNAL / OPEN-CLOSE)

**Audit:** both packages read their recipient-node lists with `U1DINT`
(`gwf2QRT8u.f:1057` `CALL U1DINT(NodQRT(IRT),...)`; `gwf2drt8u.f:833`/`:1005`
`CALL U1DINT(NodDRT(IRTSTRT),...)`). `U1DINT` (utl7u1.f) reads a control record
that may be `CONSTANT` (`LOCAT=0`, all = ICNSTNT), `INTERNAL` (`LOCAT=IN`,
inline), `EXTERNAL <unit>` (`LOCAT=unit`), or `OPEN/CLOSE <fname>`; after the
keyword it reads `ICNSTNT FMTIN IPRN`, then `JJ` integers from `LOCAT`, and (if
`ICNSTNT != 0`) multiplies the array by `ICNSTNT`. So `EXTERNAL`/`OPEN-CLOSE`
apply to the recipient list exactly as to any `U1DINT` array. (For DRT this is
the spreading block, `NR<0`; the single inline recipient `NR>0` is a field on
the data line, not a `U1DINT` block.)

**Implementation:** `read_u1dint_list` (`_usgt_returnflow.py`) now resolves
`EXTERNAL <unit>` via `ext_unit_dict` (+ `model.model_ws` for relative paths)
and `OPEN/CLOSE <fname>` (quote-aware, including filenames with spaces) by
reusing `_usgt_list._resolve_external_filename` and a shared
`_usgt_list.parse_open_close`; it applies the `ICNSTNT` multiplier consistently
with `INTERNAL`. DRT and QRT both call this one helper (no duplicated logic),
threading `model`/`ext_unit_dict` through `_parse_drain_tokens` /
`_read_sink_rows`. Internal node ids stay 0-based, the file 1-based.

**Out of scope (explicit failure):** an `EXTERNAL` unit absent from
`ext_unit_dict` raises an actionable `NotImplementedError` (inline the nodes
instead) — never a raw `ValueError` or a silent mis-read.

**Write:** unchanged — recipients are always expanded inline as
`INTERNAL (FREE)` 1-based (`Expanded valid write`); `EXTERNAL`/`OPEN-CLOSE` and
fixed formats are never produced on output.

### Polish — control-tail FMTIN hardening (review follow-up)

The first cut read every list free-format and ignored `FMTIN`. `U1DINT`
(utl7u1.f, `LOCAT>0`) reads `JJ` integers *with* `FMTIN`, list-directed only
when `FMTIN == (FREE)`. `read_u1dint_list` now parses the
`ICNSTNT FMTIN IPRN` tail explicitly (via `_parse_u1dint_tail`) for `INTERNAL`,
`EXTERNAL`, and `OPEN/CLOSE`, and:

- **supports `(FREE)`** (list-directed), the only format FloPy writes;
- **keeps the abbreviated keyword-only form** (`EXTERNAL <unit>` /
  `OPEN/CLOSE <fname>` with no tail) as a read convenience — there is no
  `FMTIN`, so it is treated as `(FREE)`. (Decision: keep it; the package docs
  say recipient U1DINT is *free-format only*.)
- **rejects a fixed Fortran `FMTIN`** (e.g. `(10I8)`) with an actionable
  `NotImplementedError` rather than free-parsing it silently;
- applies the `ICNSTNT` multiplier (non-zero) consistently across all three.

`FMTIN` is validated **before** any `EXTERNAL`/`OPEN-CLOSE` file is opened, and
the read of an opened file is wrapped in `try/finally`, so a bad format or an
unexpected EOF never leaks a file descriptor.

**Tests (`autotest/test_usg_transport.py`):** `test_mfusgqrt_recipient_u1dint_external`,
`test_mfusgdrt_recipient_u1dint_external`,
`test_mfusgqrt_recipient_u1dint_open_close_quoted` (double-quoted name with a
space), `test_mfusgdrt_recipient_u1dint_open_close`, and
`test_mfusgqrt_recipient_u1dint_external_unresolved_fails`. The first two also
assert write_file re-emits the recipients inline 1-based with no
`EXTERNAL`/`OPEN-CLOSE`. Polish adds four more:
`test_mfusgqrt_recipient_u1dint_external_full_control` (`EXTERNAL <u> 1 (FREE) -1`),
`test_mfusgdrt_recipient_u1dint_open_close_full_control` (`... 1 (FREE) -1`),
`test_mfusgqrt_recipient_u1dint_multiplier` (`ICNSTNT=10` scales file `1,2` ->
1-based `10,20` -> 0-based `9,19`), and
`test_mfusgqrt_recipient_u1dint_fixed_format_fails` (`(10I8)` ->
`NotImplementedError`). Results: `-k "mfusgdrt or mfusgqrt or usgt_recipient"`
**64 passed**; focused **238**; exe **4**; combined **242** under the USG-T 2.7
ARM binary.

## Fortran Sources

- `gwf2QRT8u.f`
- `gwf2drt8u.f`
- `utl7u1.f` for `U1DINT`

## Required Work

1. Audit `TRANSIENTQ` in QRT:
   - records,
   - time-series semantics,
   - stress-period interaction,
   - write format.
2. Audit recipient-node `U1DINT` external controls:
   - `INTERNAL`,
   - `CONSTANT`,
   - `EXTERNAL`,
   - `OPEN/CLOSE`.
3. Implement support where reasonable.
4. If preserving external recipient files is too much, at least load and expand
   them with an honest `Expanded valid write` label.

## Required Tests

- QRT `TRANSIENTQ` explicit supported case or explicit documented deferral.
- DRT recipient list via `OPEN/CLOSE`.
- QRT recipient list via `EXTERNAL` resolved through `ext_unit_dict`.
- Expanded write/reload if controls are not preserved.
- No regression in existing inline/constant recipient tests.

## Documentation Updates

Update:

- `USGT_roadmap.md`,
- `USGT_PACKAGE_BACKLOG.md`,
- `USGT_improvements.md`.

## Agent Prompt

```text
Goal: implement the remaining QRT/DRT rare list controls for USG-T 2.7, starting with recipient U1DINT EXTERNAL/OPEN/CLOSE.

Read first:
- USGT_STAGE4_MASTER_PLAN.md
- USGT_STAGE4_05_QRT_DRT_RARE_CONTROLS.md
- flopy/mfusg/mfusgqrt.py
- flopy/mfusg/mfusgdrt.py
- flopy/mfusg/_usgt_returnflow.py
- /Users/roh.sgm/Documents/SGM Co/GDM/01_Software/GSI/USG-Transport 2.7.0 (2026.03.21)/USGT_V_2-7-0_Source_Code/gwf2QRT8u.f
- /Users/roh.sgm/Documents/SGM Co/GDM/01_Software/GSI/USG-Transport 2.7.0 (2026.03.21)/USGT_V_2-7-0_Source_Code/gwf2drt8u.f

Task:
Audit and implement recipient-node U1DINT EXTERNAL/OPEN/CLOSE support for QRT/DRT. Then audit QRT TRANSIENTQ and either implement it or leave a precise explicit unsupported decision with tests. Prefer load-and-expanded-write if full preservation is too large.

Validation:
python -m pytest autotest/test_usg_transport.py -q
python -m pytest autotest/test_usg_transport_exe.py -q
USGT_EXE=/Users/roh.sgm/Documents/GitHub_Projects/GW-Software-Compiled-via-Claude/gfortran/usgt_2.7/usgt_270_arm python -m pytest autotest/test_usg_transport.py autotest/test_usg_transport_exe.py -q
git diff --check
git status --short

Constraints:
- Push only to origin, never upstream.
- Do not mix parameter preservation into this task.
```
