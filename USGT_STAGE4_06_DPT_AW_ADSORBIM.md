# Stage 4.6 - DPT A-W_ADSORBIM

Goal: implement DPT immobile-domain air-water adsorption (`A-W_ADSORBIM`) or
formally keep it unsupported with a complete spec and stronger tests.

## Current State

`MfUsgDpt.load` fails explicitly when it sees `A-W_ADSORBIM`, before later reads
can shift. This is safe but incomplete.

## Fortran Sources

- `gwt2dptu1.f`
- `dpt2aw_adsorb.f`

## Known Complexity

The option can require:

- `IAREA_FNIM`,
- `IKAWI_FNIM`,
- `NAZONESIM`,
- `NATABROWSIM`,
- `IAWIZONMAPIM`,
- `AWI_AREA_TABIM`,
- `AWAMAXIM`,
- `AWAREA_X2IM/X1IM/X0IM`,
- `ROG_SIGMAIM`,
- `ALANGAWIM`,
- `BLANGAWIM`,
- per-stress-period reads.

## Recommended Approach

1. Write a complete Fortran-derived spec first.
2. Implement the smallest complete semantic subset only if the Fortran
   condition matrix is clear.
3. Add tests for the bare keyword, indexed keyword, tabular mode, and at least
   one non-tabular mode.
4. Keep explicit failure for modes not implemented.

## Required Tests If Implemented

- load/write/reload with `IAREA_FNIM==5` or `IKAWI_FNIM==4`,
- zone map,
- tabular area arrays,
- Langmuir arrays per species,
- non-tabular constant/array mode,
- interaction with `IDPF` and `IDISPIM`,
- explicit failure for unsupported combinations.

## Agent Prompt

```text
Goal: produce a complete Fortran-derived DPT A-W_ADSORBIM implementation plan and implement the first safe semantic subset if feasible.

Read first:
- USGT_STAGE4_MASTER_PLAN.md
- USGT_STAGE4_06_DPT_AW_ADSORBIM.md
- flopy/mfusg/mfusgdpt.py
- /Users/roh.sgm/Documents/SGM Co/GDM/01_Software/GSI/USG-Transport 2.7.0 (2026.03.21)/USGT_V_2-7-0_Source_Code/gwt2dptu1.f
- /Users/roh.sgm/Documents/SGM Co/GDM/01_Software/GSI/USG-Transport 2.7.0 (2026.03.21)/USGT_V_2-7-0_Source_Code/dpt2aw_adsorb.f

Task:
Audit A-W_ADSORBIM completely. If implementation is tractable in one card, implement semantic load/write/authoring for the first complete subset and add tests. If not, update docs with the complete spec and keep explicit NotImplementedError tests for every unsupported branch.

Validation:
python -m pytest autotest/test_usg_transport.py -q
python -m pytest autotest/test_usg_transport_exe.py -q
USGT_EXE=/Users/roh.sgm/Documents/GitHub_Projects/GW-Software-Compiled-via-Claude/gfortran/usgt_2.7/usgt_270_arm python -m pytest autotest/test_usg_transport.py autotest/test_usg_transport_exe.py -q
git diff --check
git status --short

Constraints:
- Push only to origin, never upstream.
- Do not silently parse past A-W_ADSORBIM if any required arrays are unsupported.
```
