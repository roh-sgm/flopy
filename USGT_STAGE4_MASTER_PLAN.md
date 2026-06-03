# USG-T 2.7 Stage 4 Master Plan

Date: 2026-05-31

> **Status (post-Stage 4.6G):** this plan has been **executed** — Stages
> 4.6A–4.6G are done. Per-section "currently …" / "Target state …" notes below
> describe the state *when the plan was written* and are now historical. The
> authoritative final state is `USGT_STAGE4_09_FINAL_GAP_AUDIT.md` (and the
> `USGT_roadmap.md` coverage table). Notably DPT `A-W_ADSORBIM` is no longer
> "explicit unsupported" — Stage 4.6E implemented the array-only branches (the
> scalar/tabular branches stay deferred).

Purpose: define the remaining work to move the fork from strong USG-T 2.7
support to a more complete feature surface. Stage 3 closed the core
upstream-readiness concerns around honest status, executable validation,
`USGT_EXE`, and test isolation. Stage 4 now prioritizes the features the project
actually uses first, then broadens toward literal package completeness.

Primary source of truth:

`/Users/roh.sgm/Documents/SGM Co/GDM/01_Software/GSI/USG-Transport 2.7.0 (2026.03.21)/USGT_V_2-7-0_Source_Code`

Secondary references:

- `USGT_roadmap.md`
- `USGT_PACKAGE_BACKLOG.md`
- `USGT_improvements.md`
- `USGT_UPSTREAM_INFRA.md`
- `USGT_STAGE3_REREVIEW.md`

## Strategy

Work in small, reviewable stages. Do not turn every remaining package into a
giant PR. The primary Stage 4 objective is **programmatic authoring**: users must
be able to construct USG-T packages from Python/numpy inputs and write valid
USG-T 2.7 files without first loading an existing model. Round-trip support is
still required, but it is not sufficient by itself.

Each task should finish with:

- a Fortran-derived I/O note,
- semantic authoring APIs where claimed,
- valid write support from scratch,
- load/write round-trip tests,
- from-scratch authoring tests,
- explicit unsupported-mode failures where support is intentionally deferred,
- updated status docs.

## Priority Order

### Stage 4.1 - TIB Semantic Support

Plan: `USGT_STAGE4_01_TIB_SEMANTIC.md`

Reason: the project uses TIB heavily and it is currently only a raw/text
round-tripper. This is the first authoring gap to close.

Target state: promote from `Raw/text round-trip` to `Full semantic` if the parser
and constructor cover the Fortran grammar, including transport concentration
blocks. If any TIB submode is deferred, keep an honest partial label.

### Stage 4.2 - GSF Semantic Support

Plan: `USGT_STAGE4_02_GSF_SEMANTIC.md`

Reason: the project uses GSF heavily. It currently preserves text and can build
an `UnstructuredGrid`, but it lacks a semantic constructor/editor.

Target state: semantic read/write/authoring for GSF grid data, with exact 1-based
file and 0-based internal vertex/node conventions documented.

### Stage 4.3 - OC / EVT / MDT / LAK Fullness Pass

Plan: `USGT_STAGE4_03_OC_EVT_MDT_LAK.md`

Reason: these are implemented enough to work but are intentionally not `Full`.
This stage promotes them only if authoring and Fortran audit are strong enough.

Recommended package order:

1. `OC`
2. `EVT`
3. `MDT`
4. `LAK`

### Stage 4.4 - MODFLOW Parameter Preservation

Plan: `USGT_STAGE4_04_PARAMETERS.md`

Reason: parameters are currently handled honestly but incompletely. ETS expands
parameters to arrays; HFB/SGB/QRT/DRT fail explicitly.

Target state: either full parameter preservation for the scoped packages, or a
shared implementation plan with clear non-support labels. Because the user wants
complete USG-T support, this stage should aim for implementation.

### Stage 4.5 - QRT/DRT Rare Controls

Plan: `USGT_STAGE4_05_QRT_DRT_RARE_CONTROLS.md`

Reason: QRT/DRT are mostly strong, but rare controls remain unsupported:
`TRANSIENTQ` and recipient-node `U1DINT` `EXTERNAL`/`OPEN/CLOSE`.

Target state: implement or clearly close each remaining rare control.

### Stage 4.6 - DPT Immobile Air-Water Adsorption

Plan: `USGT_STAGE4_06_DPT_AW_ADSORBIM.md`

Reason: `A-W_ADSORBIM` is currently an explicit unsupported mode. It is large
and rare, so it follows project-critical features and parameter work.

Target state: semantic support for DPT immobile air-water adsorption, or a
formal decision to remain unsupported.

### Stage 4.7 - Compatibility-Only Packages

Plan: `USGT_STAGE4_07_COMPAT_PACKAGES.md`

Reason: SFR/STR/GAGE/FHB/SUB/SWT are compatibility-only. They should be audited
after project-used packages are complete, because CLN is the preferred coupling
for current work.

Target state: either implement USG-T-specific variants where needed or document
per-package deferral with stronger evidence.

### Stage 4.8 - Upstream Infrastructure

Plan: `USGT_STAGE4_08_UPSTREAM_INFRA.md`

Reason: full package support is not enough for upstream acceptance. USG-T needs a
reproducible executable source/release/CI story.

Target state: a tagged USG-T 2.7 source/release path with executable CI and a
clear optional FloPy integration path.

## Review Rhythm

For each stage:

1. Agent implements only the scoped task.
2. Agent pushes to `origin`, never upstream.
3. Reviewer checks diff, tests, docs, and status labels.
4. Reviewer issues a new prompt only if follow-up is needed.
5. When accepted, move to the next Stage 4 file.

## Standard Validation

Run at least:

```bash
python -m pytest autotest/test_usg_transport.py -q
python -m pytest autotest/test_usg_transport_exe.py -q
USGT_EXE=/Users/roh.sgm/Documents/GitHub_Projects/GW-Software-Compiled-via-Claude/gfortran/usgt_2.7/usgt_270_arm python -m pytest autotest/test_usg_transport.py autotest/test_usg_transport_exe.py -q
git diff --check
git status --short
```

For feature-specific work, add targeted tests and include them in the agent
summary.

## Current Next Task

Start with `USGT_STAGE4_01_TIB_SEMANTIC.md`.
