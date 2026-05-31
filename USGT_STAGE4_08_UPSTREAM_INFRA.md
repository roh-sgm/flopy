# Stage 4.8 - Upstream Infrastructure

Goal: prepare the non-FloPy infrastructure needed for credible upstream USG-T
2.7 support: reproducible source, tagged releases, executable CI, and optional
FloPy test integration.

## Current State

FloPy package support is improving, but upstream acceptance also depends on a
stable executable source/release story.

Current executable validation uses:

```bash
USGT_EXE=/Users/roh.sgm/Documents/GitHub_Projects/GW-Software-Compiled-via-Claude/gfortran/usgt_2.7/usgt_270_arm
```

This is valid locally but not reproducible CI infrastructure.

## Work Streams

1. USG-T source repository:
   - source tree under version control,
   - tagged releases,
   - build instructions,
   - release assets.
2. Executable CI:
   - Linux,
   - macOS,
   - Windows if feasible,
   - smoke tests.
3. FloPy optional integration:
   - keep tests optional,
   - use `USGT_EXE`,
   - later integrate with `get-modflow` or a similar installer only after the
     executable project is reproducible.

## Agent Prompt

```text
Goal: design the external USG-T 2.7 release/CI path needed before upstream FloPy adoption.

Read first:
- USGT_STAGE4_MASTER_PLAN.md
- USGT_STAGE4_08_UPSTREAM_INFRA.md
- USGT_UPSTREAM_INFRA.md
- the November 2024 email notes captured in the project discussion

Task:
Produce an actionable infrastructure plan for a USG-T 2.7 source repository with tagged releases, executable CI, smoke tests, and an optional FloPy integration path. Do not change FloPy package code in this task unless only docs are being clarified.

Validation:
git diff --check
git status --short

Constraints:
- Push only to origin, never upstream.
- Keep executable dependency optional for FloPy.
- Do not claim get-modflow support until there is a reproducible tagged release source.
```
