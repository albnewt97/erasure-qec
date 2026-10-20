# CLAUDE.md — standing rules for erasure-qec

## What this repo is
Rotated surface-code Z-memory under two-qubit-gate erasure conversion
(stim 1.16, PyMatching 2.4, sinter 1.16). A herald-conditioned MWPM decoder
sets heralded edges to weight 0; a blind decoder ignores heralds. Two noise
conventions: the repo's constant non-identity budget (data/*.csv) and wu2022
(data/wu/*.csv, validation against Wu et al., Nat. Commun. 13, 4657).

## Evidence rules (non-negotiable)
- Every number in README.md, docs/METHODS.md, docs/VALIDATION.md must come
  from a committed script or results/*.json. README result tables are
  rendered by scripts/render_readme_tables.py; never hand-edit them.
- Never drop a caveat that qualifies a number that is still published.
  docs/archive/remediation/CAVEAT_LEDGER.md is the checklist.
- Report observations before diagnoses. When you propose a mechanism, say
  what result would falsify it. Do not claim a cause you did not test.
- Do not describe results with evaluative adjectives (robust, honest,
  rigorous) in docs.
- Pre-registered definitions and gates live in the remediation plan the
  maintainer gives you; do not change an estimator, threshold, cap or gate to
  make a check pass. If a check fails, stop and report.

## Frozen things
- The circuits of the existing model are pinned by
  tests/test_circuit_golden.py. If that test fails, you changed data validity:
  revert, do not update digests.
- Never edit data/*.csv or data/wu/*.csv by hand. Only the collection scripts
  append to them, and only the maintainer runs collections.
- Do not start any command expected to run longer than 10 minutes. Print it
  (with nohup and a log path) for the maintainer instead.

## Workflow
- One phase per branch; small conventional commits; never push or merge.
- Before every commit:
  uv run ruff check . && uv run ruff format --check . && uv run mypy &&
  uv run pytest -m "not slow" -q
- Run the gate exactly as written and check its exit status directly. Never pipe it
  (`| tail`, `| tee`): a pipeline reports the last command's exit status, not the gate's.
- Do not delete tests. If behaviour changes intentionally, change the call or
  add a new test, and say so in the report.
- End every session with the report format the prompt specifies, then stop.

## Style for code and docs
- Docstrings state what a function computes and under what assumptions;
  history belongs in CHANGELOG.md or docs/archive/.
- One claim per sentence; numbers carry an interval and a method.
