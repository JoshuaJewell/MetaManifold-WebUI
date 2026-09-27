# SPDX-License-Identifier: CC-BY-SA-4.0
# SPDX-FileCopyrightText: 2026 Jonathan D.A. Jewell (hyperpolymath) <j.d.a.jewell@open.ac.uk>

# The KYAML pilot

**Status: in progress, by owner ruling of 2026-09-26.** This repository is the pilot for
migrating the estate's YAML to KYAML. The ruling is recorded here, the switch is `just
use-kyaml` / `just use-yaml`, and the escape hatch is one `git revert`.

Authority: `hyperpolymath/standards`, `3-practice/YAML-POLICY.adoc` (rules Y-1, Y-2, Y-3,
adoption order §5). This document is the operating manual for the pilot; the policy stays
the authority, and where the two disagree the policy wins and this file is wrong. The
upstream application touchpoints and the low-risk adoption sequence are mapped in
[`docs/integration/kyaml-upstream-adoption.md`](../integration/kyaml-upstream-adoption.md).

## 1. Why this is worth doing, in one paragraph

YAML's implicit typing and whitespace sensitivity are not a style question here, they are a
defect class the estate has measured: a bare `no` that a reader takes for a string and a
loader takes for `false` (the Norway problem), a version pin `1.0` that loses its trailing
zero, indentation that makes text patching unsafe. KYAML is a *strict subset* of YAML — every
YAML reader already accepts it — that removes the ambiguity at the source: flow style
throughout, every string value double-quoted, keys unquoted only where they cannot be
misread, trailing commas, explicit document header, two-space nesting. KEP-5295's decisive
property is that adopting it needs no new parser anywhere. That is what makes it a
*reversible* change rather than a migration to a new format, and that is the condition under
which the owner ruled it in.

## 2. What the ruling says

* **All YAML in this repository is to be migrated to KYAML.** The pilot covers every
  estate-authored `.yml`/`.yaml` file, workflows included.
* **YAML is deprecated here but not removed.** `just use-yaml` returns the tree to block
  style, and the pilot lands as a commit that can be reverted byte-exactly. Nothing is
  deleted and no reader has to change: KYAML is YAML.
* **Scripts inside workflows are extracted**, not quoted. A KYAML flow scalar cannot be a
  block scalar, so a 30-line `run:` block would become one enormous quoted value. Instead
  each such step becomes `run: "bash scripts/ci/<name>.sh"` and the shell lives in a file
  that shellcheck and review can see. That is better engineering independently of KYAML.
* **Bot drift is accepted in writing.** Dependabot and `gh actions-lock` rewrite `uses:`
  pins in block style. `config/kyaml/drift.txt` lists those two files as *converted but not
  gated*, with the reason, and the reconciliation is one command. YAML-POLICY §5 step 6
  requires exactly this written acceptance before workflow conversion starts; this is it.

## 3. The switch

| Command | What it does |
| --- | --- |
| `just use-kyaml` | Rewrites every estate-owned YAML file as canonical KYAML. |
| `just use-yaml` | Rewrites it back as block-style YAML, comments preserved. |
| `just check-kyaml` | Gate: every non-exempt file is byte-for-byte what the emitter writes. |
| `just kyaml-report` | Nothing is written; per-file decisions are printed (see §5). |

Exit codes: `0` clean, `1` a file is not canonical, `2` a file was **refused** — and a
refusal means nothing was written at all, so a refused tree is never half-converted.

`just check-kyaml` runs inside `just hygiene` and `just ci`: a gate that cannot run is not a
gate, so it is wired into the lanes rather than documented as a suggestion.

## 4. What the switch guarantees, and what it refuses

Guaranteed, and tested in `test/unit/test_kyaml.jl`:

* **Comments survive with their association.** An own-line comment stays on its own line above
  the same entry; an end-of-line comment stays on the same line as the same entry. A comment
  that cannot be placed losslessly is a refusal, never a silent drop.
* **Idempotence.** `check` compares the file against the emitter's own output, so a second
  run cannot give a second answer.
* **Recoverable bytes.** YAML → KYAML → YAML reproduces the canonical block form of the same
  document, and `git revert` of the pilot commit reproduces the pre-pilot file exactly.
* **A mutant dies.** `test/unit/test_kyaml.jl` deletes one comment from a converted file and
  asserts the check goes red: a check that has never failed is not a check (policy §2.2).

Refused, by name and line number, leaving the file untouched:

| Construct | Why it is refused |
| --- | --- |
| Anchors, aliases, tags, merge keys | They carry identity between nodes; a rewrite can silently change which nodes share a value. |
| Multiple documents, directives, `...` | Out of scope for a single-document pilot; the corpus has none. |
| Tabs in indentation | Illegal in YAML; guessing the intent is exactly the ambiguity KYAML removes. |
| Duplicate keys in one mapping | Last-wins is a loader detail; re-emitting would pick a winner silently. |
| Multi-line plain scalars | Folded into one line by YAML rules, so the source bytes are not recoverable. |
| An end-of-line comment on a block scalar | The scalar owns the rest of the line; there is nowhere lossless to put the comment. |

This repository's corpus needs none of those: a census of the 17 tracked YAML files found
**zero** anchors, aliases, tags, multi-document streams or directives; block scalars in two
files (`ci.yml`, 22 of them; `.github/ISSUE_TEMPLATE/bug_report.yml`, 2); 12 `~` nulls; and
comment-bearing lines concentrated in `ci.yml` (313), `config/defaults/pipeline.yml` (97) and
`config/defaults/tool_versions.yml` (64). That census is why the refusal list can be short
and honest instead of pretending to be a general YAML implementation — a general one is
standards#1022's job, not this tool's.

## 5. Decisions the switch makes for you, and prints

KYAML removes implicitness, which means the switch must sometimes choose. Every choice is
counted and printed by `just kyaml-report`, so a reviewer sees them instead of trusting them:

* `~` and an empty value become `null` (one spelling of the empty value, not three).
* A plain scalar that is not a canonical integer, float, `true` or `false` is double-quoted —
  so a bare `no` stops being a boolean by accident. This is the Norway fix, and it is a
  *change in meaning* for a YAML 1.1 reader; the report counts it rather than hiding it.
* A key that is schema-ambiguous (`on`, `off`, `yes`, `no`, `y`, `n`, `true`, `false`,
  `null`, `~`) is quoted. In GitHub workflows the canonical key is `on`, and `"on":` is the
  same key after parsing — the quote is for the readers that are not GitHub's.
* An end-of-line comment on a key whose value is a collection moves to its own line above the
  key, because a flow collection ends several lines later.
* A blank line that sat between a comment and the entry it precedes moves above the comment:
  the emitter writes blanks, then comments, then the entry. Comments and their entries stay
  together; one blank line's position can change, once, and never again.

## 6. Proof obligations and where each stands

YAML-POLICY §5 fixes an order; this is the pilot's position in it.

| Step | Issue | State here |
| --- | --- | --- |
| 1. GitHub parses a KYAML workflow | standards#1020 (closed) | discharged estate-wide by a two-arm probe. Re-probed **in this repository** by the pilot branch before merge; the arm and its result are recorded in the PR. |
| 2. Comment-preservation proof | standards#1021 | implemented here as property tests plus a dropped-comment mutant (`test/unit/test_kyaml.jl`). Contributes evidence to #1021; does not close it. |
| 3. A formatter/linter for arbitrary YAML | standards#1022 | **NOT** done estate-wide. This tool is repository-scoped and refuses what it cannot emit losslessly. That is the boundary, stated rather than blurred. |
| 4. Owner ruling on scope | standards#1023 | ruled 2026-09-26: this repository is the pilot, all YAML in scope. |
| 5. Migrate non-bot YAML | standards#1024 | this branch. |
| 6. Workflows | standards#1025 | ruled IN. Written acceptance of bot drift: §2 and `config/kyaml/drift.txt`. |

## 7. Reverting

Three levels, cheapest first:

1. `just use-yaml` — back to block style from the same document model. Comments intact.
2. `git revert <pilot commit>` — the pre-pilot bytes, because that is what a revert is.
3. Delete `scripts/kyaml/`, `config/kyaml/` and this file, drop `check-kyaml` from the
   Justfile lanes: nothing else in the tree depends on the pilot.

## 8. Follow-ups this pilot does not do

* Extract the shell out of `.github/workflows/ci.yml` into `scripts/ci/*.sh` (§2) — the
  conversion commit's largest diff, and the reason the workflow files stay reviewable.
* A CI step running `just check-kyaml` — it lands in the same commit as the conversion, so
  that the gate is never red for a reason that has nothing to do with the change under
  review. The Agda proofs lane is already wired (`proofs` job in ci.yml, gated on
  `vars.STAPELN_AGDA_IMAGE`, with the reason printed by the hygiene job when it is unset).
* `gh actions-lock` emitting KYAML, or an explicit reconciliation note per bot PR.
* The estate-wide formatter (standards#1022) should absorb this tool's parser and its
  refusal list rather than grow a second implementation.
