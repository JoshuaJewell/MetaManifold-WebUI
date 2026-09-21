// SPDX-License-Identifier: MPL-2.0
// SPDX-FileCopyrightText: 2026 Jonathan D.A. Jewell <j.d.a.jewell@open.ac.uk>
// Coupling / drift test — every toolchain pin copy must agree with its
// source of truth.
//
// Category (standards/testing-and-benchmarking/TESTING-TAXONOMY.adoc):
//   Coupling / drift. Pin web under test (see docs/reproducibility.md):
//     mise.toml          — dev-toolchain source of truth (estate convention)
//     .bun-version       — GENERATED from mise.toml (`just sync-pins`);
//                          consumed by CI via bun-version-file
//     config/defaults/tool_versions.yml — upstream pipeline-pin SOT holding
//                          overlapping julia/bun copies
//     .github/workflows/ci.yml — hardcoded julia version on the setup-julia
//                          step. It was a 1x1 `strategy.matrix` until 2026-09-21;
//                          the matrix was removed because GitHub appends the
//                          matrix combination to the posted check name, which
//                          renamed the required status check on every bump.
//   A bump in any one copy without the others is silent CI/local divergence;
//   this test makes it a build-time failure instead.
import { describe, test, expect } from 'bun:test'
import { readFileSync } from 'node:fs'
import { join, resolve } from 'node:path'

const ROOT = resolve(import.meta.dir, '../../..')
const read = (p: string) => readFileSync(join(ROOT, p), 'utf8')

function misePin(tool: string): string {
  const m = read('mise.toml').match(new RegExp(`^${tool}\\s*=\\s*"([^"]+)"`, 'm'))
  if (!m) throw new Error(`mise.toml has no ${tool} pin`)
  return m[1]
}

function toolVersionsPin(tool: string): string {
  const m = read('config/defaults/tool_versions.yml').match(
    new RegExp(`^  ${tool}:\\n    version:\\s*"([^"]+)"`, 'm'))
  if (!m) throw new Error(`tool_versions.yml has no ${tool} pin`)
  return m[1]
}

describe('coupling/drift: toolchain pins agree across all copies', () => {
  test('bun: mise.toml == .bun-version == tool_versions.yml', () => {
    const pin = misePin('bun')
    expect(read('.bun-version').trim()).toBe(pin)
    expect(toolVersionsPin('bun')).toBe(pin)
  })

  test('julia: mise.toml == tool_versions.yml == CI setup-julia step', () => {
    const pin = misePin('julia')
    expect(toolVersionsPin('julia')).toBe(pin)
    // Anchored on the action rather than on a bare `version:` key, so this
    // cannot silently latch onto some other step's version and pass for the
    // wrong reason. A missing anchor yields undefined, which fails loudly.
    const ci = read('.github/workflows/ci.yml')
      .match(/julia-actions\/setup-julia@[0-9a-f]{40}[\s\S]*?version:\s*"([^"]+)"/)
    expect(ci?.[1]).toBe(pin)
  })

  test('node is single-sourced (mise.toml only; no second copy to drift)', () => {
    // Guard for the day someone adds .nvmrc/.node-version: the drift test
    // must be extended, not the pin duplicated silently.
    expect(misePin('node')).toBe('20.20.2')
  })

  test('guard: the pin files this test reads actually exist and are parseable', () => {
    // Extraction-collapse alarms (a vacuous pass on empty matches is the
    // classic drift-test failure mode).
    expect(misePin('just')).toBe('1.43.1')
    expect(toolVersionsPin('r')).toBe('4.5.0')
  })
})
