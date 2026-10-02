// SPDX-License-Identifier: AGPL-3.0-only
// SPDX-FileCopyrightText: 2026 Jonathan D.A. Jewell (hyperpolymath) <j.d.a.jewell@open.ac.uk>
import { describe, expect, test } from 'bun:test'
import type { DifferentialResult } from '../api/types'
import { configLine, csvField, differentialCsv, effectValue, formatStat } from './differentialTable'

const result: DifferentialResult = {
  status: 'partial',
  method: 'NB GLM',
  rank: 'Genus',
  table: 'merged',
  groups: { reference: 'run_A', contrast: 'run_B' },
  n_samples: { run_A: 2, run_B: 2 },
  effect: { key: 'log2_fold_change', label: 'log2 fold change' },
  size_factors: [],
  config: { method: 'nb_glm', offset: 'tss', min_prevalence: 0 },
  diagnostics: { n_taxa: 2, n_tested: 1, n_failed: 1, n_boundary: 0, n_filtered: 0 },
  rows: [
    { taxon: 'Bacteroides', status: 'ok', note: '', estimate: 1.3862943611198906,
      log2_fold_change: 2, standard_error: 0.25, statistic: 5.545, pvalue: 2.9e-8, padj: 5.8e-8,
      dispersion_theta: 12.5, prevalence: 1 },
    { taxon: 'Odd, "quoted"', status: 'failed', note: 'the counts are constant across samples: there is nothing to estimate',
      estimate: null, log2_fold_change: null, standard_error: null, statistic: null, pvalue: null, padj: null,
      dispersion_theta: null, prevalence: 0.5 },
  ],
  figure: {},
}

describe('differential CSV', () => {
  test('quotes only fields that need it', () => {
    expect(csvField('plain')).toBe('plain')
    expect(csvField('a,b')).toBe('"a,b"')
    expect(csvField('say "hi"')).toBe('"say ""hi"""')
    expect(csvField(null)).toBe('')
    expect(csvField(0)).toBe('0')
  })

  test('keeps every taxon, full precision, and leaves missing statistics empty', () => {
    const lines = differentialCsv(result).trimEnd().split('\r\n')
    expect(lines[0]).toContain('run_B vs run_A (reference)')
    expect(lines[3]).toBe('taxon,status,log2_fold_change,estimate,standard_error,statistic,pvalue,padj,dispersion_theta,prevalence,note')
    expect(lines[4]).toBe('Bacteroides,ok,2,1.3862943611198906,0.25,5.545,2.9e-8,5.8e-8,12.5,1,')
    expect(lines[5].startsWith('"Odd, ""quoted""",failed,,,,,,,,0.5,')).toBe(true)
    expect(lines).toHaveLength(6)
  })
})

describe('formatStat', () => {
  test('shows a dash for an absent value, never zero', () => {
    expect(formatStat(null)).toBe('–')
    expect(formatStat(0)).toBe('0')
    expect(formatStat(0.04321)).toBe('0.0432')
    expect(formatStat(2.9e-8)).toBe('2.90e-8')
  })
})

describe('clr_lm results', () => {
  const clr: DifferentialResult = {
    ...result,
    method: 'CLR LM',
    effect: { key: 'estimate', label: 'CLR difference' },
    config: { method: 'clr_lm', min_prevalence: 0, replacement_delta: 0.65 },
    rows: [{ ...result.rows[0], estimate: -0.75, log2_fold_change: null, dispersion_theta: null }],
  }

  test('the effect comes from the field the result names, not log2_fold_change', () => {
    expect(effectValue(clr, clr.rows[0])).toBe(-0.75)
    expect(effectValue(result, result.rows[0])).toBe(2)
  })

  test('the settings line names the delta and says the estimate is not a fold change', () => {
    expect(configLine(clr)).toBe('method clr_lm, replacement_delta 0.65, min_prevalence 0; estimate is a CLR difference, not a fold change')
    expect(configLine(result)).toBe('method nb_glm, offset tss, min_prevalence 0')
    expect(differentialCsv(clr).split('\r\n')[2]).toBe(`# ${configLine(clr)}`)
  })
})
