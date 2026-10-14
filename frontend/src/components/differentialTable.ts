// SPDX-License-Identifier: AGPL-3.0-only
// SPDX-FileCopyrightText: 2026 Jonathan D.A. Jewell (hyperpolymath) <j.d.a.jewell@open.ac.uk>
import type { DifferentialResult, DifferentialRow } from '../api/types'

const COLUMNS = [
  'taxon', 'status', 'log2_fold_change', 'estimate', 'standard_error', 'statistic',
  'pvalue', 'padj', 'dispersion_theta', 'prevalence', 'note',
] as const

/** Quote one CSV field when it holds a comma, quote or line break (RFC 4180). */
export function csvField(value: string | number | null): string {
  if (value === null) return ''
  const s = String(value)
  return /[",\r\n]/.test(s) ? `"${s.replace(/"/g, '""')}"` : s
}

/**
 * The full results table as CSV, every taxon including failed and filtered
 * ones. Numbers keep full precision; a missing statistic is an empty field,
 * never 0 or 1. Leading `#` lines record the comparison and method (read.csv
 * skips them with `comment.char = "#"`).
 */
export function differentialCsv(result: DifferentialResult): string {
  const { reference, contrast } = result.groups
  const header = [
    `# Differential abundance: ${contrast} vs ${reference} (reference), rank ${result.rank}, table ${result.table}`,
    `# ${result.method}`,
    `# ${configLine(result)}`,
  ].map(l => l.replace(/[\r\n]+/g, ' '))
  const rows = result.rows.map(r => COLUMNS.map(c => csvField(r[c])).join(','))
  return [...header, COLUMNS.join(','), ...rows].join('\r\n') + '\r\n'
}

/** Format a statistic for display: three significant figures, or a dash when absent. */
export function formatStat(value: number | null): string {
  if (value === null) return '–'
  if (value === 0) return '0'
  const a = Math.abs(value)
  return a < 1e-3 || a >= 1e4 ? value.toExponential(2) : value.toPrecision(3)
}

/**
 * The settings line of the CSV: the method and the parameters it used. The
 * offset belongs to nb_glm and the replacement delta to clr_welch, so each line
 * names only its own method's parameters.
 */
export function configLine(result: DifferentialResult): string {
  const c = result.config
  return c.method === 'clr_welch'
    ? `method clr_welch, replacement_delta ${c.replacement_delta}, min_prevalence ${c.min_prevalence}; estimate is a CLR difference, not a fold change`
    : `method nb_glm, offset ${c.offset}, min_prevalence ${c.min_prevalence}`
}

/** The effect size of one row, from the field the result names in `effect`. */
export function effectValue(result: DifferentialResult, row: DifferentialRow): number | null {
  return row[result.effect.key]
}
