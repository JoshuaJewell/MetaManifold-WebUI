// SPDX-FileCopyrightText: 2026 Jonathan D.A. Jewell (hyperpolymath) <j.d.a.jewell@open.ac.uk>
// SPDX-License-Identifier: MPL-2.0

export type ExactInteger = number | string
export interface DescriptiveSample {
  label: string
  status: string
  total: ExactInteger
  counts?: ExactInteger[]
  features?: string[]
  proportions?: (string | null)[]
  richness?: number
  warnings: string[]
}
export interface DescriptiveResult {
  mode: 'exact_counts'
  status: string
  claim: string
  samples: DescriptiveSample[]
  analysis_id: string
  input_sha256: string
  review_status: string
}

/** Preserve exact strings from the server; never round through a JS Number. */
export function exactIntegerText(value: ExactInteger): string {
  if (typeof value === 'number') {
    if (!Number.isSafeInteger(value) || value < 0) throw new Error('The response contains an unsafe count.')
    return String(value)
  }
  if (!/^(0|[1-9][0-9]*)$/.test(value)) throw new Error('The response contains an invalid exact count.')
  return value
}

export function proportionText(value: string | null): string {
  if (value === null) return 'Undefined (zero total)'
  if (!/^[0-9]+\/[1-9][0-9]*$/.test(value)) throw new Error('The response contains an invalid exact proportion.')
  return value
}
