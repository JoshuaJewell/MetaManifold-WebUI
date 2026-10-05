// SPDX-License-Identifier: AGPL-3.0-only
// SPDX-FileCopyrightText: 2026 Jonathan D.A. Jewell (hyperpolymath) <j.d.a.jewell@open.ac.uk>
import type { ComparisonRunSpec } from '../api/types'
import { uniqueRuns } from '../api/types'

/**
 * The taxonomy ranks offered by every run in `runs`, in the first run's order.
 * Rejects when any run's ranks cannot be read: a failed request must not read
 * as "no rank shared by the selected runs".
 */
export async function sharedRanks(
  runs: ComparisonRunSpec[],
  ranksOf: (run: ComparisonRunSpec) => Promise<string[]>,
): Promise<string[]> {
  const results = await Promise.all(uniqueRuns(runs).map(ranksOf))
  return results.reduce<string[]>((acc, cur) => acc.filter(r => cur.includes(r)), results[0] ?? [])
}
