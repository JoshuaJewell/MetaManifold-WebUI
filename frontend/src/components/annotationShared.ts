// © 2026 Joshua Benjamin Jewell. All rights reserved.
// Licensed under the GNU Affero General Public License version 3 (AGPLv3).
import { useEffect, useState } from 'react'
import type { ComparisonRunSpec, ConfigSource, TableMeta } from '../api/types'
import { uniqueRuns } from '../api/types'
import { api } from '../api/client'
import { errorMessage } from '../api/errorMessage'

export const SOURCE_COLORS: Record<ConfigSource, string> = {
  default: 'var(--color-muted-fg)',
  study:   '#e67700',
  group:   '#5c940d',
  run:     'var(--color-primary)',
  placement: 'var(--color-primary)',
  tree:    'var(--color-primary)',
}

export interface AnalysisOption {
  key: string
  table: string
  label: string
}

/**
 * Discover the analysable result tables for a single run (merged, merged_otu,
 * and any saved sub-tables) from the run's results database. Analysis reads the
 * merged results table directly. Rejects when the tables cannot be listed.
 */
export async function discoverResultsTables(
  study: string, run: string, group?: string | null,
): Promise<AnalysisOption[]> {
  const tables: TableMeta[] = await api.results.runTables(study, run, group)
  return tables
    .map(meta => ({ key: meta.id, table: meta.id, label: meta.label }))
    .sort((a, b) => a.label.localeCompare(b.label))
}

/**
 * The options present in every list of `results` (keyed by `key`), sorted by
 * label. No lists means no options.
 */
export function sharedOptions(results: AnalysisOption[][]): AnalysisOption[] {
  const shared = new Map<string, AnalysisOption>()
  results.forEach((opts, index) => {
    const keys = new Set(opts.map(o => o.key))
    if (index === 0) {
      for (const o of opts) shared.set(o.key, o)
    } else {
      for (const key of [...shared.keys()]) {
        if (!keys.has(key)) shared.delete(key)
      }
    }
  })
  return [...shared.values()].sort((a, b) => a.label.localeCompare(b.label))
}

/**
 * Result tables present in every run of `runs`, sorted by label, with the
 * message of a failed lookup in `error` (and no options) so that a failure is
 * not read as "no table shared".
 */
export function useSharedResultsTablesState(
  study: string, runs: ComparisonRunSpec[],
): { options: AnalysisOption[]; error: string | null } {
  const [state, setState] = useState<{ options: AnalysisOption[]; error: string | null }>(
    { options: [], error: null })
  useEffect(() => {
    let cancelled = false
    Promise.all(uniqueRuns(runs).map(r => discoverResultsTables(study, r.run, r.group)))
      .then(results => { if (!cancelled) setState({ options: sharedOptions(results), error: null }) })
      .catch(err => { if (!cancelled) setState({ options: [], error: errorMessage(err) }) })
    return () => { cancelled = true }
  }, [study, runs])
  return state
}

/** Result tables present in every run of `runs`, sorted by label; none when the lookup fails. */
export function useSharedResultsTables(study: string, runs: ComparisonRunSpec[]): AnalysisOption[] {
  return useSharedResultsTablesState(study, runs).options
}
