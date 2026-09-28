// © 2026 Joshua Benjamin Jewell. All rights reserved.
// Licensed under the GNU Affero General Public License version 3 (AGPLv3).
import { useEffect, useState } from 'react'
import type { ComparisonRunSpec, ConfigSource, TableMeta } from '../api/types'
import { uniqueRuns } from '../api/types'
import { api } from '../api/client'

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
 * merged results table directly.
 */
export async function discoverResultsTables(
  study: string, run: string, group?: string | null,
): Promise<AnalysisOption[]> {
  const tables = await api.results.runTables(study, run, group).catch(() => [] as TableMeta[])
  return tables
    .map(meta => ({ key: meta.id, table: meta.id, label: meta.label }))
    .sort((a, b) => a.label.localeCompare(b.label))
}

/** Result tables present in every run of `runs`, sorted by label. */
export function useSharedResultsTables(study: string, runs: ComparisonRunSpec[]): AnalysisOption[] {
  const [options, setOptions] = useState<AnalysisOption[]>([])
  useEffect(() => {
    let cancelled = false
    Promise.all(uniqueRuns(runs).map(r => discoverResultsTables(study, r.run, r.group)))
      .then(results => {
        if (cancelled) return
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
        setOptions([...shared.values()].sort((a, b) => a.label.localeCompare(b.label)))
      })
    return () => { cancelled = true }
  }, [study, runs])
  return options
}
