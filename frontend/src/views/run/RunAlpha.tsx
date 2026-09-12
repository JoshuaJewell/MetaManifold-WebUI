// © 2026 Joshua Benjamin Jewell. All rights reserved.
// Licensed under the GNU Affero General Public License version 3 (AGPLv3).
import { useAnalysis } from '../../hooks/useAnalysis'
import { AnalysisControls } from '../../components/AnalysisControls'
import type { AnnotationSource, ColFilter } from '../../api/types'

export function RunAlpha({ study, run, group, table, filters, source }: {
  study: string; run: string; group?: string; table: string
  filters: Record<string, ColFilter>
  /** Classifier whose categories drive the diversity exclusions. */
  source?: AnnotationSource
}) {
  const analysis = useAnalysis({ study, run, group, table, colFilters: filters, enabled: false })
  const body = { table, colFilters: filters, ...(source ? { source } : {}) }
  const nFilters = Object.keys(filters).length

  return (
    <AnalysisControls {...analysis} body={body}>
      <code>{table}</code>
      {nFilters > 0 && <> with {nFilters} column {nFilters === 1 ? 'filter' : 'filters'} from the Tables tab</>}
    </AnalysisControls>
  )
}
