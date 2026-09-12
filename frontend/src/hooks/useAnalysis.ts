// © 2026 Joshua Benjamin Jewell. All rights reserved.
// Licensed under the GNU Affero General Public License version 3 (AGPLv3).
import { useCallback, useEffect, useMemo, useState } from 'react'
import { api } from '../api/client'
import { errorMessage } from '../api/errorMessage'
import { useToast } from '../components/Toast'
import type { AnalysisRequest, AnnotationSource, ColFilter } from '../api/types'

interface UseAnalysisOpts {
  study: string
  run: string
  group?: string | null
  table: string | null
  source?: AnnotationSource
  colFilters?: Record<string, ColFilter>
  prefix?: string | null
  /** When false, rank fetching is suppressed (e.g. table not yet ready). Default true. */
  enabled?: boolean
}

export interface UseAnalysisResult {
  ranks: string[]
  alphaFig: unknown
  loading: boolean
  runAlpha: () => Promise<void>
  resetFigures: () => void
}

export function useAnalysis(opts: UseAnalysisOpts): UseAnalysisResult {
  const { study, run, group, table, source, colFilters, prefix, enabled = true } = opts
  const toast = useToast()

  const [ranks, setRanks]       = useState<string[]>([])
  const [alphaFig, setAlphaFig] = useState<unknown>(null)
  const [loading, setLoading]   = useState(false)

  useEffect(() => {
    if (!table || !enabled) {
      setRanks([])
      return
    }
    let cancelled = false
    api.analysis.ranks(study, run, { table, group: group ?? undefined, source })
      .then(r => { if (!cancelled) setRanks(r) })
      .catch(() => { if (!cancelled) setRanks([]) })
    return () => { cancelled = true }
  }, [study, run, group, source, table, enabled])

  const body = useMemo((): AnalysisRequest | null => {
    if (!table) return null
    return { table, source, colFilters, prefix: prefix || undefined }
  }, [table, source, colFilters, prefix])

  const runAlpha = useCallback(async () => {
    if (!body) return
    setLoading(true)
    try {
      setAlphaFig(await api.analysis.alpha(study, run, body, group))
    } catch (err) {
      toast.error('Alpha analysis failed: ' + (errorMessage(err)))
    } finally { setLoading(false) }
  }, [body, group, run, study, toast])

  const resetFigures = useCallback(() => {
    setAlphaFig(null)
  }, [])

  return { ranks, alphaFig, loading, runAlpha, resetFigures }
}
