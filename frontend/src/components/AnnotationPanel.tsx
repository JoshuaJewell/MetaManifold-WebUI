// © 2026 Joshua Benjamin Jewell. All rights reserved.
// Licensed under the GNU Affero General Public License version 3 (AGPLv3).

import { useCallback, useEffect, useMemo, useState } from 'react'
import type { ReactNode } from 'react'
import { api } from '../api/client'
import { DataTable } from './DataTable'
import { useToast } from './Toast'
import { AddFuncdbModal, ContamFilterConfig } from './AnnotationPanelControls'
import { CONTAM_STYLE, prefillFromRow, RANK_COL, SOURCES, type ContamStatus } from './annotationShared'
import type { AnnotationMeta, AnnotationSource, ColFilter, ContaminationStats, TableQuery } from '../api/types'

type UIStatus = AnnotationMeta['status'] | 'generating' | 'error'

interface ContamFilterState {
  blacklist: Record<string, string>
  whitelist: Record<string, string>
}

function emptyContamConfig(): ContamFilterState {
  return { blacklist: {}, whitelist: {} }
}

function parseConfigList(cfg: Record<string, { value: unknown }>, key: string): string {
  const value = cfg[key]?.value
  return Array.isArray(value) ? value.map(String).join('\n') : ''
}

function buildContamConfig(cfg: Record<string, { value: unknown }>): ContamFilterState {
  return {
    blacklist: {
      function: parseConfigList(cfg, 'annotation.contamination.blacklist.function'),
      detailed_function: parseConfigList(cfg, 'annotation.contamination.blacklist.detailed_function'),
      associated_organism: parseConfigList(cfg, 'annotation.contamination.blacklist.associated_organism'),
      associated_material: parseConfigList(cfg, 'annotation.contamination.blacklist.associated_material'),
      environment: parseConfigList(cfg, 'annotation.contamination.blacklist.environment'),
    },
    whitelist: {
      function: parseConfigList(cfg, 'annotation.contamination.whitelist.function'),
      detailed_function: parseConfigList(cfg, 'annotation.contamination.whitelist.detailed_function'),
      associated_organism: parseConfigList(cfg, 'annotation.contamination.whitelist.associated_organism'),
      associated_material: parseConfigList(cfg, 'annotation.contamination.whitelist.associated_material'),
      environment: parseConfigList(cfg, 'annotation.contamination.whitelist.environment'),
    },
  }
}

function parseFilterLines(value: string): string[] {
  return value.split('\n').map(line => line.trim()).filter(Boolean)
}

function buildFilterPayload(config: ContamFilterState) {
  return {
    blacklist: Object.fromEntries(
      Object.entries(config.blacklist).map(([key, value]) => [key, parseFilterLines(value)]),
    ),
    whitelist: Object.fromEntries(
      Object.entries(config.whitelist).map(([key, value]) => [key, parseFilterLines(value)]),
    ),
  }
}

function AnnotationStatusBadge({ status }: { status: UIStatus }) {
  const statusStyles: Record<UIStatus, { bg: string; fg: string; label: string }> = {
    missing:    { bg: '#6b7280', fg: '#fff', label: 'Missing' },
    fresh:      { bg: '#10b981', fg: '#fff', label: 'Fresh' },
    stale:      { bg: '#f59e0b', fg: '#000', label: 'Stale' },
    generating: { bg: '#3b82f6', fg: '#fff', label: 'Generating...' },
    error:      { bg: '#ef4444', fg: '#fff', label: 'Error' },
  }

  const style = statusStyles[status]
  return (
    <span style={{
      display: 'inline-block', padding: '2px 8px', borderRadius: 4,
      fontSize: '.75rem', fontWeight: 600, background: style.bg, color: style.fg,
      marginLeft: 8, verticalAlign: 'middle',
    }}>
      {style.label}
    </span>
  )
}

export function AnnotationPanel({ study, run, group }: { study: string; run: string; group?: string }) {
  const toast = useToast()
  const [source, setSource] = useState<AnnotationSource>('VSEARCH')
  const [listing, setListing] = useState<AnnotationMeta[]>([])
  const [selected, setSelected] = useState<string | null>(null)
  const [uiStatus, setUiStatus] = useState<UIStatus>('missing')
  const [errorMsg, setErrorMsg] = useState<string | null>(null)
  const [loading, setLoading] = useState(false)
  const [contamOverrides, setContamOverrides] = useState<Record<string, ContamStatus>>({})
  const [addFuncdbPrefill, setAddFuncdbPrefill] = useState<Record<string, string> | null>(null)
  const [defaultModifiedBy, setDefaultModifiedBy] = useState('')
  const [contamConfig, setContamConfig] = useState<ContamFilterState>(emptyContamConfig)
  const [contamStats, setContamStats] = useState<ContaminationStats | null>(null)
  const [applying, setApplying] = useState(false)
  const [tableKey, setTableKey] = useState(0)

  const selectedMeta = useMemo(
    () => listing.find(item => item.table === selected) ?? null,
    [listing, selected],
  )

  useEffect(() => {
    api.config.getDefault().then(cfg => {
      const value = cfg.your_name?.value
      if (typeof value === 'string') setDefaultModifiedBy(value)
    }).catch(() => {})
  }, [])

  useEffect(() => {
    if (!study || !run) return
    api.config.getRun(study, run, group ?? null).then(cfg => {
      setContamConfig(buildContamConfig(cfg as Record<string, { value: unknown }>))
    }).catch(() => {})
  }, [study, run, group])

  const fetchListing = useCallback(async () => {
    try {
      setLoading(true)
      const items = await api.annotations.list(study, run, source, group)
      setListing(items)
      const preferred = (selected ? items.find(item => item.table === selected) : null)
        ?? items.find(item => item.status === 'fresh' || item.status === 'stale')
        ?? items[0]

      setSelected(preferred?.table ?? null)
      setUiStatus(preferred?.status ?? 'missing')
      setErrorMsg(null)
    } catch (err) {
      toast.error('Failed to load annotation listing')
      setErrorMsg(err instanceof Error ? err.message : String(err))
      setUiStatus('error')
    } finally {
      setLoading(false)
    }
  }, [study, run, source, group, selected, toast])

  useEffect(() => { void fetchListing() }, [fetchListing])

  const fetchStats = useCallback(async () => {
    if (!selected || (uiStatus !== 'fresh' && uiStatus !== 'stale')) {
      setContamStats(null)
      return
    }
    try {
      setContamStats(await api.annotations.contaminationStats(study, run, source, selected, group))
    } catch {
      setContamStats(null)
    }
  }, [study, run, source, selected, group, uiStatus])

  useEffect(() => { void fetchStats() }, [fetchStats])

  useEffect(() => {
    setContamOverrides({})
    setContamStats(null)
  }, [source, selected])

  const handleSelect = useCallback((table: string) => {
    setSelected(table)
    const meta = listing.find(item => item.table === table)
    if (meta) {
      setUiStatus(meta.status)
      setErrorMsg(null)
    }
  }, [listing])

  const handleGenerate = useCallback(async () => {
    if (!selected) return
    setUiStatus('generating')
    setErrorMsg(null)
    try {
      await api.annotations.generate(study, run, source, selected, group)
      await fetchListing()
    } catch (err) {
      setUiStatus('error')
      setErrorMsg(err instanceof Error ? err.message : String(err))
      toast.error('Annotation generation failed')
    }
  }, [fetchListing, group, run, selected, source, study, toast])

  const handleContaminationChange = useCallback(async (
    row: Record<string, unknown>,
    newStatus: ContamStatus,
  ) => {
    if (!selected) return
    const rank = String(row.funcdb_match_rank ?? '')
    if (!rank || rank === 'unmatched') return

    const taxonCol = RANK_COL[rank]?.[source]
    const taxon = taxonCol ? String(row[taxonCol] ?? '') : ''
    if (!taxon) return

    const overrideKey = `${rank}:${taxon.toLowerCase().trim()}`
    setContamOverrides(current => ({ ...current, [overrideKey]: newStatus }))

    try {
      const response = await api.annotations.updateContamination(
        study, run, source, selected, rank, taxon, newStatus, group,
      )
      if (response.rows_affected > 1) toast.success(`${response.rows_affected} rows updated`)
    } catch (err) {
      setContamOverrides(current => {
        const next = { ...current }
        delete next[overrideKey]
        return next
      })
      toast.error(err instanceof Error ? err.message : 'Update failed')
    }
  }, [group, run, selected, source, study, toast])

  const handleApplyFilter = useCallback(async () => {
    if (!selected) return
    setApplying(true)
    try {
      const response = await api.annotations.applyContaminationFilter(
        study, run, source, selected, buildFilterPayload(contamConfig), group,
      )
      setContamOverrides({})
      setContamStats(response.contamination_stats)
      setUiStatus(response.status)
      setErrorMsg(null)
      setListing(current => current.map(item =>
        item.table === selected
          ? { ...item, status: response.status, rows: response.rows, generated_at: response.generated_at }
          : item,
      ))
      setTableKey(current => current + 1)
      await fetchListing()
      toast.success('Contamination filter applied')
    } catch (err) {
      setErrorMsg(err instanceof Error ? err.message : 'Apply failed')
      toast.error(err instanceof Error ? err.message : 'Apply failed')
    } finally {
      setApplying(false)
    }
  }, [contamConfig, fetchListing, group, run, selected, source, study, toast])

  const fetcher = useCallback(
    (query: TableQuery) => api.annotations.query(study, run, source, selected!, query, group),
    [group, run, selected, source, study],
  )

  const distinctFetcher = useCallback(
    (column: string, activeFilters?: Record<string, ColFilter>) =>
      api.annotations.distinct(study, run, source, selected!, column, activeFilters, group),
    [group, run, selected, source, study],
  )

  const cellRenderer = useCallback((
    column: string,
    value: string,
    row: Record<string, unknown>,
  ): ReactNode | null => {
    if (column !== 'Contamination') return null

    const rank = String(row.funcdb_match_rank ?? '')
    const taxonCol = RANK_COL[rank]?.[source]
    const taxon = taxonCol ? String(row[taxonCol] ?? '').toLowerCase().trim() : ''
    const overrideKey = `${rank}:${taxon}`
    const resolved = (contamOverrides[overrideKey] ?? (value || 'unassigned')) as ContamStatus

    return (
      <select
        value={resolved}
        onChange={event => {
          event.stopPropagation()
          void handleContaminationChange(row, event.target.value as ContamStatus)
        }}
        onClick={event => event.stopPropagation()}
        style={{
          border: 'none',
          background: 'transparent',
          cursor: 'pointer',
          fontSize: 'inherit',
          padding: 0,
          ...CONTAM_STYLE[resolved],
        }}
      >
        <option value="unassigned">unassigned</option>
        <option value="yes">yes</option>
        <option value="no">no</option>
      </select>
    )
  }, [contamOverrides, handleContaminationChange, source])

  const extraRowActions = useCallback((row: Record<string, unknown>) => (
    <button
      className="btn"
      style={{ padding: '1px 6px', fontSize: '.75rem', marginLeft: 4 }}
      title="Add this taxon to FuncDB"
      onClick={event => {
        event.stopPropagation()
        setAddFuncdbPrefill(prefillFromRow(row, source))
      }}
    >
      +
    </button>
  ), [source])

  const showTable = uiStatus === 'fresh' || uiStatus === 'stale'
    || (uiStatus === 'generating' && selectedMeta?.status !== 'missing')

  return (
    <div className="card">
      <div className="tabs" style={{ marginBottom: 12 }}>
        {SOURCES.map(item => (
          <button
            key={item}
            className={`tab ${source === item ? 'active' : ''}`}
            onClick={() => setSource(item)}
          >
            {item}
          </button>
        ))}
      </div>

      {loading && <p style={{ color: 'var(--color-muted-fg)' }}>Loading...</p>}

      {!loading && listing.length === 0 && (
        <p style={{ color: 'var(--color-muted-fg)' }}>No annotation tables available for {source}.</p>
      )}

      {!loading && listing.length > 0 && (
        <>
          <div style={{ display: 'flex', gap: 8, alignItems: 'center', marginBottom: 12, flexWrap: 'wrap' }}>
            <label style={{ fontWeight: 600, fontSize: '.85rem' }}>Table:</label>
            <select
              value={selected ?? ''}
              onChange={event => handleSelect(event.target.value)}
              disabled={uiStatus === 'generating'}
              style={{ padding: '4px 8px', borderRadius: 4, border: '1px solid var(--color-border)' }}
            >
              {listing.map(item => (
                <option key={item.table} value={item.table}>
                  {item.table} {item.rows != null ? `(${item.rows} rows)` : ''}
                </option>
              ))}
            </select>

            {selected && <AnnotationStatusBadge status={uiStatus} />}

            {selected && uiStatus !== 'generating' && (
              <button
                className={`btn ${uiStatus === 'stale' || uiStatus === 'missing' || uiStatus === 'error' ? 'btn-primary' : ''}`}
                onClick={handleGenerate}
              >
                {uiStatus === 'missing' ? 'Generate' : uiStatus === 'error' ? 'Retry' : 'Regenerate'}
              </button>
            )}

            {uiStatus === 'generating' && (
              <span style={{ fontSize: '.85rem', color: 'var(--color-muted-fg)' }}>Generating...</span>
            )}
          </div>

          {uiStatus === 'stale' && (
            <div style={{
              padding: '8px 12px', marginBottom: 12, borderRadius: 4,
              background: '#fef3c7', color: '#92400e', fontSize: '.85rem',
              border: '1px solid #fcd34d',
            }}>
              Source data has changed - regenerate to update
            </div>
          )}

          {uiStatus === 'error' && errorMsg && (
            <div style={{
              padding: '8px 12px', marginBottom: 12, borderRadius: 4,
              background: '#fee2e2', color: '#991b1b', fontSize: '.85rem',
              border: '1px solid #fca5a5',
            }}>
              {errorMsg}
            </div>
          )}

          {uiStatus === 'generating' && (
            <div style={{
              padding: '8px 12px', marginBottom: 12, borderRadius: 4,
              background: '#dbeafe', color: '#1e40af', fontSize: '.85rem',
              border: '1px solid #93c5fd',
            }}>
              Generating...
            </div>
          )}

          {selected && showTable && (
            <ContamFilterConfig
              config={contamConfig}
              onChange={setContamConfig}
              onApply={handleApplyFilter}
              applying={applying}
              stats={contamStats}
            />
          )}

          {selected && showTable && (
            <DataTable
              key={`${source}-${selected}-${tableKey}`}
              fetcher={fetcher}
              distinctFetcher={distinctFetcher}
              cellRenderer={cellRenderer}
              extraRowActions={extraRowActions}
            />
          )}
        </>
      )}

      {addFuncdbPrefill !== null && (
        <AddFuncdbModal
          prefill={addFuncdbPrefill}
          defaultModifiedBy={defaultModifiedBy}
          onClose={() => setAddFuncdbPrefill(null)}
        />
      )}
    </div>
  )
}
