// © 2026 Joshua Benjamin Jewell. All rights reserved.
// Licensed under the GNU Affero General Public License version 3 (AGPLv3).
import { useState, useCallback, useEffect, useRef } from 'react'
import { api } from '../../api/client'
import { errorMessage } from '../../api/errorMessage'
import { DataTable } from '../../components/DataTable'
import { useToast } from '../../components/Toast'
import { NameDialog } from '../../components/NameDialog'
import type { TableMeta, TableQuery, ColFilter, FilterPreset, TableDisplay } from '../../api/types'
import type { RowPopupData, TableStats } from '../../components/DataTable'
import { parseFilterYaml } from './parseFilterYaml'

export function TablesPanel({ study, run, group, tables, onTablesChanged, cacheKey, selected, setSelected, filters, setFilters }: {
  study: string; run: string; group?: string; tables: TableMeta[]; onTablesChanged: () => void; cacheKey?: string | null
  selected: string | null; setSelected: React.Dispatch<React.SetStateAction<string | null>>
  filters: Record<string, ColFilter>; setFilters: React.Dispatch<React.SetStateAction<Record<string, ColFilter>>>
}) {
  const [sortBy, setSortBy]     = useState<string | null>(null)
  const [sortDir, setSortDir]   = useState<'asc' | 'desc'>('asc')
  const [display, setDisplay]   = useState<TableDisplay | undefined>(undefined)
  const [filterKey, setFilterKey] = useState(0)
  const [liveStats, setLiveStats] = useState<TableStats | null>(null)
  const [presets, setPresets]     = useState<FilterPreset[]>([])
  const [loadingPreset, setLoadingPreset] = useState(false)
  const [saving, setSaving]       = useState(false)
  const [nameDialog, setNameDialog] = useState<'filters' | 'table' | null>(null)
  const [managing, setManaging]   = useState(false)
  const toast = useToast()
  const fileInputRef = useRef<HTMLInputElement>(null)

  const reloadPresets = useCallback(() => { api.presets.list().then(setPresets).catch(() => {}) }, [])
  useEffect(reloadPresets, [reloadPresets])

  const fetcher = useCallback(
    (q: TableQuery) => api.results.runTable(study, run, selected!, q, group),
    [study, run, selected, group],
  )

  const distinctFetcher = useCallback(
    (column: string, activeFilters?: Record<string, ColFilter>, keywordFilter?: string) =>
      api.results.distinctValues(study, run, selected!, column, activeFilters, group, keywordFilter),
    [study, run, selected, group],
  )

  const otuPopupFetcher = useCallback(
    async (row: Record<string, unknown>): Promise<RowPopupData | null> => {
      const seqName = String(row['SeqName'] ?? '')
      if (!seqName) return null
      try {
        const data = await api.results.otuMembers(study, run, seqName, group)
        return data.rows.length > 0 ? { columns: data.columns, rows: data.rows } : null
      } catch { return null }
    },
    [study, run, group],
  )

  const [otuCellLabels, setOtuCellLabels] = useState<Record<string, Record<string, string>>>({})
  const [mergedColumns, setMergedColumns] = useState<string[]>([])
  useEffect(() => {
    if (selected !== 'merged_otu') { setOtuCellLabels({}); setMergedColumns([]); return }
    let cancelled = false
    api.results.otuCounts(study, run, group).then(({ counts }) => {
      if (cancelled) return
      const labels: Record<string, string> = {}
      for (const [otu, n] of Object.entries(counts)) {
        labels[otu] = `${otu} (${n})`
      }
      setOtuCellLabels({ SeqName: labels })
    }).catch(() => { if (!cancelled) setOtuCellLabels({}) })
    api.results.runTable(study, run, 'merged', { page: 1, perPage: 1 }, group)
      .then(d => { if (!cancelled) setMergedColumns(d.columns) })
      .catch(() => { if (!cancelled) setMergedColumns([]) })
  }, [study, run, selected, group])

  const importFilters = (e: React.ChangeEvent<HTMLInputElement>) => {
    const file = e.target.files?.[0]
    if (!file) return
    const reader = new FileReader()
    reader.onload = () => {
      try {
        const text = reader.result as string
        const parsed = parseFilterYaml(text)
        setFilters(parsed)
        setFilterKey(k => k + 1)
        toast.success(`Imported filters for ${Object.keys(parsed).length} column${Object.keys(parsed).length === 1 ? '' : 's'}`)
      } catch (err) {
        toast.error('Failed to parse filter YAML: ' + (errorMessage(err)))
      }
    }
    reader.readAsText(file)
    e.target.value = ''
  }

  const applyPreset = async (preset: FilterPreset) => {
    if (!selected) return
    setLoadingPreset(true)
    try {
      const result = await api.presets.apply(study, run, selected, preset.file, group)
      const newFilters: Record<string, ColFilter> = {}
      for (const [col, f] of Object.entries(result.filters)) {
        const cf: ColFilter = {}
        if (f.include) cf.include = f.include
        if (f.exclude) cf.exclude = f.exclude
        if (f.min != null) cf.min = f.min
        if (f.max != null) cf.max = f.max
        if (f.basis) cf.basis = f.basis
        newFilters[col] = cf
      }
      setFilters(newFilters)
      setFilterKey(k => k + 1)
      toast.success(`Applied ${preset.label}: ${result.rows_after.toLocaleString()} of ${result.rows_before.toLocaleString()} rows`)
    } catch (err) {
      toast.error('Failed to apply preset: ' + (errorMessage(err)))
    } finally {
      setLoadingPreset(false)
    }
  }

  // Errors are thrown so NameDialog shows them beside the name.
  const saveFilters = async (name: string) => {
    if (presets.some(p => p.name === name) && !confirm(`Replace the preset "${name}"?`))
      throw new Error(`"${name}" is already a preset.`)
    setSaving(true)
    try {
      await api.presets.save(name, filters, `Saved from table ${selected}`)
      reloadPresets()
      setNameDialog(null)
      toast.success(`Saved filters as ${name}`)
    } finally {
      setSaving(false)
    }
  }

  const saveTable = async (name: string) => {
    if (tables.some(t => t.id === name)) throw new Error(`A table named "${name}" already exists.`)
    setSaving(true)
    try {
      const result = await api.results.saveTable(study, run, selected!, name,
        filters, sortBy ?? undefined, sortDir, group)
      onTablesChanged()
      setNameDialog(null)
      toast.success(`Saved ${result.rows} rows to ${result.name}.csv`)
    } finally {
      setSaving(false)
    }
  }

  const deleteTable = async (id: string) => {
    if (!confirm(`Delete table "${id}"?`)) return
    try {
      await api.results.deleteTable(study, run, id, group)
      if (selected === id) setSelected(tables.find(t => t.id !== id)?.id ?? null)
      onTablesChanged()
      toast.success(`Deleted table ${id}`)
    } catch (err) {
      toast.error('Failed to delete table: ' + (errorMessage(err)))
    }
  }

  const deletePreset = async (preset: FilterPreset) => {
    if (!confirm(`Delete filter preset "${preset.label}"?`)) return
    try {
      await api.presets.delete(preset.file)
      reloadPresets()
      toast.success(`Deleted preset ${preset.label}`)
    } catch (err) {
      toast.error('Failed to delete preset: ' + (errorMessage(err)))
    }
  }

  const exportTable = async () => {
    try {
      await api.results.exportTable(study, run, selected!,
        filters, sortBy ?? undefined, sortDir, group, display)
    } catch (err) {
      toast.error('Failed to export table: ' + (errorMessage(err)))
    }
  }

  if (tables.length === 0) return <div className="empty-state">No tables generated yet.</div>

  return (
    <>
      <div style={{ display: 'flex', gap: 8, flexWrap: 'wrap', marginBottom: 8 }}>
        {tables.map(t => (
          <span key={t.id} style={{ display: 'inline-flex', alignItems: 'center' }}>
            <button
              className={`btn ${selected === t.id ? 'btn-primary' : ''}`}
              onClick={() => { setSelected(t.id); setFilters({}); setSortBy(null); setSortDir('asc'); setFilterKey(k => k + 1); setLiveStats(null) }}
            >
              {t.label} ({t.rows})
            </button>
            {t.id !== 'merged' && t.id !== 'merged_otu' && (
              <button type="button" aria-label={`Delete ${t.label}`} title={`Delete ${t.id}`}
                style={{ marginLeft: 2, border: 'none', background: 'none', cursor: 'pointer', color: 'var(--color-muted-fg)', font: 'inherit', padding: '0 4px' }}
                onClick={() => deleteTable(t.id)}
              >&times;</button>
            )}
          </span>
        ))}
      </div>

      {selected && (() => {
        const meta = tables.find(t => t.id === selected)
        if (!meta) return null
        const copy = (v: string | number) => () => {
          navigator.clipboard?.writeText(String(v)).then(() => toast.info('Copied')).catch(() => {})
        }
        const N = ({ v, raw }: { v: string | number; raw?: string | number }) => (
          <strong onClick={copy(raw ?? v)} style={{ color: 'var(--color-fg)', cursor: 'pointer' }} title="Click to copy">
            {typeof v === 'number' ? v.toLocaleString() : v}
          </strong>
        )
        const fmt = (n: number) => n.toLocaleString()
        const s = liveStats
        const isFiltered = s != null && (s.total !== s.total_unfiltered || s.samples !== s.samples_unfiltered)
        const hasReads   = s != null && s.total_reads_unfiltered > 0
        const share      = (part: number, whole: number) => whole > 0 ? Math.round(part / whole * 100) : null
        const pct        = hasReads && isFiltered ? share(s.total_reads, s.total_reads_unfiltered) : null
        const Pct = ({ p }: { p: number | null }) => p != null ? <> (<N v={`${p}%`} raw={p} />)</> : null
        const avg        = hasReads && meta.n_samples > 0
          ? Math.round((isFiltered ? s.total_reads : s!.total_reads_unfiltered) / meta.n_samples)
          : null

        const parts: React.ReactNode[] = []

        if (s != null && isFiltered && s.samples_unfiltered > 0)
          parts.push(<span key="samples"><N v={s.samples} /> / <N v={s.samples_unfiltered} /> samples<Pct p={share(s.samples, s.samples_unfiltered)} /></span>)
        else if (meta.n_samples > 0)
          parts.push(<span key="samples"><N v={meta.n_samples} /> {meta.n_samples === 1 ? 'sample' : 'samples'}</span>)

        if (s != null) {
          parts.push(
            isFiltered
              ? <span key="rows"><N v={s.total} /> / <N v={s.total_unfiltered} /> rows<Pct p={share(s.total, s.total_unfiltered)} /></span>
              : <span key="rows"><N v={s.total_unfiltered} /> rows</span>
          )
          if (hasReads) {
            parts.push(
              isFiltered
                ? <span key="reads"><N v={fmt(s.total_reads)} raw={s.total_reads} /> / <N v={fmt(s.total_reads_unfiltered)} raw={s.total_reads_unfiltered} /> reads<Pct p={pct} /></span>
                : <span key="reads"><N v={fmt(s.total_reads_unfiltered)} raw={s.total_reads_unfiltered} /> reads</span>
            )
          }
          if (avg != null)
            parts.push(<span key="avg"><N v={fmt(avg)} raw={avg} /> avg reads/sample</span>)
        }

        return parts.length > 0 ? (
          <div style={{ marginBottom: 10, fontSize: '.82rem', color: 'var(--color-muted-fg)', display: 'flex', gap: 16, flexWrap: 'wrap', alignItems: 'center' }}>
            {parts}
          </div>
        ) : null
      })()}

      {selected && (
        <div style={{ display: 'flex', gap: 8, alignItems: 'center', flexWrap: 'wrap', marginBottom: 12, fontSize: '.82rem' }}>
          {presets.length > 0 && (
            <select
              aria-label="Apply a filter preset"
              style={{ padding: '4px 8px', borderRadius: 4, border: '1px solid var(--color-border)', fontSize: '.82rem', background: 'var(--color-bg)' }}
              value=""
              onChange={e => {
                const p = presets.find(p => p.file === e.target.value)
                if (p) applyPreset(p)
              }}
              disabled={loadingPreset}
            >
              <option value="">{loadingPreset ? 'Applying…' : 'Apply preset…'}</option>
              {presets.map(p => (
                <option key={p.file} value={p.file} title={p.description}>{p.label}</option>
              ))}
            </select>
          )}
          {presets.length > 0 && (
            <span style={{ position: 'relative' }}>
              <button className="btn btn-sm" aria-expanded={managing} onClick={() => setManaging(m => !m)}>Manage presets</button>
              {managing && (
                <div role="dialog" aria-label="Filter presets" style={{
                  position: 'absolute', top: '100%', left: 0, zIndex: 20, marginTop: 4, minWidth: 240,
                  background: 'var(--color-bg)', border: '1px solid var(--color-border)', borderRadius: 6,
                  boxShadow: '0 4px 16px rgba(0,0,0,.12)', padding: 8,
                }} onKeyDown={e => { if (e.key === 'Escape') setManaging(false) }}>
                  {presets.map(p => (
                    <div key={p.file} style={{ display: 'flex', alignItems: 'center', justifyContent: 'space-between', gap: 8, padding: '2px 0' }}>
                      <span title={p.description}>{p.label}</span>
                      <button className="btn btn-sm btn-danger" aria-label={`Delete preset ${p.label}`}
                        onClick={() => deletePreset(p)}>Delete</button>
                    </div>
                  ))}
                </div>
              )}
            </span>
          )}
          <button className="btn btn-sm"
            onClick={() => setNameDialog('filters')} disabled={saving}>Save filters</button>
          <button className="btn btn-sm"
            onClick={() => fileInputRef.current?.click()}>Import filters</button>
          <input ref={fileInputRef} type="file" accept=".yml,.yaml" style={{ display: 'none' }}
            onChange={importFilters} />
          <button className="btn btn-sm"
            onClick={() => setNameDialog('table')} disabled={saving}>Save table</button>
          <button className="btn btn-sm"
            onClick={exportTable}>Export table</button>
        </div>
      )}

      {nameDialog === 'filters' && (
        <NameDialog title="Save filters as preset" placeholder="eukaryotes_only"
          onConfirm={saveFilters} onClose={() => setNameDialog(null)} />
      )}
      {nameDialog === 'table' && (
        <NameDialog title="Save filtered table" placeholder="table_name"
          onConfirm={saveTable} onClose={() => setNameDialog(null)} />
      )}

      {selected && (
        <>
          <DataTable
            key={`tbl:${study}/${run}/${group ?? ''}/${selected}:${filterKey}`}
            storageKey={`tbl:${study}/${run}/${group ?? ''}/${selected}`}
            fetcher={fetcher}
            refreshKey={cacheKey}
            distinctFetcher={distinctFetcher}
            rowPopupFetcher={selected === 'merged_otu' ? otuPopupFetcher : undefined}
            popupColumns={selected === 'merged_otu' ? mergedColumns : undefined}
            cellLabels={selected === 'merged_otu' ? otuCellLabels : undefined}
            initialFilters={filters}
            onFiltersChange={setFilters}
            onSortChange={(sb, sd) => { setSortBy(sb); setSortDir(sd) }}
            onStatsChange={setLiveStats}
            onDisplayChange={setDisplay}
            showTaxonomyPresets
          />
        </>
      )}
    </>
  )
}
