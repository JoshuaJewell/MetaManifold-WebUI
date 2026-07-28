// © 2026 Joshua Benjamin Jewell. All rights reserved.
// Licensed under the GNU Affero General Public License version 3 (AGPLv3).
import { useState, useCallback, useEffect, useRef } from 'react'
import { useApi } from '../hooks/useApi'
import type { TablePage, TableQuery, ColFilter, DistinctInfo, HeatmapMode, TableDisplay } from '../api/types'
import { SAMPLE_READS_FILTER_KEY } from '../api/types'
import styles from './DataTable.module.css'
import { ColumnDropdown } from './dataTable/ColumnDropdown'
import { SampleReadsControl } from './dataTable/SampleReadsControl'
import { RowPopup } from './dataTable/RowPopup'
import { blastUrl, flashCopy } from './dataTable/copy'

// Same palette as the server's heatmap.jl: white to #6FA8DC over 0 to max.
const HEAT_HIGH = [0x6f, 0xa8, 0xdc]
function heatColour(v: number, max: number): string | undefined {
  if (!(v !== 0 && max > 0 && Number.isFinite(v))) return undefined
  const t = Math.min(Math.abs(v) / max, 1)
  return '#' + HEAT_HIGH.map(c => Math.round(255 + (c - 255) * t).toString(16).padStart(2, '0')).join('')
}

// Star marker for the row-highlight toggle; filled when the row is highlighted,
// outlined otherwise. Uses currentColor so it follows the button's theme colour.
const StarIcon = ({ filled }: { filled: boolean }) => (
  <svg viewBox="0 0 24 24" width="13" height="13" aria-hidden="true"
    fill={filled ? 'currentColor' : 'none'} stroke="currentColor" strokeWidth="1.5">
    <path d="M12 2.6l2.7 5.9 6.4.6-4.8 4.2 1.4 6.3L12 16.9 6.3 19.6l1.4-6.3L2.9 9.1l6.4-.6z" />
  </svg>
)

export interface RowPopupData {
  columns: string[]
  rows: Record<string, unknown>[]
}

export interface TableStats {
  total:                  number
  total_unfiltered:       number
  total_reads:            number
  total_reads_unfiltered: number
  /** Sample columns with at least one read in the filtered rows. */
  samples:                number
  samples_unfiltered:     number
}

interface Props {
  fetcher: (q: TableQuery) => Promise<TablePage>
  refreshKey?: string | number | null
  /** Stable key for persisting column visibility, sort, and filters to sessionStorage. */
  storageKey?: string
  distinctFetcher?: (column: string, activeFilters?: Record<string, ColFilter>, keywordFilter?: string) => Promise<DistinctInfo>
  rowPopupFetcher?: (row: Record<string, unknown>) => Promise<RowPopupData | null>
  /** Extra columns available only in the popup (e.g. merged table columns not in merged_otu). */
  popupColumns?: string[]
  /** Map from cell value to display label, keyed by column name. E.g. { SeqName: { otu1: 'otu1 (3)' } }. */
  cellLabels?: Record<string, Record<string, string>>
  /** Override cell rendering for specific columns. Return null to fall back to default. */
  cellRenderer?: (column: string, value: string, row: Record<string, unknown>) => React.ReactNode | null
  /** Extra actions rendered in the action cell (same cell as BLAST, or its own cell if no sequence col). */
  extraRowActions?: (row: Record<string, unknown>) => React.ReactNode
  /** Show VSEARCH / DADA2 taxonomy column preset buttons (Tables tab). The All button is always shown when _dada2 cols are present. */
  showTaxonomyPresets?: boolean
  perPage?: number
  initialFilters?: Record<string, ColFilter>
  onFiltersChange?: (filters: Record<string, ColFilter>) => void
  onSortChange?: (sortBy: string | null, sortDir: 'asc' | 'desc') => void
  onStatsChange?: (stats: TableStats | null) => void
  /** Count-cell display options, for passing on to the table export. */
  onDisplayChange?: (display: TableDisplay) => void
}

interface PersistedTableState {
  hiddenCols?: string[]
  stickyCols?: string[]
  sortBy?: string | null
  sortDir?: SortDir
  colFilters?: Record<string, ColFilter>
  activePreset?: 'vsearch' | 'dada2' | null
  highlighted?: string[]
  heatmap?: HeatmapMode
  hideZeros?: boolean
}

function loadPersistedState(key: string): PersistedTableState | null {
  try {
    const raw = sessionStorage.getItem(`dt:${key}`)
    return raw ? JSON.parse(raw) : null
  } catch { return null }
}

function savePersistedState(key: string, state: PersistedTableState) {
  try {
    sessionStorage.setItem(`dt:${key}`, JSON.stringify(state))
  } catch { /* quota exceeded - ignore */ }
}

type SortDir = 'asc' | 'desc'

const isActiveFilter = (f: ColFilter) =>
  f.include != null || f.min != null || f.max != null || (f.exclude?.length ?? 0) > 0

export function DataTable({ fetcher, refreshKey, storageKey, distinctFetcher, rowPopupFetcher, popupColumns, cellLabels, cellRenderer, extraRowActions, showTaxonomyPresets = false, perPage = 100, initialFilters, onFiltersChange, onSortChange, onStatsChange, onDisplayChange }: Props) {
  const [persisted] = useState(() => storageKey ? loadPersistedState(storageKey) : null)
  const [page, setPage]             = useState(1)
  const [filter, setFilter]         = useState('')
  const [sortBy, setSortBy]         = useState<string | null>(persisted?.sortBy ?? null)
  const [sortDir, setSortDir]       = useState<SortDir>(persisted?.sortDir ?? 'asc')
  // A preset or import remounts the table with non-empty initialFilters, which win over the session copy.
  const [colFilters, _setColFilters] = useState<Record<string, ColFilter>>(
    initialFilters && Object.keys(initialFilters).length ? initialFilters : persisted?.colFilters ?? {},
  )
  const [openDropdown, setOpenDropdown] = useState<string | null>(null)
  const [hiddenCols, setHiddenCols] = useState<Set<string>>(new Set(persisted?.hiddenCols))
  const [stickyCols, setStickyCols] = useState<Set<string>>(new Set(persisted?.stickyCols))
  const [showColPicker, setShowColPicker] = useState(false)
  const [activePreset, setActivePreset] = useState<'vsearch' | 'dada2' | null>(persisted?.activePreset ?? null)
  // Manually highlighted rows, keyed by stable row identity so a highlight
  // survives filtering, sorting, and paging, and reappears when a filter clears.
  const [highlighted, setHighlighted] = useState<Set<string>>(new Set(persisted?.highlighted))
  const [heatmapOn, setHeatmapOn]       = useState(persisted?.heatmap != null && persisted.heatmap !== 'none')
  const [heatmapWhole, setHeatmapWhole] = useState(persisted?.heatmap === 'table')
  const [hideZeros, setHideZeros]       = useState(persisted?.hideZeros ?? false)
  const heatmap: HeatmapMode = heatmapOn ? (heatmapWhole ? 'table' : 'column') : 'none'
  const colPickerRef = useRef<HTMLDivElement>(null)

  const [popupData, setPopupData]       = useState<RowPopupData | null>(null)
  const [popupLoading, setPopupLoading] = useState(false)
  const [popupPos, setPopupPos]         = useState<{ x: number; y: number }>({ x: 0, y: 0 })
  const [popupRowIdx, setPopupRowIdx]   = useState<number | null>(null)
  const popupTimer = useRef<ReturnType<typeof setTimeout> | null>(null)
  const popupRef   = useRef<HTMLDivElement>(null)
  // Counts popup requests so a slow response for an earlier row is ignored.
  const popupSeq   = useRef(0)

  const startPopup = useCallback((row: Record<string, unknown>, rowIdx: number, e: React.MouseEvent) => {
    if (!rowPopupFetcher) return
    if (popupTimer.current) clearTimeout(popupTimer.current)
    const rect = (e.currentTarget as HTMLElement).getBoundingClientRect()
    popupTimer.current = setTimeout(() => {
      setPopupPos({ x: rect.left, y: Math.min(rect.bottom + 4, window.innerHeight - 340) })
      setPopupRowIdx(rowIdx)
      setPopupLoading(true)
      setPopupData(null)
      const seq = ++popupSeq.current
      rowPopupFetcher(row).then(data => {
        if (seq !== popupSeq.current) return
        setPopupData(data)
        setPopupLoading(false)
      }).catch(() => { if (seq === popupSeq.current) { setPopupData(null); setPopupLoading(false) } })
    }, 300)
  }, [rowPopupFetcher])

  const cancelTimer = useRef<ReturnType<typeof setTimeout> | null>(null)

  const cancelPopup = useCallback(() => {
    if (popupTimer.current) { clearTimeout(popupTimer.current); popupTimer.current = null }
    if (cancelTimer.current) clearTimeout(cancelTimer.current)
    cancelTimer.current = setTimeout(() => {
      popupSeq.current++
      setPopupData(null)
      setPopupLoading(false)
      setPopupRowIdx(null)
    }, 150)
  }, [])

  const keepPopup = useCallback(() => {
    if (cancelTimer.current) { clearTimeout(cancelTimer.current); cancelTimer.current = null }
  }, [])

  const setColFilters = useCallback((updater: Record<string, ColFilter> | ((prev: Record<string, ColFilter>) => Record<string, ColFilter>)) => {
    _setColFilters(prev => typeof updater === 'function' ? updater(prev) : updater)
  }, [])

  useEffect(() => {
    onFiltersChange?.(colFilters)
  }, [colFilters]) // eslint-disable-line react-hooks/exhaustive-deps

  // Notify parent of restored sort state on mount.
  useEffect(() => {
    if (persisted?.sortBy !== undefined) onSortChange?.(sortBy, sortDir)
  }, []) // eslint-disable-line react-hooks/exhaustive-deps

  useEffect(() => {
    if (initialFilters && !persisted?.colFilters) {
      _setColFilters(initialFilters)
      onFiltersChange?.(initialFilters)
      setPage(1)
    }
  }, [initialFilters]) // eslint-disable-line react-hooks/exhaustive-deps

  // Persist table UI state to sessionStorage on change.
  useEffect(() => {
    if (!storageKey) return
    savePersistedState(storageKey, {
      hiddenCols: [...hiddenCols],
      stickyCols: [...stickyCols],
      sortBy,
      sortDir,
      colFilters,
      activePreset,
      highlighted: [...highlighted],
      heatmap,
      hideZeros,
    })
  }, [storageKey, hiddenCols, stickyCols, sortBy, sortDir, colFilters, activePreset, highlighted, heatmap, hideZeros])

  useEffect(() => {
    onDisplayChange?.({ heatmap, hide_zeros: hideZeros })
  }, [heatmap, hideZeros]) // eslint-disable-line react-hooks/exhaustive-deps

  const bound = useCallback(() => {
    const q: TableQuery = { page, perPage }
    if (filter) q.filter = filter
    if (sortBy) { q.sortBy = sortBy; q.sortDir = sortDir }
    const active: Record<string, ColFilter> = {}
    for (const [col, f] of Object.entries(colFilters)) {
      if (isActiveFilter(f)) active[col] = f
    }
    if (Object.keys(active).length > 0) q.colFilters = active
    return fetcher(q)
  }, [fetcher, page, perPage, filter, sortBy, sortDir, colFilters])

  const { data, loading, error, refetch } = useApi(bound)

  const prevRefreshKey = useRef(refreshKey)
  useEffect(() => {
    if (refreshKey == null || refreshKey === prevRefreshKey.current) return
    prevRefreshKey.current = refreshKey
    refetch()
  }, [refreshKey, refetch])

  useEffect(() => {
    if (!data) { onStatsChange?.(null); return }
    onStatsChange?.({
      total: data.total, total_unfiltered: data.total_unfiltered,
      total_reads: data.total_reads, total_reads_unfiltered: data.total_reads_unfiltered,
      samples: Object.values(data.count_max ?? {}).filter(v => v > 0).length,
      samples_unfiltered: data.sample_count_columns.length + (data.excluded_samples?.length ?? 0),
    })
  }, [data]) // eslint-disable-line react-hooks/exhaustive-deps

  const rows  = data && Array.isArray(data.rows) ? data.rows : []
  const allCols = data && Array.isArray(data.columns) ? data.columns
              : (rows.length > 0 && rows[0] ? Object.keys(rows[0]) : [])
  const sampleCountColSet = new Set(data?.sample_count_columns ?? [])
  const allCountsHidden = sampleCountColSet.size > 0 && [...sampleCountColSet].every(c => hiddenCols.has(c))
  const countMax = data?.count_max ?? {}
  const wholeMax = Math.max(0, ...Object.values(countMax))
  const countCellStyle = (c: string, text: string): React.CSSProperties | undefined => {
    if (heatmap === 'none' || !sampleCountColSet.has(c) || text === '') return undefined
    const bg = heatColour(Number(text), heatmap === 'table' ? wholeMax : countMax[c] ?? 0)
    return bg ? { background: bg, color: '#111' } : undefined
  }
  const isHiddenZero = (c: string, text: string) =>
    hideZeros && sampleCountColSet.has(c) && text !== '' && Number(text) === 0
  const cols = allCols.filter(c => !hiddenCols.has(c))
  const hasSequenceCol = allCols.includes('sequence')
  const pages = data ? Math.ceil(data.total / perPage) : 0

  const tableColSet = new Set(allCols)
  const extraPopupCols = (popupColumns ?? []).filter(c => !tableColSet.has(c))
  const pickerCols = [...allCols, ...extraPopupCols]
  const popupOnlySet = new Set(extraPopupCols)
  // Count visible columns from the current set directly. Deriving it by
  // subtraction (pickerCols.length - hiddenCols.size) underflows when hiddenCols
  // retains columns absent from the present view, e.g. after switching runs.
  const visibleColCount = pickerCols.filter(c => !hiddenCols.has(c)).length

  useEffect(() => {
    if (!showColPicker) return
    const handler = (e: MouseEvent) => {
      if (colPickerRef.current && !colPickerRef.current.contains(e.target as Node)) setShowColPicker(false)
    }
    document.addEventListener('mousedown', handler)
    return () => document.removeEventListener('mousedown', handler)
  }, [showColPicker])

  const handleSort = (col: string) => {
    let newSortBy: string | null
    let newSortDir: SortDir
    if (sortBy === col) {
      if (sortDir === 'asc') { newSortBy = col; newSortDir = 'desc' }
      else { newSortBy = null; newSortDir = 'asc' }
    } else {
      newSortBy = col; newSortDir = 'asc'
    }
    setSortBy(newSortBy)
    setSortDir(newSortDir)
    onSortChange?.(newSortBy, newSortDir)
    setPage(1)
  }

  const updateColFilter = (col: string, patch: ColFilter | undefined) => {
    setColFilters(prev => {
      const next = { ...prev }
      if (!patch) delete next[col]
      else next[col] = patch
      return next
    })
    setPage(1)
  }

  const clearAllFilters = () => { setColFilters({}); setFilter(''); setPage(1) }

  const sortIndicator = (col: string) => {
    if (sortBy !== col) return <span className={styles.sortIcon}> +</span>
    return <span className={styles.sortIconActive}>{sortDir === 'asc' ? ' ^' : ' v'}</span>
  }

  const hasAnyFilter = !!filter || Object.keys(colFilters).length > 0
  const colIsFiltered = (col: string) => {
    const f = colFilters[col]
    return !!f && isActiveFilter(f)
  }

  const toggleColVisibility = (col: string) => {
    setHiddenCols(prev => {
      const next = new Set(prev)
      if (next.has(col)) next.delete(col); else next.add(col)
      return next
    })
  }

  const toggleSticky = (col: string) => {
    setStickyCols(prev => {
      const next = new Set(prev)
      if (next.has(col)) next.delete(col); else next.add(col)
      return next
    })
  }


  const applyColumnPreset = (preset: 'vsearch' | 'dada2' | 'all') => {
    if (preset === 'all') {
      setHiddenCols(new Set())
      setColFilters({})
      setActivePreset(null)
      setPage(1)
      return
    }
    const vsearchOnlyCols = ['Pident', 'Accession', 'rRNA', 'Organellum', 'specimen']
    const vsearchTaxCols = pickerCols.filter(c => pickerCols.includes(c + '_dada2'))
    const dada2Cols = pickerCols.filter(c => c.endsWith('_dada2') || c.endsWith('_boot'))
    if (preset === 'vsearch') {
      setHiddenCols(new Set(dada2Cols))
      const filters: Record<string, ColFilter> = {}
      if (pickerCols.includes('Pident')) filters['Pident'] = { min: 0 }
      setColFilters(filters)
    } else {
      setHiddenCols(new Set([...vsearchOnlyCols.filter(c => pickerCols.includes(c)), ...vsearchTaxCols]))
      const filters: Record<string, ColFilter> = {}
      const bootCol = pickerCols.find(c => c.endsWith('_boot'))
      if (bootCol) filters[bootCol] = { min: 0 }
      setColFilters(filters)
    }
    setActivePreset(preset)
    setPage(1)
  }

  // Stable identity for a row. A natural key keeps highlights on the same row
  // through filtering and sorting.
  const rowKey = (row: Record<string, unknown>): string => {
    const id = row['SeqName'] ?? row['sequence']
    return id != null ? String(id) : JSON.stringify(row)
  }
  const toggleHighlight = (key: string) => setHighlighted(prev => {
    const next = new Set(prev)
    if (next.has(key)) next.delete(key); else next.add(key)
    return next
  })

  // The leading highlight handle is a sticky column at left:0; data sticky
  // columns begin after its width.
  const HANDLE_W = 34
  const stickyOffsets = new Map<string, number>()
  {
    let offset = HANDLE_W
    for (const c of cols) {
      if (stickyCols.has(c)) {
        stickyOffsets.set(c, offset)
        offset += 140
      }
    }
  }

  return (
    <div className={styles.wrapper}>
      <div className={styles.toolbar}>
        <input
          className={styles.search}
          placeholder="Global filter…"
          aria-label="Filter all columns"
          value={filter}
          onChange={e => { setFilter(e.target.value); setPage(1) }}
        />
        <SampleReadsControl
          current={colFilters[SAMPLE_READS_FILTER_KEY]}
          excluded={data?.excluded_samples ?? []}
          onApply={f => updateColFilter(SAMPLE_READS_FILTER_KEY, f)}
        />
        {hasAnyFilter && (
          <button className="btn" style={{ fontSize: '.78rem', padding: '3px 8px' }}
            onClick={clearAllFilters}>Clear all filters</button>
        )}
        {pickerCols.length > 0 && (
          <div style={{ position: 'relative' }}>
            <button className="btn" style={{ fontSize: '.78rem', padding: '3px 8px' }}
              onClick={() => setShowColPicker(v => !v)}>
              Columns{visibleColCount < pickerCols.length ? ` (${visibleColCount}/${pickerCols.length})` : ''}
            </button>
            {showColPicker && (
              <div ref={colPickerRef} className={styles.dropdown}
                style={{ right: 0, left: 'auto', maxHeight: 320, overflowY: 'auto' }}>
                <div className={styles.dropdownActions}>
                  <button onClick={() => setHiddenCols(new Set())}>Show all</button>
                  <button onClick={() => setHiddenCols(new Set(pickerCols))}>Hide all</button>
                </div>
                <div className={styles.dropdownList}>
                  {pickerCols.map(c => (
                    <label key={c} className={styles.dropdownItem}>
                      <input type="checkbox" checked={!hiddenCols.has(c)}
                        onChange={() => toggleColVisibility(c)} />
                      <span>{c}{popupOnlySet.has(c) ? <span className={styles.popupOnlyTag}> (ASV)</span> : ''}</span>
                    </label>
                  ))}
                </div>
              </div>
            )}
          </div>
        )}
        {pickerCols.some(c => c.endsWith('_dada2')) && (
          <>
            {showTaxonomyPresets && (
              <>
                <button className={`btn${activePreset === 'vsearch' ? ' btn-primary' : ''}`} style={{ fontSize: '.78rem', padding: '3px 8px' }}
                  onClick={() => applyColumnPreset('vsearch')}>VSEARCH</button>
                <button className={`btn${activePreset === 'dada2' ? ' btn-primary' : ''}`} style={{ fontSize: '.78rem', padding: '3px 8px' }}
                  onClick={() => applyColumnPreset('dada2')}>DADA2</button>
              </>
            )}
            <button className={`btn${activePreset === null ? ' btn-primary' : ''}`} style={{ fontSize: '.78rem', padding: '3px 8px' }}
              onClick={() => applyColumnPreset('all')}>All</button>
          </>
        )}
        {sampleCountColSet.size > 0 && (
            <button
              className={`btn${allCountsHidden ? ' btn-primary' : ''}`}
              style={{ fontSize: '.78rem', padding: '3px 8px' }}
              onClick={() => setHiddenCols(prev => {
                const next = new Set(prev)
                if (allCountsHidden) sampleCountColSet.forEach(c => next.delete(c))
                else sampleCountColSet.forEach(c => next.add(c))
                return next
              })}
              title={allCountsHidden ? 'Show sample counts' : 'Hide sample counts'}
            >
              {allCountsHidden ? 'Show counts' : 'Hide counts'}
            </button>
        )}
        {sampleCountColSet.size > 0 && (
          <button
            className={`btn${hideZeros ? ' btn-primary' : ''}`}
            style={{ fontSize: '.78rem', padding: '3px 8px' }}
            aria-pressed={hideZeros}
            onClick={() => setHideZeros(h => !h)}
            title="Blank zero counts"
          >Hide zeros</button>
        )}
        {sampleCountColSet.size > 0 && (
          <button
            className={`btn${heatmapOn ? ' btn-primary' : ''}`}
            style={{ fontSize: '.78rem', padding: '3px 8px' }}
            aria-pressed={heatmapOn}
            onClick={() => setHeatmapOn(h => !h)}
            title="Shade counts from zero to the largest value"
          >Heatmap</button>
        )}
        {sampleCountColSet.size > 0 && heatmapOn && (
          <button
            className="btn"
            style={{ fontSize: '.78rem', padding: '3px 8px' }}
            onClick={() => setHeatmapWhole(w => !w)}
            title={heatmapWhole ? 'One scale across all sample columns' : 'Each sample column on its own scale'}
          >{heatmapWhole ? 'Scale: whole table' : 'Scale: per column'}</button>
        )}
        {highlighted.size > 0 && (
          <button
            className="btn"
            style={{ fontSize: '.78rem', padding: '3px 8px' }}
            onClick={() => setHighlighted(new Set())}
            title="Clear all highlighted rows"
          >Clear highlights ({highlighted.size})</button>
        )}
      </div>

      {error && <p className={styles.error}>{error}</p>}
      {loading && cols.length === 0 && <p className={styles.msg}>Loading…</p>}
      {!loading && !error && cols.length === 0 && <p className={styles.msg}>No data.</p>}

      {cols.length > 0 && (
        <>
          <div className={styles.scroll}>
            <table className={styles.table}>
              <thead>
                <tr>
                  <th className={styles.handleCol}
                    style={{ position: 'sticky', left: 0, zIndex: 12, width: HANDLE_W, background: 'var(--color-surface)' }}
                    title="Highlight rows"
                  ></th>
                  {cols.map(c => {
                    const isSticky = stickyCols.has(c)
                    const stickyStyle: React.CSSProperties | undefined = isSticky ? {
                      position: 'sticky',
                      left: stickyOffsets.get(c) ?? 0,
                      zIndex: 11,
                      background: 'var(--color-surface)',
                    } : undefined
                    return (
                      <th key={c} className={styles.sortable}
                        aria-sort={sortBy === c ? (sortDir === 'asc' ? 'ascending' : 'descending') : undefined}
                        onClick={() => handleSort(c)}
                        style={stickyStyle}>
                        <button type="button" className={styles.headerLabel}
                          onClick={e => { e.stopPropagation(); handleSort(c) }}>
                          {c}{sortIndicator(c)}
                        </button>
                        {distinctFetcher && (
                          <button
                            className={`${styles.dropdownBtn} ${colIsFiltered(c) ? styles.dropdownBtnActive : ''}`}
                            onClick={e => { e.stopPropagation(); setOpenDropdown(openDropdown === c ? null : c) }}
                            title="Filter values"
                            aria-label={`Filter ${c}`}
                          >▾</button>
                        )}
                        {openDropdown === c && distinctFetcher && (
                          <ColumnDropdown
                            column={c}
                            distinctFetcher={distinctFetcher}
                            activeFilters={colFilters}
                            keywordFilter={filter}
                            current={colFilters[c]}
                            isSticky={isSticky}
                            onToggleSticky={() => toggleSticky(c)}
                            onApply={f => { updateColFilter(c, f); setOpenDropdown(null) }}
                            onClose={() => setOpenDropdown(null)}
                          />
                        )}
                      </th>
                    )
                  })}
                  {(hasSequenceCol || extraRowActions) && <th style={{ width: hasSequenceCol && extraRowActions ? 90 : 60 }}></th>}
                </tr>
              </thead>
              <tbody>
                {loading && (
                  <tr><td colSpan={cols.length + 1 + (hasSequenceCol || extraRowActions ? 1 : 0)} className={styles.msg} style={{ textAlign: 'center' }}>Loading…</td></tr>
                )}
                {!loading && rows.length === 0 && (
                  <tr><td colSpan={cols.length + 1 + (hasSequenceCol || extraRowActions ? 1 : 0)} className={styles.msg} style={{ textAlign: 'center' }}>
                    {hasAnyFilter ? 'No matching rows.' : 'No data.'}
                  </td></tr>
                )}
                {!loading && rows.map((row, i) => {
                  const rk = rowKey(row)
                  const isHighlighted = highlighted.has(rk)
                  return (
                  <tr key={i}
                    onMouseEnter={rowPopupFetcher ? e => { keepPopup(); startPopup(row, i, e) } : undefined}
                    onMouseLeave={rowPopupFetcher ? cancelPopup : undefined}
                    className={[popupRowIdx === i ? styles.popupActiveRow : '', isHighlighted ? styles.highlightRow : ''].filter(Boolean).join(' ') || undefined}
                  >
                    <td className={styles.handleCol}
                      style={{ position: 'sticky', left: 0, zIndex: 1, width: HANDLE_W, background: 'var(--color-bg)' }}
                    >
                      <button
                        className={`${styles.highlightBtn}${isHighlighted ? ' ' + styles.highlightBtnOn : ''}`}
                        onClick={e => { e.stopPropagation(); toggleHighlight(rk) }}
                        title={isHighlighted ? 'Remove highlight' : 'Highlight row'}
                      ><StarIcon filled={isHighlighted} /></button>
                    </td>
                    {cols.map(c => {
                      const isSticky = stickyCols.has(c)
                      const stickyStyle: React.CSSProperties | undefined = isSticky ? {
                        position: 'sticky',
                        left: stickyOffsets.get(c) ?? 0,
                        zIndex: 1,
                        background: 'var(--color-bg)',
                      } : undefined
                      const text = String(row[c] ?? '')
                      const custom = cellRenderer?.(c, text, row)
                      if (custom !== null && custom !== undefined) {
                        return <td key={c} style={stickyStyle}>{custom}</td>
                      }
                      const display = isHiddenZero(c, text) ? '' : cellLabels?.[c]?.[text] ?? text
                      const heat = countCellStyle(c, text)
                      return (
                        <td key={c} style={heat ? { ...stickyStyle, ...heat } : stickyStyle}
                          onClick={flashCopy(text)}
                          title="Click to copy"
                        >{display}</td>
                      )
                    })}
                    {(hasSequenceCol || extraRowActions) && (
                      <td className={styles.blastCell}>
                        {hasSequenceCol && (
                          <a
                            href={blastUrl(String(row['sequence'] ?? ''))}
                            target="_blank"
                            rel="noopener noreferrer"
                            className={styles.blastLink}
                            title="Search this sequence on NCBI BLAST"
                          >BLAST</a>
                        )}
                        {extraRowActions?.(row)}
                      </td>
                    )}
                  </tr>
                  )
                })}
              </tbody>
            </table>
          </div>
          {popupRowIdx !== null && (popupLoading || (popupData && popupData.rows.length > 0)) && (
            <RowPopup ref={popupRef} data={popupData} loading={popupLoading} hiddenCols={hiddenCols}
              pos={popupPos} onMouseEnter={keepPopup} onMouseLeave={cancelPopup} />
          )}
          {pages > 1 && (
            <div className={styles.pager}>
              <button disabled={page <= 1} onClick={() => setPage(p => p - 1)}>{'< Prev'}</button>
              <span>Page {page} / {pages}</span>
              <button disabled={page >= pages} onClick={() => setPage(p => p + 1)}>{'Next >'}</button>
            </div>
          )}
        </>
      )}
    </div>
  )
}
