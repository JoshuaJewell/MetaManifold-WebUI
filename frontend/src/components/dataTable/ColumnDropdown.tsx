// © 2026 Joshua Benjamin Jewell. All rights reserved.
// Licensed under the GNU Affero General Public License version 3 (AGPLv3).
import { useState, useEffect, useRef } from 'react'
import type { ColFilter, DistinctInfo } from '../../api/types'
import styles from '../DataTable.module.css'
import { TextFilter } from './TextFilter'
import { NumericFilter } from './NumericFilter'

export function ColumnDropdown({ column, distinctFetcher, activeFilters, keywordFilter, current, isSticky, onToggleSticky, onApply, onClose }: {
  column:          string
  distinctFetcher: (col: string, activeFilters?: Record<string, ColFilter>, keywordFilter?: string) => Promise<DistinctInfo>
  activeFilters:   Record<string, ColFilter>
  keywordFilter?:  string
  current?:        ColFilter
  isSticky:        boolean
  onToggleSticky:  () => void
  onApply:         (f: ColFilter | undefined) => void
  onClose:         () => void
}) {
  const [info, setInfo]         = useState<DistinctInfo | null>(null)
  const [loadError, setLoadError] = useState<string | null>(null)
  const ref = useRef<HTMLDivElement>(null)

  useEffect(() => {
    let cancelled = false
    // Pass all filters except this column's so the dropdown shows contextual values.
    const otherFilters: Record<string, ColFilter> = {}
    for (const [col, f] of Object.entries(activeFilters)) {
      if (col !== column) otherFilters[col] = f
    }
    const hasOther = Object.keys(otherFilters).length > 0
    distinctFetcher(column, hasOther ? otherFilters : undefined, keywordFilter || undefined)
      .then(d => { if (!cancelled) setInfo(d) })
      .catch(e => { if (!cancelled) setLoadError(e.message ?? 'Failed to load') })
    return () => { cancelled = true }
  }, [column, distinctFetcher, activeFilters, keywordFilter])

  useEffect(() => {
    const handler = (e: MouseEvent) => {
      if (ref.current && !ref.current.contains(e.target as Node)) onClose()
    }
    document.addEventListener('mousedown', handler)
    return () => document.removeEventListener('mousedown', handler)
  }, [onClose])

  return (
    <div ref={ref} className={styles.dropdown} onClick={e => e.stopPropagation()}>
      <label className={styles.dropdownItem} style={{ borderBottom: '1px solid var(--color-border)', paddingTop: 6, paddingBottom: 6 }}>
        <input type="checkbox" checked={isSticky} onChange={onToggleSticky} />
        <span style={{ fontWeight: 600, fontSize: '.78rem' }}>Sticky column</span>
      </label>
      {loadError && <div className={styles.dropdownError}>{loadError}</div>}
      {!info && !loadError && <div className={styles.dropdownLoading}>Loading…</div>}
      {info?.type === 'text' && (
        <TextFilter
          values={info.values}
          current={current?.include}
          onApply={vals => onApply(vals ? { include: vals } : undefined)}
        />
      )}
      {info?.type === 'numeric' && (
        <NumericFilter
          dataMin={info.min}
          dataMax={info.max}
          sum={info.sum}
          mean={info.mean}
          median={info.median}
          q1={info.q1}
          q3={info.q3}
          currentMin={current?.min}
          currentMax={current?.max}
          onApply={(min, max) => {
            if (min == null && max == null) onApply(undefined)
            else onApply({ min: min ?? undefined, max: max ?? undefined })
          }}
        />
      )}
    </div>
  )
}
