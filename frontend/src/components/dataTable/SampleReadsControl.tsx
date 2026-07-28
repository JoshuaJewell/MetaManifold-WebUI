// © 2026 Joshua Benjamin Jewell. All rights reserved.
// Licensed under the GNU Affero General Public License version 3 (AGPLv3).
import { useState, useEffect } from 'react'
import type { ColFilter, SampleReadsBasis } from '../../api/types'

// Sample total-read filter. Samples are columns, so this drops whole sample
// columns. The bound travels in colFilters under a reserved key so charts, save
// and export honour it too.
export function SampleReadsControl({ current, excluded, onApply }: {
  current?:  ColFilter
  excluded:  string[]
  onApply:   (f: ColFilter | undefined) => void
}) {
  const [min, setMin] = useState(current?.min != null ? String(current.min) : '')
  const [basis, setBasis] = useState<SampleReadsBasis>(current?.basis ?? 'filtered')

  // Follow external changes (Clear all filters, a restored session, a preset).
  useEffect(() => {
    setMin(current?.min != null ? String(current.min) : '')
    setBasis(current?.basis ?? 'filtered')
  }, [current?.min, current?.basis])

  const commit = (nextMin: string, nextBasis: SampleReadsBasis) => {
    const trimmed = nextMin.trim()
    const parsed = trimmed === '' ? null : Number(trimmed)
    if (parsed == null || !Number.isFinite(parsed)) { onApply(undefined); return }
    if (parsed === current?.min && nextBasis === (current?.basis ?? 'filtered')) return
    onApply({ min: parsed, basis: nextBasis })
  }

  return (
    <span style={{ display: 'inline-flex', alignItems: 'center', gap: 4, fontSize: '.78rem' }}
      title="Exclude samples with fewer reads than this.">
      <label htmlFor="dt-sample-min-reads">Min sample reads</label>
      <input
        id="dt-sample-min-reads"
        type="number"
        min={0}
        step={1}
        placeholder="none"
        value={min}
        style={{ width: 70, fontSize: '.78rem', padding: '2px 4px' }}
        onChange={e => setMin(e.target.value)}
        onBlur={() => commit(min, basis)}
        onKeyDown={e => { if (e.key === 'Enter') commit(min, basis) }}
      />
      <select
        aria-label="Sample read count basis"
        value={basis}
        style={{ fontSize: '.78rem', padding: '2px 2px' }}
        onChange={e => {
          const b = e.target.value as SampleReadsBasis
          setBasis(b)
          if (min.trim() !== '') commit(min, b)
        }}
      >
        <option value="filtered">after filters</option>
        <option value="raw">raw</option>
      </select>
      {excluded.length > 0 && (
        <span style={{ color: 'var(--color-muted-fg)' }} title={excluded.join(', ')}>
          ({excluded.length} sample{excluded.length === 1 ? '' : 's'} excluded)
        </span>
      )}
    </span>
  )
}
