// © 2026 Joshua Benjamin Jewell. All rights reserved.
// Licensed under the GNU Affero General Public License version 3 (AGPLv3).
import { useState } from 'react'
import styles from '../DataTable.module.css'
import { flashCopy } from './copy'

// A single numeric statistic; clicking it copies the raw (unformatted) value.
function Stat({ value, fmt }: { value: number; fmt: (v: number) => string }) {
  return (
    <span className={styles.statValue} onClick={flashCopy(String(value))} title="Click to copy">
      {fmt(value)}
    </span>
  )
}

export function NumericFilter({ dataMin, dataMax, sum, mean, median, q1, q3, currentMin, currentMax, onApply }: {
  dataMin:     number
  dataMax:     number
  sum?:        number
  mean?:       number
  median?:     number
  q1?:         number
  q3?:         number
  currentMin?: number
  currentMax?: number
  onApply:     (min: number | null, max: number | null) => void
}) {
  const [minVal, setMinVal] = useState(currentMin != null ? String(currentMin) : '')
  const [maxVal, setMaxVal] = useState(currentMax != null ? String(currentMax) : '')

  const fmt = (v: number) => Number.isInteger(v) ? v.toLocaleString() : v.toLocaleString(undefined, { maximumFractionDigits: 2 })

  const apply = () => {
    const mn = minVal !== '' ? Number(minVal) : null
    const mx = maxVal !== '' ? Number(maxVal) : null
    onApply(
      mn != null && !isNaN(mn) ? mn : null,
      mx != null && !isNaN(mx) ? mx : null,
    )
  }

  return (
    <>
      {(sum != null || mean != null || median != null || q1 != null) && (
        <div className={styles.numericInfo} style={{ fontSize: '.75rem', color: 'var(--color-muted-fg)' }}>
          {sum != null && <div>Sum: <Stat value={sum} fmt={fmt} /></div>}
          <div>Range: <Stat value={dataMin} fmt={fmt} /> - <Stat value={dataMax} fmt={fmt} /></div>
          {q1 != null && q3 != null && <div>IQR: <Stat value={q1} fmt={fmt} /> - <Stat value={q3} fmt={fmt} /></div>}
          {mean != null && <div>Mean: <Stat value={mean} fmt={fmt} /></div>}
          {median != null && <div>Median: <Stat value={median} fmt={fmt} /></div>}
        </div>
      )}
      <div className={styles.numericInputs}>
        <label>
          <span>Min</span>
          <input
            type="number"
            className={styles.numericInput}
            placeholder="0"
            value={minVal}
            onChange={e => setMinVal(e.target.value)}
            step="any"
            autoFocus
          />
        </label>
        <label>
          <span>Max</span>
          <input
            type="number"
            className={styles.numericInput}
            placeholder="100"
            value={maxVal}
            onChange={e => setMaxVal(e.target.value)}
            step="any"
          />
        </label>
      </div>
      <div className={styles.dropdownFooter}>
        <button className="btn btn-primary" style={{ fontSize: '.78rem', padding: '3px 10px' }}
          onClick={apply}>Apply</button>
        <button className="btn" style={{ fontSize: '.78rem', padding: '3px 10px' }}
          onClick={() => onApply(null, null)}>Clear</button>
      </div>
    </>
  )
}
