// © 2026 Joshua Benjamin Jewell. All rights reserved.
// Licensed under the GNU Affero General Public License version 3 (AGPLv3).
import { useState } from 'react'
import styles from '../DataTable.module.css'

export function TextFilter({ values, current, onApply }: {
  values:   string[]
  current?: string[]
  onApply:  (vals: string[] | undefined) => void
}) {
  const [search, setSearch]   = useState('')
  const [checked, setChecked] = useState<Set<string>>(() => new Set(current ?? values))

  const visible = search
    ? values.filter(v => v.toLowerCase().includes(search.toLowerCase()))
    : values

  const allTicked = checked.size === values.length

  const toggle = (val: string) => {
    setChecked(prev => {
      const next = new Set(prev)
      if (next.has(val)) next.delete(val); else next.add(val)
      return next
    })
  }

  return (
    <>
      <input
        className={styles.dropdownSearch}
        placeholder="Search…"
        aria-label="Search values"
        value={search}
        onChange={e => setSearch(e.target.value)}
        autoFocus
      />
      <div className={styles.dropdownActions}>
        <button onClick={() => setChecked(new Set(values))}>All</button>
        <button onClick={() => setChecked(new Set())}>None</button>
        <span className={styles.dropdownCount}>{checked.size}/{values.length}</span>
      </div>
      <div className={styles.dropdownList}>
        {visible.map(v => (
          <label key={v} className={styles.dropdownItem}>
            <input type="checkbox" checked={checked.has(v)} onChange={() => toggle(v)} />
            <span>{v || '(empty)'}</span>
          </label>
        ))}
        {visible.length === 0 && <div className={styles.dropdownEmpty}>No matching values</div>}
      </div>
      <div className={styles.dropdownFooter}>
        <button className="btn btn-primary" style={{ fontSize: '.78rem', padding: '3px 10px' }}
          onClick={() => onApply(allTicked ? undefined : [...checked])}>Apply</button>
        <button className="btn" style={{ fontSize: '.78rem', padding: '3px 10px' }}
          onClick={() => onApply(undefined)}>Clear</button>
      </div>
    </>
  )
}
