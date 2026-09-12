// © 2026 Joshua Benjamin Jewell. All rights reserved.
// Licensed under the GNU Affero General Public License version 3 (AGPLv3).
import { useEffect, useState } from 'react'
import { indexTree, nodeKey, type DNode, type NodeStyle } from './model'
import styles from './Tree.module.css'
import { fmt } from './ui'
import { StyleEditor } from './controls'

export function Inspector({ d, ix, names, placements, onRename, onStyle, onRotate, onToggle, onReroot, onParent, onQuery }: {
  d: DNode
  ix: ReturnType<typeof indexTree>
  names: Map<string, string>
  placements: { name: string; lwr: number; rank: number; query: number }[]
  onRename: (key: string, name: string | null) => void
  onStyle: (style: NodeStyle | null) => void
  onRotate: () => void
  onToggle: () => void
  onReroot: () => void
  onParent: () => void
  onQuery: (q: number) => void
}) {
  const key = nodeKey(d)
  const original = d.tipKey !== null ? ix.src.nodes[d.src].name : ''
  const [value, setValue] = useState(d.name)
  useEffect(() => setValue(d.name), [d.name, d.clade])
  const commit = () => {
    if (!key) return
    const v = value.trim()
    const next = v === '' || v === original ? null : v
    if (next !== (names.get(key) ?? null)) onRename(key, next)
  }
  const isTip = d.tipKey !== null
  return (
    <div className="card">
      <div className="card-title">{isTip ? 'Tip' : d.parent ? 'Clade' : 'Root'}</div>
      {key && (
        <label className={styles.field}>
          <span>Name</span>
          <input value={value} placeholder={isTip ? original : 'Unnamed clade'}
            onChange={e => setValue(e.target.value)}
            onKeyDown={e => { if (e.key === 'Enter') commit() }}
            onBlur={commit} />
        </label>
      )}
      {key && names.has(key) && (
        <p className={styles.muted}>
          {isTip ? <>In file: {original} </> : null}
          <button className={styles.linkBtn} onClick={() => onRename(key, null)}>Reset name</button>
        </p>
      )}
      <dl className={styles.facts}>
        {!isTip && <><dt>Tips</dt><dd>{d.tips}</dd></>}
        {d.parent && <><dt>Branch length</dt><dd>{fmt(d.len, 6)}</dd></>}
        {d.support !== null && <><dt>Support</dt><dd>{d.support}</dd></>}
        {d.split && <><dt>Branch id</dt><dd className={styles.mono}>{d.split}</dd></>}
      </dl>
      <div className={styles.actions}>
        {!isTip && d.parent && <button className="btn" onClick={onToggle}>{d.collapsed ? 'Expand' : 'Collapse'}</button>}
        {!isTip && !d.collapsed && <button className="btn" onClick={onRotate} title="O">Rotate</button>}
        {d.parent && <button className="btn" onClick={onReroot}>Root on this branch</button>}
        {d.parent && <button className="btn" onClick={onParent}>Select parent</button>}
      </div>
      <div className="card-title" style={{ marginTop: 12 }}>Style</div>
      <StyleEditor style={d.ownStyle} onChange={onStyle} />
      {placements.length > 0 && <>
        <div className="card-title" style={{ marginTop: 12 }}>Placements here</div>
        <ul className={styles.placeList}>
          {placements.map((p, k) => (
            <li key={k}><button className={styles.linkBtn} onClick={() => onQuery(p.query)}>{p.name}</button> LWR {fmt(p.lwr, 3)}{p.rank ? ` (candidate ${p.rank + 1})` : ''}</li>
          ))}
        </ul>
      </>}
    </div>
  )
}
