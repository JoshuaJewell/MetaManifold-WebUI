// © 2026 Joshua Benjamin Jewell. All rights reserved.
// Licensed under the GNU Affero General Public License version 3 (AGPLv3).
import { useEffect, useState } from 'react'
import type { TreeFile } from '../api/types'
import { indexTree, mapSupport, matchOps, parseTreeFile, type TreeOp } from './model'
import styles from './Tree.module.css'
import { readDoc } from './viewDoc'

export function ImportPanel({ load, file, ix, candidates, onClose, onApply }: {
  load: (file: string) => Promise<TreeFile>
  file: string
  ix: ReturnType<typeof indexTree>
  candidates: string[]
  onClose: () => void
  onApply: (from: string, ops: TreeOp[], settings: unknown | null) => void
}) {
  const [from, setFrom] = useState(candidates[0] ?? '')
  const [take, setTake] = useState({ edits: true, support: true, settings: false })
  const [preview, setPreview] = useState<{
    from: string; matched: TreeOp[]; unmatched: TreeOp[]; settings: unknown
    support: [string, string][]; supportTotal: number
  } | null>(null)
  const [error, setError] = useState<string | null>(null)

  useEffect(() => {
    if (!from) return
    let cancelled = false
    setPreview(null)
    setError(null)
    load(from).then(t => {
      if (cancelled) return
      const doc = readDoc(t.view)
      const sup = mapSupport(ix, parseTreeFile(t.format, t.content))
      setPreview({
        from, ...matchOps(ix, doc.ops ?? []),
        settings: (t.view as { settings?: unknown } | null)?.settings ?? null,
        support: [...sup.support], supportTotal: sup.total,
      })
    }).catch(e => { if (!cancelled) setError(String(e.message ?? e)) })
    return () => { cancelled = true }
  }, [from, load, ix])

  const tally = (ops: TreeOp[]) => {
    const n = (k: TreeOp['op']) => ops.filter(o => o.op === k).length
    const sup = ops.find(o => o.op === 'support')
    const parts = [
      n('rename') && `${n('rename')} name${n('rename') === 1 ? '' : 's'}`,
      n('collapse') && `${n('collapse')} collapsed clade${n('collapse') === 1 ? '' : 's'}`,
      n('style') && `${n('style')} style${n('style') === 1 ? '' : 's'}`,
      n('reroot') && 'the rooting',
      sup && sup.op === 'support' && `${sup.values.length} imported support values`,
    ].filter(Boolean)
    return parts.length ? parts.join(', ') : 'none'
  }

  const ready = preview && preview.from === from
  const hasEdits = !!ready && preview.matched.length > 0
  const hasSupport = !!ready && preview.support.length > 0
  const hasSettings = !!ready && preview.settings !== null
  const chosen: TreeOp[] = !ready ? [] : [
    ...(take.edits && hasEdits ? preview.matched : []),
    ...(take.support && hasSupport ? [{ op: 'support' as const, from, values: preview.support }] : []),
  ]
  const opt = (k: keyof typeof take, enabled: boolean, label: React.ReactNode) => (
    <label style={{ opacity: enabled ? 1 : 0.5 }}>
      <input type="checkbox" disabled={!enabled} checked={enabled && take[k]}
        onChange={e => setTake(t => ({ ...t, [k]: e.target.checked }))} /> {label}
    </label>
  )

  return (
    <div className={`card ${styles.toolbar}`}>
      <div className={styles.row}>
        <span className={styles.label}>Import into {file} from</span>
        {candidates.length === 0
          ? <span className={styles.muted}>No other trees to import from.</span>
          : <select value={from} onChange={e => setFrom(e.target.value)}>
              {candidates.map(c => <option key={c} value={c}>{c}</option>)}
            </select>}
        <span className={styles.spacer} />
        <button className="btn btn-primary"
          disabled={!ready || (!chosen.length && !(take.settings && hasSettings))}
          onClick={() => ready && onApply(from, chosen, take.settings && hasSettings ? preview.settings : null)}>Import</button>
        <button className="btn" onClick={onClose}>Cancel</button>
      </div>
      {error && <p className="error-msg">{error}</p>}
      {!ready && !error && from && <p className={styles.muted}>Reading {from}…</p>}
      {ready && (
        <div className={styles.row}>
          {opt('edits', hasEdits, <>Edits: {tally(preview.matched)}</>)}
          {opt('support', hasSupport, <>Support values: {preview.support.length} of {preview.supportTotal} branches match</>)}
          {opt('settings', hasSettings, <>Display settings</>)}
        </div>
      )}
      {ready && preview.unmatched.length > 0 && (
        <p className={styles.muted}>Left out, no matching branch or tip: {tally(preview.unmatched)}.</p>
      )}
      {ready && (
        <p className={styles.muted}>
          Imported edits and support undo as one step.{take.settings && hasSettings && ' Display settings replace the current ones and are not undoable.'}
        </p>
      )}
    </div>
  )
}
