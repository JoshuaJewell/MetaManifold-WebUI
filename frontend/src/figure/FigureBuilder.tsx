// © 2026 Joshua Benjamin Jewell. All rights reserved.
// Licensed under the GNU Affero General Public License version 3 (AGPLv3).
import { useCallback, useEffect, useRef, useState } from 'react'
import { api } from '../api/client'
import { errorMessage } from '../api/errorMessage'
import type { AnnotationSource, ComparisonRunSpec } from '../api/types'
import { NameDialog } from '../components/NameDialog'
import { useToast } from '../components/Toast'
import { FigureEditor } from './FigureEditor'
import { newFigure, type FigureDoc } from './types'
import styles from './Figure.module.css'

type Listing = { id: string; title: string; modified: string }

/** Multi-panel figures for the study: pick one, lay it out, export it. */
export function FigureBuilder({ study, runs, source }: {
  study: string
  runs: ComparisonRunSpec[]
  source: AnnotationSource
}) {
  const toast = useToast()
  const [list, setList] = useState<Listing[] | null>(null)
  const [doc, setDoc] = useState<FigureDoc | null>(null)
  const [dialog, setDialog] = useState<'new' | 'rename' | null>(null)
  const pending = useRef<FigureDoc | null>(null)
  const timer = useRef<number | undefined>(undefined)

  const flush = useCallback(() => {
    if (timer.current) window.clearTimeout(timer.current)
    timer.current = undefined
    const d = pending.current
    pending.current = null
    if (d) api.figures.save(study, d).catch(e => toast.error(`Figure not saved: ${errorMessage(e)}`))
  }, [study, toast])
  useEffect(() => flush, [flush])

  const refresh = useCallback(() => api.figures.list(study).then(setList).catch(() => setList([])), [study])
  useEffect(() => { refresh() }, [refresh])

  // Only the most recent choice is shown, however the requests finish.
  const openRequest = useRef(0)
  const open = async (id: string) => {
    flush()
    const req = ++openRequest.current
    if (!id) { setDoc(null); return }
    try {
      const d = await api.figures.get(study, id)
      if (req === openRequest.current) setDoc(d)
    } catch (e) { toast.error(`Could not open the figure: ${errorMessage(e)}`) }
  }
  useEffect(() => {
    if (list && list.length && !doc && openRequest.current === 0) open(list[0].id)
  }, [list]) // eslint-disable-line react-hooks/exhaustive-deps

  const change = (next: FigureDoc) => {
    setDoc(next)
    pending.current = next
    if (timer.current) window.clearTimeout(timer.current)
    timer.current = window.setTimeout(flush, 600)
  }

  const remove = async () => {
    if (!doc || !window.confirm(`Delete the figure "${doc.title}"? This cannot be undone.`)) return
    pending.current = null
    try {
      await api.figures.remove(study, doc.id)
      setDoc(null)
      const next = await api.figures.list(study)
      setList(next)
      if (next.length) open(next[0].id)
    } catch (e) { toast.error(`Delete failed: ${errorMessage(e)}`) }
  }

  return (
    <div className="card">
      <div className={styles.bar}>
        <div className="card-title" style={{ marginBottom: 0 }}>Figures</div>
        {list && list.length > 0 && (
          <select aria-label="Figure" value={doc?.id ?? ''} onChange={e => open(e.target.value)}>
            {list.map(f => <option key={f.id} value={f.id}>{f.id === doc?.id ? doc.title : f.title}</option>)}
          </select>
        )}
        <button className="btn" onClick={() => setDialog('new')}>New figure</button>
        {doc && <button className="btn" onClick={() => setDialog('rename')}>Rename</button>}
        {doc && <button className="btn btn-danger" onClick={remove}>Delete</button>}
      </div>

      {list && list.length === 0 && !doc && (
        <p className="empty-state">No figures yet. A figure is a page of lettered groups, each a grid of charts drawn from the results.</p>
      )}
      {doc && <FigureEditor key={doc.id} study={study} runs={runs} source={source} doc={doc} onChange={change} />}

      {dialog && (
        <NameDialog title={dialog === 'new' ? 'New figure' : 'Rename figure'} placeholder="Title"
          initialValue={dialog === 'rename' ? doc?.title ?? '' : `Figure ${(list?.length ?? 0) + 1}`}
          onClose={() => setDialog(null)}
          onConfirm={async title => {
            if (dialog === 'rename' && doc) {
              change({ ...doc, title })
              flush()
              setList(l => l?.map(f => f.id === doc.id ? { ...f, title } : f) ?? l)
            } else {
              flush()
              const created = await api.figures.create(study, newFigure(title))
              setDoc(created)
              await refresh()
            }
            setDialog(null)
          }} />
      )}
    </div>
  )
}
