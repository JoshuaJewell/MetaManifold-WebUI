// SPDX-FileCopyrightText: 2026 Joshua Benjamin Jewell
// SPDX-License-Identifier: MPL-2.0

import { useCallback, useEffect, useRef, useState } from 'react'
import { api } from '../api/client'
import { errorMessage } from '../api/errorMessage'
import type { AnnotationSource, ComparisonRunSpec } from '../api/types'
import { NameDialog } from '../components/NameDialog'
import { useToast } from '../components/Toast'
import { FigureEditor } from './FigureEditor'
import { emptyHistory, record, redo, undo, undoKey, type History } from './history'
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
  // Edits to the open figure, for undo and redo; cleared when another figure opens.
  const [history, setHistory] = useState<History<FigureDoc>>(emptyHistory)

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
    setHistory(emptyHistory())
    if (!id) { setDoc(null); return }
    try {
      const d = await api.figures.get(study, id)
      if (req === openRequest.current) setDoc(d)
    } catch (e) { toast.error(`Could not open the figure: ${errorMessage(e)}`) }
  }
  useEffect(() => {
    if (list && list.length && !doc && openRequest.current === 0) open(list[0].id)
  }, [list]) // eslint-disable-line react-hooks/exhaustive-deps

  const show = (next: FigureDoc) => {
    setDoc(next)
    pending.current = next
    if (timer.current) window.clearTimeout(timer.current)
    timer.current = window.setTimeout(flush, 600)
  }
  const change = (next: FigureDoc) => {
    if (doc) setHistory(h => record(h, doc, Date.now()))
    show(next)
  }
  const step = (dir: 'undo' | 'redo') => {
    if (!doc) return
    const r = (dir === 'undo' ? undo : redo)(history, doc)
    if (!r) return
    setHistory(r.history)
    show(r.doc)
  }

  // Ctrl+Z undoes and Ctrl+Shift+Z or Ctrl+Y redoes, anywhere on the page, figure edits
  // included in text boxes; not while a dialog is open.
  const stepRef = useRef(step)
  stepRef.current = step
  useEffect(() => {
    const onKey = (e: KeyboardEvent) => {
      const dir = undoKey(e)
      if (!dir || dialog) return
      e.preventDefault()
      stepRef.current(dir)
    }
    window.addEventListener('keydown', onKey)
    return () => window.removeEventListener('keydown', onKey)
  }, [dialog])

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
        {doc && <button className="btn" disabled={!history.past.length} title="Undo (Ctrl+Z)" onClick={() => step('undo')}>Undo</button>}
        {doc && <button className="btn" disabled={!history.future.length} title="Redo (Ctrl+Shift+Z or Ctrl+Y)" onClick={() => step('redo')}>Redo</button>}
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
