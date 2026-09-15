// © 2026 Joshua Benjamin Jewell. All rights reserved.
// Licensed under the GNU Affero General Public License version 3 (AGPLv3).
import { useCallback, useEffect, useState } from 'react'
import { api } from '../api/client'
import type { ReportItem } from '../api/types'
import { useApi } from '../hooks/useApi'
import { useToast } from './Toast'

const isImage = (file: string) => /\.(svg|png)$/i.test(file)

export function ReportPanel({ study }: { study: string }) {
  const toast = useToast()
  const fetcher = useCallback(() => api.report.list(study), [study])
  const { data, loading, error } = useApi(fetcher)
  const [items, setItems] = useState<ReportItem[]>([])
  useEffect(() => { if (data) setItems(data) }, [data])

  const move = async (i: number, step: number) => {
    const j = i + step
    if (j < 0 || j >= items.length) return
    const next = [...items]
    ;[next[i], next[j]] = [next[j], next[i]]
    setItems(next)
    try { setItems(await api.report.reorder(study, next.map(it => it.id))) }
    catch (e) { toast.error(`Reorder failed: ${(e as Error).message}`) }
  }

  const rename = async (it: ReportItem, title: string) => {
    if (!title.trim() || title === it.title) return
    try {
      const saved = await api.report.rename(study, it.id, title.trim())
      setItems(xs => xs.map(x => x.id === it.id ? saved : x))
    } catch (e) { toast.error(`Rename failed: ${(e as Error).message}`) }
  }

  const remove = async (it: ReportItem) => {
    if (!window.confirm(`Remove "${it.title}" from the report?`)) return
    try {
      await api.report.remove(study, it.id)
      setItems(xs => xs.filter(x => x.id !== it.id))
    } catch (e) { toast.error(`Remove failed: ${(e as Error).message}`) }
  }

  if (loading && !data) return <p className="loading">Loading…</p>
  if (error) return <p className="error-msg">{error}</p>

  return (
    <div className="card">
      <div style={{ display: 'flex', alignItems: 'center', justifyContent: 'space-between', marginBottom: 8 }}>
        <div className="card-title" style={{ marginBottom: 0 }}>Report</div>
        <a className={`btn ${items.length ? '' : 'disabled'}`} href={items.length ? api.report.exportUrl(study) : undefined}
          aria-disabled={!items.length}>Download ZIP</a>
      </div>
      {items.length === 0 && (
        <p className="empty-state">Nothing here yet. Use Add to report on a chart, the publication tables or a tree.</p>
      )}
      <ol style={{ listStyle: 'none', display: 'flex', flexDirection: 'column', gap: 12 }}>
        {items.map((it, i) => (
          <li key={it.id} style={{ border: '1px solid var(--color-border)', borderRadius: 6, padding: 10 }}>
            <div style={{ display: 'flex', alignItems: 'center', gap: 8 }}>
              <strong style={{ minWidth: 24 }}>{i + 1}.</strong>
              <input defaultValue={it.title} aria-label={`Caption for item ${i + 1}`}
                style={{ flex: 1, font: 'inherit', padding: '2px 6px' }}
                onBlur={e => rename(it, e.target.value)}
                onKeyDown={e => { if (e.key === 'Enter') (e.target as HTMLInputElement).blur() }} />
              <span style={{ color: 'var(--color-muted-fg)', fontSize: '.8rem' }}>{it.kind}</span>
              <button className="btn btn-sm" aria-label="Move up" disabled={i === 0} onClick={() => move(i, -1)}>↑</button>
              <button className="btn btn-sm" aria-label="Move down" disabled={i === items.length - 1} onClick={() => move(i, 1)}>↓</button>
              <a className="btn btn-sm" href={api.report.fileUrl(study, it.id)} download={it.file}>Download</a>
              <button className="btn btn-sm btn-danger" onClick={() => remove(it)}>Remove</button>
            </div>
            {isImage(it.file) && (
              <img src={api.report.fileUrl(study, it.id)} alt={it.title}
                style={{ display: 'block', maxWidth: '100%', maxHeight: 360, marginTop: 8, background: '#fff' }} />
            )}
          </li>
        ))}
      </ol>
    </div>
  )
}
