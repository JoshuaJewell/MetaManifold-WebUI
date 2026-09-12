// © 2026 Joshua Benjamin Jewell. All rights reserved.
// Licensed under the GNU Affero General Public License version 3 (AGPLv3).
import { useCallback, useRef, useState } from 'react'
import { Link } from 'react-router-dom'
import { api } from '../api/client'
import { useApi } from '../hooks/useApi'
import { useToast } from './Toast'
import styles from '../tree/Tree.module.css'

const ACCEPT = '.jplace,.nwk,.newick,.tre,.tree,.treefile'

const kb = (n: number) => n < 1024 * 1024 ? `${Math.max(1, Math.round(n / 1024))} KB` : `${(n / 1024 / 1024).toFixed(1)} MB`

export function TreesPanel({ study }: { study: string }) {
  const toast = useToast()
  const fetcher = useCallback(() => api.trees.list(study), [study])
  const { data: trees, loading, error, refetch } = useApi(fetcher)
  const input = useRef<HTMLInputElement>(null)
  const [busy, setBusy] = useState(false)

  const upload = async (files: FileList | null) => {
    if (!files?.length) return
    setBusy(true)
    const existing = new Set((trees ?? []).map(t => t.file))
    for (const f of Array.from(files)) {
      const name = f.name.replace(/[^A-Za-z0-9_.-]/g, '_').replace(/^\.+/, '')
      const overwrite = existing.has(name)
      if (overwrite && !window.confirm(`Replace ${name}? Its saved view is kept and replayed on the new file.`)) continue
      try {
        await api.trees.upload(study, name, await f.text(), overwrite)
        toast.success(`Added ${name}`)
      } catch (e) {
        toast.error(`${name}: ${(e as Error).message}`)
      }
    }
    setBusy(false)
    if (input.current) input.current.value = ''
    refetch()
  }

  const remove = async (file: string) => {
    if (!window.confirm(`Delete ${file} and its saved view? This cannot be undone.`)) return
    try {
      await api.trees.delete(study, file)
      toast.success(`Deleted ${file}`)
    } catch (e) {
      toast.error(`${file}: ${(e as Error).message}`)
    }
    refetch()
  }

  return (
    <div className="card">
      <div style={{ display: 'flex', alignItems: 'center', justifyContent: 'space-between', marginBottom: 8 }}>
        <div className="card-title" style={{ marginBottom: 0 }}>Trees</div>
        <button className="btn" disabled={busy} onClick={() => input.current?.click()}>
          {busy ? 'Uploading…' : 'Add tree files'}
        </button>
        <input ref={input} type="file" multiple accept={ACCEPT} hidden onChange={e => upload(e.target.files)} />
      </div>
      {loading && <p className="loading">Loading…</p>}
      {error && <p className="error-msg">{error}</p>}
      {trees && trees.length === 0 && (
        <p className="empty-state">No trees yet. Add Newick or jplace files (EPA-ng, pplacer, RAxML -f v).</p>
      )}
      {trees && trees.length > 0 && (
        <table className={styles.fileTable}>
          <thead><tr><th>File</th><th>Format</th><th>Size</th><th>Modified</th><th>View</th><th /></tr></thead>
          <tbody>
            {trees.map(t => (
              <tr key={t.file}>
                <td><Link to={`/${study}/trees/${encodeURIComponent(t.file)}`}>{t.file}</Link></td>
                <td>{t.format === 'jplace' ? 'jplace' : 'Newick'}</td>
                <td>{kb(t.size)}</td>
                <td>{t.modified.replace('T', ' ').slice(0, 16)}</td>
                <td>{t.has_view ? 'Saved' : ''}</td>
                <td style={{ textAlign: 'right' }}>
                  <button className={styles.linkBtn} style={{ color: 'var(--color-danger)' }} onClick={() => remove(t.file)}>Delete</button>
                </td>
              </tr>
            ))}
          </tbody>
        </table>
      )}
    </div>
  )
}
