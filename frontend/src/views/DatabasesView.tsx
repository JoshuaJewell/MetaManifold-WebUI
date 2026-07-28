// © 2026 Joshua Benjamin Jewell. All rights reserved.
// Licensed under the GNU Affero General Public License version 3 (AGPLv3).
import { useCallback, useEffect, useMemo, useState } from 'react'
import { useApi } from '../hooks/useApi'
import { api } from '../api/client'
import { errorMessage } from '../api/errorMessage'
import { useToast } from '../components/Toast'
import { Skeleton } from '../components/Skeleton'
import { inputStyle, fieldLabelStyle, fieldHintStyle, patchAt, removeAt, nameIssues } from '../components/EditorCard'
import { DatabaseEditor, emptyDatabaseRow, nextId, toDatabaseEntries, toDatabaseRows } from '../components/DatabaseEditor'
import type { DatabaseRow } from '../components/DatabaseEditor'
import type { DatabaseDocument, DatabaseEntry, DatabaseWarning } from '../api/types'
import { useUnsavedGuard } from '../hooks/useUnsavedGuard'

// A version token the server could not read out of a URI's basename arrives as
// an empty string, which would otherwise render as a gap in the sentence.
const versionOr = (v: string | undefined) => (v === undefined || v === '' ? 'no version in its URI' : v)

// One readable sentence per warning kind. Every one of them reports a save that
// already succeeded, so each opens by saying so.
function warningMessage(w: DatabaseWarning): string {
  const usedBy = (w.used_by ?? []).join(', ')
  switch (w.kind) {
    case 'database_removed':
      return `Saved. Database "${w.database}" is no longer in the configuration, whether removed or renamed, but these projects still resolve to it: ${usedBy}.`
    case 'levels_changed':
      return `Saved. The taxonomy levels of database "${w.database}" changed, which shifts its taxonomy columns for these projects: ${usedBy}.`
    case 'release_mismatch':
      return `Saved. Database "${w.database}" draws its two formats from different releases: DADA2 ${versionOr(w.dada2_version)}, VSEARCH ${versionOr(w.vsearch_version)}. The consensus rank compares the two labels for equality, so a mixed pair scores genuine agreements as disagreements.`
  }
}

export function DatabasesView() {
  const toast = useToast()

  const listFetcher = useCallback(() => api.databases.list(), [])
  const { data: dbs, loading: listLoading, error: listError, refetch: refetchList } = useApi(listFetcher)

  const docFetcher = useCallback(() => api.databases.document(), [])
  const { data: doc, loading: docLoading, error: docError } = useApi(docFetcher)

  const [dir, setDir]   = useState('')
  const [rows, setRows] = useState<DatabaseRow[]>([])
  const [busy, setBusy] = useState(false)
  const [dirty, setDirty] = useState(false)
  useUnsavedGuard(dirty)

  const load = useCallback((d: DatabaseDocument) => {
    setDir(d.dir)
    setRows(toDatabaseRows(d.databases))
    setDirty(false)
  }, [])

  useEffect(() => { doc && load(doc) }, [doc, load])

  const changeDir = (next: string) => { setDir(next); setDirty(true) }
  const changeRow = (i: number, next: DatabaseRow) => { setRows(patchAt(rows, i, next)); setDirty(true) }
  const removeRow = (i: number) => { setRows(removeAt(rows, i)); setDirty(true) }
  const addRow    = () => { setRows([...rows, emptyDatabaseRow(nextId(rows))]); setDirty(true) }

  // Save is gated on the same array the editor uses to mark rows.
  const keyProblems = useMemo(() => {
    const issues = nameIssues(rows.map(r => r.key))
    return rows.map((r, i) =>
      issues[i] === 'blank'     ? 'A name is required.'
      : issues[i] === 'duplicate' ? 'Two databases share this name.'
      // `dir` names the cache directory in databases.yml.
      : r.key.trim() === 'dir'    ? 'The name "dir" is reserved for the cache directory.'
      : null)
  }, [rows])

  // Save is gated on the same `nameIssues` rule the levels editor marks.
  const levelProblem = useMemo(() => {
    for (const row of rows) {
      const issues = nameIssues(row.levels.map(l => l.name))
      if (issues.includes('blank')) return `A taxonomy level in database "${row.key}" has no name.`
      if (issues.includes('duplicate')) return `Two taxonomy levels in database "${row.key}" share a name.`
    }
    return null
  }, [rows])

  // A correction's `from` becomes a map key on save, so blanks and duplicates are rejected.
  const correctionProblem = useMemo(() => {
    for (const row of rows) {
      for (const correction of row.corrections) {
        const issues = nameIssues(correction.values.map(v => v.from))
        if (issues.includes('blank')) return `A correction value in database "${row.key}" has no source label.`
        if (issues.includes('duplicate')) return `Two correction values in database "${row.key}" share a source label.`
      }
    }
    return null
  }, [rows])

  const saveProblem = useMemo(
    () => keyProblems.find(p => p !== null) ?? levelProblem ?? correctionProblem,
    [keyProblems, levelProblem, correctionProblem],
  )

  const handleSave = async () => {
    if (saveProblem) return
    setBusy(true)
    try {
      const res = await api.databases.save({ dir, databases: toDatabaseEntries(rows) })
      load(res.document)
      toast.success('Databases saved')
      // The save succeeded, so warnings are shown as info.
      for (const w of res.warnings) toast.info(warningMessage(w))
      // The badges below describe the saved config, which has just changed.
      refetchList()
    } catch (err) {
      toast.error(errorMessage(err, 'Failed to save databases'))
    } finally {
      setBusy(false)
    }
  }

  const download = async (key: string) => {
    try {
      await api.databases.download(key)
      toast.info(`Download of "${key}" queued. Watch it on the Jobs page.`)
      refetchList()
    } catch (err) {
      toast.error(errorMessage(err, `Failed to start the download of "${key}"`))
    }
  }

  return (
    <>
      <div className="page-header">
        <h1>Databases</h1>
        <p>
          Reference databases for taxonomy assignment, saved to config/databases.yml. Saving lists
          any studies affected by a rename or removal.
        </p>
      </div>

      {(docLoading || listLoading) && <Skeleton lines={4} />}
      {docError && <p className="error-msg">{docError}</p>}

      {doc && (
        <>
          <label style={{ ...fieldLabelStyle, maxWidth: 480, marginBottom: 20 }}>
            Cache directory
            <input
              value={dir}
              onChange={e => changeDir(e.target.value)}
              placeholder="./databases"
              style={{ ...inputStyle, fontFamily: 'monospace', fontSize: '.8rem' }}
              disabled={busy}
            />
            <span style={fieldHintStyle}>
              Where downloads are cached, shared by every database. Absolute, or relative to the
              working directory.
            </span>
          </label>

          <div style={{ display: 'flex', alignItems: 'center', justifyContent: 'space-between', marginBottom: 10 }}>
            <h2 style={{ fontSize: '1rem', fontWeight: 700 }}>Configured databases</h2>
            <button className="btn" onClick={addRow} disabled={busy}>+ Database</button>
          </div>

          {rows.length === 0 && (
            <div className="empty-state">No databases configured.</div>
          )}

          {rows.map((row, i) => (
            <DatabaseEditor
              key={row.id}
              row={row}
              keyProblem={keyProblems[i]}
              onChange={next => changeRow(i, next)}
              onRemove={() => removeRow(i)}
              disabled={busy}
            />
          ))}

          <div style={{ display: 'flex', gap: 12, alignItems: 'center', marginTop: 16 }}>
            <button className="btn btn-primary" onClick={handleSave} disabled={busy || saveProblem !== null}>
              {busy ? 'Saving…' : 'Save'}
            </button>
            {saveProblem && <span className="error-msg" style={{ margin: 0 }}>{saveProblem}</span>}
            {!saveProblem && dirty && (
              <span style={{ fontSize: '.8rem', color: 'var(--color-muted-fg)' }}>Unsaved changes.</span>
            )}
          </div>
        </>
      )}

      <section style={{ marginTop: 32 }}>
        <h2 style={{ fontSize: '1rem', fontWeight: 700, marginBottom: 4 }}>Availability</h2>
        {/* Badges come from the saved config, so they are greyed and downloads
            withheld while the editor has unsaved edits. */}
        <p style={{ fontSize: '.82rem', color: 'var(--color-muted-fg)', marginBottom: 10 }}>
          {dirty
            ? 'These reflect the saved configuration. Save to refresh them.'
            : 'Files present in the cache directory.'}
        </p>

        {listError && <p className="error-msg">{listError}</p>}

        {dbs && dbs.length === 0 && (
          <div className="empty-state">Nothing configured to download.</div>
        )}

        {dbs && dbs.length > 0 && (
          <div style={{ display: 'flex', flexDirection: 'column', gap: 8, opacity: dirty ? .5 : 1 }}>
            {dbs.map((db: DatabaseEntry) => (
              <div key={db.key} className="card" style={{ display: 'grid', gridTemplateColumns: '1fr auto', alignItems: 'center', gap: 12 }}>
                <div>
                  <strong>{db.label}</strong>
                  <div style={{ fontSize: '.82rem', color: 'var(--color-muted-fg)' }}>
                    DADA2: {db.dada2_available ? 'available' : 'not downloaded'}
                    {' · '}
                    VSEARCH: {db.vsearch_available ? 'available' : 'not downloaded'}
                  </div>
                </div>
                {(!db.dada2_available || !db.vsearch_available) && (
                  <button
                    className="btn btn-primary"
                    onClick={() => download(db.key)}
                    disabled={dirty || busy}
                    title={dirty ? 'Save your changes first: a download reads the saved configuration.' : undefined}
                  >
                    Download
                  </button>
                )}
              </div>
            ))}
          </div>
        )}
      </section>
    </>
  )
}
