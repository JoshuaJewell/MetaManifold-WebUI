import { useCallback } from 'react'
import { useApi } from '../hooks/useApi'
import { api } from '../api/client'
import { Skeleton } from '../components/Skeleton'
import type { DatabaseEntry } from '../api/types'

export function DatabasesView() {
  const fetcher = useCallback(() => api.databases.list(), [])
  const { data: dbs, loading, error, refetch } = useApi(fetcher)

  const download = async (key: string) => {
    await api.databases.download(key)
    refetch()
  }

  return (
    <>
      <div className="page-header">
        <h1>Databases</h1>
        <p>Reference databases for taxonomy assignment.</p>
      </div>

      {loading && <Skeleton lines={2} />}
      {error   && <p className="error-msg">{error}</p>}

      {dbs && dbs.length === 0 && (
        <div className="empty-state">No databases configured.</div>
      )}

      {dbs && dbs.length > 0 && (
        <div style={{ display: 'flex', flexDirection: 'column', gap: 8 }}>
          {dbs.map((db: DatabaseEntry) => (
            <div key={db.key} className="card" style={{ display: 'grid', gridTemplateColumns: '1fr auto', alignItems: 'center', gap: 12 }}>
              <div>
                <strong>{db.label}</strong>
                <div style={{ fontSize: '.82rem', color: 'var(--color-muted-fg)' }}>
                  DADA2: {db.dada2_available ? 'available' : 'not downloaded'}
                  {' - '}
                  vsearch: {db.vsearch_available ? 'available' : 'not downloaded'}
                </div>
              </div>
              {(!db.dada2_available || !db.vsearch_available) && (
                <button className="btn btn-primary" onClick={() => download(db.key)}>Download</button>
              )}
            </div>
          ))}
        </div>
      )}
    </>
  )
}
