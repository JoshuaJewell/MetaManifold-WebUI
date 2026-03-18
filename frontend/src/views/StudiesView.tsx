import { Link, useNavigate } from 'react-router-dom'
import { useState } from 'react'
import { useApi } from '../hooks/useApi'
import { api } from '../api/client'
import { Skeleton } from '../components/Skeleton'
import { NameDialog } from '../components/NameDialog'
import { useToast } from '../components/Toast'

type Dialog = { mode: 'create' } | { mode: 'rename'; name: string }

export function StudiesView() {
  const { data: studies, loading, error, refetch } = useApi(api.studies.list)
  const [dialog, setDialog] = useState<Dialog | null>(null)
  const toast = useToast()
  const navigate = useNavigate()

  const handleCreate = async (name: string) => {
    await api.studies.create(name)
    setDialog(null)
    toast.success(`Study '${name}' created`)
    refetch()
  }

  const handleRename = async (oldName: string, newName: string) => {
    await api.studies.rename(oldName, newName)
    setDialog(null)
    toast.success(`Renamed to '${newName}'`)
    navigate(`/${newName}`)
  }

  const handleDelete = async (name: string) => {
    if (!window.confirm(`Delete study '${name}'? This cannot be undone.`)) return
    await api.studies.delete(name)
    toast.success(`Study '${name}' deleted`)
    refetch()
  }

  return (
    <>
      <div className="page-header" style={{ display: 'flex', alignItems: 'center', justifyContent: 'space-between', marginBottom: 20 }}>
        <div>
          <h1>Studies</h1>
          <p>Select a study to view runs, pipeline status, and results.</p>
        </div>
        <button className="btn btn-primary" onClick={() => setDialog({ mode: 'create' })}>
          New Study
        </button>
      </div>

      {loading && <Skeleton lines={3} />}
      {error   && <p className="error-msg">{error}</p>}

      {studies && studies.length === 0 && (
        <div className="empty-state">
          <p>No studies found.</p>
          <p>Create a study or add FASTQ data to <code>data/</code> and restart the server.</p>
        </div>
      )}

      {studies && studies.length > 0 && (
        <div className="card-grid">
          {studies.map(s => (
            <div key={s.name} className="study-card" style={{ display: 'flex', flexDirection: 'column' }}>
              <Link to={`/${s.name}`} style={{ textDecoration: 'none', color: 'inherit', flex: 1 }}>
                <h3>{s.name}</h3>
                <div className="meta">
                  {s.run_count} run{s.run_count !== 1 ? 's' : ''}
                  {' - '}
                  {s.group_count} group{s.group_count !== 1 ? 's' : ''}
                  {s.active_job_count > 0 && ` - ${s.active_job_count} active`}
                </div>
              </Link>
              <CardActions
                onRename={() => setDialog({ mode: 'rename', name: s.name })}
                onDelete={() => handleDelete(s.name)}
              />
            </div>
          ))}
        </div>
      )}

      {dialog?.mode === 'create' && (
        <NameDialog
          title="New Study"
          placeholder="study-name"
          onConfirm={handleCreate}
          onClose={() => setDialog(null)}
        />
      )}
      {dialog?.mode === 'rename' && (
        <NameDialog
          title="Rename Study"
          initialValue={dialog.name}
          placeholder="study-name"
          onConfirm={name => handleRename(dialog.name, name)}
          onClose={() => setDialog(null)}
        />
      )}
    </>
  )
}

function CardActions({ onRename, onDelete }: { onRename: () => void; onDelete: () => void }) {
  return (
    <div style={{
      display: 'flex', gap: 6, justifyContent: 'flex-end',
      marginTop: 10, paddingTop: 8,
      borderTop: '1px solid var(--color-border-light)',
    }}>
      <button
        className="btn"
        style={{ padding: '2px 8px', fontSize: '.78rem' }}
        onClick={e => { e.preventDefault(); e.stopPropagation(); onRename() }}
      >
        Rename
      </button>
      <button
        className="btn"
        style={{ padding: '2px 8px', fontSize: '.78rem', color: '#c92a2a', borderColor: '#ffc9c9' }}
        onClick={e => { e.preventDefault(); e.stopPropagation(); onDelete() }}
      >
        Delete
      </button>
    </div>
  )
}
