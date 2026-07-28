import { Link } from 'react-router-dom'
import type { RunStages } from '../api/types'

export function RunCard({ name, to, sampleCount, stages, onRename, onDelete }: {
  name: string
  to: string
  sampleCount: number
  stages: RunStages | Record<string, { status: string }>
  onRename: () => void
  onDelete: () => void
}) {
  const stageList = Object.values(stages ?? {}).filter(s => s.status !== 'disabled')
  const done    = stageList.filter(s => s.status === 'complete').length
  const running = stageList.filter(s => s.status === 'running').length
  const stale   = stageList.filter(s => s.status === 'stale').length

  return (
    <div className="study-card" style={{ display: 'flex', flexDirection: 'column' }}>
      <Link to={to} style={{ textDecoration: 'none', color: 'inherit', flex: 1 }}>
        <h3>{name}</h3>
        <div className="meta">
          {sampleCount} sample{sampleCount !== 1 ? 's' : ''}
          {' · '}
          {done}/{stageList.length} stages
          {running > 0 && (
            <span style={{ color: 'var(--color-primary)', fontWeight: 600 }}>
              {' · '}{running} running
            </span>
          )}
          {stale > 0 && (
            <span style={{ color: 'var(--color-warn)', fontWeight: 600 }}>
              {' · '}{stale} stale
            </span>
          )}
        </div>
      </Link>
      <CardActions onRename={onRename} onDelete={onDelete} />
    </div>
  )
}

export function CardActions({ onRename, onDelete }: { onRename: () => void; onDelete: () => void }) {
  return (
    <div style={{
      display: 'flex', gap: 6, justifyContent: 'flex-end',
      marginTop: 10, paddingTop: 8,
      borderTop: '1px solid var(--color-border-light)',
    }}>
      <button
        className="btn btn-sm"
        onClick={e => { e.preventDefault(); e.stopPropagation(); onRename() }}
      >
        Rename
      </button>
      <button
        className="btn btn-sm btn-danger"
        onClick={e => { e.preventDefault(); e.stopPropagation(); onDelete() }}
      >
        Delete
      </button>
    </div>
  )
}
