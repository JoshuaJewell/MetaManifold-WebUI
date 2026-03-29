import { Link, useParams, useNavigate } from 'react-router-dom'
import { useCallback, useState } from 'react'
import { useApi } from '../hooks/useApi'
import { api } from '../api/client'
import { Skeleton } from '../components/Skeleton'
import { NameDialog } from '../components/NameDialog'
import { useToast } from '../components/Toast'
import { ComparisonPanel } from '../components/ComparisonPanel'
import { expandRunSpecs } from '../api/types'
import {
  STAGE_CONFIG_PREFIXES,
  STAGE_LABELS,
  StageConfig,
} from '../components/PipelineStages'

const VISIBLE_STAGES = Object.keys(STAGE_CONFIG_PREFIXES) as (keyof typeof STAGE_CONFIG_PREFIXES)[]

type Dialog =
  | { mode: 'rename-group' }
  | { mode: 'new-run' }
  | { mode: 'rename-run'; name: string }

export function GroupView({ groupName }: { groupName?: string } = {}) {
  const { study, group: groupParam } = useParams<{ study: string; group?: string }>()
  const group = groupName ?? groupParam
  const navigate = useNavigate()
  const toast    = useToast()
  const [dialog, setDialog]   = useState<Dialog | null>(null)
  const [expanded, setExpanded] = useState<string | null>(null)

  const runsFetcher = useCallback(() => api.runs.listGroup(study!, group!), [study, group])
  const { data: runs, loading, error, refetch } = useApi(runsFetcher)

  const configFetcher = useCallback(() => api.config.getGroup(study!, group!), [study, group])
  const { data: configMap, refetch: refetchConfig } = useApi(configFetcher)

  const overridesFetcher = useCallback(() => api.config.groupOverrides(study!, group!), [study, group])
  const { data: overrides } = useApi(overridesFetcher)

  const patchFn = useCallback(
    (_study: string, _group: string, body: Record<string, unknown>) =>
      api.config.patchGroup(study!, group!, body),
    [study, group],
  )

  const deleteFn = useCallback(
    (_study: string, _group: string, key: string) =>
      api.config.deleteGroup(study!, group!, key),
    [study, group],
  )

  const handleRenameGroup = async (newName: string) => {
    await api.groups.rename(study!, group!, newName)
    setDialog(null)
    toast.success(`Renamed to '${newName}'`)
    navigate(`/${study}/${newName}`)
  }

  const handleDeleteGroup = async () => {
    if (!window.confirm(`Delete group '${group}' and all its runs? This cannot be undone.`)) return
    await api.groups.delete(study!, group!)
    toast.success(`Group '${group}' deleted`)
    navigate(`/${study}`)
  }

  const handleNewRun = async (name: string) => {
    await api.runs.create(study!, name)
    setDialog(null)
    toast.success(`Run '${name}' created`)
    refetch()
  }

  const handleRenameRun = async (oldName: string, newName: string) => {
    await api.runs.rename(study!, oldName, newName)
    setDialog(null)
    toast.success(`Renamed to '${newName}'`)
    refetch()
  }

  const handleDeleteRun = async (name: string) => {
    if (!window.confirm(`Delete run '${name}'? This cannot be undone.`)) return
    await api.runs.delete(study!, name)
    toast.success(`Run '${name}' deleted`)
    refetch()
  }

  return (
    <>
      <div className="page-header" style={{ display: 'flex', alignItems: 'flex-start', justifyContent: 'space-between', marginBottom: 20 }}>
        <div>
          <h1>{group}</h1>
          {runs && (
            <p>
              {runs.length} run{runs.length !== 1 ? 's' : ''}
              {' - '}
              <Link to={`/${study}`} style={{ color: 'var(--color-muted-fg)' }}>{study}</Link>
            </p>
          )}
        </div>
        <div style={{ display: 'flex', gap: 8, flexShrink: 0 }}>
          <button className="btn" onClick={() => setDialog({ mode: 'rename-group' })}>Rename</button>
          <button className="btn" style={{ color: '#c92a2a', borderColor: '#ffc9c9' }} onClick={handleDeleteGroup}>Delete</button>
        </div>
      </div>

      <h2 style={{ fontSize: '1rem', marginBottom: 8 }}>Group Config</h2>
      {configMap && (
        <div style={{ marginBottom: 24 }}>
          {VISIBLE_STAGES.map(stage => {
            const prefixes = STAGE_CONFIG_PREFIXES[stage]
            const hasKeys = prefixes.some(p => Object.keys(configMap).some(k => k.startsWith(p)))
            if (!hasKeys) return null
            const isExpanded = expanded === stage
            return (
              <div key={stage} style={{ marginBottom: 4 }}>
                <div
                  style={{ cursor: 'pointer', fontWeight: 600, fontSize: '.85rem', padding: '4px 0' }}
                  onClick={() => setExpanded(isExpanded ? null : stage)}
                >
                  <span style={{ fontSize: '.8rem', marginRight: 6, opacity: .65 }}>{isExpanded ? '▾' : '▸'}</span>
                  {STAGE_LABELS[stage]}
                </div>
                {isExpanded && (
                  <StageConfig
                    configMap={configMap}
                    prefixes={prefixes}
                    study={study!}
                    run={group!}
                    onConfigChanged={refetchConfig}
                    patchFn={patchFn}
                    deleteFn={deleteFn}
                    sourceLevel="group"
                    overrides={overrides}
                  />
                )}
              </div>
            )
          })}
        </div>
      )}

      <div style={{ display: 'flex', alignItems: 'center', justifyContent: 'space-between', marginBottom: 8 }}>
        <h2 style={{ fontSize: '1rem' }}>Runs</h2>
        <button className="btn" style={{ padding: '3px 10px', fontSize: '.82rem' }} onClick={() => setDialog({ mode: 'new-run' })}>+ New Run</button>
      </div>

      {loading && <Skeleton lines={3} />}
      {error && <p className="error-msg">{error}</p>}

      {runs && (
        <div className="card-grid">
          {runs.map(run => {
            const stages  = Object.values(run.stages ?? {}).filter(s => s.status !== 'disabled')
            const done    = stages.filter(s => s.status === 'complete').length
            const running = stages.filter(s => s.status === 'running').length
            return (
              <div key={run.name} className="study-card" style={{ display: 'flex', flexDirection: 'column' }}>
                <Link
                  to={`/${study}/${group}/${run.name}`}
                  style={{ textDecoration: 'none', color: 'inherit', flex: 1 }}
                >
                  <h3>{run.name}</h3>
                  <div className="meta">
                    {run.sample_count} sample{run.sample_count !== 1 ? 's' : ''}
                    {' - '}
                    {done}/{stages.length} stages
                    {running > 0 && ` - ${running} running`}
                  </div>
                </Link>
                <CardActions
                  onRename={() => setDialog({ mode: 'rename-run', name: run.name })}
                  onDelete={() => handleDeleteRun(run.name)}
                />
              </div>
            )
          })}
        </div>
      )}

      {runs && runs.length >= 2 && (
        <ComparisonPanel
          study={study!}
          runs={expandRunSpecs(runs, group)}
        />
      )}

      {dialog?.mode === 'rename-group' && (
        <NameDialog
          title="Rename Group"
          initialValue={group}
          placeholder="group-name"
          onConfirm={handleRenameGroup}
          onClose={() => setDialog(null)}
        />
      )}
      {dialog?.mode === 'new-run' && (
        <NameDialog
          title="New Run"
          placeholder="run-name"
          onConfirm={handleNewRun}
          onClose={() => setDialog(null)}
        />
      )}
      {dialog?.mode === 'rename-run' && (
        <NameDialog
          title="Rename Run"
          initialValue={dialog.name}
          placeholder="run-name"
          onConfirm={name => handleRenameRun(dialog.name, name)}
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
