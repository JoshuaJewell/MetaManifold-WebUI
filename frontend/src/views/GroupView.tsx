import { Link, useParams, useNavigate } from 'react-router-dom'
import { useCallback, useMemo, useState } from 'react'
import { useApi } from '../hooks/useApi'
import { api } from '../api/client'
import { Skeleton } from '../components/Skeleton'
import { NameDialog } from '../components/NameDialog'
import { RunCard } from '../components/CardActions'
import { useToast } from '../components/Toast'
import { AnalysisWorkspace } from '../components/AnalysisWorkspace'
import { useTabParam } from '../hooks/useTabParam'
import { expandRunSpecs } from '../api/types'
import { ALL_SECTIONS, ConfigAccordion } from '../components/ConfigAccordion'
import { useNavRefresh } from '../hooks/useNavRefresh'

// Placements read the study's settings, so phylogeny is set there.
const GROUP_SECTIONS = ALL_SECTIONS.filter(s => s !== 'phylogeny')

const VIEWS = ['runs', 'analysis-dada2', 'analysis-vsearch', 'config'] as const
type GroupTab = (typeof VIEWS)[number]

type Dialog =
  | { mode: 'rename-group' }
  | { mode: 'new-run' }
  | { mode: 'rename-run'; name: string }

export function GroupView({ groupName }: { groupName?: string } = {}) {
  const { study, group: groupParam } = useParams<{ study: string; group?: string }>()
  const group = groupName ?? groupParam
  const navigate = useNavigate()
  const toast    = useToast()
  const refreshNav = useNavRefresh()
  const [dialog, setDialog]   = useState<Dialog | null>(null)
  const [view, setView] = useTabParam<GroupTab>('view', VIEWS, 'runs')

  const runsFetcher = useCallback(() => api.runs.listGroup(study!, group!), [study, group])
  const { data: runs, loading, error, refetch } = useApi(runsFetcher)

  // Stable specs, so the workspace's table discovery does not refire each render.
  const comparisonRuns = useMemo(() => expandRunSpecs(runs ?? [], group), [runs, group])

  const configFetcher = useCallback(() => api.config.getGroup(study!, group!), [study, group])
  const { data: configMap, loading: configLoading, error: configError, refetch: refetchConfig } = useApi(configFetcher)

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
    refreshNav()
    navigate(`/${study}/${newName}`)
  }

  const handleDeleteGroup = async () => {
    if (!window.confirm(`Delete group '${group}' and all its runs? This cannot be undone.`)) return
    try {
      await api.groups.delete(study!, group!)
      toast.success(`Group '${group}' deleted`)
      refreshNav()
      navigate(`/${study}`)
    } catch (e) {
      toast.error(`Delete failed: ${(e as Error).message}`)
    }
  }

  const handleNewRun = async (name: string) => {
    await api.runs.create(study!, name, group)
    setDialog(null)
    toast.success(`Run '${name}' created`)
    refreshNav()
    refetch()
  }

  const handleRenameRun = async (oldName: string, newName: string) => {
    await api.runs.rename(study!, oldName, newName, group)
    setDialog(null)
    toast.success(`Renamed to '${newName}'`)
    refreshNav()
    refetch()
  }

  const handleDeleteRun = async (name: string) => {
    if (!window.confirm(`Delete run '${name}'? This cannot be undone.`)) return
    try {
      await api.runs.delete(study!, name, group)
      toast.success(`Run '${name}' deleted`)
      refreshNav()
      refetch()
    } catch (e) {
      toast.error(`Delete failed: ${(e as Error).message}`)
    }
  }

  return (
    <>
      <div className="page-header" style={{ display: 'flex', alignItems: 'flex-start', justifyContent: 'space-between', marginBottom: 20 }}>
        <div>
          <h1>{group}</h1>
          {runs && (
            <p>
              {runs.length} run{runs.length !== 1 ? 's' : ''}
              {' · '}
              <Link to={`/${study}`} style={{ color: 'var(--color-muted-fg)' }}>{study}</Link>
            </p>
          )}
        </div>
        <div style={{ display: 'flex', gap: 8, flexShrink: 0 }}>
          <button className="btn" onClick={() => setDialog({ mode: 'rename-group' })}>Rename</button>
          <button className="btn btn-danger" onClick={handleDeleteGroup}>Delete</button>
        </div>
      </div>

      <div className="tabs" role="tablist">
        {VIEWS.map(v => (
          <button key={v} role="tab" aria-selected={view === v} className={`tab ${view === v ? 'active' : ''}`} onClick={() => setView(v)}>
            {v === 'runs' ? 'Runs' : v === 'analysis-dada2' ? 'DADA2 analysis' : v === 'analysis-vsearch' ? 'VSEARCH analysis' : 'Config'}
          </button>
        ))}
      </div>

      {view === 'config' && configLoading && !configMap && <Skeleton lines={4} />}
      {view === 'config' && configError && <p className="error-msg">{configError}</p>}
      {view === 'config' && configMap && (
        <div style={{ marginBottom: 24 }}>
          <ConfigAccordion
            configMap={configMap}
            study={study!}
            run={group!}
            onConfigChanged={refetchConfig}
            patchFn={patchFn}
            deleteFn={deleteFn}
            sourceLevel="group"
            overrides={overrides}
            sections={GROUP_SECTIONS}
          />
        </div>
      )}

      {view === 'runs' && <>
      <div style={{ display: 'flex', alignItems: 'center', justifyContent: 'space-between', marginBottom: 8 }}>
        <h2 style={{ fontSize: '1rem' }}>Runs</h2>
        <button className="btn btn-sm" onClick={() => setDialog({ mode: 'new-run' })}>+ New Run</button>
      </div>

      {loading && <Skeleton lines={3} />}
      {error && <p className="error-msg">{error}</p>}

      {runs && runs.length === 0 && (
        <p style={{ color: 'var(--color-muted-fg)', fontSize: '.88rem' }}>No runs in this group.</p>
      )}

      {runs && runs.length > 0 && (
        <div className="card-grid">
          {runs.map(run => (
            <RunCard
              key={run.name}
              name={run.name}
              to={`/${study}/${group}/${run.name}`}
              sampleCount={run.sample_count}
              stages={run.stages}
              onRename={() => setDialog({ mode: 'rename-run', name: run.name })}
              onDelete={() => handleDeleteRun(run.name)}
            />
          ))}
        </div>
      )}

      </>}

      {view.startsWith('analysis') && (comparisonRuns.length > 0
        ? <AnalysisWorkspace key={view} study={study!} runs={comparisonRuns}
            source={view === 'analysis-dada2' ? 'DADA2' : 'VSEARCH'} />
        : <p className="empty-state">No runs to analyse yet.</p>)}

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

