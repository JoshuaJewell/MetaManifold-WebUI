import { Link, useParams, useNavigate } from 'react-router-dom'
import { useCallback, useEffect, useMemo, useState } from 'react'
import { useApi } from '../hooks/useApi'
import { useJobRefetch } from '../hooks/useJobEvents'
import { api } from '../api/client'
import { Skeleton } from '../components/Skeleton'
import { NameDialog } from '../components/NameDialog'
import { CardActions, RunCard } from '../components/CardActions'
import { useToast } from '../components/Toast'
import { AnalysisWorkspace } from '../components/AnalysisWorkspace'
import { useTabParam } from '../hooks/useTabParam'
import type { AnnotationSource, ComparisonRunSpec, Run } from '../api/types'
import { carriesSource, expandRunSpecs } from '../api/types'
import { ConfigAccordion } from '../components/ConfigAccordion'
import { TreesPanel } from '../components/TreesPanel'
import { PlacementPanel } from '../components/PlacementPanel'
import { ReportPanel } from '../components/ReportPanel'
import { useNavRefresh } from '../hooks/useNavRefresh'

type GroupedRun = Run & { group?: string | null }

const VIEWS = ['runs', 'analysis-dada2', 'analysis-vsearch', 'trees', 'report', 'config'] as const
type StudyTab = (typeof VIEWS)[number]

type Dialog =
  | { mode: 'rename-study' }
  | { mode: 'new-group' }
  | { mode: 'rename-group'; name: string }
  | { mode: 'new-run' }
  | { mode: 'rename-run'; name: string }

export function StudyView() {
  const { study } = useParams<{ study: string }>()
  const navigate  = useNavigate()
  const toast     = useToast()
  const refreshNav = useNavRefresh()
  const [dialog, setDialog] = useState<Dialog | null>(null)
  const [view, setView] = useTabParam<StudyTab>('view', VIEWS, 'runs')

  const fetcher = useCallback(() => api.runs.list(study!), [study])
  const { data: runs, loading, error, refetch } = useApi(fetcher)

  const studyFetcher = useCallback(() => api.studies.get(study!), [study])
  const { data: detail, refetch: refetchDetail } = useApi(studyFetcher)

  const configFetcher = useCallback(() => api.config.getStudy(study!), [study])
  const { data: configMap, loading: configLoading, error: configError, refetch: refetchConfig } = useApi(configFetcher)

  const overridesFetcher = useCallback(() => api.config.studyOverrides(study!), [study])
  const { data: overrides } = useApi(overridesFetcher)

  const patchFn = useCallback(
    (_study: string, _run: string, body: Record<string, unknown>) =>
      api.config.patchStudy(study!, body),
    [study],
  )

  const deleteFn = useCallback(
    (_study: string, _run: string, key: string) =>
      api.config.deleteStudy(study!, key),
    [study],
  )

  // Fetch runs from each group for the analysis workspace
  const [groupRuns, setGroupRuns] = useState<GroupedRun[]>([])
  useEffect(() => {
    if (!detail?.groups?.length) { setGroupRuns([]); return }
    let cancelled = false
    Promise.all(
      detail.groups.map(async (g: string) => {
        const gRuns = await api.runs.listGroup(study!, g)
        return gRuns.map(r => ({ ...r, group: g }))
      })
    ).then(results => {
      if (!cancelled) setGroupRuns(results.flat())
    }).catch(e => {
      if (!cancelled) toast.error(`Could not list group runs: ${e.message ?? e}`)
    })
    return () => { cancelled = true }
  }, [detail?.groups, study])

  const analysisSource: AnnotationSource = view === 'analysis-dada2' ? 'DADA2' : 'VSEARCH'
  const comparisonRunsFor = useCallback((source: AnnotationSource): ComparisonRunSpec[] => [
    ...expandRunSpecs((runs ?? []).filter(r => carriesSource(r, source))),
    ...groupRuns.filter(r => carriesSource(r, source)).flatMap(r => expandRunSpecs([r], r.group)),
  ], [runs, groupRuns])
  const allComparisonRuns = useMemo(() => comparisonRunsFor(analysisSource), [comparisonRunsFor, analysisSource])
  const hasSource = (source: AnnotationSource) => runs == null || comparisonRunsFor(source).length > 0

  const placementRuns = useMemo(() => [
    ...(runs ?? []).map(r => ({ run: r.name, group: null, subgroups: r.pooled ? r.subgroups : [] })),
    ...groupRuns.map(r => ({ run: r.name, group: r.group ?? null, subgroups: r.pooled ? r.subgroups : [] })),
  ], [runs, groupRuns])
  const [treesVersion, setTreesVersion] = useState(0)

  const jobFilter = useMemo(() => ({ study: study! }), [study])
  useJobRefetch(refetch, jobFilter)

  const refetchAll = useCallback(() => { refetch(); refetchDetail() }, [refetch, refetchDetail])

  const runPipeline = async () => {
    if (!window.confirm(`Run the full pipeline for every run in '${study}'?`)) return
    try {
      await api.pipeline.runStudy(study!)
      toast.success('Pipeline queued')
    } catch (e) {
      toast.error(`Could not start the pipeline: ${(e as Error).message}`)
    }
    refetch()
  }

  const handleRenameStudy = async (newName: string) => {
    await api.studies.rename(study!, newName)
    setDialog(null)
    toast.success(`Renamed to '${newName}'`)
    refreshNav()
    navigate(`/${newName}`)
  }

  const handleDeleteStudy = async () => {
    if (!window.confirm(`Delete study '${study}'? This cannot be undone.`)) return
    try {
      await api.studies.delete(study!)
      toast.success(`Study '${study}' deleted`)
      refreshNav()
      navigate('/studies')
    } catch (e) {
      toast.error(`Delete failed: ${(e as Error).message}`)
    }
  }

  const handleNewGroup = async (name: string) => {
    await api.groups.create(study!, name)
    setDialog(null)
    toast.success(`Group '${name}' created`)
    refreshNav()
    refetchAll()
  }

  const handleRenameGroup = async (oldName: string, newName: string) => {
    await api.groups.rename(study!, oldName, newName)
    setDialog(null)
    toast.success(`Renamed to '${newName}'`)
    refreshNav()
    refetchDetail()
  }

  const handleDeleteGroup = async (name: string) => {
    if (!window.confirm(`Delete group '${name}' and all its runs? This cannot be undone.`)) return
    try {
      await api.groups.delete(study!, name)
      toast.success(`Group '${name}' deleted`)
      refreshNav()
      refetchDetail()
    } catch (e) {
      toast.error(`Delete failed: ${(e as Error).message}`)
    }
  }

  const handleNewRun = async (name: string) => {
    await api.runs.create(study!, name)
    setDialog(null)
    toast.success(`Run '${name}' created`)
    refreshNav()
    refetch()
  }

  const handleRenameRun = async (oldName: string, newName: string) => {
    await api.runs.rename(study!, oldName, newName)
    setDialog(null)
    toast.success(`Renamed to '${newName}'`)
    refreshNav()
    refetch()
  }

  const handleDeleteRun = async (name: string) => {
    if (!window.confirm(`Delete run '${name}'? This cannot be undone.`)) return
    try {
      await api.runs.delete(study!, name)
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
          <h1>{study}</h1>
          {detail && (
            <p>
              {detail.run_count} run{detail.run_count !== 1 ? 's' : ''}
              {' · '}
              {detail.group_count} group{detail.group_count !== 1 ? 's' : ''}
            </p>
          )}
        </div>
        <div style={{ display: 'flex', gap: 8, flexShrink: 0 }}>
          <button className="btn" onClick={() => setDialog({ mode: 'rename-study' })}>Rename</button>
          <button className="btn btn-danger" onClick={handleDeleteStudy}>Delete</button>
          <button className="btn btn-primary" onClick={runPipeline}>Run full pipeline</button>
        </div>
      </div>

      <div className="tabs" role="tablist">
        {VIEWS.map(v => (
          <button key={v} role="tab" aria-selected={view === v} className={`tab ${view === v ? 'active' : ''}`} onClick={() => setView(v)}
            disabled={v.startsWith('analysis') && !hasSource(v === 'analysis-dada2' ? 'DADA2' : 'VSEARCH')}
            title={v.startsWith('analysis') && !hasSource(v === 'analysis-dada2' ? 'DADA2' : 'VSEARCH')
              ? `No run has a results table with ${v === 'analysis-dada2' ? 'DADA2' : 'VSEARCH'} taxonomy yet` : undefined}>
            {v === 'runs' ? 'Runs & Groups' : v === 'analysis-dada2' ? 'DADA2 analysis' : v === 'analysis-vsearch' ? 'VSEARCH analysis' : v === 'trees' ? 'Trees' : v === 'report' ? 'Report' : 'Config'}
          </button>
        ))}
      </div>

      {loading && <Skeleton lines={3} />}
      {error   && <p className="error-msg">{error}</p>}

      {view === 'config' && configLoading && !configMap && <Skeleton lines={4} />}
      {view === 'config' && configError && <p className="error-msg">{configError}</p>}
      {view === 'config' && configMap && (
        <div style={{ marginBottom: 24 }}>
          <ConfigAccordion
            configMap={configMap}
            study={study!}
            run={study!}
            onConfigChanged={refetchConfig}
            patchFn={patchFn}
            deleteFn={deleteFn}
            sourceLevel="study"
            overrides={overrides}
          />
        </div>
      )}

      {view === 'runs' && <>
      <div style={{ display: 'flex', alignItems: 'center', justifyContent: 'space-between', marginBottom: 8 }}>
        <h2 style={{ fontSize: '1rem' }}>Groups</h2>
        <button className="btn btn-sm" onClick={() => setDialog({ mode: 'new-group' })}>+ New Group</button>
      </div>
      {detail && (!detail.groups || detail.groups.length === 0) && (
        <p style={{ color: 'var(--color-muted-fg)', fontSize: '.88rem' }}>No groups yet.</p>
      )}
      {detail && detail.groups && detail.groups.length > 0 && (
        <div className="card-grid">
          {detail.groups.map((group: string) => (
            <div key={group} className="study-card" style={{ display: 'flex', flexDirection: 'column' }}>
              <Link
                to={`/${study}/${group}`}
                style={{ textDecoration: 'none', color: 'inherit', flex: 1 }}
              >
                <h3>{group}</h3>
              </Link>
              <CardActions
                onRename={() => setDialog({ mode: 'rename-group', name: group })}
                onDelete={() => handleDeleteGroup(group)}
              />
            </div>
          ))}
        </div>
      )}

      <div style={{ display: 'flex', alignItems: 'center', justifyContent: 'space-between', marginTop: 24, marginBottom: 8 }}>
        <h2 style={{ fontSize: '1rem' }}>Runs</h2>
        <button className="btn btn-sm" onClick={() => setDialog({ mode: 'new-run' })}>+ New Run</button>
      </div>

      {runs && runs.length > 0 && (
        <div className="card-grid">
          {runs.map(run => (
            <RunCard
              key={run.name}
              name={run.name}
              to={`/${study}/${run.name}`}
              sampleCount={run.sample_count}
              stages={run.stages}
              onRename={() => setDialog({ mode: 'rename-run', name: run.name })}
              onDelete={() => handleDeleteRun(run.name)}
            />
          ))}
        </div>
      )}

      {runs && runs.length === 0 && (
        <p style={{ color: 'var(--color-muted-fg)', fontSize: '.88rem' }}>No runs yet. Add FASTQ files or create a run above.</p>
      )}

      </>}

      {view.startsWith('analysis') && (allComparisonRuns.length > 0
        ? <AnalysisWorkspace key={view} study={study!} runs={allComparisonRuns}
            source={analysisSource} />
        : <p className="empty-state">No run has a results table with {analysisSource} taxonomy yet.</p>)}

      {view === 'trees' && <>
        <PlacementPanel study={study!} runs={placementRuns} onPublished={() => setTreesVersion(v => v + 1)} />
        <TreesPanel key={treesVersion} study={study!} />
      </>}
      {view === 'report' && <ReportPanel study={study!} />}

      {dialog?.mode === 'rename-study' && (
        <NameDialog
          title="Rename Study"
          initialValue={study}
          placeholder="study-name"
          onConfirm={handleRenameStudy}
          onClose={() => setDialog(null)}
        />
      )}
      {dialog?.mode === 'new-group' && (
        <NameDialog
          title="New Group"
          placeholder="group-name"
          onConfirm={handleNewGroup}
          onClose={() => setDialog(null)}
        />
      )}
      {dialog?.mode === 'rename-group' && (
        <NameDialog
          title="Rename Group"
          initialValue={dialog.name}
          placeholder="group-name"
          onConfirm={name => handleRenameGroup(dialog.name, name)}
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

