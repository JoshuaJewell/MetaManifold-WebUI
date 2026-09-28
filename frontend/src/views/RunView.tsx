import { useState, useCallback, useEffect, useMemo } from 'react'
import { useParams, useNavigate } from 'react-router-dom'
import { useApi } from '../hooks/useApi'
import { useJobRefetch, useSSEConnected } from '../hooks/useJobEvents'
import { api } from '../api/client'
import { errorMessage } from '../api/errorMessage'
import { PipelineStages } from '../components/PipelineStages'
import { Skeleton } from '../components/Skeleton'
import { useToast } from '../components/Toast'
import { NameDialog } from '../components/NameDialog'
import { AnalysisWorkspace } from '../components/AnalysisWorkspace'
import { useTabParam } from '../hooks/useTabParam'
import { CompositionPanel } from '../components/CompositionPanel'
import { TaxaCompositionChart } from '../components/TaxaCompositionChart'
import type { ColFilter, AnnotationSource } from '../api/types'
import { carriesSource } from '../api/types'
import { ConfigAccordion } from '../components/ConfigAccordion'
import { QCPanel } from './run/QCPanel'
import { DADA2Panel } from './run/DADA2Panel'
import { TablesPanel } from './run/TablesPanel'
import { RunAlpha } from './run/RunAlpha'
import { ReadFunnel } from './run/ReadFunnel'
import { useNavRefresh } from '../hooks/useNavRefresh'

const TABS = ['pipeline', 'qc', 'dada2', 'tables', 'composition', 'analysis-dada2', 'analysis-vsearch'] as const
type Tab = (typeof TABS)[number]

// Maps tabs to the backend stages whose staleness should show a yellow dot
const TAB_STALE_STAGES: Partial<Record<Tab, string[]>> = {
  qc:    ['fastqc'],
  dada2: ['dada2_denoise', 'dada2_classify'],
}

export function RunView({ runName }: { runName?: string } = {}) {
  const { study, run: runParam, group } = useParams<{ study: string; run?: string; group?: string }>()
  const run = runName ?? runParam
  const navigate = useNavigate()
  const toast    = useToast()
  const refreshNav = useNavRefresh()
  const [tab, setTab]  = useTabParam<Tab>('view', TABS, 'pipeline')
  // Analysis tabs stay mounted once opened so their computed charts survive tab switches.
  const [analysisOpened, setAnalysisOpened] = useState<Set<Tab>>(new Set())
  useEffect(() => {
    if (tab.startsWith('analysis')) setAnalysisOpened(s => s.has(tab) ? s : new Set(s).add(tab))
  }, [tab])
  const [renaming, setRenaming] = useState(false)

  const runKey = `${study}/${group ?? ''}/${run}`
  const runFetcher = useCallback(() => api.runs.get(study!, run!, group), [study, run, group])
  const { data: runData, loading, error, refetch } = useApi(runFetcher, runKey)

  const tablesFetcher = useCallback(() => api.results.runTables(study!, run!, group), [study, run, group])
  const { data: tables, refetch: refetchTables } = useApi(tablesFetcher, runKey)

  const qcFetcher = useCallback(() => api.results.qcOutputs(study!, run!, group), [study, run, group])
  const { data: qcData, refetch: refetchQc } = useApi(qcFetcher, runKey)

  const dada2Fetcher = useCallback(() => api.results.dada2Outputs(study!, run!, group), [study, run, group])
  const { data: dada2Data, refetch: refetchDada2 } = useApi(dada2Fetcher, runKey)

  const configFetcher = useCallback(() => api.config.getRun(study!, run!, group), [study, run, group])
  const { data: configMap, refetch: refetchConfig } = useApi(configFetcher, runKey)

  const refetchOutputs = useCallback(() => {
    refetch()
    refetchTables()
    refetchQc()
    refetchDada2()
  }, [refetch, refetchTables, refetchQc, refetchDada2])

  const handleConfigChanged = useCallback(() => {
    refetchConfig()
    refetch()
  }, [refetchConfig, refetch])

  const jobFilter = useMemo(() => ({ study: study!, run: run! }), [study, run])
  useJobRefetch(refetchOutputs, jobFilter)

  // Polling fallback only when SSE is disconnected and a stage is running
  const sseConnected = useSSEConnected()
  const hasRunning = runData
    ? Object.values(runData.stages ?? {}).some(s => (s as { status: string }).status === 'running')
    : false
  useEffect(() => {
    if (!hasRunning || sseConnected) return
    const id = setInterval(refetchOutputs, 3000)
    return () => clearInterval(id)
  }, [hasRunning, sseConnected, refetchOutputs])

  // Rejects after toasting, so PipelineStages can clear its pending marker.
  const handleRunStage = async (stage: string): Promise<void> => {
    try {
      await api.pipeline.runStage(study!, run!, stage, group)
      toast.success(`${stage} queued`)
    } catch (err) {
      toast.error(`Could not start ${stage}: ${errorMessage(err)}`)
      throw err
    } finally {
      refetchOutputs()
    }
  }
  // Panels that only fire a stage do not need the rejection.
  const runStageQuietly = (stage: string) => { handleRunStage(stage).catch(() => {}) }

  const [runningAll, setRunningAll] = useState(false)
  const handleRunAll = async () => {
    setRunningAll(true)
    try {
      await api.pipeline.runRun(study!, run!, group)
      toast.success('Pipeline queued')
    } catch (err) {
      toast.error(`Could not start the pipeline: ${errorMessage(err)}`)
    } finally {
      setRunningAll(false)
      refetchOutputs()
    }
  }

  const handleRename = async (newName: string) => {
    await api.runs.rename(study!, run!, newName, group)
    setRenaming(false)
    toast.success(`Renamed to '${newName}'`)
    refreshNav()
    const basePath = group ? `/${study}/${group}/${newName}` : `/${study}/${newName}`
    navigate(basePath)
  }

  const handleDelete = async () => {
    if (!window.confirm(`Delete run '${run}'? This cannot be undone.`)) return
    try {
      await api.runs.delete(study!, run!, group)
      toast.success(`Run '${run}' deleted`)
      refreshNav()
      navigate(`/${study}`)
    } catch (err) {
      toast.error(`Could not delete run: ${errorMessage(err)}`)
    }
  }

  const isTabStale = (t: Tab): boolean => {
    const stages = TAB_STALE_STAGES[t]
    if (!stages || !runData) return false
    return stages.some(s => {
      const info = runData.stages[s as keyof typeof runData.stages]
      return info?.status === 'stale'
    })
  }

  const tablesCacheKey = useMemo(() => {
    if (!runData) return null
    const stages = [
      runData.stages.dada2_denoise,
      runData.stages.dada2_classify,
      runData.stages.swarm,
      runData.stages.vsearch,
      runData.stages.merge_taxa,
    ]
    return stages.map(s => `${s.status}:${s.last_run ?? ''}`).join('|')
  }, [runData])

  // Stable run specs for the analysis workspace (one per sub-group of a pooled
  // run, else the run itself); a fresh array each render would re-fire its
  // table-discovery effect needlessly.
  const comparisonRuns = useMemo(
    () => runData?.pooled && runData.subgroups.length > 0
      ? runData.subgroups.map(sg => ({ run: run!, group: group ?? null, prefix: sg }))
      : [{ run: run!, group: group ?? null }],
    [runData?.pooled, runData?.subgroups, run, group],
  )

  // The Tables tab chooses the table and its filters; per-run analysis uses the same selection.
  const [tableSel, setTableSel] = useState<string | null>(null)
  const [tableFilters, setTableFilters] = useState<Record<string, ColFilter>>({})
  useEffect(() => {
    if (tables && tables.length > 0 && !tables.some(t => t.id === tableSel)) {
      setTableSel(tables[0].id)
      setTableFilters({})
    }
  }, [tables, tableSel])

  return (
    <>
      <div className="page-header" style={{ display: 'flex', alignItems: 'flex-start', justifyContent: 'space-between', marginBottom: 20 }}>
        <div>
          <h1>{run}</h1>
          <p>
            {study}
            {runData && ` · ${runData.sample_count} sample${runData.sample_count !== 1 ? 's' : ''}`}
            {runData?.pooled && (
              <span style={{ color: 'var(--color-muted-fg)' }}>
                {' · '}pooling {runData.subgroups.length} sub-group{runData.subgroups.length !== 1 ? 's' : ''}
              </span>
            )}
          </p>
        </div>
        <div style={{ display: 'flex', gap: 8, flexShrink: 0 }}>
          <button className="btn" onClick={() => setRenaming(true)}>Rename</button>
          <button className="btn btn-danger" onClick={handleDelete}>Delete</button>
        </div>
      </div>

      {renaming && (
        <NameDialog
          title="Rename Run"
          initialValue={run}
          placeholder="run-name"
          onConfirm={handleRename}
          onClose={() => setRenaming(false)}
        />
      )}

      <div className="tabs">
        {TABS.map(t => {
          const label = t === 'pipeline' ? 'Pipeline'
            : t === 'qc' ? 'QC'
            : t === 'dada2' ? 'DADA2'
            : t === 'tables' ? `Tables (${tables?.length ?? 0})`
            : t === 'composition' ? 'Composition'
            : t === 'analysis-dada2' ? 'DADA2 analysis'
            : 'VSEARCH analysis'
          const stale = isTabStale(t)
          const source: AnnotationSource | null = t === 'analysis-dada2' ? 'DADA2' : t === 'analysis-vsearch' ? 'VSEARCH' : null
          const missing = source != null && runData != null && !carriesSource(runData, source)
          return (
            <button key={t} className={`tab ${tab === t ? 'active' : ''}`} onClick={() => setTab(t)} disabled={missing}
              title={missing ? `This run has no results table with ${source} taxonomy yet` : undefined}>
              {label}
              {stale && <span style={{ display: 'inline-block', width: 7, height: 7, borderRadius: '50%', background: '#f59e0b', marginLeft: 6, verticalAlign: 'middle' }} title="Config changed - re-run to update outputs" />}
            </button>
          )
        })}
      </div>

      {loading && <Skeleton lines={4} />}
      {error   && <p className="error-msg">{error}</p>}

      {tab === 'pipeline' && runData && (
        <>
          <div className="card">
            <div style={{ display: 'flex', justifyContent: 'space-between', alignItems: 'center', marginBottom: 12 }}>
              <div className="card-title">Pipeline Stages</div>
              <button className="btn btn-primary" onClick={handleRunAll} disabled={runningAll}>{runningAll ? 'Starting…' : 'Run all'}</button>
            </div>
            <PipelineStages stages={runData.stages} onRun={handleRunStage} configMap={configMap} study={study!} run={run!} group={group} onConfigChanged={handleConfigChanged} />
          </div>
          {/* Settings shared by several stages, such as the remote block. */}
          {configMap && (
            <div className="card">
              <div className="card-title" style={{ marginBottom: 12 }}>Run Settings</div>
              <ConfigAccordion
                configMap={configMap}
                study={study!}
                run={run!}
                group={group}
                onConfigChanged={handleConfigChanged}
                sourceLevel="run"
                sections={['global', 'remote']}
              />
            </div>
          )}
          <ReadFunnel study={study!} run={run!} group={group} />
        </>
      )}

      {tab === 'qc' && <QCPanel qcData={qcData ?? null} stages={runData?.stages ?? null} onRunStage={runStageQuietly} configMap={configMap} cacheKey={runData?.stages?.fastqc?.last_run ?? null} />}
      {tab === 'dada2' && <DADA2Panel study={study!} run={run!} group={group} dada2Data={dada2Data ?? null} configMap={configMap ?? null} onConfigChanged={handleConfigChanged} stages={runData?.stages ?? null} onRunStage={runStageQuietly} cacheKey={runData?.stages?.dada2_denoise?.last_run ?? null} />}
      {tab === 'tables' && (
        <TablesPanel study={study!} run={run!} group={group} tables={tables ?? []}
          onTablesChanged={refetchTables} cacheKey={tablesCacheKey}
          selected={tableSel} setSelected={setTableSel}
          filters={tableFilters} setFilters={setTableFilters} />
      )}
      {tab === 'composition' && (
        <CompositionPanel
          study={study!}
          run={run!}
          group={group}
          subgroups={runData?.subgroups}
          source={(configMap?.['tagging.source']?.value as AnnotationSource | undefined) ?? 'VSEARCH'}
        />
      )}

      {runData && (['DADA2', 'VSEARCH'] as const).map(source => {
        const t: Tab = source === 'DADA2' ? 'analysis-dada2' : 'analysis-vsearch'
        if (!analysisOpened.has(t)) return null
        if (!carriesSource(runData, source)) return tab === t
          ? <p key={source} className="empty-state">This run has no results table with {source} taxonomy yet.</p>
          : null
        return (
          <div key={source} hidden={tab !== t}>
            <AnalysisWorkspace study={study!} runs={comparisonRuns} source={source} perRun={{
              diversity: tableSel
                ? <RunAlpha study={study!} run={run!} group={group} table={tableSel} filters={tableFilters} source={source} />
                : <p className="empty-state">No results tables yet.</p>,
              composition: (
                <TaxaCompositionChart study={study!} run={run!} group={group}
                  subgroups={runData.subgroups} defaultTag="rank" table={tableSel ?? undefined} />
              ),
            }} />
          </div>
        )
      })}
    </>
  )
}
