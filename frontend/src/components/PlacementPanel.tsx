// © 2026 Joshua Benjamin Jewell. All rights reserved.
// Licensed under the GNU Affero General Public License version 3 (AGPLv3).
import { useCallback, useEffect, useMemo, useState } from 'react'
import { Link } from 'react-router-dom'
import { api } from '../api/client'
import { errorMessage } from '../api/errorMessage'
import {
  PLACEMENT_STEPS,
  type ComparisonRunSpec, type ConfigMap, type PhyloSettings, type PhyloStep, type PlacementDetail, type PlacementQueries,
  type PlacementRunSpec, type PlacementSummary, type ReferenceTreeSummary,
} from '../api/types'
import { NameDialog } from './NameDialog'
import { useToast } from './Toast'
import { useSharedResultsTables } from './annotationShared'
import { CascadeSettings, FastaButtons, StepsCard, getAt, useBatchedSave, useWorkflow } from '../phylogeny/common'
import { QCPanel } from '../phylogeny/qc'
import { TrimPreviewButton, TrimPreviewPanel, type TrimPreview } from '../phylogeny/TrimEditor'
import styles from '../phylogeny/Phylo.module.css'

const LABELS: Partial<Record<PhyloStep, string>> = {
  align:      'Add queries (MAFFT --addfragments)',
  trim:       'Trim (trimAl)',
  place:      'Place (RAxML EPA)',
  accumulate: 'Accumulate (gappa)',
}

const SETTINGS = ['phylogeny.placement.align.', 'phylogeny.placement.trim.', 'phylogeny.placement.place.',
  'phylogeny.placement.accumulate.']
const SETTING_LABELS: Record<string, string> = {
  'phylogeny.placement.align.':      'Add queries (MAFFT)',
  'phylogeny.placement.trim.':       'Trim (trimAl)',
  'phylogeny.placement.place.':      'Place (RAxML EPA)',
  'phylogeny.placement.accumulate.': 'Accumulate (gappa)',
}

const ID_COLUMNS = new Set(['SeqName', 'Sequence', 'sequence', 'OTU', 'ASV'])

const runKey = (r: { run: string; group?: string | null }) => r.group ? `${r.group}/${r.run}` : r.run

/** Places a study's sequences on a reference tree from the library; the results appear under Trees. */
export function PlacementPanel({ study, runs, onPublished }: {
  study: string
  runs: PlacementRunSpec[]
  onPublished: () => void
}) {
  const toast = useToast()
  const [list, setList] = useState<PlacementSummary[] | null>(null)
  const [refs, setRefs] = useState<ReferenceTreeSummary[]>([])
  const [dialog, setDialog] = useState<'new' | 'rename' | null>(null)
  const [qcStep, setQcStep] = useState<PhyloStep | null>(null)
  const [preview, setPreview] = useState<TrimPreview | null>(null)

  const refresh = useCallback(() => api.placements.list(study).then(setList).catch(() => setList([])), [study])
  useEffect(() => { refresh() }, [refresh])
  useEffect(() => { api.referenceTrees.list().then(setRefs).catch(() => setRefs([])) }, [])
  const [studyConfig, setStudyConfig] = useState<ConfigMap | null>(null)
  useEffect(() => { api.config.getStudy(study).then(setStudyConfig).catch(() => setStudyConfig(null)) }, [study])

  const load = useCallback((id: string) => api.placements.get(study, id), [study])
  const onSettled = useCallback((d: PlacementDetail) => {
    refresh()
    onPublished()
    if (d.state === 'failed') toast.error(`${d.doc.name}: the placement failed`)
    else toast.success(`${d.doc.name}: placement finished`)
  }, [refresh, onPublished, toast])
  const { detail, setDetail, open, running, opened } = useWorkflow(load, onSettled)
  const sender = useCallback((id: string, patch: Parameters<typeof api.placements.save>[2]) =>
    api.placements.save(study, id, patch), [study])
  const save = useBatchedSave(sender)

  const choose = async (id: string | null) => {
    await save.flush()
    setQcStep(null)
    setPreview(null)
    open(id)
  }
  useEffect(() => {
    if (list?.length && !detail && opened.current === 0) choose(list[0].id)
  }, [list]) // eslint-disable-line react-hooks/exhaustive-deps

  const doc = detail?.doc
  const edit = (patch: Parameters<typeof api.placements.save>[2], delay = 600) => {
    if (!doc) return
    setDetail(d => d && { ...d, doc: { ...d.doc, ...patch } })
    save.queue(doc.id, patch, delay)
  }
  const setOverrides = (overrides: PhyloSettings) => edit({ overrides })

  const upload = async (text: string) => {
    if (!doc) return
    try {
      const r = await api.placements.uploadQueries(study, doc.id, text)
      toast.success(`${r.count} query sequences` + (r.renamed ? `; ${r.renamed} name${r.renamed === 1 ? '' : 's'} changed for the tree tools` : ''))
      open(doc.id); refresh()
    } catch (e) { toast.error(errorMessage(e)) }
  }

  const run = async () => {
    if (!doc) return
    await save.flush()
    try { await api.placements.run(study, doc.id); await open(doc.id) } catch (e) { toast.error(errorMessage(e)) }
  }

  const remove = async () => {
    if (!doc || !window.confirm(`Delete the placement "${doc.name}" and its working files? Trees already in the list below are kept.`)) return
    try { await api.placements.remove(study, doc.id); await refresh(); choose(null) } catch (e) { toast.error(errorMessage(e)) }
  }

  const loadLog = useCallback((step: PhyloStep) => doc ? api.placements.log(study, doc.id, step) : Promise.resolve(null), [study, doc?.id]) // eslint-disable-line react-hooks/exhaustive-deps
  const loadQc = useCallback((step: PhyloStep) => api.placements.qc(study, doc!.id, step), [study, doc?.id]) // eslint-disable-line react-hooks/exhaustive-deps
  const loadAlignment = useCallback((which: 'raw' | 'trimmed') => api.placements.alignment(study, doc!.id, which), [study, doc?.id]) // eslint-disable-line react-hooks/exhaustive-deps

  const ref = refs.find(r => r.id === doc?.reference)
  const blocker = !doc ? null
    : !doc.reference ? 'Choose a reference tree'
    : !ref ? 'The reference tree no longer exists'
    : ref.state === 'running' ? 'The reference tree is being built'
    : !ref.built ? 'Build the reference tree first'
    : doc.queries.source === 'fasta' && detail!.summary.queries === 0 ? 'Upload the query sequences'
    : null
  const trimGt = detail && (getAt(detail.settings, ['placement', 'trim', 'method']) === 'manual'
    ? getAt(detail.settings, ['placement', 'trim', 'gap_threshold']) as number : null)

  return (
    <div className="card">
      <div className={styles.bar}>
        <div className="card-title" style={{ marginBottom: 0 }}>Phylogenetic placement</div>
        {list && list.length > 0 && (
          <select aria-label="Placement" value={doc?.id ?? ''} onChange={e => choose(e.target.value)}>
            {list.map(p => <option key={p.id} value={p.id}>{p.id === doc?.id ? doc.name : p.name}</option>)}
          </select>
        )}
        <button className="btn" onClick={() => setDialog('new')}>New placement</button>
        {doc && <button className="btn" onClick={() => setDialog('rename')}>Rename</button>}
        {doc && <button className="btn btn-danger" disabled={running} onClick={remove}>Delete</button>}
        <span className={styles.spacer} />
        <Link to="/reference-trees">Reference trees</Link>
      </div>

      {list && list.length === 0 && !doc && (
        <p className="empty-state">No placements yet. A placement puts ASVs on a tree from the Reference trees library.</p>
      )}

      {detail && doc && (
        <div className={styles.grid}>
          <div className={styles.section}>
            <div className={styles.heading}>Reference tree</div>
            <div className={styles.row}>
              <select aria-label="Reference tree" value={doc.reference ?? ''}
                onChange={e => edit({ reference: e.target.value || null }, 0)}>
                <option value="">-</option>
                {refs.map(r => <option key={r.id} value={r.id}>{r.name}{r.built ? '' : ' (not built)'}</option>)}
              </select>
              {doc.reference && <Link to={`/reference-trees?tree=${doc.reference}`}>Open</Link>}
            </div>

            <div className={styles.heading} style={{ marginTop: 6 }}>
              Queries
              <span className={styles.spacer} />
              <span className={styles.muted}>{detail.summary.queries} sequences</span>
            </div>
            <div className={styles.row}>
              <label><input type="radio" checked={doc.queries.source === 'taxon'}
                onChange={() => edit({ queries: { ...doc.queries, source: 'taxon' } }, 0)} /> ASVs by taxon</label>
              <label><input type="radio" checked={doc.queries.source === 'fasta'}
                onChange={() => edit({ queries: { ...doc.queries, source: 'fasta' } }, 0)} /> FASTA file</label>
            </div>
            {doc.queries.source === 'taxon'
              ? <TaxonQueries study={study} runs={runs} queries={doc.queries} onChange={q => edit({ queries: q })} />
              : <FastaButtons has={detail.summary.queries > 0} url={api.placements.queriesUrl(study, doc.id)} onFile={upload} />}

            {doc.published && doc.published.length > 0 && (
              <div>
                <div className={styles.muted}>In the tree list:</div>
                {doc.published.map(f => <div key={f}><Link to={`/${study}/trees/${encodeURIComponent(f)}`}>{f}</Link></div>)}
              </div>
            )}
          </div>

          <div className={styles.section}>
            <div className={styles.heading}>Settings</div>
            <CascadeSettings prefixes={SETTINGS} labels={SETTING_LABELS} overrides={doc.overrides}
              inherited={detail.inherited} sources={studyConfig} level="placement" onChange={setOverrides} />
            <TrimPreviewButton workflow="placement" overrides={doc.overrides} inherited={detail.inherited}
              preview={t => api.placements.trimPreview(study, doc.id, t)} aligned={detail.qc.includes('align')}
              onPreview={setPreview} />
          </div>

          <StepsCard steps={PLACEMENT_STEPS} labels={LABELS} detail={detail} runLabel="Run placement"
            canRun={blocker} onRun={run} loadLog={loadLog} qcStep={qcStep} onQc={setQcStep} />

          {preview && <TrimPreviewPanel preview={preview} loadAlignment={() => loadAlignment('raw')} onClose={() => setPreview(null)} />}
          {qcStep && (
            <QCPanel step={qcStep} label={LABELS[qcStep] ?? qcStep} loadQc={loadQc} loadAlignment={loadAlignment}
              threshold={trimGt} version={detail.status?.steps?.[qcStep]?.finished ?? ''} />
          )}
        </div>
      )}

      {dialog && (
        <NameDialog title={dialog === 'new' ? 'New placement' : 'Rename placement'} placeholder="Name, e.g. Parabasalia in the caecum"
          initialValue={dialog === 'rename' ? doc?.name ?? '' : ''}
          onClose={() => setDialog(null)}
          onConfirm={async name => {
            try {
              if (dialog === 'rename' && doc) {
                await api.placements.save(study, doc.id, { name })
                setDetail(d => d && { ...d, doc: { ...d.doc, name } })
              } else {
                const created = await api.placements.create(study, name, refs.length === 1 ? refs[0].id : null)
                await choose(created.id)
              }
              await refresh()
            } catch (e) { toast.error(errorMessage(e)) }
            setDialog(null)
          }} />
      )}
    </div>
  )
}

function TaxonQueries({ study, runs, queries, onChange }: {
  study: string
  runs: PlacementRunSpec[]
  queries: PlacementQueries
  onChange: (q: PlacementQueries) => void
}) {
  const chosen = useMemo<ComparisonRunSpec[]>(
    () => queries.runs.map(r => ({ run: r.run, group: r.group ?? null })), [queries.runs])
  const tables = useSharedResultsTables(study, chosen)
  const [ranks, setRanks] = useState<string[]>([])
  const [values, setValues] = useState<string[]>([])
  const [preview, setPreview] = useState<string | null>(null)
  const first = chosen[0]
  const [text, setText] = useState(queries.values.join(', '))
  useEffect(() => { setText(queries.values.join(', ')) }, [queries.values])

  useEffect(() => {
    if (!first || !queries.table) { setRanks([]); return }
    let cancelled = false
    api.results.runTable(study, first.run, queries.table, { page: 1, perPage: 1 }, first.group).then(p => {
      if (cancelled) return
      const counts = new Set(p.sample_count_columns)
      setRanks(p.columns.filter(c => !counts.has(c) && !ID_COLUMNS.has(c)))
    }).catch(() => { if (!cancelled) setRanks([]) })
    return () => { cancelled = true }
  }, [study, first?.run, first?.group, queries.table]) // eslint-disable-line react-hooks/exhaustive-deps

  useEffect(() => {
    if (!first || !queries.table || !queries.rank) { setValues([]); return }
    let cancelled = false
    api.results.distinctValues(study, first.run, queries.table, queries.rank, undefined, first.group).then(d => {
      if (!cancelled) setValues(d.type === 'text' ? d.values : [])
    }).catch(() => { if (!cancelled) setValues([]) })
    return () => { cancelled = true }
  }, [study, first?.run, first?.group, queries.table, queries.rank]) // eslint-disable-line react-hooks/exhaustive-deps

  useEffect(() => {
    if (!queries.runs.length || !queries.table || !queries.rank || !queries.values.length) { setPreview(null); return }
    let cancelled = false
    const t = window.setTimeout(() => {
      api.placements.preview(study, queries).then(p => {
        if (cancelled) return
        const parts = p.per_run.map(r => `${runKey(r)}${r.subgroups.length ? ` (${r.subgroups.join(', ')})` : ''} ${r.count}`)
        setPreview(`${p.count} ASV${p.count === 1 ? '' : 's'}` + (parts.length > 1 ? `: ${parts.join(', ')}` : ''))
      }).catch(e => { if (!cancelled) setPreview(errorMessage(e)) })
    }, 400)
    return () => { cancelled = true; window.clearTimeout(t) }
  }, [study, queries])

  const spec = (r: PlacementRunSpec) => queries.runs.find(x => runKey(x) === runKey(r))
  const toggleRun = (r: PlacementRunSpec) => onChange({
    ...queries,
    runs: spec(r) ? queries.runs.filter(x => runKey(x) !== runKey(r)) : [...queries.runs, { run: r.run, group: r.group ?? null, subgroups: [] }],
  })
  const toggleSub = (r: PlacementRunSpec, sg: string) => {
    const s = spec(r)
    if (!s) return
    const subs = s.subgroups ?? []
    const next = subs.includes(sg) ? subs.filter(x => x !== sg) : [...subs, sg]
    onChange({ ...queries, runs: queries.runs.map(x => runKey(x) === runKey(r) ? { ...x, subgroups: next } : x) })
  }

  return (
    <>
      <div className={styles.row}>
        <label>Runs</label>
        <div className={styles.runList}>
          {runs.map(r => {
            const s = spec(r)
            return (
              <div key={runKey(r)}>
                <label><input type="checkbox" checked={!!s} onChange={() => toggleRun(r)} /> {runKey(r)}</label>
                {s && (r.subgroups?.length ?? 0) > 0 && (
                  <div className={styles.subgroups}>
                    {r.subgroups!.map(sg => (
                      <label key={sg}><input type="checkbox" checked={(s.subgroups ?? []).includes(sg)}
                        onChange={() => toggleSub(r, sg)} /> {sg.replace(/_/g, ' ')}</label>
                    ))}
                    {!(s.subgroups ?? []).length && <span className={styles.muted}>all subgroups</span>}
                  </div>
                )}
              </div>
            )
          })}
        </div>
      </div>
      <div className={styles.row}>
        <label>Table</label>
        <select value={queries.table} onChange={e => onChange({ ...queries, table: e.target.value })}>
          {!tables.some(t => t.table === queries.table) && <option value={queries.table}>{queries.table || '-'}</option>}
          {tables.map(t => <option key={t.key} value={t.table}>{t.label}</option>)}
        </select>
      </div>
      <div className={styles.row}>
        <label>Rank</label>
        <select value={queries.rank} onChange={e => onChange({ ...queries, rank: e.target.value, values: [] })}>
          <option value="">-</option>
          {!ranks.includes(queries.rank) && queries.rank && <option value={queries.rank}>{queries.rank}</option>}
          {ranks.map(r => <option key={r} value={r}>{r}</option>)}
        </select>
      </div>
      <div className={styles.row}>
        <label>Taxa</label>
        <input type="text" list="placement-taxa" value={text} placeholder="Comma-separated, e.g. Parabasalia"
          onChange={e => setText(e.target.value)}
          onBlur={() => onChange({ ...queries, values: text.split(',').map(v => v.trim()).filter(Boolean) })} />
        <datalist id="placement-taxa">{values.map(v => <option key={v} value={v} />)}</datalist>
      </div>
      <div className={styles.row}>
        <label>Min. reads</label>
        <input type="number" min={0} value={queries.min_reads}
          onChange={e => onChange({ ...queries, min_reads: Math.max(0, Number(e.target.value) || 0) })} />
      </div>
      {preview && <div className={styles.muted}>{preview}. Reads are counted in the chosen subgroups; the ASVs are taken from the tables each run.</div>}
    </>
  )
}
