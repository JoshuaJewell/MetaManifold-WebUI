// © 2026 Joshua Benjamin Jewell. All rights reserved.
// Licensed under the GNU Affero General Public License version 3 (AGPLv3).
import { useCallback, useEffect, useState } from 'react'
import { Link, useSearchParams } from 'react-router-dom'
import { api } from '../api/client'
import { errorMessage } from '../api/errorMessage'
import { REFERENCE_STEPS, type PhyloSettings, type PhyloStep, type ReferenceTreeDetail, type ReferenceTreeSummary } from '../api/types'
import { NameDialog } from '../components/NameDialog'
import { useToast } from '../components/Toast'
import { CascadeSettings, FastaButtons, StepsCard, getAt, useBatchedSave, useWorkflow } from '../phylogeny/common'
import { QCPanel } from '../phylogeny/qc'
import { TrimPreviewButton, TrimPreviewPanel, type TrimPreview } from '../phylogeny/TrimEditor'
import styles from '../phylogeny/Phylo.module.css'

const LABELS: Partial<Record<PhyloStep, string>> = {
  align: 'Align (MAFFT)', trim: 'Trim (trimAl)', tree: 'Tree (IQ-TREE)',
}

const SETTINGS = ['phylogeny.reference.align.', 'phylogeny.reference.trim.', 'phylogeny.reference.tree.']
const SETTING_LABELS: Record<string, string> = {
  'phylogeny.reference.align.': 'Align (MAFFT)',
  'phylogeny.reference.trim.':  'Trim (trimAl)',
  'phylogeny.reference.tree.':  'Tree (IQ-TREE)',
}

/** The library of reference trees that any study's placements can use. */
export function ReferenceTreesView() {
  const toast = useToast()
  const [params, setParams] = useSearchParams()
  const [list, setList] = useState<ReferenceTreeSummary[] | null>(null)
  const [dialog, setDialog] = useState<'new' | 'rename' | null>(null)
  const [qcStep, setQcStep] = useState<PhyloStep | null>(null)
  const [preview, setPreview] = useState<TrimPreview | null>(null)

  const refresh = useCallback(() => api.referenceTrees.list().then(setList).catch(() => setList([])), [])
  useEffect(() => { refresh() }, [refresh])

  const onSettled = useCallback((d: ReferenceTreeDetail) => {
    refresh()
    if (d.state === 'failed') toast.error(`${d.doc.name}: the build failed`)
    else toast.success(`${d.doc.name}: reference tree built`)
  }, [refresh, toast])
  const { detail, setDetail, open, running } = useWorkflow(api.referenceTrees.get, onSettled)
  const save = useBatchedSave(api.referenceTrees.save)

  const selected = params.get('tree')
  useEffect(() => {
    if (!list) return
    const id = selected && list.some(r => r.id === selected) ? selected : list[0]?.id ?? null
    save.flush().then(() => open(id))
    setQcStep(null)
    setPreview(null)
  }, [selected, list === null]) // eslint-disable-line react-hooks/exhaustive-deps
  const choose = (id: string) => setParams(id ? { tree: id } : {})

  const doc = detail?.doc
  const edit = (patch: Partial<Pick<NonNullable<typeof doc>, 'description' | 'overrides'>>) => {
    if (!doc) return
    setDetail(d => d && { ...d, doc: { ...d.doc, ...patch } })
    save.queue(doc.id, patch)
  }
  const setOverrides = (overrides: PhyloSettings) => edit({ overrides })

  const upload = async (text: string) => {
    if (!doc) return
    try {
      const r = await api.referenceTrees.uploadFasta(doc.id, text)
      toast.success(`${r.count} reference sequences` + (r.renamed ? `; ${r.renamed} name${r.renamed === 1 ? '' : 's'} changed for the tree tools` : ''))
      open(doc.id); refresh()
    } catch (e) { toast.error(errorMessage(e)) }
  }

  const run = async () => {
    if (!doc) return
    await save.flush()
    try { await api.referenceTrees.run(doc.id); await open(doc.id) } catch (e) { toast.error(errorMessage(e)) }
  }

  const remove = async () => {
    if (!doc || !window.confirm(`Delete the reference tree "${doc.name}" and everything built from its sequences?`)) return
    try { await api.referenceTrees.remove(doc.id); choose(''); await refresh(); open(null) } catch (e) { toast.error(errorMessage(e)) }
  }

  const loadLog = useCallback((step: PhyloStep) => doc ? api.referenceTrees.log(doc.id, step) : Promise.resolve(null), [doc?.id]) // eslint-disable-line react-hooks/exhaustive-deps
  const loadQc = useCallback((step: PhyloStep) => api.referenceTrees.qc(doc!.id, step), [doc?.id]) // eslint-disable-line react-hooks/exhaustive-deps
  const loadAlignment = useCallback((which: 'raw' | 'trimmed') => api.referenceTrees.alignment(doc!.id, which), [doc?.id]) // eslint-disable-line react-hooks/exhaustive-deps

  const trimGt = detail && (getAt(detail.settings, ['reference', 'trim', 'method']) === 'manual'
    ? getAt(detail.settings, ['reference', 'trim', 'gap_threshold']) as number : null)

  return (
    <>
      <div className="page-header">
        <h1>Reference trees</h1>
        <p>Reference alignments and trees, built once and used by placements in any study.</p>
      </div>
      <div className="card">
        <div className={styles.bar}>
          {list && list.length > 0 && (
            <select aria-label="Reference tree" value={doc?.id ?? ''} onChange={e => choose(e.target.value)}>
              {list.map(r => <option key={r.id} value={r.id}>{r.id === doc?.id ? doc.name : r.name}{r.built ? '' : ' (not built)'}</option>)}
            </select>
          )}
          <button className="btn" onClick={() => setDialog('new')}>New reference tree</button>
          {doc && <button className="btn" onClick={() => setDialog('rename')}>Rename</button>}
          {doc && <button className="btn btn-danger" disabled={running} onClick={remove}>Delete</button>}
        </div>

        {list && list.length === 0 && !doc && (
          <p className="empty-state">No reference trees yet. Upload reference sequences, then align, trim and build the tree.</p>
        )}

        {detail && doc && (
          <div className={styles.grid}>
            <div className={styles.section}>
              <div className={styles.heading}>
                Sequences
                <span className={styles.spacer} />
                <span className={styles.muted}>{detail.summary.references} references</span>
              </div>
              <FastaButtons has={detail.summary.references > 0} url={api.referenceTrees.fastaUrl(doc.id)} onFile={upload} />
              <label className={styles.muted} htmlFor="ref-desc">Description</label>
              <textarea id="ref-desc" rows={3} value={doc.description ?? ''}
                onChange={e => edit({ description: e.target.value })} />
              {detail.summary.built && (
                <div className={styles.row}>
                  <Link className="btn btn-sm" to={`/reference-trees/${doc.id}/view`}>View</Link>
                  <a href={api.referenceTrees.treefileUrl(doc.id)} download={`${doc.name}.treefile`}>Download the tree</a>
                </div>
              )}
              {detail.used_by.length > 0 && (
                <div>
                  <div className={styles.muted}>Used by</div>
                  {detail.used_by.map(u => <div key={`${u.study}/${u.id}`}><Link to={`/${u.study}?view=trees`}>{u.study}</Link>: {u.name}</div>)}
                </div>
              )}
            </div>

            <div className={styles.section}>
              <div className={styles.heading}>Settings</div>
              <div className={styles.muted}>DEFAULT values come from Default Config.</div>
              <CascadeSettings prefixes={SETTINGS} labels={SETTING_LABELS} overrides={doc.overrides}
                inherited={detail.inherited} level="tree" onChange={setOverrides} />
              <TrimPreviewButton workflow="reference" overrides={doc.overrides} inherited={detail.inherited}
                preview={t => api.referenceTrees.trimPreview(doc.id, t)} aligned={detail.qc.includes('align')}
                onPreview={setPreview} />
            </div>

            <StepsCard steps={REFERENCE_STEPS} labels={LABELS} detail={detail} runLabel="Build"
              canRun={detail.summary.references < 4 ? 'Upload at least 4 reference sequences' : null}
              onRun={run} loadLog={loadLog} qcStep={qcStep} onQc={setQcStep} />

            {preview && <TrimPreviewPanel preview={preview} loadAlignment={() => loadAlignment('raw')} onClose={() => setPreview(null)} />}
            {qcStep && (
              <QCPanel step={qcStep} label={LABELS[qcStep] ?? qcStep} loadQc={loadQc} loadAlignment={loadAlignment}
                threshold={trimGt} version={detail.status?.steps?.[qcStep]?.finished ?? ''} />
            )}
          </div>
        )}
      </div>

      {dialog && (
        <NameDialog title={dialog === 'new' ? 'New reference tree' : 'Rename reference tree'} placeholder="Name, e.g. Parabasalia"
          initialValue={dialog === 'rename' ? doc?.name ?? '' : ''}
          onClose={() => setDialog(null)}
          onConfirm={async name => {
            try {
              if (dialog === 'rename' && doc) {
                await api.referenceTrees.save(doc.id, { name })
                setDetail(d => d && { ...d, doc: { ...d.doc, name } })
              } else {
                const created = await api.referenceTrees.create(name)
                await refresh()
                choose(created.id)
                open(created.id)
              }
              await refresh()
            } catch (e) { toast.error(errorMessage(e)) }
            setDialog(null)
          }} />
      )}
    </>
  )
}
