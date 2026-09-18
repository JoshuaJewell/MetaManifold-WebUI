// © 2026 Joshua Benjamin Jewell. All rights reserved.
// Licensed under the GNU Affero General Public License version 3 (AGPLv3).
import { useCallback, useEffect, useMemo, useRef, useState } from 'react'
import { errorMessage } from '../api/errorMessage'
import type { ConfigMap, PhyloSettings, PhyloStep, WorkflowDetail } from '../api/types'
import { StageConfig } from '../components/PipelineStages'
import { useToast } from '../components/Toast'
import styles from './Phylo.module.css'

/**
 * A reference tree's or placement's detail, reloaded every 3 s while its job
 * runs. `onSettled` fires when a run ends.
 */
export function useWorkflow<T extends WorkflowDetail & { doc: { id: string } }>(
  load: (id: string) => Promise<T>,
  onSettled?: (d: T) => void,
) {
  const toast = useToast()
  const [detail, setDetail] = useState<T | null>(null)
  const request = useRef(0)
  const open = useCallback(async (id: string | null) => {
    const req = ++request.current
    if (!id) { setDetail(null); return }
    try {
      const d = await load(id)
      if (req === request.current) setDetail(d)
    } catch (e) { toast.error(errorMessage(e)) }
  }, [load, toast])

  const running = detail?.state === 'running'
  const was = useRef(false)
  useEffect(() => {
    if (!detail) return
    if (was.current && !running) onSettled?.(detail)
    was.current = running
    if (!running) return
    const id = detail.doc.id
    const t = window.setInterval(() => {
      load(id).then(d => setDetail(prev => prev?.doc.id === id ? d : prev)).catch(() => {})
    }, 3000)
    return () => window.clearInterval(t)
  }, [running, detail?.doc.id]) // eslint-disable-line react-hooks/exhaustive-deps

  return { detail, setDetail, open, running, opened: request }
}

/** Edits made within `delay` go out as one save; `flush` sends them now. */
export function useBatchedSave<P extends object>(send: (id: string, patch: P) => Promise<unknown>) {
  const toast = useToast()
  const timer = useRef<number | undefined>(undefined)
  const pending = useRef<{ id: string; patch: P } | null>(null)
  const flush = useCallback(() => {
    if (timer.current) window.clearTimeout(timer.current)
    timer.current = undefined
    const p = pending.current
    pending.current = null
    return p ? send(p.id, p.patch).then(() => {}, e => toast.error(`Not saved: ${errorMessage(e)}`)) : Promise.resolve()
  }, [send, toast])
  useEffect(() => () => { flush() }, [flush])
  const queue = (id: string, patch: P, delay = 600) => {
    if (pending.current && pending.current.id !== id) flush()
    pending.current = { id, patch: { ...pending.current?.patch, ...patch } as P }
    if (timer.current) window.clearTimeout(timer.current)
    timer.current = window.setTimeout(flush, delay)
  }
  return { queue, flush }
}

export type Path = [keyof PhyloSettings & string, string, string]

type Loose = Record<string, unknown>

export const getAt = (s: PhyloSettings | undefined, [a, b, c]: Path): unknown =>
  ((((s as Loose | undefined)?.[a] as Loose | undefined)?.[b]) as Loose | undefined)?.[c]

/** `s` with one key set, or removed when `v` is undefined; empty sections are dropped. */
export function setAt(s: PhyloSettings, [a, b, c]: Path, v: unknown): PhyloSettings {
  const top = { ...((s as Loose)[a] as Loose | undefined) }
  const sec = { ...(top[b] as Loose | undefined) }
  if (v === undefined) delete sec[c]
  else sec[c] = v
  if (Object.keys(sec).length) top[b] = sec
  else delete top[b]
  const out = { ...s } as Loose
  if (Object.keys(top).length) out[a] = top
  else delete out[a]
  return out as PhyloSettings
}

const PREFIX = 'phylogeny.'
const pathOf = (key: string) => key.slice(PREFIX.length).split('.') as Path

function flatten(s: PhyloSettings, into: (key: string, value: unknown) => void) {
  for (const [a, top] of Object.entries(s as Loose)) {
    if (!top || typeof top !== 'object') continue
    for (const [b, sec] of Object.entries(top as Loose)) {
      if (!sec || typeof sec !== 'object') continue
      for (const [c, v] of Object.entries(sec as Loose)) into(`${PREFIX}${a}.${b}.${c}`, v)
    }
  }
}

/**
 * A reference tree's or placement's settings as config rows, each showing the
 * level its value comes from. Overrides are the document's own; the rest are
 * inherited through pipeline.yml, with `sources` (the study's config) telling
 * default from study values.
 */
export function CascadeSettings({ prefixes, labels, overrides, inherited, sources, level, onChange }: {
  prefixes: string[]
  labels: Record<string, string>
  overrides: PhyloSettings
  inherited: PhyloSettings
  sources?: ConfigMap | null
  level: 'placement' | 'tree'
  onChange: (next: PhyloSettings) => void
}) {
  const configMap = useMemo(() => {
    const out: ConfigMap = {}
    flatten(inherited, (k, value) => { out[k] = { value, source: sources?.[k]?.source ?? 'default' } })
    flatten(overrides, (k, value) => { out[k] = { value, source: level } })
    return out
  }, [inherited, overrides, sources, level])
  const latest = useRef(overrides)
  latest.current = overrides
  const patch = async (_s: string, _r: string, body: Record<string, unknown>) => {
    let next = latest.current
    for (const [k, v] of Object.entries(body)) next = setAt(next, pathOf(k), v)
    onChange(next)
    return configMap
  }
  const remove = async (_s: string, _r: string, key: string) => {
    onChange(setAt(latest.current, pathOf(key), undefined))
    return configMap
  }
  return (
    <StageConfig configMap={configMap} prefixes={prefixes} labels={labels} study="" run=""
      onConfigChanged={() => {}} patchFn={patch} deleteFn={remove} sourceLevel={level} />
  )
}

export function FastaButtons({ has, url, onFile }: { has: boolean; url: string; onFile: (text: string) => void }) {
  const input = useRef<HTMLInputElement>(null)
  return (
    <div className={styles.row}>
      <button className="btn btn-sm" onClick={() => input.current?.click()}>{has ? 'Replace FASTA' : 'Upload FASTA'}</button>
      {has && <a href={url}>Download</a>}
      <input ref={input} type="file" accept=".fasta,.fa,.fas,.fna,.txt" hidden
        onChange={async e => { const f = e.target.files?.[0]; e.target.value = ''; if (f) onFile(await f.text()) }} />
    </div>
  )
}

const duration = (a?: string, b?: string) => {
  if (!a || !b) return ''
  const s = Math.max(0, Math.round((Date.parse(b) - Date.parse(a)) / 1000))
  return s < 60 ? `${s} s` : s < 3600 ? `${Math.floor(s / 60)} min ${s % 60} s` : `${(s / 3600).toFixed(1)} h`
}

/** Each step with where it runs, its state, and buttons for its log and QC. */
export function StepsCard({ steps, labels, detail, runLabel, canRun, onRun, loadLog, qcStep, onQc }: {
  steps: readonly PhyloStep[]
  labels: Partial<Record<PhyloStep, string>>
  detail: WorkflowDetail
  runLabel: string
  canRun: string | null
  onRun: () => void
  loadLog: (step: PhyloStep) => Promise<string | null>
  qcStep: PhyloStep | null
  onQc: (step: PhyloStep | null) => void
}) {
  const running = detail.state === 'running'
  const [logStep, setLogStep] = useState<PhyloStep | null>(null)
  const [log, setLog] = useState<string | null>(null)
  useEffect(() => {
    if (!logStep) { setLog(null); return }
    let cancelled = false
    const get = () => loadLog(logStep).then(t => { if (!cancelled) setLog(t) })
    get()
    const t = running ? window.setInterval(get, 3000) : undefined
    return () => { cancelled = true; if (t) window.clearInterval(t) }
  }, [logStep, running, loadLog])

  const st = detail.status?.steps ?? {}
  return (
    <div className={styles.section}>
      <div className={styles.heading}>
        Run
        <span className={styles.spacer} />
        <button className="btn btn-primary btn-sm" disabled={running || !!canRun} title={canRun ?? undefined}
          onClick={onRun}>{running ? 'Running…' : runLabel}</button>
      </div>
      <table className={styles.steps}>
        <tbody>
          {steps.map(step => {
            const s = st[step]
            const where = s?.where ?? detail.remote[step] ?? 'local'
            const state = s?.state ?? (detail.status ? 'pending' : null)
            return (
              <tr key={step}>
                <td>
                  {labels[step]}
                  <div className={styles.muted}>{where === 'local' ? 'this machine' : where}</div>
                </td>
                <td className={state ? styles[`state_${state}`] : undefined}>
                  {state === 'current' ? 'up to date' : state ?? ''}
                  {s?.state === 'done' && <div className={styles.muted}>{duration(s.started, s.finished)}</div>}
                </td>
                <td className={styles.actions}>
                  {detail.qc.includes(step) && (
                    <button className={`${styles.linkBtn} ${qcStep === step ? styles.on : ''}`}
                      onClick={() => onQc(qcStep === step ? null : step)}>QC</button>
                  )}
                  {s?.where && (
                    <button className={`${styles.linkBtn} ${logStep === step ? styles.on : ''}`}
                      onClick={() => setLogStep(logStep === step ? null : step)}>Log</button>
                  )}
                </td>
              </tr>
            )
          })}
        </tbody>
      </table>
      {canRun && !running && <div className={styles.muted}>{canRun}</div>}
      {detail.status?.state === 'failed' && detail.status.error && (
        <pre className={styles.log} style={{ color: 'var(--color-danger)' }}>{detail.status.error}</pre>
      )}
      {logStep && <pre className={styles.log}>{log ?? 'No log yet.'}</pre>}
    </div>
  )
}
