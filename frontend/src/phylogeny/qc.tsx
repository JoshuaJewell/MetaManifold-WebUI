// © 2026 Joshua Benjamin Jewell. All rights reserved.
// Licensed under the GNU Affero General Public License version 3 (AGPLv3).
import { useEffect, useMemo, useRef, useState } from 'react'
import { errorMessage } from '../api/errorMessage'
import type { AlignmentQC, PhyloStep } from '../api/types'
import styles from './Phylo.module.css'

const pct = (x: number) => `${(100 * x).toFixed(1)}%`

/** Runs of consecutive kept columns, as [start, end) pairs. */
function runsOf(kept: number[]): [number, number][] {
  const out: [number, number][] = []
  for (const c of [...kept].sort((a, b) => a - b)) {
    const last = out[out.length - 1]
    if (last && last[1] === c) last[1] = c + 1
    else out.push([c, c + 1])
  }
  return out
}

const COLORS = { all: '#495057', references: '#1c7ed6', queries: '#e8590c', kept: 'rgba(64, 192, 87, .18)' }

/** Residue occupancy along the alignment, with the columns trimming keeps shaded. */
function OccupancyChart({ qc, threshold }: { qc: AlignmentQC; threshold?: number | null }) {
  const W = 800, H = 130, pad = { l: 30, r: 6, t: 6, b: 18 }
  const L = qc.columns || 1
  const x = (c: number) => pad.l + (c / L) * (W - pad.l - pad.r)
  const y = (v: number) => pad.t + (1 - v) * (H - pad.t - pad.b)
  const line = (v: number[]) => v.map((o, i) => `${i ? 'L' : 'M'}${x(i + 0.5).toFixed(1)},${y(o).toFixed(1)}`).join('')
  const series: [keyof typeof COLORS, number[]][] = qc.occupancy_queries
    ? [['references', qc.occupancy_references ?? []], ['queries', qc.occupancy_queries]]
    : [['all', qc.occupancy]]
  return (
    <>
      <svg className={styles.chart} viewBox={`0 0 ${W} ${H}`} role="img" aria-label="Column occupancy">
        {qc.kept_columns && runsOf(qc.kept_columns).map(([a, b]) => (
          <rect key={a} x={x(a)} y={pad.t} width={Math.max(0.5, x(b) - x(a))} height={H - pad.t - pad.b} fill={COLORS.kept} />
        ))}
        {[0, 0.5, 1].map(v => (
          <g key={v}>
            <line x1={pad.l} x2={W - pad.r} y1={y(v)} y2={y(v)} stroke="var(--color-border)" strokeWidth={0.5} />
            <text x={pad.l - 4} y={y(v) + 3} textAnchor="end">{v}</text>
          </g>
        ))}
        {threshold != null && (
          <line x1={pad.l} x2={W - pad.r} y1={y(threshold)} y2={y(threshold)} stroke="#c92a2a" strokeDasharray="4 3" strokeWidth={1} />
        )}
        {series.map(([k, v]) => <path key={k} d={line(v)} fill="none" stroke={COLORS[k]} strokeWidth={1} />)}
        <text x={pad.l} y={H - 4}>1</text>
        <text x={W - pad.r} y={H - 4} textAnchor="end">{qc.columns}</text>
      </svg>
      <div className={styles.legend}>
        {series.map(([k]) => <span key={k}><span className={styles.swatch} style={{ background: COLORS[k] }} />{k === 'all' ? 'residues per column' : k}</span>)}
        {qc.kept_columns && <span><span className={styles.swatch} style={{ background: 'rgba(64, 192, 87, .45)' }} />kept by trimming</span>}
        {threshold != null && <span><span className={styles.swatch} style={{ background: '#c92a2a' }} />-gt {threshold}</span>}
      </div>
    </>
  )
}

const BASE_COLORS: Record<string, string> = { A: '#2f9e44', C: '#1971c2', G: '#f08c00', T: '#e03131', U: '#e03131' }

/** Every residue as a coloured cell; trimmed columns are faded. */
function AlignmentViewer({ fasta, kept, queries }: { fasta: string; kept?: number[]; queries: Set<string> }) {
  const rows = useMemo(() => {
    const out: { name: string; seq: string }[] = []
    for (const block of fasta.split('>').slice(1)) {
      const nl = block.indexOf('\n')
      out.push({ name: block.slice(0, nl).trim(), seq: block.slice(nl + 1).replace(/\s/g, '').toUpperCase() })
    }
    return out
  }, [fasta])
  const L = Math.max(0, ...rows.map(r => r.seq.length))
  const [cellW, setCellW] = useState(L > 4000 ? 1 : 2)
  const cellH = rows.length > 200 ? 4 : 10
  const canvas = useRef<HTMLCanvasElement>(null)
  const [hover, setHover] = useState<string | null>(null)
  const keptSet = useMemo(() => kept ? new Set(kept) : null, [kept])

  useEffect(() => {
    const c = canvas.current
    const ctx = c?.getContext('2d')
    if (!c || !ctx) return
    c.width = Math.min(32000, L * cellW)
    c.height = rows.length * cellH
    ctx.clearRect(0, 0, c.width, c.height)
    rows.forEach((r, i) => {
      for (let j = 0; j < r.seq.length; j++) {
        const ch = r.seq[j]
        if (ch === '-' || ch === '.') continue
        ctx.fillStyle = BASE_COLORS[ch] ?? '#868e96'
        ctx.fillRect(j * cellW, i * cellH, cellW, cellH - (cellH > 4 ? 1 : 0))
      }
    })
    if (keptSet) {
      ctx.fillStyle = 'rgba(255, 255, 255, .72)'
      for (let j = 0; j < L; j++) if (!keptSet.has(j)) ctx.fillRect(j * cellW, 0, cellW, c.height)
    }
  }, [rows, L, cellW, cellH, keptSet])

  const move = (e: React.MouseEvent<HTMLCanvasElement>) => {
    const rect = e.currentTarget.getBoundingClientRect()
    const j = Math.floor((e.clientX - rect.left) / cellW), i = Math.floor((e.clientY - rect.top) / cellH)
    const r = rows[i]
    setHover(r && j < L ? `${r.name} · column ${j + 1} · ${r.seq[j] ?? '-'}${keptSet && !keptSet.has(j) ? ' · trimmed' : ''}` : null)
  }

  return (
    <div>
      <div className={styles.row}>
        <label>Column width</label>
        {[1, 2, 4, 8].map(w => (
          <label key={w}><input type="radio" checked={cellW === w} onChange={() => setCellW(w)} /> {w} px</label>
        ))}
        <span className={styles.muted}>{hover ?? 'Hover for the sequence and column.'}</span>
      </div>
      <div className={styles.aligner}>
        {cellH >= 10 && (
          <div className={styles.alignNames}>
            {rows.map(r => <div key={r.name} style={{ height: cellH, lineHeight: `${cellH}px` }}
              className={queries.has(r.name) ? styles.query : undefined}>{r.name}</div>)}
          </div>
        )}
        <canvas ref={canvas} onMouseMove={move} onMouseLeave={() => setHover(null)} />
      </div>
      <div className={styles.legend}>
        {Object.entries({ A: 'A', C: 'C', G: 'G', T: 'T/U' }).map(([b, l]) =>
          <span key={b}><span className={styles.swatch} style={{ background: BASE_COLORS[b] }} />{l}</span>)}
        <span><span className={styles.swatch} style={{ background: '#868e96' }} />ambiguous</span>
        {keptSet && <span>faded: removed by trimming</span>}
      </div>
    </div>
  )
}

type SeqSort = 'kept' | 'residues' | 'name'

export function AlignmentQCView({ qc, threshold, loadAlignment }: {
  qc: AlignmentQC
  threshold?: number | null
  loadAlignment?: () => Promise<string | null>
}) {
  const [sort, setSort] = useState<SeqSort>(qc.kept_columns ? 'kept' : 'residues')
  const [all, setAll] = useState(false)
  const [fasta, setFasta] = useState<string | null | undefined>(undefined)
  const removed = useMemo(() => new Set(qc.removed ?? []), [qc.removed])
  const queries = useMemo(() => new Set(qc.per_sequence.filter(p => p.query).map(p => p.name)), [qc.per_sequence])
  const frac = (p: AlignmentQC['per_sequence'][number]) => p.kept_residues == null || !p.residues ? 1 : p.kept_residues / p.residues
  const rows = useMemo(() => [...qc.per_sequence].sort((a, b) =>
    sort === 'name' ? a.name.localeCompare(b.name) : sort === 'residues' ? a.residues - b.residues : frac(a) - frac(b)),
  [qc.per_sequence, sort])
  const shown = all ? rows : rows.slice(0, 100)
  const kept = qc.kept_columns?.length

  return (
    <div className={styles.qcBody}>
      <div className={styles.stats}>
        <span><b>{qc.sequences}</b> sequences</span>
        <span><b>{qc.columns}</b> columns</span>
        {kept != null && <span><b>{kept}</b> kept ({pct(kept / Math.max(1, qc.columns))})</span>}
        <span><b>{pct(qc.gap_fraction)}</b> gaps</span>
        {qc.removed && <span className={qc.removed.length ? styles.state_failed : undefined}>
          <b>{qc.removed.length}</b> sequence{qc.removed.length === 1 ? '' : 's'} removed by trimming</span>}
      </div>
      <OccupancyChart qc={qc} threshold={threshold} />
      <div className={styles.scroll}>
        <table className={styles.table}>
          <thead>
            <tr>
              <th onClick={() => setSort('name')}>Sequence</th>
              <th className={styles.num} onClick={() => setSort('residues')}>Residues</th>
              {qc.kept_columns && <th className={styles.num} onClick={() => setSort('kept')}>Kept</th>}
              <th className={styles.num}>Span</th>
            </tr>
          </thead>
          <tbody>
            {shown.map(p => {
              const f = frac(p)
              const warn = removed.has(p.name) || (qc.kept_columns && f < 0.5)
              return (
                <tr key={p.name} className={warn ? styles.warn : undefined}>
                  <td>{p.name}{p.query ? ' (query)' : ''}{removed.has(p.name) ? ' - removed' : ''}</td>
                  <td className={styles.num}>{p.residues}</td>
                  {qc.kept_columns && <td className={styles.num}>{p.kept_residues} ({pct(f)})</td>}
                  <td className={styles.num}>{p.span ? `${p.span[0] + 1}-${p.span[1] + 1}` : '-'}</td>
                </tr>
              )
            })}
          </tbody>
        </table>
      </div>
      {rows.length > 100 && (
        <button className={styles.linkBtn} onClick={() => setAll(!all)}>{all ? 'Show the first 100' : `Show all ${rows.length}`}</button>
      )}
      {loadAlignment && (fasta === undefined
        ? <div><button className="btn btn-sm" onClick={() => { setFasta(null); loadAlignment().then(t => setFasta(t ?? '')) }}>Show alignment</button></div>
        : fasta === null ? <p className="loading">Loading the alignment…</p>
        : <AlignmentViewer fasta={fasta} kept={qc.kept_columns} queries={queries} />)}
    </div>
  )
}

/** The QC of the align or trim step, loaded when it opens. */
export function QCPanel({ step, label, loadQc, loadAlignment, threshold, version }: {
  step: PhyloStep
  label: string
  loadQc: (step: PhyloStep) => Promise<AlignmentQC>
  loadAlignment: (which: 'raw' | 'trimmed') => Promise<string | null>
  threshold?: number | null
  version: string
}) {
  const [qc, setQc] = useState<AlignmentQC | null>(null)
  const [error, setError] = useState<string | null>(null)
  useEffect(() => {
    let cancelled = false
    setQc(null); setError(null)
    loadQc(step).then(q => { if (!cancelled) setQc(q) }).catch(e => { if (!cancelled) setError(errorMessage(e)) })
    return () => { cancelled = true }
  }, [step, loadQc, version])

  return (
    <div className={`${styles.section} ${styles.wide}`}>
      <div className={styles.heading}>QC: {label}</div>
      {error && <p className="error-msg">{error}</p>}
      {!qc && !error && <p className="loading">Loading…</p>}
      {qc && (
        <AlignmentQCView qc={qc} threshold={step === 'trim' ? threshold : null}
          loadAlignment={() => loadAlignment('raw')} />
      )}
    </div>
  )
}
