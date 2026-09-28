// © 2026 Joshua Benjamin Jewell. All rights reserved.
// Licensed under the GNU Affero General Public License version 3 (AGPLv3).
import { useEffect, useMemo, useState } from 'react'
import { AddToReport } from './AddToReport'
import { api } from '../api/client'
import { errorMessage } from '../api/errorMessage'
import type {
  CategorySet, ComparisonRunSpec, PublicationCellKind, PublicationHeatmap, PublicationTable,
  PublicationTableRequest, PublicationValue,
} from '../api/types'
import { useSharedResultsTables } from './annotationShared'
import { useToast } from './Toast'
import styles from './PublicationTables.module.css'

type Cell = string | number | null

export function formatCell(value: Cell, kind: PublicationCellKind): string {
  if (value === null) return '–'
  if (kind === 'label' || typeof value === 'string') return String(value)
  if (kind === 'int') return value.toLocaleString('en-GB')
  if (kind === 'pct') return value > 0 && value < 0.005 ? '<0.01' : value.toFixed(2)
  const r = Math.round(value * 100) / 100
  return r === 0 ? '0.00' : `${r > 0 ? '+' : '−'}${Math.abs(r).toFixed(2)}`
}

const escapeHtml = (s: string) =>
  s.replace(/&/g, '&amp;').replace(/</g, '&lt;').replace(/>/g, '&gt;')

// Word keeps inline styles on paste, so the copied table carries its own rules.
function tableHtml(t: PublicationTable): string {
  const heavy = 'border-top:1.5pt solid black;'
  const light = '0.75pt solid black'
  const cell = 'padding:1pt 6pt;font-family:"Times New Roman";font-size:10pt;'
  const align = (kind: PublicationCellKind) => kind === 'label' ? 'text-align:left;' : 'text-align:right;'
  const head = t.header_rows.map((row, i) =>
    `<tr>${row.map(h => `<th colspan="${h.span}" style="${cell}text-align:center;${i === 0 ? heavy : ''}` +
      `${h.label ? `border-bottom:${light};` : ''}">${escapeHtml(h.label)}</th>`).join('')}</tr>`).join('')
  const labels = `<tr>${t.columns.map(c =>
    `<th style="${cell}${align(c.kind)}border-bottom:${light};${t.header_rows.length === 0 ? heavy : ''}">` +
    `${escapeHtml(c.label)}</th>`).join('')}</tr>`
  const all = [...t.rows, ...t.footer]
  const body = all.map((row, i) => `<tr>${row.map((v, j) => {
    const kind = t.columns[j].kind
    const fill = t.fills?.[i]?.[j]
    const rules = (i === t.rows.length ? `border-top:${light};` : '') +
                  (i === all.length - 1 ? 'border-bottom:1.5pt solid black;' : '') +
                  (fill ? `background-color:${fill};` : '')
    return `<td style="${cell}${align(kind)}${rules}">${escapeHtml(formatCell(v, kind))}</td>`
  }).join('')}</tr>`).join('')
  const caption = `<p style="font-family:'Times New Roman';font-size:11pt;">${escapeHtml(t.title)}</p>`
  const notes = t.notes.map(n =>
    `<p style="font-family:'Times New Roman';font-size:9pt;font-style:italic;">${escapeHtml(n)}</p>`).join('')
  return `${caption}<table style="border-collapse:collapse;">${head}${labels}${body}</table>${notes}`
}

function tableText(t: PublicationTable): string {
  const lines = [t.title, t.columns.map(c => c.label).join('\t')]
  for (const row of [...t.rows, ...t.footer]) {
    lines.push(row.map((v, j) => formatCell(v, t.columns[j].kind)).join('\t'))
  }
  return [...lines, ...t.notes].join('\n')
}

async function copyTable(t: PublicationTable) {
  const html = tableHtml(t)
  const text = tableText(t)
  if (typeof ClipboardItem !== 'undefined' && navigator.clipboard?.write) {
    await navigator.clipboard.write([new ClipboardItem({
      'text/html':  new Blob([html], { type: 'text/html' }),
      'text/plain': new Blob([text], { type: 'text/plain' }),
    })])
  } else {
    await navigator.clipboard.writeText(text)
  }
}

function TableView({ table }: { table: PublicationTable }) {
  const all = [...table.rows]
  return (
    <div className={styles.scroll}>
      <table className={styles.table}>
        <thead>
          {table.header_rows.map((row, i) => (
            <tr key={i}>
              {row.map((h, j) => (
                <th key={j} colSpan={h.span} className={h.label ? styles.group : undefined}>{h.label}</th>
              ))}
            </tr>
          ))}
          <tr>
            {table.columns.map((c, j) => (
              <th key={j} className={c.kind === 'label' ? styles.left : styles.num}>{c.label}</th>
            ))}
          </tr>
        </thead>
        <tbody>
          {all.map((row, i) => (
            <tr key={i}>
              {row.map((v, j) => {
                const fill = table.fills?.[i]?.[j]
                return (
                  <td key={j} className={`${table.columns[j].kind === 'label' ? styles.left : styles.num}${fill ? ` ${styles.heat}` : ''}`}
                    style={fill ? { backgroundColor: fill } : undefined}>
                    {formatCell(v, table.columns[j].kind)}
                  </td>
                )
              })}
            </tr>
          ))}
        </tbody>
        <tbody className={styles.footer}>
          {table.footer.map((row, i) => (
            <tr key={i}>
              {row.map((v, j) => (
                <td key={j} className={table.columns[j].kind === 'label' ? styles.left : styles.num}>
                  {formatCell(v, table.columns[j].kind)}
                </td>
              ))}
            </tr>
          ))}
        </tbody>
      </table>
      {table.notes.map((n, i) => <p key={i} className={styles.note}>{n}</p>)}
    </div>
  )
}

const VALUES: { key: PublicationValue; label: string }[] = [
  { key: 'asvs', label: 'ASVs' }, { key: 'reads', label: 'Reads' }, { key: 'pct', label: '%' },
]

/**
 * One table for a manuscript, built to order: rows of taxa at a rank or of
 * categories from a set, columns of runs, sub-groups or samples, and the values
 * each column shows. It can be copied for Word, downloaded as .xlsx or CSV, or
 * added to the report.
 */
export function PublicationTablesPanel({ study, runs }: {
  study: string
  runs: ComparisonRunSpec[]
}) {
  const toast = useToast()
  const options = useSharedResultsTables(study, runs)
  const [sets, setSets] = useState<CategorySet[]>([])
  const [table, setTable] = useState<string | null>(null)
  const [rowsBy, setRowsBy] = useState<'rank' | 'category'>('rank')
  const [ranks, setRanks] = useState<string[]>([])
  const [rank, setRank] = useState('Genus')
  const [setName, setSetName] = useState('default')
  const [columns, setColumns] = useState<'run' | 'subgroup' | 'sample'>('run')
  const [values, setValues] = useState<PublicationValue[]>(['reads', 'pct'])
  const [difference, setDifference] = useState(false)
  const [title, setTitle] = useState('')
  const [heatmapOn, setHeatmapOn] = useState(false)
  const [wholeTable, setWholeTable] = useState(false)
  const [hideZeros, setHideZeros] = useState(false)
  const [built, setBuilt] = useState<PublicationTable | null>(null)
  const [busy, setBusy] = useState<string | null>(null)

  useEffect(() => {
    api.composition.categorySets().then(setSets).catch(() => setSets([]))
  }, [])

  // Keep the table choice valid as the shared tables are discovered.
  useEffect(() => {
    const keys = options.map(o => o.key)
    setTable(current => current && keys.includes(current) ? current : keys.includes('merged') ? 'merged' : keys[0] ?? null)
  }, [options])

  const firstRun = runs[0]
  useEffect(() => {
    if (!table || !firstRun) return
    let cancelled = false
    api.analysis.ranks(study, firstRun.run, { table, group: firstRun.group, ...(firstRun.source ? { source: firstRun.source } : {}) })
      .then(rs => {
        if (cancelled) return
        setRanks(rs)
        setRank(current => rs.includes(current) ? current : rs.includes('Genus') ? 'Genus' : rs[rs.length - 1] ?? current)
      })
      .catch(() => { if (!cancelled) setRanks([]) })
    return () => { cancelled = true }
  }, [study, firstRun, table])

  const hasSubgroups = runs.some(r => r.prefix)
  const request = useMemo<PublicationTableRequest | null>(() => {
    if (!table || values.length === 0) return null
    return { runs, table, rows: rowsBy, rank, category_set: setName, columns, values,
             difference: columns === 'subgroup' && difference, title }
  }, [runs, table, rowsBy, rank, setName, columns, values, difference, title])

  // A preview built from other settings would no longer match the downloads.
  useEffect(() => { setBuilt(null) }, [request])

  const heatmap: PublicationHeatmap = heatmapOn ? (wholeTable ? 'table' : 'column') : 'none'
  const presented = useMemo(
    () => request && { ...request, heatmap, hide_zeros: hideZeros },
    [request, heatmap, hideZeros])

  // Heatmap and zero-hiding only restyle the preview, so rebuild it in place.
  const hasTable = built !== null
  useEffect(() => {
    if (!hasTable || !presented) return
    let cancelled = false
    api.analysis.publicationTable(study, presented)
      .then(r => { if (!cancelled) setBuilt(r.table) })
      .catch(err => { if (!cancelled) toast.error(`Publication table: ${errorMessage(err)}`) })
    return () => { cancelled = true }
    // Only presentation changes trigger this; data changes clear the preview instead.
    // eslint-disable-next-line react-hooks/exhaustive-deps
  }, [heatmap, hideZeros])

  const act = async (key: string, f: () => Promise<void>) => {
    setBusy(key)
    try {
      await f()
    } catch (err) {
      toast.error(`Publication table: ${errorMessage(err)}`)
    } finally {
      setBusy(null)
    }
  }

  const toggleValue = (v: PublicationValue) =>
    setValues(vs => vs.includes(v) ? vs.filter(x => x !== v) : VALUES.map(o => o.key).filter(k => k === v || vs.includes(k)))

  if (runs.length === 0) return null

  return (
    <div style={{ marginTop: 8 }}>
      <div className={styles.controls}>
        <label>
          Table
          <select value={table ?? ''} onChange={e => setTable(e.target.value)}>
            {options.map(o => <option key={o.key} value={o.key}>{o.label}</option>)}
          </select>
        </label>
        <label>
          Rows
          <select value={rowsBy} onChange={e => setRowsBy(e.target.value as 'rank' | 'category')}>
            <option value="rank">Taxa at a rank</option>
            <option value="category">Categories</option>
          </select>
        </label>
        {rowsBy === 'rank' ? (
          <label>
            Rank
            <select value={rank} onChange={e => setRank(e.target.value)}>
              {(ranks.length ? ranks : [rank]).map(r => <option key={r} value={r}>{r}</option>)}
            </select>
          </label>
        ) : (
          <label title="The Compositions set that sorts ASVs into categories">
            Category set
            <select value={setName} onChange={e => setSetName(e.target.value)}>
              {sets.length === 0 && <option value={setName}>{setName}</option>}
              {sets.map(s => <option key={s.name} value={s.name}>{s.label || s.name}</option>)}
            </select>
          </label>
        )}
        <label>
          Columns
          <select value={columns} onChange={e => setColumns(e.target.value as 'run' | 'subgroup' | 'sample')}>
            <option value="run">Runs</option>
            <option value="subgroup" disabled={!hasSubgroups}>Sub-groups</option>
            <option value="sample">Samples</option>
          </select>
        </label>
        <span className={styles.values}>
          Values
          {VALUES.map(v => (
            <label key={v.key}>
              <input type="checkbox" checked={values.includes(v.key)} onChange={() => toggleValue(v.key)} />
              {v.label}
            </label>
          ))}
        </span>
        {columns === 'subgroup' && (
          <label title="Percentage points, second sub-group minus first, for runs with exactly two">
            <input type="checkbox" checked={difference} disabled={!values.includes('pct')}
              onChange={e => setDifference(e.target.checked)} />
            Δ between two sub-groups
          </label>
        )}
      </div>
      <div className={styles.controls}>
        <label style={{ flex: 1 }}>
          Title
          <input type="text" value={title} placeholder="Describes the table when empty"
            style={{ flex: 1, font: 'inherit', padding: '1px 4px' }} onChange={e => setTitle(e.target.value)} />
        </label>
        <label>
          <input type="checkbox" checked={hideZeros} onChange={e => setHideZeros(e.target.checked)} />
          Hide zeros
        </label>
        <label>
          <input type="checkbox" checked={heatmapOn} onChange={e => setHeatmapOn(e.target.checked)} />
          Heatmap
        </label>
        {heatmapOn && (
          <button className={`btn ${styles.small}`} aria-pressed={wholeTable}
            title={wholeTable ? 'One scale per column type across the table' : 'One scale per column'}
            onClick={() => setWholeTable(w => !w)}>
            {wholeTable ? 'Scale: whole table' : 'Scale: per column'}
          </button>
        )}
      </div>
      <div style={{ display: 'flex', gap: 8, marginBottom: 12 }}>
        <button className="btn" disabled={!presented || busy !== null}
          onClick={() => presented && act('build', async () => {
            setBuilt((await api.analysis.publicationTable(study, presented)).table)
          })}>
          {busy === 'build' ? 'Building…' : 'Build table'}
        </button>
        <button className="btn" disabled={!presented || busy !== null}
          onClick={() => presented && act('xlsx', () => api.analysis.downloadPublicationTable(study, presented, 'xlsx'))}>
          {busy === 'xlsx' ? 'Preparing…' : 'Download .xlsx'}
        </button>
        <button className="btn" disabled={!request || busy !== null}
          onClick={() => request && act('csv', () => api.analysis.downloadPublicationTable(study, request, 'csv'))}>
          {busy === 'csv' ? 'Preparing…' : 'Download CSV'}
        </button>
        {built && (
          <button className="btn" disabled={busy !== null}
            onClick={() => act('copy', async () => { await copyTable(built); toast.success('Table copied') })}>
            Copy for Word
          </button>
        )}
        {presented && (
          <AddToReport study={study} kind="table" className="btn" defaultTitle={built?.title || title || 'Table'}
            make={async () => ({ blob: await api.analysis.publicationTableBlob(study, presented), ext: '.xlsx' })} />
        )}
      </div>

      {built && (
        <div className={styles.block}>
          <div className={styles.blockHead}>
            <span className={styles.caption}>{built.title}</span>
          </div>
          <TableView table={built} />
        </div>
      )}
    </div>
  )
}
