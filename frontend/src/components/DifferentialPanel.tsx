// SPDX-License-Identifier: AGPL-3.0-only
// SPDX-FileCopyrightText: 2026 Jonathan D.A. Jewell (hyperpolymath) <j.d.a.jewell@open.ac.uk>
import { useEffect, useRef, useState } from 'react'
import { api } from '../api/client'
import { errorMessage } from '../api/errorMessage'
import { useToast } from './Toast'
import type { ComparisonRunSpec, DifferentialResult } from '../api/types'
import { sharedRanks } from './sharedRanks'
import type { AnalysisOption } from './annotationShared'
import { AnalysisChart } from './AnalysisChart'
import { AddToReport } from './AddToReport'
import { saveBlob } from '../utils/download'
import { differentialCsv, formatStat } from './differentialTable'

const cell = { padding: '3px 10px', textAlign: 'right' } as const
const head = { padding: '4px 10px', textAlign: 'right' } as const

/** The CSV as a blob, ready to download or file in the report. */
function csvBlob(result: DifferentialResult): Blob {
  return new Blob([differentialCsv(result)], { type: 'text/csv;charset=utf-8' })
}

/**
 * Differential abundance between exactly two runs or sub-groups: a
 * negative-binomial model per taxon at the chosen rank, with BH-adjusted
 * p-values, shown as a volcano plot and a table that can be downloaded or
 * added to the report. The first selected condition is the reference.
 */
export function DifferentialPanel({ study, runs, option, aggregate }: {
  study: string
  runs: ComparisonRunSpec[]
  option: AnalysisOption | null
  aggregate: boolean
}) {
  const toast = useToast()
  const [ranks, setRanks] = useState<string[]>([])
  const [rank, setRank] = useState<string>('')
  const [ranksReady, setRanksReady] = useState(false)
  const [ranksError, setRanksError] = useState<string | null>(null)
  const [result, setResult] = useState<DifferentialResult | null>(null)
  const [loading, setLoading] = useState(false)
  // Only the response to the latest request, for the current inputs, is shown.
  const request = useRef(0)

  // The ranks shared by every selected run.
  useEffect(() => {
    setRanksError(null)
    if (!option || runs.length === 0) { setRanks([]); setRank(''); setRanksReady(false); return }
    let cancelled = false
    setRanksReady(false)
    sharedRanks(runs, r => api.analysis.ranks(study, r.run, {
      table: option.table,
      group: r.group,
      ...(r.source ? { source: r.source } : {}),
    })).then(shared => {
      if (cancelled) return
      setRanks(shared)
      setRanksReady(true)
      setRank(current => shared.includes(current) ? current : (shared[shared.length - 1] ?? ''))
    }).catch(err => {
      if (cancelled) return
      setRanks([])
      setRank('')
      setRanksError(errorMessage(err))
    })
    return () => { cancelled = true }
  }, [study, runs, option])

  // A result answers the inputs it was computed from; new inputs discard it,
  // along with any request still in flight for the old ones.
  useEffect(() => {
    request.current++
    setResult(null)
    setLoading(false)
  }, [study, runs, option?.table, rank, aggregate])

  /** Fit the model for the current inputs and show it, unless the inputs changed meanwhile. */
  const run = async () => {
    if (!option || !rank) return
    const req = ++request.current
    setLoading(true)
    try {
      const res = await api.analysis.differential(study, { runs, table: option.table, rank, aggregate })
      if (req === request.current) setResult(res)
    } catch (err) {
      if (req !== request.current) return
      setResult(null)
      toast.error(`Differential abundance failed: ${errorMessage(err)}`)
    } finally {
      if (req === request.current) setLoading(false)
    }
  }

  const title = result
    ? `Differential abundance: ${result.groups.contrast} vs ${result.groups.reference} (${result.rank})`
    : ''

  return (
    <div style={{ marginTop: 16 }}>
      <div style={{ display: 'flex', gap: 8, alignItems: 'center', flexWrap: 'wrap', marginBottom: 12 }}>
        <button className="btn" onClick={run} disabled={loading || !option || !rank || runs.length !== 2}>
          {loading ? 'Fitting…' : 'Differential Abundance'}
        </button>
        {ranks.length > 0 && (
          <select aria-label="Rank" value={rank} onChange={e => setRank(e.target.value)}
            style={{ font: 'inherit', padding: '1px 4px', verticalAlign: 'baseline' }}>
            {ranks.map(r => <option key={r} value={r}>{r}</option>)}
          </select>
        )}
        {ranksReady && ranks.length === 0 && (
          <span style={{ fontSize: '.82rem', color: 'var(--color-muted-fg)' }}>No taxonomy rank shared by the selected runs.</span>
        )}
        {ranksError && (
          <span className="error-msg" style={{ margin: 0 }}>Could not read the taxonomy ranks: {ranksError}</span>
        )}
        {runs.length === 2 && (
          <span style={{ fontSize: '.82rem', color: 'var(--color-muted-fg)' }}>
            Reference: {runs[0].prefix ?? runs[0].run}; positive fold changes mean more abundant in {runs[1].prefix ?? runs[1].run}.
          </span>
        )}
      </div>

      {result && (
        <>
          <p style={{ fontSize: '.82rem', color: 'var(--color-muted-fg)', margin: '0 0 8px' }}>
            {result.method}. {result.diagnostics.n_tested} of {result.diagnostics.n_taxa} taxa tested
            {result.diagnostics.n_failed > 0 && `, ${result.diagnostics.n_failed} could not be fitted`}
            {result.diagnostics.n_boundary > 0 && `, ${result.diagnostics.n_boundary} at a dispersion boundary`}
            {result.diagnostics.n_filtered > 0 && `, ${result.diagnostics.n_filtered} below the prevalence threshold`}
            . Samples: {Object.entries(result.n_samples).map(([g, n]) => `${g} ${n}`).join(', ')}.
          </p>

          <AnalysisChart study={study} figure={result.figure} />

          <div style={{ display: 'flex', gap: 6, justifyContent: 'flex-end', margin: '8px 0' }}>
            <button className="btn" onClick={() => saveBlob(csvBlob(result), 'differential_abundance.csv')}>
              Download CSV
            </button>
            <AddToReport study={study} kind="table" className="btn" defaultTitle={title}
              make={async () => ({ blob: csvBlob(result), ext: '.csv' })} />
          </div>

          <div style={{ overflowX: 'auto' }}>
            <table style={{ borderCollapse: 'collapse', fontSize: '.82rem', width: '100%' }}>
              <thead>
                <tr style={{ borderBottom: '2px solid var(--color-border)' }}>
                  <th style={{ ...head, textAlign: 'left' }}>{result.rank}</th>
                  <th style={head}>log2 FC</th>
                  <th style={head}>Estimate (ln)</th>
                  <th style={head}>SE</th>
                  <th style={head}>p</th>
                  <th style={head}>padj (BH)</th>
                  <th style={{ ...head, textAlign: 'left' }}>Status</th>
                </tr>
              </thead>
              <tbody>
                {result.rows.map(r => (
                  <tr key={r.taxon} style={{ borderBottom: '1px solid var(--color-border)' }}>
                    <td style={{ ...cell, textAlign: 'left' }}>{r.taxon}</td>
                    <td style={cell}>{formatStat(r.log2_fold_change)}</td>
                    <td style={cell}>{formatStat(r.estimate)}</td>
                    <td style={cell}>{formatStat(r.standard_error)}</td>
                    <td style={cell}>{formatStat(r.pvalue)}</td>
                    <td style={cell}>{formatStat(r.padj)}</td>
                    <td style={{ ...cell, textAlign: 'left' }} title={r.note || undefined}>
                      {r.status}{r.note && <span style={{ color: 'var(--color-muted-fg)' }}> – {r.note}</span>}
                    </td>
                  </tr>
                ))}
              </tbody>
            </table>
          </div>
        </>
      )}
    </div>
  )
}
