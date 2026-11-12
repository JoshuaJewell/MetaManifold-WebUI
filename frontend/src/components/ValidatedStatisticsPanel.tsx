// SPDX-FileCopyrightText: 2026 Jonathan D.A. Jewell (hyperpolymath) <j.d.a.jewell@open.ac.uk>
// SPDX-License-Identifier: MPL-2.0

import { useEffect, useRef, useState } from 'react'
import { api } from '../api/client'
import type { ComparisonRunSpec } from '../api/types'
import type { AnalysisOption } from './annotationShared'
import { saveBlob } from '../utils/download'
import { exactIntegerText, proportionText, type DescriptiveResult } from './validatedStatistics'

export function ValidatedStatisticsPanel({ study, runs, option }: {
  study: string; runs: ComparisonRunSpec[]; option: AnalysisOption | null
}) {
  const [runIndex, setRunIndex] = useState(0)
  const [result, setResult] = useState<DescriptiveResult | null>(null)
  const [sampleIndex, setSampleIndex] = useState(0)
  const [page, setPage] = useState(0)
  const [error, setError] = useState<string | null>(null)
  const [busy, setBusy] = useState(false)
  const sequence = useRef(0)
  const run = runs[runIndex]
  const scope = JSON.stringify([study, run, option?.table])
  useEffect(() => {
    sequence.current++
    setResult(null); setError(null); setBusy(false); setSampleIndex(0); setPage(0)
    return () => { sequence.current++ }
  }, [scope])
  useEffect(() => { if (runIndex >= runs.length) setRunIndex(0) }, [runIndex, runs.length])

  async function calculate() {
    if (!run || !option) return
    const id = ++sequence.current
    setBusy(true); setError(null); setResult(null)
    try {
      const response = await api.statistics.run(study, run.run, {
        table: option.table, prefix: run.prefix ?? null,
        numeric_policy: { mode: 'exact_counts' },
        time_limit_seconds: 60,
      }, run.group)
      // Refuse a malformed count before it can appear as a rounded integer.
      response.samples.forEach(sample => {
        exactIntegerText(sample.total)
        sample.counts?.forEach(exactIntegerText)
        sample.proportions?.forEach(proportionText)
      })
      if (id === sequence.current) { setResult(response); setSampleIndex(0); setPage(0) }
    } catch (e) {
      if (id === sequence.current) setError(e instanceof Error ? e.message : String(e))
    } finally {
      if (id === sequence.current) setBusy(false)
    }
  }
  const sample = result?.samples[sampleIndex]
  const features = sample?.features ?? []
  return <details>
    <summary>Exact descriptive summaries</summary>
    <p>Summarise one run or sub-group from its stored integer counts. These summaries use all rows in the selected table; ordinary diversity normalisation and category exclusions do not apply. They do not test differences between samples.</p>
    <p>Independent scientific review and release approval are pending.</p>
    <div style={{ display: 'flex', gap: 12, flexWrap: 'wrap', alignItems: 'center' }}>
      <label>Run <select value={runIndex} onChange={e => setRunIndex(Number(e.target.value))}>
        {runs.map((r, i) => <option key={JSON.stringify([r.group, r.run, r.prefix])} value={i}>{r.group ? `${r.group} / ` : ''}{r.run}{r.prefix ? ` / ${r.prefix}` : ''}</option>)}
      </select></label>
      <button type="button" disabled={busy || !option || !run} onClick={calculate}>{busy ? 'Computing…' : 'Compute summary'}</button>
    </div>
    {error && <p role="alert" className="error-msg">{error}</p>}
    {result && <>
      <p>{result.claim}</p>
      <button type="button" onClick={() => saveBlob(new Blob([JSON.stringify(result, null, 2)], { type: 'application/json' }), `descriptive-${result.analysis_id}.json`)}>Download complete result</button>
      <p><label>Sample <select value={sampleIndex} onChange={e => { setSampleIndex(Number(e.target.value)); setPage(0) }}>
        {result.samples.map((s, i) => <option key={s.label} value={i}>{s.label}</option>)}
      </select></label></p>
      {sample && <>
        <p>Total reads: <strong>{exactIntegerText(sample.total)}</strong>{sample.richness !== undefined && <> · Observed richness: {sample.richness}</>}</p>
        {sample.warnings.map(w => <p key={w}>{w}</p>)}
        <table><thead><tr><th>Feature</th><th>Exact count</th><th>Exact proportion</th></tr></thead>
          <tbody>{features.slice(page * 25, (page + 1) * 25).map((feature, offset) => {
            const i = page * 25 + offset
            return <tr key={feature}><td style={{ overflowWrap: 'anywhere' }}>{feature}</td><td>{exactIntegerText(sample.counts![i])}</td><td>{proportionText(sample.proportions![i])}</td></tr>
          })}</tbody>
        </table>
        {features.length > 25 && <p><button disabled={page === 0} onClick={() => setPage(page - 1)}>Previous</button> Page {page + 1} of {Math.ceil(features.length / 25)} <button disabled={(page + 1) * 25 >= features.length} onClick={() => setPage(page + 1)}>Next</button></p>}
      </>}
    </>}
  </details>
}
