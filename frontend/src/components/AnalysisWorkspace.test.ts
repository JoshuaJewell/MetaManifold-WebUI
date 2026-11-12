// SPDX-FileCopyrightText: 2026 Jonathan D.A. Jewell (hyperpolymath) <j.d.a.jewell@open.ac.uk>
// SPDX-License-Identifier: MPL-2.0

import { expect, mock, spyOn, test } from 'bun:test'
import { createElement } from 'react'
import { act, create } from 'react-test-renderer'
import { MemoryRouter } from 'react-router-dom'
import { api } from '../api/client'
import * as download from '../utils/download'
import type { DescriptiveResult } from './validatedStatistics'

// Exercise the real workspace, scope controls, table discovery and descriptive
// panel. Unrelated chart panels need a browser canvas and are outside this test.
mock.module('./DiversityPanel', () => ({ DiversityPanel: () => null }))
mock.module('./TaxaCompositionChart', () => ({ TaxaCompositionChart: () => null }))
mock.module('./VennPanel', () => ({ VennPanel: () => null }))
mock.module('./DifferentialPanel', () => ({ DifferentialPanel: () => null }))
mock.module('./PublicationTablesPanel', () => ({ PublicationTablesPanel: () => null }))
mock.module('../figure/FigureBuilder', () => ({ FigureBuilder: () => null }))
const { AnalysisWorkspace } = await import('./AnalysisWorkspace')

function resultFor(prefix: string | null): DescriptiveResult {
  const total = prefix === 'A' ? 3 : 303
  return {
    mode: 'exact_counts', status: 'success', review_status: 'pending_review',
    claim: 'Exact totals and proportions', analysis_id: prefix ?? 'whole-run', input_sha256: 'fixture',
    samples: [{ label: prefix ?? 'all', total, counts: [total], features: ['f'],
      proportions: ['1/1'], warnings: [], status: 'success' }],
  }
}

test('aggregation controls the request, visible summary and downloaded result', async () => {
  const listTables = spyOn(api.results, 'runTables').mockResolvedValue([{ id: 'merged', label: 'ASVs' }])
  const calculate = spyOn(api.statistics, 'run').mockImplementation(async (_study, _run, body) => resultFor(body.prefix))
  const save = spyOn(download, 'saveBlob').mockImplementation(() => {})
  let renderer: ReturnType<typeof create> | undefined
  try {
    await act(async () => {
      renderer = create(createElement(MemoryRouter, { future: { v7_startTransition: true, v7_relativeSplatPath: true } },
        createElement(AnalysisWorkspace, { study: 'study', source: 'VSEARCH', runs: [
          { group: 'group', run: 'run1', prefix: 'A' },
          { group: 'group', run: 'run1', prefix: 'AB' },
        ] })))
    })
    const root = renderer!.root
    const button = (label: string) => root.findAllByType('button').find(b => b.children.includes(label))!
    const runSelect = () => root.findAllByType('select')[1]
    const toggle = async (checked: boolean) => {
      await act(async () => root.findByProps({ type: 'checkbox' }).props.onChange({ target: { checked } }))
    }
    const compute = async (prefix: string | null, total: string) => {
      await act(async () => button('Compute summary').props.onClick())
      expect(calculate.mock.lastCall).toEqual(['study', 'run1', {
        table: 'merged', prefix, numeric_policy: { mode: 'exact_counts' }, time_limit_seconds: 60,
      }, 'group'])
      expect(root.findByType('strong').children).toEqual([total])
      await act(async () => button('Download complete result').props.onClick())
      const [blob, filename] = save.mock.lastCall!
      expect(JSON.parse(await blob.text())).toEqual(resultFor(prefix))
      expect(filename).toBe(`descriptive-${prefix ?? 'whole-run'}.json`)
    }

    expect(runSelect().findAllByType('option')).toHaveLength(2)
    await compute('A', '3')
    await toggle(true)
    expect(runSelect().findAllByType('option')).toHaveLength(1)
    expect(button('Download complete result')).toBeUndefined()
    await compute(null, '303')
    await toggle(false)
    expect(runSelect().findAllByType('option')).toHaveLength(2)
    expect(button('Download complete result')).toBeUndefined()
    await compute('A', '3')

    // A response started under the old subgroup cannot reappear after pooling.
    let finish!: (result: DescriptiveResult) => void
    calculate.mockImplementationOnce(() => new Promise(resolve => { finish = resolve }))
    await act(async () => { void button('Compute summary').props.onClick() })
    await toggle(true)
    await act(async () => finish(resultFor('A')))
    expect(button('Download complete result')).toBeUndefined()
    expect(root.findAllByType('strong')).toHaveLength(0)
    await compute(null, '303')
  } finally {
    await act(async () => renderer?.unmount())
    listTables.mockRestore(); calculate.mockRestore(); save.mockRestore()
  }
})
