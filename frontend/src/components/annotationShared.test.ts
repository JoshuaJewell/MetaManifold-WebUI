// SPDX-License-Identifier: AGPL-3.0-only
// SPDX-FileCopyrightText: 2026 Jonathan D.A. Jewell (hyperpolymath) <j.d.a.jewell@open.ac.uk>
import { afterEach, describe, expect, spyOn, test } from 'bun:test'
import { api } from '../api/client'
import type { TableMeta } from '../api/types'
import { discoverResultsTables, sharedOptions } from './annotationShared'

const meta = (id: string, label: string) => ({ id, label }) as TableMeta

describe('discoverResultsTables', () => {
  let spy: ReturnType<typeof spyOn> | null = null
  afterEach(() => { spy?.mockRestore(); spy = null })

  test('lists the run\'s tables sorted by label', async () => {
    spy = spyOn(api.results, 'runTables').mockResolvedValue([meta('merged_otu', 'OTUs'), meta('merged', 'Merged')])
    expect(await discoverResultsTables('s', 'r')).toEqual([
      { key: 'merged', table: 'merged', label: 'Merged' },
      { key: 'merged_otu', table: 'merged_otu', label: 'OTUs' },
    ])
  })

  test('a failed lookup rejects instead of reading as no tables', async () => {
    spy = spyOn(api.results, 'runTables').mockRejectedValue(new Error('HTTP 503'))
    await expect(discoverResultsTables('s', 'r')).rejects.toThrow('HTTP 503')
  })
})

describe('sharedOptions', () => {
  const o = (key: string, label = key) => ({ key, table: key, label })

  test('keeps only the tables every run has', () => {
    expect(sharedOptions([[o('b'), o('a'), o('c')], [o('c'), o('a')]])).toEqual([o('a'), o('c')])
  })

  test('no runs means no tables', () => {
    expect(sharedOptions([])).toEqual([])
  })
})
