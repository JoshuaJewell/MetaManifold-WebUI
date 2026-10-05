// SPDX-License-Identifier: AGPL-3.0-only
// SPDX-FileCopyrightText: 2026 Jonathan D.A. Jewell (hyperpolymath) <j.d.a.jewell@open.ac.uk>
import { describe, expect, test } from 'bun:test'
import type { ComparisonRunSpec } from '../api/types'
import { sharedRanks } from './sharedRanks'

const runs: ComparisonRunSpec[] = [
  { run: 'A', group: null, prefix: 'x' },
  { run: 'A', group: null, prefix: 'y' },
  { run: 'B', group: 'g', prefix: null },
]
const RANKS: Record<string, string[]> = {
  A: ['Kingdom', 'Phylum', 'Class', 'Genus'],
  B: ['Kingdom', 'Phylum', 'Genus', 'Species'],
}

describe('sharedRanks', () => {
  test('intersects the ranks of each distinct run, in the first run\'s order', async () => {
    const asked: string[] = []
    const ranks = await sharedRanks(runs, async r => { asked.push(r.run); return RANKS[r.run] })
    expect(ranks).toEqual(['Kingdom', 'Phylum', 'Genus'])
    expect(asked).toEqual(['A', 'B'])
  })

  test('no shared rank is an empty list', async () => {
    expect(await sharedRanks(runs, async r => r.run === 'A' ? ['Genus'] : ['Species'])).toEqual([])
  })

  test('a failed request rejects instead of reading as no shared rank', async () => {
    const failing = sharedRanks(runs, async r => {
      if (r.run === 'B') throw new Error('HTTP 500: R error')
      return RANKS[r.run]
    })
    await expect(failing).rejects.toThrow('HTTP 500: R error')
  })
})
