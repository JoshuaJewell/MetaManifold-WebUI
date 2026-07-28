// © 2026 Joshua Benjamin Jewell. All rights reserved.
// Licensed under the GNU Affero General Public License version 3 (AGPLv3).
import type { ColFilter } from '../../api/types'
import { SAMPLE_READS_FILTER_KEY } from '../../api/types'

// Reads a saved filter preset: a `filters:` list of {column, type, values | value}
// entries, plus an optional {type: sample_reads, min, max, basis} entry.
export function parseFilterYaml(text: string): Record<string, ColFilter> {
  const scalar = (raw: string): string => {
    const v = raw.trim()
    if ((v.startsWith('"') && v.endsWith('"')) || (v.startsWith("'") && v.endsWith("'"))) return v.slice(1, -1)
    return v
  }
  const flowList = (raw: string): string[] | null => {
    const v = raw.trim()
    if (!v.startsWith('[') || !v.endsWith(']')) return null
    const inner = v.slice(1, -1).trim()
    return inner === '' ? [] : inner.split(',').map(scalar)
  }

  type Entry = { column?: string; type?: string; values?: string[]; value?: string; min?: string; max?: string; basis?: string }
  const entries: Entry[] = []
  let cur: Entry | null = null
  let listKey: 'values' | null = null

  for (const rawLine of text.split(/\r?\n/)) {
    const line = rawLine.replace(/\s+#.*$/, '')
    const trimmed = line.trim()
    if (trimmed === '' || trimmed.startsWith('#') || /^filters:\s*(\[\])?$/.test(trimmed)) continue

    const item = trimmed.match(/^-\s+(\w+):\s*(.*)$/)
    const bare = trimmed.match(/^-\s+(.+)$/)
    const field = trimmed.match(/^(\w+):\s*(.*)$/)

    if (item && line.search(/\S/) <= 2) {
      cur = {}
      entries.push(cur)
      listKey = null
      ;(cur as Record<string, unknown>)[item[1]] = scalar(item[2])
    } else if (bare && listKey && cur) {
      cur.values!.push(scalar(bare[1]))
    } else if (field && cur) {
      const [, key, value] = field
      if (key === 'values') {
        const flow = flowList(value)
        cur.values = flow ?? []
        listKey = flow ? null : 'values'
      } else {
        listKey = null
        ;(cur as Record<string, unknown>)[key] = scalar(value)
      }
    } else {
      throw new Error(`Unrecognised line: ${trimmed}`)
    }
  }

  const num = (v?: string) => {
    if (v == null || v === '') return undefined
    const n = Number(v)
    if (!Number.isFinite(n)) throw new Error(`Not a number: ${v}`)
    return n
  }
  const result: Record<string, ColFilter> = {}
  for (const e of entries) {
    if (e.type === 'sample_reads') {
      const f: ColFilter = {}
      const min = num(e.min), max = num(e.max)
      if (min != null) f.min = min
      if (max != null) f.max = max
      if (e.basis === 'raw' || e.basis === 'filtered') f.basis = e.basis
      result[SAMPLE_READS_FILTER_KEY] = f
      continue
    }
    if (!e.column) throw new Error('Filter entry without a column')
    const f = result[e.column] ?? (result[e.column] = {})
    switch (e.type) {
      case 'include': f.include = [...(f.include ?? []), ...(e.values ?? [])]; break
      case 'exclude': f.exclude = [...(f.exclude ?? []), ...(e.values ?? [])]; break
      case 'min':     f.min = num(e.value); break
      case 'max':     f.max = num(e.value); break
      default:        throw new Error(`Unknown filter type: ${e.type}`)
    }
  }
  return result
}
