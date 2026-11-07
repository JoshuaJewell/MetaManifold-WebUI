// SPDX-FileCopyrightText: 2026 Joshua Benjamin Jewell
// SPDX-FileCopyrightText: 2026 Jonathan D.A. Jewell (hyperpolymath) <j.d.a.jewell@open.ac.uk>
// SPDX-License-Identifier: MPL-2.0

import type { TreeOp } from './model'
import { DEFAULT_SETTINGS, type TreeSettings } from './layout'

export interface TreeViewDoc {
  version:   1
  source:    { file: string; sha256: string }
  ops:       TreeOp[]
  settings:  TreeSettings
  scroll?:   { left: number; top: number }
  selected?: string | null
}

export function mergeSettings(raw: unknown): TreeSettings {
  const saved = (raw && typeof raw === 'object' ? raw : {}) as Omit<Partial<TreeSettings>, 'placements'> & { placements?: string }
  const sym = (saved.symbols ?? {}) as Partial<TreeSettings['symbols']>
  const out: TreeSettings = {
    ...DEFAULT_SETTINGS, ...saved,
    placements: saved.placements === 'all' ? 'all' : 'best',
    bestLwrThreshold: typeof saved.bestLwrThreshold === 'number' && Number.isFinite(saved.bestLwrThreshold)
      ? Math.min(1, Math.max(0, saved.bestLwrThreshold)) : DEFAULT_SETTINGS.bestLwrThreshold,
    symbols: {
      placements: { ...DEFAULT_SETTINGS.symbols.placements, ...(sym.placements ?? {}) },
      support:    { ...DEFAULT_SETTINGS.symbols.support, ...(sym.support ?? {}) },
    },
  }
  // Views saved before symbols existed hid placements with placements: 'none'.
  if (saved.placements === 'none' && !sym.placements) out.symbols.placements.show = false
  return out
}

export function readDoc(raw: unknown): Partial<TreeViewDoc> {
  if (!raw || typeof raw !== 'object') return {}
  const d = raw as Partial<TreeViewDoc>
  return {
    ops:      Array.isArray(d.ops) ? d.ops : [],
    settings: mergeSettings(d.settings),
    scroll:   d.scroll,
    selected: d.selected ?? null,
    source:   d.source,
  }
}
