// © 2026 Joshua Benjamin Jewell. All rights reserved.
// Licensed under the GNU Affero General Public License version 3 (AGPLv3).
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
