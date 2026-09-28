// © 2026 Joshua Benjamin Jewell. All rights reserved.
// Licensed under the GNU Affero General Public License version 3 (AGPLv3).
import { useState } from 'react'
import {
  STAGE_CONFIG_PREFIXES,
  STAGE_LABELS,
  StageConfig,
} from './PipelineStages'
import type { ConfigMap, ConfigSource } from '../api/types'
import type { ConfigSection } from './PipelineStages'

export const ALL_SECTIONS = Object.keys(STAGE_CONFIG_PREFIXES) as ConfigSection[]

export function ConfigAccordion({ configMap, study, run, group, onConfigChanged, patchFn, deleteFn, sourceLevel, overrides, sections }: {
  configMap: ConfigMap
  study: string
  run: string
  group?: string
  onConfigChanged: () => void
  patchFn?: (study: string, run: string, body: Record<string, unknown>, group?: string) => Promise<ConfigMap>
  deleteFn?: (study: string, run: string, key: string, group?: string) => Promise<ConfigMap>
  sourceLevel?: ConfigSource
  overrides?: Record<string, string[]> | null
  // Which sections to render; every one by default. The run page passes the
  // handful that are not pipeline stages, because the stage sections are
  // already on that page as the runnable rows.
  sections?: ConfigSection[]
}) {
  const [expanded, setExpanded] = useState<string | null>(null)

  return (
    <div>
      {(sections ?? ALL_SECTIONS).map(stage => {
        const prefixes = STAGE_CONFIG_PREFIXES[stage]
        const hasKeys = prefixes.some(p => Object.keys(configMap).some(k => k.startsWith(p)))
        if (!hasKeys) return null
        const isExpanded = expanded === stage
        return (
          <div key={stage} style={{ marginBottom: 4 }}>
            <button
              type="button"
              className="toggle-btn"
              aria-expanded={isExpanded}
              style={{ fontWeight: 600, fontSize: '.85rem', padding: '4px 0', gap: 6 }}
              onClick={() => setExpanded(isExpanded ? null : stage)}
            >
              <span aria-hidden="true" style={{ fontSize: '.8rem', opacity: .65 }}>{isExpanded ? '▾' : '▸'}</span>
              {STAGE_LABELS[stage]}
            </button>
            {isExpanded && (
              <StageConfig
                configMap={configMap}
                prefixes={prefixes}
                study={study}
                run={run}
                group={group}
                onConfigChanged={onConfigChanged}
                patchFn={patchFn}
                deleteFn={deleteFn}
                sourceLevel={sourceLevel}
                overrides={overrides}
              />
            )}
          </div>
        )
      })}
    </div>
  )
}
