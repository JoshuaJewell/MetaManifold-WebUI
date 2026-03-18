import { useCallback, useState } from 'react'
import { useApi } from '../hooks/useApi'
import { api } from '../api/client'
import {
  STAGE_CONFIG_PREFIXES,
  STAGE_LABELS,
  StageConfig,
} from '../components/PipelineStages'

const VISIBLE_STAGES = Object.keys(STAGE_CONFIG_PREFIXES) as (keyof typeof STAGE_CONFIG_PREFIXES)[]

export function DefaultConfigView() {
  const [expanded, setExpanded] = useState<string | null>(null)

  const configFetcher = useCallback(() => api.config.getDefault(), [])
  const { data: configMap, refetch: refetchConfig } = useApi(configFetcher)

  const patchFn = useCallback(
    (_study: string, _run: string, body: Record<string, unknown>) => api.config.patchDefault(body),
    []
  )

  const deleteFn = useCallback(
    (_study: string, _run: string, key: string) => api.config.deleteDefault(key),
    []
  )

  return (
    <>
      <div className="page-header" style={{ marginBottom: 20 }}>
        <h1>Default Config</h1>
        <p>Global pipeline defaults. These apply to all studies unless overridden at the study, group, or run level.</p>
      </div>

      {configMap && (
        <div>
          {VISIBLE_STAGES.map(stage => {
            const prefixes = STAGE_CONFIG_PREFIXES[stage]
            const hasKeys = prefixes.some(p => Object.keys(configMap).some(k => k.startsWith(p)))
            if (!hasKeys) return null
            const isExpanded = expanded === stage
            return (
              <div key={stage} style={{ marginBottom: 4 }}>
                <div
                  style={{ cursor: 'pointer', fontWeight: 600, fontSize: '.85rem', padding: '4px 0' }}
                  onClick={() => setExpanded(isExpanded ? null : stage)}
                >
                  <span style={{ fontSize: '.8rem', marginRight: 6, opacity: .65 }}>{isExpanded ? '▾' : '▸'}</span>
                  {STAGE_LABELS[stage]}
                </div>
                {isExpanded && (
                  <StageConfig
                    configMap={configMap}
                    prefixes={prefixes}
                    study=""
                    run=""
                    onConfigChanged={refetchConfig}
                    patchFn={patchFn}
                    deleteFn={deleteFn}
                    sourceLevel="default"
                    overrides={undefined}
                  />
                )}
              </div>
            )
          })}
        </div>
      )}
    </>
  )
}
