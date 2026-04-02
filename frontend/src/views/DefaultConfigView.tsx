import { useCallback } from 'react'
import { useApi } from '../hooks/useApi'
import { api } from '../api/client'
import { ConfigAccordion } from '../components/ConfigAccordion'

export function DefaultConfigView() {
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
        <ConfigAccordion
          configMap={configMap}
          study=""
          run=""
          onConfigChanged={refetchConfig}
          patchFn={patchFn}
          deleteFn={deleteFn}
          sourceLevel="default"
        />
      )}
    </>
  )
}
