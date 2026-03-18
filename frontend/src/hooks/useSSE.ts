import { useEffect, useRef } from 'react'
import { openEventStream, type SSEHandler } from '../api/events'

/**
 * Subscribe to the global SSE stream. Pass a stable handlers object
 * (e.g. via useMemo) to avoid reconnecting on every render.
 */
export function useSSE(handlers: SSEHandler) {
  const ref = useRef(handlers)
  ref.current = handlers

  useEffect(() => {
    return openEventStream({
      onJobUpdate:   (job)  => ref.current.onJobUpdate?.(job),
      onStageUpdate: (ev)   => ref.current.onStageUpdate?.(ev),
    })
  }, [])
}
