import { useState, useEffect, useCallback, useRef } from 'react'

export interface State<T> {
  data:    T | null
  loading: boolean
  error:   string | null
}

/**
 * The state after a request fails. The previous data survives only when
 * `sameRequest` says it came from this very request (a refetch of unchanged
 * inputs); data fetched for other inputs, such as other filters, is dropped so
 * it cannot be shown as if it answered the new request.
 */
export function failedState<T>(prev: State<T>, error: string, sameRequest: boolean): State<T> {
  return { data: sameRequest ? prev.data : null, loading: false, error }
}

/**
 * Minimal data-fetching hook. Re-fetches when `fetcher` reference changes.
 * Data is kept across refetches unless `resetKey` changes, and is cleared when
 * a request for different inputs (a new `fetcher`) fails.
 * Returns `{ data, loading, error, refetch }`.
 */
export function useApi<T>(fetcher: () => Promise<T>, resetKey?: unknown): State<T> & { refetch: () => void } {
  const [state, setState] = useState<State<T>>({ data: null, loading: true, error: null })
  const [tick, setTick]   = useState(0)

  const refetch = useCallback(() => setTick(t => t + 1), [])
  const lastKey = useRef(resetKey)
  // The fetcher whose response is in `state.data`, or null when there is none.
  const dataFrom = useRef<(() => Promise<T>) | null>(null)

  useEffect(() => {
    let cancelled = false
    // A new resetKey names a different resource, so its old data is cleared while loading.
    const changed = lastKey.current !== resetKey
    lastKey.current = resetKey
    if (changed) dataFrom.current = null
    setState(s => ({ data: changed ? null : s.data, loading: true, error: null }))
    fetcher()
      .then(data  => {
        if (cancelled) return
        dataFrom.current = fetcher
        setState({ data, loading: false, error: null })
      })
      .catch(err  => {
        if (cancelled) return
        const same = dataFrom.current === fetcher
        if (!same) dataFrom.current = null
        setState(s => failedState(s, String(err.message ?? err), same))
      })
    return () => { cancelled = true }
  }, [tick, fetcher, resetKey])

  return { ...state, refetch }
}
