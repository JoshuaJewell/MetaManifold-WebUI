import { useState, useEffect, useCallback, useRef } from 'react'

interface State<T> {
  data:    T | null
  loading: boolean
  error:   string | null
}

/**
 * Minimal data-fetching hook. Re-fetches when `fetcher` reference changes.
 * Data is kept across refetches unless `resetKey` changes.
 * Returns `{ data, loading, error, refetch }`.
 */
export function useApi<T>(fetcher: () => Promise<T>, resetKey?: unknown): State<T> & { refetch: () => void } {
  const [state, setState] = useState<State<T>>({ data: null, loading: true, error: null })
  const [tick, setTick]   = useState(0)

  const refetch = useCallback(() => setTick(t => t + 1), [])
  const lastKey = useRef(resetKey)

  useEffect(() => {
    let cancelled = false
    // A new resetKey names a different resource, so its old data is cleared while loading.
    const changed = lastKey.current !== resetKey
    lastKey.current = resetKey
    setState(s => ({ data: changed ? null : s.data, loading: true, error: null }))
    fetcher()
      .then(data  => { if (!cancelled) setState({ data, loading: false, error: null }) })
      .catch(err  => { if (!cancelled) setState(s => ({ ...s, loading: false, error: String(err.message ?? err) })) })
    return () => { cancelled = true }
  }, [tick, fetcher, resetKey])

  return { ...state, refetch }
}
