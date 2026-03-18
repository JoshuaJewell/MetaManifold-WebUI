import { useState, useEffect, useCallback } from 'react'

interface State<T> {
  data:    T | null
  loading: boolean
  error:   string | null
}

/**
 * Minimal data-fetching hook. Re-fetches when `fetcher` reference changes.
 * Returns `{ data, loading, error, refetch }`.
 */
export function useApi<T>(fetcher: () => Promise<T>): State<T> & { refetch: () => void } {
  const [state, setState] = useState<State<T>>({ data: null, loading: true, error: null })
  const [tick, setTick]   = useState(0)

  const refetch = useCallback(() => setTick(t => t + 1), [])

  useEffect(() => {
    let cancelled = false
    setState(s => ({ ...s, loading: true, error: null }))
    fetcher()
      .then(data  => { if (!cancelled) setState({ data, loading: false, error: null }) })
      .catch(err  => { if (!cancelled) setState(s => ({ ...s, loading: false, error: String(err.message ?? err) })) })
    return () => { cancelled = true }
  // eslint-disable-next-line react-hooks/exhaustive-deps
  }, [tick, fetcher])

  return { ...state, refetch }
}
