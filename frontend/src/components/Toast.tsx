import { createContext, useContext, useState, useCallback, useRef, type ReactNode } from 'react'

type ToastVariant = 'success' | 'error' | 'info'

interface Toast {
  id: number
  message: string
  variant: ToastVariant
}

interface ToastAPI {
  success: (msg: string) => void
  error:   (msg: string) => void
  info:    (msg: string) => void
}

const ToastContext = createContext<ToastAPI | null>(null)

export function useToast(): ToastAPI {
  const ctx = useContext(ToastContext)
  if (!ctx) throw new Error('useToast must be inside ToastProvider')
  return ctx
}

export function ToastProvider({ children }: { children: ReactNode }) {
  const [toasts, setToasts] = useState<Toast[]>([])
  const nextId = useRef(0)

  const add = useCallback((message: string, variant: ToastVariant) => {
    const id = nextId.current++
    setToasts(t => [...t, { id, message, variant }])
    if (variant !== 'error') {
      setTimeout(() => setToasts(t => t.filter(x => x.id !== id)), 4000)
    }
  }, [])

  const api = useRef<ToastAPI>({
    success: (msg) => add(msg, 'success'),
    error:   (msg) => add(msg, 'error'),
    info:    (msg) => add(msg, 'info'),
  })
  api.current.success = (msg) => add(msg, 'success')
  api.current.error   = (msg) => add(msg, 'error')
  api.current.info    = (msg) => add(msg, 'info')

  return (
    <ToastContext.Provider value={api.current}>
      {children}
      <div style={{
        position: 'fixed', bottom: 16, right: 16, zIndex: 9999,
        display: 'flex', flexDirection: 'column-reverse', gap: 8,
        pointerEvents: 'none',
      }}>
        {toasts.map(t => (
          <div key={t.id} style={{
            pointerEvents: 'auto',
            padding: '10px 18px',
            borderRadius: 8,
            fontSize: '.88rem',
            fontWeight: 500,
            boxShadow: '0 4px 12px rgba(0,0,0,.15)',
            animation: 'toast-in .2s ease-out',
            ...VARIANT_STYLES[t.variant],
          }}
            onClick={() => setToasts(ts => ts.filter(x => x.id !== t.id))}
          >
            {t.message}
          </div>
        ))}
      </div>
    </ToastContext.Provider>
  )
}

const VARIANT_STYLES: Record<ToastVariant, React.CSSProperties> = {
  success: { background: '#2b8a3e', color: '#fff' },
  error:   { background: '#c92a2a', color: '#fff' },
  info:    { background: '#1864ab', color: '#fff' },
}
