import { Link } from 'react-router-dom'

export function NotFoundView() {
  return (
    <div style={{ textAlign: 'center', padding: '80px 20px' }}>
      <div style={{ fontSize: '4rem', fontWeight: 700, color: 'var(--color-muted-fg)', lineHeight: 1 }}>404</div>
      <p style={{ fontSize: '1.1rem', color: 'var(--color-muted-fg)', margin: '12px 0 24px' }}>
        This page doesn't exist.
      </p>
      <Link to="/studies" className="btn btn-primary" style={{ textDecoration: 'none' }}>
        Back to Studies
      </Link>
    </div>
  )
}
