import { Link, useParams } from 'react-router-dom'

export function Breadcrumb() {
  const { study, group, run, slug } = useParams<{ study?: string; group?: string; run?: string; slug?: string }>()

  const crumbs: { label: string; to: string }[] = [
    { label: 'Studies', to: '/studies' },
  ]

  if (study) {
    crumbs.push({ label: study, to: `/${study}` })

    if (group) {
      crumbs.push({ label: group, to: `/${study}/${group}` })
    }

    if (run) {
      crumbs.push({ label: run, to: group ? `/${study}/${group}/${run}` : `/${study}/${run}` })
    }

    if (slug && !group && !run) {
      crumbs.push({ label: slug, to: `/${study}/${slug}` })
    }
  }

  if (crumbs.length <= 1) return null

  return (
    <nav style={{
      fontSize: '.82rem',
      color: 'var(--color-muted-fg)',
      marginBottom: 12,
      display: 'flex',
      alignItems: 'center',
      gap: 6,
    }}>
      {crumbs.map((c, i) => {
        const isLast = i === crumbs.length - 1
        return (
          <span key={c.to} style={{ display: 'flex', alignItems: 'center', gap: 6 }}>
            {i > 0 && <span style={{ opacity: .4 }}>/</span>}
            {isLast
              ? <span style={{ color: 'var(--color-fg)', fontWeight: 500 }}>{c.label}</span>
              : <Link to={c.to} style={{ color: 'var(--color-muted-fg)' }}>{c.label}</Link>
            }
          </span>
        )
      })}
    </nav>
  )
}
