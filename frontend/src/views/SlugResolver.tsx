import { useParams } from 'react-router-dom'
import { useCallback } from 'react'
import { useApi } from '../hooks/useApi'
import { api } from '../api/client'
import { Skeleton } from '../components/Skeleton'
import { GroupView } from './GroupView'
import { RunView } from './RunView'
import { NotFoundView } from './NotFoundView'

export function SlugResolver() {
  const { study, slug } = useParams<{ study: string; slug: string }>()
  const fetcher = useCallback(() => api.studies.get(study!), [study])
  const { data: detail, loading } = useApi(fetcher)

  if (loading || !detail) return <Skeleton lines={4} />

  if (detail.groups?.includes(slug!)) {
    return <GroupView groupName={slug!} />
  }

  if (detail.runs?.includes(slug!)) {
    return <RunView runName={slug!} />
  }

  return <NotFoundView />
}
