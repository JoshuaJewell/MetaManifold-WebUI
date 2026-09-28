import { useCallback, useEffect } from 'react'
import { Link } from 'react-router-dom'
import { useApi } from '../hooks/useApi'
import { useJobRefetch } from '../hooks/useJobEvents'
import { api } from '../api/client'
import { JobBadge } from '../components/JobBadge'
import { Skeleton } from '../components/Skeleton'
import { parseUtc, timeAgo } from '../utils/timeago'
import type { Job } from '../api/types'
import { useToast } from '../components/Toast'

export function JobsView() {
  const fetcher = useCallback(() => api.jobs.list(), [])
  const { data: jobs, loading, error, refetch } = useApi(fetcher)
  const toast = useToast()

  useJobRefetch(refetch, {})

  // Refetch on window focus (catches jobs completed while tab was backgrounded)
  useEffect(() => {
    const onFocus = () => refetch()
    window.addEventListener('focus', onFocus)
    return () => window.removeEventListener('focus', onFocus)
  }, [refetch])

  const cancel = async (id: string) => {
    if (!window.confirm(`Cancel job ${id}?`)) return
    try {
      await api.jobs.cancel(id)
      toast.success(`Cancelled ${id}`)
    } catch (e) {
      toast.error(`Cancel failed: ${(e as Error).message}`)
    }
    refetch()
  }

  return (
    <>
      <div className="page-header">
        <h1>Jobs</h1>
        <p>Pipeline execution queue.</p>
      </div>

      {loading && <Skeleton lines={3} />}
      {error   && <p className="error-msg">{error}</p>}

      {jobs && jobs.length === 0 && (
        <div className="empty-state">No jobs have been submitted yet.</div>
      )}

      {jobs && jobs.length > 0 && (
        <div style={{ display: 'flex', flexDirection: 'column', gap: 4 }}>
          {jobs.map((job: Job) => (
            <div key={job.id} className="card" style={{ display: 'grid', gridTemplateColumns: '1fr auto auto', alignItems: 'center', gap: 12 }}>
              <div>
                <strong style={{ fontSize: '.9rem' }}>{job.id}</strong>
                <div style={{ fontSize: '.82rem', color: 'var(--color-muted-fg)' }}>
                  {job.type}
                  {job.study && <>{' · '}<Link to={`/${job.study}`}>{job.study}</Link></>}
                  {job.study && job.run && <>{' / '}<Link to={job.group ? `/${job.study}/${job.group}/${job.run}` : `/${job.study}/${job.run}`}>{job.group ? `${job.group}/${job.run}` : job.run}</Link></>}
                  {job.stage && ` · ${job.stage}`}
                  {job.message && ` · ${job.message}`}
                </div>
                <div style={{ fontSize: '.75rem', color: 'var(--color-muted-fg)' }}>
                  <span title={parseUtc(job.created_at).toLocaleString()}>
                    {timeAgo(job.created_at)}
                  </span>
                  {job.finished_at && (
                    <>
                      {' · finished '}
                      <span title={parseUtc(job.finished_at).toLocaleString()}>
                        {timeAgo(job.finished_at)}
                      </span>
                    </>
                  )}
                </div>
              </div>
              <JobBadge status={job.status} />
              {(job.status === 'queued' || job.status === 'running') && (
                <button className="btn" onClick={() => cancel(job.id)}>Cancel</button>
              )}
            </div>
          ))}
        </div>
      )}
    </>
  )
}
