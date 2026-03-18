import { NavLink, Outlet, useParams } from 'react-router-dom'
import { useCallback, useMemo, useState } from 'react'
import { useApi } from '../hooks/useApi'
import { useSSE } from '../hooks/useSSE'
import { api } from '../api/client'
import { createJobEventBus, JobEventContext } from '../hooks/useJobEvents'
import { Breadcrumb } from '../components/Breadcrumb'
import type { Job, StudySummary } from '../api/types'

export function Layout() {
  const { study, group: groupParam, slug } = useParams<{ study?: string; group?: string; slug?: string }>()
  const { data: studies } = useApi(api.studies.list)
  const [activeJobs, setActiveJobs] = useState(0)

  const jobBus = useMemo(() => createJobEventBus(), [])

  const sseHandlers = useMemo(() => ({
    onJobUpdate: (job: Job) => {
      if (job.status === 'running')
        setActiveJobs(n => n + 1)
      if (job.status === 'complete' || job.status === 'failed' || job.status === 'cancelled')
        setActiveJobs(n => Math.max(0, n - 1))
      jobBus.emit(job)
    },
  }), [jobBus])

  useSSE(sseHandlers)

  const studyFetcher = useCallback(
    () => study ? api.studies.get(study) : Promise.resolve(null),
    [study]
  )
  const { data: studyDetail } = useApi(studyFetcher)

  const activeGroup = groupParam ?? (studyDetail?.groups?.includes(slug!) ? slug : undefined)

  const groupRunsFetcher = useCallback(
    () => study && activeGroup ? api.runs.listGroup(study, activeGroup) : Promise.resolve(null),
    [study, activeGroup]
  )
  const { data: groupRuns } = useApi(groupRunsFetcher)

  return (
    <div className="app-layout">
      <nav className="sidebar">
        <div className="sidebar-brand">MetaManifold</div>

        <NavLink to="/studies" end className={({ isActive }) => `sidebar-section sidebar-section-link ${isActive ? 'active' : ''}`}>
          Studies
        </NavLink>
        {studies?.map((s: StudySummary) => (
          <NavLink
            key={s.name}
            to={`/${s.name}`}
            className={({ isActive }) => `sidebar-link sidebar-link-indented ${isActive ? 'active' : ''}`}
          >
            {s.name}
            {s.active_job_count > 0 && <span> ({s.active_job_count})</span>}
          </NavLink>
        ))}

        {study && studyDetail && (
          <>
            {studyDetail.groups && studyDetail.groups.length > 0 && (
              <div className="sidebar-sub-label">Groups</div>
            )}
            {studyDetail.groups?.map((g: string) => (
              <>
                <NavLink
                  key={g}
                  to={`/${study}/${g}`}
                  className={({ isActive }) => `sidebar-link sidebar-link-indented2 ${isActive ? 'active' : ''}`}
                >
                  {g}
                </NavLink>
                {activeGroup === g && groupRuns?.map(r => (
                  <NavLink
                    key={r.name}
                    to={`/${study}/${g}/${r.name}`}
                    className={({ isActive }) => `sidebar-link sidebar-link-indented3 ${isActive ? 'active' : ''}`}
                  >
                    {r.name}
                  </NavLink>
                ))}
              </>
            ))}
            {studyDetail.runs && studyDetail.runs.length > 0 && (
              <div className="sidebar-sub-label">Runs</div>
            )}
            {studyDetail.runs?.map((r: string) => (
              <NavLink
                key={r}
                to={`/${study}/${r}`}
                className={({ isActive }) => `sidebar-link sidebar-link-indented2 ${isActive ? 'active' : ''}`}
              >
                {r}
              </NavLink>
            ))}
          </>
        )}

        <div className="sidebar-spacer" />

        <div className="sidebar-section">System</div>
        <NavLink to="/jobs" className={({ isActive }) => `sidebar-link ${isActive ? 'active' : ''}`}>
          Jobs{activeJobs > 0 && ` (${activeJobs})`}
        </NavLink>
        <NavLink to="/databases" className={({ isActive }) => `sidebar-link ${isActive ? 'active' : ''}`}>Databases</NavLink>
        <NavLink to="/config" className={({ isActive }) => `sidebar-link ${isActive ? 'active' : ''}`}>Default Config</NavLink>
        <div style={{ height: 12 }} />
      </nav>

      <main className="main-content">
        <JobEventContext.Provider value={jobBus}>
          <Breadcrumb />
          <Outlet />
        </JobEventContext.Provider>
      </main>
    </div>
  )
}
