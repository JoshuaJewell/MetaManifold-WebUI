import React, { useCallback, useEffect, useMemo, useState } from 'react'
import { NavLink, Outlet, useLocation, useParams } from 'react-router-dom'
import { useApi } from '../hooks/useApi'
import { useSSE } from '../hooks/useSSE'
import { api } from '../api/client'
import { createJobEventBus, JobEventContext, SSEConnectedContext } from '../hooks/useJobEvents'
import { Breadcrumb } from '../components/Breadcrumb'
import { ErrorBoundary } from '../components/ErrorBoundary'
import type { Job, StudySummary } from '../api/types'
import { COPYRIGHT, LICENSE_URL, SOURCE_URL } from '../about'
import { NavRefreshContext } from '../hooks/useNavRefresh'

export function Layout() {
  const { study, group: groupParam, slug } = useParams<{ study?: string; group?: string; slug?: string }>()
  const { pathname } = useLocation()
  const { data: studies, refetch: refetchStudies } = useApi(api.studies.list)
  const [navOpen, setNavOpen] = useState(false)
  useEffect(() => { setNavOpen(false) }, [pathname])
  const [runningJobIds, setRunningJobIds] = useState<Set<string>>(new Set())

  const jobBus = useMemo(() => createJobEventBus(), [])

  useEffect(() => {
    api.jobs.list().then(jobs => {
      setRunningJobIds(new Set(jobs.filter((j: Job) => j.status === 'running').map((j: Job) => j.id)))
    }).catch(() => {})
  }, [])

  const sseHandlers = useMemo(() => ({
    onJobUpdate: (job: Job) => {
      setRunningJobIds(prev => {
        const next = new Set(prev)
        if (job.status === 'running') next.add(job.id)
        else next.delete(job.id)
        return next
      })
      jobBus.emit(job)
    },
  }), [jobBus])

  const sseConnected = useSSE(sseHandlers)

  const studyFetcher = useCallback(
    () => study ? api.studies.get(study) : Promise.resolve(null),
    [study]
  )
  const { data: studyDetail, refetch: refetchDetail } = useApi(studyFetcher)

  const activeGroup = groupParam ?? (studyDetail?.groups?.includes(slug!) ? slug : undefined)

  const groupRunsFetcher = useCallback(
    () => study && activeGroup ? api.runs.listGroup(study, activeGroup) : Promise.resolve(null),
    [study, activeGroup]
  )
  const { data: groupRuns, refetch: refetchGroupRuns } = useApi(groupRunsFetcher)

  const refreshNav = useCallback(() => {
    refetchStudies()
    refetchDetail()
    refetchGroupRuns()
  }, [refetchStudies, refetchDetail, refetchGroupRuns])

  return (
    <div className={`app-layout ${navOpen ? 'nav-open' : ''}`}>
      <nav className="sidebar" id="sidebar">
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
              <React.Fragment key={g}>
                <NavLink
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
              </React.Fragment>
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
          Jobs{runningJobIds.size > 0 && ` (${runningJobIds.size})`}
        </NavLink>
        <NavLink to="/databases" className={({ isActive }) => `sidebar-link ${isActive ? 'active' : ''}`}>Databases</NavLink>
        <NavLink to="/config" className={({ isActive }) => `sidebar-link ${isActive ? 'active' : ''}`}>Default Config</NavLink>
        <NavLink to="/primers" className={({ isActive }) => `sidebar-link ${isActive ? 'active' : ''}`}>Primers</NavLink>
        <NavLink to="/compositions" className={({ isActive }) => `sidebar-link ${isActive ? 'active' : ''}`}>Compositions</NavLink>
        <NavLink to="/reference-trees" className={({ isActive }) => `sidebar-link ${isActive ? 'active' : ''}`}>Reference trees</NavLink>
        <NavLink to="/about" className={({ isActive }) => `sidebar-link ${isActive ? 'active' : ''}`}>About</NavLink>
        <div className="sidebar-legal">
          {COPYRIGHT} · <a href={LICENSE_URL} target="_blank" rel="noreferrer">AGPL-3.0</a> ·{' '}
          <a href={SOURCE_URL} target="_blank" rel="noreferrer">Source and Documentation</a>
        </div>
      </nav>

      {navOpen && <div className="sidebar-backdrop" onClick={() => setNavOpen(false)} />}

      <main className="main-content">
        <button className="btn btn-sm nav-toggle" aria-controls="sidebar" aria-expanded={navOpen}
          onClick={() => setNavOpen(o => !o)}>☰ Menu</button>
        <NavRefreshContext.Provider value={refreshNav}>
        <JobEventContext.Provider value={jobBus}>
          <SSEConnectedContext.Provider value={sseConnected}>
            <Breadcrumb />
            <ErrorBoundary key={pathname}>
              <Outlet />
            </ErrorBoundary>
          </SSEConnectedContext.Provider>
        </JobEventContext.Provider>
        </NavRefreshContext.Provider>
      </main>
    </div>
  )
}
