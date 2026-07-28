// © 2026 Joshua Benjamin Jewell. All rights reserved.
// Licensed under the GNU Affero General Public License version 3 (AGPLv3).
import { apiUrl } from '../../api/client'
import type { ConfigMap, RunStages } from '../../api/types'
import { StaleKeysBadge } from './DADA2Panel'

interface QCOutput { has_report: boolean; report_url: string | null }

export function QCPanel({ qcData, stages, onRunStage, configMap, cacheKey }: { qcData: QCOutput | null; stages: RunStages | null; onRunStage: (stage: string) => void; configMap?: ConfigMap | null; cacheKey?: string | null }) {
  const fastqcStatus = stages?.fastqc?.status
  const isStale = fastqcStatus === 'stale'
  const isRunning = fastqcStatus === 'running'

  if (!qcData || !qcData.has_report) {
    return (
      <div className="empty-state">
        <p>No QC report generated yet. Run QC to generate FastQC/MultiQC output on raw reads.</p>
        <button className="btn btn-primary" style={{ marginTop: 12 }} onClick={() => onRunStage('fastqc')} disabled={isRunning}>
          {isRunning ? 'Running…' : 'Run QC'}
        </button>
      </div>
    )
  }

  return (
    <div className="card" style={{ padding: 0, overflow: 'hidden' }}>
      <div style={{ padding: '12px 16px', borderBottom: '1px solid var(--color-border)', display: 'flex', alignItems: 'center', justifyContent: 'space-between' }}>
        <div className="card-title" style={{ margin: 0 }}>MultiQC Report</div>
        <div style={{ display: 'flex', gap: 8, alignItems: 'center' }}>
          {isStale && <StaleKeysBadge staleKeys={stages?.fastqc?.stale_keys ?? []} configMap={configMap} />}
          <button className="btn" style={{ fontSize: '.78rem', padding: '3px 10px' }} onClick={() => onRunStage('fastqc')} disabled={isRunning}>
            {isRunning ? 'Running…' : isStale ? 'Re-run QC' : 'Run QC'}
          </button>
          <a
            href={apiUrl(qcData.report_url!)}
            target="_blank"
            rel="noopener noreferrer"
            className="btn"
            style={{ fontSize: '.78rem', padding: '3px 10px' }}
          >
            Open in new tab
          </a>
        </div>
      </div>
      <iframe
        src={apiUrl(qcData.report_url!) + (cacheKey ? `?v=${encodeURIComponent(cacheKey)}` : '')}
        style={{ width: '100%', height: 'calc(100vh - 240px)', minHeight: 500, border: 'none' }}
        title="MultiQC Report"
      />
    </div>
  )
}
