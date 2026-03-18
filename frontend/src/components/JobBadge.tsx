import type { JobStatus } from '../api/types'
import styles from './JobBadge.module.css'

const LABELS: Record<JobStatus, string> = {
  queued:    'Queued',
  running:   'Running',
  complete:  'Complete',
  failed:    'Failed',
  cancelled: 'Cancelled',
}

interface Props { status: JobStatus }

export function JobBadge({ status }: Props) {
  return <span className={`${styles.badge} ${styles[status]}`}>{LABELS[status]}</span>
}
