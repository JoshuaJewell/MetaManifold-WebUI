// © 2026 Joshua Benjamin Jewell. All rights reserved.
// Licensed under the GNU Affero General Public License version 3 (AGPLv3).
import { useState } from 'react'
import { api } from '../api/client'
import type { ReportKind } from '../api/types'
import { NameDialog } from './NameDialog'
import { useToast } from './Toast'

/** Adds the item `make` produces to the study's report, under a caption the user gives. */
export function AddToReport({ study, kind, defaultTitle, make, className = 'btn btn-sm' }: {
  study: string
  kind: ReportKind
  defaultTitle: string
  make: () => Promise<{ blob: Blob; ext: string }>
  className?: string
}) {
  const toast = useToast()
  const [open, setOpen] = useState(false)
  return (
    <>
      <button type="button" className={className} onClick={() => setOpen(true)}>Add to report</button>
      {open && (
        <NameDialog title="Add to report" initialValue={defaultTitle} placeholder="Caption"
          onClose={() => setOpen(false)}
          onConfirm={async title => {
            const { blob, ext } = await make()
            await api.report.add(study, kind, title, ext, blob)
            setOpen(false)
            toast.success(`Added "${title}" to the report`)
          }} />
      )}
    </>
  )
}
