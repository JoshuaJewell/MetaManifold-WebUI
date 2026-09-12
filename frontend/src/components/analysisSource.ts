// © 2026 Joshua Benjamin Jewell. All rights reserved.
// Licensed under the GNU Affero General Public License version 3 (AGPLv3).
import { createContext, useContext } from 'react'
import type { AnnotationSource } from '../api/types'

/** The classifier an analysis tab reads taxonomy and categories from. */
export const AnalysisSourceContext = createContext<AnnotationSource | null>(null)

export const useAnalysisSource = () => useContext(AnalysisSourceContext)
