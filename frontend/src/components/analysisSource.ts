// SPDX-FileCopyrightText: 2026 Joshua Benjamin Jewell
// SPDX-License-Identifier: MPL-2.0

import { createContext, useContext } from 'react'
import type { AnnotationSource } from '../api/types'

/** The classifier an analysis tab reads taxonomy and categories from. */
export const AnalysisSourceContext = createContext<AnnotationSource | null>(null)

export const useAnalysisSource = () => useContext(AnalysisSourceContext)
