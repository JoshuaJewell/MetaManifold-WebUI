// SPDX-FileCopyrightText: 2026 Joshua Benjamin Jewell
// SPDX-License-Identifier: MPL-2.0

import { createContext, useContext } from 'react'

/** Reloads the sidebar's studies, groups and runs after one is created, renamed or deleted. */
export const NavRefreshContext = createContext<() => void>(() => {})

export const useNavRefresh = () => useContext(NavRefreshContext)
