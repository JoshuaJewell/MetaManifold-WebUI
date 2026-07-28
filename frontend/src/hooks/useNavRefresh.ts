// © 2026 Joshua Benjamin Jewell. All rights reserved.
// Licensed under the GNU Affero General Public License version 3 (AGPLv3).
import { createContext, useContext } from 'react'

/** Reloads the sidebar's studies, groups and runs after one is created, renamed or deleted. */
export const NavRefreshContext = createContext<() => void>(() => {})

export const useNavRefresh = () => useContext(NavRefreshContext)
