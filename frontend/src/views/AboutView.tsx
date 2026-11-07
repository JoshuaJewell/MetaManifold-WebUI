// SPDX-FileCopyrightText: 2026 Joshua Benjamin Jewell
// SPDX-License-Identifier: MPL-2.0

import { COPYRIGHT, LICENSE_URL, SOURCE_URL, ZENODO_DOI } from '../about'

export function AboutView() {
  return (
    <>
      <div className="page-header">
        <h1>About</h1>
      </div>
      <div className="card" style={{ maxWidth: '70ch', lineHeight: 1.55, fontSize: '.9rem' }}>
        <div className="card-title">MetaManifold</div>
        <p>{COPYRIGHT}.</p>
        <p style={{ marginTop: 8 }}>
          This program is free software: you can redistribute it and/or modify it under the terms of
          the GNU Affero General Public License as published by the Free Software Foundation, either
          version 3 of the License, or (at your option) any later version. Individual files may also
          be available under other licences, as stated in their SPDX headers.
        </p>
        <p style={{ marginTop: 8 }}>
          This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
          without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
          See the <a href={LICENSE_URL} target="_blank" rel="noreferrer">GNU Affero General Public License</a> for
          more details.
        </p>
        <dl style={{ marginTop: 12, display: 'grid', gridTemplateColumns: 'max-content 1fr', gap: '4px 12px' }}>
          <dt>Source and documentation</dt>
          <dd><a href={SOURCE_URL} target="_blank" rel="noreferrer">{SOURCE_URL}</a></dd>
          <dt>Licence</dt>
          <dd><a href={LICENSE_URL} target="_blank" rel="noreferrer">AGPL-3.0-or-later</a></dd>
          <dt>Cite</dt>
          <dd>
            {ZENODO_DOI
              ? <a href={`https://doi.org/${ZENODO_DOI}`} target="_blank" rel="noreferrer">doi:{ZENODO_DOI}</a>
              : <span style={{ color: 'var(--color-muted-fg)' }}>Zenodo DOI to be assigned</span>}
          </dd>
        </dl>
      </div>
    </>
  )
}
