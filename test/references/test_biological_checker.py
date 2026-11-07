# SPDX-FileCopyrightText: 2026 Jonathan D.A. Jewell (hyperpolymath) <j.d.a.jewell@open.ac.uk>
# SPDX-License-Identifier: MPL-2.0
"""Negative controls for the independent CSV checker (stdlib unittest)."""
from copy import deepcopy
from decimal import Decimal, localcontext
import hashlib
import json
from pathlib import Path
import tempfile
import unittest

from check_biological_reference import check


class IndependentCheckerControls(unittest.TestCase):
    def test_accepts_known_matrix_and_rejects_corruption(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            csv_path, result_path = root / 'counts.csv', root / 'result.json'
            raw = b'feature,A,empty\nx,1,0\ny,1,0\n'
            csv_path.write_bytes(raw)
            provenance_path = root / 'provenance.json'
            provenance_bytes = b'{"source":"known test fixture"}'
            provenance_path.write_bytes(provenance_bytes)
            with localcontext() as context:
                context.prec = 100
                ln2 = str(Decimal(2).ln())
            metric = lambda value: {'kind': 'approximate', 'precision_bits': 256, 'value': value}
            result = {
                'csv_sha256': hashlib.sha256(raw).hexdigest(),
                'provenance_sha256': hashlib.sha256(provenance_bytes).hexdigest(),
                'provenance': json.loads(provenance_bytes),
                'exact': {'mode': 'exact_counts', 'input_status': 'exact', 'samples': [
                    {'label': 'A', 'features': ['x', 'y'], 'total': 2, 'counts': [1, 1],
                     'proportions': ['1/2', '1/2'], 'status': 'success'},
                    {'label': 'empty', 'features': ['x', 'y'], 'total': 0, 'counts': [0, 0],
                     'proportions': [None, None], 'status': 'undefined_zero_total'},
                ]},
                'high_precision': {'achieved_precision_bits': 256, 'samples': [
                    {'label': 'A', 'total': 2, 'richness': 2, 'shannon': metric(ln2),
                     'gini_simpson': metric('0.5'), 'status': 'success'},
                    {'label': 'empty', 'total': 0, 'richness': 0, 'shannon': None,
                     'gini_simpson': None, 'status': 'undefined_zero_total'},
                ]},
            }
            result_path.write_text(json.dumps(result))
            self.assertEqual(check(csv_path, result_path, provenance_path)['exact_proportions_checked'], 4)
            corruptions = [
                lambda r: r['exact']['samples'][0].update(total=3),
                lambda r: r['exact']['samples'][0].update(total=2.5),
                lambda r: r['exact']['samples'][0].update(counts=[True, 1]),
                lambda r: r['exact']['samples'][0].update(proportions=['1/3', '2/3']),
                lambda r: r['exact']['samples'][0].update(features=['y', 'x']),
                lambda r: r['high_precision']['samples'][0].update(shannon=metric('0')),
                lambda r: r['high_precision']['samples'][0].update(gini_simpson=metric(0.5)),
                lambda r: r['high_precision']['samples'][0].update(status='undefined_zero_total'),
                lambda r: r['high_precision']['samples'][1].update(shannon=metric('0')),
                lambda r: r.update(csv_sha256='0' * 64),
                lambda r: r.update(provenance_sha256='0' * 64),
                lambda r: r.update(provenance={'source': 'different fixture'}),
            ]
            for corrupt in corruptions:
                candidate = deepcopy(result)
                corrupt(candidate)
                result_path.write_text(json.dumps(candidate))
                with self.assertRaises(AssertionError):
                    check(csv_path, result_path, provenance_path)


if __name__ == '__main__':
    unittest.main()
