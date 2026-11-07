#!/usr/bin/env python3
# SPDX-FileCopyrightText: 2026 Jonathan D.A. Jewell (hyperpolymath) <j.d.a.jewell@open.ac.uk>
# SPDX-License-Identifier: MPL-2.0
"""Check a real counts CSV independently with Python int/Fraction/Decimal.

The CSV, not Julia's reconstructed counts, is the reference. No Julia imports,
Float64 conversions or package output are used to calculate expected values.
Usage: check_biological_reference.py COUNTS.csv RESULT.json CHECK.json [PROVENANCE.json]
"""
import csv
from decimal import Decimal, localcontext
from fractions import Fraction
import hashlib
import json
from pathlib import Path
import re
import sys


def exact_integer(value):
    if type(value) is int:
        assert 0 <= value <= 2**53 - 1, 'unsafe JSON integer'
        return value
    assert isinstance(value, str) and re.fullmatch(r'(?:0|[1-9][0-9]*)', value), 'inexact integer representation'
    return int(value)


def check(csv_path, result_path, provenance_path=None):
    if not __debug__:
        raise RuntimeError('validation must run without Python -O')
    raw = Path(csv_path).read_bytes()
    with Path(csv_path).open(newline='', encoding='utf-8-sig') as stream:
        rows = list(csv.reader(stream))
    assert len(rows) > 1 and len(rows[0]) > 1, 'empty counts table'
    labels = rows[0][1:]
    features = [row[0] for row in rows[1:]]
    assert len(set(labels)) == len(labels) and all(s.strip() for s in labels), 'invalid samples'
    assert len(set(features)) == len(features) and all(s.strip() for s in features), 'invalid features'
    assert all(len(row) == len(rows[0]) for row in rows), 'ragged CSV'
    assert all(re.fullmatch(r'[+]?[0-9]+', cell) for row in rows[1:] for cell in row[1:]), 'noninteger count'
    matrix = [[int(cell) for cell in row[1:]] for row in rows[1:]]
    result_bytes = Path(result_path).read_bytes()
    result = json.loads(result_bytes)
    provenance_hash = None
    if provenance_path is not None:
        provenance_bytes = Path(provenance_path).read_bytes()
        provenance_hash = hashlib.sha256(provenance_bytes).hexdigest()
        assert result['provenance_sha256'] == provenance_hash, 'provenance hash mismatch'
        assert result['provenance'] == json.loads(provenance_bytes), 'provenance content mismatch'
    assert result['csv_sha256'] == hashlib.sha256(raw).hexdigest(), 'CSV hash mismatch'
    assert result['exact']['input_status'] == 'exact'
    assert result['exact']['mode'] == 'exact_counts'
    assert result['high_precision']['achieved_precision_bits'] == 256
    assert len(result['exact']['samples']) == len(labels)
    assert len(result['high_precision']['samples']) == len(labels)
    errors = []
    checked = 0
    # 100 decimal digits in the independent oracle; a fixed 1e-65 absolute
    # tolerance for the 256-bit producer. This is a test threshold, not a proof.
    tolerance = Decimal('1e-65')
    for j, label in enumerate(labels):
        column = [row[j] for row in matrix]
        total = sum(column)
        sample = result['exact']['samples'][j]
        high = result['high_precision']['samples'][j]
        assert sample['label'] == high['label'] == label, 'sample alignment mismatch'
        assert sample['features'] == features, 'feature alignment mismatch'
        assert exact_integer(sample['total']) == exact_integer(high['total']) == total, 'total mismatch'
        assert [exact_integer(c) for c in sample['counts']] == column, 'count mismatch'
        assert len(sample['proportions']) == len(column), 'proportion length mismatch'
        expected = [None if total == 0 else Fraction(c, total) for c in column]
        for observed, fraction in zip(sample['proportions'], expected):
            assert observed == (None if fraction is None else f'{fraction.numerator}/{fraction.denominator}'), 'proportion mismatch'
            checked += 1
        assert exact_integer(high['richness']) == sum(c > 0 for c in column), 'richness mismatch'
        if total == 0:
            assert sample['status'] == high['status'] == 'undefined_zero_total'
            assert high['shannon'] is None and high['gini_simpson'] is None
            continue
        assert sum(expected, Fraction()) == 1, 'independent proportions do not sum to one'
        assert sample['status'] == high['status'] == 'success', 'incorrect success state'
        with localcontext() as context:
            context.prec = 100
            proportions = [Decimal(c) / Decimal(total) for c in column if c]
            shannon = -sum((p * p.ln() for p in proportions), Decimal())
            gini = Fraction(total * total - sum(c * c for c in column), total * total)
            simpson = Decimal(gini.numerator) / Decimal(gini.denominator)
            for name, reference in [('shannon', shannon), ('gini_simpson', simpson)]:
                metric = high[name]
                assert metric['kind'] == 'approximate' and metric['precision_bits'] == 256
                assert isinstance(metric['value'], str), 'metric must preserve its decimal text'
                observed = Decimal(metric['value'])
                assert observed.is_finite(), f'nonfinite {name}'
                error = abs(observed - reference)
                assert error <= tolerance, f'{label} {name} error {error} > {tolerance}'
                errors.append(str(error))
    return {
        'schema': 'metamanifold.independent-descriptive-check.v1', 'status': 'passed',
        'csv_sha256': hashlib.sha256(raw).hexdigest(),
        'result_sha256': hashlib.sha256(result_bytes).hexdigest(),
        'provenance_sha256': provenance_hash,
        'samples': labels, 'features': len(features), 'exact_proportions_checked': checked,
        'exact_tolerance': '0', 'metric_absolute_tolerance': str(tolerance),
        'maximum_metric_error': str(max(map(Decimal, errors), default=Decimal())),
        'python': sys.version, 'oracle': 'stdlib int, fractions.Fraction, decimal.Decimal (100 digits)',
        'claim': 'agreement on this supplied matrix only; no inference, taxonomy validation or scientific sign-off',
    }


if __name__ == '__main__':
    if len(sys.argv) not in (4, 5):
        raise SystemExit('usage: check_biological_reference.py COUNTS.csv RESULT.json CHECK.json [PROVENANCE.json]')
    # No success report is written until all assertions have passed. Python -O
    # would remove those assertions, so explicitly refuse that execution mode.
    if not __debug__:
        raise SystemExit('validation must run without Python -O')
    report = check(sys.argv[1], sys.argv[2], sys.argv[4] if len(sys.argv) == 5 else None)
    Path(sys.argv[3]).write_text(json.dumps(report, indent=2) + '\n', encoding='utf-8')
    print(f"PASS: {report['exact_proportions_checked']} exact proportions; {len(report['samples'])} samples")
