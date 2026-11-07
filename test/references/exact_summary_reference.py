# SPDX-FileCopyrightText: 2026 Jonathan D.A. Jewell (hyperpolymath) <j.d.a.jewell@open.ac.uk>
# SPDX-License-Identifier: MPL-2.0

#!/usr/bin/env python3
"""Independent exact-summary reference using Python fractions.Fraction.

Usage: python3 exact_summary_reference.py INPUT.json
Input/output counts are decimal strings; zero-total proportions are JSON null.
This program never imports or reads output from the Julia implementation.
"""
from __future__ import annotations

from collections import OrderedDict
from fractions import Fraction
import json
from pathlib import Path
import sys


def proportion_text(count: int, total: int) -> str | None:
    if total == 0:
        return None
    value = Fraction(count, total)
    return f"{value.numerator}/{value.denominator}"


def summarise(label: str, columns: list[list[int]], features: list[str]) -> dict:
    counts = [sum(column[i] for column in columns) for i in range(len(features))]
    total = sum(counts)
    return {
        "label": label,
        "total": str(total),
        "counts": [str(count) for count in counts],
        "proportions": [proportion_text(count, total) for count in counts],
    }


def main() -> int:
    if len(sys.argv) != 2:
        raise SystemExit("usage: exact_summary_reference.py INPUT.json")
    payload = json.loads(Path(sys.argv[1]).read_text(encoding="utf-8"))
    labels = payload["sample_labels"]
    features = payload["feature_labels"]
    rows = payload["counts"]
    assert len(rows) == len(features)
    columns = [[int(row[j]) for row in rows] for j in range(len(labels))]
    result = {
        "samples": [summarise(label, [column], features)
                    for label, column in zip(labels, columns)],
        "groups": [],
    }
    groups = payload.get("groups")
    if groups is not None:
        assert len(groups) == len(labels)
        by_group = OrderedDict()
        for label, column in zip(groups, columns):
            by_group.setdefault(label, []).append(column)
        result["groups"] = [summarise(label, members, features)
                            for label, members in by_group.items()]
    print(json.dumps(result, separators=(",", ":")))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
