# SPDX-FileCopyrightText: 2026 Jonathan D.A. Jewell (hyperpolymath) <j.d.a.jewell@open.ac.uk>
# SPDX-License-Identifier: MPL-2.0

#!/usr/bin/env python3
"""Reproduce the independent high-precision known answers used by Julia tests.

Python's decimal module is deliberately independent of Julia/MPFR. No result
produced by the implementation is read here. With --check (the default when a
path is supplied), verify the constants embedded in test_high_precision.jl.
"""
from __future__ import annotations

import argparse
from decimal import Decimal, localcontext
import json
from pathlib import Path
import re
import sys

PRECISION = 460


def references(precision: int = PRECISION) -> dict[str, str]:
    with localcontext() as context:
        context.prec = precision
        counts = (1, 2, 3, 4)
        total = sum(counts)
        h = sum(
            -(Decimal(count) / total) * (Decimal(count) / total).ln()
            for count in counts
        )
        return {
            "REF_LN2": str(Decimal(2).ln()),
            "REF_LN3": str(Decimal(3).ln()),
            "REF_TWO_THIRDS": str(Decimal(2) / Decimal(3)),
            "REF_H_1234": str(h),
        }


def check(path: Path, expected: dict[str, str]) -> list[str]:
    text = path.read_text(encoding="utf-8")
    failures = []
    for name in expected:
        match = re.search(rf'^const {name} = "([^"]+)"$', text, re.MULTILINE)
        if not match:
            failures.append(f"{name}: constant not found")
            continue
        embedded = match.group(1)
        # REF_LN2 deliberately carries enough digits for the 1024-bit test;
        # the other references currently use 220. Verify each value at the
        # precision actually embedded instead of rejecting extra correct digits.
        precision = len(Decimal(embedded).as_tuple().digits)
        if embedded != references(precision)[name]:
            failures.append(f"{name}: embedded reference differs from Python decimal")
    return failures


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("julia_test", nargs="?", type=Path)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    expected = references()
    if args.check or args.julia_test:
        if not args.julia_test:
            parser.error("--check requires the Julia test path")
        failures = check(args.julia_test, expected)
        if failures:
            print("\n".join(failures), file=sys.stderr)
            return 1
        print("Independent decimal references agree at their embedded precisions.")
        return 0
    print(json.dumps({"decimal_precision": PRECISION, "references": expected}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
