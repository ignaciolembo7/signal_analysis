#!/usr/bin/env python3
"""Expand sparse phantom PGSE result tables to the directional table layout.

The affected acquisitions contain two b0 rows, three nonzero b-value rows that
were acquired once, and one final b-value step acquired in every direction.
The signal-analysis table contract requires every nonzero b-value step to have
one row per direction. This utility repeats each singleton row to that width.
"""

from __future__ import annotations

import argparse
import os
import shutil
from pathlib import Path

import numpy as np
import pandas as pd


SHEETS = ("avg", "std", "med", "mad", "mode")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("results_root", type=Path, help="Results directory containing the sparse PGSE workbooks.")
    parser.add_argument("--glob", default="*acq-hz000*_results.xlsx", help="Workbook glob inside results_root.")
    parser.add_argument("--ndirs", type=int, default=6, help="Required number of directional rows per b-value step.")
    parser.add_argument("--b0-rows", type=int, default=2, help="Number of leading b0 rows to preserve.")
    parser.add_argument("--singleton-steps", type=int, default=3, help="Number of singleton b-value rows to expand.")
    parser.add_argument("--dry-run", action="store_true", help="Validate and report without writing files.")
    return parser.parse_args()


def read_and_validate(path: Path, *, b0_rows: int, singleton_steps: int, ndirs: int) -> dict[str, pd.DataFrame]:
    book = pd.ExcelFile(path)
    if tuple(book.sheet_names) != SHEETS:
        raise ValueError(f"{path}: expected sheets {list(SHEETS)}, got {book.sheet_names}")

    tables = {sheet: pd.read_excel(path, sheet_name=sheet) for sheet in SHEETS}
    reference = tables[SHEETS[0]]
    expected_sparse_rows = b0_rows + singleton_steps + ndirs
    expected_expanded_rows = b0_rows + (singleton_steps + 1) * ndirs

    if len(reference) == expected_expanded_rows:
        return tables
    if len(reference) != expected_sparse_rows:
        raise ValueError(
            f"{path}: expected {expected_sparse_rows} sparse rows or "
            f"{expected_expanded_rows} expanded rows, got {len(reference)}"
        )
    if reference.columns.empty or str(reference.columns[0]).strip().lower() != "bvalues":
        raise ValueError(f"{path}: first column must be 'bvalues'")

    reference_bvalues = pd.to_numeric(reference.iloc[:, 0], errors="raise").to_numpy(float)
    if not np.allclose(reference_bvalues[:b0_rows], 0.0, rtol=0.0, atol=1e-9):
        raise ValueError(f"{path}: the first {b0_rows} rows are not all b0 rows")
    if np.any(np.isclose(reference_bvalues[b0_rows:], 0.0, rtol=0.0, atol=1e-9)):
        raise ValueError(f"{path}: unexpected b0 row after the leading b0 block")

    for sheet, table in tables.items():
        if table.columns.tolist() != reference.columns.tolist():
            raise ValueError(f"{path}/{sheet}: columns differ from the avg sheet")
        if len(table) != len(reference):
            raise ValueError(f"{path}/{sheet}: row count differs from the avg sheet")
        sheet_bvalues = pd.to_numeric(table.iloc[:, 0], errors="raise").to_numpy(float)
        if not np.allclose(sheet_bvalues, reference_bvalues, rtol=0.0, atol=1e-9):
            raise ValueError(f"{path}/{sheet}: b-values differ from the avg sheet")

    return tables


def expand_table(table: pd.DataFrame, *, b0_rows: int, singleton_steps: int, ndirs: int) -> pd.DataFrame:
    expected_expanded_rows = b0_rows + (singleton_steps + 1) * ndirs
    if len(table) == expected_expanded_rows:
        return table.copy()

    blocks = [table.iloc[:b0_rows]]
    for row_index in range(b0_rows, b0_rows + singleton_steps):
        blocks.append(pd.concat([table.iloc[[row_index]]] * ndirs, ignore_index=True))
    blocks.append(table.iloc[b0_rows + singleton_steps :])
    return pd.concat(blocks, ignore_index=True)


def main() -> None:
    args = parse_args()
    if args.ndirs <= 0 or args.b0_rows < 0 or args.singleton_steps <= 0:
        raise SystemExit("[STOP] --ndirs and --singleton-steps must be positive; --b0-rows must be non-negative.")

    results_root = args.results_root.resolve()
    if not results_root.is_dir():
        raise SystemExit(f"[STOP] Results directory not found: {results_root}")

    paths = sorted(results_root.glob(args.glob))
    if not paths:
        raise SystemExit(f"[STOP] No workbooks matched {args.glob!r} in {results_root}")

    backup_root = results_root / "original_sparse_pgse"
    changed = 0
    skipped = 0

    for path in paths:
        tables = read_and_validate(
            path,
            b0_rows=args.b0_rows,
            singleton_steps=args.singleton_steps,
            ndirs=args.ndirs,
        )
        old_rows = len(tables[SHEETS[0]])
        expanded = {
            sheet: expand_table(
                table,
                b0_rows=args.b0_rows,
                singleton_steps=args.singleton_steps,
                ndirs=args.ndirs,
            )
            for sheet, table in tables.items()
        }
        new_rows = len(expanded[SHEETS[0]])

        if old_rows == new_rows:
            print(f"[SKIP] Already expanded ({new_rows} rows): {path.name}")
            skipped += 1
            continue

        print(f"[EXPAND] {path.name}: {old_rows} -> {new_rows} rows per sheet")
        if args.dry_run:
            changed += 1
            continue

        backup_root.mkdir(parents=True, exist_ok=True)
        backup_path = backup_root / path.name
        if not backup_path.exists():
            shutil.copy2(path, backup_path)

        temp_path = path.with_name(f".{path.stem}.expanding.xlsx")
        try:
            with pd.ExcelWriter(temp_path, engine="openpyxl") as writer:
                for sheet in SHEETS:
                    expanded[sheet].to_excel(writer, sheet_name=sheet, index=False)
            os.replace(temp_path, path)
        finally:
            if temp_path.exists():
                temp_path.unlink()
        changed += 1

    print(f"Done. Matched={len(paths)} expanded={changed} already_expanded={skipped} dry_run={args.dry_run}")
    if not args.dry_run and changed:
        print(f"Original workbooks: {backup_root}")


if __name__ == "__main__":
    main()
