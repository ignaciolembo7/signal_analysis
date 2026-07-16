from __future__ import annotations

import argparse
from pathlib import Path

import repo_bootstrap  # noqa: F401

from data_processing.master_table import load_master_table


def build_parser() -> argparse.ArgumentParser:
    ap = argparse.ArgumentParser(
        description=(
            "Export the canonical master.long.parquet table to an inspection "
            "workbook. This script is intentionally separate from the pipeline "
            "so append steps never rewrite a large Excel file."
        )
    )
    ap.add_argument("master_parquet", type=Path, help="Input master.long.parquet.")
    ap.add_argument(
        "--out-xlsx",
        type=Path,
        default=None,
        help="Output .xlsx path. Default: input path with .xlsx suffix.",
    )
    ap.add_argument(
        "--row-kind",
        action="append",
        default=[],
        help="Optional row_kind filter. Can be repeated, e.g. --row-kind signal --row-kind signal_rotated.",
    )
    ap.add_argument(
        "--head",
        type=int,
        default=None,
        help="Write only the first N rows after filtering, useful for quick checks.",
    )
    return ap


def main() -> None:
    args = build_parser().parse_args()
    out_xlsx = args.out_xlsx or args.master_parquet.with_suffix(".xlsx")

    df = load_master_table(args.master_parquet)
    if args.row_kind:
        wanted = {str(v) for v in args.row_kind}
        df = df[df["row_kind"].astype(str).isin(wanted)].copy()
    if args.head is not None:
        df = df.head(args.head).copy()

    out_xlsx.parent.mkdir(parents=True, exist_ok=True)
    df.to_excel(out_xlsx, index=False)
    print(f"Saved: {out_xlsx}")
    print(f"Rows: {len(df)}")
    print(f"Columns: {len(df.columns)}")


if __name__ == "__main__":
    main()
