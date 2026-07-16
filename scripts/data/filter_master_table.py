from __future__ import annotations

import argparse
from pathlib import Path

import repo_bootstrap  # noqa: F401

import numpy as np
import pandas as pd

from data_processing.io import write_table_outputs
from data_processing.master_table import load_master_table


def parse_last_points_spec(spec: str) -> dict[tuple[float, float | None], int | None]:
    """Parse a last-points spec string into a mapping of (td_ms, N) -> n_points.

    Format: comma- or space-separated tokens of the form ``TD=POINTS`` or
    ``TD:N=POINTS``.  When only ``TD`` is given the rule applies to all N
    values at that td_ms.  POINTS may be ``ALL`` to keep every point.

    Returns a dict keyed by ``(td_ms, N_or_None)``.
    """
    out: dict[tuple[float, float | None], int | None] = {}
    for token in str(spec).replace(",", " ").split():
        if not token:
            continue
        if "=" not in token:
            raise ValueError(f"Expected TD=POINTS or TD:N=POINTS, received {token!r}.")
        lhs, raw_points = token.split("=", 1)
        if ":" in lhs:
            raw_td, raw_n = lhs.split(":", 1)
            td_ms = float(raw_td)
            n_val: float | None = float(raw_n)
        else:
            td_ms = float(lhs)
            n_val = None
        points_text = raw_points.strip()
        if points_text.upper() == "ALL":
            out[(td_ms, n_val)] = None
            continue
        points = int(points_text)
        if points <= 0:
            raise ValueError(f"Last-points limit for td_ms={td_ms:g} must be > 0.")
        out[(td_ms, n_val)] = points
    return out


def filter_last_points(df: pd.DataFrame, spec: str) -> pd.DataFrame:
    """Keep only the last *n* b_step values for each (td_ms[, N]) group.

    Rules with a N qualifier (``TD:N=POINTS``) are matched on both td_ms and N.
    Rules without a N qualifier (``TD=POINTS``) match all N values at that td_ms.
    More-specific ``TD:N`` rules take precedence over broader ``TD`` rules when
    both apply to the same row.
    """
    limits = parse_last_points_spec(spec)
    if not limits:
        return df.copy().reset_index(drop=True)
    if "td_ms" not in df.columns:
        raise KeyError("Cannot filter last points because column 'td_ms' is missing.")
    if "b_step" not in df.columns:
        raise KeyError("Cannot filter last points because column 'b_step' is missing.")

    has_n_col = "N" in df.columns
    out = df.copy()
    keep = pd.Series(True, index=out.index)
    td_values = pd.to_numeric(out["td_ms"], errors="coerce")
    point_values = pd.to_numeric(out["b_step"], errors="coerce")
    n_values = pd.to_numeric(out["N"], errors="coerce") if has_n_col else None

    # Apply rules from least specific (td-only) to most specific (td+N) so that
    # the more-specific rule wins on any overlap.
    td_only = {k: v for k, v in limits.items() if k[1] is None}
    td_and_n = {k: v for k, v in limits.items() if k[1] is not None}

    for (td_ms, n_val), point_limit in {**td_only, **td_and_n}.items():
        if point_limit is None:
            continue
        td_mask = np.isclose(td_values, float(td_ms), atol=1e-6, equal_nan=False)
        if n_val is not None:
            if not has_n_col:
                raise KeyError("Cannot filter by N because column 'N' is missing.")
            n_mask = np.isclose(n_values, float(n_val), atol=1e-6, equal_nan=False)
            td_mask = td_mask & n_mask
        if not td_mask.any():
            continue
        point_order = np.sort(point_values[td_mask].dropna().unique())
        allowed = set(point_order[-int(point_limit):])
        keep.loc[td_mask] = point_values[td_mask].isin(allowed) | point_values[td_mask].isna()

    return out.loc[keep].reset_index(drop=True)


def build_parser() -> argparse.ArgumentParser:
    ap = argparse.ArgumentParser(description="Write a filtered master table for downstream analysis steps.")
    ap.add_argument("master_parquet", type=Path, help="Input master.long.parquet.")
    ap.add_argument("--out-parquet", type=Path, required=True, help="Output filtered master parquet.")
    ap.add_argument(
        "--last-points-by-td-n",
        required=True,
        help=(
            "Comma- or space-separated rules of the form TD=POINTS or TD:N=POINTS. "
            "Each rule keeps the last POINTS b_step values for that td_ms (and optionally N). "
            "Example: 120:8=6,120:4=4,210=8. Unlisted td_ms values keep all points."
        ),
    )
    return ap


def main() -> None:
    args = build_parser().parse_args()
    df = load_master_table(args.master_parquet)
    out = filter_last_points(df, args.last_points_by_td_n)
    write_table_outputs(out, args.out_parquet)
    print(f"Saved filtered master: {args.out_parquet}")
    print(f"Rows: {len(df)} -> {len(out)}")


if __name__ == "__main__":
    main()
