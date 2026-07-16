from __future__ import annotations
from pathlib import Path
import numpy as np
import pandas as pd

def to_long(
    stats: dict[str, pd.DataFrame],
    *,
    ndirs: int,
    nbvals: int,
    bcol: str = "bvalues",
    out_col: str = "bvalue",
    b0_reps: int = 2,
    source_file: str | None = None,
) -> pd.DataFrame:
    """
    Convert interleaved tables (stride=ndirs) to long/tidy format.
    Replicate one averaged b0 row per direction, matching the notebook convention.
    """
    if not stats:
        raise ValueError("stats is empty.")

    any_df = next(iter(stats.values()))
    if bcol not in any_df.columns:
        for alt in ["bvalue", "bval", "b"]:
            if alt in any_df.columns:
                bcol = alt
                break
    if bcol not in any_df.columns:
        raise ValueError(f"Could not find column {bcol!r}. Columns: {list(any_df.columns)}")

    rois = [c for c in any_df.columns if c != bcol]

    rows_expected = ndirs * nbvals
    long_parts: list[pd.DataFrame] = []

    for stat, df in stats.items():
        # Average b0 per ROI.
        b0 = df.loc[: b0_reps - 1, rois].mean(axis=0)

        # Data rows without the repeated b0 measurements.
        data = df.loc[b0_reps:, [bcol] + rois].reset_index(drop=True)
        if len(data) != rows_expected:
            raise ValueError(
                f"[{stat}] Expected {rows_expected} rows (ndirs*nbvals), got {len(data)}. "
                "This usually means ndirs/nbvals do not match the file."
            )

        idx = np.arange(len(data))
        direction = (idx % ndirs) + 1          # 1..ndirs
        b_step = (idx // ndirs) + 1            # 1..nbvals

        data = data.assign(direction=direction, b_step=b_step)

        # Long table without b0 rows.
        long = data.melt(
            id_vars=[bcol, "direction", "b_step"],
            value_vars=rois,
            var_name="roi",
            value_name="value",
        )
        long = long.rename(columns={bcol: out_col})
        long["stat"] = stat

        # Duplicate b0 for each direction (b_step=0).
        b0_long = pd.DataFrame({
            out_col: 0.0,
            "direction": np.repeat(np.arange(1, ndirs + 1), len(rois)),
            "b_step": 0,
            "roi": np.tile(np.array(rois), ndirs),
            "value": np.tile(b0.to_numpy(), ndirs),
            "stat": stat,
        })

        long_parts.append(pd.concat([b0_long, long], ignore_index=True))

    out = pd.concat(long_parts, ignore_index=True)
    out["source_file"] = source_file if source_file is not None else ""

    # Final order: ROI, direction, then b_step within each curve (0 first).
    sort_cols = [c for c in ["stat", "roi", "direction", "b_step"] if c in out.columns]
    out = out.sort_values(sort_cols, kind="stable").reset_index(drop=True)

    return out
