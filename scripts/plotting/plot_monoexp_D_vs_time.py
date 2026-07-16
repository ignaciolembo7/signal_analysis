from __future__ import annotations

import repo_bootstrap  # noqa: F401

import argparse
from pathlib import Path

from data_processing.io import write_xlsx_csv_outputs
from monoexp_fitting.plot_monoexp_D_vs_time import (
    aggregate_monoexp_by_x,
    load_monoexp_fit_measurements,
    plot_compare_N_within_sheet,
    plot_compare_direction,
    plot_compare_roi,
)


def _clear_compare_n_pngs(out_dir: Path, xcol: str) -> None:
    compare_dir = out_dir / xcol / "compare_N"
    if not compare_dir.exists():
        return
    for png in compare_dir.glob("*.png"):
        png.unlink()


def main() -> None:
    ap = argparse.ArgumentParser(
        description="Build monoexp D vs td_ms and Delta_app_ms plots from fit_params*.parquet."
    )
    ap.add_argument("--fits-root", required=True, help="Root folder with monoexp fit_params*.parquet files.")
    ap.add_argument("--out-dir", required=True, help="Output folder for combined tables and PNGs.")
    ap.add_argument("--pattern", default="**/fit_params*.parquet", help="Relative glob inside fits-root.")
    ap.add_argument("--subjs", nargs="+", default=None, help="Subjects/phantoms to include.")
    ap.add_argument("--rois", nargs="+", default=None, help="ROIs to include.")
    ap.add_argument("--dirs", nargs="+", default=None, help="Directions to include.")
    ap.add_argument("--Ns", "--ns", nargs="+", type=float, default=None, help="N values to include.")
    ap.add_argument("--stat", default="avg", help="Stat to keep from the monoexp fit.")
    ap.add_argument("--ycol", default="value_norm", help="ycol to keep from the monoexp fit.")
    args = ap.parse_args()

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    raw = load_monoexp_fit_measurements(
        args.fits_root,
        pattern=args.pattern,
        subjs=args.subjs,
        rois=args.rois,
        directions=args.dirs,
        Ns=args.Ns,
        stat=args.stat,
        ycol=args.ycol,
    )

    write_xlsx_csv_outputs(raw, out_dir / "monoexp_D.raw.xlsx", csv_path=out_dir / "monoexp_D.raw.csv")

    for xcol in ("td_ms", "Delta_app_ms"):
        avg = aggregate_monoexp_by_x(raw, xcol=xcol)
        write_xlsx_csv_outputs(
            avg,
            out_dir / f"monoexp_D_vs_{xcol}.combined.xlsx",
            csv_path=out_dir / f"monoexp_D_vs_{xcol}.combined.csv",
        )
        plot_compare_roi(avg, xcol=xcol, out_dir=out_dir)
        plot_compare_direction(avg, xcol=xcol, out_dir=out_dir)
        _clear_compare_n_pngs(out_dir, xcol)
        plot_compare_N_within_sheet(avg, xcol=xcol, out_dir=out_dir)

    print(f"[OK] Monoexp D summary plots + tables in: {out_dir}")


if __name__ == "__main__":
    main()
