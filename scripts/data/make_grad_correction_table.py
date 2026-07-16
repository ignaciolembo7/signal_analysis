from __future__ import annotations

import repo_bootstrap  # noqa: F401

import argparse
from pathlib import Path

import pandas as pd

from data_processing.io import write_xlsx_sheets
from ogse_fitting.make_grad_correction_table import (
    make_grad_correction_from_manifest,
    _fill_missing_correction_factors_by_avg,
    _update_master_correction_factors,
)


def main() -> None:
    ap = argparse.ArgumentParser(
        description=(
            'Build gradient-correction factors from a syringe (or water) manifest.\n'
            '\n'
            'For each curve listed in the manifest the script fits both the NOGSE free model\n'
            'and a monoexp model to the raw signal in the master parquet, then computes\n'
            '  correction_factor = sqrt(D0_nogse / D0_monoexp)\n'
            'and writes it back to the master table for ALL rois that share the same\n'
            '(subj, sheet, direction, td_ms, N) acquisition parameters.'
        )
    )
    ap.add_argument(
        '--manifest', required=True,
        help='CSV manifest with columns subj,sheet,roi,direction,td_ms,N[,Hz,model]. '
             'Typically bash_template/manifests/<dataset>/grad_correction.csv.',
    )
    ap.add_argument(
        '--master-parquet', type=Path, required=True,
        help='Path to master.long.parquet.',
    )
    ap.add_argument('--out-xlsx', required=True, help='Output .xlsx path.')
    ap.add_argument('--out-csv', default=None, help='Optional .csv output.')
    ap.add_argument(
        '--plot-dir', type=Path, default=None,
        help='Directory to save per-curve comparison plots (NOGSE fit vs monoexp fit).',
    )
    ap.add_argument(
        '--roi', default=None,
        help='Override the manifest roi column for every row, so all curves are '
             'matched against this reference ROI in the master table instead of '
             "whatever the manifest's roi column says (e.g. 'Syringe' for brains, "
             "'Water1' for phantoms). Matching against the master table is "
             'case-insensitive. Default: use each row\'s roi as given in the manifest.',
    )

    ap.add_argument(
        '--stat', default='avg',
        help='Statistic row to use from master (default: avg).',
    )
    ap.add_argument(
        '--row-kind', default='signal_rotated',
        help='row_kind to load from master (default: signal_rotated).',
    )
    ap.add_argument(
        '--gbase', default='g_lin_max',
        help='Gradient column for NOGSE free fit (default: g_lin_max).',
    )
    ap.add_argument(
        '--bbase', default='bvalue_thorsten',
        help='B-value column for monoexp fit (default: bvalue_thorsten).',
    )
    ap.add_argument(
        '--ycol', default='value_norm',
        help='Signal column to fit (default: value_norm). Use value for raw signal.',
    )
    ap.add_argument(
        '--D0-init', type=float, default=2.3e-12,
        help='D0 seed for NOGSE free fit in m2/ms (default: 2.3e-12).',
    )
    ap.add_argument(
        '--tol-ms', type=float, default=1e-3,
        help='Tolerance in ms when matching td_ms values (default: 1e-3).',
    )

    ap.add_argument(
        '--avg-N',
        nargs='*',
        type=int,
        default=None,
        metavar='N',
        help=(
            'Average D0_monoexp across N values (in addition to directions) when '
            'computing the correction factor. '
            'No values: average over ALL N values. '
            'Specific values: average only over those N (e.g. --avg-N 4 8). '
            'Omit entirely (default): no averaging over N, one D0_monoexp per '
            '(subj, sheet, roi, td_ms, N).'
        ),
    )

    ap.add_argument(
        '--no-fill-missing',
        action='store_true',
        help=(
            'Skip the cross-subject fill step. By default, after writing manifest-derived '
            'factors, any signal rows still missing a factor receive the mean factor from '
            'other subjects at the same (direction, td_ms, N).'
        ),
    )

    m0_group = ap.add_mutually_exclusive_group()
    m0_group.add_argument(
        '--fix-M0', type=float, default=1.0,
        help='Fix M0 to this value in both fits (default: 1.0).',
    )
    m0_group.add_argument(
        '--free-M0', nargs='?', const=1.0, type=float, default=None,
        help='Allow M0 to vary with optional seed.',
    )

    args = ap.parse_args()

    M0_vary = args.free_M0 is not None
    M0_value = float(args.free_M0) if args.free_M0 is not None else float(args.fix_M0)

    out = make_grad_correction_from_manifest(
        master_parquet=args.master_parquet,
        manifest=args.manifest,
        stat_keep=args.stat,
        row_kind=args.row_kind,
        gbase=args.gbase,
        bbase=args.bbase,
        ycol=args.ycol,
        M0_vary=M0_vary,
        M0_value=M0_value,
        D0_init=args.D0_init,
        tol_ms=args.tol_ms,
        avg_N=args.avg_N,
        plot_dir=args.plot_dir,
        roi_override=args.roi,
    )

    out_xlsx = Path(args.out_xlsx)
    out_xlsx.parent.mkdir(parents=True, exist_ok=True)
    write_xlsx_sheets({'grad_correction': out}, out_xlsx)
    print('OK:', out_xlsx)

    if args.out_csv:
        out_csv = Path(args.out_csv)
        out_csv.parent.mkdir(parents=True, exist_ok=True)
        out.to_csv(out_csv, index=False)
        print('OK:', out_csv)

    _update_master_correction_factors(args.master_parquet, out, tol_ms=args.tol_ms)

    if not args.no_fill_missing:
        _fill_missing_correction_factors_by_avg(args.master_parquet)

    n_ok = int(out['correction_factor'].notna().sum())
    n_total = len(out)
    print(f'Correction factors computed: {n_ok}/{n_total}')
    if n_ok > 0:
        cf = out.loc[out['correction_factor'].notna(), 'correction_factor']
        print(f'correction_factor range: [{cf.min():.4f}, {cf.max():.4f}], mean={cf.mean():.4f}')


if __name__ == '__main__':
    main()
