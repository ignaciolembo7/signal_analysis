from __future__ import annotations

from pathlib import Path
import sys
import tempfile
import unittest

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

REPO_ROOT = Path(__file__).resolve().parents[1]
SRC_ROOT = REPO_ROOT / "src"
if str(SRC_ROOT) not in sys.path:
    sys.path.insert(0, str(SRC_ROOT))

from presentation_analysis.abstract_figures import (  # noqa: E402
    GAMMA_RAD_MS_MT,
    _derived_peak_axes,
    _resolve_example,
    _tc_pseudohuber,
    build_contrast_table,
    build_resampled_contrasts,
    discover_analysis_runs,
    export_all_contrast_lcf_panels,
    export_all_signal_contrast_panels,
    fit_tc_pseudohuber,
    load_alpha_summaries,
    normalize_alpha_to_internal_reference,
    plot_contrast_lcf_grid,
    summarize_contrast_shape,
)
from models.model_fitting import M_ogse_rest_rician  # noqa: E402


class AbstractFigureAnalysisTests(unittest.TestCase):
    def test_derived_peak_axes_use_capiglioni_filter_length(self) -> None:
        td_ms = 143.4
        gradient_mtm = 35.0
        d0_m2_ms = 3.2e-12

        lcf_um, tau_f_ms = _derived_peak_axes(td_ms, gradient_mtm, d0_m2_ms)
        expected_lcf_m = (
            1.5 * d0_m2_ms / (GAMMA_RAD_MS_MT**2 * gradient_mtm**2 * td_ms)
        ) ** 0.25

        self.assertAlmostEqual(lcf_um, expected_lcf_m * 1e6)
        self.assertAlmostEqual(tau_f_ms, expected_lcf_m**2 / d0_m2_ms)

    def test_old_phantom_diffusion_times_use_colored_palette(self) -> None:
        diffusion_times = [75.1, 97.1, 119.1, 142.5, 209.1]
        contrast_rows = []
        fit_rows = []
        tc_rows = []
        for td_ms in diffusion_times:
            fit_rows.append(
                {
                    "variant": "manual",
                    "subj": "P1",
                    "roi": "fiber1",
                    "direction": "long",
                    "td_ms": td_ms,
                    "ok": True,
                }
            )
            tc_rows.append(
                {
                    "variant": "manual",
                    "subj": "P1",
                    "roi": "fiber1",
                    "direction": "long",
                    "td_ms": td_ms,
                    "tc_peak_ms": 5.0 + td_ms / 100.0,
                }
            )
            for lcf_um in (3.5, 5.0, 8.0):
                contrast_rows.append(
                    {
                        "variant": "manual",
                        "subj": "P1",
                        "roi": "fiber1",
                        "direction": "long",
                        "td_ms": td_ms,
                        "lcf_um": lcf_um,
                        "contrast": 0.2,
                    }
                )

        figure = plot_contrast_lcf_grid(
            pd.DataFrame(contrast_rows),
            pd.DataFrame(fit_rows),
            pd.DataFrame(tc_rows),
            pd.DataFrame(
                columns=["variant", "subj", "roi", "direction", "ok", "c_ms", "delta_ms", "alpha_macro"]
            ),
            variants=["manual"],
            subject="P1",
            roi="fiber1",
            directions=["long"],
        )
        curve_colors = [line.get_color() for line in figure.axes[0].lines[: len(diffusion_times)]]
        plt.close(figure)

        self.assertNotIn("#555555", curve_colors)
        self.assertEqual(len(set(curve_colors)), len(diffusion_times))

    def test_alpha_loader_harmonizes_reference_and_prefers_direct_direction(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_directory:
            alpha_path = Path(temporary_directory) / "summary.xlsx"
            pd.DataFrame(
                {
                    "subj": ["P1", "P1"],
                    "sheet": ["S1", "S1"],
                    "roi": ["fiber1", "fiber1"],
                    "direction": ["long", "long"],
                    "direction_kind": ["raw", "derived"],
                    "D0_mean_mm2_s": [0.0010, 0.00115],
                    "D0_std_mm2_s": [0.0001, 0.0001],
                    "alpha_macro": [0.3125, 0.359375],
                    "alpha_macro_error": [0.01, 0.01],
                }
            ).to_excel(alpha_path, index=False)
            runs = pd.DataFrame(
                {
                    "variant": ["manual"],
                    "dataset": ["phantoms"],
                    "alpha_path": [alpha_path],
                }
            )

            loaded = load_alpha_summaries(
                runs,
                reference_d0_m2_ms=2.3e-12,
            )

        self.assertEqual(len(loaded), 1)
        self.assertEqual(loaded.iloc[0]["direction_kind"], "raw")
        self.assertAlmostEqual(float(loaded.iloc[0]["alpha_macro"]), 0.0010 / 0.0023)
        self.assertTrue(bool(loaded.iloc[0]["alpha_reference_harmonized"]))

    def test_internal_alpha_normalizes_against_same_scan_water(self) -> None:
        table = pd.DataFrame(
            {
                "variant": ["manual"] * 3,
                "subj": ["P1"] * 3,
                "sheet": ["S1"] * 3,
                "roi": ["fiber1", "water1", "water2"],
                "direction": ["tra"] * 3,
                "D0_mean_mm2_s": [0.0009, 0.0018, 0.0020],
                "D0_std_mm2_s": [0.0001, 0.0001, 0.0001],
            }
        )

        normalized = normalize_alpha_to_internal_reference(
            table, reference_rois=["water1", "water2"]
        )

        fiber = normalized[normalized["roi"].eq("fiber1")].iloc[0]
        self.assertAlmostEqual(float(fiber["internal_reference_D0_mm2_s"]), 0.0019)
        self.assertAlmostEqual(float(fiber["alpha_internal"]), 0.0009 / 0.0019)

    def test_discovery_finds_future_variant_names(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_directory:
            root = Path(temporary_directory)
            for variant in ["plain", "future-method"]:
                path = (
                    root
                    / "Data-BIDS"
                    / "derivatives"
                    / "signal_analysis"
                    / f"den_gr-topup--{variant}"
                    / "brains"
                    / "ogse_experiments"
                    / "master.long.parquet"
                )
                path.parent.mkdir(parents=True)
                path.touch()

            runs = discover_analysis_runs(root, dataset="brains", dwi_level="den_gr-topup")

        self.assertEqual(set(runs["variant"]), {"plain", "future-method"})

    def test_discovery_supports_phantom_master_layout(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_directory:
            root = Path(temporary_directory)
            path = (
                root
                / "Data-BIDS"
                / "derivatives"
                / "signal_analysis"
                / "raw--manual"
                / "phantoms"
                / "ogse_experiments"
                / "master.long.parquet"
            )
            path.parent.mkdir(parents=True)
            path.touch()

            runs = discover_analysis_runs(root, dataset="phantoms", dwi_level="raw")

        self.assertEqual(runs.iloc[0]["variant"], "manual")
        self.assertEqual(runs.iloc[0]["dataset"], "phantoms")

    def test_build_contrast_pairs_high_and_low_n(self) -> None:
        rows = []
        for stat, values in [("avg", {8: 0.8, 4: 0.6}), ("std", {8: 0.08, 4: 0.06})]:
            for n_value in [8, 4]:
                rows.append(
                    {
                        "variant": "plain",
                        "row_kind": "signal_rotated",
                        "subj": "S1",
                        "sheet": "sheet1",
                        "roi": "PostCC",
                        "direction": "long",
                        "stat": stat,
                        "td_ms": 100.0,
                        "b_step": 1,
                        "N": n_value,
                        "Hz": float(n_value),
                        "sequence": n_value,
                        "source_file": f"n{n_value}.xlsx",
                        "value_norm": values[n_value],
                        "g": float(n_value),
                        "g_thorsten": float(n_value) * 100.0,
                        "grad_correction_factor": 0.5,
                        "delta_ms": 10.0,
                        "Delta_app_ms": 9.0,
                    }
                )

        contrast = build_contrast_table(pd.DataFrame(rows), rois=["PostCC"])

        self.assertEqual(len(contrast), 1)
        self.assertAlmostEqual(float(contrast.iloc[0]["contrast"]), 0.2)
        self.assertAlmostEqual(float(contrast.iloc[0]["contrast_std"]), 0.1)
        self.assertAlmostEqual(float(contrast.iloc[0]["g_high_corr"]), 4.0)
        self.assertAlmostEqual(float(contrast.iloc[0]["g_low_corr"]), 2.0)

    def test_auto_example_prefers_combination_shared_across_variants(self) -> None:
        table = pd.DataFrame(
            {
                "variant": ["manual", "manual", "future-method"],
                "subj": ["S0", "S1", "S1"],
                "roi": ["sample-z", "sample-a", "sample-a"],
                "direction": ["long", "long", "long"],
                "td_ms": [90.0, 120.0, 120.0],
            }
        )

        resolved = _resolve_example(table, subject=None, roi=None, direction=None, td_ms=None)

        self.assertEqual(resolved, ("S1", "sample-a", "long", 120.0))

    def test_resampled_contrast_uses_one_common_corrected_gradient(self) -> None:
        td_ms = 143.4
        d0_m2_ms = 3.2e-12
        c_value = 0.2
        branches = {
            8: (np.linspace(0.0, 70.0, 10), 8.0),
            4: (np.linspace(0.0, 55.0, 10), 5.0),
        }
        rows = []
        for n_value, (gradients, tc_ms) in branches.items():
            signals = M_ogse_rest_rician(
                td_ms,
                gradients,
                n_value,
                td_ms / n_value,
                tc_ms,
                1.0,
                d0_m2_ms,
                c_value,
            )
            for b_step, (gradient, signal) in enumerate(zip(gradients, signals)):
                rows.append(
                    {
                        "variant": "manual",
                        "row_kind": "signal_rotated",
                        "subj": "P1",
                        "sheet": "sheet1",
                        "roi": "fiber1",
                        "direction": "long",
                        "stat": "avg",
                        "td_ms": td_ms,
                        "b_step": b_step,
                        "N": n_value,
                        "source_file": f"n{n_value}.xlsx",
                        "value_norm": signal,
                        "g": gradient,
                        "g_thorsten": gradient * 10.0,
                        "grad_correction_factor": 1.0,
                    }
                )

        masters = pd.DataFrame(rows)
        contrast, summary = build_resampled_contrasts(
            masters,
            d0_m2_ms=d0_m2_ms,
            directions=["long"],
            rois=["fiber1"],
            grid_size=256,
        )

        self.assertTrue(bool(summary.iloc[0]["ok"]))
        self.assertEqual(len(contrast), 256)
        np.testing.assert_allclose(contrast["g_high_corr"], contrast["g_low_corr"])
        self.assertAlmostEqual(float(contrast["g_resampled_corr"].max()), 55.0)
        self.assertAlmostEqual(float(contrast.iloc[0]["contrast"]), 0.0, places=10)
        self.assertAlmostEqual(float(summary.iloc[0]["tc_high_ms"]), 8.0, places=4)
        self.assertAlmostEqual(float(summary.iloc[0]["tc_low_ms"]), 5.0, places=4)
        self.assertAlmostEqual(float(summary.iloc[0]["C"]), c_value, places=4)

        with tempfile.TemporaryDirectory() as temporary_directory:
            output_dir = Path(temporary_directory) / "all_panels"
            legacy_directory = output_dir / "sub-old" / "sheet-old" / "roi-old" / "dir-old"
            legacy_directory.mkdir(parents=True)
            (legacy_directory / "td-90ms.png").touch()
            (legacy_directory / "td-90ms.pdf").touch()
            (legacy_directory / "Thumbs.db").touch()
            (output_dir / "manifest.csv").touch()
            failed_variant = summary.iloc[[0]].copy()
            failed_variant["variant"] = "future-method"
            failed_variant["ok"] = False
            failed_variant["message"] = "Synthetic missing branch"
            batch_summary = pd.concat([summary, failed_variant], ignore_index=True)
            manifest = export_all_signal_contrast_panels(
                masters,
                contrast,
                batch_summary,
                variants=["manual", "future-method"],
                output_dir=output_dir,
                dpi=72,
                progress_every=0,
            )
            figure_path = output_dir / manifest.iloc[0]["figure_path"]
            self.assertEqual(len(manifest), 1)
            self.assertFalse(bool(manifest.iloc[0]["all_selected_variants_successful"]))
            self.assertEqual(manifest.iloc[0]["missing_or_failed_variants"], "future-method")
            self.assertEqual(
                Path(manifest.iloc[0]["figure_path"]).parts,
                (
                    "sub-P1",
                    "sheet-sheet1",
                    "roi-fiber1__dir-long__td-143.4ms.png",
                ),
            )
            self.assertTrue(figure_path.is_file())
            self.assertTrue(figure_path.with_suffix(".pdf").is_file())
            self.assertTrue((output_dir / "manifest.csv").is_file())
            self.assertFalse((output_dir / "sub-old").exists())

            lcf_output_dir = Path(temporary_directory) / "all_lcf_panels"
            stale_directory = lcf_output_dir / "sub-old"
            stale_directory.mkdir(parents=True)
            (stale_directory / "old.png").touch()
            (stale_directory / "old.pdf").touch()
            tc_data = summary.loc[:, ["variant", "subj", "roi", "direction", "td_ms", "tc_peak_ms"]].copy()
            tc_summary = pd.DataFrame(
                columns=["variant", "subj", "roi", "direction", "ok", "c_ms", "delta_ms", "alpha_macro"]
            )
            lcf_manifest = export_all_contrast_lcf_panels(
                contrast,
                batch_summary,
                tc_data,
                tc_summary,
                variants=["manual", "future-method"],
                directions=["long"],
                output_dir=lcf_output_dir,
                dpi=72,
                progress_every=0,
            )
            lcf_figure_path = lcf_output_dir / lcf_manifest.iloc[0]["figure_path"]
            self.assertEqual(len(lcf_manifest), 1)
            self.assertEqual(lcf_manifest.iloc[0]["observed_directions"], "long")
            self.assertEqual(
                Path(lcf_manifest.iloc[0]["figure_path"]).parts,
                ("sub-P1", "roi-fiber1__contrast-vs-filtered-length.png"),
            )
            self.assertTrue(lcf_figure_path.is_file())
            self.assertTrue(lcf_figure_path.with_suffix(".pdf").is_file())
            self.assertTrue((lcf_output_dir / "manifest.csv").is_file())
            self.assertFalse((lcf_output_dir / "sub-old").exists())

    def test_shared_tc_mode_uses_one_physical_correlation_time(self) -> None:
        td_ms = 120.0
        d0_m2_ms = 2.3e-12
        expected_tc_ms = 7.0
        rows = []
        for n_value, maximum_gradient in [(8, 70.0), (4, 55.0)]:
            gradients = np.linspace(0.0, maximum_gradient, 10)
            signals = M_ogse_rest_rician(
                td_ms,
                gradients,
                n_value,
                td_ms / n_value,
                expected_tc_ms,
                1.0,
                d0_m2_ms,
                0.15,
            )
            for b_step, (gradient, signal) in enumerate(zip(gradients, signals)):
                rows.append(
                    {
                        "variant": "manual",
                        "row_kind": "signal_rotated",
                        "subj": "P1",
                        "sheet": "sheet1",
                        "roi": "fiber1",
                        "direction": "tra",
                        "stat": "avg",
                        "td_ms": td_ms,
                        "b_step": b_step,
                        "N": n_value,
                        "source_file": f"n{n_value}.xlsx",
                        "value_norm": signal,
                        "g": gradient,
                        "g_thorsten": gradient * 10.0,
                        "grad_correction_factor": 1.0,
                    }
                )

        contrast, summary = build_resampled_contrasts(
            pd.DataFrame(rows),
            d0_m2_ms=d0_m2_ms,
            directions=["tra"],
            rois=["fiber1"],
            grid_size=128,
            tc_mode="shared",
        )
        shape = summarize_contrast_shape(contrast)

        self.assertEqual(summary.iloc[0]["tc_model"], "shared")
        self.assertEqual(int(summary.iloc[0]["n_parameters"]), 2)
        self.assertAlmostEqual(float(summary.iloc[0]["tc_shared_ms"]), expected_tc_ms, places=4)
        self.assertAlmostEqual(float(summary.iloc[0]["tc_high_ms"]), float(summary.iloc[0]["tc_low_ms"]))
        self.assertTrue(bool(shape.iloc[0]["shape_ok"]))
        self.assertGreater(float(shape.iloc[0]["contrast_auc_normalized_g"]), 0.0)
        self.assertGreater(float(shape.iloc[0]["lcf_centroid_um"]), 0.0)

    def test_pseudohuber_recovers_delta_with_fixed_alpha(self) -> None:
        td = np.array([90.0, 120.0, 143.4, 210.0])
        alpha = 0.25
        expected_delta = 350.0
        tc_peak = _tc_pseudohuber(td, 5.0, expected_delta, alpha)
        contrast_fits = pd.DataFrame(
            {
                "variant": "plain",
                "subj": "S1",
                "sheet": "sheet1",
                "roi": "PostCC",
                "direction": "long",
                "td_ms": td,
                "tc_peak_ms": tc_peak,
                "ok": True,
            }
        )
        alpha_table = pd.DataFrame(
            {
                "variant": ["plain"],
                "subj": ["S1"],
                "roi": ["PostCC"],
                "direction": ["long"],
                "alpha_macro": [alpha],
                "alpha_macro_error": [0.01],
            }
        )

        _, summary = fit_tc_pseudohuber(contrast_fits, alpha_table, exclude_td_ms=())

        self.assertTrue(bool(summary.iloc[0]["ok"]))
        self.assertAlmostEqual(float(summary.iloc[0]["delta_ms"]), expected_delta, places=3)


if __name__ == "__main__":
    unittest.main()
