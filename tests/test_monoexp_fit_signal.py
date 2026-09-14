from __future__ import annotations

from pathlib import Path
import sys
import unittest

import numpy as np
import pandas as pd

REPO_ROOT = Path(__file__).resolve().parents[1]
SRC_ROOT = REPO_ROOT / "src"
if str(SRC_ROOT) not in sys.path:
    sys.path.insert(0, str(SRC_ROOT))

from monoexp_fitting.fit_signal import fit_signal_monoexp, monoexp, select_monoexp_fit_result  # noqa: E402


class MonoexpFitSignalTests(unittest.TestCase):
    def test_auto_fit_points_keeps_largest_prefix_within_tolerance(self) -> None:
        b = np.array([0.0, 100.0, 200.0, 300.0, 400.0])
        y = monoexp(b, 1.0, 0.0023)
        y[-1] = 0.95

        result = select_monoexp_fit_result(
            b,
            y,
            auto_fit_points=True,
            auto_fit_min_points=3,
            auto_fit_max_points=5,
            auto_fit_rel_tol=0.05,
            auto_fit_err_floor=0.005,
            fix_M0=1.0,
            D0_init=0.0023,
        )

        self.assertTrue(result["ok"], result.get("msg"))
        self.assertEqual(result["fit_points"], 4)
        self.assertEqual(result["n_fit"], 4)
        self.assertIn("First rejected prefix: k=5", result["msg"])

    def test_auto_fit_points_is_not_anchored_to_an_artificially_perfect_short_prefix(self) -> None:
        b = np.array([0.0, 100.0, 200.0, 300.0, 400.0, 500.0])
        y = monoexp(b, 1.0, 0.0023)
        y[3:5] *= np.array([1.025, 0.975])
        y[5] = 0.9

        result = select_monoexp_fit_result(
            b,
            y,
            auto_fit_points=True,
            auto_fit_min_points=3,
            auto_fit_max_points=6,
            auto_fit_rel_tol=0.05,
            auto_fit_err_floor=0.005,
            fix_M0=1.0,
            D0_init=0.0023,
        )

        self.assertTrue(result["ok"], result.get("msg"))
        self.assertEqual(result["fit_points"], 5)
        self.assertEqual(result["n_fit"], 5)

    def test_fixed_fit_uses_requested_leading_points(self) -> None:
        b = np.array([0.0, 100.0, 200.0, 300.0, 400.0])
        y = monoexp(b, 1.0, 0.0023)

        result = select_monoexp_fit_result(
            b,
            y,
            fit_points=5,
            auto_fit_points=False,
            fix_M0=1.0,
            D0_init=0.0023,
        )

        self.assertTrue(result["ok"], result.get("msg"))
        self.assertEqual(result["fit_strategy"], "fixed")
        self.assertEqual(result["fit_points"], 5)
        self.assertEqual(result["n_fit"], 5)

    def test_fit_signal_monoexp_marks_used_points(self) -> None:
        b = np.array([0.0, 100.0, 200.0, 300.0, 400.0])
        y = monoexp(b, 1.0, 0.0023)
        y[-1] = 0.95
        df = pd.DataFrame(
            {
                "subj": ["S1"] * len(b),
                "sheet": ["S1"] * len(b),
                "td_ms": [90.0] * len(b),
                "N": [4] * len(b),
                "Hz": [25.0] * len(b),
                "roi": ["ROI_A"] * len(b),
                "direction": ["long"] * len(b),
                "stat": ["avg"] * len(b),
                "b_step": np.arange(len(b)),
                "bvalue_thorsten": b,
                "value_norm": y,
            }
        )

        outputs = fit_signal_monoexp(
            df,
            b_axis="bvalue_thorsten",
            auto_fit_points=True,
            auto_fit_min_points=3,
            auto_fit_max_points=5,
            auto_fit_rel_tol=0.05,
            auto_fit_err_floor=0.005,
            fix_M0=1.0,
            D0_init=0.0023,
        )

        self.assertEqual(len(outputs.fit_params), 1)
        self.assertEqual(int(outputs.fit_params.loc[0, "fit_points"]), 4)
        self.assertEqual(int(outputs.fit_points["used_for_fit"].sum()), 4)

    def test_fit_signal_uses_group_gradient_correction_factor(self) -> None:
        b = np.array([0.0, 100.0, 200.0, 300.0])
        y = monoexp(b * 4.0, 1.0, 0.0023)
        df = pd.DataFrame(
            {
                "subj": ["S1"] * len(b),
                "sheet": ["S1"] * len(b),
                "td_ms": [90.0] * len(b),
                "N": [4] * len(b),
                "Hz": [25.0] * len(b),
                "roi": ["ROI_A"] * len(b),
                "direction": ["long"] * len(b),
                "stat": ["avg"] * len(b),
                "b_step": np.arange(len(b)),
                "bvalue_thorsten": b,
                "value_norm": y,
                "grad_correction_factor": [2.0] * len(b),
            }
        )

        outputs = fit_signal_monoexp(
            df,
            b_axis="bvalue_thorsten",
            auto_fit_points=False,
            fit_points=4,
            correction_factor_col="grad_correction_factor",
            fix_M0=1.0,
            D0_init=0.0023,
        )

        self.assertTrue(bool(outputs.fit_params.loc[0, "ok"]))
        self.assertAlmostEqual(float(outputs.fit_params.loc[0, "D0_mm2_s"]), 0.0023, places=6)
        self.assertEqual(outputs.fit_points["bvalue_used"].tolist(), (b * 4.0).tolist())

    def test_fit_signal_marks_curve_failed_when_gradient_correction_factor_is_missing(self) -> None:
        b = np.array([0.0, 100.0, 200.0, 300.0])
        y = monoexp(b, 1.0, 0.0023)
        df = pd.DataFrame(
            {
                "subj": ["S1"] * len(b),
                "sheet": ["S1"] * len(b),
                "td_ms": [90.0] * len(b),
                "N": [4] * len(b),
                "Hz": [25.0] * len(b),
                "roi": ["ROI_A"] * len(b),
                "direction": ["x"] * len(b),
                "stat": ["avg"] * len(b),
                "b_step": np.arange(len(b)),
                "bvalue_thorsten": b,
                "value_norm": y,
                "grad_correction_factor": [np.nan] * len(b),
            }
        )

        outputs = fit_signal_monoexp(
            df,
            b_axis="bvalue_thorsten",
            auto_fit_points=True,
            correction_factor_col="grad_correction_factor",
            fix_M0=1.0,
            D0_init=0.0023,
        )

        self.assertFalse(bool(outputs.fit_params.loc[0, "ok"]))
        self.assertIn("Missing grad_correction_factor", outputs.fit_params.loc[0, "msg"])
        self.assertEqual(int(outputs.fit_points["used_for_fit"].sum()), 0)


if __name__ == "__main__":
    unittest.main()
