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

from models.model_fitting import M_ogse_free  # noqa: E402
from ogse_fitting.fit_ogse_free_signal import fit_ogse_free_signal  # noqa: E402


class OgseFreeSignalFitTests(unittest.TestCase):
    def setUp(self) -> None:
        self.td_ms = 90.0
        self.n_value = 4
        self.true_d0_m2_ms = 2.3e-12
        self.factor = 2.0
        gradients = np.array([0.0, 5.0, 10.0, 15.0, 20.0])
        signal = M_ogse_free(
            self.td_ms,
            gradients * self.factor,
            self.n_value,
            self.td_ms / self.n_value,
            1.0,
            self.true_d0_m2_ms,
        )
        self.frame = pd.DataFrame(
            {
                "subj": ["S1"] * len(gradients),
                "sheet": ["SHEET"] * len(gradients),
                "td_ms": [self.td_ms] * len(gradients),
                "N": [self.n_value] * len(gradients),
                "Hz": [20.0] * len(gradients),
                "roi": ["ROI"] * len(gradients),
                "direction": ["long"] * len(gradients),
                "stat": ["avg"] * len(gradients),
                "b_step": np.arange(len(gradients)),
                "g": gradients,
                "value_norm": signal,
                "grad_correction_factor": [self.factor] * len(gradients),
            }
        )

    def test_raw_gradient_fit_absorbs_squared_gradient_error_into_d0(self) -> None:
        outputs = fit_ogse_free_signal(self.frame, fix_M0=1.0)

        fitted = float(outputs.fit_params.loc[0, "D0_m2_ms"])
        self.assertAlmostEqual(fitted / self.true_d0_m2_ms, self.factor**2, places=5)
        self.assertAlmostEqual(float(outputs.fit_params.loc[0, "f_corr"]), 1.0)
        self.assertEqual(outputs.fit_params.loc[0, "model"], "ogse_free")

    def test_corrected_gradient_fit_recovers_true_d0(self) -> None:
        outputs = fit_ogse_free_signal(
            self.frame,
            fix_M0=1.0,
            correction_factor_col="grad_correction_factor",
        )

        fitted = float(outputs.fit_params.loc[0, "D0_m2_ms"])
        self.assertAlmostEqual(fitted / self.true_d0_m2_ms, 1.0, places=5)
        self.assertAlmostEqual(float(outputs.fit_params.loc[0, "f_corr"]), self.factor)
        self.assertTrue(outputs.fit_points["used_for_fit"].all())
        self.assertTrue((outputs.fit_points["fit_kind"] == "ogse_free_signal").all())


if __name__ == "__main__":
    unittest.main()
