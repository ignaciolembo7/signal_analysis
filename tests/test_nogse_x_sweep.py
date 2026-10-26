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

from data_processing.reshape import to_long  # noqa: E402
from data_processing.result_signals import (  # noqa: E402
    CleanSequenceParams,
    add_s0_and_value_norm,
    add_x_sweep_derivatives,
    is_x_sweep_params,
)
from fitting.b_from_g import nogse_trapezoid_b_per_g2, nogse_x_from_b_per_g2  # noqa: E402

# x sweep of nogse_v2003_fov220_m128_sl01_R1_autox/T050_g030_3dir20x_AP_rep01 (sequence version v2.00.3):
# T = 50 ms, N = 2, ramps 0.32 ms, G = 30 mT/m; the protocol x (ms) and the .bval written by the Pulseq generator.
X_MS = [0.0, 1.32, 2.64, 3.94, 5.26, 6.58, 7.90, 9.22, 10.52, 11.84, 13.16, 14.48, 15.78, 17.10, 18.42, 19.74,
        21.06, 22.36, 23.68, 25.0]
B_PER_X = [657.996813, 606.927261, 558.626234, 513.764062, 470.960372, 430.926332, 393.661943, 359.167204,
           327.902135, 298.904733, 272.676981, 249.21888, 228.823227, 210.862463, 195.671349, 183.249886,
           173.598073, 166.799522, 162.645046, 161.260221]


def autox_params(**changes) -> CleanSequenceParams:
    values = dict(
        sheet="PHANTOM", subj="phantom", protocol="T050_g030_3dir20x_AP_rep01", sequence=5, group=None,
        G=30.0, TN=50.0, x=np.nan, y=np.nan, type="AUTOX", Hz=np.nan, bmax=np.nan, N=2,
        delta_ms=np.nan, Delta_app_ms=np.nan, max_dur_ms=np.nan, tm_ms=np.nan, td_ms=50.0,
        TE=152.86, TR=1000.0, g_thorsten=None, ramp_time_ms=0.32,
    )
    values.update(changes)
    return CleanSequenceParams(**values)


class NogseXSweepTests(unittest.TestCase):
    def test_x_from_b_inverts_the_trapezoidal_waveform(self) -> None:
        for x, b in zip(X_MS, B_PER_X):
            self.assertAlmostEqual(
                nogse_x_from_b_per_g2(b / 30.0**2, N=2, TN_ms=50.0, ramp_time_ms=0.32), x, places=5
            )
        # Between the pure Hahn limit and x = 4 ramps no waveform can be played.
        k0 = nogse_trapezoid_b_per_g2(N=2, TN_ms=50.0, x_ms=0.0, ramp_time_ms=0.32)
        k_min_x = nogse_trapezoid_b_per_g2(N=2, TN_ms=50.0, x_ms=1.28, ramp_time_ms=0.32)
        with self.assertRaises(ValueError):
            nogse_x_from_b_per_g2(0.5 * (k0 + k_min_x), N=2, TN_ms=50.0, ramp_time_ms=0.32)

    def test_x_sweep_is_recognised_by_type_or_name(self) -> None:
        self.assertTrue(is_x_sweep_params(autox_params()))
        self.assertTrue(is_x_sweep_params(autox_params(type="NOGSE")))
        self.assertFalse(is_x_sweep_params(autox_params(type="HAHN", protocol="T050_xhahn_3dir10b_AP_rep01")))

    def test_results_table_gets_x_per_volume_and_fixed_g(self) -> None:
        # Results table as signal_extraction writes it: 2 b0, then the triples x, y, z of every x.
        bvalues = [0.0, 0.0] + [b for b in B_PER_X for _ in range(3)]
        signal = 600.0 * np.exp(-np.asarray(bvalues) * 2.15e-3)
        stats = {stat: pd.DataFrame({"bvalues": bvalues, "water": signal}) for stat in ("mean", "std")}
        long = to_long(stats, ndirs=3, nbvals=20, source_file="autox_results.xlsx")
        out = add_s0_and_value_norm(add_x_sweep_derivatives(long, autox_params()))
        weighted = out.loc[out["b_step"] > 0]
        for step, x in enumerate(X_MS, start=1):
            np.testing.assert_allclose(weighted.loc[weighted["b_step"] == step, "x"], x, atol=1e-5)
        np.testing.assert_allclose(weighted["y"], 50.0 - weighted["x"], atol=1e-12)
        np.testing.assert_allclose(weighted["g"], 30.0)
        np.testing.assert_allclose(out.loc[out["b_step"] == 0, "g"], 0.0)
        self.assertTrue(out.loc[out["b_step"] == 0, "x"].isna().all())
        np.testing.assert_allclose(out["value_norm"], np.exp(-out["bvalue"].to_numpy(float) * 2.15e-3), rtol=1e-12)

    def test_missing_ramp_is_rejected(self) -> None:
        long = to_long({"mean": pd.DataFrame({"bvalues": [0.0, 0.0] + [B_PER_X[0]] * 3, "water": 1.0})},
                       ndirs=3, nbvals=1)
        with self.assertRaises(ValueError):
            add_x_sweep_derivatives(long, autox_params(ramp_time_ms=np.nan))


if __name__ == "__main__":
    unittest.main()
