from __future__ import annotations

import json
from pathlib import Path
import sys
import tempfile
import unittest

import numpy as np
import pandas as pd

REPO_ROOT = Path(__file__).resolve().parents[1]
SRC_ROOT = REPO_ROOT / "src"
SCRIPT_DATA_ROOT = REPO_ROOT / "scripts" / "data"
if str(SRC_ROOT) not in sys.path:
    sys.path.insert(0, str(SRC_ROOT))
if str(SCRIPT_DATA_ROOT) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DATA_ROOT))

from data_processing.result_signals import (  # noqa: E402
    CleanSequenceParams,
    add_g_and_b_derivatives,
    is_nogse_params,
    resolve_nogse_b_per_g2,
)
from fitting.b_from_g import (  # noqa: E402
    b_from_g,
    g_from_b_nogse,
    nogse_rect_b_per_g2,
    nogse_trapezoid_b_per_g2,
)
from models.model_fitting import _GYRO, M_nogse_free  # noqa: E402

GAMMA = 267.5221900  # 1/(ms*mT), default of process_one_results.py


def nogse_params(**changes) -> CleanSequenceParams:
    values = dict(
        sheet="PULSEQ_NOGSE", subj="phantom", protocol="NOGSE_TN50", sequence=1, group=None,
        G=78.0, TN=50.0, x=25.0, y=25.0, type="NOGSE", Hz=np.nan, bmax=np.nan, N=2,
        delta_ms=np.nan, Delta_app_ms=np.nan, max_dur_ms=np.nan, tm_ms=np.nan, td_ms=np.nan,
        TE=np.nan, TR=np.nan, g_thorsten=None, ramp_time_ms=1.2,
    )
    values.update(changes)
    return CleanSequenceParams(**values)


def b_long(bvalues) -> pd.DataFrame:
    b = np.asarray(bvalues, dtype=float)
    steps = np.searchsorted(np.unique(b), b)
    return pd.DataFrame({"stat": "avg", "roi": "ROI", "direction": "1", "b_step": steps, "bvalue": b, "value": 1.0})


class NogseBToGTests(unittest.TestCase):
    def test_rectangular_b_matches_the_nogse_free_model(self) -> None:
        # ln M = -b D0 with b in s/mm^2 and D0 in mm^2/s; the model uses D0 in m^2/ms.
        D0_mm2_s = 2.0e-3
        for N, TN, x in ((2, 50.0, 25.0), (2, 80.0, 2.5), (4, 40.0, 8.0)):
            b = nogse_rect_b_per_g2(N=N, TN_ms=TN, x_ms=x, gamma=_GYRO) * 60.0**2  # gamma of the model
            M = M_nogse_free(TN, 60.0, N, x, 1.0, D0_mm2_s * 1e-9)
            self.assertAlmostEqual(-np.log(M), b * D0_mm2_s, places=9)

    def test_g_from_b_inverts_b(self) -> None:
        b_per_g2 = nogse_rect_b_per_g2(N=2, TN_ms=50.0, x_ms=25.0, gamma=GAMMA)
        g = np.array([0.0, 7.8, 39.0, 78.0])
        np.testing.assert_allclose(g_from_b_nogse(b_per_g2 * g**2, b_per_g2=b_per_g2), g, rtol=1e-12)
        with self.assertRaises(ValueError):
            g_from_b_nogse(g, b_per_g2=0.0)
        with self.assertRaises(ValueError):
            nogse_rect_b_per_g2(N=2, TN_ms=10.0, x_ms=25.0, gamma=GAMMA)  # y < 0

    def test_trapezoidal_formula_matches_played_phantom_sequences(self) -> None:
        expected = {
            (30.0, 15.0): 0.0304675246546,
            (30.0, 5.0): 0.0798460984326,
            (40.0, 20.0): 0.0780307798506,
            (40.0, 5.0): 0.229386842953,
            (50.0, 25.0): 0.159228715736,
            (50.0, 5.0): 0.499869311654,
            (70.0, 35.0): 0.458310204778,
            (70.0, 5.0): 1.546785722396,
        }
        for (TN, x), coefficient in expected.items():
            actual = nogse_trapezoid_b_per_g2(N=2, TN_ms=TN, x_ms=x, ramp_time_ms=1.2)
            self.assertAlmostEqual(actual, coefficient, places=10)

    def test_pure_hahn_limit_matches_played_v2_00_2_sequences(self) -> None:
        # x = 0 (sequence version v2.00.2): two Hahn lobes of TN/2; b/g^2 of the _diffusion.json sidecars.
        expected = {30.0: 0.141503363056, 40.0: 0.347042272766, 50.0: 0.691413694259, 65.0: 1.546553981841}
        for TN, coefficient in expected.items():
            actual = nogse_trapezoid_b_per_g2(N=2, TN_ms=TN, x_ms=0.0, ramp_time_ms=1.2)
            self.assertAlmostEqual(actual, coefficient, places=10)
        # Without ramps, the rectangular Hahn limit gamma^2 TN^3 / 12.
        rect = nogse_rect_b_per_g2(N=2, TN_ms=50.0, x_ms=0.0, gamma=GAMMA)
        self.assertAlmostEqual(
            nogse_trapezoid_b_per_g2(N=2, TN_ms=50.0, x_ms=0.0, ramp_time_ms=1e-6, gamma=GAMMA), rect, places=6
        )
        with self.assertRaises(ValueError):
            nogse_trapezoid_b_per_g2(N=2, TN_ms=50.0, x_ms=-1.0, ramp_time_ms=1.2)

    def test_nogse_rows_are_detected(self) -> None:
        self.assertTrue(is_nogse_params(nogse_params()))
        self.assertTrue(is_nogse_params(nogse_params(type=None, sheet="S", protocol=None)))  # TN and x only
        ogse = nogse_params(type="OGSE", sheet="OGSE_d40", protocol="OGSE", TN=np.nan, x=np.nan,
                            delta_ms=10.0, Delta_app_ms=20.0)
        self.assertFalse(is_nogse_params(ogse))

    def test_waveform_formula_priority_and_legacy_fallbacks(self) -> None:
        rect = nogse_rect_b_per_g2(N=2, TN_ms=50.0, x_ms=25.0, gamma=GAMMA)
        exact = nogse_trapezoid_b_per_g2(N=2, TN_ms=50.0, x_ms=25.0, ramp_time_ms=1.2)
        self.assertEqual(
            resolve_nogse_b_per_g2(nogse_params(), gamma=GAMMA),
            (exact, "trapezoidal_waveform(N,TN,x,ramp)"),
        )
        with tempfile.TemporaryDirectory() as tmp:
            summary = Path(tmp) / "seq_diffusion.json"
            summary.write_text(json.dumps({"b_per_g2_s_per_mm2_per_mt2_m2": 0.9 * rect}), encoding="utf-8")
            value, source = resolve_nogse_b_per_g2(
                nogse_params(ramp_time_ms=np.nan), gamma=GAMMA, diffusion_summary=summary
            )
            self.assertAlmostEqual(value, 0.9 * rect)
            self.assertEqual(source, str(summary))
            # Older workbooks retain their previous priority when no ramp is known.
            value, source = resolve_nogse_b_per_g2(
                nogse_params(ramp_time_ms=np.nan, b_per_g2=1.5), gamma=GAMMA, diffusion_summary=summary
            )
            self.assertEqual((value, source), (1.5, "sequence_parameters"))

            with self.assertRaisesRegex(ValueError, "disagrees"):
                resolve_nogse_b_per_g2(nogse_params(), gamma=GAMMA, diffusion_summary=summary)

        self.assertEqual(
            resolve_nogse_b_per_g2(nogse_params(ramp_time_ms=np.nan), gamma=GAMMA),
            (rect, "rectangular_nogse"),
        )

    def test_bvalue_axis_gives_the_played_nogse_g(self) -> None:
        # Exact (trapezoid) b of the Pulseq sequence: b_per_g2 smaller than the rectangular value.
        b_per_g2 = 0.97 * nogse_rect_b_per_g2(N=2, TN_ms=50.0, x_ms=25.0, gamma=GAMMA)
        g_played = np.array([0.0, 0.0, 7.8, 15.6, 78.0])
        out = add_g_and_b_derivatives(
            b_long(b_per_g2 * g_played**2), gamma=GAMMA, N=2, delta_ms=np.nan, Delta_app_ms=np.nan,
            g_thorsten=None, nogse_b_per_g2=b_per_g2,
        )
        np.testing.assert_allclose(out["g"], g_played, rtol=1e-12, atol=1e-12)
        np.testing.assert_allclose(out["bvalue_g"], out["bvalue"], rtol=1e-12, atol=1e-12)
        np.testing.assert_allclose(out["g_max"], 78.0)

    def test_ogse_conversion_is_unchanged(self) -> None:
        g = np.array([0.0, 20.0, 40.0])
        b = b_from_g(g, N=2, gamma=GAMMA, delta_ms=10.0, delta_app_ms=20.0, g_type="g")
        out = add_g_and_b_derivatives(b_long(b), gamma=GAMMA, N=2, delta_ms=10.0, Delta_app_ms=20.0, g_thorsten=None)
        np.testing.assert_allclose(out["g"], g, rtol=1e-12)


if __name__ == "__main__":
    unittest.main()
