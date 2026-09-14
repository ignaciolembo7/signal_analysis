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

from plotting.ogse.signal_vs_g import _prepare_avg_std, filter_plot_n_values  # noqa: E402


class OgseSignalPlottingTests(unittest.TestCase):
    def test_default_plot_filter_keeps_only_supported_n_values(self) -> None:
        frame = pd.DataFrame({"N": [1, 2, 4, 8, 12, 16], "value": range(6)})

        filtered = filter_plot_n_values(frame)

        self.assertEqual(filtered["N"].tolist(), [1, 4, 8, 12])

    def test_normalized_error_uses_raw_std_over_mean_s0(self) -> None:
        common = {
            "subj": "S1",
            "roi": "ROI",
            "N": 4,
            "td_ms": 90.0,
            "direction": "long",
            "b_step": 1,
            "g_thorsten": 20.0,
        }
        frame = pd.DataFrame(
            [
                {**common, "stat": "avg", "value": 50.0, "value_norm": 0.5, "S0": 100.0},
                {**common, "stat": "std", "value": 10.0, "value_norm": 1.0, "S0": 10.0},
            ]
        )

        plotted = _prepare_avg_std(
            frame,
            xcol="g_thorsten",
            ycol="value_norm",
            stat="avg",
        )

        self.assertAlmostEqual(float(plotted.loc[0, "y_mean"]), 0.5)
        self.assertAlmostEqual(float(plotted.loc[0, "y_std"]), 0.1)

    def test_raw_signal_error_remains_in_raw_units(self) -> None:
        common = {
            "subj": "S1",
            "roi": "ROI",
            "N": 4,
            "td_ms": 90.0,
            "direction": "long",
            "b_step": 1,
            "g_thorsten": 20.0,
        }
        frame = pd.DataFrame(
            [
                {**common, "stat": "avg", "value": 50.0, "value_norm": 0.5, "S0": 100.0},
                {**common, "stat": "std", "value": 10.0, "value_norm": 1.0, "S0": 10.0},
            ]
        )

        plotted = _prepare_avg_std(frame, xcol="g_thorsten", ycol="value", stat="avg")

        self.assertAlmostEqual(float(plotted.loc[0, "y_mean"]), 50.0)
        self.assertAlmostEqual(float(plotted.loc[0, "y_std"]), 10.0)


if __name__ == "__main__":
    unittest.main()
