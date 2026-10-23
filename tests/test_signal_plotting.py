from __future__ import annotations

from pathlib import Path
import sys
import tempfile
import unittest

import numpy as np
import pandas as pd

REPO_ROOT = Path(__file__).resolve().parents[1]
SRC_ROOT = REPO_ROOT / "src"
if str(SRC_ROOT) not in sys.path:
    sys.path.insert(0, str(SRC_ROOT))

from plotting.ogse.signal_vs_g import _prepare_avg_std, filter_plot_n_values  # noqa: E402
from plotting.nogse.signal_vs_g import (  # noqa: E402
    _prepare_avg_std as _prepare_nogse_avg_std,
    plot_nogse_signal_table,
)


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


class NogseSignalPlottingTests(unittest.TestCase):
    @staticmethod
    def _mixed_type_frame() -> pd.DataFrame:
        common = {
            "subj": "PHANTOM-1",
            "roi": "water",
            "N": 2,
            "td_ms": 30.0,
            "direction": "1",
            "b_step": 1,
            "g": 20.0,
            "S0": 100.0,
        }
        return pd.DataFrame(
            [
                {**common, "type": "CPMG", "stat": "avg", "value": 80.0, "value_norm": 0.8},
                {**common, "type": "CPMG", "stat": "std", "value": 8.0, "value_norm": 0.08},
                {**common, "type": "HAHN", "stat": "avg", "value": 60.0, "value_norm": 0.6},
                {**common, "type": "HAHN", "stat": "std", "value": 6.0, "value_norm": 0.06},
            ]
        )

    def test_avg_std_join_keeps_hahn_and_cpmg_independent(self) -> None:
        plotted = _prepare_nogse_avg_std(
            self._mixed_type_frame(),
            xcol="g",
            ycol="value_norm",
            stat="avg",
        ).sort_values("type")

        self.assertEqual(plotted["type"].tolist(), ["CPMG", "HAHN"])
        self.assertEqual(plotted["y_mean"].tolist(), [0.8, 0.6])
        self.assertEqual(plotted["y_std"].tolist(), [0.08, 0.06])

    def test_plot_filenames_are_split_by_acquisition_type(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            paths = plot_nogse_signal_table(
                self._mixed_type_frame(),
                out_root=Path(tmp),
                analysis_id="",
                xcol="g",
                ycol="value_norm",
                stat="avg",
                rois=None,
                directions=None,
            )

        names = [path.name for path in paths]
        self.assertEqual(len(names), 4)
        self.assertTrue(any("type-CPMG" in name for name in names))
        self.assertTrue(any("type-HAHN" in name for name in names))


if __name__ == "__main__":
    unittest.main()
