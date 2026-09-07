from __future__ import annotations

from pathlib import Path
import sys
import unittest

import pandas as pd

REPO_ROOT = Path(__file__).resolve().parents[1]
SRC_ROOT = REPO_ROOT / "src"
if str(SRC_ROOT) not in sys.path:
    sys.path.insert(0, str(SRC_ROOT))

from tc_fittings.alpha_macro_summary import compute_alpha_macro_summary  # noqa: E402


def _dproj_table() -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    for roi in ["Left-Lateral-Ventricle", "Syringe"]:
        for bvalue in [100.0, 500.0, 980.0, 1280.0, 2000.0]:
            for delta in [20.0, 40.0]:
                rows.append(
                    {
                        "subj": "MBBL",
                        "sheets": "20230630_MBBL-3",
                        "roi": roi,
                        "direction": "long",
                        "bvalue": bvalue,
                        "Delta_app_ms": delta,
                        "D_mean_mm2_s": bvalue / 1_000_000.0,
                    }
                )
    return pd.DataFrame(rows)


class AlphaMacroSummaryTests(unittest.TestCase):
    def test_bvalmax_selects_candidate_bstep_from_plot_bvalues(self) -> None:
        _, summary = compute_alpha_macro_summary(
            _dproj_table(),
            selected_bstep=4,
            candidate_bvalues=[500.0, 980.0, 1280.0, 2000.0],
            direction_aliases={},
        )

        rows = summary[summary["roi"] == "Left-Lateral-Ventricle"]
        self.assertEqual(float(rows.iloc[0]["selected_bvalue"]), 2000.0)
        self.assertEqual(int(rows.iloc[0]["selected_bstep"]), 5)

    def test_roi_bvalmax_can_select_bvalue_from_plot_bvalues(self) -> None:
        _, summary = compute_alpha_macro_summary(
            _dproj_table(),
            selected_bstep=2000,
            roi_selected_bsteps={"Left-Lateral-Ventricle": 500.0, "Syringe": 1280.0},
            candidate_bvalues=[500.0, 980.0, 1280.0, 2000.0],
            direction_aliases={},
        )

        by_roi = summary.set_index("roi")
        self.assertEqual(float(by_roi.loc["Left-Lateral-Ventricle", "selected_bvalue"]), 500.0)
        self.assertEqual(int(by_roi.loc["Left-Lateral-Ventricle", "selected_bstep"]), 2)
        self.assertEqual(float(by_roi.loc["Syringe", "selected_bvalue"]), 1280.0)
        self.assertEqual(int(by_roi.loc["Syringe", "selected_bstep"]), 4)


if __name__ == "__main__":
    unittest.main()
