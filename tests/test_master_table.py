from __future__ import annotations

from pathlib import Path
import sys
import tempfile
import unittest

import pandas as pd

REPO_ROOT = Path(__file__).resolve().parents[1]
SRC_ROOT = REPO_ROOT / "src"
if str(SRC_ROOT) not in sys.path:
    sys.path.insert(0, str(SRC_ROOT))

from data_processing.master_table import (  # noqa: E402
    append_master_rows,
    build_analysis_id_from_columns,
    filter_table_rows,
    load_master_table,
    select_alpha_macro,
    select_contrast_pair,
    select_fit_params,
    select_signal,
    write_master_table,
)


def _signal_rows() -> pd.DataFrame:
    return pd.DataFrame(
        {
            "stat": ["avg", "avg"],
            "roi": ["ROI_A", "ROI_A"],
            "direction": ["long", "long"],
            "b_step": [0, 1],
            "bvalue": [0.0, 1000.0],
            "value": [10.0, 7.0],
            "td_ms": [40.0, 40.0],
            "N": [2, 2],
            "Hz": [50.0, 50.0],
            "subj": ["S1", "S1"],
            "sheet": ["S1", "S1"],
            "source_file": ["legacy.long.parquet", "legacy.long.parquet"],
        }
    )


def _contrast_rows() -> pd.DataFrame:
    return pd.DataFrame(
        {
            "row_kind": ["contrast", "contrast"],
            "analysis_id": ["S1_N4-N2_td40", "S1_N4-N2_td40"],
            "stat": ["avg", "avg"],
            "roi": ["ROI_A", "ROI_A"],
            "direction": ["long", "long"],
            "b_step": [0, 1],
            "value": [0.0, 0.25],
            "value_norm": [0.0, 0.25],
            "td_ms": [40.0, 40.0],
            "subj": ["S1", "S1"],
            "sheet": ["S1", "S1"],
        }
    )


class MasterTableTests(unittest.TestCase):
    def test_append_and_select_signal(self) -> None:
        master = append_master_rows(None, _signal_rows(), row_kind="signal_rotated", analysis_id="S1_td40_N2")

        selected = select_signal(master, rotated=True, subj="S1", roi="ROI_A", direction="long")

        self.assertEqual(len(selected), 2)
        self.assertEqual(selected["row_kind"].unique().tolist(), ["signal_rotated"])
        self.assertIn("value_norm", selected.columns)
        self.assertEqual(float(selected.loc[selected["b_step"] == 0, "value_norm"].iloc[0]), 1.0)

    def test_select_contrast_pair(self) -> None:
        master = append_master_rows(None, _signal_rows(), row_kind="signal", analysis_id="S1_td40_N2")
        master = append_master_rows(master, _contrast_rows())

        selected = select_contrast_pair(master, subj="S1", roi="ROI_A", direction="long")

        self.assertEqual(len(selected), 2)
        self.assertEqual(selected["analysis_id"].unique().tolist(), ["S1_N4-N2_td40"])

    def test_write_and_load_master_table(self) -> None:
        master = append_master_rows(None, _signal_rows(), row_kind="signal", analysis_id="S1_td40_N2")
        with tempfile.TemporaryDirectory() as tmp:
            out_path = Path(tmp) / "master.long.parquet"
            write_master_table(master, out_path)
            loaded = load_master_table(out_path)

        self.assertEqual(len(loaded), len(master))
        self.assertIn("row_kind", loaded.columns)

    def test_numeric_direction_is_stored_as_selectable_label(self) -> None:
        rows = _signal_rows().copy()
        rows["direction"] = [1, 1]
        master = append_master_rows(None, rows, row_kind="signal", analysis_id="S1_td40_N2")

        selected = select_signal(master, rotated=False, subj="S1", roi="ROI_A", direction=1)

        self.assertEqual(len(selected), 2)
        self.assertEqual(selected["direction"].unique().tolist(), ["1"])

    def test_build_analysis_id_from_columns_requires_unique_values(self) -> None:
        analysis_id = build_analysis_id_from_columns(
            _signal_rows(),
            columns=("subj", "td_ms", "N", "direction"),
            prefix="signal",
        )

        self.assertEqual(analysis_id, "signal_subj-S1_td_ms-40_N-2_direction-long")

        mixed = _signal_rows().copy()
        mixed.loc[1, "N"] = 4
        with self.assertRaises(ValueError):
            build_analysis_id_from_columns(mixed, columns=("subj", "N"))

    def test_select_fit_params_from_cumulative_table(self) -> None:
        table = pd.DataFrame(
            {
                "fit_kind": ["ogse_contrast", "ogse_signal"],
                "model": ["ogse_free", "ogse_free"],
                "subj": ["S1", "S1"],
                "roi": ["ROI_A", "ROI_A"],
                "direction": ["long", "tra"],
                "td_ms": [40.0, 40.0],
                "ok": [True, True],
            }
        )

        selected = select_fit_params(table, fit_kind="ogse_contrast", subj="S1", direction="long", ok_only=True)

        self.assertEqual(len(selected), 1)
        self.assertEqual(selected["roi"].iloc[0], "ROI_A")

    def test_select_alpha_macro_from_cumulative_table(self) -> None:
        table = pd.DataFrame(
            {
                "fit_kind": ["alpha_macro_summary", "ogse_contrast"],
                "subj": ["S1", "S1"],
                "roi": ["ROI_A", "ROI_A"],
                "direction": ["long", "long"],
                "td_ms": [40.0, 40.0],
                "alpha_macro": [0.22, float("nan")],
                "alpha_macro_se": [0.01, float("nan")],
            }
        )

        selected = select_alpha_macro(table, subj="S1", roi="ROI_A", direction="long")

        self.assertEqual(len(selected), 1)
        self.assertAlmostEqual(float(selected["alpha_macro"].iloc[0]), 0.22)
        self.assertIn("alpha_macro_error", selected.columns)

    def test_filter_table_rows_numeric_selector(self) -> None:
        table = pd.DataFrame({"td_ms": [40.0, 90.0], "subj": ["S1", "S2"]})

        selected = filter_table_rows(table, td_ms=40)

        self.assertEqual(selected["subj"].tolist(), ["S1"])


if __name__ == "__main__":
    unittest.main()
