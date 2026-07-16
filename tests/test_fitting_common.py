from __future__ import annotations

import argparse
from pathlib import Path
import sys
import tempfile
import unittest

import pandas as pd

REPO_ROOT = Path(__file__).resolve().parents[1]
SRC_ROOT = REPO_ROOT / "src"
if str(SRC_ROOT) not in sys.path:
    sys.path.insert(0, str(SRC_ROOT))

from data_processing.master_table import append_master_rows, write_master_table  # noqa: E402
from fitting.cli_common import build_common_parameter_plan, load_master_input  # noqa: E402


def _signal_rows() -> pd.DataFrame:
    return pd.DataFrame(
        {
            "stat": ["avg", "avg"],
            "roi": ["ROI_A", "ROI_A"],
            "direction": ["long", "long"],
            "b_step": [0, 1],
            "bvalue": [0.0, 1000.0],
            "value": [1.0, 0.7],
            "td_ms": [90.0, 90.0],
            "N": [4, 4],
            "Hz": [25.0, 25.0],
            "subj": ["S1", "S1"],
            "sheet": ["S1", "S1"],
            "source_file": ["legacy.long.parquet", "legacy.long.parquet"],
        }
    )


class FittingCommonTests(unittest.TestCase):
    def test_common_parameter_plan_parses_modes_values_and_bounds(self) -> None:
        args = argparse.Namespace(
            param_mode=["M0=fixed", "tc=global_td"],
            param_init=["D0=3.2e-12", "tc=7"],
            param_fixed=["M0=1"],
            param_bounds=["tc=0.2:500"],
        )
        plan = build_common_parameter_plan(
            args,
            param_names=("M0", "D0_m2_ms", "tc_ms"),
            default_modes={"M0": "free", "D0_m2_ms": "fixed", "tc_ms": "free"},
            default_inits={"M0": 0.9, "D0_m2_ms": 2.3e-12, "tc_ms": 5.0},
            default_bounds={"M0": (0.0, 5.0), "D0_m2_ms": (1e-16, 1e-10), "tc_ms": (0.1, 1000.0)},
            log_params=("D0_m2_ms", "tc_ms"),
        )

        self.assertEqual(plan.mode("M0"), "fixed")
        self.assertEqual(plan.fixed("M0"), 1.0)
        self.assertEqual(plan.mode("tc_ms"), "global_td")
        self.assertEqual(plan.init("D0_m2_ms", 0.0), 3.2e-12)
        self.assertEqual(plan.bounds("tc_ms", (0, 1)), (0.2, 500.0))
        self.assertEqual(plan.global_params(), ["tc_ms"])

    def test_load_master_input_selects_signal_rotated(self) -> None:
        master = append_master_rows(None, _signal_rows(), row_kind="signal_rotated", analysis_id="S1")
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "master.long.parquet"
            write_master_table(master, path)
            args = argparse.Namespace(
                master_parquet=path,
                row_kind="signal_rotated",
                analysis_id=None,
                subj=["S1"],
                sheet=None,
                roi=["ROI_A"],
                direction=["long"],
                stat="avg",
                source_file=None,
                td_ms=90.0,
                N=4.0,
                Hz=25.0,
            )
            selected = load_master_input(args, default_row_kind="signal_rotated", signal_rotated=True)

        self.assertEqual(len(selected), 2)
        self.assertEqual(selected["row_kind"].unique().tolist(), ["signal_rotated"])


if __name__ == "__main__":
    unittest.main()
