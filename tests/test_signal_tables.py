from __future__ import annotations

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

from data_processing.io import read_table_file  # noqa: E402
from data_processing.match_params import parse_results_filename  # noqa: E402
from process_one_results import resolve_table_layout  # noqa: E402


class SignalTableTests(unittest.TestCase):
    def test_read_table_file_reads_csv_tables(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "table.csv"
            pd.DataFrame({"roi": ["ROI"], "value": [1.0]}).to_csv(path, index=False)

            out = read_table_file(path)

        self.assertEqual(out.to_dict(orient="records"), [{"roi": "ROI", "value": 1.0}])

    def test_bids_results_filename_metadata_is_parsed(self) -> None:
        meta = parse_results_filename("sub-BRAIN-1_ses-T0_acq-hz025d40b1075s6_results.xlsx")

        self.assertEqual(meta.subject_id, "BRAIN-1")
        self.assertEqual(meta.seq, 6)
        self.assertEqual(meta.Hz, 25.0)
        self.assertEqual(meta.bmax, 1075.0)
        self.assertEqual(meta.d_ms, 40.0)

    def test_layout_can_come_from_sequence_parameters(self) -> None:
        stats = {"avg": pd.DataFrame({"bvalues": [0, 0] + list(range(12)), "roi": np.arange(14)})}
        params_row = pd.Series({"ndirs": 3, "nbvals": 4})

        ndirs, nbvals = resolve_table_layout(
            stats=stats,
            results_file=Path("sub-TEST_ses-T0_acq-hz025d40b100s1_results.xlsx"),
            params_row=params_row,
            gradient_col="bvalues",
        )

        self.assertEqual((ndirs, nbvals), (3, 4))

    def test_bids_signal_extraction_table_shape_falls_back_to_six_by_ten(self) -> None:
        bvalues = [0, 0] + [0] * 6 + [5] * 12 + [10] * 6 + [15] * 6 + [20] * 6
        bvalues += [25] * 6 + [30] * 6 + [40] * 6 + [45] * 6
        stats = {"avg": pd.DataFrame({"bvalues": bvalues, "roi": np.arange(len(bvalues))})}

        ndirs, nbvals = resolve_table_layout(
            stats=stats,
            results_file=Path("sub-BRAIN-1_ses-T0_acq-hz100d40b0045s9_results.xlsx"),
            params_row=None,
            gradient_col="bvalues",
        )

        self.assertEqual((ndirs, nbvals), (6, 10))


if __name__ == "__main__":
    unittest.main()
