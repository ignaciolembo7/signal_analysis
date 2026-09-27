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
from data_processing.match_params import parse_results_filename, select_params_row  # noqa: E402
from data_processing.result_signals import CleanSequenceParams, grouped_direct_g_output_stem  # noqa: E402
from process_one_results import _drop_existing_direct_g_curve, resolve_table_layout  # noqa: E402


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

    def test_bids_nogse_filename_parses_type_and_sequence_before_dwi(self) -> None:
        meta = parse_results_filename(
            "sub-PHANTOM-1_ses-20260505_acq-002_NOGSE_CPMG_N2_TN50_G04_3_dwi_results.xlsx"
        )

        self.assertEqual(meta.encoding, "NOGSE")
        self.assertEqual(meta.sequence_type, "CPMG")
        self.assertEqual(meta.group, 2)
        self.assertEqual(meta.G, 4.0)
        self.assertEqual(meta.TN, 50.0)
        self.assertEqual(meta.N, 2.0)
        self.assertEqual(meta.seq, 3)

    def test_nogse_type_disambiguates_cpmg_and_hahn_parameter_rows(self) -> None:
        meta = parse_results_filename(
            "sub-PHANTOM-1_ses-20260505_acq-002_NOGSE_CPMG_N2_TN50_G04_dwi_results.xlsx"
        )
        params = pd.DataFrame(
            {
                "sheet": ["20260505-PHANTOM_FIBER", "20260505-PHANTOM_FIBER"],
                "subj": ["PHANTOM-1", "PHANTOM-1"],
                "group": [2, 2],
                "G": [4, 4],
                "TN": [50, 50],
                "N": [2, 2],
                "type": ["CPMG", "HAHN"],
            }
        )

        selected = select_params_row(params, meta)

        self.assertIsNotNone(selected)
        self.assertEqual(selected["type"], "CPMG")

    def test_direct_g_curve_stem_excludes_g_and_sequence(self) -> None:
        common = dict(
            sheet="20260505-PHANTOM_FIBER",
            subj="PHANTOM-1",
            protocol="NOGSE_CPMG_group002_TN50",
            group=2,
            TN=50.0,
            x=np.nan,
            y=np.nan,
            type="CPMG",
            Hz=np.nan,
            bmax=np.nan,
            N=2,
            delta_ms=10.0,
            Delta_app_ms=50.0,
            max_dur_ms=np.nan,
            tm_ms=np.nan,
            td_ms=50.0,
            TE=100.0,
            TR=2000.0,
            g_thorsten=None,
        )
        g0 = CleanSequenceParams(sequence=2, G=0.0, **common)
        g4 = CleanSequenceParams(sequence=3, G=4.0, **common)

        self.assertEqual(grouped_direct_g_output_stem(g0), grouped_direct_g_output_stem(g4))
        self.assertNotIn("G-", grouped_direct_g_output_stem(g0))

    def test_reingest_replaces_only_the_matching_direct_g_curve(self) -> None:
        common = {
            "row_kind": "signal",
            "sheet": "20260505-PHANTOM_FIBER",
            "subj": "PHANTOM-1",
            "protocol": "NOGSE_group002_TN50",
            "group": 2,
            "TN": 50,
            "N": 2,
        }
        master = pd.DataFrame(
            [
                {**common, "type": "CPMG", "G": 0, "value": 100},
                {**common, "type": "HAHN", "G": 0, "value": 90},
            ]
        )
        incoming = pd.DataFrame(
            [
                {**common, "type": "CPMG", "G": 0, "value": 100},
                {**common, "type": "CPMG", "G": 4, "value": 95},
            ]
        )

        remaining = _drop_existing_direct_g_curve(master, incoming)

        self.assertEqual(len(remaining), 1)
        self.assertEqual(remaining.iloc[0]["type"], "HAHN")

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
