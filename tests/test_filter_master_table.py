from __future__ import annotations

from pathlib import Path
import sys
import unittest

import pandas as pd

REPO_ROOT = Path(__file__).resolve().parents[1]
SRC_ROOT = REPO_ROOT / "src"
SCRIPT_ROOT = REPO_ROOT / "scripts" / "data"
for path in (SRC_ROOT, SCRIPT_ROOT):
    if str(path) not in sys.path:
        sys.path.insert(0, str(path))

from filter_master_table import filter_last_points, parse_last_points_spec  # noqa: E402


class ParseLastPointsSpecTests(unittest.TestCase):
    def test_td_only(self) -> None:
        parsed = parse_last_points_spec("120=8,210=6")
        self.assertEqual(parsed, {(120.0, None): 8, (210.0, None): 6})

    def test_td_and_n(self) -> None:
        parsed = parse_last_points_spec("120:8=6,120:4=4")
        self.assertEqual(parsed, {(120.0, 8.0): 6, (120.0, 4.0): 4})

    def test_mixed_and_all(self) -> None:
        parsed = parse_last_points_spec("120:8=6,210=8,90=ALL")
        self.assertEqual(parsed, {(120.0, 8.0): 6, (210.0, None): 8, (90.0, None): None})


class FilterLastPointsTests(unittest.TestCase):
    def _make_df(self, td_n_pairs: list[tuple[float, float]], b_steps: list[int]) -> pd.DataFrame:
        rows = []
        for td_ms, n_val in td_n_pairs:
            for b in b_steps:
                rows.append({"td_ms": td_ms, "N": n_val, "b_step": b, "value": 1.0})
        return pd.DataFrame(rows)

    def test_keeps_last_n_by_td(self) -> None:
        df = self._make_df([(120.0, 8.0), (120.0, 4.0)], [0, 1, 2, 3])
        result = filter_last_points(df, "120=2")
        # last 2 b_steps are 2 and 3
        self.assertEqual(sorted(result["b_step"].unique()), [2, 3])

    def test_keeps_unlisted_td(self) -> None:
        df = self._make_df([(120.0, 8.0), (90.0, 8.0)], [0, 1, 2, 3])
        result = filter_last_points(df, "120=2")
        td90 = result[result["td_ms"] == 90.0]
        self.assertEqual(sorted(td90["b_step"].tolist()), [0, 1, 2, 3])

    def test_td_and_n_rule(self) -> None:
        df = self._make_df([(120.0, 8.0), (120.0, 4.0)], [0, 1, 2, 3, 4, 5])
        # last 3 for N=8, last 2 for N=4
        result = filter_last_points(df, "120:8=3,120:4=2")
        n8 = sorted(result[result["N"] == 8.0]["b_step"].tolist())
        n4 = sorted(result[result["N"] == 4.0]["b_step"].tolist())
        self.assertEqual(n8, [3, 4, 5])
        self.assertEqual(n4, [4, 5])

    def test_specific_n_overrides_td_only(self) -> None:
        # TD-only rule keeps last 4; TD:N rule for N=8 keeps last 2
        df = self._make_df([(120.0, 8.0), (120.0, 4.0)], [0, 1, 2, 3, 4, 5])
        result = filter_last_points(df, "120=4,120:8=2")
        n8 = sorted(result[result["N"] == 8.0]["b_step"].tolist())
        n4 = sorted(result[result["N"] == 4.0]["b_step"].tolist())
        self.assertEqual(n8, [4, 5])   # overridden to last 2
        self.assertEqual(n4, [2, 3, 4, 5])  # TD-only keeps last 4

    def test_all_keyword_keeps_everything(self) -> None:
        df = self._make_df([(120.0, 8.0)], [0, 1, 2, 3, 4, 5])
        result = filter_last_points(df, "120=ALL")
        self.assertEqual(len(result), len(df))

    def test_empty_spec(self) -> None:
        df = self._make_df([(120.0, 8.0)], [0, 1, 2])
        result = filter_last_points(df, "")
        self.assertEqual(len(result), len(df))


if __name__ == "__main__":
    unittest.main()
