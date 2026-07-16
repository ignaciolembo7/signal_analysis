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

from data_processing.io import read_table_file  # noqa: E402


class SignalTableTests(unittest.TestCase):
    def test_read_table_file_reads_csv_tables(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "table.csv"
            pd.DataFrame({"roi": ["ROI"], "value": [1.0]}).to_csv(path, index=False)

            out = read_table_file(path)

        self.assertEqual(out.to_dict(orient="records"), [{"roi": "ROI", "value": 1.0}])


if __name__ == "__main__":
    unittest.main()
