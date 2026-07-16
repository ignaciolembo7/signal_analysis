from __future__ import annotations

from argparse import Namespace
from pathlib import Path
import tempfile
import unittest

import pandas as pd

from pipeline.recipe import selected_rows_or_legacy_dataframe, selected_rows_or_legacy_table


class PipelineRecipeTests(unittest.TestCase):
    def test_legacy_dataframe_selection_reads_one_table(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "signal.long.parquet"
            pd.DataFrame({"value": [1.0, 2.0]}).to_parquet(path, index=False)
            args = Namespace(master_parquet=None)

            selected = selected_rows_or_legacy_dataframe(
                args,
                legacy_path=path,
                default_row_kind="signal",
            )

            self.assertEqual(selected.source, "legacy")
            self.assertEqual(selected.paths, [path])
            self.assertEqual(selected.df["value"].tolist(), [1.0, 2.0])

    def test_legacy_path_selection_keeps_paths_without_temp_files(self) -> None:
        path = Path("signal.long.parquet")
        args = Namespace(master_parquet=None)

        selected = selected_rows_or_legacy_table(
            args,
            legacy_path=path,
            default_row_kind="signal",
            temp_prefix="recipe_test_",
        )

        self.assertEqual(selected.source, "legacy")
        self.assertEqual(selected.paths, [path])
        self.assertIsNone(selected.df)


if __name__ == "__main__":
    unittest.main()
