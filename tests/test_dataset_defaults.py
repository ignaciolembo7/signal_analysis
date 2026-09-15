from __future__ import annotations

import os
from pathlib import Path
import subprocess
import unittest


REPO_ROOT = Path(__file__).resolve().parents[1]
HELPER = REPO_ROOT / "helpers" / "master_table_common.sh"


class DatasetDefaultTests(unittest.TestCase):
    def _alpha_reference_defaults(self, dataset: str) -> tuple[str, str]:
        environment = os.environ.copy()
        environment.pop("ALPHA_REFERENCE_D0_MM2_S", None)
        environment.pop("ALPHA_REFERENCE_D0_ERROR_MM2_S", None)
        command = (
            f'source "{HELPER}"; '
            "pipeline_setup_common; "
            f'pipeline_set_dataset_defaults "{dataset}" ogse; '
            'printf "%s|%s" "$ALPHA_REFERENCE_D0_MM2_S" "$ALPHA_REFERENCE_D0_ERROR_MM2_S"'
        )
        completed = subprocess.run(
            ["bash", "-c", command],
            cwd=REPO_ROOT,
            env=environment,
            check=True,
            capture_output=True,
            text=True,
        )
        return tuple(completed.stdout.split("|", maxsplit=1))

    def test_brain_alpha_reference_defaults(self) -> None:
        self.assertEqual(self._alpha_reference_defaults("brain"), ("0.0032", "0.0000283512"))

    def test_phantom_alpha_reference_defaults(self) -> None:
        self.assertEqual(self._alpha_reference_defaults("phantom"), ("0.0023", "0.0"))


if __name__ == "__main__":
    unittest.main()
