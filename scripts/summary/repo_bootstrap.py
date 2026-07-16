"""Import-path guard for scripts/summary command-line tools."""

from __future__ import annotations

import sys
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[2]
SRC_ROOT = REPO_ROOT / "src"

if SRC_ROOT.is_dir():
    src = str(SRC_ROOT)
    if src in sys.path:
        sys.path.remove(src)
    sys.path.insert(0, src)
