from __future__ import annotations
from pathlib import Path
import numpy as np


def _normalize_rows(v: np.ndarray) -> np.ndarray:
    n = np.linalg.norm(v, axis=1, keepdims=True)
    n[n == 0] = 1.0
    return v / n


def repo_root_from_src_layout() -> Path:
    return Path(__file__).resolve().parents[2]


def _resolve_dirs_path(path: str | Path) -> Path:
    p = Path(path)
    suffix = p.suffix.lower()

    if suffix == ".csv":
        fallback = p.with_suffix(".txt")
        if fallback.is_file():
            return fallback
        raise ValueError(
            f"Only .txt direction files are supported now. Got: {p}."
        )

    if suffix != ".txt":
        raise ValueError(f"Unsupported directions file format: {p}. Expected .txt")

    return p


def load_dirs_txt(path: str | Path) -> np.ndarray:
    p = _resolve_dirs_path(path)
    arr = np.loadtxt(p, dtype=float)
    if arr.ndim == 1:
        if arr.size == 3:
            arr = arr.reshape(1, 3)
        else:
            raise ValueError(f"Invalid directions table: expected Nx3, got shape {arr.shape}")
    if arr.ndim != 2 or arr.shape[1] != 3:
        raise ValueError(f"Invalid directions table: expected Nx3, got shape {arr.shape}")
    return _normalize_rows(arr)


def load_dirs_csv(path: str | Path) -> np.ndarray:
    return load_dirs_txt(path)


def load_default_dirs(ndirs: int) -> np.ndarray:
    p = repo_root_from_src_layout() / "assets" / "dirs" / f"dirs_{ndirs}.txt"
    if not p.is_file():
        raise FileNotFoundError(
            f"Missing {p}. For ndirs={ndirs}, create this .txt table "
            f"(one row per direction, in direction=1..{ndirs} order)."
        )
    dirs = load_dirs_txt(p)
    if dirs.shape[0] != ndirs:
        raise ValueError(f"{p} tiene {dirs.shape[0]} filas, pero ndirs={ndirs}.")
    return dirs
