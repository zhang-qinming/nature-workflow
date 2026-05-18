from __future__ import annotations

from pathlib import Path


def script_path(relative_path: str) -> Path:
    path = Path(__file__).resolve().parent / "scripts" / relative_path
    if not path.exists():
        raise FileNotFoundError(path)
    return path
