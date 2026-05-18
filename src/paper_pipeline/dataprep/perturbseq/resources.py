from __future__ import annotations

from pathlib import Path


def r_script_path(name: str) -> Path:
    path = Path(__file__).resolve().parent / "r_scripts" / name
    if not path.exists():
        raise FileNotFoundError(path)
    return path
