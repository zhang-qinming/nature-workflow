from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Callable


Action = Callable[[], None]


@dataclass
class Task:
    name: str
    preview: str
    command: list[str] | None = None
    cwd: Path | None = None
    action: Action | None = None

    def validate(self) -> None:
        if self.command is None and self.action is None:
            raise ValueError(f"Task '{self.name}' has neither command nor action")
        if self.command is not None and self.action is not None:
            raise ValueError(f"Task '{self.name}' cannot have both command and action")
