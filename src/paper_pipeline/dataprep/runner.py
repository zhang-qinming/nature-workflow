from __future__ import annotations

import subprocess
from pathlib import Path
from typing import Iterable

from .tasks import Task


def _cmdline(command: list[str]) -> str:
    return subprocess.list2cmdline([str(item) for item in command])


def print_tasks(tasks: Iterable[Task]) -> None:
    for index, task in enumerate(tasks, start=1):
        print(f"{index}. {task.name}")
        print(f"   {task.preview}")


def run_tasks(tasks: list[Task], dry_run: bool = False) -> None:
    for index, task in enumerate(tasks, start=1):
        task.validate()
        print(f"[{index}/{len(tasks)}] {task.name}")
        print(f"  {task.preview}")

        if dry_run:
            continue

        if task.command is not None:
            subprocess.run(
                [str(part) for part in task.command],
                cwd=str(task.cwd) if task.cwd else None,
                check=True,
            )
        elif task.action is not None:
            task.action()
        else:
            raise RuntimeError(f"Task '{task.name}' has no executable payload")


def preview_command(command: list[str], cwd: Path | None = None) -> str:
    if cwd is None:
        return _cmdline(command)
    return f"(cd {cwd} && {_cmdline(command)})"
