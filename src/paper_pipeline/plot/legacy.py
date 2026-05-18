from __future__ import annotations

import os
import shutil
import subprocess
from pathlib import Path
from typing import Iterable, Sequence

from ..dataprep.config import ConfigError
from ..dataprep.runner import preview_command


def ensure_directories(paths: Iterable[Path]) -> None:
    for path in paths:
        path.mkdir(parents=True, exist_ok=True)


def _same_target(path: Path, source: Path) -> bool:
    if not path.exists() and not path.is_symlink():
        return False
    try:
        return path.resolve() == source.resolve()
    except OSError:
        return False


def link_path(source: Path, target: Path) -> None:
    if not source.exists():
        raise ConfigError(f"Legacy plotting input does not exist: {source}")

    ensure_directories([target.parent])
    if _same_target(target, source):
        return
    if target.is_symlink() or target.is_file():
        target.unlink()
    elif target.is_dir():
        return

    try:
        target.symlink_to(source, target_is_directory=source.is_dir())
    except OSError:
        if source.is_dir():
            shutil.copytree(source, target, dirs_exist_ok=True)
        else:
            shutil.copy2(source, target)


def populate_directory(
    source_dir: Path,
    target_dir: Path,
    *,
    exclude_names: Sequence[str] = (),
    include_names: Sequence[str] | None = None,
) -> None:
    if not source_dir.exists():
        raise ConfigError(f"Legacy plotting input directory does not exist: {source_dir}")
    if not source_dir.is_dir():
        raise ConfigError(f"Expected directory for legacy plotting input: {source_dir}")

    ensure_directories([target_dir])
    excluded = set(exclude_names)
    included = set(include_names) if include_names is not None else None

    for child in source_dir.iterdir():
        if child.name in excluded:
            continue
        if included is not None and child.name not in included:
            continue
        link_path(child, target_dir / child.name)


def write_text(path: Path, content: str) -> None:
    ensure_directories([path.parent])
    path.write_text(content, encoding="utf-8")


def run_command(command: Sequence[str], *, cwd: Path, env: dict[str, str] | None = None) -> None:
    subprocess.run(
        [str(part) for part in command],
        cwd=str(cwd),
        check=True,
        env=None if env is None else {**os.environ, **env},
    )


def run_commands(commands: Iterable[Sequence[str]], *, cwd: Path, env: dict[str, str] | None = None) -> None:
    for command in commands:
        run_command(command, cwd=cwd, env=env)


def format_command_preview(command: Sequence[str], *, cwd: Path) -> str:
    return preview_command([str(part) for part in command], cwd)


def link_alias(source: Path, target: Path) -> None:
    link_path(source, target)


def link_output_path(target: Path, alias: Path, *, is_directory: bool) -> None:
    ensure_directories([alias.parent])
    if is_directory:
        ensure_directories([target])
    else:
        ensure_directories([target.parent])
        target.touch(exist_ok=True)

    if _same_target(alias, target):
        return
    replace_tree(alias)

    try:
        alias.symlink_to(target, target_is_directory=is_directory)
    except OSError:
        if is_directory:
            ensure_directories([alias])
        else:
            alias.write_text(target.read_text(encoding="utf-8"), encoding="utf-8")


def replace_tree(path: Path) -> None:
    if path.is_symlink() or path.is_file():
        path.unlink()
    elif path.is_dir():
        shutil.rmtree(path)
