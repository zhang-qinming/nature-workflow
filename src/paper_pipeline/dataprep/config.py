from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Any

try:
    import yaml
except ModuleNotFoundError:  # pragma: no cover - depends on local environment
    yaml = None


class ConfigError(ValueError):
    pass


def _resolve_path(base_dir: Path, value: str | Path | None) -> Path | None:
    if value is None:
        return None
    path = Path(value)
    if path.is_absolute():
        return path
    return (base_dir / path).resolve()


@dataclass
class LoadedConfig:
    path: Path
    raw: dict[str, Any]

    @property
    def base_dir(self) -> Path:
        return self.path.parent.resolve()

    @property
    def project_root(self) -> Path:
        return _resolve_path(self.base_dir, self.raw.get("project_root", "."))  # type: ignore[return-value]

    @property
    def artifact_root(self) -> Path:
        raw_root = (
            self.raw.get("artifact_root")
            or self.raw.get("output_root")
            or self.raw.get("outputs_root")
            or "outputs"
        )
        return _resolve_path(self.base_dir, raw_root)  # type: ignore[return-value]

    def resolve_path(self, value: str | Path | None) -> Path | None:
        return _resolve_path(self.base_dir, value)

    def artifact_path(self, *parts: str | Path) -> Path:
        path = self.artifact_root
        for part in parts:
            path = path / Path(part)
        return path.resolve()

    def resolve_path_or_artifact(
        self,
        value: str | Path | None,
        *artifact_parts: str | Path,
    ) -> Path:
        resolved = self.resolve_path(value)
        if resolved is not None:
            return resolved
        if not artifact_parts:
            return self.artifact_root.resolve()
        return self.artifact_path(*artifact_parts)

    def executable_command(
        self,
        name: str,
        default: str | list[str] | tuple[str, ...] | None = None,
    ) -> list[str]:
        raw = self.raw.get("executables", {}).get(name, default)
        if raw is None:
            return []
        if isinstance(raw, str):
            return [raw]
        if isinstance(raw, (list, tuple)) and all(isinstance(value, str) for value in raw):
            return [str(value) for value in raw]
        raise ConfigError(
            f"executables.{name} must be a string or list of strings, got {type(raw).__name__}"
        )

    def executable(self, name: str, default: str) -> str:
        command = self.executable_command(name, default)
        if len(command) != 1:
            raise ConfigError(
                f"executables.{name} must resolve to a single executable string for this call site"
            )
        return command[0]

    def workflow(self, name: str) -> dict[str, Any]:
        workflows = self.raw.get("workflows", {})
        if name not in workflows:
            raise ConfigError(f"Workflow '{name}' is not defined in {self.path}")
        value = workflows[name]
        if not isinstance(value, dict):
            raise ConfigError(f"Workflow '{name}' must be a mapping")
        return value


def load_config(path: str | Path) -> LoadedConfig:
    if yaml is None:
        raise ConfigError(
            "PyYAML is required to load the paper-pipeline config. "
            "Install project dependencies with `pip install -e .`."
        )

    config_path = Path(path).resolve()
    if not config_path.exists():
        raise ConfigError(f"Config file does not exist: {config_path}")

    data = yaml.safe_load(config_path.read_text(encoding="utf-8")) or {}
    if not isinstance(data, dict):
        raise ConfigError(f"Config root must be a mapping: {config_path}")

    return LoadedConfig(path=config_path, raw=data)
