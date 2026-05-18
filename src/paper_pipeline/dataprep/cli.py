from __future__ import annotations

import argparse
import subprocess
import sys

from ..workflows import TOP_LEVEL_WORKFLOWS, WORKFLOW_BUILDERS
from .config import ConfigError, load_config
from .runner import print_tasks, run_tasks


ALL_WORKFLOW_NAME = "all"


def _workflow_choices() -> list[str]:
    return sorted(list(WORKFLOW_BUILDERS) + [ALL_WORKFLOW_NAME])


def _workflow_help() -> str:
    top_level = ", ".join(TOP_LEVEL_WORKFLOWS)
    return (
        "Workflow name to plan or run. "
        f"`{ALL_WORKFLOW_NAME}` expands only the configured top-level aggregate workflows: {top_level}."
    )


def _collect_tasks(config_path: str, workflow_name: str):
    config = load_config(config_path)
    if workflow_name == ALL_WORKFLOW_NAME:
        tasks = []
        configured_workflows = config.raw.get("workflows", {})
        if not isinstance(configured_workflows, dict):
            raise ConfigError(f"Config key `workflows` must be a mapping in {config.path}")
        configured_top_level = []
        for builder_name in TOP_LEVEL_WORKFLOWS:
            if builder_name not in configured_workflows:
                continue
            configured_top_level.append(builder_name)
            tasks.extend(WORKFLOW_BUILDERS[builder_name](config))
        if not configured_top_level:
            top_level = ", ".join(TOP_LEVEL_WORKFLOWS)
            raise ConfigError(
                f"`{ALL_WORKFLOW_NAME}` requires at least one configured top-level workflow in "
                f"`workflows`: {top_level}"
            )
        return tasks
    return WORKFLOW_BUILDERS[workflow_name](config)


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Config-driven entrypoint for paper-pipeline workflows."
    )
    subparsers = parser.add_subparsers(dest="command", required=True)

    plan_parser = subparsers.add_parser(
        "plan",
        help="Print the configured workflow plan",
        description=(
            "Print the configured workflow plan. "
            f"`{ALL_WORKFLOW_NAME}` includes only configured top-level aggregate workflows."
        ),
    )
    plan_parser.add_argument("--config", required=True, help="Path to YAML config")
    plan_parser.add_argument("workflow", choices=_workflow_choices(), help=_workflow_help())

    run_parser = subparsers.add_parser(
        "run",
        help="Run a configured workflow",
        description=(
            "Run a configured workflow. "
            f"`{ALL_WORKFLOW_NAME}` includes only configured top-level aggregate workflows."
        ),
    )
    run_parser.add_argument("--config", required=True, help="Path to YAML config")
    run_parser.add_argument("workflow", choices=_workflow_choices(), help=_workflow_help())
    run_parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Print the steps without executing them",
    )

    args = parser.parse_args()

    try:
        tasks = _collect_tasks(args.config, args.workflow)
        if args.command == "plan":
            print_tasks(tasks)
        else:
            run_tasks(tasks, dry_run=bool(args.dry_run))
    except ConfigError as exc:
        print(f"Configuration error: {exc}", file=sys.stderr)
        raise SystemExit(1) from exc
    except ValueError as exc:
        print(f"Task validation error: {exc}", file=sys.stderr)
        raise SystemExit(1) from exc
    except subprocess.CalledProcessError as exc:
        print(f"Command failed (exit code {exc.returncode}): {exc}", file=sys.stderr)
        raise SystemExit(exc.returncode) from exc


if __name__ == "__main__":
    main()
