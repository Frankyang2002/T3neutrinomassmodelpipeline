"""Persistent numerical-config and benchmark bookkeeping for pipeline runs.

This module owns the JSON/file-management details associated with model-specific
numerical configurations.  The central ``pipeline.py`` only decides when these
operations occur.
"""

from __future__ import annotations

from copy import deepcopy
import json
from pathlib import Path
from typing import Any

from common.RunRecords import RunRecord
from Numerical.SobolBenchmarkSearch import (
    _parameters_from_config,
    apply_benchmark_parameters,
    model_benchmark_key,
    retarget_payload_to_record,
)


PROJECT_ROOT = Path(__file__).resolve().parents[1]
GENERATED_NUMERICAL_CONFIG_DIR = PROJECT_ROOT / "configs" / "generated_models"


def _load_json_object(path: Path) -> dict[str, Any]:
    payload = json.loads(Path(path).read_text(encoding="utf-8"))
    if not isinstance(payload, dict):
        raise ValueError(f"{path} must contain a JSON object.")
    return payload


def _atomic_write_json(path: Path, payload: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + ".tmp")
    temporary.write_text(
        json.dumps(payload, indent=2) + "\n",
        encoding="utf-8",
    )
    temporary.replace(path)


def _generated_numerical_config_path(record: RunRecord) -> Path:
    return GENERATED_NUMERICAL_CONFIG_DIR / f"{model_benchmark_key(record)}.json"


def prepare_model_numerical_config(
    record: RunRecord,
    template_path: Path,
    *,
    reset: bool,
) -> Path:
    """Return the persistent model-specific numerical config for one record.

    Existing generated configs are reused by default.  With ``reset=True`` the
    supplied numerical template is retargeted to the current model and replaces
    the generated config before the numerical run.
    """

    generated_path = _generated_numerical_config_path(record)

    if generated_path.is_file() and not reset:
        print(
            f"  {record.name}: using saved numerical config -> {generated_path}",
            flush=True,
        )
        return generated_path

    template = _load_json_object(Path(template_path))
    retargeted = retarget_payload_to_record(record, template)
    retargeted["generated_config"] = {
        "model_key": model_benchmark_key(record),
        "source_template": str(Path(template_path).resolve()),
        "resolved_benchmark": False,
    }
    _atomic_write_json(generated_path, retargeted)
    print(
        f"  {record.name}: saved model numerical config -> {generated_path}",
        flush=True,
    )
    return generated_path


def freeze_resolved_benchmark_config(
    record: RunRecord,
    generated_path: Path,
) -> None:
    """Freeze the benchmark actually used into the persistent model config.

    The integrated numerical run writes benchmark provenance into
    ``running_diagnostics.json``.  When that provenance points to a Sobol cache,
    apply the cache's best parameters to the generated config and disable
    ``benchmark_search`` so future runs reproduce that exact point.
    """

    diagnostics_relative = record.summary.get("RunningDiagnosticsFile")
    if not diagnostics_relative:
        return

    diagnostics_path = record.output_dir / str(diagnostics_relative)
    if not diagnostics_path.is_file():
        return

    diagnostics = _load_json_object(diagnostics_path)
    metadata = diagnostics.get("benchmark_search")
    if not isinstance(metadata, dict):
        return

    cache_value = metadata.get("cache_path")
    if not cache_value:
        return

    cache_path = Path(str(cache_value))
    if not cache_path.is_absolute():
        cache_path = PROJECT_ROOT / cache_path
    if not cache_path.is_file():
        return

    cache = _load_json_object(cache_path)
    best = cache.get("best")
    if not isinstance(best, dict) or not isinstance(best.get("parameters"), dict):
        return

    payload = _load_json_object(generated_path)
    definitions = _parameters_from_config(payload)
    resolved = apply_benchmark_parameters(
        payload,
        best["parameters"],
        definitions,
    )

    search = resolved.get("benchmark_search")
    if isinstance(search, dict):
        search = deepcopy(search)
        search["enabled"] = False
        resolved["benchmark_search"] = search

    resolved["generated_config"] = {
        "model_key": model_benchmark_key(record),
        "source_template": payload.get("generated_config", {}).get(
            "source_template"
        ),
        "resolved_benchmark": True,
        "benchmark_cache": str(cache_path.resolve()),
        "benchmark_status": metadata.get("benchmark_status"),
        "benchmark_chi2": metadata.get("chi2"),
        "benchmark_method": metadata.get("method"),
    }
    _atomic_write_json(generated_path, resolved)
    print(
        f"  {record.name}: froze resolved benchmark config -> {generated_path}",
        flush=True,
    )
