"""Within-model dashboard using the *same renderer* as cross-model comparison.

The four scenario diagnostics are staged under scenario names and then passed to
InteractiveModelComparison.write_dashboard. This is presentation-only; no
matching, RG evolution or physical quantities are recalculated. The temporary
staging directory is deleted after the HTML has been generated.
"""
from __future__ import annotations

from copy import deepcopy
import json
from pathlib import Path
from tempfile import TemporaryDirectory
from typing import Mapping

from Numerical.plotting.InteractiveModelComparison import collect, write_dashboard
from Numerical.plotting.RunningTrajectoryData import SCENARIOS, SavedRun

LABELS = {
    "smallY_smallL": "Y=0.005, lambda(ref)=0.1",
    "smallY_largeL": "Y=0.005, lambda(ref)=1.0",
    "largeY_smallL": "Y=0.5, lambda(ref)=0.1",
    "largeY_largeL": "Y=0.5, lambda(ref)=1.0",
}


def _stage_scenarios(runs: Mapping[str, SavedRun], root: Path) -> None:
    """Supply cross-model `collect` with scenario-labelled diagnostics.

    Only the diagnostic's *display identity* is changed, never its physics
    arrays, matching contributions, scales, or benchmark metadata.
    """
    for scenario in SCENARIOS:
        if scenario not in runs:
            continue
        run = runs[scenario]
        payload = deepcopy(run.payload)
        search = payload.get("benchmark_search")
        if not isinstance(search, dict):
            search = {}
            payload["benchmark_search"] = search
        model = search.get("model")
        if not isinstance(model, dict):
            model = {}
            search["model"] = model
        # InteractiveModelComparison._model_name prioritises this field.
        model["model_key"] = scenario
        path = root / scenario / "data" / "running_diagnostics.json"
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(json.dumps(payload, allow_nan=False), encoding="utf-8")


def build_interactive_data(runs: Mapping[str, SavedRun]) -> dict:
    """Return the exact structured input used by the shared dashboard."""
    if not runs:
        raise ValueError("At least one saved comparison scenario is required")
    with TemporaryDirectory(prefix="t3_within_dashboard_") as dirname:
        root = Path(dirname)
        _stage_scenarios(runs, root)
        return collect(root)


def write_interactive_comparison(
    runs: Mapping[str, SavedRun], output_path: Path
) -> tuple[Path, list[str]]:
    """Create a within-model HTML using the canonical cross-model renderer."""
    if not runs:
        raise ValueError("At least one saved comparison scenario is required")
    model_key = next(iter(runs.values())).model_key
    if any(run.model_key != model_key for run in runs.values()):
        raise ValueError("A within-model comparison cannot mix different models")
    unknown = set(runs).difference(SCENARIOS)
    if unknown:
        raise ValueError(f"Unexpected scenarios: {sorted(unknown)}")

    output_path = Path(output_path)
    with TemporaryDirectory(prefix="t3_within_dashboard_") as dirname:
        root = Path(dirname)
        _stage_scenarios(runs, root)
        # Reuse ALL of the cross-model front-end: plotting canvas, EFT bands,
        # scale zoom presets, component overlays, axes and model toggles.
        write_dashboard(
            root,
            output_path,
            f"{model_key} — four-scenario within-model comparison",
        )
        data = collect(root)

    # Change only the displayed grouping nouns. The rendering code is identical.
    page = output_path.read_text(encoding="utf-8")
    page = page.replace("Show all models", "Show all scenarios")
    page = page.replace("Hide all models", "Hide all scenarios")
    page = page.replace('<label class="caption">Models</label>',
                        '<label class="caption">Scenarios</label>')
    page = page.replace("visible models / current zoom", "visible scenarios / current zoom")
    page = page.replace("all models / all scales", "all scenarios / all scales")
    page = page.replace("when models are shown or hidden", "when scenarios are shown or hidden")
    output_path.write_text(page, encoding="utf-8")

    warnings = [f"{scenario}: {reason}" for scenario, reason in
                sorted(data["continuation_issues"].items())]
    return output_path, warnings
