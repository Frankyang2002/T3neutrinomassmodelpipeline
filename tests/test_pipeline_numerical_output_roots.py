"""Prevent recurrence of the Numerical/output vs repository/output mistake."""
from pathlib import Path
from types import SimpleNamespace

from Numerical.orchestration import PipelineNumericalResults as pipeline


def test_repo_output_roots():
    repository = Path(__file__).resolve().parents[1]
    assert pipeline.PROJECT_ROOT.resolve() == repository.resolve()
    assert pipeline.RAW_OUTPUT_ROOT.resolve() == (repository / "output").resolve()
    assert pipeline.REPORT_OUTPUT_ROOT.resolve() == (repository / "Reports" / "output").resolve()


def test_report_figures_mirror_raw_tree():
    case = pipeline.RAW_OUTPUT_ROOT / "full" / "comparison" / "scenario" / "model" / "run"
    expected = pipeline.REPORT_OUTPUT_ROOT / "full" / "comparison" / "scenario" / "model" / "run" / "figures"
    assert pipeline._report_figure_dir(SimpleNamespace(output_dir=case)) == expected
