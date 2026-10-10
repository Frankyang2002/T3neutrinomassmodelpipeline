from pathlib import Path
from Numerical.orchestration.ModelReportLayout import migrate

MODEL = "T3_dS1_2_dS2_2_dF_1_alpha_m1"


def report(root: Path, scenario: str, name: str, content: bytes) -> Path:
    path = root / "full" / "comparison" / scenario / MODEL / "RGE" / name
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_bytes(content)
    return path


def test_dry_run_does_not_move(tmp_path):
    src = report(tmp_path, "smallY_smallL", "UV.pdf", b"a")
    summary = migrate(tmp_path)
    assert summary["source_files"] == 1 and src.exists()
    assert not (tmp_path / "models").exists()


def test_identical_reports_move_to_one_destination(tmp_path):
    a = report(tmp_path, "smallY_smallL", "UV.pdf", b"a")
    b = report(tmp_path, "largeY_largeL", "UV.pdf", b"a")
    summary = migrate(tmp_path, execute=True)
    assert summary["status"] == "Complete"
    assert summary["removed_source_files"] == 2
    assert not a.exists() and not b.exists()
    assert (tmp_path / "models" / MODEL / "RGE" / "UV.pdf").read_bytes() == b"a"
    assert (tmp_path / "models" / "report_migration_manifest.json").exists()


def test_differing_reports_preserved_as_variant(tmp_path):
    report(tmp_path, "smallY_smallL", "UV.pdf", b"canonical")
    report(tmp_path, "largeY_largeL", "UV.pdf", b"variant")
    summary = migrate(tmp_path, execute=True)
    assert len(summary["different_content_variants"]) == 1
    assert (tmp_path / "models" / MODEL / "RGE" / "UV.pdf").read_bytes() == b"canonical"
    assert (tmp_path / "models" / MODEL / "_scenario_variants" / "largeY_largeL" / "RGE" / "UV.pdf").read_bytes() == b"variant"


def test_numerical_figure_untouched(tmp_path):
    report(tmp_path, "smallY_smallL", "UV.tex", b"text")
    figure = tmp_path / "full" / "comparison" / "smallY_smallL" / MODEL / "figures" / "plot.png"
    figure.parent.mkdir(parents=True, exist_ok=True)
    figure.write_bytes(b"png")
    migrate(tmp_path, execute=True)
    assert figure.read_bytes() == b"png"
