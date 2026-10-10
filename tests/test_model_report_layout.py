"""Isolated filesystem tests; never call matching or numerical evolution."""
from pathlib import Path

from Numerical.orchestration.ModelReportLayout import migrate

MODEL = 'T3_dS1_3_dS2_3_dF_2_alpha_p0'


def _report(root: Path, scenario: str, name: str, content: bytes) -> Path:
    path = root / 'full' / 'comparison' / scenario / MODEL / 'RGE' / name
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_bytes(content)
    return path


def test_dry_run_and_copy(tmp_path: Path) -> None:
    a = _report(tmp_path, 'smallY_smallL', 'UV.pdf', b'valid PDF placeholder')
    b = _report(tmp_path, 'largeY_largeL', 'UV.pdf', a.read_bytes())
    assert migrate(tmp_path)['unique_files'] == 1
    assert not (tmp_path / 'models').exists()
    result = migrate(tmp_path, execute=True)
    assert result['status'] == 'Complete' and result['copied'] == 1
    assert a.exists() and b.exists()
    assert (tmp_path / 'models' / MODEL / 'RGE' / 'UV.pdf').read_bytes() == a.read_bytes()
    assert (tmp_path / 'models' / 'report_migration_manifest.json').exists()


def test_conflicts_block_all_writes(tmp_path: Path) -> None:
    _report(tmp_path, 'smallY_smallL', 'UV.pdf', b'A')
    _report(tmp_path, 'largeY_largeL', 'UV.pdf', b'B')
    result = migrate(tmp_path, execute=True)
    assert result['status'] == 'Conflict' and result['conflicts']
    assert not (tmp_path / 'models').exists()


def test_prune_verified_duplicates_only(tmp_path: Path) -> None:
    a = _report(tmp_path, 'smallY_smallL', 'UV.tex', b'symbolic')
    b = _report(tmp_path, 'smallY_largeL', 'UV.tex', b'symbolic')
    migrate(tmp_path, execute=True, prune=True)
    assert not a.exists() and not b.exists()
    assert (tmp_path / 'models' / MODEL / 'RGE' / 'UV.tex').read_bytes() == b'symbolic'


def test_ignores_numerical_figures(tmp_path: Path) -> None:
    _report(tmp_path, 'smallY_smallL', 'UV.tex', b'model')
    figure = tmp_path / 'full' / 'comparison' / 'smallY_smallL' / MODEL / 'figures' / 'plot.png'
    figure.parent.mkdir(parents=True, exist_ok=True)
    figure.write_bytes(b'numerical')
    migrate(tmp_path, execute=True, prune=True)
    assert figure.read_bytes() == b'numerical'
