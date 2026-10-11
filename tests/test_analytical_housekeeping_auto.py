from pathlib import Path
import tempfile

from scripts.AnalyticalReportHousekeeping import clean, classify


def test_classification():
    assert classify(Path('x.aux')) == 'latex_build'
    assert classify(Path('x.tex.bak')) == 'backup'
    assert classify(Path('C5_original_matchete.pdf')) == 'superseded_presentation'
    assert classify(Path('C5.pdf')) is None
    assert classify(Path('GF_UV.tex')) is None


def test_auto_clean_preserves_reports_and_raw():
    with tempfile.TemporaryDirectory() as d:
        repo = Path(d)
        root = repo / 'Reports' / 'output' / 'analytical'
        root.mkdir(parents=True)
        for name in ('GF_UV.tex', 'GF_UV.pdf', 'GF_UV.aux', 'GF_UV.tex.bak', 'C5_original_matchete.pdf'):
            (root / name).write_text(name)
        raw = repo / 'output' / 'full'
        raw.mkdir(parents=True)
        (raw / 'data.json').write_text('important')
        r = clean(root)
        assert r['archived'] == 3
        assert (root / 'GF_UV.tex').exists()
        assert (root / 'GF_UV.pdf').exists()
        assert (raw / 'data.json').read_text() == 'important'
        assert (Path(r['archive']) / 'GF_UV.aux').exists()
