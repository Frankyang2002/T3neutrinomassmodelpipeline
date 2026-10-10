from pathlib import Path
from scripts.BuildCompactC5Report import parse_ordered, build_tex

def test_uploaded_matchete_c5():
    source = Path(__file__).resolve().parent / 'C5_original_matchete_sample.tex'
    if not source.exists():
        return
    rows, kernel = parse_ordered(source.read_text(encoding='utf-8'))
    assert len(rows) == 16
    assert kernel.count(r'\hbar') == 18
    assert r'C_{5,m}^{ij}' in build_tex(rows)
