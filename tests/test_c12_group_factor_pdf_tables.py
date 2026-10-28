"""Current C12 PDF report table API regression tests.

The retired HTML/TeX parsing API (_extract_c12_sections and
transform_group_factor_tex) is not part of the current PDF generator.
"""
from __future__ import annotations

import pytest
import sympy as sp

from Reports import C12GroupFactorPDF as pdf


def _item(model, alpha, factor):
    return {'model': model, 'status': 'Projected',
            'representation': {'dS1': 2, 'dS2': 2, 'dF': 1, 'alpha': alpha},
            'terms': {'g2**2': sp.Integer(factor)}, 'total_collinear': True,
            'data_dir': None}


def test_current_table_grouping_keeps_representations_distinct():
    models = [_item('T3-B-m3', -3, 3), _item('T3-B-m1', -1, 5),
              _item('T3-B-p1', 1, 7)]
    blocks = pdf._run_blocks(models)
    assert len(blocks) == 1
    assert [x['model'] for x in blocks[0]] == [x['model'] for x in models]
    assert [x['terms']['g2**2'] for x in blocks[0]] == [3, 5, 7]


def test_current_table_generation_preserves_independent_factors():
    models = [_item('T3-B-m3', -3, 3), _item('T3-B-m1', -1, 5)]
    rendered = '\n'.join(pdf._comparison_table(models, 1))
    assert 'Run block 1' in rendered
    assert 'T3-B-m3' in rendered and 'T3-B-m1' in rendered
    assert '$3$' in rendered and '$5$' in rendered


def test_c12_pdf_section_replacement_preserves_existing_content():
    before = r'\documentclass{article}' + '\n' + r'\begin{document}' + '\n'
    before += 'Existing C5 content.\n' + r'\end{document}' + '\n'
    generated = pdf.BEGIN + '\n' + 'New C12 content.\n' + pdf.END + '\n'
    once = pdf.replace_pdf_tex_section(before, generated)
    twice = pdf.replace_pdf_tex_section(once, generated)
    assert once == twice
    assert 'Existing C5 content.' in twice
    assert twice.count('New C12 content.') == 1
    assert twice.strip().endswith(r'\end{document}')


def test_mismatched_section_markers_fail_closed():
    with pytest.raises(ValueError, match='markers'):
        pdf.replace_pdf_tex_section(pdf.BEGIN + r'\end{document}', 'replacement')
