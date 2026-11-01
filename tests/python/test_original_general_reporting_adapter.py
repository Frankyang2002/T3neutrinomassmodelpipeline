"""Functionality-only tests: original writers are used, not a new PDF layout."""
import ast
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]


def test_original_report_writers_are_called():
    src = (ROOT / 'Reports/GeneralSymbolicAdapter.py').read_text(encoding='utf-8')
    ast.parse(src)
    for name in ('write_bsm_uv_field_table(', 'write_rge_comparison(',
                 'write_eft1_rge_comparison(', '_saved_rge_beta_sections('):
        assert name in src
    assert 'WriteT3GeneralSymbolicPDFs' not in src


def test_presentation_module_has_short_quartic_labels():
    src = (ROOT / 'Reports/GeneralSymbolicPresentation.py').read_text(encoding='utf-8')
    assert r'\lambda_{X_{' in src
    assert 'GENERAL_QUARTIC' not in src
    ast.parse(src)


def test_legacy_conversion_module_has_no_report_layout():
    src = (ROOT / 'Reports/GeneralSymbolicAdapterLegacyConversions.py').read_text(encoding='utf-8')
    assert 'def matchete_terms(' in src
    assert 'def obtain_beta_latex(' in src
    assert 'def render(' not in src
