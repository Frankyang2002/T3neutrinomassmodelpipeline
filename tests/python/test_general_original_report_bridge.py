from pathlib import Path

def test_main_pipeline_routes_to_original_report_adapter():
    root = Path(__file__).resolve().parents[2]
    p = (root / 'pipeline.py').read_text(encoding='utf-8')
    assert 'from Reports.GeneralSymbolicAdapter import render' in p
    assert 'from scripts.WriteT3GeneralSymbolicPDFs import' not in p

def test_adapter_reuses_original_report_modules():
    root = Path(__file__).resolve().parents[2]
    text = (root / 'Reports/GeneralSymbolicAdapter.py').read_text(encoding='utf-8')
    assert 'from Reports.ReportGeneration import compile_latex_document' in text
    assert 'from Reports.ReportGeneration import paper_notation_key_lines' in text
    assert 'from Reports.GroupFactorReports import _saved_rge_term_table' in text
    assert 'from Reports.RGEComparison import normalise_rgbeta_latex' in text
    assert 'report_betas' in text and 'bsm_lagrangian_inputform_file' in text
