"""Adapt saved general-interaction outputs to the ORIGINAL T3 report writers.

No new report layout: calls write_bsm_uv_field_table, write_rge_comparison,
write_eft1_rge_comparison and the original saved-RGE group-factor tables.
Only renormalisable EFT1 RGEs are available; no general matching is claimed.
"""
from __future__ import annotations

import csv
import json
from pathlib import Path
from types import SimpleNamespace

from Reports.ReportGeneration import (
    PAPER_SYMBOL_LATEX, compile_latex_document, latex_escape_text,
    paper_notation_key_lines as original_notation,
    write_bsm_uv_field_table,
)
from Reports.RGEComparison import write_rge_comparison, write_eft1_rge_comparison
from Reports.GroupFactorReports import (
    _group_factor_document_header, _representation_table,
    _saved_rge_beta_sections,
)
from Reports import ReportGeneration as RG
from Reports import RGEComparison as RC
from Reports import GroupFactorReports as GF
from Numerical.experimental.GeneralT3ScalarQuadraticMasses import REP
from scripts.RunT3GeneralSymbolicStudy import DEFAULT_DATA_ROOT, DEFAULT_REPORT_ROOT

# Reuse the existing conversion/notation functions, but not its document writers.
from Reports.GeneralSymbolicPresentation import (
    generated_symbol, format_general_latex, field_content, notation_rows,
    matchete_terms, obtain_beta_latex,
)


def _json(path: Path):
    if not path.is_file():
        raise FileNotFoundError(path)
    return json.loads(path.read_text(encoding='utf-8-sig'))


def _model_tag(model):
    return model.replace('+', 'p').replace('-', 'm')


def _shared_notation(models):
    result = original_notation(heading='Coupling and interaction notation')
    rows = notation_rows(models)
    result += [r'\section*{Additional general-interaction couplings}',
        r'\begin{longtable}{@{}p{0.20\linewidth}p{0.34\linewidth}p{0.38\linewidth}@{}}',
        r'\toprule', r'Symbol & Field content & Stored coupling / contraction \\', r'\midrule']
    for symbol, content, detail in rows:
        result.append('$' + symbol + '$ & $' + content + '$ & ' + detail + r' \\')
    result += [r'\bottomrule', r'\end{longtable}']
    return result


def _compile(path):
    compile_latex_document(path)
    pdf = path.with_suffix('.pdf')
    if not pdf.is_file() or pdf.stat().st_size < 1000:
        raise RuntimeError(f'LaTeX compilation did not produce {pdf}; see compile log')
    return pdf


def render(*, study='full16_general_symbolic', run_id='full16_general_v1', wolframscript='wolframscript'):
    import re
    if not all(re.fullmatch(r'[A-Za-z0-9][A-Za-z0-9_-]{0,63}', x) for x in (study, run_id)):
        raise ValueError('Invalid study or run ID')
    base = DEFAULT_REPORT_ROOT / study / 'GeneralSymbolic'
    root = DEFAULT_DATA_ROOT / run_id
    summary = _json(base / 'summary.json')
    if not summary.get('all_requested_models_validated'):
        raise RuntimeError('General-symbolic study is incomplete')
    with (base / 'summary.csv').open(encoding='utf-8-sig', newline='') as f:
        models = [row['model'] for row in csv.DictReader(f)]
    if not models or len(models) != int(summary['model_count_validated']):
        raise RuntimeError('General-symbolic model inventory is incomplete')

    # Presentation-only staged JSON files allow the original loaders to run
    # unchanged, using their normal RunRecord.summary keys.
    staging = base / '_original_report_inputs'
    staging.mkdir(parents=True, exist_ok=True)
    records = []
    all_names = {'UV': set(), 'EFT1': set()}
    for model in models:
        kind, alpha_string = model.split('_')
        d1, d2, df = REP[kind]
        folder = root / _model_tag(model)
        lag = _json(folder / 'lagrangian.json')
        if not lag.get('general_scalar_quadratics_requested') or not lag.get('extra_yukawas_requested'):
            raise RuntimeError(f'{model}: restricted Lagrangian encountered')
        lag_path = Path(lag['bsm_lagrangian_inputform_file'])
        terms = matchete_terms(lag_path, staging / f'{_model_tag(model)}_bsm_tex.json', wolframscript)
        # The original field-content sorter handles these terms, instead of a
        # model-specific independent Lagrangian layout.
        lag_tex = ' + '.join(format_general_latex(term) for term in terms)
        saved_dir = staging / _model_tag(model)
        saved_dir.mkdir(exist_ok=True)
        values = {}
        for stage in ('UV', 'EFT1'):
            source = folder / (stage.lower() + '_rgbeta.json')
            payload = _json(source)
            betas = payload.get('report_betas')
            if not isinstance(betas, dict) or not betas:
                raise RuntimeError(f'{model}: {stage} beta functions unavailable')
            latex = obtain_beta_latex(payload, source, saved_dir / stage, wolframscript)
            if set(latex) < set(betas):
                raise RuntimeError(f'{model}: missing {stage} beta LaTeX')
            all_names[stage].update(betas)
            # Data compatibility only. Original RGE writers read the file via
            # UVRGEFile/EFT1RenormalisableRGEFile and status=Success.
            prepared = dict(payload)
            prepared['status'] = 'Success'
            prepared['report_beta_latex'] = {n: format_general_latex(tex) for n, tex in latex.items()}
            name = stage.lower() + '.json'
            (saved_dir / name).write_text(json.dumps(prepared, indent=2), encoding='utf-8')
            values[stage] = name
        records.append(SimpleNamespace(
            name=kind, alpha=int(alpha_string), d_s1=d1, d_s2=d2, d_f=df,
            shared_scalar=False, output_dir=saved_dir, return_code=0, eft_stages=[],
            summary={
                'BuildStatus': 'Success',
                'BSMUVLagrangianLaTeX': lag_tex,
                'UVRGEStatus': 'Success', 'UVRGEFile': values['UV'],
                'EFT1RenormalisableRGEStatus': 'Success',
                'EFT1RenormalisableRGEFile': values['EFT1'],
                'EFT1WilsonRGEStatus': 'NotComputed',
                'MatchingStatus': 'NotPerformed',
            },
        ))

    # Make the existing shared notation key include the model general couplings.
    # No separate original-writer layout is created.
    original_map = dict(PAPER_SYMBOL_LATEX)
    original_notation_fn = RG.paper_notation_key_lines
    original_rc_notation = RC.paper_notation_key_lines
    original_gf_notation = GF.paper_notation_key_lines
    for name in all_names['UV'] | all_names['EFT1']:
        symbol = generated_symbol(name)
        if symbol:
            PAPER_SYMBOL_LATEX[name] = '{' + symbol + '}'

    def notation(*, heading='Notation key'):
        return _shared_notation(models)

    RG.paper_notation_key_lines = notation
    RC.paper_notation_key_lines = notation
    GF.paper_notation_key_lines = notation
    outputs = []
    try:
        # THESE ARE THE ORIGINAL WRITER FUNCTIONS, UNALTERED.
        uv = write_bsm_uv_field_table(records, report_root=base / 'Lagrangian')
        outputs.append(_compile(uv))
        uv_rge = write_rge_comparison(records, report_root=base / 'RGE')
        outputs.append(_compile(uv_rge))
        eft_rge = write_eft1_rge_comparison(records, report_root=base / 'RGE')
        outputs.append(_compile(eft_rge))

        for stage, label, loader in (
            ('UV', 'GF_UV', RC.load_uv_rge_payload),
            ('EFT1', 'GF_EFT1', RC.load_eft1_renormalisable_rge_payload),
        ):
            # The SAME existing GF RGE structure/prefactor table generator.
            lines = _group_factor_document_header(
                f'T3 group factors: general {stage}',
                r'One-loop renormalisable RGBeta structure prefactors across T3 models.',
            )
            lines += notation()
            lines += _representation_table(records, include_f=(stage == 'UV'))
            lines += _saved_rge_beta_sections(
                records, loader=loader, couplings=sorted(all_names[stage]),
                title='Complete general-interaction RGE prefactors',
            )
            lines += [r'\end{document}', '']
            tex = base / 'GroupFactors' / stage / (label + '.tex')
            tex.parent.mkdir(parents=True, exist_ok=True)
            tex.write_text('\n'.join(lines), encoding='utf-8')
            outputs.append(_compile(tex))
    finally:
        PAPER_SYMBOL_LATEX.clear()
        PAPER_SYMBOL_LATEX.update(original_map)
        RG.paper_notation_key_lines = original_notation_fn
        RC.paper_notation_key_lines = original_rc_notation
        GF.paper_notation_key_lines = original_gf_notation

    index = ['# Original T3 reports with the general BSM interaction data', '',
        'These PDFs are generated by the existing field-content, RGE comparison and',
        'group-factor table functions. No general threshold matching was performed.', '']
    for path in outputs:
        index.append(f'- [{path.name}]({path.relative_to(base).as_posix()})')
    (base / 'PDF_REPORTS.md').write_text('\n'.join(index) + '\n', encoding='utf-8')
    return outputs
