"""No physics checks: file routing and PDF emission only."""
from pathlib import Path
import csv
import json
import importlib.util
import shutil

SCRIPT = Path(__file__).resolve().parents[2] / 'scripts' / 'WriteT3GeneralSymbolicPDFs.py'
spec = importlib.util.spec_from_file_location('general_pdf_writer', SCRIPT)
module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)


def test_pdf_routing(tmp_path, monkeypatch):
    monkeypatch.setattr(module, 'DATA_ROOT', tmp_path / 'data')
    monkeypatch.setattr(module, 'REPORT_ROOT', tmp_path / 'reports')
    study = 'trial'
    report = module.REPORT_ROOT / study / 'GeneralSymbolic'
    folder = module.DATA_ROOT / 'full16_general_v1' / 'A_m2'
    report.mkdir(parents=True)
    folder.mkdir(parents=True)
    (report / 'summary.json').write_text(json.dumps({'all_requested_models_validated': True, 'model_count_validated': 1}))
    with (report / 'summary.csv').open('w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=['model','status','extra_quartics','extra_bilinears','extra_yukawa','mass_basis'])
        w.writeheader()
        w.writerow(dict(model='A_-2',status='SymbolicRGExportedMatchingUnverified', extra_quartics=5,extra_bilinears=1,extra_yukawa='y2Conj',mass_basis='CheckedNotMatched'))
    with (report / 'group_invariants.csv').open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=['model','field','d_SU2','Y','C2_SU2','T_SU2'])
        w.writeheader()
        for field in ['S1','S2','F']:
            w.writerow(dict(model='A_-2',field=field,d_SU2=2,Y=0,C2_SU2='3/4',T_SU2='1/2'))
    bsm = folder / 'bsm.txt'
    bsm.write_text('Plus[NewScalar1[],y2Conj[],lambdaT3GeneralQ000004I1[]]\n')
    (folder / 'lagrangian.json').write_text(json.dumps({'general_scalar_quadratics_requested': True, 'extra_yukawas_requested': True, 'bsm_lagrangian_inputform_file': str(bsm), 'interactions': ['Yukawa2ConjugateScalar_GenericCG','T3ExtraQ_X']}))
    for st in ['uv','eft1']:
        (folder / f'{st}_rgbeta.json').write_text(json.dumps({'report_betas': {'g2':'-19/6*g2^3','lambdaT3GeneralQ000004I1':'lambdaT3GeneralQ000004I1^2', **({'y2Conj':'y2Conj^3'} if st=='uv' else {})}}))
    if shutil.which('pdflatex'):
        files = module.render(study=study)
        assert {p.name for p in files} == {'UV.pdf','EFT_1_after_F.pdf','GroupFactors.pdf'}
        assert all(p.stat().st_size>1000 for p in files)
    else:
        assert module.tag('A_-2')=='A_m2'
