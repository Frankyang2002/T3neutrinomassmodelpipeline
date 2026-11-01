"""Functional test: PDF export works when a TeX engine is absent or fails."""
from __future__ import annotations
import csv
import importlib.util
import json
from pathlib import Path


def test_native_fallback(tmp_path):
    source = Path(__file__).resolve().parents[2] / 'scripts' / 'WriteT3GeneralSymbolicPDFs.py'
    spec = importlib.util.spec_from_file_location('t3_general_pdf_fallback_under_test', source)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    out = tmp_path / 'Reports' / 'output' / 'pilot' / 'GeneralSymbolic'
    native = out / 'PDFNative'
    native.mkdir(parents=True)
    (out / 'summary.json').write_text(json.dumps({'all_requested_models_validated': True,
                                                    'model_count_validated': 1}), encoding='utf-8')
    with (out / 'summary.csv').open('w', newline='', encoding='utf-8') as f:
        w = csv.DictWriter(f, fieldnames=['model','status','extra_quartics','extra_bilinears','extra_yukawa'])
        w.writeheader()
        w.writerow({'model':'A_-2','status':'SymbolicRGExportedMatchingUnverified',
                    'extra_quartics':5,'extra_bilinears':1,'extra_yukawa':'y2Conj'})
    with (out / 'group_invariants.csv').open('w', newline='', encoding='utf-8') as f:
        w = csv.DictWriter(f, fieldnames=['model','field','d_SU2','Y','C2_SU2','T_SU2'])
        w.writeheader()
        for field, dim, y, c2, t in [('S1',1,'-1','0','0'),('S2',3,'0','2','2'),('F',2,'-1/2','3/4','1/2')]:
            w.writerow({'model':'A_-2','field':field,'d_SU2':dim,'Y':y,'C2_SU2':c2,'T_SU2':t})
    (native/'A_m2_bsm.txt').write_text('MassTerm[muT3GeneralS2Sq] + ExtraQuartic[T3A] + y2Conj[]\n',encoding='utf-8')
    for st in ('uv','eft1'):
        (native/f'A_m2_{st}.txt').write_text('beta_lambdaT3GeneralQ000004I1 =\ng2^2 * lambdaT3GeneralQ000004I1\n',encoding='utf-8')
    for stem in ('UV', 'EFT_1_after_F','GroupFactors'):
        p = out / (stem+'.tex')
        p.write_text('fake tex', encoding='utf-8')
        pdf = module.compile_tex(p, pdflatex='compiler_that_does_not_exist_t3')
        assert pdf.is_file() and pdf.stat().st_size > 1000
    import fitz
    uv = '\n'.join(p.get_text() for p in fitz.open(out/'UV.pdf'))
    assert 'y2Conj' in uv and 'MassTerm' in uv and 'lambdaT3General' in uv
    eft = '\n'.join(p.get_text() for p in fitz.open(out/'EFT_1_after_F.pdf'))
    assert 'lambdaT3General' in eft
    group = '\n'.join(p.get_text() for p in fitz.open(out/'GroupFactors.pdf'))
    assert 'A_-2' in group and 'S2' in group
