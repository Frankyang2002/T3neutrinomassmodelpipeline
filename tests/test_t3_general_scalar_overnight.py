import csv
from pathlib import Path
import pytest
from scripts import RunT3GeneralScalarOvernight as overnight

def test_identifiers():
    assert overnight.check_id('night_1') == 'night_1'
    with pytest.raises(ValueError): overnight.check_id('../oops')

def test_cases():
    assert len(overnight.CASES)==3 and len(overnight.MODELS)==16

def test_report_complete(tmp_path):
    f=tmp_path/'comparison.csv'
    fieldnames=['case','model','yukawa_mode','status','chi2','donor_status',*overnight.FIELDS]
    with f.open('w',newline='',encoding='utf-8') as stream:
        w=csv.DictWriter(stream,fieldnames=fieldnames);w.writeheader()
        for c in overnight.CASES:
            for m in overnight.MODELS:
                for mode in ('optimised','shared_a_m2'):
                    w.writerow(dict(case=c,model=m,yukawa_mode=mode,status='Pending'))
    report=tmp_path/'index.html'
    result=overnight.render_report(f,report)
    assert sum(result.values())==96
    assert 'Full sequential Wilson' not in report.read_text()
    assert 'conditional candidate' in report.read_text().lower()

def test_report_rejects_missing(tmp_path):
    f=tmp_path/'empty.csv';f.write_text('case,model,yukawa_mode\n',encoding='utf-8')
    with pytest.raises(ValueError):overnight.render_report(f,tmp_path/'report.html')

def test_prepare_only_parser_rejects_non_power_of_two():
    with pytest.raises(SystemExit):overnight.main(['--run-id','x','--sobol-points','7'])
