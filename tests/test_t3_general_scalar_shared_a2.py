"""Offline schema and bookkeeping tests; no physics integration is claimed."""
import importlib.util
import json
from pathlib import Path

import pytest

SCRIPT=Path(__file__).resolve().parents[1]/'scripts'/'CompareT3GeneralScalarSharedA2.py'
spec=importlib.util.spec_from_file_location('t3_shared_a2',SCRIPT)
assert spec and spec.loader
code=importlib.util.module_from_spec(spec)
spec.loader.exec_module(code)


def fake_rows():
    return {(c,m):{'study':c,'model':m,'status':'NoMassSafeNeutrinoFit',
                   'chi2':None} for c in code.CASES for m in code.MODELS}


def test_exact_six_cases():
    assert len(code.MODELS)==16
    assert len(code.CASES)==3
    assert len({(x,y) for x in code.CASES for y in ('optimised','shared_a_m2')})==6


def test_reads_exact_48_and_fails_closed(tmp_path):
    doc={'format':'t3-48-model-general-quartic-conditional-campaign-v1',
         'conditional_only':True,'no_sequential_wilson_matching':True,
         'rows':list(fake_rows().values())}
    (tmp_path/'campaign.json').write_text(json.dumps(doc),encoding='utf-8')
    found=code.campaign_rows(tmp_path)
    assert len(found)==48
    doc['rows']=doc['rows'][:-1]
    (tmp_path/'campaign.json').write_text(json.dumps(doc),encoding='utf-8')
    with pytest.raises(ValueError,match='48'):
        code.campaign_rows(tmp_path)


def test_all_96_rows_and_missing_not_filled():
    original=fake_rows()
    shared=[{'case':c,'model':m,'status':'UnavailableDonor','donor_status':'NoMassSafeA2Donor'}
            for c in code.CASES for m in code.MODELS]
    rows=code.summarize(original,shared)
    assert len(rows)==96
    assert all(r['conditional_matching_only'] for r in rows)
    assert all(r['chi2'] is None for r in rows)
    assert sum(r['yukawa_mode']=='optimised' for r in rows)==48
    assert sum(r['yukawa_mode']=='shared_a_m2' for r in rows)==48


def test_missing_donor_is_never_replaced_with_historical(tmp_path):
    for case in code.CASES:
        d=code.donor(tmp_path,case,allow_poor=False)
        assert d['status']=='MissingDonorFit'


def test_atomic_output(tmp_path):
    p=tmp_path/'data'/'x.json'
    code.atomic(p,{'x':5})
    assert code.read(p)=={'x':5}
    assert not p.with_name(p.name+'.tmp').exists()
