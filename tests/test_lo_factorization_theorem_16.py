"""Batch 7: analytic series, source-audit and 16-model no-overwrite validation."""
from __future__ import annotations

import json
import math
import os
from pathlib import Path
import numpy as np
import pytest

from Numerical.lo_patterns_16 import factorization_theorem as t7
from Numerical.lo_patterns_16 import kernel_factorization as k6
from Numerical.lo_patterns_16 import texture_synthesis as t5
from Numerical.lo_patterns_16.analytic_explanation import (
    b0, b0_derivative, divdiff, scalar_kernels
)
from Numerical.lo_certification.direct_uv.native_amplitude_contract import model_dictionary
from Numerical.lo_certification.direct_uv.certification import charge_bases


@pytest.mark.parametrize('m2,a,b',[(3e4,2e5,3.2e5),(9.8e8,1.15e6,8e5),
                                   (2.1e6,4e6,1.3e6)])
def test_partial_derivatives_against_symmetric_finite_difference(m2,a,b):
    D,da,db=t7.derivative_D(m2,a,b)
    h=5e-5
    fa=(divdiff(m2,a*(1+h),b)-divdiff(m2,a*(1-h),b))/(2*h*a)
    fb=(divdiff(m2,a,b*(1+h))-divdiff(m2,a,b*(1-h)))/(2*h*b)
    assert np.isclose(da,fa,rtol=2e-5,atol=1e-17)
    assert np.isclose(db,fb,rtol=2e-5,atol=1e-17)
    assert D==pytest.approx(divdiff(m2,a,b),rel=1e-13)


@pytest.mark.parametrize('m2,a,b',[(9e4,3e5,8e5),(5e7,8e5,1.4e6)])
def test_exact_scalar_kernel_cubic_perturbation_converges(m2,a,b):
    errors=[]
    for d in [12000.,6000.,3000.]:
        exact=scalar_kernels(m2,a,b,complex(d,-d/3))[0]
        lin,cub=t7.kernel_series(m2,a,b,complex(d,-d/3))
        assert abs(cub-exact)<abs(lin-exact)
        errors.append(abs(cub-exact))
    assert errors[2]<errors[1]/15
    assert errors[1]<errors[0]/15


def test_complex_scalar_phase_carries_through_cubic_kernel():
    m2,a,b=1e5,2e5,4e5
    lin,cub=t7.kernel_series(m2,a,b,350*(1+1j))
    lin2,cub2=t7.kernel_series(m2,a,b,350*(1-1j))
    assert lin2==pytest.approx(lin.conjugate(),rel=1e-14)
    assert cub2==pytest.approx(cub.conjugate(),rel=1e-14)


def test_cubic_only_odd_mixing_and_sign_change():
    lin,cub=t7.kernel_series(1e5,2e5,4e5,250)
    lin2,cub2=t7.kernel_series(1e5,2e5,4e5,-250)
    assert lin2==-lin and cub2==-cub


def test_zero_mixing_gives_zero_mass_cross_kernel():
    lin,cub=t7.kernel_series(1e5,2e5,4e5,0j)
    assert lin==cub==0


def test_non_degenerate_derivative_guard():
    with pytest.raises(t7.TheoremError,match='Nondegenerate'):
        t7.derivative_D(1e5,2e5,2e5)


def test_complex_tensor_factorization_and_bound_for_arbitrary_yukawas():
    rng=np.random.default_rng(33)
    a=rng.normal(size=(3,3,3))+1j*rng.normal(size=(3,3,3))
    e=1e-3*(rng.normal(size=(3,3,3))+1j*rng.normal(size=(3,3,3)))
    k=1.4-0.9j
    b=k*a+e
    fit=t7.compare_coefficient_arrays(a,b)
    assert fit['kappa']['real']==pytest.approx(k.real,abs=.001)
    assert fit['kappa']['imag']==pytest.approx(k.imag,abs=.001)
    y1=rng.normal(size=(3,3))+1j*rng.normal(size=(3,3))
    y2=rng.normal(size=(3,3))+1j*rng.normal(size=(3,3))
    def mat(t):
        z=sum((y1.conj()@g@y2.conj().T for g in t),np.zeros((3,3),complex))
        return k6.LOOP*(z+z.T)
    measured=k6.norm(mat(b)-k*mat(a))
    bound=2*k6.LOOP*k6.norm(y1)*k6.norm(y2)*sum(k6.norm(g) for g in e)
    assert measured<=bound*(1+1e-12)


def _toy_source_audit(models,ks=None):
    if ks is None:ks={m['model']:1. for m in models}
    rows=[]
    for item in models:
        model=item['model'];kind=model[0];alpha=int(model[2:]);dim=item['SU2_dimensions']
        charge=charge_bases(kind,alpha)
        set1=charge['scalar'][:dim['S1']]
        set2=charge['scalar'][dim['S1']:]
        ferm=[(dim['F']-1)/2-f+(alpha+1)/2 for f in range(dim['F'])]
        pairs=[(i,j) for i,x in enumerate(set1) for j,y in enumerate(set2) if abs(x-y)<1e-10]
        n=2 if kind=='E' else 1
        pairs=pairs[:n]
        assert len(pairs)==n
        channels=[]
        for i,j in pairs:
            q=set1[i];f=next(f for f,v in enumerate(ferm) if abs(q-v)<1e-10)
            w=-2*math.sqrt(3)/3*ks[model]/n
            channels.append({'charge':str(q),'indices_weight_basis_one_based':
                             {'F':f+1,'S1':i+1,'S2':j+1},
                             'y1':str(w),'y2':'1','mixing':'1',
                             'complex_weight':str(w),'nonzero':True})
        rows.append({'model':{'dS1':dim['S1'],'dS2':dim['S2'],
                              'dF':dim['F'],'alpha':alpha},'class':kind,
                     'channels':channels,
                     'total_group_coefficient':str(-2*math.sqrt(3)/3*ks[model])})
    return {'status':'charge_resolved_group_coefficients','count':16,
            'contraction':'Y1 * Y2 * conjugate(T3MixCG) using same physical internal weight indices',
            'models':rows}


def test_all16_source_CG_audit_validator_on_synthetic_consistent_table():
    models=model_dictionary();w=t7.expected_model_weights(_toy_source_audit(models),[x['model'] for x in models])
    assert len(w)==16
    assert all(x['active_charge_count']==(2 if x['family']=='E' else 1) for x in w.values())
    assert all(abs(x['magnitude_relative_to_ref']-1)<1e-12 for x in w.values())


def test_corrupt_CG_weight_and_contraction_rejected():
    models=model_dictionary();data=_toy_source_audit(models)
    data['models'][0]['channels'][0]['complex_weight']='42'
    with pytest.raises(Exception,match='CG product'):
        t7.expected_model_weights(data,[x['model'] for x in models])
    data=_toy_source_audit(models);data['contraction']='Y1*Y2*mix'
    with pytest.raises(Exception,match='contraction'):
        t7.expected_model_weights(data,[x['model'] for x in models])


def test_invalid_charge_index_rejected():
    models=model_dictionary();data=_toy_source_audit(models)
    data['models'][0]['channels'][0]['indices_weight_basis_one_based']['F']=999
    with pytest.raises(Exception,match='CG component index'):
        t7.expected_model_weights(data,[x['model'] for x in models])


def test_output_refuses_overwrite_and_outside_dir(tmp_path):
    target=tmp_path/t7.TARGET
    target.mkdir(parents=True);(target/'keep.txt').write_text('preserve')
    with pytest.raises(t7.TheoremError,match='already populated'):
        t7.analyse(tmp_path)
    assert (target/'keep.txt').read_text()=='preserve'
    with pytest.raises(t7.TheoremError,match='inside Chat A'):
        t7.analyse(tmp_path,output_dir=tmp_path/'different')


def _synthetic_result(tmp_path):
    from test_lo_texture_synthesis_16 import fixture_16
    root=fixture_16(tmp_path)
    t5.analyse(root,make_plots=False)
    data=k6.analyse(root,make_plots=False)
    table={r['model']:r['absolute_scale'] for r in data['coefficient_comparisons']
           if r['scenario']=='baseline'}
    audit=_toy_source_audit(model_dictionary(),table)
    path=root/t7.AUDIT;path.parent.mkdir(parents=True,exist_ok=True)
    path.write_text(json.dumps(audit))
    return root


def test_all16_synthetic_source_to_theorem_and_no_original_mutation(tmp_path):
    root=_synthetic_result(tmp_path)
    before={str(p):t7.digest(p) for p in root.rglob('*') if p.is_file()}
    p=t7.analyse(root)
    assert p['model_count']==16
    assert len(p['cg_and_factorization'])==16
    assert len(p['correction_tests'])==16
    assert len(p['scenario_summary'])==80
    assert all(row['cubic_improves_over_linear'] for row in p['correction_tests'])
    for file,sha in before.items():assert t7.digest(Path(file))==sha
    out=root/t7.TARGET
    for name in ('BATCH07_FACTORISATION_THEOREM.md','batch07_theorem.json',
                 'batch07_cg_weight_table.csv','batch07_first_corrections.csv',
                 'batch07_scenario_residuals.csv','SOURCE_MANIFEST.json','HANDOFF.json'):
        assert (out/name).exists(),name
    report=(out/'BATCH07_FACTORISATION_THEOREM.md').read_text()
    assert '0/16' in report and 'C_-1' in report and 'source' in report.lower()


def test_source_input_tampering_prevents_report(tmp_path):
    root=_synthetic_result(tmp_path)
    path=root/'output/general_observable/lo_direct_one_loop/provisional_inputs/C_+1.json'
    path.write_text(path.read_text()+'\n')
    with pytest.raises(t7.TheoremError,match='Source SHA drift'):
        t7.analyse(root)
    assert not (root/t7.TARGET).exists()


def test_false_physical_certification_prevents_report(tmp_path):
    root=_synthetic_result(tmp_path)
    path=root/t7.BATCH06/'batch06_factorization.json';obj=t7.load(path)
    obj['physical_normalisation_certified']=True
    path.write_text(json.dumps(obj))
    with pytest.raises(t7.TheoremError,match='physical certification'):
        t7.analyse(root)
    assert not (root/t7.TARGET).exists()


def test_modified_batch6_manifest_prevents_report(tmp_path):
    root=_synthetic_result(tmp_path)
    path=root/t7.BATCH06/'SOURCE_MANIFEST.json';obj=t7.load(path)
    obj['models']['A_-2']['input_sha256']='0'*64
    path.write_text(json.dumps(obj))
    with pytest.raises(t7.TheoremError,match='manifest hash chain'):
        t7.analyse(root)
    assert not (root/t7.TARGET).exists()


@pytest.mark.skipif(os.environ.get('RUN_LO_THEOREM_INTEGRATION')!='1',
                    reason='Requires actual local frozen 16-model source and a new analysis output directory')
def test_opt_in_real_checkout(tmp_path):
    root=Path(__file__).resolve().parents[1]
    p=t7.analyse(root,output_dir=root/t7.BASE/'analysis_batch07_test_integration_01')
    assert p['model_count']==16 and p['physical_source_certified_count']==0
