"""Chat A Batch 6: coefficient proof, stress counterexamples and source-integrity tests.

Pure algebra unit tests run in every pytest invocation. Optional genuine
integration requires RUN_LO_KERNEL_INTEGRATION=1 and an unoccupied output dir.
"""
from __future__ import annotations
import json
import os
from hashlib import sha256
from pathlib import Path

import numpy as np
import pytest

from Numerical.lo_patterns_16 import kernel_factorization as k
from Numerical.lo_patterns_16 import texture_synthesis as t5
from Numerical.lo_patterns_16.analytic_explanation import reconstruct
from Numerical.lo_certification.direct_uv.native_amplitude_contract import model_dictionary


def rand_y(seed=6):
    rng=np.random.default_rng(seed)
    return (np.eye(3)+rng.normal(size=(3,3))*.1+1j*rng.normal(size=(3,3))*.2)*.001,\
           (np.eye(3)+rng.normal(size=(3,3))*.1+1j*rng.normal(size=(3,3))*.2)*.001


def base_toy(model='C_+1',seed=8):
    from test_lo_texture_synthesis_16 import _synthetic_model_data
    item=next(r for r in model_dictionary() if r['model']==model)
    inp,frozen=_synthetic_model_data(item,seed)
    return inp,t5.extract_source_tensors(inp,frozen,model),item['SU2_dimensions']['S1']


def sha(p): return sha256(p.read_bytes()).hexdigest()


def test_finite_dimensional_factorization_exact_for_arbitrary_complex_yukawas():
    rng=np.random.default_rng(12)
    a=[rng.normal(size=(3,3))+1j*rng.normal(size=(3,3)) for _ in range(3)]
    b=[{ 'm2': float(10**(i+2)), 'B':x,'B_leading':x } for i,x in enumerate(a)]
    one={'groups':b}
    two={'groups':[{'m2':g['m2'],'B':(2-1.3j)*g['B'],
                    'B_leading':(2-1.3j)*g['B_leading']} for g in b]}
    cmp=k.compare_coeff(one,two)
    assert cmp['matrix_residual']<1e-14
    assert cmp['max_group_scale_fractional_spread']<1e-14
    for seed in [1,2,3,18,33]:
        y1,y2=rand_y(seed)
        assert k.rel(k.assemble(two,y1,y2)['matrix'],(2-1.3j)*k.assemble(one,y1,y2)['matrix'])<1e-14


def test_group_by_group_is_not_sufficient_if_coefficients_differ():
    a=[np.diag([1.+i,2.+i,3.+i]).astype(complex) for i in range(3)]
    x={'groups':[{'m2':float(i+1),'B':a[i]} for i in range(3)]}
    y={'groups':[{'m2':float(i+1),'B':a[i]*c} for i,c in enumerate([1,2,3])]}
    result=k.compare_coeff(x,y)
    assert result['max_group_shape_residual']<1e-14
    assert result['matrix_residual']>.1
    assert result['max_group_scale_fractional_spread']>.1


def test_spectral_weyl_mass_tilt_preserves_symmetry_and_spectrum():
    inp,tensor,d1=base_toy('B_-1')
    M=t5.decode(inp['fermion_mass_gev'],'M')
    out=k.deform_heavy_generations(inp,(.71,1.19,1.64))
    N=t5.decode(out['fermion_mass_gev'],'N')
    assert k.rel(N,N.T)<1e-12
    assert k.rel(M,t5.decode(inp['fermion_mass_gev'],'M'))==0
    got=np.unique(np.round(np.linalg.svd(N,compute_uv=False),8))
    assert np.allclose(np.sort(got),np.sort(np.array([120.,175.,260.])*[.71,1.19,1.64]))


def test_spectral_tilt_rejects_invalid_factors():
    inp,_,_=base_toy()
    for x in [(1.,2.), (1.,0.,2.), (-1.,1.,1.), (1.,float('inf'),2.)]:
        with pytest.raises(k.FactorizationError,match='positive factor'):
            k.deform_heavy_generations(inp,x)


@pytest.mark.parametrize('model',['A_-2','B_-1','C_+1','C_-1','D_+0','E_-4'])
def test_extracted_coefficient_tensor_matches_independent_b3_kernel(model):
    inp,ts,d1=base_toy(model)
    coeff=k.coefficient_groups(inp,ts,model,d1)
    yy1,yy2=rand_y(123)
    new=t5.substitute_textures(inp,ts,yy1,yy2)
    lhs=reconstruct(new,model[0],int(model[2:]),d1)
    rhs=k.assemble(coeff,yy1,yy2)
    assert k.rel(lhs['exact'],rhs['matrix'])<1e-9
    assert all(k.rel(g['matrix'],h)<1e-9 for g,h in zip(lhs['heavy_groups'],rhs['groups']))


@pytest.mark.parametrize('scenario',k.SCENARIOS)
def test_all_scenarios_preserve_kernel_formula(scenario):
    inp,ts,d1=base_toy('E_-2')
    changed=k.deform(inp,d1,scenario)
    coeff=k.coefficient_groups(changed,ts,'E_-2',d1)
    y1,y2=rand_y(19)
    expected=reconstruct(t5.substitute_textures(changed,ts,y1,y2),'E',-2,d1)
    assert k.rel(expected['exact'],k.assemble(coeff,y1,y2)['matrix'])<1e-9


def test_charge_splitting_breaks_equal_E_channel_kernel_norms():
    inp,ts,d1=base_toy('E_-2')
    original=k.coefficient_groups(inp,ts,'E_-2',d1)
    broken=k.coefficient_groups(k.deform(inp,d1,'charge_split_diagnostic'),ts,'E_-2',d1)
    def two(x):
        return [np.sqrt(sum(g['coefficient_norm']**2 for g in cs)) for cs in x['charges'].values()]
    orig=two(original);new=two(broken)
    assert len(orig)==len(new)==2
    # The synthetic two sectors may differ already (different source CG weights),
    # so assert relative *change* instead of a physical equality claim.
    assert abs(new[0]/new[1]-orig[0]/orig[1])>1e-3


def test_scalar_mass_negative_diagonal_rejected():
    inp,ts,d1=base_toy()
    A=t5.decode(inp['scalar_A_gev2'],'A')
    A[0,0]=-10
    inp['scalar_A_gev2']=t5.encode(A)
    with pytest.raises(k.FactorizationError,match='diagonal tachyon'):
        k.deform(inp,d1,'baseline')


def test_fermion_nonsymmetric_mass_rejected():
    inp,ts,d1=base_toy()
    M=t5.decode(inp['fermion_mass_gev'],'M');M[0,1]=5j
    inp['fermion_mass_gev']=t5.encode(M)
    with pytest.raises(k.FactorizationError,match='symmetric'):
        k.coefficient_groups(inp,ts,'C_+1',d1)


def test_refuse_existing_output_and_no_mutation(tmp_path):
    root=tmp_path.resolve();target=root/k.TARGET;target.mkdir(parents=True)
    p=target/'existing.json';p.write_text('{"keep":1}')
    with pytest.raises(k.FactorizationError,match='already populated'):
        k.analyse(root,make_plots=False)
    assert p.read_text()=='{"keep":1}'


def test_refuse_output_outside_project(tmp_path):
    with pytest.raises(k.FactorizationError,match='distinct and contained'):
        k.analyse(tmp_path.resolve(),output_dir=tmp_path/'another',make_plots=False)


def _fixture(tmp_path):
    from test_lo_texture_synthesis_16 import fixture_16
    root=fixture_16(tmp_path.resolve())
    t5.analyse(root,make_plots=False)
    return root


def test_synthetic_complete_all16_and_source_preservation(tmp_path):
    root=_fixture(tmp_path)
    before={str(p):sha(p) for p in root.rglob('*') if p.is_file()}
    results=k.analyse(root,make_plots=False)
    assert results['model_count']==16
    assert results['physical_normalisation_certified'] is False
    assert len(results['coefficient_comparisons'])==80
    assert len(results['flavour_matrix_controls'])==160
    assert len(results['tensor_reconstruction_checks'])==80
    assert len(results['E_channel_diagnostics'])==12
    assert all(sha(Path(p))==h for p,h in before.items())
    out=root/k.TARGET
    assert (out/'BATCH06_ANALYTICAL_FACTORISATION.md').exists()
    assert (out/'batch06_factorization.json').exists()
    assert (out/'batch06_coefficient_tests.csv').exists()


def test_synthetic_tampered_frozen_input_rejected(tmp_path):
    root=_fixture(tmp_path)
    p=root/'output/general_observable/lo_direct_one_loop/provisional_inputs/C_+1.json'
    z=json.loads(p.read_text());z['model']='C_-1';p.write_text(json.dumps(z))
    with pytest.raises(k.FactorizationError,match='SHA'):
        k.analyse(root,make_plots=False)
    assert not (root/k.TARGET).exists()


def test_synthetic_false_certificate_rejected(tmp_path):
    root=_fixture(tmp_path)
    p=root/k.BATCH05/'batch05_synthesis.json'
    z=json.loads(p.read_text());z['physical_normalisation_certified']=True
    p.write_text(json.dumps(z))
    with pytest.raises(k.FactorizationError,match='physical certification'):
        k.analyse(root,make_plots=False)
    assert not (root/k.TARGET).exists()


def test_synthetic_batch05_manifest_chain_refused(tmp_path):
    root=_fixture(tmp_path)
    p=root/k.BATCH05/'SOURCE_MANIFEST.json'
    z=json.loads(p.read_text());z['study_sha256']='0'*64
    p.write_text(json.dumps(z))
    with pytest.raises(k.FactorizationError,match='hash chain'):
        k.analyse(root,make_plots=False)
    assert not (root/k.TARGET).exists()


@pytest.mark.skipif(os.environ.get('RUN_LO_KERNEL_INTEGRATION')!='1',
                    reason='Requires real frozen Windows checkout and new test output directory')
def test_real_frozen_all16_integration():
    root=Path(__file__).resolve().parents[1]
    result=k.analyse(root,output_dir=root/k.BASE/'analysis_batch06_pytest_integration_01',make_plots=False)
    assert result['model_count']==16
    assert all(x['exact_tensor_reconstruction_error']<5e-8 for x in result['tensor_reconstruction_checks'])
