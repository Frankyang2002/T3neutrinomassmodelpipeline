"""Offline tests; no Wolfram or repository source copies required."""
from __future__ import annotations

import json
from pathlib import Path
import zipfile

import numpy as np
import pytest
import sympy as sp

from Numerical.running.T3OrdinaryHardFlavor import HardFlavorError
from Numerical.running.T3OrdinarySixFixedOrder import (
    MODELS, SYM, Y1,Y2,LAM,MF,MS1,MS2,MU2,EPS,HBAR,
    parse_one_generation_coeff, inherited_prefactor, audit_inherited,
    verify_native_hard, evaluate, matrices, compact,
)
from scripts.RunT3GeneralOrdinaryScalarHardBatch5 import (
    SUPPORTED, alpha_argument, select_models, preflight, EXPECTED_WILSON_TERMS,
)
from scripts.RunT3OrdinarySixFixedOrder import select_models as select_fixed_models

TEST_ROOT=Path(__file__).resolve().parents[1]
SOURCE_FIXTURE=TEST_ROOT/'fixtures/t3_ordinary_hard_c5_candidate.txt'


def pole_text(family: str) -> str:
    v={'B':'(2*hbar*','C':'(-2*hbar*','E':'(-2*Sqrt[2/3]*hbar*'}[family]
    numerator=('Bar[Coupling[y1, {Index[i, Flavor], Index[r, NFlavor]}, 0]]*'
               'Bar[Coupling[y2, {Index[j, Flavor], Index[r, NFlavor]}, 0]]*'
               'Coupling[lambdaT3, {}, 0]')
    denom={'B':')/(Sqrt[3]*Coupling[MF, {}, 0])',
           'C':')/(3*Coupling[MF, {}, 0])',
           'E':')/Coupling[MF, {}, 0]'}[family]
    return v+numerator+denom


def pole_report(family: str) -> dict:
    return {'Status':'Success','ConsistencyValidated':True,
      'DirectLogEqualsTwicePoleResidueOneGeneration':True,
      'DirectLogToPoleRatioOneGenerationInputForm':'2',
      'HardPoleResidueOneGenerationInputForm':pole_text(family)}


def hard_expression(family: str) -> str:
    base=SOURCE_FIXTURE.read_text(encoding='utf-8-sig')
    multiplier={'B':'1','C':'-1/Sqrt[3]','E':'-Sqrt[2]'}[family]
    return f'({multiplier}*({base}))' if multiplier!='1' else base


def create_inherited(tmp:Path,model:str)->Path:
    """Six exact synthetic source-format triplets, not real Matchete output."""
    fam=model[0]
    raw=(
      'Bar[Coupling[y1, {Index[i, Flavor], Index[r, NFlavor]}, 0]]*'
      'Bar[Coupling[y2, {Index[j, Flavor], Index[r, NFlavor]}, 0]]*'
      'Coupling[lambdaT3, {}, 0]*Field[H, Scalar, {}, {}]*'
      'Field[H, Scalar, {}, {}]*NCM[Field[l, Fermion, {}, {}],'
      ' DiracProduct[GammaCC], Field[l, Fermion, {}, {}]]')
    hc=(
      'Coupling[y1, {Index[i, Flavor], Index[r, NFlavor]}, 0]*'
      'Coupling[y2, {Index[j, Flavor], Index[r, NFlavor]}, 0]*'
      'Bar[Coupling[lambdaT3, {}, 0]]*Field[H, Scalar, {}, {}]*'
      'Field[H, Scalar, {}, {}]*NCM[Field[l, Fermion, {}, {}],'
      ' DiracProduct[GammaCC], Field[l, Fermion, {}, {}]]')
    parts=[]
    for term in (raw,hc):
        for suffix in ('','ϵ*','Log[μbar2/Coupling[MF, {}, 0]^2]*'):
            sig={'B':'+ (2*hbar*','C':'- (2*hbar*','E':'- (2*Sqrt[2/3]*hbar*'}[fam]
            factor=('Sqrt[3]*' if fam=='B' else '3*' if fam=='C' else '')
            # The log multiplies the numerator; epsilon divides the denominator.
            if suffix.startswith('Log['):
                left=sig+term.replace('*Field[H,', '*'+suffix+'Field[H,', 1)
                right=')/('+factor+'Coupling[MF, {}, 0])' if factor else ')/Coupling[MF, {}, 0]'
            else:
                left=sig+term
                right=')/('+factor+suffix+'Coupling[MF, {}, 0])' if factor or suffix else ')/Coupling[MF, {}, 0]'
            parts.append(left+right)
    path=tmp/'inherited.txt';path.write_text('\n'.join(parts),encoding='utf-8')
    return path


def create_hard(tmp:Path,model:str)->Path:
    f=tmp/'hard';f.mkdir()
    fam=model[0]
    (f/'hard_c5_candidate.txt').write_text(hard_expression(fam),encoding='utf-8')
    (f/'hard_scalar_report.json').write_text(json.dumps({
      'Status':'HardScalarCandidateExported','Model':model,
      'HardOneLoopMatched':True,'PerHeavyMassEigenvalueTemplate':True,
      'FermionMatchingRecomputed':False}),encoding='utf-8')
    n={'B':2,'C':4,'E':4}[fam]
    (f/'hard_scalar_input_audit.json').write_text(json.dumps({
      'Status':'Passed','Model':model,'CheckLagrangianPassed':True,
      'TreeWilsonTermCount':n,'IndexedFermionMassFactorsInOriginalTree':n,
      'AdditionalQuarticsRetained':{'B':1,'C':1,'E':0}[fam],
      'DiagonalNFlavorMassIndexEliminatedForTemplate':True}),encoding='utf-8')
    (f/'provenance.json').write_text(json.dumps({
      'model':model,'status':'Ready','tree_wilson_terms':n,
      'additional_quartics_retained':{'B':1,'C':1,'E':0}[fam]}),encoding='utf-8')
    return f


def running_doc(model:str)->dict:
    # Illustrative three-fermion flavour input; no native fit claimed.
    m=lambda a:[[{'real':float(a[i] if i==j else 0),'imag':0} for j in range(3)]
                 for i in range(3)]
    return {'model':model,'status':'NumericalRenormalisableRunningOnly',
      'scales_gev':{'UV':1e7,'F':1e5,'S':1000},
      'uv_at_fermion_threshold':{
        'y1':m([0.0004,0.0003,0.0005]),'y2':m([0.0003,0.0006,0.0004]),
        'MF':m([110000,130000,140000])},
      'eft1_at_fermion_threshold':{'lambdaT3':0.01,
        'mS1Sq':1300**2,'mS2Sq':1500**2},
      'eft1_at_scalar_threshold':{'lambdaT3':0.012,
        'mS1Sq':1400**2,'mS2Sq':1600**2}}


def test_source_parser_accepts_original_D_template():
    expr=parse_one_generation_coeff(SOURCE_FIXTURE.read_text())
    kernel=1/EPS+1+(MS1**2*sp.log(MU2/MS1**2)-MS2**2*sp.log(MU2/MS2**2))/(MS1**2-MS2**2)
    assert sp.cancel(sp.together(expr+2/sp.sqrt(3)*HBAR*Y1*Y2*LAM/MF*kernel))==0


@pytest.mark.parametrize('family,expected',[('B',2/sp.sqrt(3)),('C',-sp.Rational(2,3)),('E',-2*sp.sqrt(sp.Rational(2,3)))])
def test_model_specific_inherited_prefactors(family,expected):
    label=next(x for x in MODELS if x[0]==family)
    assert sp.simplify(inherited_prefactor(pole_report(family),label)-expected)==0


@pytest.mark.parametrize('model',MODELS)
def test_native_hard_symbolic_parity_and_inherited_triplets(tmp_path,model):
    hard=create_hard(tmp_path,model)
    got=verify_native_hard(hard,model,pole_report(model[0]))
    assert got['mass_kernel_proven_by_symbolic_difference']
    archive=create_inherited(tmp_path,model)
    audited=audit_inherited(archive,model,inherited_prefactor(pole_report(model[0]),model))
    assert audited['pureSM_twoL_twoH_terms']==6


@pytest.mark.parametrize('model',MODELS)
def test_numerical_common_coupling_closure(tmp_path,model):
    hard=create_hard(tmp_path,model)
    archive=create_inherited(tmp_path,model)
    out=evaluate(model,hard,running_doc(model),pole_report(model[0]),archive)
    assert out['status']=='InheritedHardAndCommonCouplingPoleCancelled_Conditional'
    assert out['checks']['matching_scale_variation_relative_residual']<3e-12
    assert out['checks']['common_coupling_pole_cancellation_residual']<1e-13
    assert len(out['conditional_masses']['conditional_takagi_masses_ev_ascending'])==3


def test_unknown_extra_quartic_hard_expression_is_rejected(tmp_path):
    folder=create_hard(tmp_path,'B_-3')
    f=folder/'hard_c5_candidate.txt'
    f.write_text(f.read_text().replace('Coupling[lambdaT3, {}, 0]',
      'Coupling[lambdaT3ExtraQ020002I1, {}, 0]'),encoding='utf-8')
    with pytest.raises(HardFlavorError,match='Additional operators|couplings'):
        verify_native_hard(folder,'B_-3',pole_report('B'))


def test_wrong_representation_coefficient_rejected(tmp_path):
    folder=create_hard(tmp_path,'C_-3')
    (folder/'hard_c5_candidate.txt').write_text(hard_expression('B'),encoding='utf-8')
    with pytest.raises(HardFlavorError,match='not the assumed one-loop'):
        verify_native_hard(folder,'C_-3',pole_report('C'))


def test_incorrect_inherited_triplet_blocked(tmp_path):
    f=create_inherited(tmp_path,'E_-4')
    f.write_text(f.read_text().replace('ϵ*','',1),encoding='utf-8')
    with pytest.raises(HardFlavorError):
        audit_inherited(f,'E_-4',inherited_prefactor(pole_report('E'),'E_-4'))


def test_new_model_selection_and_alphas():
    assert tuple(SUPPORTED)==MODELS
    assert alpha_argument('C_-3')=='m3'
    assert alpha_argument('E_+2')=='p2'
    assert select_models('B_-3,C_+1')==('B_-3','C_+1')
    assert select_fixed_models('B_-3,C_+1')==('B_-3','C_+1')
    with pytest.raises(Exception): select_models('B_-3,B_-3')
    with pytest.raises(Exception): select_fixed_models('D_+2')


def test_batch5_fixed_order_cli_dryrun_resume(tmp_path):
    from scripts.RunT3OrdinarySixFixedOrder import main
    model='B_-3'; tag='B_m3'
    hard_root=tmp_path/'hardroot';hard_root.mkdir()
    example=create_hard(tmp_path,model)
    example.rename(hard_root/tag)
    running=tmp_path/'running'/tag;running.mkdir(parents=True)
    (running/'running.json').write_text(json.dumps(running_doc(model)),encoding='utf-8')
    archive=tmp_path/'archive'/'T3_B_alpha_m3'/'data';archive.mkdir(parents=True)
    p=create_inherited(tmp_path,model);p.rename(archive/'eft1_after_F_one_loop_inputform.txt')
    (archive/'c5_pole_rge_consistency.json').write_text(json.dumps(pole_report('B')))
    output=tmp_path/'results'
    args=['--models',model,'--hard-root',str(hard_root),
          '--running-root',str(tmp_path/'running'),
          '--archive-root',str(tmp_path/'archive'),
          '--reference-root',str(tmp_path/'missing-reference'),
          '--output-root',str(output),'--run-id','smoke']
    assert main(args+['--check-only'])==0
    assert main(args)==0
    result=output/'smoke'/tag/'result.json'
    assert json.loads(result.read_text())['status']=='InheritedHardAndCommonCouplingPoleCancelled_Conditional'
    assert main(args+['--resume'])==0


@pytest.mark.parametrize('model',MODELS)
def test_matching_preflight_respects_exact_wilson_and_extra_quartic_inventory(tmp_path,model):
    from scripts.RunT3GeneralOrdinaryScalarHardBatch5 import preflight
    tag=model.replace('+','p').replace('-','m')
    root=tmp_path/'Ffirst'
    base=root/tag
    (base/'data').mkdir(parents=True)
    n={'B':2,'C':4,'E':4}[model[0]]
    extra={'B':1,'C':1,'E':0}[model[0]]
    (base/'native_status.json').write_text(json.dumps({
       'model':model,'status':'FThresholdTreeWilsonExported','integrated_fields':['F'],
       'light_bsm_fields':['S1','S2'],'heavy_flavour_count':3,'loop_order':0,
       'operator_basis_is_full_general_requested':True,
       'full_general_extra_yukawa_count':0,'full_general_quartic_count':extra,
       'full_general_quadratic_count':0}))
    (base/'status.json').write_text(json.dumps({'model':model,
       'status':'FermionThresholdTreeWilsonExported'}))
    common={'IntegratedFields':['F'],'ActiveHeavyFields':['S1','S2']}
    terms=[{'Chirality':'PL' if k%2==0 else 'PR'} for k in range(n)]
    (base/'data/eft1_after_F_wilson_seed.json').write_text(json.dumps(dict(common,
       TreeWilsonTermCount=n,TreeWilsonTerms=terms)))
    (base/'data/eft1_after_F_scalar_quartic_seed.json').write_text(json.dumps(dict(common,
       ScalarQuarticTermCount=1)))
    (base/'data/eft1_after_F_rge_export.json').write_text(json.dumps(common))
    mass='Coupling[MF, {Index[x, NFlavor]}, 0]'
    tree=('Field[NewScalar1, Scalar, {}, {}]*Field[NewScalar2, Scalar, {}, {}]+'
          'Coupling[MS1, {}, 1]+Coupling[MS2, {}, 1]+'
          '+'.join(mass for _ in range(n)))
    (base/'data/eft1_after_F_tree_inputform.txt').write_text(tree)
    wl=tmp_path/'batch5.wl';wl.write_text('dummy stage source')
    good=preflight(model,root,wl)
    assert good['tree_wilson_terms']==n
    assert good['additional_quartics_retained']==extra
    assert good['indexed_fermion_mass_factors']==n
    native=json.loads((base/'native_status.json').read_text())
    native['full_general_quartic_count']=99
    (base/'native_status.json').write_text(json.dumps(native))
    with pytest.raises(Exception,match='inventory'):
        preflight(model,root,wl)
