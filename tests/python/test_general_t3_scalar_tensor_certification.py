"""Stage 22 algebra regression tests (run without Mathematica / RGBeta)."""
from __future__ import annotations

from copy import deepcopy
import pytest
from sympy import Matrix

from scripts.RGBetaTensorAlgebra import (TensorCertificateError, exact,
                                          solve_basis, potential_map, same_species_symmetry)
from scripts.RunGeneralT3ScalarTensorCertification import certify


def t(symbol: str, values: list[str], counts=None, dims=None):
    return dict(symbol=symbol,components=values,
                counts=counts or [0,0,0,0,0,4], dimensions=dims or [1,1,1,1],
                potential_prefactor='1', species=['S2dag']*4)


def test_exact_scalars_and_safe_parser():
    assert exact('Sqrt[2]/2')**2==exact('1/2')
    assert exact('Complex[0, 1]/2')==exact('I/2')
    assert exact('Rational[-3,4]')==exact('-3/4')
    for illegal in ('__import__("os")','x','Float[0.3]', '1.2', 'open("secret")'):
        with pytest.raises(TensorCertificateError): exact(illegal)


def test_one_dimensional_normalisation_and_sign():
    block=solve_basis([t('m',['1'])],[t('r',['-2'])])
    assert block['r_in_m']==[['-2']]
    p=potential_map(block,['1'],['-1/2'])
    assert p['rgbeta_from_matchete']==[['1']]
    assert p['zero_coupling_recovery']


def test_two_invariant_mixing_exact_reconstruction():
    # 4 distinguishable bosons, 16 tensor components. The 2D subspace is
    # mixed with a nontrivial integer GL(2,Q) matrix.
    counts=[1,1,1,1,0,0]
    dim=[2,2,2,2]
    def sample(n,v):
        return dict(symbol=n,counts=counts,dimensions=dim,
                    components=[str(x) for x in v],potential_prefactor='1')
    e1=[1]+[0]*15
    e2=[0,1]+[0]*14
    M=[sample('m1',e1),sample('m2',e2)]
    R=[sample('r1',[a+2*b for a,b in zip(e1,e2)]),
       sample('r2',[3*a+7*b for a,b in zip(e1,e2)])]
    block=solve_basis(M,R)
    assert block['r_in_m']==[['1','3'],['2','7']]
    assert block['m_in_r']==[['7','-3'],['-2','1']]
    assert block['maximum_exact_component_residual']=='0'
    p=potential_map(block,['1','2'],['3','4'], 'conjugate')
    assert p['orientation']=='conjugate'
    assert p['rgbeta_from_matchete']==[['7/3','-2'],['-1/2','1/2']]


def test_rank_deficiency_and_mismatch_fail_closed():
    m1=t('m',['1']); m2=t('m2',['1'])
    with pytest.raises(TensorCertificateError,match='Rank deficient'):
        solve_basis([m1,m2],[t('r',['1']),t('r2',['2'])])
    with pytest.raises(TensorCertificateError,match='Incompatible'):
        solve_basis([m1],[t('r',['1'],counts=[0,0,0,4,0,0])])
    with pytest.raises(TensorCertificateError,match='expected 1 components'):
        solve_basis([m1],[t('r',['1','1'])])
    with pytest.raises(TensorCertificateError,match='Missing any potential prefactor'):
        potential_map(solve_basis([m1],[t('r',['1'])]),[],['1'])


def test_bose_permutation_catches_asymmetry():
    record={'symbol':'m','counts':[0,2,0,2,0,0],
            'dimensions':[2,2,2,2], 'components':['0']*16}
    record['components'][1]='1'
    with pytest.raises(TensorCertificateError,match='nonsymmetric'):
        same_species_symmetry(record)


def test_missing_native_16_model_data_does_not_map_to_identity():
    # No native records: all 1-slot test couplings explicitly uncertified.
    dummy={k:{'new_rgbeta_symbols':[{'name':'r1','counts':[0,0,0,0,0,4]}],
              'preexisting_rgbeta_symbols':[]} for k in ('A_+0','B_-1')}
    reg={'models':dummy}
    out,ok=certify(reg,{}, {})
    assert not ok and out['uncertified_slots']==2
    assert all(x['status']=='Uncertified' for m in out['models'].values() for x in m['blocks'])


def test_fake_candidate_rgbeta_provenance_rejected():
    from scripts.RunGeneralT3ScalarTensorCertification import tensor_index
    datum=t('r',['1'])
    candidate={'format':'stage22-exact-tensors-v1','basis':'RGBeta',
               'source':'GeneratedCandidateBasisNotWolframVerified',
               'models':{'A_+0':{'tensors':[datum]}}}
    with pytest.raises(TensorCertificateError,match='cannot use candidate'):
        tensor_index(candidate,'A_+0','RGBeta',True)


def test_tensor_equal_but_potential_prefactor_missing():
    from scripts.RunGeneralT3ScalarTensorCertification import tensor_index
    rec=t('m',['1']);rec.pop('potential_prefactor')
    data={'format':'stage22-exact-tensors-v1','source':'MatcheteInvariantTensors',
          'basis':'Matchete','models':{'A_+0':{'tensors':[rec]}}}
    with pytest.raises(TensorCertificateError,match='prefactor'):
        tensor_index(data,'A_+0','Matchete',False)


def test_native_two_channel_model_certification_and_sign_provenance():
    counts=[1,1,1,1,0,0]
    species=['H','Hdag','S1','S1dag']
    def rec(symbol,vec,factor):
        return {'symbol':symbol,'counts':counts,'species':species,
                'dimensions':[2,2,2,2], 'components':[str(x) for x in vec],
                'potential_prefactor':factor}
    e1=[1]+[0]*15
    e2=[0,1]+[0]*14
    m=[rec('lambdaT3ExtraQ111100I1',e1,'1'),
       rec('lambdaT3ExtraQ111100I2',e2,'1')]
    r=[rec('r1',[a+2*b for a,b in zip(e1,e2)],'-1'),
       rec('r2',[3*a+7*b for a,b in zip(e1,e2)],'-1')]
    common={'format':'stage22-exact-tensors-v1','models':{
        'B_-1':{'holomorphic_orientation':'direct'}}}
    mat=deepcopy(common);mat['basis']='Matchete';mat['source']='MatcheteInvariantTensors'
    mat['models']['B_-1']['tensors']=m
    rgb=deepcopy(common);rgb['basis']='RGBeta';rgb['source']='RGBetaNativeQuartic'
    rgb.update(native_invariant_definitions_verified=True,
               native_potential_prefactors_verified=True,
               field_component_alignment_verified=True)
    rgb['models']['B_-1']['tensors']=r
    reg={'models':{'B_-1':{'new_rgbeta_symbols':[
        {'name':'r1','counts':counts},{'name':'r2','counts':counts}],
        'preexisting_rgbeta_symbols':[]}}}
    report,ok=certify(reg,mat,rgb)
    assert ok and report['uncertified_slots']==0 and report['exact_algebra_slots']==2
    block=report['models']['B_-1']['blocks'][0]
    assert block['tensor_map']['r_in_m']==[['1','3'],['2','7']]
    assert block['potential_map']['rgbeta_from_matchete']==[['-7','3'],['2','-1']]
    # Changing a native tensor without updating its partner cannot silently pass.
    rgb['models']['B_-1']['tensors'][0]['components'][2]='1'
    bad,ok=certify(reg,mat,rgb)
    assert not ok and bad['uncertified_slots']==2


def test_singlet_axes_are_explicit_singletons_and_not_su2_indices():
    """A_+0 starts with a 1x1x3x3 tensor, but Matchete stores 3x3."""
    from scripts.RunGeneralT3ScalarTensorCertification import tensor_index
    counts=[0,0,0,2,1,1]
    rec={'symbol':'lambdaT3ExtraQ000211I1',
         'species':['S1dag','S1dag','S2','S2dag'],
         'counts':counts,'dimensions':[1,1,3,3],
         'components':['1']+['0']*8,
         'potential_prefactor':'1',
         'native_non_singlet_dimensions':[3,3],
         'native_tensor_dimensions':[3,3]}
    payload={'format':'stage22-exact-tensors-v1',
             'basis':'Matchete','source':'MatcheteInvariantTensors',
             'models':{'A_+0':{'tensors':[rec]}}}
    indexed=tensor_index(payload,'A_+0','Matchete',False)
    assert len(indexed[tuple(counts)])==1
    assert solve_basis([rec],[dict(rec, symbol='rgbeta')])['r_in_m']==[['1']]


def test_incomplete_matchete_export_refused(tmp_path):
    import json
    from scripts.RunGeneralT3ScalarTensorCertification import audit_native_matchete
    fp=tmp_path/'matchete.json'
    fp.write_text(json.dumps({'format':'stage22-exact-tensors-v1',
        'source':'MatcheteInvariantTensors','basis':'Matchete',
        'models':{'A_+0':{'tensors':[]}}}))
    with pytest.raises(TensorCertificateError, match='exactly the 16'):
        audit_native_matchete({'models':{'A_+0':{},'B_-1':{}}}, fp)


def test_rgbeta_raw_expressions_are_diagnostic_only(tmp_path):
    import json
    from scripts.RunGeneralT3ScalarTensorCertification import audit_native_rgbeta_invariants
    reg={'models':{}}
    result={'format':'stage22-native-rgbeta-invariants-v1',
        'source':'RGBetaRegisteredQuarticInvariantFunctions',
        'status':'NativeInvariantExpressionsNotCommonBasisTensors',
        'native_invariant_component_alignment_verified':False,
        'native_potential_prefactors_verified':False,'models':{}}
    for i in range(16):
        model=f'M{i}'
        symbols=[{'counts':[0,0,0,0,0,4], 'name':f'r{i}_{j}'}
                 for j in range(4 if i<2 else 3)]
        reg['models'][model]={'new_rgbeta_symbols':symbols,
                              'preexisting_rgbeta_symbols':[]}
        result['models'][model]={'invariants':[
            {'symbol':v['name'], 'native_fields_inputform':'{S2,S2,S2,S2}',
             'native_invariant_function_inputform':'Function[1]',
             'native_invariant_expression_inputform':'1',
             'native_coupling_function_inputform':'Function[6*lambda]',
             'native_unique_arrangements_inputform':'{}',
             'self_conjugate':False,'component_tensor_exported':False,
             'verified_potential_prefactor':None}
            for v in symbols]}
    fp=tmp_path/'rgbeta_expressions.json'
    fp.write_text(json.dumps(result))
    report=audit_native_rgbeta_invariants(reg,fp)
    assert report['quartic_slots_checked']==50
    assert report['component_tensors_exported'] is False
    result['models']['M0']['invariants'][0]['verified_potential_prefactor']='1'
    fp.write_text(json.dumps(result))
    with pytest.raises(TensorCertificateError, match='fail-closed'):
        audit_native_rgbeta_invariants(reg,fp)


def test_zero_coupling_recovery_not_a_proof_of_native_normalisation():
    # A zero initial coupling maps to zero for every invertible change of basis;
    # this must not be mistaken for evidence of native CG normalisation.
    block=solve_basis([t('m',['2'])],[t('r',['5'])])
    assert block['zero_coupling_recovery']
    assert block['tensor_equality']
    assert block['r_in_m']==[['5/2']]


def test_missing_key_records_in_uploaded_rgbeta_evidence_are_rejected(tmp_path):
    """Observed Windows export: unresolved Missing[] must not pass as metadata."""
    import json
    from scripts.RunGeneralT3ScalarTensorCertification import (audit_native_declarations,
                                               audit_native_rgbeta_invariants)
    reg={'models':{f'M{i}':{'new_rgbeta_symbols':[
        {'name':f'q{i}_{j}','counts':[0,0,0,0,0,4]}
        for j in range(4 if i<2 else 3)],
        'preexisting_rgbeta_symbols':[]} for i in range(16)}}
    missing='Missing["KeyAbsent", Invariant]'
    d={'format':'stage22-native-rgbeta-declarations-v1',
       'source':'RGBetaRegisteredQuarticAssociation',
       'status':'NativeDeclarationMetadataOnly',
       'matchete_rgbeta_normalisation_verified':False,
       'models':{}}
    inv={'format':'stage22-native-rgbeta-invariants-v1',
       'source':'RGBetaRegisteredQuarticInvariantFunctions',
       'status':'NativeInvariantExpressionsNotCommonBasisTensors',
       'native_invariant_component_alignment_verified':False,
       'native_potential_prefactors_verified':False,'models':{}}
    for model, entry in reg['models'].items():
        syms=[v['name'] for v in entry['new_rgbeta_symbols']]
        d['models'][model]={'declarations':[
            {'symbol':sym,'fields_inputform': 'Missing["KeyAbsent", Fields]',
             'invariant_inputform':missing,'coupling_factor_inputform':missing,
             'unique_arrangements_inputform':missing,'self_conjugate':False}
            for sym in syms]}
        inv['models'][model]={'invariants':[
            {'symbol':sym,'native_fields_inputform':'Missing["KeyAbsent", Fields]',
             'native_invariant_function_inputform':missing,
             'native_invariant_expression_inputform':missing+'[t3i1,t3i2,t3i3,t3i4]',
             'native_coupling_function_inputform':missing,
             'native_unique_arrangements_inputform':missing,
             'self_conjugate':False,'component_tensor_exported':False,
             'verified_potential_prefactor':None}
            for sym in syms]}
    f=tmp_path/'test.json'
    f.write_text(json.dumps(d))
    with pytest.raises(TensorCertificateError,match='unresolved native fields_inputform'):
        audit_native_declarations(reg,f)
    f.write_text(json.dumps(inv))
    with pytest.raises(TensorCertificateError,match='unresolved native native_fields_inputform'):
        audit_native_rgbeta_invariants(reg,f)


def test_native_diagnostic_metadata_must_contain_resolved_fields(tmp_path):
    """Reject arbitrary nonmissing text, not just Missing[] values."""
    import json
    from scripts.RunGeneralT3ScalarTensorCertification import audit_native_declarations
    reg={'models': {f'M{i}':{'new_rgbeta_symbols':[
        {'name':f'q{i}_{j}','counts':[0,0,0,0,0,4]}
        for j in range(4 if i<2 else 3)], 'preexisting_rgbeta_symbols': []}
        for i in range(16)}}
    data={'format':'stage22-native-rgbeta-declarations-v1',
          'source':'RGBetaRegisteredQuarticAssociation',
          'status':'NativeDeclarationMetadataOnly',
          'matchete_rgbeta_normalisation_verified':False,'models':{}}
    for model, entry in reg['models'].items():
        data['models'][model]={'declarations': [
            {'symbol':v['name'], 'fields_inputform':'some random value',
             'invariant_inputform':'1 &','coupling_factor_inputform':'6*lambda &',
             'unique_arrangements_inputform':'{}','self_conjugate':False}
            for v in entry['new_rgbeta_symbols']]}
    f=tmp_path/'badfields.json'; f.write_text(json.dumps(data))
    with pytest.raises(TensorCertificateError,match='explicit four-field list'):
        audit_native_declarations(reg,f)


def test_rgbeta_native_declaration_invariant_coherence_guard(tmp_path):
    import json
    from scripts.RunGeneralT3ScalarTensorCertification import audit_native_rgbeta_coherence
    d={'models':{}}
    i={'models':{}}
    for k in range(16):
        model=f'M{k}'
        count=4 if k<2 else 3
        d['models'][model]={'declarations':[]}
        i['models'][model]={'invariants':[]}
        for j in range(count):
            sym=f'q{k}_{j}'
            d['models'][model]['declarations'].append({
                'symbol':sym,'fields_inputform':'{S2, S2, S2, S2}',
                'invariant_inputform':'(1 &)', 'coupling_factor_inputform':'(6*q &)',
                'unique_arrangements_inputform':'{{1,2,3,4}}','self_conjugate':False})
            i['models'][model]['invariants'].append({
                'symbol':sym,'native_fields_inputform':'{S2, S2, S2, S2}',
                'native_invariant_function_inputform':'(1 &)',
                'native_coupling_function_inputform':'(6*q &)',
                'native_unique_arrangements_inputform':'{{1,2,3,4}}',
                'self_conjugate':False})
    df=tmp_path/'decl.json'; fi=tmp_path/'inv.json'
    df.write_text(json.dumps(d)); fi.write_text(json.dumps(i))
    assert audit_native_rgbeta_coherence(df,fi)['coupling_slots_checked']==50
    i['models']['M0']['invariants'][0]['native_invariant_function_inputform']='(-1 &)'
    fi.write_text(json.dumps(i))
    with pytest.raises(TensorCertificateError,match='independent native exports disagree'):
        audit_native_rgbeta_coherence(df,fi)
