"""Stage26 exact scalar-action algebra; fixtures are NOT live certification."""
from copy import deepcopy
from pathlib import Path

import pytest
import sympy as sp

from scripts.QuarticActionIdentity import (
    actions_equal, assemble_action, audit_block, hermitian_conjugate,
    scalar_action_conversion, tensor_polynomial, verify_matchete_source,
    run,
)

ROOT = Path(__file__).resolve().parents[2]


def toy(species, dims, components, name='x'):
    return {'symbol':name, 'species':species, 'dimensions':dims,
            'components':list(components), 'counts':[0]*6}


def block(A, mat, rgb, cls='DirectGeometryCandidate'):
    matrix = sp.Matrix(A)
    minus_inverse = -matrix.inv()
    return {
        'model':'toy','block':'test',
        'matchete_couplings':[x['symbol'] for x in mat],
        'rgbeta_couplings':[x['symbol'] for x in rgb],
        'geometry_class':cls,
        'r_in_m':[[str(matrix[i,j]) for j in range(matrix.cols)] for i in range(matrix.rows)],
        'lambda_rgbeta_from_lambda_matchete_if_assumptions_hold':[
            [str(minus_inverse[i,j]) for j in range(minus_inverse.cols)]
            for i in range(minus_inverse.rows)],
        'assumed_relative_potential_sign':'-1',
        'component_prefactor_one_ordered_orientation':'1/4',
        'physical_mapping_certified':False,
    }


def test_exact_normalisation_conventions():
    c = scalar_action_conversion()
    assert c['net_monomial_action_factor'] == '1'
    assert c['rgbeta_quartic_lagrangian_sign'] == '-1'


def test_identical_bosons_collect_without_factorial():
    t=toy(['H','H','S1dag','S2'],[2,2,1,1], ['1','1','1','1'])
    p=tensor_polynomial(t)
    assert len(p)==3
    assert sorted(p.values())==[1,1,2]


def test_explicit_hermitian_conjugate_includes_complex_coefficient():
    t=toy(['Hdag','Hdag','S1','S2'],[1,1,1,1],['I'])
    z=sp.Symbol('z',complex=True)
    a=assemble_action([t],[z],plus_hc=True)
    assert a[(('Hdag',0),('Hdag',0),('S1',0),('S2',0))]==sp.I*z
    assert a[(('H',0),('H',0),('S1dag',0),('S2dag',0))]==-sp.I*sp.conjugate(z)


def test_formal_conjugation_round_trip():
    t=toy(['Hdag','Hdag','S1','S2'],[1,1,1,1],['1+I'])
    assert hermitian_conjugate(hermitian_conjugate(tensor_polynomial(t)))==tensor_polynomial(t)


def test_one_dimensional_signed_candidate_identifies_full_action():
    m=toy(['Hdag','Hdag','S1','S2'],[1,1,1,1],['2'],'m')
    r=toy(['Hdag','Hdag','S1','S2'],[1,1,1,1],['2'],'r')
    import scripts.QuarticActionIdentity as stage26
    with pytest.MonkeyPatch.context() as mp:
        mp.setattr(stage26,'reconstruct',lambda native,ref,cls:native)
        result=audit_block(block([[1]],[m],[r]),{'m':m},{'r':r})
    assert result['exact_signed_action_identity_in_provisional_basis']
    assert result['physical_mapping_certified'] is False


def test_two_coupling_signed_map_has_no_factor_4():
    # Two independent symmetric doublet tensors, full component source.
    m1=toy(['Hdag','Hdag','S1','S2'],[2,2,1,1],['1','0','0','0'],'m1')
    m2=toy(['Hdag','Hdag','S1','S2'],[2,2,1,1],['0','0','0','1'],'m2')
    r1=toy(['Hdag','Hdag','S1','S2'],[2,2,1,1],['2','0','0','0'],'r1')
    r2=toy(['Hdag','Hdag','S1','S2'],[2,2,1,1],['0','0','0','3'],'r2')
    import scripts.QuarticActionIdentity as stage26
    with pytest.MonkeyPatch.context() as mp:
        mp.setattr(stage26,'reconstruct',lambda native,ref,cls:native)
        result=audit_block(block([[2,0],[0,3]],[m1,m2],[r1,r2]),
                           {m['symbol']:m for m in [m1,m2]},
                           {r['symbol']:r for r in [r1,r2]})
    assert result['quartic_couplings']==2
    assert result['full_action_monomials_matchete']==4


def test_wrong_signed_mapping_is_rejected():
    m=toy(['Hdag','Hdag','S1','S2'],[1,1,1,1],['1'],'m')
    r=toy(['Hdag','Hdag','S1','S2'],[1,1,1,1],['1'],'r')
    b=block([[1]],[m],[r]);b['lambda_rgbeta_from_lambda_matchete_if_assumptions_hold']=[['1']]
    import scripts.QuarticActionIdentity as stage26
    with pytest.MonkeyPatch.context() as mp:
        mp.setattr(stage26,'reconstruct',lambda native,ref,cls:native)
        with pytest.raises(ValueError,match='Signed block map'):
            audit_block(b,{'m':m},{'r':r})


def test_component_mutation_is_rejected():
    m=toy(['Hdag','Hdag','S1','S2'],[1,1,1,1],['1'],'m')
    r=toy(['Hdag','Hdag','S1','S2'],[1,1,1,1],['3'],'r')
    import scripts.QuarticActionIdentity as stage26
    with pytest.MonkeyPatch.context() as mp:
        mp.setattr(stage26,'reconstruct',lambda native,ref,cls:native)
        with pytest.raises(ValueError,match='action identity FAILED'):
            audit_block(block([[1]],[m],[r]),{'m':m},{'r':r})


def test_raw_quartic_multiplicity_negative_control():
    m=toy(['H','H','S1dag','S2'],[2,2,1,1],['1','1','1','1'],'m')
    p=tensor_polynomial(m)
    assert any(val==2 for val in p.values())
    # A wrong coefficient at one ordered tensor component changes an entire
    # commuting monomial coefficient; no implicit factorial may absorb it.
    v=toy(['H','H','S1dag','S2'],[2,2,1,1],['1','0','1','1'],'v')
    assert not actions_equal(assemble_action([m],[1]),assemble_action([v],[1]))


def test_unreviewed_geometry_class_refused():
    m=toy(['Hdag','Hdag','S1','S2'],[1,1,1,1],['1'],'m')
    r=toy(['Hdag','Hdag','S1','S2'],[1,1,1,1],['1'],'r')
    with pytest.raises(ValueError,match='Uncertified geometry class'):
        audit_block(block([[1]],[m],[r],cls='InventedBasis'),{'m':m},{'r':r})


def test_project_source_validation_if_present():
    if not (ROOT/'Lagrangian/model/LagrangianBuilder.wl').is_file():
        pytest.skip('Requires local project Matchete source files')
    assert verify_matchete_source(ROOT)


def test_live_stage26_optional_integration():
    path=ROOT/'output/general_t3_stage22/stage25_live_native_quartic_assembly.json'
    if not path.is_file():
        pytest.skip('Requires user-generated Stage25 live Mathematica evidence')
    result=run(ROOT/'output/general_t3_stage22',None,ROOT)
    assert result['slots_checked']==50
    assert result['blocks_checked']==46
    assert result['physical_coupling_mapping_certified'] is False
