"""Stage 27 exact native Matchete six-case triplet witness (not certification)."""
from copy import deepcopy
import json
from pathlib import Path

import pytest
import sympy as sp

from scripts.TripletIntertwiner import (
    ALL, SIMPLE, PAULI, REFERENCE, audit, entry, generators, run,
    simple_cg, tensor_for, triple_contraction, zero,
)

ROOT=Path(__file__).resolve().parents[2]
FILES=ROOT/'output/general_t3_stage22'

def evidence():
    if not (FILES/'matchete_tensors.json').is_file() or not (FILES/'rgbeta_native_invariants.json').is_file():
        pytest.skip('Requires native Stage22 tensor JSON inputs')
    folder=FILES
    return (json.loads((folder/'matchete_tensors.json').read_text(encoding='utf-8-sig')),
            json.loads((folder/'rgbeta_native_invariants.json').read_text(encoding='utf-8-sig')))


def test_reference_symmetric_basis_orthonormal():
    assert all(zero(sum(sp.conjugate(REFERENCE[a][i,j])*REFERENCE[b][i,j]
                        for i in range(2) for j in range(2)) - (a==b))
               for a in range(3) for b in range(3))


def test_spin_half_commutator():
    assert zero(PAULI[0]*PAULI[1]-PAULI[1]*PAULI[0]-sp.I*PAULI[2])


def test_native_six_triplet_witness():
    m,r=evidence();result=audit(m,r)
    assert len(result['six_triplet_checks'])==6
    assert result['one_global_triplet_intertwiner_candidate'] is True
    assert result['native_rgbeta_s2_generator_components_exported'] is False
    assert result['physical_coupling_mapping_certified'] is False
    assert all(x['all_components_exact'] for x in result['six_triplet_checks'])


def test_extracted_unitary_is_nontrivial():
    m,r=evidence();record=tensor_for(SIMPLE[0],m,r)[1]
    cg,norm=simple_cg(record);U, reference, matchete = generators(cg)
    assert zero(U*U.H-sp.eye(3))
    assert not zero(U-sp.eye(3))
    assert zero(U*U.T-sp.diag(1,1,-1))
    assert all(zero(matchete[i]*U-U*reference[i]) for i in range(3))


def test_tampered_triplet_matchete_component_rejected():
    m,r=evidence();bad=deepcopy(m)
    label='D_+0';rec=tensor_for(label,bad,r)[1]
    rec['components'][0]='2'
    with pytest.raises(ValueError,match='Matchete|CG|norm|common'):
        audit(bad,r)


def test_model_specific_triplet_rephasing_rejected():
    m,r=evidence();bad=deepcopy(m)
    rec=tensor_for('A_-2',bad,r)[1]
    rec['components']=[('-('+z+')' if z!='0' else '0') for z in rec['components']]
    with pytest.raises(ValueError,match='common'):
        audit(bad,r)


def test_native_triple_invariant_changed_rejected():
    m,r=evidence();bad=deepcopy(r)
    rec=tensor_for('E_+0',m,bad)[0]
    rec['native_invariant_expression_inputform']=rec['native_invariant_expression_inputform'].replace('eps[','del[',1)
    with pytest.raises(ValueError,match='E-series'):
        audit(m,bad)


def test_e_series_component_rejected():
    m,r=evidence();bad=deepcopy(m)
    rec=tensor_for('E_+0',bad,r)[1]
    idx=next(k for k,z in enumerate(rec['components']) if z!='0')
    rec['components'][idx]='0'
    with pytest.raises(ValueError,match='disagrees|missing|Zero'):
        audit(bad,r)


def test_reference_cg_disallows_scaled_singular_basis():
    m,r=evidence();rec=tensor_for('A_+0',m,r)[1]; cg,_=simple_cg(rec)
    broken=list(cg);broken[0]=broken[1]
    with pytest.raises(ValueError,match='unitary'):
        generators(broken)


def test_generators_reject_wrong_local_cg_sign():
    m,r=evidence();rec=tensor_for('A_+0',m,r)[1]; cg,_=simple_cg(rec)
    bad=list(cg);bad[0]=2*bad[0]
    with pytest.raises(ValueError,match='unitary'):
        generators(bad)


def test_runtime_fail_closed_without_stage26(tmp_path):
    m,r=evidence()
    (tmp_path/'matchete_tensors.json').write_text(json.dumps(m),encoding='utf-8')
    (tmp_path/'rgbeta_native_invariants.json').write_text(json.dumps(r),encoding='utf-8')
    with pytest.raises(ValueError,match='Stage26'):
        run(tmp_path,None)


def test_full_local_integration_if_available():
    if not (FILES/'stage26_conditional_quartic_action.json').is_file() or \
       not (FILES/'stage25_live_native_quartic_assembly.json').is_file():
        pytest.skip('User-generated live Stage25/26 evidence not present')
    d=run(FILES,None)
    assert d['stage26_fresh_evidence_checked'] is True
    assert d['source_native_matchete_tensor_count']==6
    assert d['physical_coupling_mapping_certified'] is False
