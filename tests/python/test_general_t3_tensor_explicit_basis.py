"""Stage 22 formal invariant comparison, not physical coupling certification."""
import json
from pathlib import Path
import pytest
import sympy as sp

from scripts.RGBetaExplicitBasis import (
    parse_fields, eval_component, reconstructed_tensor, verify_sources, compare,
)


def test_native_species_parser():
    assert parse_fields('{H, Bar[H], Bar[S1], Bar[S2]}') == ['H','Hdag','S1dag','S2dag']
    with pytest.raises(ValueError):
        parse_fields('{H, Bad[], S1, S2}')


def test_delta_epsilon_orientation():
    assert eval_component('eps[SU2L[fund], t3i1, t3i2]', (0,1,0,0)) == 1
    assert eval_component('eps[SU2L[fund], t3i1, t3i2]', (1,0,0,0)) == -1
    assert eval_component('del[SU2L[S2], t3i1, t3i2]', (2,2,0,0)) == 1
    assert eval_component('del[SU2L[S2], t3i1, t3i2]', (1,2,0,0)) == 0


def test_exact_symmetrised_doublet_projection():
    expr = ('(del[SU2L[fund], t3i1, t3i3]*del[SU2L[fund], t3i2, t3i4] + '
            'del[SU2L[fund], t3i1, t3i4]*del[SU2L[fund], t3i2, t3i3])/2')
    assert eval_component(expr,(0,1,0,1)) == sp.Rational(1,2)
    assert eval_component(expr,(0,0,0,0)) == 1


def test_triplet_bridge_fails_closed():
    with pytest.raises(ValueError, match='triplet CG'):
        eval_component('delS2[SU2L, t3i4, t3i1, t3i2]',(0,0,0,0))


def test_no_arbitrary_expression_execution():
    for expr in ('foo(t3i1)', '__import__("os").system("whoami")',
                 'del[SU2L[fund], t3i1, t3i2]+arbitraryFunction[x]'):
        with pytest.raises((ValueError, SyntaxError)):
            eval_component(expr,(0,0,0,0))


def test_end_to_end_exact_doublet_block():
    m = {'symbol':'lambdaT3ExtraQ022000I1',
         'species':['Hdag','Hdag','S1','S1'],
         'dimensions':[2,2,2,2], 'counts':[0,2,2,0,0,0],
         'components':[str(sp.Rational(int(i==k and j==l)+int(i==l and j==k),2))
                       for i in range(2) for j in range(2)
                       for k in range(2) for l in range(2)]}
    inv = {'symbol':'lambdaT3GeneralQ022000I1',
           'native_fields_inputform':'{Bar[H], Bar[H], S1, S1}',
           'native_invariant_expression_inputform':
           '(del[SU2L[fund], t3i1, t3i3]*del[SU2L[fund], t3i2, t3i4] + '
           'del[SU2L[fund], t3i1, t3i4]*del[SU2L[fund], t3i2, t3i3])/2'}
    from scripts.RGBetaTensorAlgebra import solve_basis
    assert solve_basis([m],[reconstructed_tensor(inv,m)])['r_in_m'] == [['1']]


def test_uploaded_files_when_present():
    root = Path(__file__).resolve().parents[2]
    p = root/'output/general_t3_stage22'
    required = ['matchete_tensors.json','rgbeta_native_declarations.json',
                'rgbeta_native_invariants.json']
    if not all((p/n).is_file() for n in required):
        pytest.skip('Native evidence not on this machine')
    m,d,i = (json.loads((p/n).read_text(encoding='utf-8-sig')) for n in required)
    assert verify_sources(m,i,d)
    result = compare(m,i,d)
    counts = result['totals']
    assert counts['slots'] == 50
    assert counts['blocks'] == 46
    # Corrected native B/C expressions may have zero mismatch slots.
    exact_count = counts.get('ExactGeometryCandidate',0)
    mismatch_count = counts.get('ExplicitGeometryMismatch',0)
    triplet_count = counts.get('TripletBasisPending',0)
    assert triplet_count == 6
    assert exact_count + mismatch_count + triplet_count == 50
    assert exact_count in (32,36,40,44)
    assert result['models']['E_+0'][0]['basis_change_r_in_m'] == [
        ['1/2','1'], ['sqrt(33)/6','0']]
    assert result['physical_normalisation_verified'] is False
    assert result['rgbeta_lagrangian_prefactor_verified'] is False
