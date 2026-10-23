"""Pure symbolic closure tests. Full production tensor scan runs via CLI in the repository."""
import sympy as sp
import pytest

from RGE.running.intermediate.T3BOperatorClosure import (
    canonical, canonical_values, project_onto_two_templates,
)

S, K, a, b, epsilon = sp.symbols("S K a b epsilon", nonzero=True)
C12 = (1,1,5,7)
C5 = (1,1,1,3)
EXTRA = (1,1,5,5)
COMPONENTS = [C5, EXTRA, C12]
TEMPLATES = {"C12": {C12:S}, "C5": {C5:K}}


def test_closed_two_template_operator_map():
    beta = {C12:a*S, C5:b*K}
    factors, residual = project_onto_two_templates(beta, TEMPLATES, COMPONENTS)
    assert residual == {}
    assert sp.simplify(factors["C12"]-a) == 0
    assert sp.simplify(factors["C5"]-b) == 0


def test_new_operator_is_not_hidden_by_pivot_projection():
    beta = {C12:a*S, C5:b*K, EXTRA:epsilon*S}
    factors, residual = project_onto_two_templates(beta, TEMPLATES, COMPONENTS)
    assert len(residual) == 1
    assert sp.simplify(residual[EXTRA]-epsilon*S) == 0
    assert factors == {"C12":a,"C5":b}


def test_residual_within_LLSS_is_detected():
    beta = {C12:a*S+epsilon, C5:b*K}
    factors, residual = project_onto_two_templates(beta, TEMPLATES, COMPONENTS)
    assert residual == {}
    # A pivot alone cannot distinguish a shifted coefficient from a true
    # template term; a second LLSS tensor component is needed.
    T2 = (1,1,6,8)
    comps = [C5, EXTRA, C12, T2]
    seed = {"C12": {C12:S, T2:-S}, "C5":{C5:K}}
    _, residual = project_onto_two_templates(
        {C12:a*S+epsilon, T2:-a*S, C5:b*K}, seed, comps,
    )
    assert residual[T2] != 0


def test_external_component_cannot_be_ignored():
    X=(2,3,1,3)
    _, residual = project_onto_two_templates(
        {C12:a*S,C5:b*K,X:sp.Integer(1)},TEMPLATES,[*COMPONENTS,X],
    )
    assert residual[X] == 1


def test_pair_symmetry_conflict_is_rejected():
    with pytest.raises(ValueError,match="Inconsistent tensor symmetry"):
        canonical_values({(1,2,5,7):S, (2,1,7,5):S+1})


def test_disjoint_support_required():
    with pytest.raises(ValueError,match="disjoint supports"):
        project_onto_two_templates({}, {"C12":{C12:S},"C5":{C12:K}}, COMPONENTS)


def test_unscanned_template_is_rejected():
    with pytest.raises(ValueError, match="support is not in the scanned"):
        project_onto_two_templates({}, TEMPLATES, [C12,EXTRA])
