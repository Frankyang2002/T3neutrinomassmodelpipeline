"""B/C-only conditional signed HHSS map; no physical basis certification."""
from __future__ import annotations

import ast
from pathlib import Path

import numpy as np
import pytest

from Numerical.running.GeneralT3SequentialLeadingLog import TensorBundle
from Numerical.running.GeneralT3ConditionalBCQuarticMap import (
    CHANNEL_NAMES, DIRECT_QUARTICS, CONVENTION,
    apply_conditional_bc_direct_quartic_map,
)


def candidate(model: str = 'B_-1', *, names: tuple[str, ...] | None = None):
    if names is None:
        names = ('native:lambdaT3', CHANNEL_NAMES[0],
                 'portal:lambdaH1Inv1', CHANNEL_NAMES[1],
                 'general:lambdaOther')
    a = np.zeros((len(names), 4, 4, 8, 8), dtype=float)
    for k in range(len(names)):
        for i in range(4):
            for j in range(4):
                a[k, i, j, :, :] = np.eye(8)*(k+1)*0.013*(1+(i==j))
    entries=[]
    for name in names:
        entry={'name':name.split(':',1)[-1],
               'origin': 'exact SU2 reference invariant basis' if name.startswith('general:')
               else 'native Matchete source CG'}
        entries.append(entry)
    return TensorBundle(model,np.diag(np.arange(1,9,dtype=float)),
                        a.sum(axis=0),a,names,
                        {'quartic_channels':entries,
                         'matchete_rgbeta_map_certified':False})


@pytest.mark.parametrize('model',['B_-1','C_-1'])
def test_sign_flip_only_two_selected_tensor_channels(model):
    original=candidate(model)
    initial=original.channels.copy()
    endpoint_state={'lambdaHHdagS2S2':0.01,'lambdaHHdagS1barS1bar':-0.02}
    changed=apply_conditional_bc_direct_quartic_map(model,original)
    assert changed is not original
    assert changed.channels is not original.channels
    np.testing.assert_array_equal(original.channels,initial)
    np.testing.assert_array_equal(changed.mass2,original.mass2)
    for k,name in enumerate(original.names):
        sign=-1 if name in CHANNEL_NAMES else 1
        np.testing.assert_allclose(changed.channels[k],sign*initial[k],rtol=0,atol=0)
    np.testing.assert_allclose(changed.quartic,changed.channels.sum(axis=0),rtol=0,atol=0)
    np.testing.assert_allclose(changed.quartic,
        original.quartic-2*sum((initial[original.names.index(n)] for n in CHANNEL_NAMES)),
        rtol=1e-14,atol=1e-14)
    assert endpoint_state=={'lambdaHHdagS2S2':0.01,'lambdaHHdagS1barS1bar':-0.02}
    assert 'bc_direct_quartic_convention' not in original.provenance
    assert changed.provenance['bc_direct_quartic_convention']==CONVENTION
    assert changed.provenance['bc_direct_quartics_signed']==list(DIRECT_QUARTICS)
    assert changed.provenance['bc_direct_quartic_map_affects_rges'] is False
    assert changed.provenance['bc_direct_quartic_map_physically_certified'] is False
    selected=[x for x in changed.provenance['quartic_channels']
              if x['name'] in DIRECT_QUARTICS]
    assert len(selected)==2 and all(
        x['conditional_matchete_over_rgbeta_sign']==-1 for x in selected)
    assert all('conditional_matchete_over_rgbeta_sign' not in x for x in
               original.provenance['quartic_channels'])


@pytest.mark.parametrize('model',['A_-4','C_+1','E_-2'])
def test_nonmixed_models_cannot_apply_bc_sign_map(model):
    with pytest.raises(ValueError,match='restricted to B_-1/C_-1'):
        apply_conditional_bc_direct_quartic_map(model,candidate(model))


def test_bundle_model_mismatch_rejected():
    with pytest.raises(ValueError,match='restricted to B_-1/C_-1'):
        apply_conditional_bc_direct_quartic_map('B_-1',candidate('C_-1'))


def test_missing_or_duplicate_direct_channels_rejected():
    names=('native:lambdaT3',CHANNEL_NAMES[0],'general:lambdaOther')
    with pytest.raises(ValueError,match='missing or duplicated'):
        apply_conditional_bc_direct_quartic_map('B_-1',candidate(names=names))
    names=('native:lambdaT3',CHANNEL_NAMES[0],CHANNEL_NAMES[0],CHANNEL_NAMES[1])
    with pytest.raises(ValueError,match='duplicate HHSS'):
        apply_conditional_bc_direct_quartic_map('C_-1',candidate('C_-1',names=names))


def test_map_is_not_silently_applied_twice():
    signed=apply_conditional_bc_direct_quartic_map('B_-1',candidate())
    with pytest.raises(ValueError,match='already applied'):
        apply_conditional_bc_direct_quartic_map('B_-1',signed)


def test_inconsistent_initial_quartic_total_rejected():
    from dataclasses import replace
    original=candidate()
    broken=replace(original,quartic=original.quartic+0.1)
    with pytest.raises(ValueError,match='not the sum'):
        apply_conditional_bc_direct_quartic_map('B_-1',broken)


def test_non_native_mapping_guard():
    from dataclasses import replace
    original=candidate()
    updated=dict(original.provenance)
    updated['quartic_channels']=[dict(x) for x in original.provenance['quartic_channels']]
    updated['quartic_channels'][1]['origin']='native Matchete source CG'
    with pytest.raises(ValueError,match='native source quartic'):
        apply_conditional_bc_direct_quartic_map('B_-1',replace(original,provenance=updated))


def test_controller_applies_map_once_and_only_in_mixed_branch():
    root=Path(__file__).resolve().parents[2]
    controller=root/'scripts'/'RunT3GeneralSequentialFitAblation.py'
    source=controller.read_text(encoding='utf-8')
    tree=ast.parse(source)
    functions={x.name:x for x in tree.body if isinstance(x,ast.FunctionDef)}
    assert '_mixed' in functions and '_nonmixed' in functions
    def count_calls(func):
        return sum(isinstance(n,ast.Call) and isinstance(n.func,ast.Name)
                   and n.func.id=='apply_conditional_bc_direct_quartic_map'
                   for n in ast.walk(func))
    assert count_calls(functions['_mixed'])==1
    assert count_calls(functions['_nonmixed'])==0
    assert source.index('tensor=scalar_tensors(model,tree,channels,snap,registry,geometry=geometry)',
                        source.index('def _mixed(')) < source.index(
        'tensor=apply_conditional_bc_direct_quartic_map(model,tensor)')
