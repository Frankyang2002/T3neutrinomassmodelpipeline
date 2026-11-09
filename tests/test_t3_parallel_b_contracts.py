"""Basic real-provenance, shape, EFT-order, and no-invented-zero B smoke tests.
Run from repository root with: python -m pytest -q tests/test_t3_parallel_b_contracts.py
"""
from __future__ import annotations

import hashlib
import json
from pathlib import Path

import numpy as np
import pytest

from Numerical.parallel_t3.background.contracts import (
    MODELS, MIXED, SELF, BASE, HEAVY, ContractError, check_pair,
    check_beta, hash_file, save_json, tag, verify_batch1_provenance)
from Numerical.parallel_t3.background.runner import (
    _matrix_pattern, _endpoint_fingerprint, Options)
from Numerical.parallel_t3.background.api import (
    Background, load_formal_boundary)


def pair(model='A_-4'):
    fam, alpha = model.split('_')
    dims = {'A':(1,3,2),'B':(2,2,1),'C':(2,2,3),
            'D':(3,1,2),'E':(3,3,2)}[fam]
    base = set(BASE) | {'lambdaH2Adj'}
    return tuple({
        'format':'general-t3-quadratic-rgbeta-experimental-v1',
        'status':'SuccessCandidateMassBeta',
        'metadata':{'dS1':dims[0],'dS2':dims[1],'dF':dims[2],
                    'alpha':int(alpha),'GeneralScalarModel':model,
                    'GeneralScalarStage':stage, 'SharedScalar':False,
                    'ReportBetaConvention':
                      '16*pi^2*dX/dln(mu); gauge BetaTerm divided by 2*g',
                    'GeneralScalarBasisStatus':'CandidateNativeRGBetaUnmappedToMatchete',
                    'GeneralScalarNewBetaKeys':[]},
        'report_betas':{n:'1' for n in sorted(base | (set(('y1','y2','MF')) if stage=='UV' else set()))},
        'beta_symbols': sorted(base | (set(('y1','y2','MF')) if stage=='UV' else set())),
    } for stage in ('UV','EFT1'))


def test_all_16_and_special_models():
    assert len(MODELS)==16 and len(set(MODELS))==16
    assert MIXED=={'B_-1','C_-1'} and len(SELF)==6
    assert tag('A_-4')=='A_m4' and tag('D_+2')=='D_p2'


def test_beta_inventory_and_all_heavy_removed():
    uv,eft=pair()
    a,b=check_pair('A_-4',uv,eft)
    assert set(a)-set(b)=={'MF','y1','y2'}
    eft['report_betas']['MF']='1'
    eft['beta_symbols'].append('MF')
    with pytest.raises(ContractError, match='fermion'):
        check_pair('A_-4',uv,eft)


def test_reject_wrong_stage_and_empty_beta():
    uv,eft=pair('B_-1')
    uv['metadata']['GeneralScalarModel']='A_-4'
    with pytest.raises(ContractError,match='representation'):
        check_pair('B_-1',uv,eft)
    uv,eft=pair()
    del uv['report_betas']['lambdaH']
    with pytest.raises(ContractError):
        check_pair('A_-4',uv,eft)


def test_source_and_beta_sha_fail_closed(tmp_path):
    root=tmp_path
    p=root/'Lagrangian'/'interactions'/'ScalarPotential.wl'
    p.parent.mkdir(parents=True)
    p.write_text('T3HiggsPortalHermitianPhase\n')
    batch=root/'output'/'b1'
    (batch/'models').mkdir(parents=True)
    b= batch/'models/A_m4.json'
    b.write_text('{}')
    u=batch/'uv.json';e=batch/'eft.json'
    u.write_text('{}');e.write_text('{}')
    record={'model':'A_-4',
            'local_source_sha256':{'Lagrangian/interactions/ScalarPotential.wl':hash_file(p)},
            'rgbeta':{'uv':{'status':'SuccessCandidateMassBeta','sha256':hash_file(u)},
                      'eft1':{'status':'SuccessCandidateMassBeta','sha256':hash_file(e)}},
            'f_threshold':{'a4_portal':{'inv1_survives':True,'inv2_survives':True}}}
    verify_batch1_provenance(root,batch,'A_-4',record,u,e)
    p.write_text('old cancelled portal')
    with pytest.raises(ContractError,match='changed'):
        verify_batch1_provenance(root,batch,'A_-4',record,u,e)


def test_complex_flavour_and_endpoint_hash():
    a=_matrix_pattern('y1','A_-4')
    b=_matrix_pattern('y2','A_-4')
    assert a.shape==(3,3) and np.max(np.abs(a.imag))>0
    assert np.max(np.abs(a-np.diag(np.diag(a))))>0
    assert not np.allclose(a,b)
    h=_endpoint_fingerprint({'y1':np.stack((a,b))},['y1'])
    assert len(h)==64
    assert h!=_endpoint_fingerprint({'y1':np.stack((a,a))},['y1'])


def test_scales_and_finite_guard():
    with pytest.raises(ContractError):
        Options(root=Path('.'),mu_uv=1e4,mu_f=1e5).validated()
    Options(root=Path('.'),mu_uv=1e6,mu_f=1e5,mu_s=1e3).validated()


def test_interface_interpolation_and_no_f(tmp_path):
    out=tmp_path
    uv=np.array([1e6,1e5],float)
    # Reader requires >=3 sampled scales.
    uv=np.geomspace(1e6,1e5,3)
    eft=np.geomspace(1e5,1e3,3)
    arrays={'mu_uv_gev':uv, 'mu_eft1_gev':eft,
            't_uv':np.log(uv),'t_eft1':np.log(eft)}
    names=['gY','lambdaS1','mS1Sq']
    for n in names:
        arrays['uv__'+n]=np.full((3,), .1+0.02j, dtype=np.complex128)
        arrays['eft1__'+n]=np.full((3,), .1+0.02j,dtype=np.complex128)
    arrays['uv__MF']=np.tile(np.eye(3,dtype=np.complex128),(3,1,1))*1e5
    np.savez_compressed(out/'background.npz',**arrays)
    save_json(out/'background_manifest.json',{
        'schema':'t3-parallel-B-background-v1',
        'beta_symbol_order_uv':names+['MF'], 'beta_symbol_order_eft1':names,
        'background_npz_sha256':hash_file(out/'background.npz')})
    b=Background.load(out)
    mid=b.native_eft1_at(1e4)
    assert set(mid)==set(names) and np.isclose(mid['lambdaS1'], .1+.02j)
    assert 'MF' not in mid
    assert np.isclose(b.native_scalar_ms()['mS1Sq'],.1+.02j)
    with pytest.raises(ContractError):
        b.native_eft1_at(1e2)
    save_json(out/'fermion_boundary_manifest.json',{
        'schema':'t3-parallel-B-fermion-boundary-v1',
        'formal_LLsS_tree_available':False,'hard_C5_status':'hard_C5_unavailable'})
    with pytest.raises(ContractError, match='no evaluated'):
        load_formal_boundary(out)
