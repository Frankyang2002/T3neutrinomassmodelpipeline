"""Synthetic/structural regression tests for the B/C native mixed-F tree adapter.

These do not certify physical GammaCC normalization. Run the repository's
native CG tests as well (see README) before interpreting a numerical audit.
"""
from __future__ import annotations

import ast
from dataclasses import dataclass
from hashlib import sha256
import json
from pathlib import Path
import re
import sys
from types import ModuleType

import numpy as np
import pytest

from Numerical.running import GeneralT3MixedExactFermionWilson as mix


def term(*, y1='y1', y2='y2', denominator='2*Sqrt[3]',
         mixing='T3MEU12', bar=True):
    f=lambda y,n:(f'Bar[Coupling[{y}, {{Index[d$${n}, Flavor], '
                    f'Index[d$$1, NFlavor]}}, 0]]' if bar else
                    f'Coupling[{y}, {{Index[d$${n}, Flavor], '
                    f'Index[d$$1, NFlavor]}}, 0]')
    u = f'Bar[Coupling[{mixing}, {{}}, 0]]' if bar else f'Coupling[{mixing}, {{}}, 0]'
    return (f'({u}^2*{f(y1,1)}*{f(y2,2)}*'
            f'Field[NewScalarX, Scalar, {{Index[d$$1, T3BSMd2]}}, {{}}]*'
            f'NCM[Transp[Field[l, Fermion, {{Index[d$$2, SU2L[fund]], Index[d$$2, Flavor]}}, {{}}]], '
            f'DiracProduct[GammaCC, Proj[-1]], Field[l, Fermion, '
            f'{{Index[d$$1, SU2L[fund]], Index[d$$1, Flavor]}}, {{}}]])'
            f'/({denominator}*Coupling[MF, {{}}, 0]))')


def test_mixed_c_minus1_denominator_and_barred_yukawas():
    row=mix.normalise_mixed_term(term(y1='y1',y2='y2'))
    assert [(f['name'],f['conjugated']) for f in row['factors']]==[
        ('y1',True),('y2',True)]
    assert row['flavour_matrix_axes_external_indices']==['d$$1','d$$2']
    assert 'Sqrt[3]' in row['operator_template_unit_flavor_inputform']
    assert 'Coupling[y' not in row['operator_template_unit_flavor_inputform']
    assert 'Coupling[MF' not in row['operator_template_unit_flavor_inputform']
    assert 'Proj[-1]' in row['operator_template_unit_flavor_inputform']


def test_parser_rejects_incorrect_chirality_and_missing_mass():
    with pytest.raises(ValueError, match='PL'):
        mix.normalise_mixed_term(term().replace('Proj[-1]','Proj[1]'))
    with pytest.raises(ValueError, match='MF'):
        mix.normalise_mixed_term(term().replace('Coupling[MF, {}, 0]','100000'))


def test_current_complex_scalar_U_keeps_bar_and_squared_powers():
    U=np.array([[0.9,0.2+0.3j],[-0.2+0.3j,0.9]],complex)
    raw='Bar[Coupling[T3MEU12, {}, 0]]^2*Coupling[T3MEU21, {}, 0]'
    converted,seen=mix.insert_current_mixing(raw,U)
    assert seen==('T3MEU12','T3MEU21')
    assert 'Coupling[' not in converted
    assert '+(-0.29999999999999999)*I)^2' in converted
    assert '(-0.20000000000000001+(0.29999999999999999)*I)' in converted or '(-0.2+' in converted


def test_zero_mixing_is_finite_and_not_divided_by_archived_U():
    U=np.eye(2)
    expression,labels=mix.insert_current_mixing('Coupling[T3MEU12, {}, 0]^2',U)
    assert labels==('T3MEU12',)
    assert '0' in expression
    assert 'Coupling[' not in expression


def make_native_records(tmp_path:Path, model='B_-1'):
    count=mix.EXPECTED_PL[model]
    dirname=tmp_path/mix.SOURCE_ROOT/mix.source_tag(model)
    dirname.mkdir(parents=True)
    raw=term()
    start=11 if model=='B_-1' else 13
    pl=[dict(source_index=i+start,chirality='PL',source_term_inputform=raw,
             spurion_names=['T3MEU12'],factors=[
                  {'name':'y1','conjugated':True,'external_index':'d$$1','heavy_index':'d$$1'},
                  {'name':'y2','conjugated':True,'external_index':'d$$2','heavy_index':'d$$1'}],
             flavour_matrix_axes_external_indices=['d$$1','d$$2']) for i in range(count)]
    pr=[dict(source_index=i+1,chirality='PR') for i in range(count)]
    source={'format':'t3-mixed-factorised-LLSS-source-v1',
            'model':model,'mu_F_gev':1e5,'source_status_sha256':'abc',
            'source_mass_basis_sha256':'def','llss_sha256':'ghi','operator_terms':pr+pl}
    cg={'format':'t3-mixed-matchete-native-cg-registry-v1',
        'model':model,'CGRegistry':[{'Name':str(i)} for i in range(2 if model=='B_-1' else 3)]}
    cf=dirname/'mixed_native_tree_llss_cg_registry.json'
    cf.write_text(json.dumps(cg),encoding='utf8')
    ordered={'format':'t3-mixed-native-matchete-ordered-tensor-v1',
             'model':model,'mu_F_gev':1e5,'source_status_sha256':'abc',
             'source_mass_basis_sha256':'def','source_synced_matchete_llss_sha256':'ghi',
             'source_native_cg_sha256':sha256(cf.read_bytes()).hexdigest(),
             'operators':[dict(source_index=i+start) for i in range(count)]}
    (dirname/'source_flavour.json').write_text(json.dumps(source),encoding='utf8')
    (dirname/'ordered_tensor.json').write_text(json.dumps(ordered),encoding='utf8')
    return source,pl,cg,ordered,dirname


@pytest.mark.parametrize('model', ['B_-1','C_-1'])
def test_native_source_inventory_and_provenance(tmp_path,model):
    _=make_native_records(tmp_path,model)
    source,pl,cg,ordered=mix.read_native(model,tmp_path)
    assert len(pl)==mix.EXPECTED_PL[model]
    assert [r['source_index'] for r in pl]==[
        i+(11 if model=='B_-1' else 13) for i in range(mix.EXPECTED_PL[model])]
    assert len(cg['CGRegistry']) in (2,3)


def test_native_source_checksum_mismatch_is_rejected(tmp_path):
    source,pl,cg,ordered,dirpath=make_native_records(tmp_path)
    ordered['source_native_cg_sha256']='wrong'
    (dirpath/'ordered_tensor.json').write_text(json.dumps(ordered),encoding='utf8')
    with pytest.raises(ValueError,match='CG registry checksum'):
        mix.read_native('B_-1',tmp_path)


def test_component_reference_reproduction_detects_changes(monkeypatch):
    projected=ModuleType('Numerical.running.GeneralT3OrderedTensor')
    monkeypatch.setitem(sys.modules, 'Numerical.running.GeneralT3OrderedTensor',projected)
    fake={'source_index':11,'components':[
         {'lepton_su2':[1,1],'scalar_real':[5,5],
         'ordered_flavour_coefficient_gev_inv':{
         'real':np.eye(3).tolist(),'imag':np.zeros((3,3)).tolist()}}]}
    projected.project_ordered_term=lambda *a:fake
    pl=[{'source_index':11}]
    assert mix.archived_reference_reproduction('B_-1',pl,{'CGRegistry':[]},
                {'operators':[fake]})==0
    wrong=json.loads(json.dumps(fake))
    wrong['components'][0]['ordered_flavour_coefficient_gev_inv']['real'][0][0]=2
    with pytest.raises(ValueError,match='cannot reproduce'):
        mix.archived_reference_reproduction('B_-1',pl,{'CGRegistry':[]},
                                            {'operators':[wrong]})


def test_runner_does_not_change_nonmixed_backend():
    src=Path(__file__).resolve().parents[2]/'scripts/RunT3GeneralSequentialFitAblation.py'
    text=src.read_text(encoding='utf8')
    program=ast.parse(text)
    methods={f.name:ast.get_source_segment(text,f) for f in program.body
             if isinstance(f,(ast.FunctionDef,ast.AsyncFunctionDef))}
    nonmixed=methods['_nonmixed']
    mixed=methods['_mixed']
    assert 'nonmixed_boundary(model,uv_end' in nonmixed
    assert 'mixed_exact_boundary(' not in nonmixed
    assert 'mixed_exact_boundary(' in mixed
    assert 'mixed_reference_boundary(' not in mixed
    assert "source_details[f'{int(uv_on)}{int(eft_on)}']" in mixed


def test_source_exact_preserves_independent_majorana_yukawa_channels(monkeypatch):
    # Projector and canonicalizer are isolated so this tests the dynamic native
    # Yukawa reconstruction, not synthetic Clebsch-Gordan physics.
    projected=ModuleType('Numerical.running.GeneralT3OrderedTensor')
    dense=ModuleType('Numerical.running.GeneralT3DynamicFermionWilson')
    monkeypatch.setitem(sys.modules,'Numerical.running.GeneralT3OrderedTensor',projected)
    monkeypatch.setitem(sys.modules,'Numerical.running.GeneralT3DynamicFermionWilson',dense)
    def projector(row,registry,model):
        return {'source_index':row['source_index'],'components':[
            {'lepton_su2':[1,1],'scalar_real':[5,5],
             'ordered_flavour_coefficient_gev_inv':row['ordered_flavour_matrix_gev_inv']}]}
    def canonicalizer(operators,rows,model):
        mats=[np.asarray(o['components'][0]['ordered_flavour_coefficient_gev_inv']['real'])+
              1j*np.asarray(o['components'][0]['ordered_flavour_coefficient_gev_inv']['imag'])
              for o in operators]
        out=np.zeros((len(mats),3,3,2,2,8,8),complex)
        for i,w in enumerate(mats):out[i,:,:,0,0,0,0]=w
        return out,tuple(int(x['source_index']) for x in rows),{'fake_projector':True}
    projected.project_ordered_term=projector
    dense._dense_operator=canonicalizer
    ys={
       'y1':np.array([[1,0,0],[0,2,0],[1,2,3]],complex),
       'y2':np.array([[0,3,0],[4,0,5],[0,2,1]],complex),
    }
    U=np.eye(2,dtype=complex)
    rows=[]
    for idx,pa,pb in [(13,'y1','y1'),(14,'y1','y2'),
                      (15,'y2','y2'),(23,'y1','y2')]:
        raw=term(y1=pa,y2=pb,mixing='T3MEU11')
        normalized=mix.normalise_mixed_term(raw)
        rows.append({'source_index':idx,'source_term_inputform':raw,
                     'factors':normalized['factors'],
                     'spurion_names':['T3MEU11'],
                     'flavour_matrix_axes_external_indices':
                         normalized['flavour_matrix_axes_external_indices']})
    by,indices,details,audit=mix.project_source('C_-1',{},rows,
                                {'CGRegistry':[]},np.array([1,2,4]),ys,U)
    assert indices==(13,14,15,23)
    assert by.shape==(4,3,3,2,2,8,8)
    assert np.linalg.norm(by[0]-by[1])>1e-4
    assert np.linalg.norm(by[0]-by[2])>1e-4
    assert np.linalg.norm(by[1]-by[2])>1e-4
    assert np.linalg.norm(by[3])>1e-4
    assert details[0]['yukawa_pair']==['y1','y1']
    assert details[1]['yukawa_pair']==['y1','y2']
    assert details[2]['yukawa_pair']==['y2','y2']
