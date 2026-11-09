"""v18-v20 native complex scalar source checks, rank atlas, resume and negatives."""
from __future__ import annotations
import copy
from hashlib import sha256
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch
import numpy as np
from Numerical.parallel_t3.wilson.complex_scalar_bridge_v19 import (
    native_document, historical_extra_map, topology_action, SOURCES, MAT_FORMAT,
    MAT_STATUS, STATUS, FORMAT, SPECIAL, BC, V19Blocked, require,
    ALL16, DIMENSIONS,
)
from Numerical.parallel_t3.wilson.scalar_topology_action_v12 import topology_tensor
from Numerical.parallel_t3.wilson.full_scalar_census_v20 import (
    audit as census20, FORMAT as V20_FORMAT, STATUS as V20_STATUS,
)
from Numerical.parallel_t3.wilson.full_scalar_census_v17 import (
    FORMAT as V17_FORMAT, STATUS as V17_STATUS,
)


def as_input(z):
    z=complex(z)
    if abs(z.imag)<1e-14:return repr(float(z.real))
    if abs(z.real)<1e-14:return repr(float(z.imag))+'*I'
    return '('+repr(float(z.real))+'+'+repr(float(z.imag))+'*I)'


def fake_sources(root):
    result={}
    for p in SOURCES:
        s=root/p;s.parent.mkdir(parents=True,exist_ok=True)
        s.write_text('(* synthetic native source for unit tests '+p+' *)\n',encoding='utf-8')
        result[p]=sha256(s.read_bytes()).hexdigest()
    return result


def cg_rec(symbol,t,fields,cg,**extras):
    return dict(symbol=symbol,dimensions=list(t.shape),fields=fields,cg=cg,
                components_inputform=[as_input(x) for x in t.ravel()],
                native_dimensions=[d for d in t.shape if d>1],
                assembled_lagrangian_contains_operator=True,**extras)


def fake_native(root,mode):
    hashes=fake_sources(root)
    models={}
    rng=np.random.default_rng(1921)
    names=[('110101',1,['H','Hdag','S1dag','S2dag']),
           ('110101',2,['H','Hdag','S1dag','S2dag']),
           ('000112',1,['S1dag','S2','S2dag','S2dag']),
           ('000202',1,['S1dag','S1dag','S2dag','S2dag']),
           ('001201',1,['S1','S1dag','S1dag','S2dag']),
           ('020020',1,['Hdag','Hdag','S2','S2']),
           ('020200',1,['Hdag','Hdag','S1dag','S1dag'])]
    for tag in sorted(ALL16 if mode=='TOPOLOGY' else BC):
        d1,d2,df=DIMENSIONS[tag[0]]
        if mode=='TOPOLOGY':
            t=topology_tensor(d1,d2)
            rows=[cg_rec('lambdaT3',t,['H','H','S1','Bar[S2]'],'T3MixCG')]
        else:
            rows=[]
            for key,k,fields in names:
                a=rng.normal(size=(2,2,2,2))+.2j*rng.normal(size=(2,2,2,2))
                symbol='lambdaT3ExtraQ'+key+'I'+str(k)
                rows.append(cg_rec(symbol,a,fields,'T3ExtraCGQ'+key+'I'+str(k),
                                   counts=[int(x) for x in key],operator='T3ExtraQ_'+
                                   '_'.join(key)+'_'+str(k)))
        models[tag]={'dimensions':[d1,d2,df],'alpha':int(tag.split('_')[1]),
                     'native_builder_status':'Success',
                     'matching_mass_basis_physically_certified':False,
                     'symbolic_quartic_only':True,'records':rows}
    doc={'format':MAT_FORMAT,'status':MAT_STATUS,'mode':mode,'models':models,
         'source_sha256':hashes,'models_checked':len(models),
         'native_complex_CG_slots':16 if mode=='TOPOLOGY' else 14,
         'independent_RGBeta_native_component_frame_certified':False,
         'BC_minus1_mixed_scalar_mass_diagonalized':False,
         'physical_Wilson_evolution_authorized':False}
    return doc


class TestV18Source(unittest.TestCase):
    def test_all16_topology_native_tensor_shapes_and_actions(self):
        with tempfile.TemporaryDirectory() as td:
            root=Path(td);obs=native_document(fake_native(root,'TOPOLOGY'),root,'TOPOLOGY')
            self.assertEqual(len(obs),16)
            self.assertLess(max(topology_action(v['lambdaT3'][1],samples=3)
                                for v in obs.values()),3e-14)

    def test_two_bc_seven_extra_records(self):
        with tempfile.TemporaryDirectory() as td:
            root=Path(td);obs=native_document(fake_native(root,'BCEXTRA'),root,'BCEXTRA')
            self.assertEqual(sum(len(row) for row in obs.values()),14)

    def test_source_hash_or_missing_source_fails(self):
        with tempfile.TemporaryDirectory() as td:
            root=Path(td);doc=fake_native(root,'TOPOLOGY')
            (root/SOURCES[2]).write_text('modified native source',encoding='utf-8')
            with self.assertRaises(V19Blocked):native_document(doc,root,'TOPOLOGY')

    def test_broken_topology_higgs_exchange_rejected(self):
        with tempfile.TemporaryDirectory() as td:
            root=Path(td);doc=fake_native(root,'TOPOLOGY')
            row=doc['models']['B_-1']['records'][0]
            a=topology_tensor(2,2);a[0,1,0,0]+=1.
            row['components_inputform']=[as_input(x) for x in a.ravel()]
            with self.assertRaises(V19Blocked):native_document(doc,root,'TOPOLOGY')

    def test_false_mass_basis_physical_promotion_rejected(self):
        with tempfile.TemporaryDirectory() as td:
            root=Path(td);doc=fake_native(root,'BCEXTRA')
            doc['BC_minus1_mixed_scalar_mass_diagonalized']=True
            with self.assertRaises(V19Blocked):native_document(doc,root,'BCEXTRA')

    def test_wrong_native_cg_record_or_duplicates_rejected(self):
        with tempfile.TemporaryDirectory() as td:
            root=Path(td);doc=fake_native(root,'BCEXTRA')
            doc['models']['B_-1']['records'][1]['symbol']=doc['models']['B_-1']['records'][0]['symbol']
            with self.assertRaises(V19Blocked):native_document(doc,root,'BCEXTRA')


class TestV19SixOfSeven(unittest.TestCase):
    def data(self):
        with tempfile.TemporaryDirectory() as td:
            root=Path(td);native=native_document(fake_native(root,'BCEXTRA'),root,'BCEXTRA')
        old={'models':{}};fresh={'models':{}};saved={'models':{}}
        rgb=['lambdaT3GeneralQ110101I1']+list(SPECIAL)
        a=np.eye(7)
        for tag in BC:
            old['models'][tag]={'tensors':[]}
            current=native[tag]
            for symbol,(rec,t) in current.items():
                old['models'][tag]['tensors'].append({
                    'symbol':symbol,'species':rec['fields'],
                    'dimensions':rec['dimensions'],'components':[as_input(x) for x in t.ravel()]})
            names=list(current)
            row={'matchete_symbols':names,'rgbeta_symbols':rgb,
                 'map_rgbeta_from_matchete':{'real':a.tolist(),'imag':np.zeros((7,7)).tolist()}}
            fresh['models'][tag]=copy.deepcopy(row);saved['models'][tag]=copy.deepcopy(row)
        return native,old,fresh,saved

    def test_BC_6_rank_plus_missing_one_7d_direction(self):
        native,old,fresh,saved=self.data()
        r=historical_extra_map(fresh,native,old,{},saved)
        self.assertEqual(set(r),set(BC))
        for tag in BC:
            self.assertEqual(r[tag]['baseline_map_rank'],6)
            self.assertEqual(r[tag]['baseline_kernel_dimension_in_Matchete_seven_space'],1)
            self.assertEqual(r[tag]['optional_extra_rgbeta_symbol'],'lambdaT3GeneralQ110101I1')

    def test_changed_native_current_extra_CG_rejected(self):
        native,old,fresh,saved=self.data()
        key=next(iter(native['B_-1']))
        rec,t=native['B_-1'][key];t[0,0,0,0]+=1.
        with self.assertRaises(V19Blocked):historical_extra_map(fresh,native,old,{},saved)

    def test_missing_baseline_special_symbol_rejected(self):
        native,old,fresh,saved=self.data()
        fresh['models']['C_-1']['rgbeta_symbols'][1]='not_a_baseline_special'
        with self.assertRaises(V19Blocked):historical_extra_map(fresh,native,old,{},saved)

    def test_rank_deficient_BC_baseline_rejected(self):
        native,old,fresh,saved=self.data()
        fresh['models']['B_-1']['map_rgbeta_from_matchete']['real'][1]=[0.]*7
        with self.assertRaises(V19Blocked):historical_extra_map(fresh,native,old,{},saved)

    def test_historical_name_changes_rejected(self):
        native,old,fresh,saved=self.data()
        old['models']['B_-1']['tensors'][0]['symbol']='fabricated_operator'
        with self.assertRaises(V19Blocked):historical_extra_map(fresh,native,old,{},saved)


class TestV20Census(unittest.TestCase):
    def input(self,root):
        names=['v17','v19','native14','topo','bc','v15','m50','r50','v3']
        paths=[root/n for n in names]
        for p in paths:p.write_text('{}',encoding='utf-8')
        models={};real_models={}
        for tag in ALL16:
            d1,d2,_=DIMENSIONS[tag[0]]
            portal=(1 if d1==1 else 2)+(1 if d2==1 else 2)
            nonportal=(2 if d1==3 else 1)+(2 if d2==3 else 1)+(
                3 if d1==d2==3 else 2 if d1>1 and d2>1 else 1)
            real_models[tag]={'total_real_quartic_slots':1+portal+nonportal}
            models[tag]={'baseline_complex_quartic_slots':7 if tag in BC else 1}
        self.assertEqual(sum(x['total_real_quartic_slots'] for x in real_models.values()),150)
        self.assertEqual(sum(x['baseline_complex_quartic_slots'] for x in models.values()),28)
        v19={'format':FORMAT,'status':STATUS,'input_sha256':{'x':'y'},
            'certification':{'physical_Wilson_evolution_authorized':False},
            'native_Matchete_topology_CGs':16,'native_Matchete_BC_extra_CGs_rechecked':14,
            'baseline_EFT1_complex_slots':28,
            'max_native_topology_6_over_24_real_action_relative_residual':1e-14,
            'models':models}
        v17={'format':V17_FORMAT,'status':V17_STATUS,'models':real_models,
             'native_RGBeta_real_quartic_slots':150,
             'native_RGBeta_real_projector_unit_diagonals':150,
             'native_RGBeta_portal_response_matrix_entries':110,
             'Matchete_nonportal_scalar_slots':76,'Matchete_Higgs_portal_slots':58,
             'Matchete_Higgs_quartic_slots':16,
             'certification':{'physical_Wilson_evolution_authorized':False}}
        from Numerical.parallel_t3.wilson.full_scalar_census_v20 import sha
        base=root/'output/parallel_t3/wilson';base.mkdir(parents=True)
        names={
            'Matchete_nonportal_v16':'matchete_nonportal_native_v16.json',
            'Matchete_portals_v5':'matchete_portals_native_v5.json',
            'native_RGBeta_v14':'rgbeta_quartic_projectors_native_v14.json',
            'nonportal_dictionary_v16':'nonportal_dictionary_conditional_v16.json',
            'portal_real_action_v9':'portal_real_action_conditional_v9.json',
            'native_quartic_v14_summary':'native_quartic_projectors_conditional_v14.json',
        }
        for name in names.values():(base/name).write_text('{}',encoding='utf-8')
        v17['input_sha256']={label:sha(base/name) for label,name in names.items()}
        paths[0].write_text(json.dumps(v17));paths[1].write_text(json.dumps(v19))
        return paths,v19,v17

    def invoke(self,modifier=None):
        with tempfile.TemporaryDirectory() as td:
            root=Path(td);paths,v19,v17=self.input(root)
            if modifier:modifier(v19,v17,paths)
            paths[0].write_text(json.dumps(v17));paths[1].write_text(json.dumps(v19))
            native={'source_census':{'quartic_registered':178,'real_registered':150,
                   'complex_registered':28,'real_diagonal_unit':150},
                   'complex_probe_census':{'complex_direct_unit_conj_zero':28}}
            with patch('Numerical.parallel_t3.wilson.full_scalar_census_v20.audit_v19',return_value=v19),\
                 patch('Numerical.parallel_t3.wilson.full_scalar_census_v20.audit_witness',return_value=native),\
                 patch('Numerical.parallel_t3.wilson.full_scalar_census_v20.read_v14',return_value={}):
                 return census20(root,*paths)

    def test_178_150_28_completed_conditionally(self):
        result=self.invoke()
        self.assertEqual(result['status'],V20_STATUS)
        self.assertEqual(result['total_EFT1_scalar_quartic_registration_slots'],178)
        self.assertFalse(result['certification']['physical_Wilson_evolution_authorized'])

    def test_saved_BC_complex_slot_count_tamper_fails(self):
        def change(v19,v17,paths):v19['models']['B_-1']['baseline_complex_quartic_slots']=6
        with self.assertRaises(V19Blocked):self.invoke(change)

    def test_saved_real_census_wrong_fails(self):
        def change(v19,v17,paths):v17['native_RGBeta_real_projector_unit_diagonals']=149
        with self.assertRaises(V19Blocked):self.invoke(change)

    def test_saved_false_physics_authorization_fails(self):
        def change(v19,v17,paths):v17['certification']['physical_Wilson_evolution_authorized']=True
        with self.assertRaises(V19Blocked):self.invoke(change)


class TestResumeAndSafety(unittest.TestCase):
    def test_cli_v19_resume_immutable_and_reject_changes(self):
        from scripts.RunT3ParallelAComplexBridgeV19 import resume
        with tempfile.TemporaryDirectory() as td:
            out=Path(td)/'v19.json'
            fresh={'format':FORMAT,'status':STATUS,'input_sha256':{'a':'b'},
                   'certification':{'physical_Wilson_evolution_authorized':False},
                   'baseline_EFT1_complex_slots':28,'native_Matchete_BC_extra_CGs_rechecked':14,
                   'max_native_topology_6_over_24_real_action_relative_residual':1e-15,
                   'models':{t:{'baseline_complex_quartic_slots':7 if t in BC else 1} for t in ALL16},
                   'BC_baseline_six_of_seven_conditional_source_maps':{
                       t:{'coefficients':[1,2]} for t in BC}}
            out.write_text(json.dumps(fresh),encoding='utf-8');before=out.read_bytes()
            resume(out,fresh);self.assertEqual(out.read_bytes(),before)
            bad=copy.deepcopy(fresh)
            bad['BC_baseline_six_of_seven_conditional_source_maps']['B_-1']['coefficients'][0]+=1
            with self.assertRaises(V19Blocked):resume(out,bad)

    def test_runner_restores_isolated_quartic_flags(self):
        source=(Path(__file__).parents[1]/'scripts'/'RunT3ParallelAComplexScalarBatchV18V20.ps1').read_text()
        self.assertIn('finally{',source)
        self.assertIn("SetEnvironmentVariable($f,$previous[$f],'Process')",source)
        self.assertIn('BCEXTRA',source)
        self.assertIn('--resume',source)
        self.assertNotIn('Match[]',source)

    def test_exporter_two_modes_and_fail_closed_physical(self):
        source=(Path(__file__).parents[1]/'scripts'/'RunT3ParallelAExportComplexMatcheteV18.wl').read_text()
        for k in ('TOPOLOGY','BCEXTRA','GetCGTensor','AdditionalQuarticsIncluded',
                  'physical_Wilson_evolution_authorized','source_sha256'):
            self.assertIn(k,source)
        self.assertNotIn('\nMatch[',source)

if __name__=='__main__':unittest.main()
