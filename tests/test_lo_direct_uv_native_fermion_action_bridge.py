"""Synthetic regression tests. They do not constitute live native validation."""
from __future__ import annotations
import copy
import hashlib
import json
from pathlib import Path
import sys
import tempfile
import unittest
import numpy as np

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
from Numerical.lo_certification.direct_uv import native_fermion_action_bridge as br


def mock_native(label):
    maj=label!='A_p0'; masskind='epsilon' if not maj else 'singlet' if label=='B_m1' else 'delta'
    rgb={'format':'t3-ws4-simultaneous-native-rgbeta-yukawa-mass-v1',
      'status':'NativeRGBetaYukawaAndMassRegistryVerified','model':label,
      'dimensions_alpha':{'A_p0':[1,3,2,0],'B_m1':[2,2,1,-1],
                          'C_m1':[2,2,3,-1]}[label],
      'majorana':maj,'yukawas':[
        {'name':'y1','fields_inputform':'{Bar[S1], l, F}','chirality':'Right',
         'coupling_internal_factor':'Sqrt[2]','couplingbar_internal_factor':'Sqrt[2]',
         'physical_action_coefficient_certified':False},
        {'name':'y2','fields_inputform':'{S2, l, F}' if maj else '{S2, l, Fc}',
         'chirality':'Right','coupling_internal_factor':'Sqrt[2]',
         'couplingbar_internal_factor':'Sqrt[2]',
         'physical_action_coefficient_certified':False}],
      'mass':{'fields_inputform':'{F,F}' if maj else '{F,Fc}',
              'chirality':'Right','invariant_kind':masskind,
              'mass_internal_tensor_factor':'2' if not maj else '1',
              'factor_is_physical_Lagrangian_coefficient':False},
      'four_spinor_to_weyl_phase_measured':False,
      'native_su2_numeric_component_axes_measured':False,
      'physical_Matchete_RGBeta_dictionary_certified':False}
    mat={'format':'t3-ws4-restricted-native-matchete-fermion-action-v1',
       'status':'NativeRestrictedFermionActionCaptured','model':label,
       'native_lagrangian_validation':True,
       'checks':{'y1_present':True,'y2_present':True,
                 'both_yukawa_hc_present':True,'heavy_fermion_free_mass_symbol_present':True,
                 'mass_and_yukawa_use_newfermion':True},
       'mass_free_lagrangian_inputform':'NewFermion Coupling[MF,{}]',
       'yukawa_and_mixing_interactions_inputform':
          'Coupling[y1,{}] Coupling[y2,{}] NewFermion GammaCC',
       'physical_charge_conjugation_phase_measured':False,
       'physical_Matchete_RGBeta_dictionary_certified':False}
    return rgb,mat


def dump(p:Path,obj):
    p.parent.mkdir(parents=True,exist_ok=True)
    p.write_text(json.dumps(obj,indent=2),encoding='utf-8')


class NativeRegistryTests(unittest.TestCase):
    def test_three_representatives(self):
        for name in br.REPRESENTATIVES:
            self.assertFalse(br.validate_native(*mock_native(name),name)
                             ['physical_cross_package_normalisation_certified'])

    def test_wrong_epsilon_fails(self):
        a,b=mock_native('A_p0');a['mass']['invariant_kind']='delta'
        with self.assertRaises(br.FermionBridgeError):br.validate_native(a,b,'A_p0')

    def test_wrong_majorana_mass_fields_fails(self):
        a,b=mock_native('B_m1');a['mass']['fields_inputform']='{F,Fc}'
        with self.assertRaises(br.FermionBridgeError):br.validate_native(a,b,'B_m1')

    def test_wrong_dirac_y2_orientation_fails(self):
        a,b=mock_native('A_p0');a['yukawas'][1]['fields_inputform']='{S2,l,F}'
        with self.assertRaises(br.FermionBridgeError):br.validate_native(a,b,'A_p0')

    def test_wrong_yukawa_norm_fails(self):
        a,b=mock_native('C_m1');a['yukawas'][0]['coupling_internal_factor']='2'
        with self.assertRaises(br.FermionBridgeError):br.validate_native(a,b,'C_m1')

    def test_claimed_physical_phase_fails(self):
        a,b=mock_native('A_p0');a['four_spinor_to_weyl_phase_measured']=True
        with self.assertRaises(br.FermionBridgeError):br.validate_native(a,b,'A_p0')

    def test_declared_native_certificate_fails(self):
        a,b=mock_native('B_m1');b['physical_Matchete_RGBeta_dictionary_certified']=True
        with self.assertRaises(br.FermionBridgeError):br.validate_native(a,b,'B_m1')

    def test_unverified_action_fails(self):
        a,b=mock_native('C_m1');b['checks']['y1_present']=False
        with self.assertRaises(br.FermionBridgeError):br.validate_native(a,b,'C_m1')

    def test_mass_factor_non_numeric_fails(self):
        a,b=mock_native('A_p0');a['mass']['mass_internal_tensor_factor']='MF'
        with self.assertRaises(br.FermionBridgeError):br.validate_native(a,b,'A_p0')

    def test_malformed_native_format_fails(self):
        a,b=mock_native('A_p0');a['format']='old-format'
        with self.assertRaises(br.FermionBridgeError):br.validate_native(a,b,'A_p0')


class ExactAndTensorTests(unittest.TestCase):
    def test_exact_complex_values(self):
        self.assertAlmostEqual(br.exact('Sqrt(2)') .real, 2**0.5)
        self.assertEqual(br.exact('I'),1j)

    def test_rotation_identity(self):
        v={'dimensions':[2,1], 'components':[{'indices':[1,1],'cg':'2'},
                                              {'indices':[2,1],'cg':'-3'}]}
        actual=br.rotated_tensor(v,(False,False),{1:np.eye(1),2:np.eye(2)})
        np.testing.assert_allclose(actual,np.array([[2],[-3]]))

    def test_conjugate_rotation(self):
        u=np.array([[1,1j],[1j,1]])/np.sqrt(2)
        v={'dimensions':[2,1],'components':[{'indices':[1,1],'cg':'I'}]}
        actual=br.rotated_tensor(v,(True,False),{1:np.eye(1),2:u})
        np.testing.assert_allclose(actual[:,0],1j*np.conjugate(u[0]))

    def test_bad_cg_index(self):
        v={'dimensions':[2,1],'components':[{'indices':[3,1],'cg':'1'}]}
        with self.assertRaises(br.FermionBridgeError):
            br.rotated_tensor(v,(True,False),{1:np.eye(1),2:np.eye(2)})

    def test_bad_exact_value(self):
        with self.assertRaises(br.FermionBridgeError):br.exact('q')

    def test_benchmark_names(self):
        self.assertEqual(br.benchmark_name('E_-4'),
                         'T3_dS1_3_dS2_3_dF_2_alpha_m4')
        self.assertEqual(br.benchmark_name('B_+1'),
                         'T3_dS1_2_dS2_2_dF_1_alpha_p1')


class SyntheticFullBatch(unittest.TestCase):
    def setUp(self):
        self.tmp=tempfile.TemporaryDirectory()
        self.root=Path(self.tmp.name)
        self.native=self.root/'native'
        self.audit=self.root/br.COMPONENT_AUDITS
        self.audit.mkdir(parents=True)
        # Tiny contract-matching snippets; tests of the guards are separate
        source_lines={
          'rgbeta':'''AddYukawa[y1,{Bar@S1,l,F}, AddYukawa[y2,{S2,l,secondHeavy},
           secondHeavy=If[yF===0,F,Fc]; AddFermionMass[MF,{F,F},GroupInvariant->fermionInvariant,
           AddFermionMass[MF,{F,Fc},GroupInvariant->fermionInvariant,
           dF===2,(eps[T3RGBetaRep[dF],#1,#2]&)
           dF===3,(del[T3RGBetaRep[dF],#1,#2]&)''',
          'topology':'''useCConj=which===1;scalarBar=which===2;
             y[p,r]NCM[Bar[l[i,p]],f]scalarCG[cg,labels]
             lambdaT3[]H[i]H[j]s1s2CG[T3MixCG,labels]''',
          'model':'''SelfConjugate->True Mass->T3MassSpec["F",Lookup[f,"MassSymbol",MF]]''',
          'builder':'''metric=np.asarray(us[df].T*us[df],complex)
            y1dag[:,index,s]+=np.conjugate(y1[:,r])*factor
            y2plain[:,r*df+f,ds1+s]+=np.conjugate(y2[:,r])*factor'''}
        for name,rel in br.SRC.items():
            p=self.root/rel;p.parent.mkdir(parents=True,exist_ok=True)
            p.write_text(source_lines[name],encoding='utf-8')
        for label in br.REPRESENTATIVES:
            rgb,mat=mock_native(label)
            dump(self.native/(label+'_rgbeta.json'),rgb)
            dump(self.native/(label+'_matchete.json'),mat)
        us={1:np.eye(1),2:np.eye(2),3:np.eye(3)}
        dump(self.audit/'matchete_generator_analysis.json',{
            'representations':[{'dimension':d,'status':'generator_algebra_verified',
                'U_columns_are_weight_eigenvectors':[[str(int(k==j)) for j in range(d)] for k in range(d)]}
                 for d in (2,3)]})
        cases=[];cov=[]
        for label,(d1,d2,df,alphas) in br.CLASSES.items():
            for a in alphas:
                short=f'{label}_{a:+d}'; info={'dS1':d1,'dS2':d2,'dF':df,'alpha':a}
                vtx={}
                for name,dsc in (('T3Y1CG',d1),('T3Y2CG',d2)):
                    found=[]
                    for k in range(df):
                        for s in range(dsc):
                            yf=(a+1)/2; ys=a/2 if name=='T3Y1CG' else (a+2)/2
                            if abs(br.charge(df,yf,k)-br.charge(dsc,ys,s))<1e-10:
                                found.append({'indices':[1,k+1,s+1], 'cg':'1'})
                    vtx[name]={'dimensions':[2,df,dsc],'components':found}
                case={'model':info,'vertices':vtx}
                cases.append(case);cov.append({'model':info,'all_covariant':True})
                mf=np.diag([7.,9.,11.]);y1=np.arange(9).reshape(3,3)/5+1j/7
                y2=np.arange(9).reshape(3,3)/7-1j/9
                bench={'base_state':{'ordinary':{'MF':mf.tolist(),
                    'y1':{'real':y1.real.tolist(),'imag':y1.imag.tolist()},
                    'y2':{'real':y2.real.tolist(),'imag':y2.imag.tolist()}}}}
                benchpath=self.root/br.BENCHMARK/(br.benchmark_name(short)+'.json')
                dump(benchpath,bench)
                # Independently synthesise frozen conditional tensor entries.
                nc=3*df;majorana=bool(df%2 and a==-1)
                nf=nc if majorana else 2*nc;ns=d1+d2
                m=np.zeros((nf,nf),complex);yy=np.zeros((3,nf,ns),complex)
                zz=np.zeros((3,nf,ns),complex)
                for r in range(3):
                    for t in range(3):
                        if majorana:m[r*df:(r+1)*df,t*df:(t+1)*df]=mf[r,t]*np.eye(df)
                        else:
                            for k in range(df):
                                m[r*df+k,nc+t*df+k]=mf[r,t]
                                m[nc+t*df+k,r*df+k]=mf[r,t]
                for r in range(3):
                    for c in vtx['T3Y1CG']['components']:
                        k=c['indices'][1]-1;s=c['indices'][2]-1
                        zz[:,r*df+k if majorana else nc+r*df+k,s]+=y1[:,r].conj()
                    for c in vtx['T3Y2CG']['components']:
                        k=c['indices'][1]-1;s=c['indices'][2]-1
                        yy[:,r*df+k,d1+s]+=y2[:,r].conj()
                enc=lambda x: {'real':x.real.tolist(),'imag':x.imag.tolist()}
                f=self.root/br.FROZEN/(short+'.json')
                dump(f,{'model':short,'fermion_mass_gev':enc(m),
                    'vertices':[{'name':'y1_HC','conjugate':enc(zz)},
                                {'name':'y2_HC','field':enc(yy)}],
                    'provenance':{'cg_support_sha256':'PENDING',
                                  'generator_analysis_sha256':'PENDING',
                                  'benchmark_sha256':br.sha(benchpath)}})
        dump(self.audit/'component_support_audit.json',{'cases':cases,'errors':[]})
        dump(self.audit/'cg_covariance_audit.json',{
            'status':'cg_generator_covariance_verified','models':cov})
        self.frozen_hash=br.sha(self.audit/'component_support_audit.json')
        self.gen_hash=br.sha(self.audit/'matchete_generator_analysis.json')
        for p in (self.root/br.FROZEN).glob('*.json'):
            x=br.load(p);x['provenance']['cg_support_sha256']=self.frozen_hash
            x['provenance']['generator_analysis_sha256']=self.gen_hash
            dump(p,x)

    def tearDown(self):self.tmp.cleanup()

    def test_full_sixteen(self):
        out=br.evaluate(self.root,self.native)
        self.assertEqual(out['model_count'],16)
        self.assertEqual(out['native_representative_count'],3)
        self.assertEqual(out['physical_source_to_weyl_normalisation_certified_models'],0)
        self.assertLess(max(max(r['recomputed_candidate_vertex_and_mass_residuals'].values())
                            for r in out['models']),2e-10)

    def test_wrong_hash_blocks(self):
        p=self.root/br.FROZEN/'A_+0.json';data=br.load(p)
        data['provenance']['cg_support_sha256']='deadbeef';dump(p,data)
        with self.assertRaisesRegex(br.FermionBridgeError,'provenance hash'):
            br.evaluate(self.root,self.native)

    def test_wrong_frozen_yukawa_blocks(self):
        p=self.root/br.FROZEN/'B_-1.json';data=br.load(p)
        data['vertices'][1]['field']['real'][0][0][2]+=1
        dump(p,data)
        with self.assertRaisesRegex(br.FermionBridgeError,'differs'):
            br.evaluate(self.root,self.native)

    def test_source_drift_blocks(self):
        p=self.root/br.SRC['rgbeta'];p.write_text('changed',encoding='utf-8')
        with self.assertRaisesRegex(br.FermionBridgeError,'Source contract'):
            br.evaluate(self.root,self.native)

    def test_missing_source_audit_blocks(self):
        (self.audit/'cg_covariance_audit.json').unlink()
        with self.assertRaisesRegex(br.FermionBridgeError,'Missing canonical source audit'):
            br.evaluate(self.root,self.native)

    def test_missing_native_evidence_blocks(self):
        (self.native/'C_m1_rgbeta.json').unlink()
        with self.assertRaisesRegex(br.FermionBridgeError,'Missing evidence'):
            br.evaluate(self.root,self.native)


if __name__=='__main__':unittest.main()
