"""v25-v27 integration fixtures are SYNTHETIC, never native Matchete evidence."""
from __future__ import annotations
from copy import deepcopy
from pathlib import Path
import json
import sys
import tempfile
import types
import unittest
from unittest import mock
import numpy as np

from Numerical.parallel_t3.wilson import bc_native_mass_eigen_reconcile_v25 as v25
from Numerical.parallel_t3.wilson import bc_one_scale_probe_audit_v26 as v26
from Numerical.parallel_t3.wilson import bc_threshold_order_ledger_v27 as v27
from scripts import RunT3ParallelABCMatcheteThresholdBridgeV25V27 as cli

ROOT_NAMES=('B_-1','C_-1')
BASE=Path('output/parallel_t3/wilson')
SYMBOLS=('lambdaT3ExtraQ000112I1','lambdaT3ExtraQ000202I1',
         'lambdaT3ExtraQ001201I1','lambdaT3ExtraQ020020I1',
         'lambdaT3ExtraQ020200I1','lambdaT3ExtraQ110101I1',
         'lambdaT3ExtraQ110101I2')

def write(path,obj):
    path.parent.mkdir(parents=True,exist_ok=True)
    path.write_text(json.dumps(obj,indent=2,allow_nan=False)+'\n',encoding='utf8')
    return path

def enc(a):
    a=np.asarray(a,complex)
    return {'real':a.real.tolist(),'imag':a.imag.tolist()}

def cc(z):return {'real':float(z.real),'imag':float(z.imag)}

def fixtures(root):
    v21={'format':'FAKE-V21','status':'FAKE-V21-STATUS','models_checked':2,'models':{}}
    v22={'format':'FAKE-V22','status':'FAKE-V22-STATUS','full_scalar_quartic_slots':178,
         'certification':{'physical_Wilson_evolution_authorized':False}}
    for idx,tag in enumerate(ROOT_NAMES):
        m1,m2=1.9e6+idx*1e4,2.55e6+idx*2e4
        mu=5200.+8100.j
        A=np.array([[m1,mu],[mu.conjugate(),m2]],complex)
        vals,U=np.linalg.eigh(A)
        for k in range(2):
            pivot=np.argmax(abs(U[:,k]));U[:,k]*=np.conj(U[pivot,k]/abs(U[pivot,k]))
        masses=np.sqrt(vals)
        stem=tag[0]+'_m1_final'
        path_mass=root/'output/general_t3_mass_basis'/(stem+'_mass_basis.json')
        mass={'format':v25.MASS_FORMAT,'model':tag,'endpoint':'final',
              'physical_mass_matching_performed':False,
              'MatcheteRGBetaPhysicalCouplingMapVerified':False,
              'basis':{'mass_squared_flavour_matrix_gev2':enc(A),
                       'unitary_flavour_rotation_U':enc(U),
                       'masses_squared_gev2':vals.tolist()}}
        write(path_mass,mass)
        param={'format':v25.PARAM_FORMAT,'status':'NumericMassEigenstateBoundaryPreparedMatcheteNotRun',
               'model':tag,'endpoint':'final','source_sha256':v25.digest(path_mass),
               'source_path':'output\\general_t3_mass_basis\\'+stem+'_mass_basis.json',
               'dimensions_alpha':[2,2,1 if tag.startswith('B') else 3,-1],
               'quartic_source_count':17,'yukawa_source_count':2,
               'matchete_mass_basis_verified':False,'physical_matching_performed':False,
               'matchete_rgbeta_component_map_certified':False,'coefficient_map_certified':False,
               'original_mass_squared_A_gev2':enc(A),'unitary_U':enc(U),
               'heavy_mass_squared_gev2':vals.tolist(),
               'epsilon':[[0,1],[-1,0]],
               'source_eft1_scale_gev':1000,
               'threshold_plan':{'decreasing_scale_order':['NewScalarY','NewScalarX'],
                                 'sequential_matching_performed':False,
                                 'existing_eft1_endpoint_below_both_mass_thresholds':True}}
        write(root/'output/general_t3_mass_eigen_matchete'/(stem+'_parameters.json'),param)
        mat={'format':v25.NATIVE_FORMAT,'status':v25.NATIVE_STATUS,'model':tag,
             'mass_basis_ordinary_check_passed':True,'heavy_mass_basis_detailed_check_passed':True,
             'all_old_scalar_fields_removed':True,'all_17_quartics_present':True,
             'both_yukawa_interactions_present':True,'matchete_rotation_uses_exact_symbolic_spurions':True,
             'physical_matching_performed':False,'matching_to_two_scalar_thresholds_performed':False,
             'component_basis_map_certified':False,'matchete_rgbeta_quartic_normalization_certified':False,
             'diagonal_heavy_fields':['NewScalarX','NewScalarY'],
             'source_baseline_interactions':11,'source_extra_quartics':7,'source_new_interaction_count':18,
             'source_interactions':[{'name':f'Baseline{i}','class':'Baseline'} for i in range(11)]+
                                   [{'name':'T3ExtraQ'+str(i),'class':'AdditionalQuartic'} for i in range(7)],
             'rotation_spurion_values_from_numerical_U':{'T3MEU'+str(a+1)+str(b+1):cc(U[a,b])
                                                          for a in range(2) for b in range(2)},
             'native_additional_quartic_tensors':{'Certificates':{
                 'T3ExtraQ'+str(i):{'SourceCoupling':sym,'SourceTensorIndex':1,
                                   'ExactTensorProjectionVerified':True,
                                   'NativeInvariantSpanRank':2,'CombinedBasisRank':2}
                  for i,sym in enumerate(SYMBOLS)},
                  'SourceOperatorCount':7,'RecoupledOperatorCount':7,
                  'AllProjectedTensorsInvariant':True,'PhysicalCouplingMapCertified':False},
             'native_yukawa1_tensor':{'OldToNewCGFactor':'1','ExactTensorProjectionVerified':True,
                                    'FermionSU2Dimension':1 if tag.startswith('B') else 3,
                                    'PhysicalCouplingMapCertified':False},
             'native_yukawa2_epsilon_tensor':{'OldTimesEpsilonCGFactor':'-1' if tag.startswith('B') else '1',
                                    'ExactTensorProjectionVerified':True,
                                    'FermionSU2Dimension':1 if tag.startswith('B') else 3,
                                    'PhysicalCouplingMapCertified':False}}
        for key in ('epsilon_bridge_BB','epsilon_bridge_UU'):
            mat[key]={'native_intertwiner_inputform':'{{0, 1}, {-1, 0}}',
                      'phase_relative_to_matchete_component_epsilon_verified':True,
                      'physical_RGBeta_component_map_certified':False}
        for key,count in [('native_H2_portal_epsilon_tensors',2),('native_cross_scalar_epsilon_tensors',2),
                          ('native_S2_self_epsilon_tensor',1)]:
            mat[key]={'OldInvariantCount':count,'NewInvariantCount':count,
                      'TransportedRank':count,'CombinedBasisRank':count,
                      'ExactTensorProjectionVerified':True,'PhysicalCouplingMapCertified':False}
        path_native=write(root/'output/general_t3_mass_eigen_matchete'/(stem+'_matchete.json'),mat)
        v21['models'][tag]={'native_CG_sectors':[
            {'source_matchete_coupling':'lambdaT3'}]+
            [{'source_matchete_coupling':x} for x in SYMBOLS]}
        p=deepcopy(mat)
        p.update({'status':v26.PROBE_STATUS,'one_loop_match_stage_completed':True,
                  'matching_approach':'SimultaneousAllHeavyOneScaleExploratory',
                  'matching_EFT_order':5,'matching_loop_order':1,
                  'weinberg_extraction_status':'Success','weinberg_term_count':236,
                  'weinberg_holomorphic_term_count':118,'weinberg_conjugate_term_count':118,
                  'weinberg_operator_convention':v26.CONVENTION,
                  'ordered_coefficient_is_canonical_flavour_symmetric_c5':False,
                  'production_sequential_matching_performed':False,
                  'numerical_c5_computed':False})
        for key,suffix in [('matched_eft_inputform_file','_one_loop_probe_matched_eft.inputform.txt'),
                           ('weinberg_ordered_coefficient_inputform_file','_one_loop_probe_weinberg_ordered_coefficient.inputform.txt')]:
            name=stem+suffix
            p[key]=name
            dest=root/'output/general_t3_mass_eigen_matchete/one_loop_probe'/name
            dest.parent.mkdir(parents=True,exist_ok=True)
            dest.write_text('Symbolic Mathematica InputForm placeholder, NOT native matching evidence.\n'*40,encoding='utf8')
        write(root/'output/general_t3_mass_eigen_matchete/one_loop_probe'/(stem+'_one_loop_probe.json'),p)
    write(root/BASE/'bc_native_mass_eigen_complex_conditional_v21.json',v21)
    write(root/BASE/'bc_quartic_mass_readiness_conditional_v22.json',v22)
    return v21,v22

class TestV25V27(unittest.TestCase):
    def setUp(self):
        self.tmp=tempfile.TemporaryDirectory();self.addCleanup(self.tmp.cleanup)
        self.root=Path(self.tmp.name)
        fixtures(self.root)
        # Test-only replacements for existing v21/v22 module contracts.
        import Numerical.parallel_t3.wilson as ns
        m21=types.ModuleType(ns.__name__+'.bc_mass_eigen_native_bridge_v21')
        m21.FORMAT='FAKE-V21';m21.STATUS='FAKE-V21-STATUS'
        m22=types.ModuleType(ns.__name__+'.bc_mass_bridge_gate_v22')
        m22.FORMAT='FAKE-V22';m22.STATUS='FAKE-V22-STATUS'
        m22.audit=lambda root,v20_path,v21_path: v25.load(self.root/BASE/'bc_quartic_mass_readiness_conditional_v22.json')
        self.module_patches=mock.patch.dict(sys.modules,{m21.__name__:m21,m22.__name__:m22})
        self.module_patches.start();self.addCleanup(self.module_patches.stop)
        self.p21=self.root/BASE/'bc_native_mass_eigen_complex_conditional_v21.json'
        self.p22=self.root/BASE/'bc_quartic_mass_readiness_conditional_v22.json'
        self.p25=self.root/BASE/'bc_native_matchete_mass_eigen_concordance_v25.json'
        self.p26=self.root/BASE/'bc_existing_one_scale_matchete_probe_v26.json'
        self.p27=self.root/BASE/'bc_fermion_first_scalar_threshold_ledger_v27.json'
    def get25(self):return v25.audit(self.root,self.p21,self.p22)
    def get26(self):
        write(self.p25,self.get25());return v26.audit(self.root,self.p25)
    def run_full(self):
        write(self.p25,self.get25());write(self.p26,v26.audit(self.root,self.p25))
        return v27.audit(self.root,self.p22,self.p25,self.p26)
    def alter(self,model='B_-1',suffix='_matchete.json'):
        p=self.root/'output/general_t3_mass_eigen_matchete'/(model[0]+'_m1_final'+suffix)
        return p,v25.load(p)
    def test_all_three_reports_source_conditioned(self):
        report=self.run_full()
        self.assertEqual(report['models_checked'],2)
        self.assertEqual(report['full_scalar_quartic_registration_slots'],178)
        self.assertTrue(report['certification']['physical_Wilson_evolution_authorized'] is False)
        self.assertEqual(report['evidence']['native_BC_complex_CG_mass_rotations'],16)
    def test_endpoint_below_both_thresholds(self):
        result=self.run_full()
        for model,row in result['models'].items():
            self.assertGreater(row['scalar_mass_eigenstates']['X_gev'],1000)
            self.assertGreater(row['scalar_mass_eigenstates']['Y_gev'],row['scalar_mass_eigenstates']['X_gev'])
            self.assertFalse(row['source_candidate_EFT1_state_is_at_scalar_matching_scale'])
    def test_v25_all_native_cg_counts(self):
        p=self.get25()
        self.assertEqual(p['native_17_quartic_eigenfield_lagrangians_checked'],2)
        for k,row in p['models'].items():self.assertEqual(len(row['native_cg_source_symbol_reconciliation']),7)
    def test_mutated_spurion_blocked(self):
        path,doc=self.alter();doc['rotation_spurion_values_from_numerical_U']['T3MEU21']['imag']+=0.2
        write(path,doc)
        with self.assertRaisesRegex(v25.GateError,'spurions'):self.get25()
    def test_changed_y2_epsilon_factor_blocked(self):
        path,doc=self.alter();doc['native_yukawa2_epsilon_tensor']['OldTimesEpsilonCGFactor']='1'
        write(path,doc)
        with self.assertRaisesRegex(v25.GateError,'Yukawa epsilon'):self.get25()
    def test_missing_seventh_cg_blocked(self):
        path,doc=self.alter();del doc['native_additional_quartic_tensors']['Certificates']['T3ExtraQ0']
        write(path,doc)
        with self.assertRaisesRegex(v25.GateError,'seven CG'):self.get25()
    def test_native_eigenfield_status_change_blocked(self):
        path,doc=self.alter();doc['mass_basis_ordinary_check_passed']=False
        write(path,doc)
        with self.assertRaisesRegex(v25.GateError,'native check'):self.get25()
    def test_false_physical_certification_blocked(self):
        path,doc=self.alter();doc['physical_matching_performed']=True
        write(path,doc)
        with self.assertRaisesRegex(v25.GateError,'matching claim'):self.get25()
    def test_changed_parameter_source_hash_blocked(self):
        path,doc=self.alter(suffix='_parameters.json');doc['source_sha256']='0'*64
        write(path,doc)
        with self.assertRaisesRegex(v25.GateError,'parameter provenance'):self.get25()
    def test_unsupported_provenance_path_blocked(self):
        path,doc=self.alter(suffix='_parameters.json');doc['source_path']='..\\credentials.json'
        write(path,doc)
        with self.assertRaisesRegex(v25.GateError,'unsafe provenance path'):self.get25()
    def test_changed_one_step_approach_blocked(self):
        self.get26() # creates v25
        path=self.root/'output/general_t3_mass_eigen_matchete/one_loop_probe/B_m1_final_one_loop_probe.json'
        doc=v25.load(path);doc['matching_approach']='FermionFirstPhysicalEFT1';write(path,doc)
        with self.assertRaisesRegex(v25.GateError,'Weinberg extraction schema'):v26.audit(self.root,self.p25)
    def test_symbolic_probe_c5_overclaim_blocked(self):
        self.get26()
        path=self.root/'output/general_t3_mass_eigen_matchete/one_loop_probe/C_m1_final_one_loop_probe.json'
        doc=v25.load(path);doc['numerical_c5_computed']=True;write(path,doc)
        with self.assertRaisesRegex(v25.GateError,'overclaim'):v26.audit(self.root,self.p25)
    def test_missing_inputform_probe_artifact_blocked(self):
        self.get26()
        path=self.root/'output/general_t3_mass_eigen_matchete/one_loop_probe/B_m1_final_one_loop_probe_matched_eft.inputform.txt'
        path.unlink()
        with self.assertRaisesRegex(v25.GateError,'missing/empty'):v26.audit(self.root,self.p25)
    def test_cross_evidence_hash_rejected(self):
        self.run_full()
        p=v25.load(self.p25);p['status']='TAMPERED';write(self.p25,p)
        with self.assertRaises(v25.GateError):v27.audit(self.root,self.p22,self.p25,self.p26)
    def test_resume_immutable_and_verify_only(self):
        first=cli.run(self.root,resume=True)
        second=cli.run(self.root,verify_only=True)
        self.assertEqual(first,second)
        p=v25.load(self.p27);p['full_scalar_quartic_registration_slots']=177;write(self.p27,p)
        with self.assertRaisesRegex(v25.GateError,'differs'):cli.run(self.root,resume=True)
    def test_verify_only_missing_rejected(self):
        with self.assertRaisesRegex(v25.GateError,'verify-only'):cli.run(self.root,verify_only=True)
    def test_cli_input_file_not_overwritten(self):
        before=self.p21.read_bytes();cli.run(self.root,resume=True)
        self.assertEqual(before,self.p21.read_bytes())
    def test_nonfinite_complex_rejected(self):
        with self.assertRaises(v25.GateError):v25.cm({'real':[[1,float('nan')],[0,1]],'imag':[[0,0],[0,0]]},2,'bad')
    def test_path_absolute_or_collision_rejected(self):
        for key in ('C:\\x','/root/x','a/../x','a//b','a/./b'):
            with self.subTest(key=key),self.assertRaises(v25.GateError):v25.relpath(key)

class TestLiveLocalRepository(unittest.TestCase):
    def test_current_repository_native_sources_when_present(self):
        project=Path(__file__).resolve().parents[1]
        paths=[project/BASE/'bc_native_mass_eigen_complex_conditional_v21.json',
               project/BASE/'bc_quartic_mass_readiness_conditional_v22.json',
               project/'output/general_t3_mass_eigen_matchete/B_m1_final_matchete.json',
               project/'output/general_t3_mass_eigen_matchete/C_m1_final_matchete.json']
        if not all(x.is_file() for x in paths):
            self.skipTest('Installed native Matchete + completed local v21/v22 not mounted')
        output=v25.audit(project,paths[0],paths[1])
        self.assertEqual(output['models_checked'],2)

if __name__=='__main__':unittest.main()
