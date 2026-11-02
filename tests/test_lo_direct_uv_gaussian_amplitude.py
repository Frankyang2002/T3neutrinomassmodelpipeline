"""All tests are mathematical and source-contract fixtures, NOT native proof."""
from __future__ import annotations

import copy
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
import numpy as np

ROOT=Path(__file__).resolve().parents[1]
TESTS=Path(__file__).resolve().parent
for p in (ROOT,TESTS):
    if str(p) not in sys.path:sys.path.insert(0,str(p))
from Numerical.direct_one_loop.core import encoded,_complex
from Numerical.lo_certification.direct_uv.gaussian_amplitude import (
    AmplitudeError,CHANNELS,batch,dagger,exterior,explicit_four_spinor_witness,
    gaussian_tensors,phase_invariant_source_orientation,_eft1_boundary_gate)
from Numerical.lo_certification.direct_uv.certification import CLASSES,all_models
from Numerical.lo_certification.direct_uv.phase_portal import rephase
from test_lo_direct_uv_certification import synthetic_all16,toy_B
from test_lo_direct_uv_spinor_mass import fixture


class GaussianAmplitude(unittest.TestCase):
    def test_chosen_clifford_C_and_both_yukawas(self):
        r=explicit_four_spinor_witness()
        self.assertTrue(r['negative_wrong_C_sign_rejected'])
        self.assertEqual(r['exact_identities']['BarL_F'],'dagger(l epsilon xi)')
        self.assertFalse(r['relative_Matchete_RGBeta_signed_mass_and_yukawa_dictionary_certified'])

    def test_majorana_half_and_one_dirac_mass(self):
        r=explicit_four_spinor_witness()
        self.assertTrue(r['majorana_repeated_grassmann_spinor_factor_two_checked'])
        self.assertEqual(r['exact_identities']['BarF_F'],'chi epsilon xi + h.c.')

    def test_independent_complex_mass_and_phase(self):
        for majorana in (False,True):
            r=phase_invariant_source_orientation(majorana)
            self.assertEqual(r['sum_charge'],[0,0,0,0])
            self.assertNotEqual(r['wrong_lambda_conjugation_charge'],[0,0,0,0])

    def test_exterior_anticommutation_and_no_double_spinor(self):
        a={('ch0',):1};b={('ch1',):1}
        self.assertEqual(exterior(a,a),{})
        self.assertEqual(exterior(a,b),{('ch0','ch1'):1})
        self.assertEqual(exterior(b,a),{('ch0','ch1'):-1})
        self.assertEqual(dagger(dagger(exterior(a,b))),exterior(a,b))


    def test_dirac_chirality_cross_and_no_automatic_diagonal_mass(self):
        d=synthetic_all16('A_-4')
        M=_complex(d['fermion_mass_gev'],'M')
        nc=M.shape[0]//2
        self.assertTrue(np.allclose(M[:nc,:nc],0))
        self.assertTrue(np.allclose(M[nc:,nc:],0))
        v=gaussian_tensors(d,'A',-4)
        # The mass propagator connects opposite left Weyl chiralities.
        self.assertGreater(v['heavy_fermion_state_count'],3)
        self.assertLess(v['LLSS_exchange_symmetry_max_relative_error'],1e-9)

    def test_fermion_release_scale_guard(self):
        with tempfile.TemporaryDirectory() as t:
            root=Path(t);p=root/'output/general_observable/lo_fermion_boundary/releases/source_test/B_m1/boundary.json'
            p.parent.mkdir(parents=True)
            # Boundary directory is B_m1, but JSON uses the signed model ID.
            p.write_text(json.dumps({'model':'B_-1',
              'mu_F_gev':1e5,
              'field_content':{'heavy_fermions_integrated_out':'all at mu_F',
                               'heavy_scalar_multiplets':'both retained in EFT1'},
              'spinor_normalization':{'certified':False}}))
            gate=_eft1_boundary_gate(root,'B_-1',1e7)
            self.assertFalse(gate['same_renormalization_scale'])
            self.assertFalse(gate['compared_LLSS_tensors'])
            self.assertEqual(gate['eft1_threshold_mu_F_gev'],1e5)
            gate2=_eft1_boundary_gate(root,'B_-1',1e5)
            self.assertTrue(gate2['same_renormalization_scale'])
            self.assertFalse(gate2['compared_LLSS_tensors'])
            self.assertIn('PHYSICAL',gate2['status'])

    def test_all_16_computed_all_four_scalar_sectors_and_pairs(self):
        for short,_ in all_models():
            with self.subTest(model=short):
                v=gaussian_tensors(synthetic_all16(short),short[0],int(short[2:]),include_tensors=False)
                self.assertEqual(len(v['tensor_channels']),4)
                self.assertEqual(set(v['tensor_channels']),{ch for _,_,ch in CHANNELS})
                self.assertEqual(set(v['tensor_channels']['YZ']['source_pair_tensors_gev_minus1']),
                                 {'y1_HC x y1_HC','y1_HC x y2_HC','y2_HC x y2_HC'})
                self.assertLess(v['heavy_inverse_direct_vs_takagi_relative_residual'],1e-9)
                self.assertLess(v['LLSS_exchange_symmetry_max_relative_error'],1e-9)
                self.assertLess(v['LLSS_charge_violation_max_relative_error'],1e-9)
                self.assertFalse(v['fermion_threshold_hard_one_loop_LLHH_included'])

    def test_majorana_cross_term_has_one_not_spurious_two(self):
        d=toy_B()
        v=gaussian_tensors(d,'B',-1)
        inv=_complex(d['fermion_mass_gev'],'M')
        Z=_complex(d['vertices'][0]['conjugate'],'Z')
        Y=_complex(d['vertices'][1]['field'],'Y')
        expected=np.einsum('ais,ij,bjt->abst',Z,np.linalg.inv(inv),Y)
        actual=_complex(v['tensor_channels']['ZY']['source_pair_tensors_gev_minus1']['y1_HC x y2_HC'],'K')
        # Z comes from y1_HC, Y from y2_HC: convention ZY (conjugate,field)
        self.assertLess(np.linalg.norm(actual-expected)/max(np.linalg.norm(expected),1e-100),1e-12)

    def test_gaussian_direct_schur_complement_with_complex_offdiagonal_mass(self):
        d=toy_B()
        M=np.array([[101+3j,1+2j,0],[1+2j,118-5j,2],[0,2,139+1j]],complex)
        d['fermion_mass_gev']=encoded(M)
        v=gaussian_tensors(d,'B',-1)
        self.assertLess(v['gaussian_equations_of_motion_relative_residual'],1e-11)
        self.assertLess(v['heavy_inverse_direct_vs_takagi_relative_residual'],1e-10)

    def test_scalar_charge_negative_control(self):
        d=toy_B()
        z=_complex(d['vertices'][0]['conjugate'],'Z')
        z[0,0,1]=.8
        d['vertices'][0]['conjugate']=encoded(z)
        with self.assertRaises(ValueError):gaussian_tensors(d,'B',-1)

    def test_fermion_mass_singularity_negative_control(self):
        d=toy_B()
        d['fermion_mass_gev']=encoded(np.diag([1.,2.,0.]))
        with self.assertRaises(AmplitudeError):gaussian_tensors(d,'B',-1)

    def test_pure_field_phase_covariance_all_tensor_sectors(self):
        d=toy_B()
        pf=np.exp(1j*np.array([.19,-.32,.47]))
        ps=np.exp(1j*np.array([.71,-.24,.09,-.45]))
        rephased=rephase(d,pf,ps)
        a=gaussian_tensors(d,'B',-1)
        b=gaussian_tensors(rephased,'B',-1)
        for tag1,tag2,ch in CHANNELS:
            x=_complex(a['tensor_channels'][ch]['full_tensor_gev_minus1'],'K')
            y=_complex(b['tensor_channels'][ch]['full_tensor_gev_minus1'],'Kp')
            p=ps if tag1=='field' else ps.conj()
            q=ps if tag2=='field' else ps.conj()
            expected=x*p[None,None,:,None]*q[None,None,None,:]
            self.assertLess(np.linalg.norm(y-expected)/max(np.linalg.norm(x),1e-100),1e-10)

    def test_crossing_two_boson_orientation_sectors(self):
        v=gaussian_tensors(toy_B(),'B',-1)
        zy=_complex(v['tensor_channels']['ZY']['full_tensor_gev_minus1'],'ZY')
        yz=_complex(v['tensor_channels']['YZ']['full_tensor_gev_minus1'],'YZ')
        self.assertLess(np.linalg.norm(zy-yz.transpose(1,0,3,2)),1e-16)

    def test_summary_only_keeps_source_pair_norms(self):
        r=gaussian_tensors(toy_B(),'B',-1,include_tensors=False)
        self.assertIsNone(r['tensor_channels']['YZ']['full_tensor_gev_minus1'])
        self.assertIn('norm_gev_minus1',r['tensor_channels']['YZ']['source_pair_tensors_gev_minus1']['y1_HC x y2_HC'])

    def test_source_provenance_guard_and_all_16_batch(self):
        with tempfile.TemporaryDirectory() as t:
            root=Path(t);fixture(root)
            indir=root/'output/general_observable/lo_direct_one_loop/provisional_inputs'
            indir.mkdir(parents=True)
            for short,key in all_models():
                kind,alpha=short[0],int(short[2:]);d1,d2,df,_=CLASSES[kind]
                bp=root/'Numerical/benchmarks/fitted_reference'/(key+'.json')
                bp.parent.mkdir(parents=True,exist_ok=True)
                bp.write_text(json.dumps({'representation':{'d_s1':d1,'d_s2':d2,'d_f':df,'alpha':alpha},
                         'base_state':{'mu_gev':1e7}}))
                inp=synthetic_all16(short)
                import hashlib
                inp['provenance']={'benchmark_sha256':hashlib.sha256(bp.read_bytes()).hexdigest()}
                (indir/(short+'.json')).write_text(json.dumps(inp))
            result=batch(root,include_tensors=False)
            self.assertEqual(result['tree_amplitude_computed_models'],16)
            self.assertEqual(result['blocked_models'],0)
            self.assertEqual(result['physical_certified_models'],0)
            self.assertEqual(len(result['models']),16)
            # One bad benchmark binds no non-current synthetic input.
            bp=root/'Numerical/benchmarks/fitted_reference'/(next(iter(all_models()))[1]+'.json')
            bp.write_text('{}')
            wrong=batch(root,include_tensors=False)
            self.assertEqual(wrong['blocked_models'],1)
            self.assertEqual(wrong['tree_amplitude_computed_models'],15)

    def test_successful_cli_all16_writes_owned_tree_and_csv(self):
        with tempfile.TemporaryDirectory() as t:
            root=Path(t);fixture(root)
            indir=root/'output/general_observable/lo_direct_one_loop/provisional_inputs'
            indir.mkdir(parents=True)
            for short,key in all_models():
                kind,alpha=short[0],int(short[2:]);d1,d2,df,_=CLASSES[kind]
                bp=root/'Numerical/benchmarks/fitted_reference'/(key+'.json')
                bp.parent.mkdir(parents=True,exist_ok=True)
                bp.write_text(json.dumps({'representation':{'d_s1':d1,'d_s2':d2,'d_f':df,'alpha':alpha},
                    'base_state':{'mu_gev':1e7}}))
                inp=synthetic_all16(short)
                import hashlib
                inp['provenance']={'benchmark_sha256':hashlib.sha256(bp.read_bytes()).hexdigest()}
                (indir/(short+'.json')).write_text(json.dumps(inp))
            out=root/'output/general_observable/lo_certification_direct_uv/full_01'
            script=ROOT/'scripts/RunLODirectUVAmplitudeBridge.py'
            cmd=[sys.executable,str(script),'--project-root',str(root),'--output-dir',str(out),
                 '--summary-only']
            process=subprocess.run(cmd,capture_output=True,text=True)
            self.assertEqual(process.returncode,0,process.stderr+process.stdout)
            status=json.loads((out/'all_16_gaussian_amplitude_status.json').read_text())
            self.assertEqual(status['tree_amplitude_computed_models'],16)
            self.assertEqual(status['blocked_models'],0)
            self.assertEqual(status['physical_certified_models'],0)
            self.assertTrue((out/'all_16_gaussian_amplitude_summary.csv').is_file())
            self.assertEqual(len(list((out/'models').glob('*.json'))),16)
            self.assertEqual(subprocess.run(cmd,capture_output=True,text=True).returncode,2)

    def test_cli_guards_and_no_overwrite(self):
        with tempfile.TemporaryDirectory() as t:
            root=Path(t);fixture(root)
            out=root/'output/general_observable/lo_certification_direct_uv/test_amplitude_01'
            script=ROOT/'scripts/RunLODirectUVAmplitudeBridge.py'
            cmd=[sys.executable,str(script),'--project-root',str(root),'--output-dir',str(out),'--summary-only']
            # No benchmark files: all 16 are blocked, report still written.
            a=subprocess.run(cmd,capture_output=True,text=True)
            self.assertEqual(a.returncode,2,a.stderr)
            self.assertTrue((out/'all_16_gaussian_amplitude_status.json').exists())
            self.assertEqual(subprocess.run(cmd,capture_output=True,text=True).returncode,2)
            self.assertEqual(subprocess.run([sys.executable,str(script),'--project-root',str(root),
                           '--output-dir',str(root/'outside')],capture_output=True,text=True).returncode,2)


if __name__=='__main__':unittest.main()
