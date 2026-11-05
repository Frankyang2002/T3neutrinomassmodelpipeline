"""Isolated tests for uniform-0.01 thesis running-stage heatmap campaign."""
from __future__ import annotations
import tempfile
import unittest
from dataclasses import dataclass
from pathlib import Path
import hashlib
import json
from unittest.mock import patch
import numpy as np
from scripts import RunT3Uniform001StageRunning as r


def cases():
    base={
      'delta_m21_sq_ev2':7.537e-5,'delta_m3l_sq_ev2':2.511e-3,
      'sin2_theta12':0.3088,'sin2_theta13':0.02248,'sin2_theta23':0.47,
    }
    mat={'real':[[1e-15,0.,0.],[0.,2e-15,0.],[0.,0.,3e-15]],
         'imag':[[0.,0.,0.],[0.,0.,0.],[0.,0.,0.]]}
    return {key:{'prediction':{k:v*scale for k,v in base.items()},
                 'C5_low_gev_inv':mat,'chi2':0.0,'neutrino_masses_ev':[.001,.009,.05]}
            for key,scale in [('111',1.),('011',.9),('101',1.1),('110',.99)]}


class RunTests(unittest.TestCase):
    def test_model_contract(self):
        self.assertEqual(len(r.MODELS),16)
        self.assertEqual(r.select_models('pilot'),('A_-4','B_-1'))
        self.assertEqual(r.select_models('all'),r.MODELS)
        with self.assertRaises(ValueError):r.select_models('A_-4,A_-4')
        with self.assertRaises(ValueError):r.select_models('Z_-7')

    def test_all_quartics_policy_source(self):
        self.assertEqual(r.uniform.QUARTIC,0.01)
        self.assertEqual(set(r.CASES.values()),{'011','101','110'})

    def test_frozen_state_updates_only_scale(self):
        @dataclass(frozen=True)
        class State:
            mu_gev:float
            lambdaT3:float
            def validated(self):return self
        a=State(1e7,.1)
        b=r.frozen_state(a,1e5)
        self.assertEqual(a.mu_gev,1e7)
        self.assertEqual(b.mu_gev,1e5)
        self.assertEqual(b.lambdaT3,.1)

    def test_four_case_contract(self):
        c=cases()
        self.assertEqual(len(r.stage_rows('A_-4',c)),15)
        with self.assertRaises(ValueError):r.stage_rows('A_-4',{'111':c['111']})

    def test_signed_percentage_denominator_is_off(self):
        self.assertAlmostEqual(r.percent_change(1.0,0.8),25.0)
        self.assertAlmostEqual(r.percent_change(0.8,1.0),-20.0)
        with self.assertRaises(ValueError):r.percent_change(1.0,0.0)

    def test_stage_rows_stay_separate(self):
        data=r.stage_rows('A_-4',cases())
        one=[x for x in data if x['observable']=='sin2_theta12']
        self.assertEqual([x['switch'] for x in one],list(r.STAGES))
        self.assertEqual([x['off_case'] for x in one],['011','101','110'])
        self.assertAlmostEqual(one[0]['relative_change_percent'],100/9)
        self.assertLess(one[1]['relative_change_percent'],0)
        self.assertIn('conditional-one-scale',one[2]['matching_scope'])

    def test_c5_shape_and_symmetry(self):
        matrix=np.eye(3)*1e-15
        self.assertEqual(r.require_c5(matrix).shape,(3,3))
        broken=np.eye(3)*1e-15
        broken[0,1]=1e-15
        with self.assertRaises(ValueError):r.require_c5(broken)
        with self.assertRaises(ValueError):r.require_c5(np.zeros((2,2)))

    def test_reference_baseline_check(self):
        c=cases()
        ref={'prediction':c['111']['prediction'],
             'numerical_details':{'c5_at_low_scale':c['111']['C5_low_gev_inv']}}
        checked=r.verify_baseline('A_-4',c,ref)
        self.assertEqual(checked['C5_low_relative_difference'],0)
        ref['prediction']=dict(ref['prediction'],sin2_theta13=.03)
        with self.assertRaises(ValueError):r.verify_baseline('A_-4',c,ref)

    def test_reference_c5_mismatch_blocks(self):
        c=cases()
        bad={'real':[[1e-14,0.,0.],[0.,2e-15,0.],[0.,0.,3e-15]],
             'imag':[[0.,0.,0.],[0.,0.,0.],[0.,0.,0.]]}
        ref={'prediction':c['111']['prediction'],
             'numerical_details':{'c5_at_low_scale':bad}}
        with self.assertRaisesRegex(ValueError,'C5 baseline'):
            r.verify_baseline('D_-2',c,ref)

    def test_plot_refuses_incomplete_models(self):
        with tempfile.TemporaryDirectory() as tmp:
            with self.assertRaises(ValueError):
                r.write_heatmap_csv(Path(tmp),('A_-4',))

    def test_complete_csv_preserves_schema_and_unique_rows(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp)
            for model in r.MODELS:
                r.save(root/'models'/r.tag(model)/'result.json',{
                    'schema':r.SCHEMA,'status':'Computed','model':model,
                    'extra_BSM_quartics_uv':0.01,'conjugate_yukawas_uv_zero':True,
                    'baseline_reproduction':{'C5_low_relative_difference':0.0},
                    'cases':cases()})
            csv_path=r.write_heatmap_csv(root,r.MODELS)
            self.assertTrue(csv_path.is_file())
            import csv
            with csv_path.open(encoding='utf-8',newline='') as f:
                rows=list(csv.DictReader(f))
            self.assertEqual(len(rows),240)
            self.assertEqual(len({(x['model'],x['switch'],x['observable']) for x in rows}),240)

    def test_portable_provenance_accepts_only_line_ending_conversion(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp)
            src=root/'saved'/'uv_beta.json'
            src.parent.mkdir(parents=True)
            original=b'{"lambdaT3":0.1,\r\n"mass":1000}\r\n'
            src.write_bytes(original.replace(b'\r\n',b'\n'))
            manifest=root/'merged.json'
            manifest.write_bytes(b'{"accepted":16}\n')
            row={'provenance':{
                 'source_sha256':{'uv_beta':hashlib.sha256(original).hexdigest()},
                 'source_paths':{'uv_beta':'saved/uv_beta.json'},
                 'accepted_manifest_sha256':hashlib.sha256(manifest.read_bytes()).hexdigest(),
                 }}
            with patch.object(r.campaign,'ROOT',root):
                result=r.validate_portable_source_provenance(
                    row,{'uv_beta':src},manifest)
                self.assertEqual(result,['uv_beta'])
                src.write_bytes(original)
                self.assertEqual(r.validate_portable_source_provenance(
                    row,{'uv_beta':src},manifest),[])

    def test_portable_provenance_rejects_physics_change_and_wrong_path(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp)
            src=root/'uv.json'
            src.write_bytes(b'{"mass":1000}\r\n')
            manifest=root/'merged.json'
            manifest.write_bytes(b'manifest\n')
            row={'provenance':{
                'source_sha256':{'uv_beta':hashlib.sha256(src.read_bytes()).hexdigest()},
                'source_paths':{'uv_beta':'uv.json'},
                'accepted_manifest_sha256':hashlib.sha256(manifest.read_bytes()).hexdigest()}}
            with patch.object(r.campaign,'ROOT',root):
                src.write_bytes(b'{"mass":1100}\n')
                with self.assertRaisesRegex(ValueError,'content differs'):
                    r.validate_portable_source_provenance(row,{'uv_beta':src},manifest)
                src.write_bytes(b'{"mass":1000}\n')
                row['provenance']['source_paths']={'uv_beta':'other.json'}
                with self.assertRaisesRegex(ValueError,'paths changed'):
                    r.validate_portable_source_provenance(row,{'uv_beta':src},manifest)
                row['provenance']['source_paths']={'uv_beta':'uv.json'}
                manifest.write_bytes(b'changed\n')
                with self.assertRaisesRegex(ValueError,'manifest content changed'):
                    r.validate_portable_source_provenance(row,{'uv_beta':src},manifest)

    def test_portable_reference_fallback_retains_acceptance_and_yukawa_checks(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp)
            folder=root/'B_m1';folder.mkdir()
            src=root/'beta.json'
            src.write_bytes(b'{"b":1}\n')
            old_hash=hashlib.sha256(b'{"b":1}\r\n').hexdigest()
            manifest=root/'merged.json';manifest.write_bytes(b'accepted\n')
            matrices={'real':[[0.0]*3 for _ in range(3)],
                      'imag':[[0.0]*3 for _ in range(3)]}
            cfg={'base_state':{'ordinary':{'y1':matrices,'y2':matrices,
                                          'lambdaT3':{'real':0.1,'imag':0.0}}},
                 'thesis_uniform_scalar_benchmark':{'schema':r.uniform.SCHEMA}}
            row={'schema':r.campaign.SCHEMA,'stage':'reference','model':'B_-1',
                 'status':'AcceptedApproxFit','chi2':0.01,'lambdaT3':0.1,
                 'fitted_yukawas':{'y1':matrices,'y2':matrices},
                 'provenance':{
                     'source_sha256':{'uv_beta':old_hash},
                     'source_paths':{'uv_beta':'beta.json'},
                     'accepted_manifest_sha256':hashlib.sha256(manifest.read_bytes()).hexdigest()}}
            r.save(folder/'result.json',row)
            r.save(folder/'fitted_config.json',cfg)
            kwargs={'model':'B_-1','source_paths':{'uv_beta':src},
                    'manifest_path':manifest,'original_config':{},'chi2_cutoff':5.86}
            with patch.object(r.campaign,'ROOT',root), \
                 patch.object(r.campaign,'validated_reference',
                       side_effect=ValueError('Original accepted-source fingerprints changed')), \
                 patch.object(r.campaign,'check_uniform_reference',return_value=None):
                accepted,refcfg=r.validated_reference_portable(folder,**kwargs)
                self.assertEqual(accepted['model'],'B_-1')
                self.assertEqual(refcfg['base_state']['ordinary']['y1'],matrices)
                row['status']='PoorFit'
                r.save(folder/'result.json',row)
                with self.assertRaisesRegex(ValueError,'not accepted'):
                    r.validated_reference_portable(folder,**kwargs)

    def test_invalid_benchmark_not_rendered(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp)
            for model in r.MODELS:
                r.save(root/'models'/r.tag(model)/'result.json',{
                    'schema':r.SCHEMA,'status':'Computed','model':model,
                    'extra_BSM_quartics_uv':.02,'conjugate_yukawas_uv_zero':True,
                    'baseline_reproduction':{},'cases':cases()})
            with self.assertRaisesRegex(ValueError,'not the requested'):
                r.write_heatmap_csv(root,r.MODELS)

if __name__=='__main__':unittest.main()
