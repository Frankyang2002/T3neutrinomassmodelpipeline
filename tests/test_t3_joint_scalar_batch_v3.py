"""Complete synthetic/full-SU2 parity tests; do not infer a physical GammaCC sign."""
from __future__ import annotations

from dataclasses import replace
from pathlib import Path
import json
import shutil
import tempfile
import unittest

import numpy as np

from Numerical.joint_scalar_matching.adapter import ALL_MODELS
from Numerical.joint_scalar_matching.fermion_boundary_bridge import digest,load_source,slug
from Numerical.joint_scalar_matching.large_batch import load_background
from Numerical.joint_scalar_matching.source_su2_batch import (
    FORMAT,STATUS,SUMMARY_FORMAT,SUMMARY_STATUS,
    EXPECTED_CONTRACTION, EXPECTED_WEINBERG,
    audit_signed_evidence,independent_source_beta_check,
    validate_su2_manifest,validate_su2_response)
from scripts.RunT3JointScalarBatchV3 import SAFE,main
from test_t3_joint_scalar_fermion_boundary_bridge import make_release
from test_t3_joint_scalar_integrated_batch import (
    clean_synthetic_source,make_background,write_pivot_kernel)


def complex_doc(a):
    a=np.asarray(a,dtype=complex)
    return {'real':a.real.tolist(),'imag':a.imag.tolist()}


def write_su2(path,source):
    n=len(source.labels)
    per_source=[]
    for s,idx in enumerate(source.source_indices):
        pairs=[]
        for a in range(n):
            for b in range(a,n):
                ar=source.by_source[s,...,a,b]
                if not np.any(ar!=0): continue
                pairs.append({'heavy_pair_one_based':[a+5,b+5],
                              'ordered_pair_multiplicity':1 if a==b else 2,
                              'LLSS_C_pqij_per_NGamma_GeV_inv':complex_doc(ar)})
        per_source.append({'matchete_source_index':idx,'scalar_pairs':pairs})
    data={'format':FORMAT,'status':STATUS,'model':source.model,
          'source_boundary_sha256':source.hash,'mu_F_gev':source.mu_f_gev,
          'mu_S_gev':None,'source_count':len(source.source_indices),
          'N_GammaCC_to_Weyl':None,
          'absolute_spinor_normalization_certified':False,
          'physical_C5_computed':False,
          'complete_one_loop_operator_inventory_certified':False,
          'one_loop_fermion_hard_matching_included':False,
          'finite_joint_scalar_matching_included':False,
          'full_su2_contraction_convention':EXPECTED_CONTRACTION,
          'weinberg_tensor_convention':EXPECTED_WEINBERG,
          'source_resolved_blocks':per_source}
    path.parent.mkdir(parents=True,exist_ok=True)
    path.write_text(json.dumps(data,indent=2),encoding='utf-8')
    return path


def write_manifest(root,release,models):
    rows=[]
    for model in models:
        src=load_source(release/slug(model)/'boundary.json')
        write_su2(root/slug(model)/'full_su2_response.json',src)
        rows.append({'model':model,'source_count':len(src.source_indices),
                     'source_boundary_sha256':src.hash,'status':STATUS})
    o={'format':SUMMARY_FORMAT,'status':SUMMARY_STATUS,
       'model_count':16,'fermion_release_summary_sha256':digest(release/'summary.json'),
       'quartic_input_kind':'none','physical_C5_computed':False,
       'source_spinor_normalization_certified':False,'models':rows}
    (root/'summary.json').write_text(json.dumps(o,indent=2),encoding='utf-8')
    return root


class FullSU2BatchTests(unittest.TestCase):
    def test_16_model_full_su2_parity_including_mixed_and_self_split(self):
        with tempfile.TemporaryDirectory() as td:
            root=Path(td)
            release=make_release(root)
            path=write_manifest(root/'su2',release,ALL_MODELS)
            manifest=validate_su2_manifest(path,release,models=ALL_MODELS)
            self.assertEqual(len(manifest),16)
            for model in ALL_MODELS:
                with self.subTest(model=model):
                    source=load_source(release/slug(model)/'boundary.json')
                    report=validate_su2_response(path/slug(model)/'full_su2_response.json',source)
                    self.assertTrue(report['full_two_lepton_SU2_component_parity'])
                    self.assertFalse(report['physical_GammaCC_normalization_certified'])
                    self.assertEqual(report['source_count_verified'],2)

    def test_rejects_signed_su2_component_corruption(self):
        with tempfile.TemporaryDirectory() as td:
            src,_,_=clean_synthetic_source('B_-1',Path(td))
            p=write_su2(Path(td)/'SU2.json',src)
            obj=json.loads(p.read_text())
            pair=obj['source_resolved_blocks'][0]['scalar_pairs'][0]
            pair['LLSS_C_pqij_per_NGamma_GeV_inv']['real'][0][0][0][0]+=1e-10
            p.write_text(json.dumps(obj))
            with self.assertRaisesRegex(ValueError,'differs from fermion boundary'):
                validate_su2_response(p,src)

    def test_rejects_pair_multiplicity_and_missing_pair(self):
        with tempfile.TemporaryDirectory() as td:
            src,_,_=clean_synthetic_source('E_-2',Path(td))
            p=write_su2(Path(td)/'SU2.json',src)
            obj=json.loads(p.read_text())
            obj['source_resolved_blocks'][0]['scalar_pairs'][0]['ordered_pair_multiplicity']=999
            p.write_text(json.dumps(obj))
            with self.assertRaisesRegex(ValueError,'multiplicity'):
                validate_su2_response(p,src)
            p=write_su2(p,src)
            obj=json.loads(p.read_text())
            obj['source_resolved_blocks'][0]['scalar_pairs'].pop(0)
            p.write_text(json.dumps(obj))
            with self.assertRaisesRegex(ValueError,'omitted or added scalar pairs'):
                validate_su2_response(p,src)

    def test_zero_at_benchmark_PL_source_and_inter_source_cancellation(self):
        with tempfile.TemporaryDirectory() as td:
            src,_,_=clean_synthetic_source('E_+0',Path(td))
            by=np.stack((src.raw,-src.raw,np.zeros_like(src.raw)))
            source=replace(src,raw=np.zeros_like(src.raw),by_source=by,
                           source_indices=(10,11,12))
            path=write_su2(Path(td)/'SU2.json',source)
            report=validate_su2_response(path,source)
            self.assertGreater(report['distinct_pair_count_before_cancellation'],0)
            self.assertEqual(report['distinct_pair_count_after_cancellation'],0)
            self.assertEqual(report['source_count_verified'],3)

    def test_manifest_refuses_other_source_release(self):
        with tempfile.TemporaryDirectory() as td:
            root=Path(td);rel=make_release(root)
            su2=write_manifest(root/'su2',rel,ALL_MODELS)
            doc=json.loads((su2/'summary.json').read_text())
            doc['fermion_release_summary_sha256']='0'*64
            (su2/'summary.json').write_text(json.dumps(doc))
            with self.assertRaisesRegex(ValueError,'manifest disagrees'):
                validate_su2_manifest(su2,rel,models=ALL_MODELS)

    def test_full_SU2_unordered_pair_beta_rotation_16(self):
        with tempfile.TemporaryDirectory() as td:
            root=Path(td)
            for model in ALL_MODELS:
                with self.subTest(model=model):
                    src,q,channels=clean_synthetic_source(model,root/model)
                    n=len(src.labels)
                    rot=np.eye(n)
                    rot[0,0]=0.6;rot[0,1]=0.8
                    rot[1,0]=-0.8;rot[1,1]=0.6
                    bg=load_background(make_background(root/model/'bg',src,q,channels,U=rot),src)
                    result=independent_source_beta_check(src,bg)
                    self.assertTrue(result['full_SU2_unordered_pair_beta_parity'])
                    self.assertLess(result['weinberg_projection_residual'],1e-10)

    def test_spinor_evidence_missing_does_not_certify(self):
        with tempfile.TemporaryDirectory() as td:
            root=Path(td);release=make_release(root)
            report=audit_signed_evidence(root/'absent_native.json',root/'absent_signed.json',release)
            self.assertFalse(report['physical_spinor_normalization_certified'])
            self.assertIsNone(report['absolute_matchete_to_weyl_factor'])
            self.assertGreaterEqual(len(report['blocking_items']),3)

    def test_signed_evidence_correct_source_not_a_normalization_certificate(self):
        with tempfile.TemporaryDirectory() as td:
            root=Path(td);release=make_release(root)
            signed=root/'signed.json';native=root/'native.json'
            signed.write_text(json.dumps({
                'format':'t3-signed-zero-momentum-grassmann-diagnostic-v1',
                'status':'ConditionalSignedContactOnly_NoPhysicalCertification',
                'input_release_summary_sha256':digest(release/'summary.json'),
                'assumed_GammaCC_phase_certified':False}))
            native.write_text(json.dumps({
                'format':'t3-native-matchete-fr-evidence-v1',
                'status':'Blocked','blocked_phase':'EFT1 Feynman rules',
                'reason':'No nonempty contact vertex',
                'physical_master_wilson_normalization_certified':False}))
            report=audit_signed_evidence(native,signed,release)
            self.assertFalse(report['physical_spinor_normalization_certified'])
            self.assertIn('native signed UV-to-EFT1 four-point amplitude not certified',
                          report['blocking_items'])

    def test_cli_full_16_without_background_exits_success_and_does_not_invent_C5(self):
        with tempfile.TemporaryDirectory() as td:
            root=Path(td);release=make_release(root)
            su2=write_manifest(root/'su2',release,ALL_MODELS)
            kernels=root/'pivots'
            for model in ALL_MODELS:
                src=load_source(release/slug(model)/'boundary.json')
                write_pivot_kernel(kernels/slug(model)/'quartic_pivot_kernel.json',src)
            dest=SAFE/'_batch_v3_unit_test_sandbox'
            shutil.rmtree(dest,ignore_errors=True)
            try:
                rc=main(['--release-dir',str(release),'--su2-dir',str(su2),
                          '--pivot-dir',str(kernels),'--native-evidence',str(root/'missing.json'),
                          '--signed-evidence',str(root/'missing2.json'),
                          '--output-dir',str(dest)])
                self.assertEqual(rc,0)
                o=json.loads((dest/'summary.json').read_text())
                self.assertEqual(o['models_source_pivot_su2_verified'],16)
                self.assertEqual(o['conditional_responses_computed'],0)
                self.assertFalse(o['physical_C5_certified'])
                self.assertEqual(len(o['models']),16)
                self.assertTrue(all(not r['background_at_muS_available'] for r in o['models']))
                self.assertEqual(list(dest.glob('conditional_C5_per_GammaCC/*.json')),[])
            finally:
                shutil.rmtree(dest,ignore_errors=True)

    def test_cli_blocks_when_su2_source_is_damaged(self):
        with tempfile.TemporaryDirectory() as td:
            root=Path(td);release=make_release(root)
            su2=write_manifest(root/'su2',release,ALL_MODELS)
            p=su2/'B_m1/full_su2_response.json'
            obj=json.loads(p.read_text())
            obj['source_resolved_blocks'][0]['matchete_source_index']=999
            p.write_text(json.dumps(obj))
            kern=root/'pivots'
            src=load_source(release/'B_m1/boundary.json')
            write_pivot_kernel(kern/'B_m1/quartic_pivot_kernel.json',src)
            dest=SAFE/'_batch_v3_unit_test_sandbox'
            shutil.rmtree(dest,ignore_errors=True)
            try:
                rc=main(['--release-dir',str(release),'--su2-dir',str(su2),
                         '--pivot-dir',str(kern), '--model','B_-1',
                         '--native-evidence',str(root/'no.json'),
                         '--signed-evidence',str(root/'no2.json'),
                         '--output-dir',str(dest)])
                self.assertEqual(rc,2)
                report=json.loads((dest/'summary.json').read_text())
                self.assertEqual(report['blocked_models'],1)
                self.assertIn('Matchete PL source list order',report['models'][0]['error'])
            finally:
                shutil.rmtree(dest,ignore_errors=True)


if __name__=='__main__':unittest.main()
