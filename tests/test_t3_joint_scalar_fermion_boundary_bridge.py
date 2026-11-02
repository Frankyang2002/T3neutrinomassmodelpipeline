"""Synthetic schema/regression tests. Do NOT certify physical Matchete spinor factors."""
from __future__ import annotations

from hashlib import sha256
import json
from pathlib import Path
import tempfile
import unittest

import numpy as np

from Numerical.joint_scalar_matching.adapter import ALL_MODELS,CONVENTION,load_input,save_input
from Numerical.joint_scalar_matching.fermion_boundary_bridge import (
    ENDPOINT_FORMAT, RELEASE_STATUS, SOURCE_FORMAT, SOURCE_STATUS,
    canonical, digest, load_eft1_endpoint, load_release_manifest,
    load_source, orbit,slug,write_provisional_source)
from Numerical.joint_scalar_matching.fixtures import make_fixture
from scripts.RunT3JointScalarMatching import compute
from scripts.BridgeT3JointScalarEFT1 import main as bridge_main


def _layout(model: str):
    d1,d2={'A':(1,3),'B':(2,2),'C':(2,2),'D':(3,1),'E':(3,3)}[model[0]]
    n=2*(d1+d2)
    if model in ('B_-1','C_-1'):
        return n,{'NewScalarX':5,'NewScalarY':9},{'NewScalarX':2,'NewScalarY':2}
    if model in ('A_+0','D_+0','E_+0'):
        return n,{'T3Scalar1Minus':5,'T3Scalar1Plus':5+d1,'NewScalar2':5+2*d1},\
               {'T3Scalar1Minus':d1,'T3Scalar1Plus':d1,'NewScalar2':d2}
    if model in ('A_-2','D_-2','E_-2'):
        return n,{'NewScalar1':5,'T3Scalar2Minus':5+2*d1,'T3Scalar2Plus':5+2*d1+d2},\
               {'NewScalar1':d1,'T3Scalar2Minus':d2,'T3Scalar2Plus':d2}
    return n,{'NewScalar1':5,'NewScalar2':5+2*d1},\
           {'NewScalar1':d1,'NewScalar2':d2}


def _row(key, real, imag=0):
    return {'coordinate':list(key),'coefficient_real':real,
            'coefficient_imag':imag,'orbit_size':len(orbit(key))}


def write_source(where: Path,model: str):
    n,offsets,dims=_layout(model)
    key=(1,1,2,3,5,4+n)
    assert canonical(key)==key
    r1=_row(key,1e-13,2e-14)
    r2=_row(key,-2e-14,1e-14)
    total=_row(key,8e-14,3e-14)
    sources=[{'matchete_source_index':11,'native_chirality':'PL',
              'independent_symmetric_coordinates':[r1]},
             {'matchete_source_index':12,'native_chirality':'PL',
              'independent_symmetric_coordinates':[r2]}]
    o={'format':SOURCE_FORMAT,'status':SOURCE_STATUS,'model':model,
       'mu_F_gev':100000.0,'units':'GeV^-1',
       'field_content':{'heavy_fermions_integrated_out':'all at mu_F',
                        'heavy_scalar_multiplets':'both retained in EFT1'},
       'index_order':['lepton_su2_i','lepton_flavour_p','lepton_su2_j',
                      'lepton_flavour_q','real_scalar_a','real_scalar_b'],
       'index_base':1,'real_scalar_dimension':4+n,
       'scalar_field_offsets_one_based':offsets,
       'scalar_field_complex_dimensions':dims,
       'symmetric_tensor_orbit_rule':
          'All distinct i,p<->j,q and a<->b orbit members equal representative. No extra orbit factor.',
       'formal_polynomial':'sum C[i,p,j,q,a,b] psi[i,p] psi[j,q] phi[a] phi[b]',
       'spinor_normalization':{'certified':False,'gamma_cc_to_master_weyl_factor':None,
                               'gamma_cc_to_master_weyl_sign':None},
       'completeness':{'F_tree_LLSS_included':True,
                       'distinct_PL_sources_retained':True,
                       'PR_Hermitian_conjugate_not_double_counted':True,
                       'one_loop_F_hard_matching_included':False,
                       'EFT1_running_performed':False,
                       'mu_S_scalar_matching_performed':False},
       'sources':sources,'summed_independent_symmetric_coordinates':[total]}
    where.parent.mkdir(parents=True,exist_ok=True)
    where.write_text(json.dumps(o,indent=2),encoding='utf-8')
    return where


def make_release(root:Path):
    release=root/'release'
    models=[]
    for model in ALL_MODELS:
        path=f'{slug(model)}/boundary.json'
        write_source(release/path,model)
        models.append({'model':model,'path':path,'status':SOURCE_STATUS})
    (release/'summary.json').write_text(json.dumps({
         'format':SOURCE_FORMAT,'status':RELEASE_STATUS,
         'n_models_verified':16,'absolute_spinor_normalization_certified':False,
         'physical_eft1_workers_must_reject_this_as_certified_input':True,
         'models':models},indent=2),encoding='utf-8')
    return release


def make_endpoint(root:Path, source, *, factor=-2.0):
    path=root/f'{slug(source.model)}.json'
    path.parent.mkdir(parents=True,exist_ok=True)
    fixture=make_fixture(source.model)
    n=len(source.labels)
    U=np.eye(n)
    # Deliberately make field basis nontrivial to test orientation of transform.
    U[0,0]=0; U[1,1]=0; U[0,1]=1; U[1,0]=1
    muF=np.einsum('fgijab,ak,bl->fgijkl',factor*source.raw,U,U,optimize=True)
    # Fixture muS tensors are independent synthetic algebra and cannot be
    # construed as a source-derived renormalized EFT1 trajectory.
    np.savez_compressed(path.with_suffix('.npz'),source_to_matching_real_U=U,
        C_LLSS_muF_matching=muF,C_LLSS_muS_hbar0=fixture.c_llss,
        lambda_HHSS_muS=fixture.lambda_hhss,
        lambda_HHSS_channels_muS=fixture.quartic_channels,
        mass_squared_gev2_muS=fixture.mass_squared,
        C5_EFT1_muS=fixture.c5_eft1)
    evidence=path.parent/(slug(source.model)+'_signed_amplitude.txt')
    evidence.write_text('SYNTHETIC TEST ONLY. Not an amplitude verification.\n',encoding='utf-8')
    meta={'format':ENDPOINT_FORMAT,'model':source.model,
        'source_boundary_sha256':source.hash,
        'source_scalar_real_labels':list(source.labels),
        'tensor_npz':path.with_suffix('.npz').name,
        'mu_f_gev':source.mu_f_gev,'mu_s_gev':fixture.mu_s_gev,
        'reference_scale_gev':fixture.reference_scale_gev,
        'renormalization_scheme':'MSbar','tensor_convention':CONVENTION,
        'scalar_mass_origin':'unbroken_electroweak',
        'quartic_channel_names':list(fixture.quartic_channel_names),
        'eft1_evolution':{'status':'renormalized_full_EFT1_to_muS',
             'same_muS_wilson_quartic_mass_C5':True,
             'no_separate_scalar_threshold':True,
             'C_LLSS_used_in_bubble_order':'hbar0',
             'C5_EFT1_includes_running_once':True},
        'normalization_certificate':{
             'status':'independent_signed_amplitude_verified',
             'source_boundary_sha256':source.hash,
             'master_operator_convention':CONVENTION,
             'gamma_cc_to_master_weyl_signed_factor':factor,
             'evidence_file':evidence.name,'evidence_sha256':digest(evidence)}}
    path.write_text(json.dumps(meta,indent=2),encoding='utf-8')
    return path


class FermionBridgeTests(unittest.TestCase):
    def test_all_16_source_real_layouts_and_complete_orbits(self):
        with tempfile.TemporaryDirectory() as td:
            release=make_release(Path(td))
            paths=load_release_manifest(release)
            self.assertEqual(set(paths),set(ALL_MODELS))
            for model in ALL_MODELS:
                with self.subTest(model=model):
                    s=load_source(paths[model],expected_model=model)
                    self.assertEqual(s.by_source.shape[0],2)
                    self.assertEqual(s.raw.shape[-1],8 if model[0]!='E' else 12)
                    self.assertEqual(len(s.labels),s.raw.shape[-1])
                    self.assertFalse(s.normalization_certified)
                    self.assertTrue(np.allclose(s.raw,np.sum(s.by_source,axis=0),rtol=1e-12,atol=1e-26))
                    self.assertEqual(np.count_nonzero(s.raw),4)
                    self.assertNotEqual(s.raw.imag.max(),0)

    def test_named_source_export_no_unwanted_spinor_factor(self):
        with tempfile.TemporaryDirectory() as td:
            root=Path(td)
            s=load_source(write_source(root/'B_m1/boundary.json','B_-1'))
            target=root/'result/B_m1.npz'
            write_provisional_source(target,s)
            obj=json.loads(target.with_suffix('.json').read_text())
            self.assertFalse(obj['gamma_cc_weyl_normalization_certified'])
            self.assertFalse(obj['full_physical_C5_computed'])
            with np.load(target) as z:
                self.assertTrue(np.array_equal(z['C_LLSS_formal_muF'],s.raw))
                self.assertTrue(np.array_equal(z['C_LLSS_by_PL_source_formal_muF'],s.by_source))
                self.assertEqual(z['matchete_PL_source_indices'].tolist(),[11,12])

    def test_rejects_source_orbit_multiplicity_or_double_count(self):
        with tempfile.TemporaryDirectory() as td:
            p=write_source(Path(td)/'A_m4/boundary.json','A_-4')
            o=json.loads(p.read_text())
            o['summed_independent_symmetric_coordinates'][0]['orbit_size']=2
            p.write_text(json.dumps(o))
            with self.assertRaisesRegex(ValueError,'orbit size'):
                load_source(p)
            p=write_source(p,'A_-4')
            o=json.loads(p.read_text())
            o['sources'][0]['native_chirality']='PR'
            p.write_text(json.dumps(o))
            with self.assertRaisesRegex(ValueError,'PR or unknown'):
                load_source(p)
            p=write_source(p,'A_-4')
            o=json.loads(p.read_text())
            o['summed_independent_symmetric_coordinates'][0]['coefficient_real']*=2
            p.write_text(json.dumps(o))
            with self.assertRaisesRegex(ValueError,'source sum mismatches'):
                load_source(p)

    def test_zero_at_benchmark_PL_terms_are_kept(self):
        with tempfile.TemporaryDirectory() as td:
            p=write_source(Path(td)/'C_m1/boundary.json','C_-1')
            o=json.loads(p.read_text())
            o['sources'].append({'matchete_source_index':13,'native_chirality':'PL',
                                 'zero_at_this_benchmark':True,
                                 'independent_symmetric_coordinates':[]})
            p.write_text(json.dumps(o))
            s=load_source(p)
            self.assertEqual(s.by_source.shape[0],3)
            self.assertTrue(np.array_equal(s.by_source[2],np.zeros_like(s.by_source[2])))
            self.assertEqual(s.source_indices,(11,12,13))
            self.assertTrue(np.allclose(np.sum(s.by_source,axis=0),s.raw,rtol=1e-12,atol=1e-26))

    def test_endpoints_fail_closed_for_undefined_gamma_cc(self):
        with tempfile.TemporaryDirectory() as td:
            root=Path(td)
            s=load_source(write_source(root/'src/B_m1/boundary.json','B_-1'))
            p=make_endpoint(root/'end',s)
            o=json.loads(p.read_text())
            o.pop('normalization_certificate')
            p.write_text(json.dumps(o))
            with self.assertRaisesRegex(ValueError,'GammaCC/Weyl normalization certificate missing'):
                load_eft1_endpoint(p,s)
            p=make_endpoint(root/'end',s)
            o=json.loads(p.read_text())
            o['normalization_certificate']['gamma_cc_to_master_weyl_signed_factor']=None
            p.write_text(json.dumps(o))
            with self.assertRaisesRegex(ValueError,'explicit signed'):
                load_eft1_endpoint(p,s)

    def test_endpoint_basis_rotation_and_gamma_factor_exact(self):
        with tempfile.TemporaryDirectory() as td:
            root=Path(td)
            s=load_source(write_source(root/'src/B_m1/boundary.json','B_-1'))
            p=make_endpoint(root/'end',s)
            inp=load_eft1_endpoint(p,s)
            saved=save_input(root/'normal/B_m1.json',inp)
            again=load_input(saved)
            out=compute(again)
            self.assertLess(out['weinberg_projection_residual'],1e-10)
            self.assertTrue(np.iscomplexobj(again.c_llss))
            self.assertFalse(out['physical_full_one_loop_certified'])
            with np.load(p.with_suffix('.npz'),allow_pickle=False) as z:
                arrays={k:np.array(z[k]) for k in z.files}
            arrays['C_LLSS_muF_matching']*=1.001
            np.savez_compressed(p.with_suffix('.npz'),**arrays)
            with self.assertRaisesRegex(ValueError,'muF tensor differs'):
                load_eft1_endpoint(p,s)

    def test_bridge_cli_source_only_and_missing_endpoints(self):
        with tempfile.TemporaryDirectory() as td:
            root=Path(td)
            release=make_release(root)
            # CLI writes outside project's output only if already located there;
            # exercise the CLI in a temporary subdirectory of the allowed output.
            from scripts.BridgeT3JointScalarEFT1 import SAFE_OUTPUT
            parent=SAFE_OUTPUT/'_bridge_unit_test_sandbox'
            import shutil
            shutil.rmtree(parent,ignore_errors=True)
            try:
                rc=bridge_main(['--release-dir',str(release),
                                '--output-dir',str(parent/'audit')])
                self.assertEqual(rc,0)
                summary=json.loads((parent/'audit/summary.json').read_text())
                self.assertEqual(summary['source_muF_exported'],16)
                self.assertEqual(summary['conditional_matching_candidates'],0)
                self.assertFalse(summary['physical_C5_certified'])
                rc=bridge_main(['--release-dir',str(release),'--model','B_-1',
                                '--endpoint-dir',str(root/'missing'),
                                '--output-dir',str(parent/'blocked')])
                self.assertEqual(rc,2)
                self.assertFalse((parent/'blocked/candidate_C5').exists())
            finally:
                shutil.rmtree(parent,ignore_errors=True)


if __name__=='__main__':
    unittest.main()
