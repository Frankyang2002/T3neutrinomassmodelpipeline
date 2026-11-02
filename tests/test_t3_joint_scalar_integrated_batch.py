"""Pure algebra/contract tests; no physical T3 data or signed GammaCC assumed."""
from __future__ import annotations

from dataclasses import replace
import json
from pathlib import Path
import shutil
import tempfile
import unittest

import numpy as np

from Numerical.joint_scalar_matching.adapter import ALL_MODELS, CONVENTION
from Numerical.joint_scalar_matching.core import LOOP, weinberg_tensor
from Numerical.joint_scalar_matching.fermion_boundary_bridge import (
    SourceBoundary,SOURCE_STATUS,digest,load_source,slug)
from Numerical.joint_scalar_matching.large_batch import (
    BACKGROUND_FORMAT,KERNEL_FORMAT,KERNEL_STATUS,
    ScalarBackground,encode_matrix,parse_matrix,_beta_project,
    _response,evaluate_conditional,load_background,validate_pivot_kernel)
from scripts.RunT3JointScalarIntegratedBatch import SAFE_OUTPUT,main as batch_main
from test_t3_joint_scalar_fermion_boundary_bridge import make_release


def clean_synthetic_source(model,where):
    n=8 if model[0]!='E' else 12
    T=weinberg_tensor()
    flavor=np.array([[1+0.21j,.3-.12j,-.11+.04j],
                     [.3-.12j,.8+.04j,.2-.17j],
                     [-.11+.04j,.2-.17j,.6+.3j]],complex)*1e-11
    C=np.zeros((3,3,2,2,n,n),complex)
    Q=np.zeros((4,4,n,n),float)
    chans=np.zeros((6,4,4,n,n),float)
    dim1=2*{'A':1,'B':2,'C':2,'D':3,'E':3}[model[0]]
    pairs=[(0,0),(1,1),(0,dim1),(1,dim1+1),(dim1,dim1),(dim1+1,dim1+1)]
    for t,(a,b) in enumerate(pairs):
        weak=(0,0) if t<2 else (0,1) if t<4 else (1,1)
        w=1 if t%2==0 else 1j
        cweak=np.zeros((2,2),complex)
        cweak[weak]=w;cweak[weak[::-1]]=w
        coeff=T[weak].real if t%2==0 else T[weak].imag
        val=coeff/(2 if a!=b else 1)
        C[...,a,b]+=np.einsum('fg,ij->fgij',flavor,cweak)
        if a!=b:C[...,b,a]+=np.einsum('fg,ij->fgij',flavor,cweak)
        Q[...,a,b]+=val
        if a!=b: Q[...,b,a]+=val
        chans[t,...,a,b]+=val
        if a!=b:chans[t,...,b,a]+=val
    path=where/'source.json'
    path.parent.mkdir(parents=True,exist_ok=True)
    path.write_text('synthetic source',encoding='utf-8')
    return SourceBoundary(model,path,digest(path),100000.,C,np.stack((0.6*C,0.4*C)),
                          (101,102),tuple('S1' if i<dim1 else 'S2' for i in range(n)),
                          False,SOURCE_STATUS),Q,chans


def make_background(where,source,Q,channels,*,U=None,mu_s=1000.):
    n=len(source.labels)
    U=np.eye(n) if U is None else U
    Qm=np.einsum('ABab,ak,bl->ABkl',Q,U,U,optimize=True)
    chm=np.einsum('tABab,ak,bl->tABkl',channels,U,U,optimize=True)
    where.mkdir(parents=True,exist_ok=True)
    doc={'format':BACKGROUND_FORMAT,'model':source.model,
         'source_boundary_sha256':source.hash,
         'source_scalar_real_labels':list(source.labels),
         'renormalization_scheme':'MSbar','tensor_convention':CONVENTION,
         'mu_f_gev':source.mu_f_gev,'mu_s_gev':mu_s,
         'reference_scale_gev':mu_s,
         'no_separate_scalar_threshold':True,
         'same_muS_quartic_and_mass':True,
         'scalar_mass_origin':'unbroken_electroweak',
         'tensor_npz':'background.npz',
         'quartic_channel_names':[f'quartic_{i}' for i in range(len(channels))],
         'provenance':{'source':'synthetic_unit_test_NOT_physical',
            'quartic_potential_normalization':'V=lambda_abcd*phi_a*phi_b*phi_c*phi_d/4!'}}
    path=where/'background.json'
    path.write_text(json.dumps(doc),encoding='utf-8')
    np.savez_compressed(where/'background.npz',source_to_matching_real_U=U,
        mass_squared_gev2_muS=np.eye(n)*(1200.**2),
        lambda_HHSS_muS=Qm,lambda_HHSS_channels_muS=chm)
    return path


def write_pivot_kernel(where,source):
    n=len(source.labels)
    rows=[];summed={}
    for t,idx in enumerate(source.source_indices):
        channels=[]
        for a in range(n):
            for b in range(a,n):
                mult=1 if a==b else 2
                value=4*mult*source.by_source[t,:,:,0,0,a,b]
                if np.max(abs(value))==0:continue
                channels.append({'heavy_scalar_pair_one_based':[a+5,b+5],
                                 'multiplicity_in_ordered_ab_sum':mult,
                                 'indeterminate_quartic':f'lambda_[H0R,H0R,{a+5},{b+5}]',
                                 'conditional_beta_C5_per_unit_quartic':encode_matrix(value)})
                key=(a,b)
                summed[key]=summed.get(key,np.zeros((3,3),complex))+value
        rows.append({'matchete_source_index':idx,'kernel_channels':channels})
    both=[]
    for (a,b),value in sorted(summed.items()):
        both.append({'heavy_scalar_pair_one_based':[a+5,b+5],
                    'conditional_beta_C5_per_unit_quartic':encode_matrix(value)})
    doc={'format':KERNEL_FORMAT,'status':KERNEL_STATUS,'model':source.model,
         'source_boundary_sha256':source.hash,
         'real_scalar_dimension':4+n,
         'H0R_real_scalar_index_one_based':3,
         'physical_normalization':{'certified':False,'N_GammaCC_to_Weyl':None},
         'quartic_component_couplings':'UNSUPPLIED - no numerical C5 is computed',
         'per_source':rows,'summed_kernel_channels':both}
    where.parent.mkdir(parents=True,exist_ok=True)
    where.write_text(json.dumps(doc),encoding='utf-8')
    return where


class IntegratedBatchTests(unittest.TestCase):
    def test_16_source_pivots_agree_with_independent_formula(self):
        with tempfile.TemporaryDirectory() as td:
            root=Path(td)
            for model in ALL_MODELS:
                with self.subTest(model=model):
                    source,_,_=clean_synthetic_source(model,root/model)
                    p=write_pivot_kernel(root/(model+'_pivot.json'),source)
                    audit=validate_pivot_kernel(p,source)
                    self.assertEqual(audit['source_count'],2)
                    self.assertTrue(audit['distinct_scalar_pairs']>0)
                    self.assertLess(audit['source_vs_kernel_max_abs_error'],1e-26)

    def test_source_pivot_corrupted_multiplicity_and_signed_amplitude(self):
        with tempfile.TemporaryDirectory() as td:
            root=Path(td)
            source,_,_=clean_synthetic_source('A_-4',root)
            p=write_pivot_kernel(root/'pivot.json',source)
            obj=json.loads(p.read_text())
            obj['per_source'][0]['kernel_channels'][0]['multiplicity_in_ordered_ab_sum']=7
            p.write_text(json.dumps(obj))
            with self.assertRaisesRegex(ValueError,'multiplicity'):
                validate_pivot_kernel(p,source)
            p=write_pivot_kernel(p,source)
            obj=json.loads(p.read_text())
            obj['per_source'][0]['kernel_channels'][0]['conditional_beta_C5_per_unit_quartic']['real'][0][0] *= -1
            p.write_text(json.dumps(obj))
            with self.assertRaisesRegex(ValueError,'disagrees'):
                validate_pivot_kernel(p,source)

    def test_exact_inter_source_pair_cancellation_is_retained_and_omitted(self):
        with tempfile.TemporaryDirectory() as td:
            root=Path(td)
            src,_,_=clean_synthetic_source('E_-2',root)
            # Nonzero individual PL contributions cancel in the total.  This
            # occurs in the actual published E-class pivot exports.
            source=replace(src,raw=np.zeros_like(src.raw),
                           by_source=np.stack((src.raw,-src.raw)))
            p=write_pivot_kernel(root/'pivot.json',source)
            o=json.loads(p.read_text())
            o['summed_kernel_channels']=[]
            p.write_text(json.dumps(o))
            audit=validate_pivot_kernel(p,source)
            self.assertEqual(audit['distinct_scalar_pairs'],0)
            self.assertGreater(audit['cancelled_source_pairs_omitted_by_kernel'],0)
            # If the cancellation is even slightly broken, the missing
            # summed channel must instead block the audit.
            bad=source.by_source.copy()
            bad[1] *= 0.999
            damaged=replace(source,raw=bad.sum(axis=0),by_source=bad)
            p=write_pivot_kernel(p,damaged)
            o=json.loads(p.read_text())
            o['summed_kernel_channels']=[]
            p.write_text(json.dumps(o))
            with self.assertRaisesRegex(ValueError,'omitted an uncancelled'):
                validate_pivot_kernel(p,damaged)

    def test_16_fixed_order_log_and_bubble_are_ward_consistent(self):
        with tempfile.TemporaryDirectory() as td:
            root=Path(td)
            for model in ALL_MODELS:
                with self.subTest(model=model):
                    source,q,channels=clean_synthetic_source(model,root/model)
                    p=make_background(root/model/'bg',source,q,channels)
                    background=load_background(p,source)
                    result=evaluate_conditional(source,background)
                    self.assertLess(result['weinberg_beta_projection_residual'],1e-10)
                    self.assertLess(result['weinberg_bubble_projection_residual'],1e-10)
                    self.assertLess(result['relative_matching_scale_cancellation_residual'],1e-9)
                    self.assertIsNone(result['source_to_eft_weyl_normalization_factor'])
                    self.assertFalse(result['physical_C5_certified'])
                    total=parse_matrix(result['C5_LLSS_total_one_loop_per_N_gev_inverse'])
                    self.assertTrue(np.allclose(total,total.T,atol=1e-26,rtol=1e-12))

    def test_nontrivial_basis_rotation_matches_full_response(self):
        with tempfile.TemporaryDirectory() as td:
            root=Path(td)
            source,q,ch=clean_synthetic_source('B_-1',root)
            n=len(source.labels)
            U=np.eye(n)
            U[0,0]=0; U[1,1]=0; U[0,1]=1; U[1,0]=1
            a=load_background(make_background(root/'a',source,q,ch,U=np.eye(n)),source)
            b=load_background(make_background(root/'b',source,q,ch,U=U),source)
            ra=evaluate_conditional(source,a);rb=evaluate_conditional(source,b)
            for key in ('C5_joint_scalar_finite_bubble_per_N_gev_inverse',
                        'C5_EFT1_induced_leading_log_per_N_gev_inverse',
                        'C5_LLSS_total_one_loop_per_N_gev_inverse'):
                self.assertTrue(np.allclose(parse_matrix(ra[key]),parse_matrix(rb[key]),
                                            rtol=1e-10,atol=1e-27))

    def test_log_cancels_scale_of_finite_bubble(self):
        with tempfile.TemporaryDirectory() as td:
            source,q,ch=clean_synthetic_source('E_-2',Path(td))
            C=source.raw
            m2=1200.**2*np.eye(len(source.labels))
            base=_response(C,q,m2,100000.,1000.)
            shift=_response(C,q,m2,100000.,1180.)
            self.assertTrue(np.allclose(base['sum'],shift['sum'],atol=1e-25,rtol=1e-10))
            self.assertTrue(np.allclose(shift['running']-base['running'],
               -(shift['bubble']-base['bubble']),atol=1e-25,rtol=1e-9))
            self.assertFalse(np.allclose(base['running'],base['bubble'],atol=1e-25,rtol=1e-9))

    def test_background_fail_closed_and_quartic_complete(self):
        with tempfile.TemporaryDirectory() as td:
            root=Path(td)
            source,q,ch=clean_synthetic_source('D_-2',root)
            p=make_background(root/'bg',source,q,ch)
            obj=json.loads(p.read_text());obj['source_boundary_sha256']='BAD';p.write_text(json.dumps(obj))
            with self.assertRaisesRegex(ValueError,'correct fermion'):
                load_background(p,source)
            obj['source_boundary_sha256']=source.hash
            obj['quartic_channel_names']=[];p.write_text(json.dumps(obj))
            with self.assertRaisesRegex(ValueError,'quartic channels'):
                load_background(p,source)
            obj['quartic_channel_names']=[f'quartic_{i}' for i in range(6)]
            p.write_text(json.dumps(obj))
            with np.load(root/'bg/background.npz') as z: ar={k:np.array(z[k]) for k in z.files}
            ar['lambda_HHSS_channels_muS'][0]*=1.02
            np.savez_compressed(root/'bg/background.npz',**ar)
            with self.assertRaisesRegex(ValueError,'do not sum'):
                load_background(p,source)

    def test_background_rejects_nonorthogonal_transform(self):
        with tempfile.TemporaryDirectory() as td:
            root=Path(td)
            source,q,ch=clean_synthetic_source('A_+0',root)
            n=len(source.labels)
            U=np.eye(n);U[0,1]=0.1
            p=make_background(root/'bg',source,q,ch,U=U)
            with self.assertRaisesRegex(ValueError,'orthogonal'):
                load_background(p,source)

    def test_background_missing_inputs_do_not_imply_zero_C5(self):
        with tempfile.TemporaryDirectory() as td:
            root=Path(td)
            source,q,ch=clean_synthetic_source('C_-1',root)
            p=make_background(root/'bg',source,q,ch)
            (root/'bg/background.npz').unlink()
            with self.assertRaises(FileNotFoundError):
                load_background(p,source)

    def test_real_repository_audit_cli_16_synthetic_layouts(self):
        with tempfile.TemporaryDirectory() as td:
            root=Path(td)
            release=make_release(root)
            kern=root/'kernels'
            for model in ALL_MODELS:
                source=load_source(release/slug(model)/'boundary.json')
                write_pivot_kernel(kern/slug(model)/'quartic_pivot_kernel.json',source)
            dest=SAFE_OUTPUT/'_integrated_batch_unit_test_sandbox'
            shutil.rmtree(dest,ignore_errors=True)
            try:
                rc=batch_main(['--release-dir',str(release),'--kernel-dir',str(kern),
                               '--output-dir',str(dest)])
                self.assertEqual(rc,0)
                summary=json.loads((dest/'summary.json').read_text())
                self.assertEqual(summary['models_pivot_verified'],16)
                self.assertEqual(summary['conditional_response_models'],0)
                self.assertFalse(summary['physical_C5_certified'])
                self.assertEqual(len(list((dest/'scalar_background_contracts').glob('*.json'))),16)
                rc2=batch_main(['--release-dir',str(release),'--kernel-dir',str(kern),
                               '--output-dir',str(dest/'second'),
                               '--background-dir',str(root/'missing'),
                               '--model','B_-1'])
                self.assertEqual(rc2,2)
                summary2=json.loads((dest/'second/summary.json').read_text())
                self.assertEqual(summary2['blocked_models'],1)
            finally:
                shutil.rmtree(dest,ignore_errors=True)


if __name__=='__main__':unittest.main()
