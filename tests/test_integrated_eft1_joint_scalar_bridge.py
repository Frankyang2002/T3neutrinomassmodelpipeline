"""Synthetic-only contract and algebra tests; no physical T3 couplings assumed."""
from __future__ import annotations

import copy
import json
import math
from pathlib import Path
import sys
import tempfile
import types
import unittest
from unittest.mock import patch

import numpy as np

from RGE.running.intermediate.ProvisionalSourceKernel import SOURCE_FORMAT, _orbit, MODELS, tag, SourceKernelError
from RGE.running.intermediate.SourceBackedEFT1Batch import (
    COUPLING_FORMAT, _validate_couplings,
)
from RGE.running.intermediate.GeneralEFT1C5Mixing import (
    LOOP, load_boundary,
)
from RGE.running.intermediate.FullSU2SourceResponse import weinberg_template
from RGE.running.intermediate.IntegratedEFT1JointScalarBridge import (
    FORMAT, contracts, complex_matrix, _same, _master_source_heavy,
    _project, _bubble_c5, _backrotate, evaluate,
)
from scripts import RunIntegratedEFT1JointScalarBridge as runner

# In the installed repository these are the real pre-existing joint matcher
# modules. The local isolated package tests permit a minimal algebraic stand-in
# ONLY when that upstream workstream has not been mounted in this test session.
try:
    from Numerical.joint_scalar_matching.core import scalar_bubble_matching as real_bubble
    from Numerical.joint_scalar_matching.core import project_weinberg as real_project
    from Numerical.joint_scalar_matching.fermion_boundary_bridge import load_source
    from Numerical.joint_scalar_matching.large_batch import load_background
    HAS_REAL_JOINT=True
except ImportError:
    HAS_REAL_JOINT=False


def enc(z):
    return {'real':np.asarray(z).real.tolist(),'imag':np.asarray(z).imag.tolist()}


def _setup_joint_stubs():
    """Exactly reproduce the bubble and template sign for local algebra tests."""
    if HAS_REAL_JOINT:return
    num=types.ModuleType('Numerical');num.__path__=[]
    mod=types.ModuleType('Numerical.joint_scalar_matching');mod.__path__=[]
    core=types.ModuleType('Numerical.joint_scalar_matching.core')
    fm=types.ModuleType('Numerical.joint_scalar_matching.fermion_boundary_bridge')
    lb=types.ModuleType('Numerical.joint_scalar_matching.large_batch')
    def project_joint(amp,strict=False):
        t=2*weinberg_template()
        z=np.einsum('fgijAB,ijAB->fg',amp,t.conj())/np.vdot(t,t)
        res=np.linalg.norm(amp-np.einsum('fg,ijAB->fgijAB',z,t))/max(np.linalg.norm(amp),1e-30)
        return z,float(res)
    def bub(C,q,m2,mu_s,reference_scale=None):
        from math import log,pi
        n=C.shape[-1]
        amp=np.zeros((3,3,2,2,4,4),complex)
        pairs={}
        for a in range(n):
            for b in range(a,n):
                d=(1 if a==b else 2)
                x,y=m2[a,a],m2[b,b]
                B= -log(x/(mu_s*mu_s)) if abs(x-y)<1e-12 else \
                   1-(x*log(x/(mu_s*mu_s))-y*log(y/(mu_s*mu_s)))/(x-y)
                v=-d/(2*16*pi*pi)*np.einsum('fgij,AB->fgijAB',C[...,a,b],q[...,a,b])*B
                pairs[f'e{a:02d}_e{b:02d}']=v
                amp+=v
        return {'amplitude':amp,'pair_amplitudes':pairs,
                'eigenvalues_sq_gev2':np.diag(m2),'mass_rotation_U':np.eye(n)}
    core.project_weinberg=project_joint
    core.scalar_bubble_matching=bub
    fm.load_source=lambda *a,**k:None
    lb.load_background=lambda *a,**k:None
    for name,module in [('Numerical',num),('Numerical.joint_scalar_matching',mod),
       ('Numerical.joint_scalar_matching.core',core),
       ('Numerical.joint_scalar_matching.fermion_boundary_bridge',fm),
       ('Numerical.joint_scalar_matching.large_batch',lb)]:
        sys.modules[name]=module


def synthetic_boundary(model='A_-4', *, asymmetric=False):
    """Build source C using 6 heavy diagonal scalar channels; not T3 physics."""
    rows=[]
    heavy_entries=[5,6,7,8,9,10]
    for k,(i,j) in enumerate([(1,1),(1,2),(2,2)]):
        for mode in range(2):
            a=heavy_entries[2*k+mode]
            key=(i,1,j,1,a,a)
            rows.append({'coordinate':list(key),'coefficient_real':(1. if mode==0 else 0.),
                         'coefficient_imag':(0. if mode==0 else 1.),
                         'orbit_size':len(_orbit(key))})
    if asymmetric:rows[0]['coefficient_real']=0.9
    return {'format':SOURCE_FORMAT,'model':model,
      'status':'NumericalLLSSBySourceVerified_SpinorNormalizationOpen',
      'mu_F_gev':1e5,'units':'GeV^-1','index_base':1,
      'index_order':['lepton_su2_i','lepton_flavour_p','lepton_su2_j','lepton_flavour_q',
                     'real_scalar_a','real_scalar_b'],
      'real_scalar_dimension':12,
      'scalar_field_offsets_one_based':{'NewScalar1':5,'NewScalar2':7},
      'scalar_field_complex_dimensions':{'NewScalar1':1,'NewScalar2':3},
      'spinor_normalization':{'certified':False,
          'gamma_cc_to_master_weyl_factor':None,'gamma_cc_to_master_weyl_sign':None,
          'do_not_use_as_physical_master_wilson_boundary':True},
      'completeness':{'F_tree_LLSS_included':True,
        'distinct_PL_sources_retained':True,'PR_Hermitian_conjugate_not_double_counted':True,
        'complete_physical_leading_one_loop_neutrino_mass_certified':False},
      'sources':[{'matchete_source_index':17,'native_chirality':'PL',
           'zero_at_this_benchmark':False,'independent_symmetric_coordinates':rows}],
      'summed_independent_symmetric_coordinates':copy.deepcopy(rows)}


def synthetic_couplings(boundary):
    """A *mathematical* SU2 Weinberg tensor constructed from source and q.

    Separate real/imag heavy diagonal channels encode complex Higgs polynomials.
    Avoids assuming any physical scalar quartic potential in T3.
    """
    T=weinberg_template()
    entries=[]
    for k,(i,j) in enumerate([(0,0),(0,1),(1,1)]):
        for mode in range(2):
            a=5+2*k+mode
            x=(T[i,j].real if mode==0 else T[i,j].imag)
            for H in range(4):
                for J in range(H,4):
                    if abs(x[H,J])>1e-14:
                        entries.append({'indices':[H+1,J+1,a,a], 'value':float(x[H,J])})
    return {'format':COUPLING_FORMAT,'model':boundary['model'],
        'source_boundary_sha256':'not-set','mu_F_gev':1e5,'mu_S_gev':1000.,
        'n_fermions':6,'n_real_scalars':12,'renormalization_scheme':'MSbar',
        'master_conventions':{
            'master_wilson':'C_ijab psi_i psi_j phi_a phi_b + hc',
            'weinberg':'(1/2) C5_ij (L_i.H)(L_j.H) + hc',
            'quartic':'V4=lambda_abcd phi_a phi_b phi_c phi_d/4!'},
        'coverage_assertion':{'authority':'diagnostic_only',
            'all_active_weyl_fermion_components_enumerated':True,
            'all_yukawa_components_enumerated':True,
            'all_quartic_channels_enumerated':True,
            'all_gauge_generators_enumerated_or_explicitly_toy_zero':True,
            'toy_gauge_couplings_set_to_zero':True,
            'fixed_order_quartic_identical_muF_muS':True},
        'yukawa_components':[],'gauge_sectors':[],
        'quartic_channels':[{'label':'mathematical_SU2_invariant', 'components':entries}]}


class AlgebraTests(unittest.TestCase):
    def test_joint_and_master_templates_have_exact_factor_two(self):
        _setup_joint_stubs()
        from Numerical.joint_scalar_matching.core import project_weinberg
        T=weinberg_template()
        amp=np.zeros((3,3,2,2,4,4),complex)
        amp[0,0]=T*(3+4j)
        master,_=_project(amp,'test',strict=True)
        joint, residual=project_weinberg(amp,strict=False)
        np.testing.assert_allclose(master,2*joint,atol=1e-14,rtol=1e-14)
        self.assertLess(residual,1e-13)

    def test_neutrino_pivot_cannot_fake_ward(self):
        amp=np.zeros((3,3,2,2,4,4),complex)
        amp[0,0,0,0,2,2]=1.
        with self.assertRaisesRegex(SourceKernelError,'SU2 Ward'):
            _project(amp,'bad',strict=True)

    def test_wrong_basis_fails_close(self):
        a=np.eye(4);b=a.copy();b[0,0]*=0.9
        with self.assertRaisesRegex(SourceKernelError,'differ'):_same(a,b,'quartic')

    def test_unknown_factor_and_independent_hard_matching_in_contract(self):
        b=synthetic_boundary()
        c=contracts(b,'abc')
        self.assertIsNone(c['physical_source_spinor_factor'])
        self.assertFalse(c['complete_one_loop_UV_to_SMEFT_certified'])
        self.assertEqual(c['required_joint_scalar_background']['numerical_array_shapes']['lambda_HHSS_muS'],
                         [4,4,8,8])

    def test_source_heavy_tensor_is_not_doubled_by_symmetric_orbit(self):
        b=synthetic_boundary();q=synthetic_couplings(b)
        p,_=_validate_couplings(q,b,'not-set');ctx=load_boundary(p,allow_synthetic=True)
        c=_master_source_heavy(ctx)
        self.assertEqual(c.shape,(3,3,2,2,8,8))
        self.assertEqual(c[0,0,0,1,2,2],1)
        self.assertEqual(c[0,0,1,0,2,2],1)
        self.assertEqual(c[0,0,0,0,0,0],1)

    def test_rotation_source_matching_identity(self):
        U=np.linalg.qr(np.random.default_rng(17).normal(size=(8,8)))[0]
        q=np.random.default_rng(29).normal(size=(4,4,8,8))
        q=(q+q.swapaxes(0,1)+q.swapaxes(2,3)+q.swapaxes(0,1).swapaxes(2,3))/4
        match=np.einsum('ABab,ai,bj->ABij',q,U,U)
        np.testing.assert_allclose(_backrotate(match,U),q,atol=1e-14)


class IntegratedTests(unittest.TestCase):
    def setUp(self):
        _setup_joint_stubs()

    def fixture(self, folder, *, qscale=1, wrong_label=False):
        folder=Path(folder)
        b=synthetic_boundary()
        boundary_path=folder/'boundary.json'
        boundary_path.write_text(json.dumps(b))
        from hashlib import sha256
        sha=sha256(boundary_path.read_bytes()).hexdigest()
        q=synthetic_couplings(b)
        q['source_boundary_sha256']=sha
        qpath=folder/'couplings.json';qpath.write_text(json.dumps(q))
        payload,_=_validate_couplings(q,b,sha)
        ctx=load_boundary(payload,allow_synthetic=True)
        C=_master_source_heavy(ctx)
        # SourceBoundary emulation used only when actual Chat-1 release is absent.
        source=types.SimpleNamespace(model=b['model'],hash=sha,mu_f_gev=1e5,
             raw=C,by_source=C[None,...],source_indices=(17,),
             labels=tuple(f'R{i}' for i in range(8)),normalization_certified=False)
        qt=ctx.quartic[:4,:4,4:,4:]
        if qscale!=1:qt*=qscale
        bg=types.SimpleNamespace(model=b['model'],mu_f=1e5,mu_s=1000.,
             reference_scale=1000., mass_origin='unbroken_electroweak',
             U=np.eye(8),quartic=qt,
             channels=np.array([qt]),names=('WRONG' if wrong_label else 'mathematical_SU2_invariant',),
             mass2=np.eye(8)*1200**2)
        bgpath=folder/'A_m4.json';bgpath.write_text(json.dumps({'tensor_npz':'A_m4.npz'}))
        np.savez(folder/'A_m4.npz', fake=np.eye(8))
        return boundary_path,qpath,bgpath,source,bg

    def test_full_seven_term_plus_joint_threshold_nonzero_and_log_cancel(self):
        with tempfile.TemporaryDirectory() as td:
            b,q,s,source,bg=self.fixture(td)
            with patch('Numerical.joint_scalar_matching.fermion_boundary_bridge.load_source',return_value=source),\
                 patch('Numerical.joint_scalar_matching.large_batch.load_background',return_value=bg):
                result=evaluate(b,q,s,trajectory_points=9)
            self.assertEqual(result['format'],FORMAT)
            self.assertFalse(result['physical_C5_computed'])
            self.assertIsNone(result['N_GammaCC_to_master_Weyl'])
            self.assertEqual(result['scalar_threshold_count'],1)
            self.assertLess(result['full_SU2_quartic_residual'],1e-12)
            self.assertLess(result['full_SU2_joint_bubble_residual'],1e-12)
            self.assertLess(result['quartic_muS_log_cancellation_relative_residual'],1e-10)
            beta=complex_matrix(result['C5_quartic_beta_16pi2_per_N'],'beta')
            self.assertAlmostEqual(beta[0,0].real,1.,places=11)
            self.assertGreater(abs(complex_matrix(
               result['C5_LLSS_induced_increment_to_SMEFT_per_N_gev_inv'],'sum')[0,0]),0)
            self.assertEqual(len(result['per_Matchete_PL_source']),1)
            self.assertEqual(len(result['per_independent_quartic_channel']),1)

    def test_wrong_background_quartic_blocks_without_output(self):
        with tempfile.TemporaryDirectory() as td:
            b,q,s,source,bg=self.fixture(td,qscale=1.1)
            with patch('Numerical.joint_scalar_matching.fermion_boundary_bridge.load_source',return_value=source),\
                 patch('Numerical.joint_scalar_matching.large_batch.load_background',return_value=bg):
                with self.assertRaisesRegex(SourceKernelError,'background differ'):
                    evaluate(b,q,s)

    def test_wrong_quartic_channel_fails(self):
        with tempfile.TemporaryDirectory() as td:
            b,q,s,source,bg=self.fixture(td,wrong_label=True)
            with patch('Numerical.joint_scalar_matching.fermion_boundary_bridge.load_source',return_value=source),\
                 patch('Numerical.joint_scalar_matching.large_batch.load_background',return_value=bg):
                with self.assertRaisesRegex(SourceKernelError,'quartic channel'):
                    evaluate(b,q,s)

    def test_missing_fixed_coupling_evidence_blocks(self):
        with tempfile.TemporaryDirectory() as td:
            b,q,s,source,bg=self.fixture(td)
            qdata=json.loads(q.read_text())
            del qdata['coverage_assertion']['fixed_order_quartic_identical_muF_muS']
            q.write_text(json.dumps(qdata))
            with patch('Numerical.joint_scalar_matching.fermion_boundary_bridge.load_source',return_value=source),\
                 patch('Numerical.joint_scalar_matching.large_batch.load_background',return_value=bg):
                with self.assertRaisesRegex(SourceKernelError,'frozen quartic'):
                    evaluate(b,q,s)

    def test_cli_preflight_all_16_no_overwrites(self):
        with tempfile.TemporaryDirectory() as td:
            root=Path(td)
            source_root=root/'release'
            rows=[]
            for model in MODELS:
                b=synthetic_boundary(model)
                path=source_root/tag(model)/'boundary.json';path.parent.mkdir(parents=True,exist_ok=True)
                path.write_text(json.dumps(b))
                rows.append({'model':model,'path':f'{tag(model)}/boundary.json',
                             'mu_F_gev':1e5,'PL_sources':1})
            (source_root/'summary.json').write_text(json.dumps({
                'format':'t3-lo-fermion-boundary-source-resolved-provisional-v1',
                'n_models_verified':16,'absolute_spinor_normalization_certified':False,
                'physical_eft1_workers_must_reject_this_as_certified_input':True,
                'models':rows}))
            o=root/'output/general_observable/lo_eft1_mixing/integration'
            original=runner.ROOT
            try:
                runner.ROOT=root
                args=['--release-root',str(source_root),'--output-root',str(o)]
                self.assertEqual(runner.main(args+['--check-only']),0)
                self.assertFalse(o.exists())
                self.assertEqual(runner.main(args),0)
                summary=json.loads((o/'summary.json').read_text())
                self.assertEqual(summary['model_count'],16)
                self.assertEqual(summary['models_evaluated_conditionally'],0)
                self.assertFalse(summary['physical_C5_computed'])
                with self.assertRaises(FileExistsError):runner.main(args)
                with self.assertRaises(SystemExit):runner.main(args+['--background-root',str(root)])
            finally:runner.ROOT=original

if __name__=='__main__':unittest.main()
