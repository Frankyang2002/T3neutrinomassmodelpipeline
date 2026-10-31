"""Pure numerical regression and conditional live-source smoke tests."""
from __future__ import annotations

import json
from pathlib import Path
import sys
from types import ModuleType, SimpleNamespace
import unittest
from unittest.mock import patch

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:sys.path.insert(0,str(ROOT))

from Numerical.experimental.T3GeneralMatrixBeta import (
    EXTRA_MODELS, ExtendedUVEquationSystem, encode_matrix, decode_matrix,
    trial_extra_seed,
)
from Numerical.experimental.T3BCThresholdInteractionBridge import extract_complex


class FakeLayout:
    names=('y1','y2','MF')
    size=54


def make_meta(model,extras):
    return dict(GeneralScalarModel=model,GeneralScalarStage='UV',
        ReportBetaConvention='16*pi^2*dX/dln(mu); gauge BetaTerm divided by 2*g',
        GeneralScalarQuarticBetaBasis='RGBetaFullQuarticProjection',
        GeneralScalarNewBetaKeys=list(extras),
        GeneralScalarQuadraticMassMixing='NotIncluded')


class StateTests(unittest.TestCase):
    def setUp(self):
        self.y1=np.eye(3,dtype=complex)
        self.y2=np.diag([.3,.5,.7]).astype(complex)

    def test_six_extra_yukawa_models(self):
        self.assertEqual(len(EXTRA_MODELS),6)
        self.assertEqual(sum(v=='y1Conj' for v in EXTRA_MODELS.values()),3)
        self.assertEqual(sum(v=='y2Conj' for v in EXTRA_MODELS.values()),3)

    def test_trial_is_independent_and_zero_limit(self):
        z=trial_extra_seed(self.y1,self.y2,0)
        self.assertTrue(np.array_equal(z,np.zeros((3,3))))
        mat=trial_extra_seed(self.y1,self.y2,.3)
        self.assertEqual(mat.shape,(3,3))
        self.assertGreater(abs(mat[0,1].imag),0)
        self.assertFalse(np.allclose(mat,self.y1))

    def test_matrix_round_trip(self):
        m=self.y1+1j*self.y2
        self.assertTrue(np.array_equal(m,decode_matrix(encode_matrix(m),'matrix')))
        with self.assertRaises(ValueError):decode_matrix({'real':[1],'imag':[2]},'matrix')

    def test_extended_block_packing(self):
        labels=('q0','q1')
        names=('y1','y2','MF','q0','q1','y1Conj')
        sys=ExtendedUVEquationSystem('A_+0',FakeLayout(),labels,'y1Conj',
            {name:'1' for name in names},make_meta('A_+0',labels))
        extra=self.y1+2j*self.y2
        initial=np.arange(54,dtype=float)
        v=sys.initial_vector(initial,{'q0':1-3j,'q1':-2+4j},extra)
        self.assertEqual(len(v),76)
        base,qs,y=sys.split(v)
        self.assertTrue(np.array_equal(base,initial))
        self.assertEqual(qs,{'q0':1-3j,'q1':-2+4j})
        self.assertTrue(np.array_equal(y,extra))
        with self.assertRaises(ValueError):sys.initial_vector(initial,{'q0':1},extra)

    def test_no_extra_yukawa_models_have_zero_sidecar(self):
        sys=ExtendedUVEquationSystem('B_-1',FakeLayout(),('q0',),None,
            {name:'1' for name in ('y1','y2','MF','q0')},make_meta('B_-1',('q0',)))
        self.assertEqual(sys.size,56)
        v=sys.initial_vector(np.zeros(54),{'q0':0.7+0.2j})
        self.assertIsNone(sys.split(v)[2])

    def test_missing_beta_rejected_before_integration(self):
        with self.assertRaisesRegex(ValueError,'mismatch packed UV state'):
            ExtendedUVEquationSystem('A_+0',FakeLayout(),(), 'y1Conj',
              {'y1':'1','y2':'1','MF':'1'},make_meta('A_+0',()))

    def test_wrong_physical_metadata_rejected(self):
        m=make_meta('B_-1',())
        m['GeneralScalarQuadraticMassMixing']='IncludedCandidate'
        with self.assertRaisesRegex(ValueError,'without extra quadratic'):
            ExtendedUVEquationSystem('B_-1',FakeLayout(),(),None,
                {name:'1' for name in FakeLayout.names},m)

    def test_beta_matrix_pack_using_existing_parser_contract(self):
        # Minimal mocks are only for independent module tests when the complete
        # user's repository is unavailable in the builder environment.
        beta=ModuleType('Numerical.core.BetaVector')
        def derivative(values,layout,divide_by_loop_factor):
            self.assertTrue(divide_by_loop_factor)
            return np.full(layout.size,.125)
        beta.beta_values_to_derivative=derivative
        evalmod=ModuleType('Numerical.core.RGBetaEvaluator')
        evalmod.state_environment=lambda state:{'y1':self.y1,'y2':self.y2,'MF':np.eye(3)}
        evalmod.evaluate_inputform_expression=lambda e,env:env[e] if e in env else (2j if e=='imaginary' else 3)
        stmod=ModuleType('Numerical.core.StateVector')
        stmod.unpack_uv_state=lambda arr,layout,mu_gev:SimpleNamespace(y1=self.y1,y2=self.y2)
        sys1=ExtendedUVEquationSystem('A_+0',FakeLayout(),('q0',),'y1Conj',
            {'y1':'y1','y2':'y2','MF':'MF','q0':'imaginary','y1Conj':'y1Conj'},
            make_meta('A_+0',('q0',)))
        v=sys1.initial_vector(np.ones(54),{'q0':.3},self.y2)
        with patch.dict(sys.modules,{'Numerical.core.BetaVector':beta,
                                     'Numerical.core.RGBetaEvaluator':evalmod,
                                     'Numerical.core.StateVector':stmod}):
            dy=sys1.derivative(np.log(1e7),v)
        self.assertEqual(dy.shape,(74,))
        self.assertTrue(np.all(dy[:54]==.125))
        self.assertAlmostEqual(dy[55],2/(16*np.pi**2))
        self.assertTrue(np.array_equal(dy[56:65],self.y2.real.flatten()/(16*np.pi**2)))

    def test_negative_metadata_not_relabelled_certified(self):
        self.assertEqual(extract_complex({'real':1,'imag':-2}),1-2j)
        with self.assertRaises(ValueError):extract_complex({'real':1})


class LocalSourceTests(unittest.TestCase):
    @unittest.skipUnless((ROOT/'Numerical/core/RGBetaEvaluator.py').exists(),
                         'Full repository not available in package-test environment')
    def test_actual_native_rgbeta_matrix_parser(self):
        from Numerical.core.RGBetaEvaluator import evaluate_inputform_expression
        a=np.array([[1,2j,0],[0,3,0],[0,0,4]],dtype=complex)
        expr='2*Matrix[y1, Trans[Bar[y1]], y1][gen[$i], heavy[$j]] + Tr[y1 . Trans[Bar[y1]]]*Matrix[y1][gen[$i], heavy[$j]]'
        result=evaluate_inputform_expression(expr,{'y1':a})
        expected=2*a@a.conj().T@a+np.trace(a@a.conj().T)*a
        self.assertTrue(np.allclose(result,expected))

    @unittest.skipUnless((ROOT/'scripts/CoupledUVRunning.py').exists() and
                         (ROOT/'Reports/output/workstream_a/integrated_general_yukawa_uv/summary.json').exists(),
                         'Current complete local repository with native UV exports not available')
    def test_latest_native_uv_source_inventory(self):
        from scripts.RunT3FullStackOvernight import (read_json,REGISTRY,model_specifications,preflight)
        registry=read_json(REGISTRY)
        checks=preflight(model_specifications(registry),registry)
        bad={k:v for k,v in checks.items() if v['status']=='Blocked'}
        self.assertFalse(bad,str(bad))


if __name__=='__main__':unittest.main()

class BridgeContractTests(unittest.TestCase):
    def test_source_mixed_doublet_bridge_with_mocked_native_basis(self):
        # Reuses the SOURCE interface, but supplies an independent minimal
        # numeric reference in this repository-free packaging environment.
        from Numerical.experimental.T3BCThresholdInteractionBridge import bridge_snapshot
        mod=ModuleType('Numerical.experimental.GeneralT3MixedScalarMassBasis')
        class Reference:
            def __init__(self,m1,m2,mu):
                self.A=np.array([[m1,mu],[mu.conjugate(),m2]])
                self.e,self.U=np.linalg.eigh(self.A)
                self.tachyonic=bool((self.e<0).any());self.degenerate=False
            @classmethod
            def build(cls,model,m1,m2,mu):return cls(m1,m2,mu)
            def original_fields(self,x,y):
                s1=self.U[0,0]*x+self.U[0,1]*y
                v2=self.U[1,0]*x+self.U[1,1]*y
                eps=np.array([[0,1],[-1,0]],complex)
                return s1,-eps@v2.conj()
            def eigen_fields(self,s1,s2):
                eps=np.array([[0,1],[-1,0]],complex)
                return tuple(self.U.conj().T@np.stack([s1,eps@s2.conj()]))
            def old_potential(self,s1,s2):
                eps=np.array([[0,1],[-1,0]],complex)
                z=np.stack([s1,eps@s2.conj()]);return float(np.real(np.einsum('ai,ab,bi->',z.conj(),self.A,z)))
            def new_potential(self,x,y):return float(self.e[0]*np.vdot(x,x).real+self.e[1]*np.vdot(y,y).real)
            def info(self):return {'masses_squared_gev2':self.e.tolist()}
        mod.MixedDoubletMassBasis=Reference
        mod.verify_tensor_transport=lambda fields,U,seed,trials:0.
        mod.verify_su2_covariance=lambda basis,trials:0.
        mod.quadratic_mass_tensor_transport=lambda basis:{'mass_tensor_diagonalization_residual_gev2':0.}
        m1,m2,mu=1.1e6,1.4e6,3.2e4-8.1e3j
        eigen=np.linalg.eigvalsh(np.array([[m1,mu],[mu.conjugate(),m2]]))
        snapshot={'mu_gev':1e5,'mS1Sq_gev2':m1,'mS2Sq_gev2':m2,
                  'muT3GeneralS1S2_gev2':{'real':mu.real,'imag':mu.imag},
                  'scalar_spectrum':{'mass_squared_eigenvalues_gev2':np.repeat(eigen,4).tolist()}}
        with patch.dict(sys.modules,{'Numerical.experimental.GeneralT3MixedScalarMassBasis':mod}):
            result=bridge_snapshot('B_-1',snapshot)
        self.assertTrue(result['valid_running_scalar_eigenbasis'])
        self.assertLess(result['quadratic_action_relative_residual_max'],1e-11)
        self.assertFalse(result['physical_matching_performed'])


class SyntheticCampaignTests(unittest.TestCase):
    def test_real_dop853_full_matrix_extra_yukawa_evolution(self):
        """ODE is genuinely integrated; only the source expressions are mocked."""
        import scripts.RunT3FullStackOvernight as campaign
        from Numerical.experimental.T3GeneralMatrixBeta import LOOP
        initial_matrices = {
            'y1': np.diag([.01,.02,.03]).astype(complex),
            'y2': np.diag([.04,.05,.06]).astype(complex),
            'MF': np.diag([1e5,1.1e5,1.2e5]).astype(complex),
        }
        def pack_matrix(m):
            return np.concatenate((m.real.reshape(-1),m.imag.reshape(-1)))
        baseline = np.concatenate(tuple(pack_matrix(initial_matrices[k]) for k in FakeLayout.names))
        eq = ExtendedUVEquationSystem('A_+0',FakeLayout(),('q0',),'y1Conj',
            {k:k for k in ('y1','y2','MF','q0','y1Conj')},
            make_meta('A_+0',('q0',)))
        def unpack(v,layout,mu_gev):
            coords={}
            for idx,k in enumerate(FakeLayout.names):
                a=np.asarray(v[18*idx:18*(idx+1)])
                coords[k]=a[:9].reshape(3,3)+1j*a[9:].reshape(3,3)
            return SimpleNamespace(**coords)
        evaluator=ModuleType('Numerical.core.RGBetaEvaluator')
        evaluator.state_environment=lambda state:{k:getattr(state,k) for k in FakeLayout.names}
        evaluator.evaluate_inputform_expression=lambda expr,env:env[expr]
        statevector=ModuleType('Numerical.core.StateVector')
        statevector.unpack_uv_state=unpack
        beta=ModuleType('Numerical.core.BetaVector')
        beta.beta_values_to_derivative=lambda vals,layout,divide_by_loop_factor: (
            np.concatenate([pack_matrix(vals[k])/LOOP for k in FakeLayout.names]))
        crossing=ModuleType('Numerical.experimental.GeneralT3RunningFermionThresholds')
        crossing.locate_fermion_crossings=lambda *a,**kw: {'crossings':[], 'unresolved_crossings':[]}
        crossing.mass_matrix_from_vector=lambda v,layout:unpack(v,layout,1e7).MF
        registry={'models':{'A_+0':{'alpha':0,'dimensions':[1,3,2]}}}
        with patch.dict(sys.modules,{
                'Numerical.core.RGBetaEvaluator':evaluator,
                'Numerical.core.StateVector':statevector,
                'Numerical.core.BetaVector':beta,
                'Numerical.experimental.GeneralT3RunningFermionThresholds':crossing,
             }),patch.object(campaign,'load_uv_inputs',return_value=(baseline,FakeLayout(),{'q0':.1+0j})),\
               patch.object(campaign,'source_equations',return_value=eq):
            report=campaign.solve_uv('A_+0',registry=registry,
                seed_fraction=.3,points=5,rtol=1e-9,atol=1e-11,convergence=True)
        self.assertEqual(report['status'],'CandidateFullMatrixRGBetaUVIntegrated')
        self.assertTrue(report['convergence']['passed'])
        self.assertEqual(len(report['trajectory']),5)
        self.assertGreater(report['nfev'],5)
        expected=np.exp(np.log(1e5/1e7)/LOOP)
        self.assertAlmostEqual(report['final']['extra_yukawa_frobenius_norm']/
            report['initial']['extra_yukawa_frobenius_norm'],expected,places=8)
        self.assertFalse(report['physical_neutrino_masses_computed'])
