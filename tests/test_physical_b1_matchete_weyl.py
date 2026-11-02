"""B1 independent convention and negative-control tests (no native proof)."""
from __future__ import annotations
import json
from pathlib import Path
import tempfile
import sys
import unittest
import numpy as np

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
from Numerical.lo_certification.physical_normalisation.b1_matchete_weyl import (
    B1Error,four_spinor_witness,spurion_charge,triplet_metric,construct_components,
    model_specs,relative,congruence,effective_tree_kernel,
    b0_finite_independent,diagonal_toy_self_energy)
from scripts.RunPhysicalB1WeylAudit import guard_source
from scripts.RunPhysicalB1WeylAudit import current_native

class WeylB1Test(unittest.TestCase):
    def test_exact_dirac_and_majorana_four_spinor_identity(self):
        w=four_spinor_witness()
        self.assertTrue(w['exact_chosen_basis_checked'])
        self.assertFalse(w['absolute_native_GammaCC_phase_measured'])

    def test_rephasing_spurions_and_wrong_lambda_conjugation(self):
        for majorana in (False,True):
            witness=spurion_charge(majorana)
            self.assertEqual(witness['charge'],[0,0,0,0])
            self.assertNotEqual(witness['wrong_lambda_charge'],[0,0,0,0])

    def test_catalogue_has_14_dirac_2_majorana(self):
        rows=list(model_specs())
        self.assertEqual(len(rows),16)
        self.assertEqual({name for name,(_,_,df,a) in rows if df%2 and a==-1},
                         {'B_-1','C_-1'})

    def test_nontrivial_complex_dirac_flavour_matrix_rephasing(self):
        rng=np.random.default_rng(6107)
        X=rng.normal(size=(3,3))+1j*rng.normal(size=(3,3))
        X+=3*np.eye(3)
        M=np.block([[np.zeros((3,3)),X],[X.T,np.zeros((3,3))]])
        V=rng.normal(size=(3,6))+1j*rng.normal(size=(3,6))
        D=np.diag(np.exp(1j*rng.uniform(-3,3,size=6)))
        trans,_,_=congruence(M,V[:,:,None],np.zeros((3,6,1)),D,np.eye(1))
        new=V@D
        self.assertLess(relative(effective_tree_kernel(M,V),effective_tree_kernel(trans,new)),1e-12)
        self.assertLess(relative(M,M.T),1e-14)
        # Wrong heavy flavour transpose is not a harmless convention.
        malformed=np.block([[np.zeros((3,3)),X],[X,np.zeros((3,3))]])
        self.assertGreater(np.linalg.norm(malformed-malformed.T),1.e-2)

    def test_sign_change_with_compensating_field_change_is_harmless(self):
        # Flipping chi changes y1, MF simultaneously; not a new physical theory.
        M=np.array([[0,4+1j],[4+1j,0]],complex)
        V=np.array([[1+2j,3-.5j],[.2+.7j,-.3+1.1j]],complex)
        D=np.diag([1.,-1.])
        transformed,_,_=congruence(M,V[:,:,None],np.zeros((2,2,1)),D,np.eye(1))
        self.assertLess(relative(effective_tree_kernel(M,V),effective_tree_kernel(transformed,V@D)),1e-13)
        self.assertGreater(relative(effective_tree_kernel(M,V),effective_tree_kernel(M,V@D)),1e-2)

    def test_wrong_majorana_half_changes_physical_mass(self):
        M=np.array([[120+0j]])
        V=np.array([[1+1j],[2-.2j]])
        reference=effective_tree_kernel(M,V)
        missing_half=effective_tree_kernel(2*M,V)
        self.assertLess(relative(missing_half,reference/2),1e-14)
        self.assertGreater(relative(reference,missing_half),.49)

    def test_yukawa_conjugation_changes_complex_flavour_invariants(self):
        # Wrong y1 H.c. conjugation generally changes Takagi singular values.
        u=np.array([.8+1.1j,-.3+.7j,.6-.4j])
        v=np.array([.5-1j,1.3+.1j,-.1+.6j])
        correct=np.outer(u,v)+np.outer(v,u)
        wrong=np.outer(u.conj(),v)+np.outer(v,u.conj())
        self.assertGreater(np.linalg.norm(np.sort(np.linalg.svd(correct)[1])-
                                          np.sort(np.linalg.svd(wrong)[1])),.01)

    def test_triplet_metric_and_wrong_delta_charge_pairing(self):
        root=np.sqrt(2.)
        U=np.array([[1/root,0,1/root],[1j/root,0,-1j/root],[0,1,0]],complex)
        t3=np.array([[0,-1j,0],[1j,0,0],[0,0,0]],complex)
        G=triplet_metric(U,t3)
        self.assertLess(relative(G,np.array([[0,0,1],[0,1,0],[1,0,0]])),1e-13)
        self.assertGreater(relative(G,np.eye(3)),.5)
        self.assertLess(np.linalg.norm(G-G.T),1e-14)
        with self.assertRaises(B1Error):
            triplet_metric(np.eye(3),t3)
        Uphase=U@np.diag([1,1j,1])
        self.assertAlmostEqual((triplet_metric(Uphase,t3))[1,1].real,-1.)

    def test_scalar_YZ_swap_requires_conjugation_of_scalar_data(self):
        # Under conjugating the scalar, exchanging Y,Z is permitted only with
        # A -> A.T and B -> B.conj(); it is not a universal standalone swap.
        A=np.array([[2, .3+.7j],[.3-.7j,4]],complex)
        B=np.array([[.1+.2j, .7-.3j],[.7-.3j,-.2+.1j]],complex)
        self.assertGreater(np.linalg.norm(A-A.T),.1)
        self.assertGreater(np.linalg.norm(B-B.conj()),.1)
        self.assertLess(np.linalg.norm(A.T-A.conj()),1e-15)
        # Two mixed complex scalars, fixed A: swapping only Y/Z changes
        # even the magnitude of a one-flavour Majorana self-energy.
        Aphy=np.array([[40000,14000+10000j],[14000-10000j,90000]],complex)
        h=np.block([[Aphy.real,-Aphy.imag],[Aphy.imag,Aphy.real]])
        ms2,O=np.linalg.eigh(h)
        C=(O[:2,:]+1j*O[2:,:])/np.sqrt(2)
        Y=np.array([[[.7+1.1j,.2-.3j]]])
        Z=np.array([[[.4-.8j,-.9+.4j]]])
        reference=diagonal_toy_self_energy(Y@C+Z@C.conj(),[110.],ms2,110.**2)[0,0]
        swapped=diagonal_toy_self_energy(Z@C+Y@C.conj(),[110.],ms2,110.**2)[0,0]
        self.assertGreater(abs(abs(reference)-abs(swapped)),1e-3)

    def test_loop_no_spurious_majorana_two_and_canonical_real_scalar(self):
        m=100.
        f1=b0_finite_independent(m*m,200.*200.,100.*100.)
        f2=b0_finite_independent(m*m,300.*300.,100.*100.)
        h=.2+.1j
        # Real and imaginary scalar vertices h/sqrt2, i*h/sqrt2.
        g=np.array([[[h/np.sqrt(2.),1j*h/np.sqrt(2.)]]])
        ans=diagonal_toy_self_energy(g,[m],[200.**2,300.**2],100.**2)[0,0]
        expected=h*h*m*(f1-f2)/(32*np.pi**2)
        self.assertAlmostEqual(abs(ans-expected),0,places=13)
        self.assertGreater(abs(ans-2*expected),1e-8)

    def test_unsupported_source_drift_is_rejected(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp)
            p=root/'Lagrangian/interactions/T3Topology.wl'
            p.parent.mkdir(parents=True)
            p.write_text('changed',encoding='utf-8')
            with self.assertRaises(B1Error):
                guard_source(root)

    def test_native_certificate_promotion_blocked(self):
        with tempfile.TemporaryDirectory() as tmp:
            p=Path(tmp)
            for name in ('A_p0','B_m1','C_m1'):
                (p/(name+'.json')).write_text(json.dumps({
                  'format':'t3-b1-native-restricted-action-v1',
                  'model':name,'status':'NativeActionTermsIsolated_PhysicalDictionaryPending',
                  'physical_source_certified':True,
                  'source_gammaCC_concrete_matrix_measured':False,'tests':{'action_success':True},
                  'yukawas':{'y1':{'original_term_count':1,'hc_term_count':1},
                             'y2':{'original_term_count':1,'hc_term_count':1}}
                }),encoding='utf-8')
            with self.assertRaises(B1Error):
                current_native(p)

    def test_construction_uses_distinct_weyl_blocks_complex_couplings(self):
        # Synthetic 2 x 1 x 2 source tensor. Only charge-neutral component
        # used for this test; these are not asserted Matchete CG measurements.
        def vertex(dims,entries):
            return {'dimensions':list(dims),'components':[
                {'indices':list(i),'cg':str(v)} for i,v in entries]}
        # B_-1: dF=1; qF=0. S1,Y=-1/2: weight 0 at s=0;
        # S2,Y=+1/2: weight 0 at s=1.
        case={'model':{'dS1':2,'dS2':2,'dF':1,'alpha':-1},
              'vertices':{'T3Y1CG':vertex((2,1,2),[((1,1,1),'1')]),
                          'T3Y2CG':vertex((2,1,2),[((1,1,2),'1')])}}
        ys={'real':np.eye(3).tolist(),'imag':(np.eye(3)*.17).tolist()}
        ms={'real':(np.eye(3)*110).tolist(),'imag':np.zeros((3,3)).tolist()}
        benchmark={'representation':{'d_s1':2,'d_s2':2,'d_f':1,'alpha':-1},
                   'base_state':{'ordinary':{'y1':ys,'y2':ys,'MF':ms}}}
        made=construct_components(case,benchmark,{1:np.eye(1),2:np.eye(2)},'B_-1')
        self.assertEqual(made['fermion_type'],'Majorana')
        self.assertEqual(made['M'].shape,(3,3))
        self.assertAlmostEqual(made['Z'][0,0,0],-(1-.17j))
        self.assertAlmostEqual(made['Y'][0,0,3],-(1-.17j))
        fake=construct_components(case,benchmark,{1:np.eye(1),2:np.eye(2)},'B_-1',mass_factor=2)
        self.assertGreater(relative(fake['M'],made['M']),.3)

if __name__=='__main__':unittest.main()
