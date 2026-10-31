"""Standalone mega batch contracts; source tests opt in if local files present."""
from __future__ import annotations
import json
from pathlib import Path
import sys
import tempfile
import unittest
from dataclasses import dataclass
import numpy as np

ROOT=Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:sys.path.insert(0,str(ROOT))
from Numerical.experimental import T3BCMegaPhysicalBridge as m
from Numerical.experimental import T3BCComplexUVStress as uv


@dataclass
class Block:
    name:str
    start:int
    stop:int
    kind:str='complex_matrix'
    shape:tuple=(3,3)

class Layout:
    def __init__(self):
        self.blocks={k:Block(k,18*i,18*(i+1)) for i,k in enumerate(('MF','y1','y2'))}
        self.names=tuple(self.blocks)
    def block(self,n):return self.blocks[n]

class MatrixAndActionTests(unittest.TestCase):
    def test_complex_3x3_roundtrip(self):
        x=np.arange(9).reshape(3,3)+1j*np.arange(9).reshape(3,3)/17
        np.testing.assert_array_equal(m.unpack_matrix(m.pack_matrix(x)),x)
    def test_wrong_schema(self):
        with self.assertRaises(ValueError):m.unpack_matrix({'real':[[0]],'imag':[[0]]})
    def test_complex_action_unitary_covariance(self):
        r=m.complex_spinor_action_test(trials=100)
        self.assertLess(r['max_correct_action_relative_residual'],1e-12)
        self.assertGreater(r['wrong_same_V_negative_control_max_relative'],0.1)
    def test_mass_summary_reference(self):
        K=np.diag([2e-15,1e-15,0.1e-15]).astype(complex)
        x=m.mass_summary(K,246.22)
        self.assertTrue(x['masses_descending_ev'][0]>x['masses_descending_ev'][1])
        self.assertAlmostEqual(x['masses_descending_ev'][0],0.0606242884,places=7)
    def test_nonsymmetric_neutrino_matrix_rejected(self):
        K=np.eye(3,dtype=complex);K[0,2]=0.1
        with self.assertRaises(ValueError):m.mass_summary(K,246.22)
    def test_fixed_random_complex_unitary(self):
        Q=uv.unitary(14)
        self.assertGreater(np.linalg.norm(Q.imag),0.5)
        self.assertLess(m.relative(Q.conj().T@Q,np.eye(3)),1e-12)
    def test_rge_state_basis_transform_is_invertible(self):
        v=np.arange(54,dtype=float)/50
        L=Layout();Q=uv.unitary(31)
        vt=uv.change_basis(v,L,Q)
        back=uv.change_basis(vt,L,Q.conj().T)
        self.assertLess(np.linalg.norm(back-v),1e-11)
    def test_matrices_correct_rotation(self):
        v=np.arange(54,dtype=float)/50
        L=Layout();Q=uv.unitary(34)
        original={k:uv.block_complex(v,L,k) for k in L.names}
        vt=uv.change_basis(v,L,Q)
        self.assertLess(m.relative(uv.block_complex(vt,L,'MF'),Q.T@original['MF']@Q),1e-12)
        self.assertLess(m.relative(uv.block_complex(vt,L,'y1'),original['y1']@Q.conj()),1e-12)
        self.assertLess(m.relative(uv.block_complex(vt,L,'y2'),original['y2']@Q),1e-12)
    def test_wrong_y1_transformation_detected(self):
        v=np.arange(54,dtype=float)/50;L=Layout();Q=uv.unitary(27)
        good=uv.change_basis(v,L,Q)
        wrong=uv.change_basis(v,L,Q,y1_conjugated=False)
        self.assertGreater(np.linalg.norm(good-wrong),.05)
    def test_beta_covariance_zero_under_simple_covariant_model(self):
        L=Layout();Q=uv.unitary(12)
        z=np.arange(54,dtype=float)/100
        def rhs(x):return 0.071*x
        rotated=uv.change_basis(z,L,Q)
        expected=uv.change_basis(rhs(z),L,Q)
        actual=rhs(rotated)
        self.assertLess(np.linalg.norm(expected-actual),1e-12)
    def test_no_physical_promotion(self):
        gates=m.blocked_physical_gates()
        self.assertIs(gates['physical_C5_certified'],False)
        self.assertIs(gates['same_trajectory_sequential_matching_performed'],False)

class InputGates(unittest.TestCase):
    def test_source_sha_mismatch_caught(self):
        with tempfile.TemporaryDirectory() as td:
            src=Path(td)/'source.json';one=Path(td)/'candidate.json'
            src.write_text('{}');one.write_text('{}')
            with self.assertRaises(ValueError):m.check_running_pair('B_-1',m.load(src),m.load(one),source_path=src,running_path=one)
    def test_live_topology_head_if_present(self):
        path=ROOT/'Lagrangian/interactions/T3Topology.wl'
        if not path.is_file():self.skipTest('full GitHub checkout not mounted')
        x=m.check_topology_orientation(path)
        self.assertEqual(x['y1_source_heavy_spinor'],'CConj[F]')
    def test_uploaded_reference_if_present(self):
        # The test works with user-provided benchmark JSON when present in
        # this local repository; it is optional on a clean source checkout.
        p=ROOT/'output/studies/t3_bc_running_one_scale_triangle_v1/B_m1_running_one_scale.json'
        if not p.is_file():self.skipTest('previous B/C one-scale output not mounted')
        a=m.load(p)
        c5=m.unpack_matrix(a['candidate_ordered_source_C5_3x3_gev_inverse_NOT_PHYSICAL'])
        r=m.mass_summary(c5,a['reference_EW_vev_gev'])
        np.testing.assert_allclose(r['masses_descending_ev'],
                                   a['reference_mnu_singular_ev_descending_NOT_PHYSICAL'],rtol=1e-11)
    def test_live_sm_weinberg_no_k_if_present(self):
        try:
            from Numerical.running.SMWeinbergEvolution import SMWeinbergInitialConditions,evolve_sm_weinberg
        except ImportError:self.skipTest('source SM+Weinberg RGE not mounted')
        eye=np.eye(3,dtype=complex)
        st=SMWeinbergInitialConditions(.37,.61,.88,.27,eye*.01,eye*.2,eye*.01,np.zeros((3,3),complex))
        out=evolve_sm_weinberg(st,1e4,1e3,rtol=1e-8,atol=1e-21)
        self.assertLess(np.linalg.norm(out.K),1e-30)

if __name__=='__main__':unittest.main()
