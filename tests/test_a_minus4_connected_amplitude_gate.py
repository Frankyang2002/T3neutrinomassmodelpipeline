"""Targeted Batch-1 tests for source CG and two CC-sign exchange gate.

Full integration against local original source is performed by the CLI; these
unit tests run without Matchete or saved audit status flags.
"""
from __future__ import annotations

import copy
from fractions import Fraction
import json
from pathlib import Path
import unittest

import numpy as np

from Numerical.lo_certification.spinor import a_minus4_connected_amplitude_gate as gate
from Numerical.lo_certification.spinor.spinors import G, wedge


MIX = ('SparseArray[Automatic, {2, 2, 3}, 0, '
       '{1, {{0, 3, 6}, {{1, 1}, {1, 3}, {2, 2}, {1, 2}, {2, 1}, {2, 3}}}, '
       '{3^(-1/4), I/3^(1/4), -3^(-1/4), -3^(-1/4), -3^(-1/4), I/3^(1/4)}}]')
Y1 = ('SparseArray[Automatic, {2, 2}, 0, '
      '{1, {{0, 1, 2}, {{2}, {1}}}, {1, -1}}]')
Y2 = ('SparseArray[Automatic, {2, 2, 3}, 0, '
      '{1, {{0, 3, 6}, {{1, 2}, {2, 1}, {2, 3}, {1, 1}, {1, 3}, {2, 2}}}, '
      '{3^(-1/4), 3^(-1/4), (-I)/3^(1/4), '
      '3^(-1/4), I/3^(1/4), -3^(-1/4)}}]')


def fixture() -> tuple[dict, dict]:
    records = [{'Name':nm, 'TensorInputForm':expr,
                'RepsInputForm':rep} for nm,expr,rep in (
        ('T3Y1CG',Y1,'{SU2L[fund], T3BSMd2}'),
        ('T3Y2CG',Y2,'{SU2L[fund], Bar[T3BSMd2], T3BSMd3}'),
        ('T3MixCG',MIX,'{Bar[SU2L[fund]], Bar[SU2L[fund]], T3BSMd3}'))]
    return ({'status':'Success','metadata':{'dS1':1,'dS2':3,'dF':2,'alpha':-4},
             'CGRegistry': records},
            {'CGRegistry':[records[2]]})


class NativeConnectedGate(unittest.TestCase):
    def test_native_sparse_reader_yukawa(self):
        uv,_ = fixture()
        matrix=gate.native_cg_sparse(uv['CGRegistry'][0],'T3Y1CG',(2,2))
        np.testing.assert_array_equal(matrix,np.array([[0,1],[-1,0]],complex))

    def test_yukawa_cg_contraction_equals_six_native_matched_entries(self):
        uv,seed=fixture()
        arrays,info=gate.native_cg_identity(uv,seed)
        self.assertEqual(info['UV_Y1_Y2_CG_to_Matchete_F_only_max_abs_residual'],0.0)
        self.assertEqual(info['matched_nonzero'],6)
        np.testing.assert_allclose(arrays['Gmatch'],
                                   np.einsum('ah,bhc->abc',
                                             arrays['G1'].conj(),arrays['G2'].conj()),
                                   atol=1e-15,rtol=0)

    def test_native_y2_phase_flip_detected(self):
        uv,seed=fixture()
        uv['CGRegistry'][1]['TensorInputForm']=uv['CGRegistry'][1]['TensorInputForm'].replace(
            '(-I)/3^(1/4)','I/3^(1/4)')
        with self.assertRaisesRegex(gate.GateBlocked,'does not reproduce'):
            gate.native_cg_identity(uv,seed)

    def test_matched_tensor_sign_flip_detected(self):
        uv,seed=fixture()
        seed['CGRegistry'][0]['TensorInputForm']=seed['CGRegistry'][0]['TensorInputForm'].replace(
            'I/3^(1/4)', '(-I)/3^(1/4)', 1)
        with self.assertRaises((gate.GateBlocked,ValueError)):
            gate.native_cg_identity(uv,seed)

    def test_sparse_pointer_corruption_fails(self):
        uv,_=fixture()
        uv['CGRegistry'][0]['TensorInputForm']=uv['CGRegistry'][0]['TensorInputForm'].replace(
            '{{0, 1, 2}', '{{0, 2, 1}')
        with self.assertRaisesRegex(gate.GateBlocked,'pointers'):
            gate.native_cg_sparse(uv['CGRegistry'][0],'T3Y1CG',(2,2))

    def test_exact_two_current_fermionic_saddle_both_signs(self):
        result=gate.exact_two_current_saddle()
        self.assertTrue(result['UV_EFT_saddle_equal_for_both_signs'])
        self.assertEqual(result['fixed_action_relative_dirac_tree_coefficient'],'1/M')
        self.assertFalse(result['mapping_native_Dirac_action_to_declared_Weyl_components_certified'])

    def test_wrong_dirac_mass_sign_does_not_match(self):
        m=Fraction(13)
        A=(G.var('a0'),G.var('a1'))
        B=(G.var('b0'),G.var('b1'))
        eta=tuple(-x*(1/m) for x in B)
        chi=tuple(x*(1/m) for x in A)
        wrong=+m*wedge(eta,chi)+wedge(eta,A)+wedge(chi,B)
        expect=wedge(A,B)*(1/m)
        self.assertNotEqual(wrong,expect)

    def test_complex_offdiagonal_UV_EFT_both_signs_equal_but_opposite(self):
        arrays,_=gate.native_cg_identity(*fixture())
        result=gate.sample_connected_amplitudes(arrays)
        samples=result['candidate_amplitudes']
        self.assertEqual({s['C_choice'] for s in samples},{'+','-'})
        self.assertTrue(all(s['relative_amplitude_residual']<5e-13 for s in samples))
        self.assertFalse(result['native_phase_discriminated'])
        self.assertFalse(result['physical_master_C_tensor_certified'])
        a=complex(**samples[0]['EFT_Matchete_PL_contact'])
        b=complex(**samples[1]['EFT_Matchete_PL_contact'])
        self.assertGreater(abs(a),1e-7)
        self.assertAlmostEqual(abs(a+b),0,places=13)

    def test_native_action_requires_actual_yukawa_terms(self):
        # Tampering with literal Matchete action is rejected before any
        # guessed spinor phase can enter the numerical amplitude.
        action={'native_yukawa_1_terms_inputform':['Coupling[y1]','Bar[Coupling[y1]]'],
                'native_yukawa_2_terms_inputform':['Coupling[y2]','Bar[Coupling[y2]]'],
                'native_heavy_free_inputform':'-MF Bar[F] F',
                'native_concrete_4x4_cc_matrix':None,
                'native_4x4_CC_to_Weyl_epsilon_ratio':None,
                'native_connected_UV_amplitude_calculated':False}
        with self.assertRaisesRegex(gate.GateBlocked,r'Missing holomorphic HC Yukawas|action changed'):
            gate.validate_actions(action)

    def test_unsupported_source_shape_rejected(self):
        uv,_=fixture()
        with self.assertRaisesRegex(gate.GateBlocked,'wrong native CG shape'):
            gate.native_cg_sparse(uv['CGRegistry'][1],'T3Y2CG',(3,2,3))


if __name__=='__main__':
    unittest.main()
