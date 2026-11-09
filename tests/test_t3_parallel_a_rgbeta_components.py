"""Small tests of the conditional source-expression interpreter and source guards."""
import math
import unittest

from Numerical.parallel_t3.wilson.rgbeta_components import (
    ComponentFrameError, delta, epsilon, symmetric_cg, evaluate_inputform, triple_invariant,
    validate_frame, reconstruct,
)


class TestParallelARGBetaComponents(unittest.TestCase):
    def test_diagonal_delta_is_one_not_rgbeta_einstein_dimension(self):
        self.assertEqual(evaluate_inputform('del[SU2L[S2], #1, #2] &', (1,1,1,1)),1)
        self.assertEqual(evaluate_inputform('del[SU2L[S2], #1, #2] &', (1,2,1,1)),0)
        self.assertEqual(delta('fund',1,1),1)

    def test_epsilon_orientation_explicit(self):
        self.assertEqual(epsilon('fund',1,2),1)
        self.assertEqual(epsilon('fund',2,1),-1)
        self.assertEqual(evaluate_inputform('del[SU2L[fund], #4, #1]*eps[SU2L[fund], #2, #3] &',
                                            (1,1,2,1)),1)

    def test_s2_projector_and_normalization(self):
        for value in validate_frame().values():
            self.assertLess(value,1e-14)
        self.assertAlmostEqual(symmetric_cg('SU2L',2,1,2),1/math.sqrt(2))
        self.assertEqual(evaluate_inputform('delS2[SU2L, #4, #1, #2] &', (1,2,1,2)),1/math.sqrt(2))


    def test_reused_internal_triplet_index_and_project_triple(self):
        self.assertAlmostEqual(evaluate_inputform(
            'delS2[SU2L, a, #1, #2]*delS2[SU2L, a, #3, #4] &',
            (1,2,1,2)),0.5)
        for a in (1,2,3):
            self.assertAlmostEqual(triple_invariant(a,a,a),0.0)
        # This source-defined triple tensor exists in the two E-class records.
        val=evaluate_inputform(
            'delS2[SU2L, a, #1, #2]*T3RGBetaS2TripleInvariant[a, #3, #4] &',
            (1,1,1,1))
        self.assertTrue(math.isfinite(val))

    def test_multi_delta_symmetrized_potential(self):
        s='(1/3)*(del[SU2L[S2], #1, #2]*del[SU2L[S2], #3, #4]) + (1/3)*' \
          '(del[SU2L[S2], #1, #3]*del[SU2L[S2], #2, #4]) + (1/3)*' \
          '(del[SU2L[S2], #1, #4]*del[SU2L[S2], #2, #3]) &'
        self.assertAlmostEqual(evaluate_inputform(s,(1,1,2,2)),1/3)
        self.assertAlmostEqual(evaluate_inputform(s,(1,1,1,1)),1)

    def test_safe_parser_blocks_unknown_heads_and_code(self):
        for s in ('exp[1] &','__import__("os").system("ls") &','del[SU2L[S2],#1,#2] &;Bad[]',
                  'del[SU2L[S2],#1,#2]', 'Sqrt[2] &'):
            with self.subTest(s=s), self.assertRaises(ComponentFrameError):
                evaluate_inputform(s,(1,1,1,1))

    def test_mini_inventory_and_sha(self):
        d={'format':'stage22-native-rgbeta-declarations-v1','status':'NativeDeclarationMetadataOnly',
           'models':{'A_-4':{'tensors_not_exported':True,'declarations':[{
               'symbol':'x','fields_inputform':'{H, Bar[H], S1, Bar[S1]}',
               'invariant_inputform':'del[SU2L[fund],#1,#2] &',
               'coupling_factor_inputform':'6*Matrix[x][##1]&',
               'unique_arrangements_inputform':'{{1,2,3,4}}',
               'self_conjugate':False}]}}}
        r={'models':{'A_-4':{'dimensions':[1,3,2], 'total_additional_quartics':1,
                           'new_rgbeta_symbols':[{'name':'x'}], 'preexisting_rgbeta_symbols':[]}}}
        z=reconstruct(d,r,declarations_sha256='d',registry_sha256='r',strict_all16=False)
        self.assertEqual(z['quartic_slot_count'],1)
        self.assertEqual(z['models']['A_-4']['records'][0]['dimensions'],[2,2,1,1])
        self.assertFalse(z['interpretation']['physical_matchete_to_RGBeta_dictionary_certified'])
        with self.assertRaises(ComponentFrameError):
            reconstruct(d,r,declarations_sha256='d',registry_sha256='r',strict_all16=True)


if __name__=='__main__': unittest.main()
