"""Offline SU(2) invariance contract tests for the traced bad Matchete CG."""
from pathlib import Path
import sys
import tempfile
import unittest
from fractions import Fraction

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
from scripts.DiagnoseT3A4FourDoubletCG import (
    parse_four_bar_tensor, audit_tensor, audit_trace_file, run,
)

SPARSE=(
  'HoldCompleteForm[DefineCG[cg1, '
  '{Bar[SU2L[fund]], Bar[SU2L[fund]], Bar[SU2L[fund]], '
  'Bar[SU2L[fund]]}, SparseArray[Automatic, {2, 2, 2, 2}, 0, '
  '{1, {{0, 3, 6}, {{1, 1, 1}, {2, 1, 2}, {2, 2, 1}, '
  '{1, 1, 2}, {1, 2, 1}, {2, 2, 2}}}, '
  '{2/Sqrt[3], 1/Sqrt[3], 1/Sqrt[3], 1/Sqrt[3], '
  '1/Sqrt[3], 2/Sqrt[3]}}]]]'
)

class TestFourBarInvariance(unittest.TestCase):
    def test_trace_compressed_matchete_format(self):
        a=parse_four_bar_tensor(SPARSE)
        self.assertEqual(a['name'],'cg1')
        self.assertEqual(len(a['components']),6)
        self.assertEqual(a['components'][(1,1,1,1)],Fraction(2))
        self.assertEqual(a['components'][(2,2,2,2)],Fraction(2))

    def test_maximal_torus_is_not_invariant(self):
        a=audit_tensor(parse_four_bar_tensor(SPARSE))
        self.assertFalse(a['maximal_torus_invariant'])
        self.assertEqual(a['maximal_torus_weight_violations'],['1111','2222'])

    def test_exact_singlet_projection(self):
        a=audit_tensor(parse_four_bar_tensor(SPARSE))
        self.assertEqual(a['singlet_projection_norm_fraction_exact'],'1/9')
        self.assertEqual(a['squared_tensor_norm_exact'],'4')
        self.assertEqual(a['squared_singlet_projection_norm_exact'],'4/9')
        self.assertEqual(a['squared_non_singlet_norm_exact'],'32/9')

    def test_reject_unknown_coefficient(self):
        with self.assertRaisesRegex(ValueError,'Unsupported coefficient'):
            parse_four_bar_tensor(SPARSE.replace('2/Sqrt[3]','Sqrt[2]'))

    def test_reject_changed_representation(self):
        with self.assertRaisesRegex(ValueError,'Four Bar'):
            parse_four_bar_tensor(SPARSE.replace('Bar[SU2L[fund]]','SU2L[fund]',1))

    def test_two_independent_trace_files(self):
        with tempfile.TemporaryDirectory() as td:
            p=Path(td)
            for name in ('A_one_loop_DefineCG_calls_DIAGNOSTIC.txt',
                         'A_two_heavy_DefineCG_calls_DIAGNOSTIC.txt'):
                (p/name).write_text('# Matchete traced\n1\t'+SPARSE+'\n',encoding='utf-8')
            dst=p/'audit.json'
            result=run(p,dst)
            self.assertTrue(dst.is_file())
            self.assertEqual(result['trace_tensor_count_checked'],2)
            self.assertTrue(result['non_invariant_composite_cg_confirmed'])
            self.assertEqual(result['unparsed_calls'],0)
            self.assertFalse(result['physical_C5_certified'])

    def test_probe_source_does_not_override_matchete(self):
        src=(ROOT/'Lagrangian/RunT3A4FourDoubletCGIsolatedProbe.wl').read_text(encoding='utf-8')
        self.assertIn('InvariantTensors[SU[2],t3AlgebraReps]',src)
        self.assertIn('DefineCG[T3A4BadFourBarProbeCG',src)
        self.assertIn('DefineCG[T3A4GoodFourBarProbeCG',src)
        self.assertNotIn('Unprotect[',src)
        self.assertNotIn('UpValues[DefineCG]',src)
        self.assertNotIn('Match[',src)

if __name__=='__main__':unittest.main()
