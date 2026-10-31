"""Fully standalone regression tests; real fixture replay is optional."""
import unittest
from pathlib import Path
from tempfile import TemporaryDirectory
from fractions import Fraction
import sys
import json
from hashlib import sha256

from Numerical.experimental.T3BCPoleLogConsistency import (
    ConsistencyError, polynomial, terms, compare_residue_and_log, finite_groups,
)
from scripts.AuditT3BCPoleLogConsistency import run

L='Log[μbar2/Coupling[MSY, {}, 0]^2]'
A='Bar[Coupling[lambdaT3, {}, 0]]'
F='NCM[Field[l, Fermion, {Index[i, Flavor]}, {}], GammaCC]'

class TestParser(unittest.TestCase):
    def test_signed_terms(self):
        self.assertEqual(len(terms('- 3*Coupling[a, {}, 0] + 4*Coupling[b, {}, 0]')),2)
    def test_nested_minus_not_split(self):
        self.assertEqual(len(terms('-2*(3-1)*Field[H, Scalar, {}, {}]')),1)
    def test_exact_rationals(self):
        p=polynomial(('+(2*A[foo]*B[bar])/(3*MF[baz])', '-(A[foo]*B[bar])/(6*MF[baz])'),table={})
        self.assertEqual(next(iter(p.values())),Fraction(1,2))
    def test_negative_power(self):
        p=polynomial(('2*A[foo]^(-2)',),table={})
        self.assertEqual(next(iter(p.values())),Fraction(2,1))
    def test_noncommutative_container_opaque(self):
        p=polynomial(('NCM[A,B]',),table={})
        q=polynomial(('NCM[B,A]',),table={})
        # Different opaque objects are never compared by their arbitrary local IDs;
        # compare with a SHARED atom table to confirm ordering is preserved.
        table={}
        self.assertNotEqual(polynomial(('NCM[A,B]',),table=table),
                            polynomial(('NCM[B,A]',),table=table))
    def test_reject_functions_of_logs(self):
        with self.assertRaises(ConsistencyError):
            compare_residue_and_log('A[foo]',f'A[foo]*{L}*Log[other]','Y')
    def test_reject_pole_without_log(self):
        with self.assertRaises(ConsistencyError):
            compare_residue_and_log('A[foo]','0','Y')
    def test_matching_opposite_forms(self):
        a=compare_residue_and_log('-(A[foo]*B[bar])/(2*Q[z])',
                                  f'-(2*A[foo]*B[bar]*{L})/(4*Q[z])','Y')
        self.assertTrue(a['pole_equals_log_coefficient_in_opaque_factor_algebra'])
    def test_reject_coeff_mismatch(self):
        with self.assertRaises(ConsistencyError):
            compare_residue_and_log('A[foo]/2',f'A[foo]*{L}/3','Y')
    def test_reject_mixed_quartics(self):
        with self.assertRaises(ConsistencyError):
            finite_groups('Coupling[lambdaA, {}, 0]*Coupling[lambdaB, {}, 0]','LLHH')
    def test_group_quartic(self):
        s=f'{A}*{F} + 2*{A}*{F}'
        g=finite_groups(s,'LLHH')
        self.assertEqual(g['lambdaT3']['raw_additive_terms'],2)
        self.assertEqual(g['lambdaT3']['independent_opaque_monomials_after_exact_combining'],1)
    def test_reject_scale_in_constant(self):
        with self.assertRaises(ConsistencyError):
            finite_groups(f'{A}*{L}','LLHH')
    def test_unsupported_pow(self):
        with self.assertRaises(ConsistencyError):
            polynomial(('A[foo]^(B[b])',),table={})
    def test_reject_ast_call(self):
        with self.assertRaises(ConsistencyError):
            polynomial(('__import__(3)',),table={})

class TestRealArchivedData(unittest.TestCase):
    def test_replay_if_mounted(self):
        base=Path('/mnt/data/t3_finite_fixture/output/audits/t3_bc_scalar_finite_term_extraction_v1')
        if not (base/'summary.json').is_file():
            self.skipTest('Archived fixture not present on local Windows checkout')
        with TemporaryDirectory() as d:
            out=run(base,Path(d),('B_-1','C_-1'),root=None)
            self.assertEqual(len(out['sector_results']),8)
            self.assertTrue(all(x['pole_equals_log_coefficient_in_opaque_factor_algebra'] for x in out['sector_results']))

if __name__=='__main__':unittest.main()
