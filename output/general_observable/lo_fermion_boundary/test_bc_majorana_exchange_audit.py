"""Tests for the B_-1/C_-1 signed Majorana p=0 native-source audit.
Project path: output/general_observable/lo_fermion_boundary/test_bc_majorana_exchange_audit.py
"""
from __future__ import annotations
from copy import deepcopy
from fractions import Fraction
import unittest

import bc_majorana_exchange_audit as bc


def example_term(model, row):
    ix, ys, spur, scalars, sign, denominator, numer, cg = row
    left = '- (' if sign < 0 else '('
    if numer:
        left += 'Sqrt[3]*'
    parts=[]
    for s in spur.split():
        parts.append(f'Bar[Coupling[T3MEU{s}, {{}}, 0]]')
    for y in ys:
        parts.append(f'Bar[Coupling[y{y}, {{Index[p, Flavor], Index[r, NFlavor]}}, 0]]')
    for sc in scalars:
        parts.append(f'Bar[Field[NewScalar{sc}, Scalar, {{}}, {{}}]]')
    for marker in bc.CG_MARKER[cg]:
        parts.append(marker + ', {Index[i, SU2L[fund]]}]')
    parts.append('NCM[Transp[Field[l, Fermion, {Index[i, SU2L[fund]], Index[p, Flavor]}, {}]], DiracProduct[GammaCC, Proj[-1]], Field[l, Fermion, {Index[j, SU2L[fund]], Index[q, Flavor]}, {}]]')
    raw=left+'*'.join(parts)+')'+bc.DENOMS[denominator]
    return {'source_index':ix,'chirality':'PL','source_term_inputform':raw,
            'factors':[{'name':f'y{y}'} for y in ys],
            'spurion_names':[f'T3MEU{s}' for s in spur.split()]}


class TestBCMajorana(unittest.TestCase):
    def test_exact_signed_majorana_saddle(self):
        r=bc.exact_majorana_saddle()
        self.assertEqual(r['independent_diagonal_channels'],4)
        self.assertEqual(r['independent_cross_channels'],6)

    def test_nontrivial_signs_and_fractions(self):
        for mass, weights in ((Fraction(7),(1,-2,3,-4)),
                              (Fraction(13,2),(-3,5,-7,11)),
                              (Fraction(25,3),(7,4,-3,-5))):
            self.assertTrue(bc.exact_majorana_saddle(mass,weights)['zero_momentum_saddle_verified'])

    def test_bad_mass_rejected(self):
        for mass in (0,-2):
            with self.assertRaises(ValueError):bc.exact_majorana_saddle(mass)

    def test_all_22_source_fixture_signatures(self):
        for model, rows in bc.ROWS.items():
            for row in rows:
                with self.subTest(model=model, index=row[0]):
                    self.assertEqual(bc.verify_native_term(example_term(model,row), model,row)['source_index'],row[0])

    def test_c_bridge_not_misclassified_as_half_current(self):
        for row in bc.ROWS['C_-1'][-2:]:
            self.assertEqual(row[5],'2M')
            self.assertTrue(row[6])
            self.assertEqual(row[7],'Bridge')

    def test_mutated_sign_is_blocked(self):
        row=bc.ROWS['C_-1'][1]
        term=example_term('C_-1',row)
        term['source_term_inputform']=term['source_term_inputform'][2:]
        with self.assertRaises(ValueError):bc.verify_native_term(term,'C_-1',row)

    def test_mutated_denominator_is_blocked(self):
        row=bc.ROWS['B_-1'][1]
        term=example_term('B_-1',row)
        term['source_term_inputform']=term['source_term_inputform'].replace('/(2*Coupling[MF, {}, 0])','/Coupling[MF, {}, 0]')
        with self.assertRaises(ValueError):bc.verify_native_term(term,'B_-1',row)

    def test_mutated_scalar_is_blocked(self):
        row=bc.ROWS['B_-1'][0]
        term=example_term('B_-1',row)
        term['source_term_inputform']=term['source_term_inputform'].replace('NewScalarY','NewScalarX')
        with self.assertRaises(ValueError):bc.verify_native_term(term,'B_-1',row)

    def test_mutated_spurion_is_blocked(self):
        row=bc.ROWS['C_-1'][0]
        term=example_term('C_-1',row)
        term['source_term_inputform']=term['source_term_inputform'].replace('T3MEU11','T3MEU22')
        with self.assertRaises(ValueError):bc.verify_native_term(term,'C_-1',row)

    def test_wrong_projector_is_blocked(self):
        row=bc.ROWS['B_-1'][0]
        term=example_term('B_-1',row)
        term['source_term_inputform']=term['source_term_inputform'].replace('Proj[-1]','Proj[1]')
        with self.assertRaises(ValueError):bc.verify_native_term(term,'B_-1',row)

    def test_claims_are_fail_closed(self):
        self.assertFalse(bc.exact_majorana_saddle()['native_GammaCC_to_master_Weyl_phase_certified'])
        self.assertEqual(len(bc.ROWS['B_-1']),10)
        self.assertEqual(len(bc.ROWS['C_-1']),12)


if __name__ == '__main__':
    unittest.main(verbosity=2)
