"""Exact synthetic checks of quartic Cartesian conversion and factorials."""
import unittest
import sympy as sp
from Numerical.running.T3MatcheteCartesianResponses import (
    convert_record, CartesianConversionError,
)


def rec(species, counts, dims, entries, symbol='c'):
    return {'symbol':symbol,'species':species,'counts':counts,
            'dimensions':dims,'components':entries,'potential_prefactor':'1'}


def entries(r,part):
    return {tuple(item['indices']):sp.sympify(item['lambda_abcd'])
            for item in r['responses'][part]['entries']}


class TestExactQuarticConversion(unittest.TestCase):
    def test_abs_phi_four(self):
        # S1*S1*S1dag*S1dag in A_-4 (S1 singlet),
        # L = c |phi|^4 + hc; V(real c=1)=-0.5*(r^2+i^2)^2
        x=rec(['S1','S1','S1dag','S1dag'],[0,0,2,2,0,0],
              [1,1,1,1],['1'])
        row=convert_record('A_-4',x)
        self.assertEqual(entries(row,0),{
            (4,4,4,4):sp.Integer(-12),
            (4,4,5,5):sp.Integer(-4),
            (5,5,5,5):sp.Integer(-12),
        })
        self.assertEqual(entries(row,1),{})

    def test_holomorphic_phi_four(self):
        # L=c*S1^4+hc, dV/dIm(c)= 2*(r^3 i-r i^3)
        x=rec(['S1']*4,[0,0,4,0,0,0],[1]*4,['1'])
        row=convert_record('A_-4',x)
        self.assertEqual(entries(row,1),{
            (4,4,4,5):sp.Integer(12),
            (4,5,5,5):sp.Integer(-12),
        })
        self.assertEqual(entries(row,0),{
            (4,4,4,4):sp.Integer(-12),
            (4,4,5,5):sp.Integer(12),
            (5,5,5,5):sp.Integer(-12),
        })

    def test_hhss_doublet(self):
        # H Hdag S1 S1dag: scalar doublet 2, S1 singlet 1
        x=rec(['H','Hdag','S1','S1dag'],[1,1,1,1,0,0],
              [2,2,1,1],['1','0','0','1'])
        row=convert_record('A_-4',x)
        real=entries(row,0)
        self.assertEqual(real[(0,0,4,4)],sp.Integer(-2))
        self.assertEqual(real[(2,2,4,4)],sp.Integer(-2))
        self.assertEqual(real[(0,0,5,5)],sp.Integer(-2))
        self.assertEqual(entries(row,1),{})

    def test_species_mismatch_fails(self):
        x=rec(['S1','S1dag','S1','S1dag'],[0,0,2,2,0,0],
              [1,1,1,1],['1'])
        with self.assertRaises(CartesianConversionError):
            convert_record('A_-4',x)

    def test_symmetry_fails(self):
        x=rec(['H','H','S1','S1dag'],[2,0,1,1,0,0],
              [2,2,1,1],['0','1','0','0'])
        with self.assertRaises(Exception):
            convert_record('A_-4',x)


if __name__=='__main__':
    unittest.main()
