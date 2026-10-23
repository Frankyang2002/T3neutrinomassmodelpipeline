import numpy as np
import sympy as sp
from validation.DeriveT3BCrossQuarticMap import numerical_cross_quartic_tensors,solve_map
from Numerical.running.CoupledT3BWilson import _complex,_pair

def test_rgbeta_doublet_invariants_independent():
    a,b=numerical_cross_quartic_tensors()
    assert a and b
    assert sp.Matrix([[a.get(k,0),b.get(k,0)] for k in sorted(set(a)|set(b))]).rank()==2

def test_exact_identity_map():
    basis=numerical_cross_quartic_tensors()
    m,n=solve_map(basis,basis)
    assert m==sp.eye(2)
    assert n>0

def test_complex_io():
    assert _complex([1,-2])==1-2j
    assert _pair(1-2j)==[1.,-2.]
