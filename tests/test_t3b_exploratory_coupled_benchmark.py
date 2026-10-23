import numpy as np
from Numerical.diagnostics.T3BExploratoryCoupledBenchmark import _gamma, coupled_heavy_resolved, _rhs

def test_one_generation_limit_of_diagonal_lift():
    c={'g2':0.6,'lambda12':0.1,'lambda12Adj':0.2,'lambdaH':0.3,
       'lambdaT3':1,'yu':np.diag([0,0,0.8]).astype(complex),
       'yd':np.diag([0,0,0.02]).astype(complex),
       'ye':np.diag([0,0,0.01]).astype(complex)}
    a,b,e=_gamma(c)
    assert np.isclose(a+e[2,2],-3*.6**2+2*.1-.2/2+.01**2)
    assert np.isclose(b-e[2,2],-3*.6**2+6*.8**2+6*.02**2+2*.3-.01**2)

def test_offdiagonal_exact_constant_solution():
    mus=[1e5,1e4,1e3]
    samples=[{'mu_gev':m,'g2':0.5,'lambdaH':0.2,'lambda12':0.0,
      'lambda12Adj':0.,'lambdaT3':[1.,0.],
      'yu':{'real':np.zeros((3,3)).tolist(),'imag':np.zeros((3,3)).tolist()},
      'yd':{'real':np.zeros((3,3)).tolist(),'imag':np.zeros((3,3)).tolist()},
      'ye':{'real':np.zeros((3,3)).tolist(),'imag':np.zeros((3,3)).tolist()}}
      for m in mus]
    y1=np.eye(3);y2=np.eye(3);masses=np.array([1e5]*3)
    _,c5=coupled_heavy_resolved(samples,y1,y2,masses,diagonal=False)
    from Numerical.diagnostics.T3BExploratoryCoupledBenchmark import LOOP,FWD,REV
    t=np.log(1e3/1e5)
    expected=np.sinh(np.sqrt(FWD*REV)*t/LOOP)*np.sqrt(FWD/REV)/1e5
    assert np.allclose(c5,np.eye(3)*expected,rtol=1e-7,atol=1e-14)
