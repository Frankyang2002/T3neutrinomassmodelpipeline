import numpy as np
import pytest
from scripts.DiagnoseT3FullUVMassResponse import boundary_vector, scalar_values, response
class Block:
    def __init__(self,i):self.start=i;self.size=1
class Layout:
    def block(self,n):return Block({'mS1Sq':0,'mS2Sq':2}[n])
def test_scaling_only_two_squared_masses():
    a=np.array([4.,99.,9.,77.]);b=boundary_vector(a,Layout(),.25)
    assert b.tolist()==[1.,99.,2.25,77.]
    assert a.tolist()==[4.,99.,9.,77.]
def test_nonpositive_factors_rejected():
    with pytest.raises(ValueError):boundary_vector(np.array([4.,0.,9.]),Layout(),0)
def test_scalar_values():
    assert scalar_values(np.array([4.,7.,9.]),Layout())=={'mS1Sq':4.,'mS2Sq':9.}
def test_response():
    rows=[{'status':'Computed','factor_mSsq':1.,'scalar_masses_muS_tev':[40.,41.]},
          {'status':'Computed','factor_mSsq':.25,'scalar_masses_muS_tev':[39.,41.]}]
    r=response(rows)
    assert r['comparisons_to_baseline'][1]['delta_mass_TeV']==[-1.,0.]
