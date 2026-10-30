from scripts.RebuildT3FrozenFullScalarRefitDisplay import parse_grid
from scripts.GenerateT3FrozenSensitivityScans import apply_direction
import numpy as np

def test_grid():
    assert parse_grid('0.95,1.0,1.05')==[0.95,1.0,1.05]

def test_direction_copy():
    x=np.ones(18)
    out=apply_direction(x,'T1',1.05)
    assert np.allclose(out[:3],1.05)
    assert np.allclose(out[3:],1)
    assert np.all(x==1)
