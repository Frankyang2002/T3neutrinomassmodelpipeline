from types import SimpleNamespace
import math
import numpy as np
from unittest.mock import patch
from scripts.DiagnoseT3GeneralScalarScreens import screen_details


def test_pass_and_failed_components():
    uv=SimpleNamespace(mS1Sq=1e6,mS2Sq=1e6,
                       y1=np.zeros((3,3)),y2=np.zeros((3,3)))
    eft=object()
    with patch('Numerical.running.FermionThresholdBoundary.fermion_threshold_masses',
               return_value=SimpleNamespace(masses_gev=np.array([1e5,1e5,1e5]))), \
         patch('Numerical.running.IntermediateScalarRGBetaEvaluator.intermediate_scalar_state_environment',
               return_value={'lambdaH':0.13}):
        ok=screen_details(uv,eft,[1000,1200],(),())
        assert ok['passed']
        bad=screen_details(uv,eft,[100,6000],(),())
        assert not bad['passed']
        assert 'scalar_min_ge_0p5TeV' in bad['failed_checks']
        assert 'scalar_max_le_5TeV' in bad['failed_checks']
