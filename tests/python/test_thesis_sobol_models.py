"""Pure parameter-profile tests; does not require Mathematica or fitting data."""
from importlib.util import module_from_spec, spec_from_file_location
from pathlib import Path
import pytest

SCRIPT = Path(__file__).resolve().parents[2] / 'scripts/RunThesisSobolModels.py'
spec = spec_from_file_location('thesis_sobol_models', SCRIPT)
mod = module_from_spec(spec)
spec.loader.exec_module(mod)


def fixture():
    return {
        'kind':'t3_pipeline_numerical_results_v1',
        'base_state':{
            'ordinary':{
                'lambdaT3':{'real':0.1,'imag':0.0},
                'y1':{'real':[[0.0001, 0, 0], [0,0.0001,0], [0,0,0.0001]],
                      'imag':[[0,0,0]]*3},
                'y2':{'real':[[0.0002, 0, 0], [0,0.0002,0], [0,0,0.0002]],
                      'imag':[[0,0,0]]*3},
            }
        }
    }


def test_models_16_and_reference_scalar_values():
    assert len(mod.ALL_MODELS) == 16
    assert mod.key_for('A_-2') == 'T3_dS1_1_dS2_3_dF_2_alpha_m2'
    assert mod.key_for('E_+2') == 'T3_dS1_3_dS2_3_dF_2_alpha_p2'
    assert mod.lambda_for('A_-2', .1, .2, True) == pytest.approx(.2)
    assert mod.lambda_for('E_+2', .05, .2, True) == pytest.approx(.1)
    assert mod.lambda_for('E_+2', .05, None, True) == pytest.approx(.05)


def test_fixed_scalar_and_rescaled_seed():
    cfg = mod.build_search_payload(fixture(), 'A_-2', .4, n_points=128,
                                  target_chi2=.01, seed=7)
    assert cfg['base_state']['ordinary']['lambdaT3'] == {'real':.4, 'imag':0}
    assert cfg['base_state']['ordinary']['y1']['real'][0][0] == pytest.approx(.00005)
    assert len(cfg['benchmark_search']['parameters']) == 18
    assert all(x['name'] != 'lambdaT3_real' for x in cfg['benchmark_search']['parameters'])
    assert cfg['benchmark_search']['local_refinement']['active_counts'] == [5,8,12]
    assert cfg['sensitivity']['n_points'] == 11


def test_optional_scalar_optimisation_and_no_mutation():
    ref = fixture()
    cfg = mod.build_search_payload(ref, 'A_-2', .1, n_points=64,
                                  target_chi2=.01, seed=7, fit_lambda=True)
    assert len(cfg['benchmark_search']['parameters']) == 19
    assert cfg['benchmark_search']['parameters'][-1]['binding'] == 'ordinary.lambdaT3.real'
    assert ref['base_state']['ordinary']['lambdaT3']['real'] == .1


def test_invalid_profile():
    with pytest.raises(ValueError):
        mod.key_for('A_-1')
    with pytest.raises(ValueError):
        mod.build_search_payload(fixture(),'A_-2',.1,n_points=100,target_chi2=.01,seed=1)
    with pytest.raises(ValueError):
        mod.lambda_for('A_-2',.1,-.2,True)
