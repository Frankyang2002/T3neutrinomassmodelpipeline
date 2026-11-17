"""Safe numerical evaluation of saved Matchete expressions; no live matching.

Absolute mass reference values are not frozen here: they depend on the installed
illustrative configuration and saved source coefficients. Retain mathematical
consistency and coupling-scaling checks instead of an obsolete snapshot.
"""
from __future__ import annotations

import copy
from pathlib import Path
import sys
import mpmath as mp
import pytest

ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))
from scripts import RunGeneralT3NumericalC5 as num


def basic_expr():
    y = 'Coupling[y1, {Index[d$$1, Flavor], Index[d$$1, NFlavor]}, 0]'
    z = 'Coupling[y2, {Index[d$$2, Flavor], Index[d$$1, NFlavor]}, 0]'
    return (
        f'hbar*Bar[{y}]*{z}*Coupling[lambdaT3, {{}}, 0]/Coupling[MF, {{}}, 0]'
        ' + Sqrt[3]*Log[Coupling[MSX, {}, 0]/Coupling[MSY, {}, 0]]'
    )


def test_conversion_exact_matchete_subset():
    expr = num.compile_coefficient(basic_expr())
    mp.mp.dps = 50
    val = eval(expr, {'__builtins__': {}}, {
        'conj': mp.conj, 'sqrt': mp.sqrt, 'log': mp.log,
        'hbar': mp.mpf(1)/157,
        'y1a': mp.mpc('1','2'), 'y2b': mp.mpc('3','4'),
        'lambdaT3': mp.mpc('0.1'), 'MF': mp.mpf('100000'),
        'MSX': mp.mpf('1000'), 'MSY': mp.mpf('2000'),
    })
    reference = (
        mp.conj(mp.mpc(1, 2))*mp.mpc(3, 4)*mp.mpf('0.1')/(157*100000)
        + mp.sqrt(3)*mp.log(mp.mpf(1)/2)
    )
    assert abs(val-reference) < mp.mpf('1e-45')


@pytest.mark.parametrize('invalid', [
    "__import__('os').system('echo nope')", 'Foo[1]',
    'Coupling[unknown, {}, 0]', 'hbar /. x -> y', 'Log[MF;4]', '1/0 + x',
])
def test_unknown_input_fails_closed(invalid):
    with pytest.raises(ValueError, match='fail-closed'):
        num.compile_coefficient(invalid)


def test_example_schema_and_results_if_present(tmp_path):
    source = ROOT / 'output/general_t3_mass_eigen_matchete/one_loop_probe'
    params = ROOT / 'output/general_t3_mass_eigen_matchete'
    if not (source / 'B_m1_final_one_loop_probe.json').is_file():
        pytest.skip('Real Matchete coefficient files not installed')
    config = num.read_json(ROOT / 'configs/general_t3/c5_illustrative.json')
    for model in ('B_-1', 'C_-1'):
        result, out = num.numerical_model(model, config, source, params, tmp_path)
        assert out.is_file()
        assert result['physical_c5_certified'] is False
        assert result['sequential_EFT1_performed'] is False
        assert result['symmetry_max_abs_C5_gev_inv'] == 0
        assert result['max_abs_mnu_eV'] > 0
        masses = result['takagi_singular_values_eV_descending']
        assert len(masses) == 3
        assert masses[0] > 0
        assert masses[2] < 1e-10
        # Turning off the interaction must remove the topology.
        no_topology = copy.deepcopy(config)
        no_topology['quartics']['lambdaT3']['real'] = 0.0
        zero, _ = num.numerical_model(model, no_topology, source, params, tmp_path)
        assert zero['max_abs_mnu_eV'] < 1e-30
        # Doubling each Yukawa gives quadratic scaling of the mass matrix.
        bigger = copy.deepcopy(config)
        for matrix in (bigger['y1'], bigger['y2']):
            for row in matrix:
                for value in row:
                    value['real'] *= 2
        scaled, _ = num.numerical_model(model, bigger, source, params, tmp_path)
        assert scaled['max_abs_mnu_eV'] == pytest.approx(
            4*result['max_abs_mnu_eV'], rel=1e-10,
        )


def test_no_arbitrary_source_replacements():
    assert num.ALLOWED_COUPLINGS == {
        'MF','MSX','MSY','T3MEU11','T3MEU12','T3MEU21','T3MEU22',
        'lambdaT3','lambdaT3ExtraQ020020I1','lambdaT3ExtraQ020200I1',
        'y1','y2',
    }


def test_complex_json_and_shapes():
    assert num.complex_value({'real': 2, 'imag': 3}, 'z') == mp.mpc(2, 3)
    with pytest.raises(ValueError):
        num.matrix_values([[{'real':1,'imag':0}], [{'real':1,'imag':0}]],
                          'y', rows=3)
    with pytest.raises(ValueError):
        num.complex_value({'real':'nan','imag':0}, 'z')
