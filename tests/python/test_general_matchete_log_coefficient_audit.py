"""Pure regression tests for independent one-scale log coefficient extraction.

Does not require wolframscript, expensive general EFT1 states, or saved run data.
"""
import math
from pathlib import Path
import sys
import numpy as np
import pytest

ROOT=Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:sys.path.insert(0,str(ROOT))

from Numerical.running.GeneralT3MatcheteLogCoefficientAudit import (
    MODELS,MIXED,replace_common_mass_logs,compare_coefficients,
    original_one_scale_parity,
)

MF='Coupling[MF, {}, 0]'
MS1='Coupling[MS1, {}, 0]'
MS2='Coupling[MS2, {}, 0]'
MSX='Coupling[MSX, {}, 0]'


def test_all_16_representations_present():
    assert len(MODELS)==16 and len(set(MODELS))==16
    assert MIXED=={'B_-1','C_-1'}


def test_two_regular_logs_both_shift_and_leave_nonlogs_unchanged():
    src=f'3 + 4*Log[{MF}^2/{MS1}^2] - 2*Log[{MF}^2/{MS2}^2]'
    shifted,inventory=replace_common_mass_logs(src)
    assert inventory['log_occurrences']==2
    assert inventory['scalar_mass_symbols']=={'MS1':1,'MS2':1}
    assert '3 + 4*' in shifted
    assert shifted.count('Log[2*(')==2
    assert src.count('Log[')==shifted.count('Log[')


def test_mixed_and_split_mass_logs():
    src=(f'Log[{MF}^2/{MSX}^2] + '
         f'Log[{MF}^2/Coupling[T32PlusMass, {{}}, 0]^2]')
    s,info=replace_common_mass_logs(src)
    assert 'MSX' in s and 'T32PlusMass' in s
    assert info['scalar_mass_symbols']=={'MSX':1,'T32PlusMass':1}


def test_reject_nonmass_ratio_log_to_prevent_false_spurion():
    with pytest.raises(ValueError,match='non-heavy/scalar'):
        replace_common_mass_logs('Log[Coupling[lambdaT3, {}, 0]]')
    with pytest.raises(ValueError,match='non-heavy/scalar'):
        replace_common_mass_logs(f'Log[{MS1}^2/{MS2}^2]')


def test_unbalanced_log_rejected():
    with pytest.raises(ValueError,match='unbalanced'):
        replace_common_mass_logs(f'Log[{MF}^2/{MS1}^2')


def test_no_log_rejected():
    with pytest.raises(ValueError,match='no recognised'):
        replace_common_mass_logs('1/(MF^2-MS1^2)')


def test_non_integer_shift_rejected():
    with pytest.raises(ValueError):
        replace_common_mass_logs(f'Log[{MF}^2/{MS1}^2]',multiplier=1)


def test_unity_coefficient_normalization():
    m=np.array([[1.+0.3j,2j,0.5], [2j,-1.3+0j,0],[0.5,0,1]],complex)*1e-16
    z=compare_coefficients(m,m)
    assert abs(z['best_fit_complex_factor']['real']-1)<1e-14
    assert abs(z['best_fit_complex_factor']['imag'])<1e-14
    assert z['relative_matrix_difference_at_assumed_factor_1']<1e-14
    assert z['relative_matrix_residual_after_best_fit_factor']<1e-14


def test_signed_factor_and_complex_phase_detected():
    m=np.array([[1.,0.2],[0.2,2]],complex)
    m=np.pad(m,((0,1),(0,1)))*1e-16
    factor=-2+0.4j
    report=compare_coefficients(factor*m,m)
    assert abs(report['best_fit_complex_factor']['real']+2)<1e-13
    assert abs(report['best_fit_complex_factor']['imag']-0.4)<1e-13
    assert report['relative_matrix_residual_after_best_fit_factor']<1e-14
    assert report['relative_matrix_difference_at_assumed_factor_1']>1


def test_orthogonal_matrix_mismatch_not_absorbed_by_factor():
    x=np.eye(3,dtype=complex)*1e-15
    y=x+np.diag([2,-2,0])*1e-15
    report=compare_coefficients(y,x)
    assert report['relative_matrix_residual_after_best_fit_factor']>0.3


def test_original_matchete_reproduction():
    a=np.eye(3)*1e-15
    assert original_one_scale_parity(a,a)==0.
    assert abs(original_one_scale_parity(1.01*a,a)-.01)<1e-13


def test_common_log_shift_mathematical_identity_for_linear_log():
    # Synthetic symbolic coefficient A*Log(MF^2/MS^2)+B.
    A=complex(1.23,-0.4);B=complex(2.,0.5)
    MFmass=1e5;MSmass=1e3
    original=A*math.log(MFmass**2/MSmass**2)+B
    shifted=A*math.log(2*MFmass**2/MSmass**2)+B
    assert abs((shifted-original)/math.log(2)-A)<1e-14


def test_wolfram_channel_driver_preserves_source_and_coupling_points(tmp_path, monkeypatch):
    """Synthetic end-to-end runner contract; no external Wolfram executable."""
    import json
    import subprocess
    import types
    from Numerical.running import GeneralT3MatcheteLogCoefficientAudit as audit

    src=tmp_path/'coefficient.txt';src.write_text('dummy coefficient')
    executable=tmp_path/'repo'/'Lagrangian'/'EvaluateGeneralC5Channels.wl'
    executable.parent.mkdir(parents=True)
    executable.write_text('synthetic evaluator fixture')
    fake=types.ModuleType('Numerical.matching.GeneralT3FlavourC5')
    fake.EVALUATOR_WL=executable
    fake.flavour_c5_from_channels=lambda rows,ys,channels: np.eye(3,dtype=complex)*float(rows[0]['real'])
    monkeypatch.setitem(sys.modules,'Numerical.matching.GeneralT3FlavourC5',fake)
    seen=[]

    def run(command, **kwargs):
        assert command[:3]==['fake-wolfram','-file',str(executable)]
        assert kwargs['cwd']==executable.parents[1]
        assert kwargs['check'] is False
        points=json.loads(Path(command[4]).read_text())
        assert points['hbar_value']==pytest.approx(1/(16*math.pi**2))
        assert len(points['points'])==3
        Path(command[5]).write_text(json.dumps({
            'results':[{'real':4.5}], 'channels':['y1Bar__y2Bar']}))
        seen.append(True)
        return types.SimpleNamespace(returncode=0,stdout='ok',stderr='')
    monkeypatch.setattr(subprocess,'run',run)
    got=audit._run_wolfram_channels(src,[{'index':i,'couplings':{'MF':1e5}}
        for i in range(3)],{},wolframscript='fake-wolfram',timeout_sec=15)
    assert seen==[True]
    assert np.array_equal(got,4.5*np.eye(3))
