"""Offline tests for historical-input paired B/C pilot. No physics run needed."""
from __future__ import annotations

from copy import deepcopy
import math
from pathlib import Path

import numpy as np
import pytest

from scripts import RunT3GeneralScalarBCPilot as p
from scripts import RunT3ZeroExtraBoundaryBenchmark as b7

NEW=('lambdaT3GeneralQ110101I1',)
REUSED=('lambdaHHdagS2S2','lambdaHHdagS1barS1bar',
 'lambdaHHdagS1barS2barCross','lambdaS1bar2S2bar2',
 'lambdaS1barS2S2bar2','lambdaS1S1bar2S2bar')


def snap(lam=0.5, extra=0j,mu=0j,m1=1e6,m2=1.4e6):
    return dict(mu_gev=1e3,mS1Sq_gev2=m1,mS2Sq_gev2=m2,
                lambdaT3=b7.complex_json(lam),
                quadratic_mu2_gev2=b7.complex_json(mu),
                new_quartics={k:b7.complex_json(extra) for k in NEW},
                reused_quartics={k:b7.complex_json(extra) for k in REUSED})


def traj(end=None):
    return dict(eft1_end=snap() if end is None else end)


def test_case_inventory():
    assert len(p.CASES)==3
    assert p.CASES['small_uv_running_5e-9'][1] is True
    assert p.CASES['small_uv_frozen_5e-9'][1] is False
    assert p.CASES['large_uv_frozen_0p5'][2] == 0.5
    assert p.MODELS==('B_-1','C_-1')


def test_run_id_safe():
    assert p.parse_run_id('pilot_2026-10-28')=='pilot_2026-10-28'
    for name in ('../old','/tmp/x','C:/tmp/x','','a'*65,'a b','test\\x'):
        with pytest.raises(ValueError):p.parse_run_id(name)


def test_refuse_existing_result(tmp_path):
    x=p.check_destination('small_uv_running_5e-9','r1',root=tmp_path)
    assert x.parent.parent.name=='small_uv_running_5e-9' or x.parent.name=='small_uv_running_5e-9'
    x.mkdir(parents=True)
    with pytest.raises(ValueError,match='already exists'):
        p.check_destination('small_uv_running_5e-9','r1',root=tmp_path)


def test_input_fractions():
    for v in (-0.1,float('nan'),float('inf'),1.01):
        with pytest.raises(ValueError):p.finite_fraction(v,'q',maximum=1)
    assert p.finite_fraction(.01,'mu2',maximum=.2)==.01


def test_zero_scenario_really_zero():
    q,mu=p.added_boundary(NEW,REUSED,5e-9,1e6,1.4e6,
         scenario='restricted_zero',quartic_fraction=.1,mixing_fraction=.01)
    assert len(q)==7 and all(v==0 for v in q.values()) and mu==0


def test_quartics_only_scenario():
    q,mu=p.added_boundary(NEW,REUSED,5e-9,1e6,1.4e6,
         scenario='general_quartics',quartic_fraction=.1,mixing_fraction=.01)
    assert len(q)==7 and all(math.isclose(v.real,5e-10,rel_tol=1e-12) for v in q.values())
    assert mu==0


def test_mixed_scenario():
    q,mu=p.added_boundary(NEW,REUSED,.5,1e6,1.4e6,
         scenario='general_quartics_and_mixing',quartic_fraction=.1,mixing_fraction=.01)
    assert all(v==.05 for v in q.values())
    assert mu==1e4


def test_bad_inventory_fails():
    with pytest.raises(ValueError,match='inventory'):
        p.added_boundary(NEW,REUSED[:-1],.5,1e6,1.4e6,
            scenario='restricted_zero',quartic_fraction=.1,mixing_fraction=.01)


def test_bad_masses_fail():
    with pytest.raises(ValueError,match='mass'):
        p.added_boundary(NEW,REUSED,.5,0.,1e6,
            scenario='restricted_zero',quartic_fraction=.1,mixing_fraction=.01)


def test_rg_comparison_zero_self():
    data=traj()
    row=p.compare_rg(data,deepcopy(data),'B_-1')
    assert row['lambdaT3_relative_change']==0
    assert row['scalar_diagonal_m2_relative_shifts']==[0,0]
    assert all(v['delta_abs']==0 for v in row['quartic_changes'].values())


def test_rg_changes():
    a=traj();s=snap(lam=.51,m1=1.01e6,m2=1.5e6,extra=.05)
    row=p.compare_rg(a,traj(s),'C_-1')
    assert math.isclose(row['lambdaT3_relative_change'],.02)
    assert math.isclose(row['scalar_diagonal_m2_relative_shifts'][0],.01)
    assert len(row['quartic_changes'])==7


def test_missing_quartic_fails():
    a=traj();b=deepcopy(a)
    b['eft1_end']['reused_quartics'].pop(REUSED[0])
    with pytest.raises(ValueError,match='inventory'):p.compare_rg(a,b,'B_-1')


def test_mass_mixing_diagnostic():
    sp=p.mass_spectrum('B_-1',snap(mu=1e4))
    assert sp['mixing_sin_2theta']>0
    assert sp['mu2_over_min_diagonal']==.01
    assert sp['physical_matchete_coupling_map_certified'] is False


def test_tachyon_blocked():
    with pytest.raises(ValueError):p.mass_spectrum('C_-1',snap(mu=2e6))


def test_c5_unmixed_projection_reports_generated_mu2():
    a=traj(snap(mu=0j));b=traj(snap(mu=1e4))
    info=p.unmixed_c5_projection_info(a,b)
    assert info['projected_mu2_gev2']=={'real':0.,'imag':0.}
    assert info['any_generated_mixing_omitted'] is True
    assert info['actual_mixing_contribution_to_c5_evaluated'] is False
    assert not info['endpoint_mixing']['restricted_zero']['endpoint_mixing_discarded_in_c5_projection']
    assert info['endpoint_mixing']['general_quartics']['endpoint_mixing_discarded_in_c5_projection']
    assert info['endpoint_mixing']['general_quartics']['actual_endpoint_mu2_over_min_diagonal_m2']==.01


def test_c5_unmixed_projection_with_two_generated_mu2():
    info=p.unmixed_c5_projection_info(traj(snap(mu=2e3)),traj(snap(mu=1e4)))
    assert info['any_generated_mixing_omitted']
    assert info['endpoint_mixing']['restricted_zero']['actual_endpoint_mu2_gev2']['real']==2e3
    assert info['endpoint_mixing']['general_quartics']['actual_endpoint_mu2_gev2']['real']==1e4


def test_c5_unmixed_projection_rejects_unphysical_diagonals():
    with pytest.raises(ValueError,match='diagonal mass'):
        p.unmixed_c5_projection_info(traj(),traj(snap(m1=-1)))


def test_summary_labels_noncertified():
    a=traj(); b=traj(snap(extra=.01,mu=1e4))
    compare=p.compare_rg(a,b,'B_-1')
    r={'chi2_fit':33.0,'passes_original_fit_screen':False,
       'quartic_vs_zero':compare,
       'scenarios':{'restricted_zero':a, 'general_quartics_and_mixing':
                     {'spectrum':p.mass_spectrum('B_-1',snap(mu=1e4))}}}
    md=p.format_summary({'case':'large_uv_frozen_0p5','run_id':'r1'},{'B_-1':r})
    assert 'NO' in md and 'uncertified' in md.lower() and 'sin(2 theta)' in md
    assert 'not a' in md.lower() or 'not an' in md.lower()


def test_zero_and_mixed_scenario_distinct():
    a,b=p.added_boundary(NEW,REUSED,.5,1e6,1.4e6,
        scenario='general_quartics_and_mixing',quartic_fraction=.1,mixing_fraction=.01)
    c,d=p.added_boundary(NEW,REUSED,.5,1e6,1.4e6,
        scenario='general_quartics',quartic_fraction=.1,mixing_fraction=.01)
    assert a==c and b!=d and d==0


def test_summary_c5_preview_displays_value():
    a=traj()
    comparison=p.compare_rg(a,a,'C_-1')
    row={'chi2_fit':0.1,'passes_original_fit_screen':True,
         'quartic_vs_zero':comparison,'scenarios':{},
         'C5_preview':{'epsilon_C5_frobenius_uncertified':0.0045}}
    md=p.format_summary({'case':'small_uv_frozen_5e-9','run_id':'r2'},{'C_-1':row})
    assert '| 0.0045 |' in md
    assert 'Shift abs(lambdaT3)' in md


def test_summary_distinguishes_rge_generated_mu2_without_explicit_input():
    a=traj(snap(mu=0j)); b=traj(snap(mu=1e4))
    a['spectrum']=p.mass_spectrum('B_-1',a['eft1_end'])
    b['spectrum']=p.mass_spectrum('B_-1',b['eft1_end'])
    row={'chi2_fit':1.,'passes_original_fit_screen':True,
         'quartic_vs_zero':p.compare_rg(a,b,'B_-1'),
         'scenarios':{'restricted_zero':a,'general_quartics':b}}
    md=p.format_summary({'case':'small_uv_frozen_5e-9','run_id':'r3'},{'B_-1':row})
    assert 'RGE-generated mixed quadratic' in md
    assert 'Zero quadratic *input*' in md
    assert '0.01' in md
    assert 'mu2 explicitly projected to zero' in md


def test_one_scale_c5_preview_runs_with_rge_generated_mixing(monkeypatch,tmp_path):
    """Exercises the full preview path with stubbed one-scale source services.

    This is a control-flow/contract test, NOT the real Matchete coefficient.
    """
    import sys
    import types
    source=tmp_path/'B_m1_final_one_loop_probe_weinberg_ordered_coefficient.inputform.txt'
    source.write_text('synthetic source for unit test',encoding='utf-8')
    matched=types.ModuleType('scripts.RunGeneralT3NumericalC5')
    matched.DEFAULT_SOURCE=tmp_path
    fake_expression=('(lambdaT3*y1a*y2b + '
        'lambdaT3ExtraQ020020I1*y1a*y1b + '
        'lambdaT3ExtraQ020200I1*y2a*y2b)*hbar/MF')
    matched.compile_coefficient=lambda _expr:compile(fake_expression,'<unit-test>','eval')
    monkeypatch.setitem(sys.modules,'scripts.RunGeneralT3NumericalC5',matched)
    # Only import-dependent services are substituted. The actual preview
    # arithmetic, mixing bookkeeping and JSON report creation are exercised.
    for module_name in ('Numerical','Numerical.experimental'):
        if module_name not in sys.modules:
            package=types.ModuleType(module_name)
            package.__path__=[]
            monkeypatch.setitem(sys.modules,module_name,package)
    basis_mod=types.ModuleType('Numerical.experimental.GeneralT3MixedScalarMassBasis')
    class FakeBasis:
        @staticmethod
        def build(_model,m1,m2,mu):
            assert mu==0j, 'Preview must explicitly use the unmixed projection'
            return types.SimpleNamespace(eigenvalues_sq=np.array(sorted([m1,m2])),
                                         U=np.eye(2,dtype=complex))
    basis_mod.MixedDoubletMassBasis=FakeBasis
    monkeypatch.setitem(sys.modules,basis_mod.__name__,basis_mod)
    def mat(a):
        v=np.asarray(a,dtype=float)
        return {'real':v.tolist(),'imag':np.zeros_like(v).tolist()}
    record={'cfg':{'base_state':{'ordinary':{
        'MF':mat(np.diag([1e5,1.1e5,1.2e5])),
        'y1':mat(np.diag([1,2,3])),
        'y2':mat(np.diag([.5,1,1.5]))}}}}
    a=traj(snap(lam=.1,mu=0j))
    b=traj(snap(lam=.1,extra=.01,mu=1e4))
    result=p.one_scale_c5_preview('B_-1',record,a,b)
    assert result['epsilon_C5_frobenius_uncertified']>0
    assert result['actual_endpoint_mixing_included_in_c5'] is False
    assert result['unmixed_projection_info']['any_generated_mixing_omitted'] is True
    assert result['unmixed_projection_info']['endpoint_mixing']['general_quartics'][
       'actual_endpoint_mu2_gev2']['real']==1e4
