"""Offline tests of persistent sampling and diverse local restart selection."""
import ast
from pathlib import Path
from types import SimpleNamespace
import numpy as np
import pytest
from scipy.stats import qmc

SOURCE=Path(__file__).resolve().parents[1]/'scripts'/'SearchT3FrozenPersistent.py'
TREE=ast.parse(SOURCE.read_text(encoding='utf-8'))
NODES=[n for n in TREE.body if isinstance(n,ast.FunctionDef) and n.name in {'sobol_batch','choose_seeds','parse_radii'}]
NS={'np':np,'qmc':qmc,'math':__import__('math'),
    'bulk':SimpleNamespace(NAMES=('y1_00','y2_00')),
    'eligible_candidate':lambda x:isinstance(x,dict) and x.get('chi2') is not None}
exec(compile(ast.Module(body=NODES,type_ignores=[]),str(SOURCE),'exec'),NS)

def row(chi,x,y=0):
    return {'chi2':chi,'parameters':{'y1_00':x,'y2_00':y}}

def test_deterministic_sobol():
    lo=np.array([-1.,0.]);hi=np.array([1.,2.])
    a=NS['sobol_batch'](lo,hi,16,47)
    assert a.shape==(16,2)
    assert np.array_equal(a,NS['sobol_batch'](lo,hi,16,47))
    assert not np.array_equal(a,NS['sobol_batch'](lo,hi,16,48))

def test_invalid_sobol_count():
    with pytest.raises(ValueError):NS['sobol_batch'](np.zeros(2),np.ones(2),15,4)

def test_diverse_seeds():
    rows=[row(5,0),row(6,.001),row(7,.002),row(30,.8),row(35,-.9)]
    selected=NS['choose_seeds'](rows,3,np.array([-1.,-1.]),np.array([1.,1.]),.1)
    assert [r['chi2'] for r in selected]==[5,30,35]

def test_fallback_seeds():
    rows=[row(5,0),row(6,.001),row(7,.002)]
    selected=NS['choose_seeds'](rows,3,np.array([-1.,-1.]),np.array([1.,1.]),.9)
    assert len(selected)==3
    assert selected[0]['chi2']==5

def test_radii_validation():
    assert NS['parse_radii']('0.02,0.05,0.10')==(.02,.05,.1)
    for bad in ('0','0.1,0.3','', '-0.1'):
        with pytest.raises(ValueError):NS['parse_radii'](bad)

def test_unbounded_cycles_default():
    main=next(n for n in TREE.body if isinstance(n,ast.FunctionDef) and n.name=='main')
    calls=[n for n in ast.walk(main) if isinstance(n,ast.Call) and isinstance(n.func,ast.Attribute)
           and n.func.attr=='add_argument']
    hit=[n for n in calls if n.args and isinstance(n.args[0],ast.Constant) and n.args[0].value=='--max-cycles']
    assert len(hit)==1 and any(k.arg=='default' and isinstance(k.value,ast.Constant) and k.value.value==0 for k in hit[0].keywords)

# Mass normalisation operates independently of the expensive physics evaluator.
MASS_NODES=[n for n in TREE.body if isinstance(n,ast.FunctionDef) and n.name in {'mass_scale_factor','normalize_mass_candidate'}]
MN={'np':np,'bulk':SimpleNamespace(NAMES=('y1_00','y2_00')),
    'eligible_candidate':lambda x:isinstance(x,dict) and x.get('chi2') is not None}
exec(compile(ast.Module(body=MASS_NODES,type_ignores=[]),str(SOURCE),'exec'),MN)

def mass_row(x=1., gap=.002511 * 1e16):
    return {'chi2':2e18,'parameters':{'y1_00':x,'y2_00':x},
            'prediction':{'delta_m3l_sq_ev2':gap}}

def test_mass_scaling_fourth_root():
    f=MN['mass_scale_factor']({'delta_m3l_sq_ev2':.002511*1e16},.002511)
    assert f == pytest.approx(1e-4)

def test_mass_scaling_wrong_sign_and_zero():
    for gap in (0,-.003,float('nan')):
        assert MN['mass_scale_factor']({'delta_m3l_sq_ev2':gap},.002511) is None

def test_mass_normalization_checks_authoritative_evaluator():
    start=mass_row()
    seen=[]
    def evaluate(x):
        seen.append(x.copy())
        factor=x[0]
        return {'chi2':2.,'parameters':{'y1_00':factor,'y2_00':factor},
                'prediction':{'delta_m3l_sq_ev2':.002511*1e16*factor**4}}
    result,status=MN['normalize_mass_candidate'](start,np.zeros(2),np.ones(2)*2,.002511,evaluate)
    assert status=='normalized'
    assert result['chi2']==2.
    assert len(seen)==1
    assert seen[0][0]==pytest.approx(1e-4)

def test_mass_normalization_rejects_mismatched_scaling():
    start=mass_row()
    def evaluate(x):
        return {'chi2':5.,'parameters':{'y1_00':float(x[0]),'y2_00':float(x[1])},
                'prediction':{'delta_m3l_sq_ev2':.009}}
    result,status=MN['normalize_mass_candidate'](start,np.zeros(2),np.ones(2)*2,.002511,evaluate)
    assert result is start and status=='scaling_mismatch'

def test_mass_normalization_preserves_bounds():
    start=mass_row(x=1e-6)
    result,status=MN['normalize_mass_candidate'](start,np.ones(2)*.01,np.ones(2)*2,.002511,lambda x:None)
    assert result is start and status=='outside_search_box'
