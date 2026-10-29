"""Offline tests of the new persistent search's pure scheduling helpers."""
import ast
from pathlib import Path
from types import SimpleNamespace
import numpy as np
import pytest
from scipy.stats import qmc

SOURCE = Path(__file__).resolve().parents[1] / 'scripts' / 'SearchT3FrozenPersistent.py'
TREE = ast.parse(SOURCE.read_text(encoding='utf-8'))
NODES = [node for node in TREE.body if isinstance(node, ast.FunctionDef) and node.name in {'sobol_batch','choose_seeds'}]
NS = {'np': np, 'qmc': qmc, 'math': __import__('math'),
      'bulk': SimpleNamespace(NAMES=('y1_00','y2_00')),
      'eligible_candidate': lambda x: isinstance(x,dict) and x.get('chi2') is not None}
exec(compile(ast.Module(body=NODES,type_ignores=[]),str(SOURCE),'exec'),NS)

def test_deterministic_sobol():
    lo=np.array([-1.,0.]);hi=np.array([1.,2.])
    a=NS['sobol_batch'](lo,hi,16,47)
    assert a.shape == (16,2)
    assert np.array_equal(a,NS['sobol_batch'](lo,hi,16,47))
    assert not np.array_equal(a,NS['sobol_batch'](lo,hi,16,48))
    assert np.all(a>=lo) and np.all(a<=hi)

def test_invalid_sobol_count():
    with pytest.raises(ValueError): NS['sobol_batch'](np.zeros(2),np.ones(2),15,4)

def test_seeds_sorted_and_deduplicated():
    def row(chi,x):return {'chi2':chi,'parameters':{'y1_00':x,'y2_00':0}}
    rows=[row(100,1),row(20,2),row(50,2),row(5,3)]
    assert [r['chi2'] for r in NS['choose_seeds'](rows,2)] == [5,20]

def test_loop_default_is_unbounded():
    main=next(n for n in TREE.body if isinstance(n,ast.FunctionDef) and n.name=='main')
    defaults=[n for n in ast.walk(main) if isinstance(n,ast.Call) and
              isinstance(n.func,ast.Attribute) and n.func.attr=='add_argument']
    hit=[n for n in defaults if n.args and isinstance(n.args[0],ast.Constant) and n.args[0].value=='--max-cycles']
    assert len(hit)==1 and any(k.arg=='default' and isinstance(k.value,ast.Constant) and k.value.value==0 for k in hit[0].keywords)
