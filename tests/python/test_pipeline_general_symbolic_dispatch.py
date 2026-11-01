"""Offline routing tests: no Matchete, RGBeta, or physics assertions."""
import argparse
import ast
from pathlib import Path
from types import SimpleNamespace

SOURCE = Path(__file__).resolve().parents[2] / 'pipeline.py'


def funcs():
    tree = ast.parse(SOURCE.read_text(encoding='utf8'))
    chosen = [n for n in tree.body if isinstance(n, ast.FunctionDef)
              and n.name in {'_use_general_symbolic', '_run_general_symbolic'}]
    code = compile(ast.Module(body=chosen, type_ignores=[]), str(SOURCE), 'exec')
    ns = {'argparse': argparse}
    exec(code, ns)
    return ns


def args(**updates):
    values = dict(general_symbolic=False, restricted_production=False,
                  reports_only=False, hypercharge_comparison=True,
                  numerical=None, full=False, dims=None, force=False,
                  reset_numerical_configs=False, study=None,
                  threshold=None, threshold_scale=None, debug_reports=False)
    values.update(updates)
    return argparse.Namespace(**values)


def test_16_model_default_routes_general():
    assert funcs()['_use_general_symbolic'](args(), argparse.ArgumentParser())


def test_explicit_legacy_preserved():
    assert not funcs()['_use_general_symbolic'](
        args(restricted_production=True), argparse.ArgumentParser())


def test_numerical_uses_legacy():
    assert not funcs()['_use_general_symbolic'](
        args(numerical=Path('benchmark.json')), argparse.ArgumentParser())


def test_reports_only_routes_general():
    assert funcs()['_use_general_symbolic'](
        args(reports_only=True), argparse.ArgumentParser())


def test_bad_combination_fails():
    try:
        funcs()['_use_general_symbolic'](
            args(general_symbolic=True, restricted_production=True),
            argparse.ArgumentParser())
    except SystemExit as exc:
        assert exc.code == 2
    else:
        raise AssertionError('Conflicting modes were accepted')


def test_delegation_uses_same_study_name(monkeypatch):
    import sys
    import types
    fake = types.ModuleType('scripts.RunT3GeneralSymbolicStudy')
    calls = []
    fake.main = lambda argv: calls.append(argv) or 0
    monkeypatch.setitem(sys.modules, 'scripts.RunT3GeneralSymbolicStudy', fake)
    ns = funcs()
    ns['_run_general_symbolic'].__globals__['__builtins__'] = __builtins__
    assert ns['_run_general_symbolic'](args(study='my_study', reports_only=True)) == 0
    assert calls == [['--study', 'my_study', '--reports-only']]
