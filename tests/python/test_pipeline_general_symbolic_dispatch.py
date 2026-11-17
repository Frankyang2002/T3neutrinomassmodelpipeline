"""Offline pipeline routing tests; no Matchete, RGBeta or physical matching."""
from __future__ import annotations

import argparse
import ast
from pathlib import Path
from types import SimpleNamespace, ModuleType
import sys

SOURCE = Path(__file__).resolve().parents[2] / 'pipeline.py'


def funcs():
    tree = ast.parse(SOURCE.read_text(encoding='utf8'))
    chosen = [node for node in tree.body if isinstance(node, ast.FunctionDef)
              and node.name in {'_use_general_symbolic', '_run_general_symbolic'}]
    code = compile(ast.Module(body=chosen, type_ignores=[]), str(SOURCE), 'exec')
    ns = {'argparse': argparse}
    exec(code, ns)
    return ns


def args(**updates):
    values = dict(
        general_symbolic=False, restricted_production=False,
        reports_only=False, hypercharge_comparison=True,
        numerical=None, full=False, dims=None, force=False,
        reset_numerical_configs=False, study=None,
        threshold=None, threshold_scale=None, debug_reports=False,
    )
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
    # The real pipeline delegates to both symbolic generation and PDF rendering.
    # Mock BOTH dependencies: the old test intercepted only the symbolic stage
    # and incorrectly invoked a real adapter with an incomplete fake module.
    symbolic = ModuleType('scripts.RunT3GeneralSymbolicStudy')
    reports = ModuleType('Reports.GeneralSymbolicAdapter')
    calls = []
    symbolic.main = lambda argv: calls.append(('symbolic', argv)) or 0
    reports.render = lambda *, study: calls.append(('render', study)) or []
    monkeypatch.setitem(sys.modules, symbolic.__name__, symbolic)
    monkeypatch.setitem(sys.modules, reports.__name__, reports)
    ns = funcs()
    ns['_run_general_symbolic'].__globals__['__builtins__'] = __builtins__
    assert ns['_run_general_symbolic'](
        args(study='my_study', reports_only=True)) == 0
    assert calls == [
        ('symbolic', ['--study', 'my_study', '--reports-only']),
        ('render', 'my_study'),
    ]
