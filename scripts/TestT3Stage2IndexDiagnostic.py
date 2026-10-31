#!/usr/bin/env python3
"""Regression for uploaded stage-2 tree-index diagnosis and mass-mode gating."""
from __future__ import annotations

import importlib.util
from pathlib import Path
import sys

SCRIPTS = Path(__file__).resolve().parent


def load(name: str, file: str):
    spec = importlib.util.spec_from_file_location(name, SCRIPTS / file)
    if spec is None or spec.loader is None:
        raise RuntimeError('Cannot load module ' + file)
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


def main() -> int:
    audit = load('t3_diagnose_index', 'DiagnoseT3Stage2Index.py')
    runner = load('t3_hierarchy_runner', 'RunT3ExtendedHierarchyBatch.py')
    audit.self_test()
    runner.self_test()
    a = {
        'format': 't3-extra-yukawa-hierarchical-native-v1',
        'model': 'A_p0', 'mode': 'stage1', 'variant': 'seeded',
        'heavy_flavours': 3, 'mass_mode': 'common',
        'independent_heavy_masses_certified': False,
        'mass_mode_common_is_diagnostic_only': True,
        'status': runner.SUCCESS['stage1'], 'stage_count': 1,
        'fermion_threshold_integrated_out': True,
        'scalar_threshold_integrated_out': False,
        'eft1_Wilson_RGE_evolved': False,
        'absolute_matchete_rgbeta_dictionary_certified': False,
        'physical_low_energy_prediction_certified': False,
        'full_scalar_quadratic_basis_included': False,
    }
    runner.check_native(a, 'A_p0', 'stage1', 'seeded', 'common')
    for tweak in ({'mass_mode':'indexed'},
                  {'independent_heavy_masses_certified':True},
                  {'mass_mode_common_is_diagnostic_only':False},
                  {'eft1_Wilson_RGE_evolved':True}):
        try:
            runner.check_native({**a, **tweak}, 'A_p0', 'stage1', 'seeded', 'common')
        except runner.GateError:
            pass
        else:
            raise AssertionError('Unsafe common-mass report accepted: '+str(tweak))
    native = (SCRIPTS.parent/'Lagrangian/RunT3ExtendedHierarchyNative.wl').read_text(encoding='utf-8')
    assert 'hMassMode = If[Length[hArgs] === 6' in native
    assert '"mass_mode_common_is_diagnostic_only"->(hMassMode === "common")' in native
    assert 'RunSequentialT3Matching[uv,model,plan,5,1,hOut,hModel]' in native
    assert '"eft1_Wilson_RGE_evolved"->False' in native
    assert '"T3_EXPERIMENTAL_FLAVOR_MASSES"] =!=' in native
    path = audit.DEFAULT
    if path.is_file():
        audit.assert_a0(audit.audit_tree(path))
        print('PASS: real local A(0) 43-term/4-violation EFT1 WXF witness')
    else:
        print('SKIP: real local WXF witness is not present; synthetic contract tests passed')
    print('PASS: native source/common-mass provenance gates, negative controls')
    print('NOTE: no native Wolfram execution was performed by this test')
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
