#!/usr/bin/env python3
"""Static/native-WXF evidence tests for independent-mass EFT1 Wilson tensors.

No Wolfram or Matchete execution is claimed by this test. The actual live
matching probe is RunT3EFT1WilsonTensorBatch.py --pilot --mode soft.
"""
from __future__ import annotations

from collections import Counter
from pathlib import Path
import sys

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
from RunT3EFT1WilsonTensorBatch import (  # noqa: E402
    Blocked, check_definitions, check_native, preflight, MODELS,
)


def native_report(model: str, mode: str, variant: str) -> dict:
    return {
        'format': 't3-eft1-flavour-wilson-soft-v1',
        'model': model, 'mode': mode, 'variant': variant,
        'status': ('FermionThresholdWilsonTreeContracted' if mode == 'prepare'
                   else 'ScalarThresholdWilsonSoftMatchedNotRenormalized'),
        'heavy_generations_integrated_together': 3,
        'independent_heavy_masses_preserved_in_wilson_definition': True,
        'tree_level_eft1_wilson_contraction_applied': True,
        'stage1_hard_one_loop_inherited_at_stage2': False,
        'eft1_wilson_running_included': False,
        'renormalized_hierarchical_c5_certified': False,
        'absolute_matchete_rgbeta_dictionary_certified': False,
        'physical_neutrino_prediction_certified': False,
        'scalar_threshold_native_passed': mode != 'prepare',
        'stage2_loop_order': 0 if mode == 'tree' else 1,
        'eft1_wilson_running_insertion_added': False,
        'wilson_counts': {'baseline': 2, 'extra': 2 if variant == 'seeded' else 0,
                          'unmatched_mass': 0},
    }


def test_schemas() -> None:
    for model in MODELS:
        for mode in ('prepare', 'tree', 'soft'):
            for variant in ('seeded', 'zero'):
                original = native_report(model, mode, variant)
                check_native(original, model, variant, mode)
                for field in ('stage1_hard_one_loop_inherited_at_stage2',
                              'eft1_wilson_running_included',
                              'renormalized_hierarchical_c5_certified',
                              'physical_neutrino_prediction_certified'):
                    changed = dict(original, **{field: True})
                    try:
                        check_native(changed, model, variant, mode)
                    except Blocked:
                        pass
                    else:
                        raise AssertionError('Bad certification accepted: ' + field)
                changed = dict(original, wilson_counts={
                    'baseline': 1, 'extra': 0, 'unmatched_mass': 1})
                try:
                    check_native(changed, model, variant, mode)
                except Blocked:
                    pass
                else:
                    raise AssertionError('Invalid tensor count accepted')
                extra = ('Sum[y1Conj[p,r] y2[q,r] / MF[r], {r,1,3}]'
                         if model.endswith('p0') else
                         'Sum[y1[p,r] y2Conj[q,r] / MF[r], {r,1,3}]')
                definitions = {
                    'T3EFT1W12[p,q]': 'Sum[y1[p,r] y2[q,r] / MF[r], {r,1,3}]',
                    'T3EFT1WExtra[p,q]': extra,
                    'mass_assumption': 'Three diagonal independent masses, integrated at one common F threshold'
                }
                check_definitions(definitions, variant, model)
                definitions['T3EFT1WExtra[p,q]'] = 'wrong conjugate channel'
                try:
                    check_definitions(definitions, variant, model)
                except Blocked:
                    pass
                else:
                    raise AssertionError('Wrong EFT1 flavour tensor accepted')
    print('PASS: 36 six-model/variant/mode native schema checks and negative controls')


def test_numeric_three_mass_tensor() -> None:
    rng = np.random.default_rng(20261101)
    y1 = rng.normal(size=(3, 3)) + 1j * rng.normal(size=(3, 3))
    y2 = rng.normal(size=(3, 3)) + 1j * rng.normal(size=(3, 3))
    y1c = rng.normal(size=(3, 3)) + 1j * rng.normal(size=(3, 3))
    y2c = rng.normal(size=(3, 3)) + 1j * rng.normal(size=(3, 3))
    masses = np.array([9.2e4, 1.15e5, 1.37e5])
    w = (y1 / masses) @ y2.T
    wexp = sum(np.outer(y1[:, r], y2[:, r]) / masses[r] for r in range(3))
    assert np.allclose(w, wexp, rtol=1e-14, atol=0)
    wx1 = (y1c / masses) @ y2.T
    wx2 = (y1 / masses) @ y2c.T
    assert np.linalg.norm(wx1 - wx2) > 1e-8
    assert np.allclose(np.conjugate(w),
                       (np.conjugate(y1) / masses) @ np.conjugate(y2).T)
    # A common matching scale does not replace the three mass denominators
    # with a common mass; doing so is observably a different approximation.
    w_degenerate = (y1 / masses.mean()) @ y2.T
    assert np.linalg.norm(w - w_degenerate) > 1e-7
    print('PASS: 3 independent MF[r], matrix Wilson contractions, conjugate channels')


def test_uploaded_native_tree_witness() -> None:
    # Earlier diagnostic script was installed in the user's previous batch.
    # The attached WXF is authoritative for A_p0 only; other models must be
    # validated by actual live matching, not by reusing this A_p0 input.
    path = (ROOT / 'Reports/output/workstream_a/extended_c5_hierarchy/hierarchy/'
            'A_p0/seeded/debug/fresh_kernel_stage_2/input_transition_tree.wxf')
    if not path.is_file():
        print('SKIP: archived A_p0 stage1 WXF unavailable; live probe will rebuild it')
        return
    try:
        from DiagnoseT3Stage2Index import WXFReader, Function, symbol_name, walk
    except ImportError:
        print('SKIP: previous WXF diagnostic script unavailable')
        return
    root = WXFReader(path.read_bytes()).read()
    assert isinstance(root, Function) and symbol_name(root.head) == 'Plus'
    counts = Counter()
    for term in root.args:
        if not isinstance(term, Function) or symbol_name(term.head) != 'Times':
            continue
        factors = term.args
        mf = [x for x in factors if isinstance(x, Function)
              and symbol_name(x.head) == 'Power' and len(x.args) == 2
              and x.args[1] == -1 and isinstance(x.args[0], Function)
              and symbol_name(x.args[0].head) == 'Coupling'
              and symbol_name(x.args[0].args[0]) == 'MF'
              ]
        # The standard WXF printer does not turn the 3-index inverse mass
        # into a valid Einstein product; this audit therefore counts only
        # the structural inverse-MF factors and their Yukawas.
        if not mf:
            continue
        names = [symbol_name(node.args[0]) for node in walk(term)
                 if isinstance(node, Function) and symbol_name(node.head) == 'Coupling'
                 and node.args]
        assert 'y2' in names and any(s in names for s in ('y1', 'y1Conj'))
        counts['extra' if 'y1Conj' in names else 'baseline'] += 1
    assert counts == {'baseline': 2, 'extra': 2}, counts
    print('PASS: uploaded A_p0 stage1 WXF has precisely 2 baseline and 2 extra MF[r] factors')


def main() -> int:
    assert preflight(ROOT, verify_real=False)['status'] == 'EFT1WilsonTensorSourceGatePassed'
    test_schemas()
    test_numeric_three_mass_tensor()
    test_uploaded_native_tree_witness()
    print('SYNTHETIC + SOURCE CHECK ONLY: no live scalar Matchete result yet')
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
