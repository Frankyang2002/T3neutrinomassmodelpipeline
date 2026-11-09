"""Unit tests for the explicit Matchete polynomial -> Weyl coordinate bridge.

These tests verify geometry, source symmetries, no implicit combinatorial factors,
and the fail-closed source metadata contract. They do NOT certify GammaCC.
"""
from __future__ import annotations

import json
from pathlib import Path
import sys
import tempfile
import types
import unittest
from unittest.mock import patch

import numpy as np

from Numerical.running.T3EFT1MasterSourceBridge import (
    SourceBridgeError, embed_formal_source, load_formal_source,
)


def physical_sm_fermion_labels():
    labels = []
    for species in ('Q', 'uc', 'dc', 'L', 'ec'):
        for g in (1, 2, 3):
            if species == 'Q':
                labels.extend(f'Q[{g},{col},{w}]' for col in (1,2,3)
                              for w in (1,2))
            elif species in ('uc', 'dc'):
                labels.extend(f'{species}[{g},{col}]' for col in (1,2,3))
            elif species == 'L':
                labels.extend(f'L[{g},{w}]' for w in (1,2))
            else:
                labels.append(f'ec[{g}]')
    return tuple(labels)


F = physical_sm_fermion_labels()
S = ('H[1,R]', 'H[1,I]', 'H[2,R]', 'H[2,I]') + tuple(
    f'S[{i}]' for i in range(8))


def fixture():
    a = np.zeros((3, 3, 2, 2, 8, 8), dtype=complex)
    for p, q, i, j in ((0, 1, 0, 1), (1, 0, 1, 0)):
        a[p, q, i, j, 0, 1] = 1 + 2j
        a[p, q, i, j, 1, 0] = 1 + 2j
    return a


class BridgeTests(unittest.TestCase):
    def setUp(self):
        # Isolate source geometry tests from the rest of the (large) repository.
        stub = types.ModuleType('Numerical.running.GeneralT3OrderedTensor')
        stub.scalar_basis = lambda model: types.SimpleNamespace(total_dimension=12)
        self.mock = patch.dict(sys.modules, {'Numerical.running.GeneralT3OrderedTensor': stub})
        self.mock.start()
        self.addCleanup(self.mock.stop)

    def test_lepton_embedding_and_exact_polynomial_factor(self):
        a = fixture()
        c, meta = embed_formal_source(a, model='A_-4',
                                       fermion_labels=F, scalar_labels=S)
        self.assertEqual(c.shape, (45, 45, 12, 12))
        self.assertEqual(c[36, 39, 4, 5], 1 + 2j)
        self.assertEqual(c[39, 36, 5, 4], 1 + 2j)
        self.assertEqual(c[36, 39, 4, 4], 0)
        self.assertFalse(meta['claimed_physical_master_boundary'])
        self.assertEqual(meta['scalar_coordinate_map_status'],
                         'IDENTITY_ASSUMED_FORMAL_ONLY')
        self.assertTrue(np.array_equal(c, c.swapaxes(0, 1)))
        self.assertTrue(np.array_equal(c, c.swapaxes(2, 3)))

    def test_explicit_real_scalar_rotation(self):
        a = fixture()
        R = np.eye(8)
        R[0, 0] = R[1, 1] = 0
        R[0, 1], R[1, 0] = 1, -1
        c, meta = embed_formal_source(a, model='A_-4',
                                       fermion_labels=F, scalar_labels=S,
                                       heavy_source_from_target=R,
                                       spinor_signed_factor=-2.0)
        self.assertEqual(c[36, 39, 4, 5], 2 * (1 + 2j))
        self.assertEqual(meta['spinor_normalization_status'],
                         'EXTERNAL_FACTOR_UNCERTIFIED')

    def test_rejects_exchange_asymmetry(self):
        a = fixture()
        a[0, 1, 0, 1, 1, 0] = 0
        a[1, 0, 1, 0, 1, 0] = 0
        with self.assertRaisesRegex(SourceBridgeError, 'Scalar interchange'):
            embed_formal_source(a, model='A_-4', fermion_labels=F, scalar_labels=S)

    def test_rejects_bad_basis_or_map(self):
        with self.assertRaisesRegex(SourceBridgeError, 'Higgs'):
            embed_formal_source(fixture(), model='A_-4',
                                fermion_labels=F, scalar_labels=S[::-1])
        with self.assertRaisesRegex(SourceBridgeError, 'orthogonal'):
            embed_formal_source(fixture(), model='A_-4',
                                fermion_labels=F, scalar_labels=S,
                                heavy_source_from_target=np.ones((8, 8)))

    def test_exact_formal_export_read(self):
        with tempfile.TemporaryDirectory() as td:
            f = Path(td)/'source.npz'
            np.savez_compressed(f, C_LLSS_formal_muF=fixture())
            meta = {
                'status':'FormalLLSSExported_PhysicalNormalizationUncertified',
                'gamma_cc_weyl_normalization_certified':False,
                'model':'A_-4', 'mu_f_gev':1e5, 'n_heavy_real':8,
                'source_scalar_real_labels':['s']*8,
            }
            f.with_suffix('.json').write_text(json.dumps(meta),encoding='utf8')
            inp = load_formal_source(f, model='A_-4')
            self.assertTrue(np.array_equal(inp.raw, fixture()))
            self.assertEqual(inp.mu_f_gev, 1e5)
            meta['gamma_cc_weyl_normalization_certified'] = True
            f.with_suffix('.json').write_text(json.dumps(meta),encoding='utf8')
            with self.assertRaisesRegex(SourceBridgeError, 'normalization status'):
                load_formal_source(f)


if __name__ == '__main__':
    unittest.main()
