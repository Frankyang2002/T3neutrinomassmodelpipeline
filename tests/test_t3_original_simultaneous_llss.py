"""Run with: python -m unittest tests/test_t3_original_simultaneous_llss.py"""
import unittest
import json
import tempfile
import types
import sys
from pathlib import Path
from unittest.mock import patch
from types import SimpleNamespace
from unittest.mock import Mock
import numpy as np

from scripts.RunT3OriginalSimultaneousLLSSHeatmap import (
    source_in_original_basis, LLSSIsolator, _comparison,
    _archived_symbolic_fingerprint, ALLOWED_STALE_SYMBOLIC_SOURCE,
)


class OriginalLLSSIsolationTest(unittest.TestCase):
    def setUp(self):
        self.ctx=SimpleNamespace(model='B_-1', nf=6, ns=8,
            higgs=(0,1,2,3), leptons=((0,1),(2,3),(4,5)),
            quartic=np.zeros((8,)*4))
        self.c=np.zeros((6,6,8,8),complex)
        self.c[0,0,4,5]=3
        self.c[0,0,5,4]=3
        self.ctx.quartic[2,2,4,5]=.5
        self.ctx.quartic[2,2,5,4]=.5

    def test_source_only_in_llhh(self):
        src=source_in_original_basis(self.ctx,self.c)
        self.assertAlmostEqual(src[0,0,2,2].real,3)
        self.assertEqual(src[0,0,4,5],0)
        self.assertEqual(src[2,2,2,2],0)

    def test_initial_parity_remains_unmodified(self):
        original=Mock(side_effect=lambda ctx,c:{'scalar_pair': np.ones_like(c)})
        iso=LLSSIsolator(original)
        with iso.integration(eft_on=True):
            pre=iso.beta(self.ctx,self.c)
            post=iso.beta(self.ctx,self.c)
        self.assertAlmostEqual(pre['scalar_pair'][0,0,2,2].real,1)
        self.assertAlmostEqual(post['scalar_pair'][0,0,2,2].real,-2)
        self.assertEqual(iso.per_model['B_-1']['rhs_calls'],1)
        self.assertGreater(iso.per_model['B_-1']['max_rhs_source_norm'],0)

    def test_inactive_passthrough(self):
        iso=LLSSIsolator(lambda ctx,c:{'scalar_pair':np.ones_like(c)})
        with iso.integration(eft_on=False):
            out=iso.beta(self.ctx,self.c)
        self.assertEqual(out['scalar_pair'][0,0,2,2],1)

    def test_comparison_same_frozen_denominator(self):
        fields=('delta_m21_sq_ev2','delta_m3l_sq_ev2','sin2_theta12','sin2_theta13','sin2_theta23')
        case=lambda v:{'prediction':dict.fromkeys(fields,v)}
        row=_comparison(case(10),case(2),case(3),'A_-4')
        self.assertEqual(row['delta_m21_sq_ev2_pct_running'],400)
        self.assertEqual(row['delta_m21_sq_ev2_pct_mixing_off'],50)

    def test_bad_denominator_is_blocked(self):
        fields=('delta_m21_sq_ev2','delta_m3l_sq_ev2','sin2_theta12','sin2_theta13','sin2_theta23')
        case=lambda v:{'prediction':dict.fromkeys(fields,v)}
        with self.assertRaises(ValueError):_comparison(case(1),case(0),case(2),'B_-1')

    def _check_fingerprint(self, changed, extra=None):
        with tempfile.TemporaryDirectory() as td:
            root=Path(td)
            general=types.ModuleType('scripts.RunT3GeneralAll16FromReference')
            general.OUTPUT=root/'fit'
            canonical=types.ModuleType('scripts.RunT3GeneralSymbolicStudy')
            canonical.DEFAULT_DATA_ROOT=root/'canonical'
            canonical.model_tag=lambda m:m.replace('-','m').replace('+','p')
            actual={ALLOWED_STALE_SYMBOLIC_SOURCE:'current', 'Lagrangian/model/T3Fields.wl':'same'}
            canonical.fingerprint=lambda:dict(actual)
            (general.OUTPUT/'thesis16_merged_v1').mkdir(parents=True)
            (general.OUTPUT/'thesis16_merged_v1'/'settings.json').write_text(
                json.dumps({'canonical_run_id':'saved_source'}))
            cached={ALLOWED_STALE_SYMBOLIC_SOURCE:'archived',
                    'Lagrangian/model/T3Fields.wl':changed}
            if extra:cached.update(extra)
            for model in ('A_-4','A_-2'):
                directory=canonical.DEFAULT_DATA_ROOT/'saved_source'/canonical.model_tag(model)
                directory.mkdir(parents=True)
                (directory/'provenance.json').write_text(json.dumps({
                    'model':model, 'source_sha256':cached}))
            with patch.dict(sys.modules, {
                general.__name__:general,canonical.__name__:canonical}):
                return _archived_symbolic_fingerprint(('A_-4','A_-2'),
                                                       'thesis16_merged_v1')

    def test_archived_source_allows_only_known_wolfram_file(self):
        cached,audit=self._check_fingerprint('same')
        self.assertEqual(cached[ALLOWED_STALE_SYMBOLIC_SOURCE],'archived')
        self.assertEqual(len(audit['models']),2)
        self.assertEqual(audit['models'][0]['source_mismatches'],
                         [ALLOWED_STALE_SYMBOLIC_SOURCE])

    def test_archived_source_rejects_other_source_changes(self):
        with self.assertRaisesRegex(ValueError,'unapproved symbolic-source changes'):
            self._check_fingerprint('different')

    def test_archived_source_rejects_source_inventory_drift(self):
        with self.assertRaisesRegex(ValueError,'source inventory differs'):
            self._check_fingerprint('same', {'Lagrangian/new_source.wl':'different'})


if __name__=='__main__':unittest.main()
