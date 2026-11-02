"""Cross-check imported READ-ONLY existing GitHub provisional reference snapshot.

This is not a physical-certification test: it checks numerical transcription and
source convention statuses, and deliberately refuses to promote the baseline.
"""
import json
from pathlib import Path
import unittest
import numpy as np

ROOT=Path(__file__).resolve().parents[1]
BASE=ROOT/'output/general_observable/lo_certification_direct_uv'

class ExistingSourceSnapshot(unittest.TestCase):
    def test_all16_reference_matrices_symmetric_and_takagi(self):
        data=json.loads((BASE/'github_provisional_matrix_reference.json').read_text())
        self.assertEqual(data['source_physical_certification_count'],0)
        rows=data['models'];self.assertEqual(len(rows),16)
        self.assertEqual(len({x['model'] for x in rows}),16)
        for row in rows:
            with self.subTest(model=row['model']):
                z=row['mnu_3x3_gev']
                m=np.asarray(z['real'])+1j*np.asarray(z['imag'])
                self.assertEqual(m.shape,(3,3))
                self.assertLess(np.linalg.norm(m-m.T),1e-24)
                masses=np.sort(np.linalg.svd(m,compute_uv=False))*1e9
                self.assertTrue(np.allclose(masses,row['takagi_masses_ev'],rtol=2e-12,atol=2e-13))
                self.assertIn('provisional',row['status'])

    def test_source_coverage_six_extra_yukawa_models(self):
        data=json.loads((BASE/'source_coverage_all16.json').read_text())
        self.assertEqual(data['algebraic_audit_passed'],16)
        self.assertEqual(data['physical_source_certified'],0)
        self.assertEqual(data['geometric_native_portal_invariants_representable'],32)
        self.assertEqual(data['native_extra_yukawa_source_allowed_count'],6)
        rows=data['rows'];self.assertEqual(len(rows),16)
        extra=[r['model'] for r in rows if r['additional_hypercharge_zero_scalar_yukawa_allowed']]
        self.assertEqual(extra,['A_-2','A_+0','D_-2','D_+0','E_-2','E_+0'])
        self.assertEqual(sum(r['number_of_common_charge_channels'] for r in rows),20)
        self.assertTrue(all(r['physical_source_mass_matrix_certified'] is False for r in rows))

if __name__=='__main__':unittest.main()
