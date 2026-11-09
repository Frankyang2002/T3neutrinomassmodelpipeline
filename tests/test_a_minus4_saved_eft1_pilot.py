"""Self-contained tests of native saved-tensor integrity and provenance checks."""
from __future__ import annotations

from hashlib import sha256
import json
from pathlib import Path
import tempfile
import unittest

import numpy as np

from scripts.RunAminus4SavedEFT1Pilot import (
    SavedPilotError, decode_complex, repackage_saved_source, same_numbers,
)


class SavedAminus4PilotTests(unittest.TestCase):
    def setUp(self):
        self.t = tempfile.TemporaryDirectory()
        self.addCleanup(self.t.cleanup)
        self.root = Path(self.t.name)
        self.raw = np.zeros((3,3,2,2,8,8), dtype=np.complex128)
        self.raw[0,0,0,0,0,0] = 1e-11
        self.raw[0,1,0,1,0,4] = 0.2e-11 + 0.1e-11j
        self.raw[1,0,1,0,0,4] = 0.2e-11 + 0.1e-11j
        self.raw[0,1,0,1,4,0] = 0.2e-11 + 0.1e-11j
        self.raw[1,0,1,0,4,0] = 0.2e-11 + 0.1e-11j
        self.seed = self.root/'seed.json'
        self.seed.write_text('{"fake":true}', encoding='utf-8')
        self.sig = sha256(self.seed.read_bytes()+self.raw.tobytes()).hexdigest()
        self.saved = self.root/'dynamic_matching.npz'
        np.savez_compressed(self.saved,model=np.asarray('A_-4'),
                            C_LLSS_muF=self.raw,
                            C_LLSS_by_source_muF=self.raw[None,...],
                            muF_gev=np.asarray(1e5),muS_gev=np.asarray(1e3),
                            C5_direct_log=np.zeros((3,3)))
        self.result = {
            'status':'ConditionalDynamicFToBothScalarsToSMEFTCompleted',
            'model':'A_-4', 'scales_gev':{'F':1e5,'S':1e3},
            'flags':{'UV_evolved_fresh':True,
                     'all_heavy_fermions_integrated_out_at_F':True,
                     'LLSS_F_boundary_evaluated_from_current_Yukawas':True},
            'fermion_tree':{'shape':list(self.raw.shape),'source_count':1,
                             'norm_gev_inverse':float(np.linalg.norm(self.raw)),
                             'source_hash':self.sig,
                             'source_geometry':{'sources':[{'source_index':2}]}}
        }

    def test_repackage_preserves_formal_source_and_verifies_hash(self):
        out = self.root/'formal.npz'
        info = repackage_saved_source(self.saved,self.result,self.seed,out)
        self.assertEqual(info['source_hash'],self.sig)
        with np.load(out) as d:
            np.testing.assert_array_equal(d['C_LLSS_formal_muF'],self.raw)
            np.testing.assert_array_equal(d['matchete_PL_source_indices'],[2])
        meta=json.loads(out.with_suffix('.json').read_text())
        self.assertEqual(meta['format'],'t3-joint-scalar-source-bridge-v1')
        self.assertFalse(meta['gamma_cc_weyl_normalization_certified'])

    def test_bad_hash_refuses_different_matchete_source(self):
        self.seed.write_text('different source',encoding='utf-8')
        with self.assertRaisesRegex(SavedPilotError,'hash'):
            repackage_saved_source(self.saved,self.result,self.seed,self.root/'formal.npz')

    def test_bad_matching_tensor_refuses(self):
        wrong=self.raw.copy()
        wrong[1,0,1,0,4,0] += 1e-9
        np.savez_compressed(self.saved,model=np.asarray('A_-4'),
                            C_LLSS_muF=wrong,
                            C_LLSS_by_source_muF=wrong[None,...],
                            muF_gev=np.asarray(1e5),muS_gev=np.asarray(1e3),
                            C5_direct_log=np.zeros((3,3)))
        with self.assertRaises(SavedPilotError):
            repackage_saved_source(self.saved,self.result,self.seed,self.root/'formal.npz')

    def test_matrix_nested_complex_entries_preserve_phase(self):
        inp=[[{'real':1.,'imag':2.}, {'real':3.,'imag':0.}],
             [{'real':0.,'imag':4.}, {'real':2.,'imag':-1.}]]
        np.testing.assert_array_equal(decode_complex(inp),
                                      [[1+2j,3],[4j,2-1j]])

    def test_endpoint_comparison_rejects_cross_benchmark(self):
        self.assertLess(same_numbers(.6,.60000001,rtol=1e-4,atol=1e-10,
                                     name='g2'),1e-6)
        with self.assertRaisesRegex(SavedPilotError,'inconsistent'):
            same_numbers(.6,.65,rtol=1e-4,atol=1e-10,name='g2')


if __name__=='__main__':
    unittest.main()
