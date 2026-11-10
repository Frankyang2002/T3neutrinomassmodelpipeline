"""Small interface smoke tests; do not claim to execute Wolfram or physical matching.

Run: python -m unittest discover -s tests -p test_t3_thesis_general_batch1.py -v
"""
from __future__ import annotations

from hashlib import sha256
import json
from pathlib import Path
import tempfile
import unittest

import numpy as np

from Numerical.parallel_t3.background.general_batch1_api import (
    FORMAT, NativeBatch1Trajectory, decode_native, encode_native,
)
from scripts.RunT3ThesisGeneralBatch1 import (read_benchmark, real_yukawas,
                                              source_inventory, thesis_overrides)


class TestNativeStateInterface(unittest.TestCase):
    def construct(self, root: Path, *, model='D_-2', shifted=False) -> Path:
        folder = root / model
        folder.mkdir()
        f = np.array([1e7, 1e6, 1e5])
        s = np.array([1e5, 1e4, 1e3])
        uv_names = ['gY', 'lambdaT3', 'mS1Sq', 'mS2Sq', 'y1', 'y2', 'MF']
        eft_names = ['gY', 'lambdaT3', 'mS1Sq', 'mS2Sq']
        if model == 'D_-2':
            extras = ['muT3GeneralS2Sq'] + ['lambdaT3GeneralQ'+str(i) for i in range(5)]
            uv_names += extras + ['y2Conj']
            eft_names += extras
        def values(stage, key):
            t = f if stage == 'uv' else s
            if key in ['MF', 'y1', 'y2', 'y2Conj']:
                return np.stack([np.eye(3, dtype=complex) * (i+1) for i in range(len(t))])
            if key == 'gY':
                x = np.array([0.30, 0.35, 0.4]) if stage == 'uv' else np.array([0.4,0.45,0.5])
                if shifted and stage == 'eft1': x += .01
                return x.astype(complex)
            return np.ones(len(t), dtype=complex)
        arrays = {'uv__mu_gev': f, 'eft1__mu_gev': s}
        for n in uv_names: arrays['uv__'+n] = values('uv', n)
        for n in eft_names: arrays['eft1__'+n] = values('eft1', n)
        # Enforce exact continuity for all common keys in this synthetic dataset.
        for n in eft_names:
            if n != 'gY':
                arrays['eft1__'+n][0] = arrays['uv__'+n][-1]
        p = folder / 'trajectories.npz'
        np.savez_compressed(p, **arrays)
        fin = folder/'fermion_matching_input.json'
        sin = folder/'scalar_matching_input.json'
        fin.write_text('{"model":"' + model + '"}')
        sin.write_text('{"model":"' + model + '"}')
        def digest(path): return sha256(path.read_bytes()).hexdigest()
        m = {
            'format': FORMAT, 'status': 'CandidateGeneralRenormalisableRGIntegrated',
            'model': model,
            'scales_gev': {'UV':1e7,'F':1e5,'S':1e3,'low':100},
            'coupling_order': {'uv': uv_names, 'eft1': eft_names},
            'coupling_shapes': {stage: {n:list(arrays[stage+'__'+n].shape[1:])
                                        for n in names}
                                for stage, names in [('uv',uv_names),('eft1',eft_names)]},
            'eft1_boundary_deltas': {'gY': .01} if shifted else {},
            'trajectory_sha256': digest(p),
            'fermion_matching_input_sha256': digest(fin),
            'scalar_matching_input_sha256': digest(sin),
        }
        (folder/'manifest.json').write_text(json.dumps(m))
        return folder

    def test_round_trip_complex_mass_and_matrix(self):
        src = np.array([[1+2j, 3-4j], [0.1j, 5]])
        self.assertTrue(np.array_equal(src, decode_native(encode_native(src))))

    def test_complete_D2_matrix_and_bilinear_endpoints(self):
        with tempfile.TemporaryDirectory() as td:
            f = self.construct(Path(td))
            traj = NativeBatch1Trajectory.load(f)
            self.assertIn('y2Conj', traj.at_fermion_threshold())
            self.assertIn('muT3GeneralS2Sq', traj.at_scalar_threshold())
            self.assertEqual(len(traj.at_scalar_threshold()), 10)
            self.assertAlmostEqual(traj.at('uv', 1e6)['gY'].real, .35)
            self.assertAlmostEqual(traj.at('eft1', 1e4)['gY'].real, .45)
            self.assertEqual(traj.fermion_matching_input()['model'], 'D_-2')

    def test_explicit_boundary_shift_tracked(self):
        with tempfile.TemporaryDirectory() as td:
            f = self.construct(Path(td), shifted=True)
            t = NativeBatch1Trajectory.load(f)
            self.assertAlmostEqual(t.at('eft1', 1e5)['gY'].real, .41)

    def test_checksum_tampering_rejected(self):
        with tempfile.TemporaryDirectory() as td:
            f = self.construct(Path(td), model='A_-4')
            path = f/'fermion_matching_input.json'
            path.write_text('not original')
            with self.assertRaisesRegex(ValueError, 'checksum'):
                NativeBatch1Trajectory.load(f).fermion_matching_input()

    def test_inventory_fail_closed_on_missing_conjugate_yukawa(self):
        registry = {'models': {'D_-2': {
            'total_additional_quartics': 5,
            'new_rgbeta_symbols': [{'name': 'lambdaT3GeneralQ'+str(i)} for i in range(5)]}}}
        base = {'y1','y2','MF','mS1Sq','mS2Sq','lambdaT3'}
        additions = {'lambdaT3GeneralQ'+str(i) for i in range(5)}
        uv = base | additions | {'muT3GeneralS2Sq'}
        eft = uv - {'y1','y2','MF'}
        def payload(stage, names):
            return {'format':'general-t3-quadratic-rgbeta-experimental-v1',
                    'status':'SuccessCandidateMassBeta',
                    'metadata':{
                        'GeneralScalarModel':'D_-2','GeneralScalarStage':stage,
                        'ReportBetaConvention':'16*pi^2*dX/dln(mu); gauge BetaTerm divided by 2*g',
                        'GeneralScalarQuarticBetaBasis':'RGBetaFullQuarticProjection',
                        'GeneralScalarAdditionalQuartics':5,
                        'GeneralScalarQuadraticMassMixing':'IncludedCandidate',
                        'PhysicalMatcheteRGBetaMassBasisCertified':False,
                        'GeneralScalarNewBetaKeys':sorted(additions),
                        'QuadraticMassSymbol':'muT3GeneralS2Sq'},
                    'beta_symbols':sorted(names),'report_betas':{n:'1' for n in names}}
        with self.assertRaisesRegex(ValueError, 'conjugate Yukawa'):
            source_inventory(payload('UV',uv),payload('EFT1',eft),registry,'D_-2')


if __name__ == '__main__':
    unittest.main()
