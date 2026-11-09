"""Fast numerical tests; does not require Matchete, RGBeta or Wolfram."""
import unittest
from pathlib import Path
from tempfile import TemporaryDirectory
from unittest.mock import patch
import json
import numpy as np

from Numerical.basis.T3Batch1TensorMap import (
    numeric_component, Polynomial, solve_dictionary, scan_hermitian_phases, canonical_species, TensorMapError
)
from scripts.RunT3Batch1GeneralContract import (
    MODELS, CATALOG, MIXED, SELF, native_commands, wolfram_integer_token,
    repair_rgbeta, inspect_beta
)


class TensorTests(unittest.TestCase):
    def test_models(self):
        self.assertEqual(len(MODELS),16)
        self.assertEqual(len(MIXED),2)
        self.assertEqual(len(SELF),6)
        self.assertEqual(sum(len(x[1]) for x in CATALOG.values()),16)

    def test_actual_repo_species_labels(self):
        self.assertEqual(canonical_species('Hdag'),('H',True))
        self.assertEqual(canonical_species('S2dag'),('S2',True))
        self.assertEqual(canonical_species('S1'),('S1',False))

    def test_safe_mathematica_numbers(self):
        self.assertAlmostEqual(numeric_component('1/Sqrt[2]').real,2**-.5)
        self.assertAlmostEqual(numeric_component('Rational[-3,4]').real,-.75)
        self.assertEqual(numeric_component('2 I'),2j)
        with self.assertRaises(TensorMapError):
            numeric_component('__import__("os")')

    def _rec(self,symbol,coefficient='1',add_hc=False,sp=None):
        return {'symbol':symbol, 'species':sp or ['Bar[H]','H','Bar[S1]','S1'],
                'dimensions':[2,2,1,1],
                'components':[coefficient,'0','0',coefficient],
                'potential_prefactor':'1','add_hc':add_hc}

    def test_matrix_dictionary_known_factor(self):
        # V_M = (Hdag H)(Sdag S), V_R = 2 V_M, so r=m/2.
        ms=[self._rec('lambdaM')]
        rs=[self._rec('lambdaR','2')]
        m=solve_dictionary(ms,rs,samples=64,seed=31)
        self.assertEqual(m['status'],'NumericalPolynomialMapComplete')
        self.assertAlmostEqual(m['map_rgbeta_from_matchete'][0][0],0.5,places=10)
        self.assertLess(m['holdout_relative_residual'],1e-12)

    def test_antihermitian_portal_is_rescued_by_i(self):
        # Anti-Hermitian invariant: imaginary matrix on a two-component H.
        # 1/2 PlusHc[i*O] survives; without i it cancels identically.
        t=np.zeros((2,2,1,1),complex)
        t[0,1,0,0]=1
        t[1,0,0,0]=-1
        base={'symbol':'portal','species':['Bar[H]','H','Bar[S1]','S1'],
              'dimensions':[2,2,1,1],
              'components':t.ravel().tolist(),'potential_prefactor':'1/2',
              'add_hc':True}
        h={'H':np.array([1+1j,.3+2j]),'S1':np.array([2+1j])}
        old=Polynomial.from_record({**base,'phase':'1'}).evaluate(h)
        new=Polynomial.from_record({**base,'phase':'I'}).evaluate(h)
        self.assertAlmostEqual(old,0.0,places=12)
        self.assertGreater(abs(new),1)

    def test_quartic_phase_scan_detects_antihermitian(self):
        t=np.zeros((2,2,1,1),complex)
        t[0,1,0,0]=1
        t[1,0,0,0]=-1
        record={'symbol':'anti','species':['Bar[H]','H','Bar[S1]','S1'],
                'dimensions':[2,2,1,1],
                'components':t.ravel().tolist(),'potential_prefactor':'1',
                'add_hc':True}
        self.assertEqual(scan_hermitian_phases([record])[0]['classification'],
                         'AntiHermitianRequiresI')

    def test_wolframscript_safe_negative_hypercharge(self):
        self.assertEqual(wolfram_integer_token(-4), 'm4')
        self.assertEqual(wolfram_integer_token(-3), 'm3')
        self.assertEqual(wolfram_integer_token(-2), 'm2')
        self.assertEqual(wolfram_integer_token(-1), 'm1')
        self.assertEqual(wolfram_integer_token(0), '0')
        self.assertEqual(wolfram_integer_token(2), '2')
        with TemporaryDirectory() as tmp:
            for model in MODELS:
                commands=native_commands(model,Path(tmp),0)
                self.assertEqual(len(commands), 2 if model in MIXED else 3)
                for label, command, _ in commands[:2]:
                    self.assertEqual(command[7],wolfram_integer_token(int(model.split('_')[1])))
                    self.assertFalse(command[7].startswith('-'))
                    self.assertEqual(command[8],model)

    def test_repair_only_invalid_and_missing_rgbeta(self):
        with TemporaryDirectory() as tmp:
            root=Path(tmp)
            tag='A_m4'
            folder=root/'rgbeta'/tag
            folder.mkdir(parents=True)
            good={'status':'SuccessCandidateMassBeta',
                  'metadata':{'GeneralScalarModel':'A_-4','GeneralScalarStage':'UV'},
                  'report_betas':{'g2':'some beta'}}
            (folder/'uv.json').write_text(json.dumps(good),encoding='utf-8')
            before=(folder/'uv.json').read_bytes()
            self.assertTrue(inspect_beta('A_-4',root)['uv']['candidate_ready'])
            with patch('scripts.RunT3Batch1GeneralContract.run_local',
                       return_value={'status':'Completed','returncode':0}) as run:
                records=repair_rgbeta('A_-4',root,600)
            self.assertEqual([r['status'] for r in records],
                             ['ReusedExistingCandidate','Completed'])
            self.assertEqual(run.call_count,1)
            self.assertEqual(run.call_args.args[0],'rgbeta_eft1')
            self.assertEqual(run.call_args.args[1][7],'m4')
            self.assertEqual((folder/'uv.json').read_bytes(),before)

    def test_repair_preserves_bad_log(self):
        with TemporaryDirectory() as tmp:
            root=Path(tmp)
            log=root/'logs'/'B_m3'/'rgbeta_uv.log'
            log.parent.mkdir(parents=True)
            log.write_text('Quadratic RGBeta error: Usage...',encoding='utf-8')
            with patch('scripts.RunT3Batch1GeneralContract.run_local',
                       return_value={'status':'Completed','returncode':0}):
                records=repair_rgbeta('B_-3',root,600)
            self.assertEqual(len(records),2)
            archive=Path(records[0]['archived_old_log'])
            self.assertTrue(archive.is_file())
            self.assertIn('Usage',archive.read_text(encoding='utf-8'))

    def test_incomplete_basis_blocks(self):
        ms=[self._rec('m')]
        rs=[self._rec('r','0')]
        with self.assertRaises(TensorMapError):
            solve_dictionary(ms,rs,samples=32)

    def test_bad_complex_potential_blocks(self):
        p=self._rec('unphysical', 'I', add_hc=False)
        with self.assertRaises(TensorMapError):
            Polynomial.from_record(p).evaluate({'H':np.array([1,0]),'S1':np.array([1])})

if __name__=='__main__':
    unittest.main()
