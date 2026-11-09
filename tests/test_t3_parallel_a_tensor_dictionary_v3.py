"""Workstream A conditional native-evidence extra-quartic mapping tests."""
from __future__ import annotations

import copy
from pathlib import Path
import math
import json
import subprocess
import sys
import tempfile
import unittest

import numpy as np

from Numerical.parallel_t3.wilson.matchete_rgbeta_dictionary import (
    DictionaryEvidenceError, _native_invariant, _mat_tensor,
    _target_tensor, calibrated_triplet_cg, matchete_number,
    solve_from_paths, solve_from_uploaded_evidence,
)

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT/'output/parallel_t3/wilson/matchete_native_tensors.json'
DST = ROOT/'output/parallel_t3/wilson/rgbeta_abstract_components_conditional_v2.json'


def cartoon_cg():
    """Cartesian CG with non-real phases, normalized on triplet axes."""
    g=np.zeros((2,2,3),complex)
    a=1/math.sqrt(2)
    g[0,0,0]=a; g[0,0,2]=-1j*a
    g[0,1,1]=g[1,0,1]=-a
    g[1,1,0]=-a;g[1,1,2]=-1j*a
    return g


class TestT3ParallelANativeDictionaryV3(unittest.TestCase):
    def test_parser_matchete_sqrt_quartic_and_phase(self):
        trials={
            '3^(-1/4)':3**(-.25),'(-I)/3^(1/4)':-1j/3**.25,
            '1/Sqrt[2]':1/math.sqrt(2),'-(1/Sqrt[3])':-1/math.sqrt(3),
            '-1/2*Sqrt[3/11]':-0.5*math.sqrt(3/11),
            '2*Sqrt[3/11]':2*math.sqrt(3/11),
        }
        for expr,want in trials.items():
            with self.subTest(expr=expr):
                self.assertAlmostEqual(abs(matchete_number(expr)-want),0,places=12)

    def test_reject_unsafe_numeric_or_unsupported_native_AST(self):
        for expr in ('__import__("os").system("exit 4")', 'x',
                     'Sqrt[2,3]', '[1,2]', 'I.real'):
            with self.subTest(expr=expr), self.assertRaises(DictionaryEvidenceError):
                matchete_number(expr)
        g=cartoon_cg()
        for native in ('__import__("os").system("exit 4") &',
                       'tGen[SU2L, a, #1, #2] &', 'foo(#1) &'):
            with self.subTest(native=native), self.assertRaises(DictionaryEvidenceError):
                _native_invariant(native,(1,1,1,1),('Bar[H]',)*4,g)

    def test_cartesian_triplet_cg_orthonormal(self):
        g=cartoon_cg()
        self.assertLess(np.max(abs(np.einsum('ija,ijb->ab',g.conj(),g)-np.eye(3))),1e-14)
        self.assertLess(np.max(abs(g-g.swapaxes(0,1))),1e-14)

    def test_unbarred_doublet_CG_conjugation_is_essential(self):
        g=cartoon_cg()
        species=('Bar[H]','Bar[H]','S2','S2')
        expr='delS2[SU2L, a, #1, #2]*delS2[SU2L, a, #3, #4] &'
        actual=np.array([_native_invariant(expr,tuple(idx),species,g)
                         for idx in __import__('itertools').product((1,2),repeat=4)]).reshape((2,)*4)
        expected=np.array([0.5*(int(i==k and j==l)+int(i==l and j==k))
                           for i,j,k,l in __import__('itertools').product((1,2),repeat=4)]).reshape((2,)*4)
        self.assertLess(np.max(abs(actual-expected)),1e-14)
        # Wrong (unstarred) Cartesian contraction has an O(1) error.
        wrong=np.einsum('ija,kla->ijkl',g,g)
        self.assertGreater(np.max(abs(wrong-expected)),.9)

    def test_triplet_anchor_normalization_synthetic_record(self):
        g=cartoon_cg();scale=math.sqrt(2/math.sqrt(3))
        raw=g*scale
        rec={'symbol':'lambdaT3ExtraQ021010I1',
             'species':['Hdag','Hdag','S1','S2'],
             'dimensions':[2,2,1,3],
             'potential_prefactor':'1',
             'components':[f'({complex(x).real}+({complex(x).imag})*I)' for x in raw.reshape(-1)]}
        m={'models':{'A_+0':{'tensors':[rec]}}}
        got,mag=calibrated_triplet_cg(m)
        self.assertLess(np.max(abs(got-g)),1e-14)
        self.assertAlmostEqual(mag,scale,places=12)

    def test_source_record_malformed_rejected(self):
        src={'symbol':'bad','potential_prefactor':'1',
             'dimensions':[2,2,1,3],'components':['0']*2}
        with self.assertRaises(DictionaryEvidenceError):
            _mat_tensor(src)


@unittest.skipUnless(SRC.exists() and DST.exists(),
                     'Optional full native input files not in local output folder')
class TestT3ParallelANativeEvidenceIntegration(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.md=json.loads(SRC.read_text(encoding='utf-8-sig'))
        cls.rd=json.loads(DST.read_text(encoding='utf-8-sig'))
        cls.result=solve_from_paths(SRC,DST)

    def test_50_component_spaces_and_complex_polynomial_holdouts(self):
        r=self.result
        self.assertEqual(r['status'],'All50RawExtraQuarticInvariantsMatched_ConditionalFrameOnly')
        self.assertEqual(r['additional_quartic_slots'],50)
        self.assertEqual(r['model_count'],16)
        self.assertEqual(r['n_failed_models'],0)
        self.assertLess(r['max_component_fit_residual'],1e-12)
        self.assertLess(r['max_independent_potential_holdout_residual'],1e-12)

    def test_complex_maps_have_expected_nondiagonal_blocks(self):
        for model in ('B_-1','C_-1','E_-2','E_+0'):
            data=self.result['models'][model]
            self.assertEqual(data['rank_matchete'],data['n_extra_quartics'])
            self.assertEqual(data['rank_rgbeta'],data['n_extra_quartics'])
            blocks=[g for g in data['groups'] if len(g['source_indices'])==2]
            self.assertEqual(len(blocks),1)
            A=np.array(blocks[0]['map_rgbeta_from_matchete']['real'])+1j*np.array(
                blocks[0]['map_rgbeta_from_matchete']['imag'])
            self.assertEqual(A.shape,(2,2))
            self.assertGreater(abs(A[0,1]),0.5)
            self.assertGreater(abs(A[1,0]),0.5)

    def test_not_automatically_physically_certified(self):
        r=self.result
        self.assertFalse(r['scope']['physical_scalar_dictionary_certified'])
        self.assertFalse(r['scope']['native_scalar_basis_verified'])
        self.assertFalse(r['scope']['signed_GammaCC_to_Weyl_verified'])
        self.assertFalse(r['scope']['correct_complete_one_loop_ADM_certified'])
        self.assertEqual(r['models']['A_-4']['status'],'NoAdditionalQuartics')

    def test_matchete_cartesian_calibration_not_v2_automatic_identification(self):
        cg,scale=calibrated_triplet_cg(self.md)
        self.assertAlmostEqual(scale,math.sqrt(2/math.sqrt(3)),places=12)
        a=next(x for x in self.md['models']['A_+0']['tensors']
               if x['symbol']=='lambdaT3ExtraQ021010I1')
        b=next(x for x in self.rd['models']['A_+0']['records']
               if x['symbol']=='lambdaT3GeneralQ021010I1')
        raw=_mat_tensor(a)
        naive=np.asarray(b['components'],complex).reshape(b['dimensions'])
        ratio=np.vdot(naive,raw)/np.vdot(naive,naive)
        self.assertGreater(np.linalg.norm(raw-ratio*naive)/np.linalg.norm(raw),.90)
        corrected=_target_tensor(b,cg)
        self.assertLess(np.linalg.norm(raw-scale*corrected),1e-12)

    def test_changed_native_factor_fails_closed(self):
        corrupted=copy.deepcopy(self.rd)
        corrupted['models']['A_+0']['records'][0]['native_coupling_function']='12*Matrix[x][##1] & '
        with self.assertRaises(DictionaryEvidenceError):
            solve_from_uploaded_evidence(self.md,corrupted)

    def test_no_silent_incomplete_quartic_inventory(self):
        corrupted=copy.deepcopy(self.md)
        corrupted['models']['B_-1']['tensors'].pop()
        with self.assertRaises(DictionaryEvidenceError):
            solve_from_uploaded_evidence(corrupted,self.rd)

    def test_manifest_has_actual_evidence_sha256(self):
        self.assertEqual(len(self.result['source_evidence_sha256']['matchete']),64)
        self.assertEqual(len(self.result['source_evidence_sha256']['conditional_rgbeta']),64)

    def test_cli_refuses_overwrite(self):
        with tempfile.TemporaryDirectory() as d:
            outfile=Path(d)/'fit.json'
            command=[sys.executable,'-m','scripts.RunT3ParallelAFromNativeEvidence',
                     '--matchete',str(SRC),'--rgbeta',str(DST),'--output',str(outfile)]
            first=subprocess.run(command,cwd=ROOT,capture_output=True,text=True)
            self.assertEqual(first.returncode,0,first.stderr)
            manifest=json.loads(outfile.read_text())
            self.assertEqual(manifest['additional_quartic_slots'],50)
            again=subprocess.run(command,cwd=ROOT,capture_output=True,text=True)
            self.assertNotEqual(again.returncode,0)


if __name__=='__main__':
    unittest.main()
