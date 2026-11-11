"""Offline fail-closed regression tests; do not require local Matchete/Wolfram."""
from __future__ import annotations
import json
from pathlib import Path
from tempfile import TemporaryDirectory
import unittest
from unittest.mock import patch
import numpy as np
from scripts.RunT3WorkerCNumericOverlap import (
    AUDIT_INPUTS, COMPONENTS, preflight, run, load_complex_matrix,
    isolated_numeric_evaluator, DIAGNOSTIC_MASS_LINE, LEGACY_MASS_LINE,
)


class TestWorkerCNumericOverlap(unittest.TestCase):
    def setUp(self):
        self.temp=TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.base=Path(self.temp.name)
        self.restricted=self.base/'restricted';self.restricted.mkdir()
        self.profile=self.base/'profile';self.profile.mkdir()
        self.config=self.base/'config.json'
        self.evaluator=self.base/'evaluator.wl'
        self.evaluator.write_text('(* untouched evaluator source *)\n'
            'mnuEV = N[-(vev^2/2) matrixC5 10^9, 15];\n'
            + LEGACY_MASS_LINE, encoding='utf8')
        for n in AUDIT_INPUTS:
            (self.profile/n).write_text('{}',encoding='utf8')
        for filename in COMPONENTS.values():
            (self.restricted/filename).write_bytes(b'WXF fixture')
        (self.restricted/'result.json').write_text(json.dumps({
            'status':'RestrictedNF1PoleFreeCandidateExported',
            'restricted_total_pole_zero_exact_NF1':True,
            'safe_to_add_to_worker_c_joint':False}),encoding='utf8')
        self.params={
            'MF_GeV':4000.,'MS1_GeV':3200.,'MS2_GeV':2400.,
            'AbsKappaS_GeV2':90000.,'MuScalar_GeV':3000.,'NF':1,
            'Yukawas':{'y1':[[1.],[2.],[3.]],'y2':[[2.],[3.],[4.]],
                       'y2Conj':[[1.],[2.],[3.]]}}
        self.config.write_text(json.dumps(self.params),encoding='utf8')
        self.out=self.base/'fresh'

    def pf(self):
        return preflight(self.restricted,self.profile,self.config,self.evaluator,self.out)

    def test_preflight_and_explicit_conjugate_zero(self):
        p=self.pf()
        self.assertTrue(p['independent_y2Conj_zeroed_in_isolated_benchmark'])
        self.assertEqual(p['restricted_benchmark']['Yukawas']['y2Conj'],[[0.],[0.],[0.]])
        self.assertEqual(p['restricted_benchmark']['PoleHandling'],'FinitePartDiagnostic')
        self.assertEqual(json.loads(self.config.read_text())['Yukawas']['y2Conj'][0][0],1.)

    def test_takagi_zero_mode_evaluator_isolation(self):
        path=self.base/'diagnostic_copy.wl'
        original=self.evaluator.read_bytes()
        isolated_numeric_evaluator(self.evaluator,path)
        self.assertIn(DIAGNOSTIC_MASS_LINE,path.read_text())
        self.assertNotIn(LEGACY_MASS_LINE,path.read_text())
        self.assertEqual(self.evaluator.read_bytes(),original)

    def test_changed_evaluator_is_rejected(self):
        self.evaluator.write_text('unrecognized massesEV',encoding='utf8')
        with self.assertRaisesRegex(ValueError,'source changed'):
            isolated_numeric_evaluator(self.evaluator,self.base/'copy.wl')

    def test_rank_two_weinberg_matrix_zero_mass(self):
        y1=np.array([1.,2.,3.]); y2=np.array([2.,3.,4.])
        c5=np.outer(y1,y2)+np.outer(y2,y1)
        masses=np.linalg.svd(c5,compute_uv=False)
        self.assertEqual(len(masses),3)
        self.assertLess(masses[-1],1e-12)

    def test_source_is_fail_closed(self):
        (self.restricted/'result.json').write_text(json.dumps({
            'status':'RestrictedNF1PoleFreeCandidateExported',
            'restricted_total_pole_zero_exact_NF1':False,
            'safe_to_add_to_worker_c_joint':False}))
        with self.assertRaisesRegex(ValueError,'exact NF=1 pole'):
            self.pf()

    def test_refuses_existing_output(self):
        self.out.mkdir()
        with self.assertRaisesRegex(ValueError,'overwrite'):
            self.pf()

    def test_nf_not_one_rejected(self):
        self.params['NF']=2
        self.config.write_text(json.dumps(self.params))
        with self.assertRaisesRegex(ValueError,'NF=1'):
            self.pf()

    def test_isolated_wolfram_calls_and_output_validation(self):
        original_eval=self.evaluator.read_bytes()
        original_config=self.config.read_bytes()
        calls=[]
        def fake_run(command,**kwargs):
            calls.append(command)
            self.assertNotEqual(Path(command[2]), self.evaluator)
            self.assertIn(DIAGNOSTIC_MASS_LINE, Path(command[2]).read_text())
            root=Path(command[-2]); stage=root/'data/general_d_m2_scalar_threshold/profile_tree'
            self.assertTrue((stage/'profile_c5_flavour_candidate_total_provisional.wxf').exists())
            self.assertTrue(all((stage/n).exists() for n in AUDIT_INPUTS))
            self.assertFalse(str(stage).startswith(str(self.profile)))
            entry={'Status':'ProvisionalFinitePartBenchmarkEvaluated',
                   'PoleHandling':'FinitePartDiagnostic','FinitePartOnly':True,
                   'NoNewMatchPerformed':True,'NumericalBenchmarkCertifiedPhysical':False,
                   'RenormalizationCertified':False,
                   'DimRegPolePresent':root.name=='hard',
                   'C5PoleResidueMaxAbs_GeVinv':1e-15 if root.name=='hard' else 1e-30}
            (stage/'profile_c5_rough_benchmark_result.json').write_text(json.dumps(entry))
            m=np.eye(3)* (2.0 if root.name=='hard' else 5.0)
            np.savetxt(stage/'profile_c5_rough_benchmark_C5_real_GeVinv.csv',m,delimiter=',')
            np.savetxt(stage/'profile_c5_rough_benchmark_C5_imag_GeVinv.csv',np.zeros((3,3)),delimiter=',')
            class Proc:returncode=0
            return Proc()
        with patch('scripts.RunT3WorkerCNumericOverlap.subprocess.run',side_effect=fake_run):
            report_path=run(restricted=self.restricted,profile=self.profile,config=self.config,
                            evaluator=self.evaluator,output=self.out)
        report=json.loads(report_path.read_text())
        self.assertEqual(report['status'],'RestrictedHardAndTotalFinitePartsEvaluatedDiagnosticOnly')
        self.assertEqual(len(calls),2)
        self.assertEqual(report['overlap']['status'],'NotAttemptedNoCanonicalTensorInput')
        self.assertEqual(report['components']['hard']['finite_c5_gev_inv']['real'][0][0],2.)
        self.assertEqual(report['total_minus_hard_not_nonbubble_correction']['real'][0][0],3.)
        self.assertFalse(report['safe_to_add_to_worker_c_joint'])
        self.assertTrue(report['components']['hard']['pole_reported'])
        self.assertTrue(report['total_numeric_pole_consistent_with_exact_symbolic_zero'])
        self.assertTrue(all(c['finite_part_only'] for c in report['components'].values()))
        self.assertEqual(self.evaluator.read_bytes(),original_eval)
        self.assertTrue(report['isolated_takagi_zero_mode_workaround'])
        self.assertTrue((self.out/'worker_c_benchmark_evaluator.wl').is_file())
        self.assertEqual(self.config.read_bytes(),original_config)

    def test_rejects_large_total_pole(self):
        def fake_run(command, **kwargs):
            root=Path(command[-2]); stage=root/'data/general_d_m2_scalar_threshold/profile_tree'
            entry={'Status':'ProvisionalFinitePartBenchmarkEvaluated',
                   'PoleHandling':'FinitePartDiagnostic','FinitePartOnly':True,
                   'NoNewMatchPerformed':True,'NumericalBenchmarkCertifiedPhysical':False,
                   'RenormalizationCertified':False,'DimRegPolePresent':True,
                   'C5PoleResidueMaxAbs_GeVinv':1e-15}
            (stage/'profile_c5_rough_benchmark_result.json').write_text(json.dumps(entry))
            np.savetxt(stage/'profile_c5_rough_benchmark_C5_real_GeVinv.csv',np.eye(3),delimiter=',')
            np.savetxt(stage/'profile_c5_rough_benchmark_C5_imag_GeVinv.csv',np.zeros((3,3)),delimiter=',')
            class Proc:returncode=0
            return Proc()
        with patch('scripts.RunT3WorkerCNumericOverlap.subprocess.run',side_effect=fake_run):
            with self.assertRaisesRegex(RuntimeError,'exceeds diagnostic noise bound'):
                run(restricted=self.restricted,profile=self.profile,config=self.config,
                    evaluator=self.evaluator,output=self.out)
        self.assertEqual(json.loads((self.out/'result.json').read_text())['status'],'Blocked')

    def test_exit_zero_without_json_is_error(self):
        class Proc:returncode=0
        with patch('scripts.RunT3WorkerCNumericOverlap.subprocess.run',return_value=Proc()):
            with self.assertRaisesRegex(RuntimeError,'result is missing'):
                run(restricted=self.restricted,profile=self.profile,config=self.config,
                    evaluator=self.evaluator,output=self.out)
        report=json.loads((self.out/'result.json').read_text())
        self.assertEqual(report['status'],'Blocked')

    def test_invalid_matrix_rejected(self):
        real=self.base/'real.csv'; imag=self.base/'imag.csv'
        np.savetxt(real,np.zeros((2,2)),delimiter=',')
        np.savetxt(imag,np.zeros((2,2)),delimiter=',')
        with self.assertRaisesRegex(ValueError,'Invalid 3x3'):
            load_complex_matrix(real,imag)


if __name__=='__main__':unittest.main()
