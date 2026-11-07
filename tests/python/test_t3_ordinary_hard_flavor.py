"""Offline tests for the source-locked nondegenerate hard scalar flavor lift."""
from __future__ import annotations

import json
import math
from pathlib import Path
import shutil
import tempfile
import unittest

import numpy as np

from Numerical.running.T3OrdinaryHardFlavor import (
    HardFlavorError, compare_archived_v6, decode_complex, evaluate,
    mass_kernel, running_parameters, source_pole_matrix, verify_hard_stage,
)

ROOT = Path(__file__).resolve().parents[2]
FIXTURE = ROOT/'tests/fixtures/t3_ordinary_hard_c5_candidate.txt'


def cmatrix(a):
    return [[{'real':float(complex(z).real),'imag':float(complex(z).imag)}
             for z in row] for row in np.asarray(a)]


def example_running(model='D_+2', *, different=False):
    yy = np.diag([.00042,.00043,.00044]).astype(complex)
    y2 = np.diag([.00047,.00046,.00045]).astype(complex)
    if different:
        yy[0,1] = 1.5e-4 + .3e-4j
        y2[2,1] = 2e-4 - .2e-4j
    return {'model':model,'status':'NumericalRenormalisableRunningOnly',
            'scales_gev':{'F':1e5,'S':1e3},
            'uv_at_fermion_threshold':{'MF':cmatrix(np.diag([110000,150000,190000])),
                                       'y1':cmatrix(yy),'y2':cmatrix(y2)},
            'eft1_at_fermion_threshold':{'mS1Sq':1.8e6,'mS2Sq':2.5e6,
                                        'lambdaT3':{'real':.009,'imag':0}},
            'eft1_at_scalar_threshold':{'mS1Sq':2.1e6,'mS2Sq':2.8e6,
                                       'lambdaT3':{'real':.012,'imag':0}}}


def stage(root, model):
    folder=root/model.replace('+','p').replace('-','m')
    folder.mkdir(parents=True)
    shutil.copyfile(FIXTURE, folder/'hard_c5_candidate.txt')
    (folder/'hard_scalar_report.json').write_text(json.dumps({
        'Status':'HardScalarCandidateExported','Model':model,
        'HardOneLoopMatched':True,'HardHbarSectorExtracted':True,
        'PerHeavyMassEigenvalueTemplate':True,'FermionMatchingRecomputed':False,
        'CompleteC5Assembled':False,'SU2NormalizationIndependentlyCertified':False,
        'ThreeHeavyMassNumericalSumImplemented':False}),encoding='utf-8')
    (folder/'hard_scalar_input_audit.json').write_text(json.dumps({
        'Status':'Passed','Model':model,'CheckLagrangianPassed':True,
        'DiagonalNFlavorMassIndexEliminatedForTemplate':True,
        'PerEigenvalueHeavyMassSumRequired':True,'TreeWilsonTermCount':2,
        'OddHeavyScalarDegreeAbsent':True}),encoding='utf-8')
    (folder/'provenance.json').write_text(json.dumps({
        'model':model,'status':'Ready','three_heavy_mass_numerical_sum_required':True,
        'indexed_fermion_mass_factors':2}),encoding='utf-8')
    return folder


class T3HardFlavorTests(unittest.TestCase):
    def setUp(self):
        self.tmp=tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.root=Path(self.tmp.name)

    def test_stable_scalar_kernel_for_equal_mass(self):
        a=2.2e6
        self.assertAlmostEqual(mass_kernel(a,a,1e3),math.log(1e6/a),places=12)
        self.assertAlmostEqual(mass_kernel(a,a*(1+1e-12),1e3),
                               mass_kernel(a,a,1e3),places=10)

    def test_scalar_kernel_vs_direct_expression(self):
        a,b,mu=2e6,2.5e6,1e3
        expected=1+(a*math.log(mu*mu/a)-b*math.log(mu*mu/b))/(a-b)
        self.assertAlmostEqual(mass_kernel(a,b,mu),expected,places=12)

    def test_exact_source_locked(self):
        stage_dir=stage(self.root,'D_+2')
        obj=verify_hard_stage(stage_dir,'D_+2')
        self.assertEqual(obj['status'],'VerifiedOneEigenvalueTemplate')

    def test_source_mutation_is_rejected(self):
        stage_dir=stage(self.root,'D_+2')
        with (stage_dir/'hard_c5_candidate.txt').open('a',encoding='utf-8') as file:
            file.write(' + 123')
        with self.assertRaisesRegex(HardFlavorError,'unexpected Mathematica'):
            verify_hard_stage(stage_dir,'D_+2')

    def test_wrong_model_fails(self):
        stage_dir=stage(self.root,'D_+2')
        with self.assertRaises(HardFlavorError):
            verify_hard_stage(stage_dir,'B_-1')

    def test_nondegenerate_three_flavour_benchmark(self):
        stage_dir=stage(self.root,'D_+2')
        doc=example_running()
        result=evaluate(model='D_+2',stage_folder=stage_dir,running_doc=doc)
        self.assertEqual(result['status'],'HardThreeFlavorAndDirectLogComputed_NotFullC5')
        self.assertEqual(len(result['per_heavy_eigenstate']),3)
        self.assertGreater(result['checks']['pole_residue_norm_gev_inverse'],0)
        self.assertLess(result['checks']['hard_plus_direct_log_matching_scale_variation_residual'],1e-12)
        self.assertFalse(result['checks']['complete_sequential_C5_certified'])
        self.assertFalse(result['checks']['independent_MSbar_counterterm_proved'])
        self.assertEqual(result['heavy_fermion_mass_eigenvalues_gev'],[110000.,150000.,190000.])
        self.assertEqual(len(result['hard_plus_running_subset_metrics']['conditional_takagi_masses_ev_ascending']),3)

    def test_flavour_components_add_to_pole(self):
        stage_dir=stage(self.root,'A_-4')
        r=evaluate(model='A_-4',stage_folder=stage_dir,running_doc=example_running('A_-4',different=True))
        summed=sum((decode_complex(k['hard_pole_residue_gev_inverse'],'x')
                     for k in r['per_heavy_eigenstate']),np.zeros((3,3),complex))
        pole=decode_complex(r['raw_unsubtracted_hard_1_over_epsilon_residue_gev_inverse'],'pole')
        self.assertTrue(np.allclose(summed,pole,rtol=1e-12,atol=0))
        self.assertTrue(np.allclose(pole,pole.T,rtol=1e-12,atol=0))

    def test_complex_bars_and_lambda(self):
        y1=np.array([[1+2j,0,0],[0,1,0],[0,0,1]],complex)
        y2=np.eye(3,dtype=complex)
        pole,rows=source_pole_matrix(y1,y2,np.array([1e5,2e5,3e5]),.012+.03j)
        expected=(-2/math.sqrt(3)/(16*math.pi**2))*(.012+.03j)*(1-2j)/1e5
        self.assertAlmostEqual(pole[0,0].real,expected.real,places=20)
        self.assertAlmostEqual(pole[0,0].imag,expected.imag,places=20)

    def test_non_diagonal_heavy_mass_rejected(self):
        d=example_running()
        d['uv_at_fermion_threshold']['MF'][0][1]['real']=99.
        with self.assertRaisesRegex(HardFlavorError,'not diagonal'):
            running_parameters(d,'D_+2')

    def test_negative_heavy_mass_rejected(self):
        d=example_running()
        d['uv_at_fermion_threshold']['MF'][0][0]['real']=-1
        with self.assertRaisesRegex(HardFlavorError,'positive real'):
            running_parameters(d,'D_+2')

    def test_matching_scale_identity_non_diagonal(self):
        d=example_running(different=True)
        r=evaluate(model='D_+2',stage_folder=stage(self.root,'D_+2'),running_doc=d)
        hard=decode_complex(r['hard_finite_scalar_gev_inverse'],'h')
        running=decode_complex(r['direct_eft1_llss_log_gev_inverse'],'r')
        combined=decode_complex(r['c5_hard_finite_plus_direct_log_subset_gev_inverse'],'c')
        self.assertTrue(np.allclose(combined,hard+running,rtol=1e-12,atol=0))
        self.assertLess(r['checks']['hard_plus_direct_log_matching_scale_variation_residual'],1e-12)

    def test_independent_archived_llss_parity(self):
        d=example_running()
        r=evaluate(model='D_+2',stage_folder=stage(self.root,'D_+2'),running_doc=d)
        responses={'model':'D_+2',
                   'C5_joint_scalar_finite_bubble_per_N_gev_inverse':r['hard_finite_scalar_gev_inverse'],
                   'C5_EFT1_induced_leading_log_per_N_gev_inverse':r['direct_eft1_llss_log_gev_inverse'],
                   'C5_LLSS_total_one_loop_per_N_gev_inverse':r['c5_hard_finite_plus_direct_log_subset_gev_inverse']}
        path=self.root/'reference.json'
        path.write_text(json.dumps({'model':'D_+2','complete_one_loop_C5_certified':False,
                                    'response':responses}),encoding='utf-8')
        parity=compare_archived_v6(r,path)
        self.assertEqual(parity['status'],'SourceBackedLLSSKernelAgreement')

    def test_reference_mismatch_not_silenced(self):
        d=example_running()
        r=evaluate(model='D_+2',stage_folder=stage(self.root,'D_+2'),running_doc=d)
        zero={'real':np.zeros((3,3)).tolist(),'imag':np.zeros((3,3)).tolist()}
        responses={'model':'D_+2',
                   'C5_joint_scalar_finite_bubble_per_N_gev_inverse':zero,
                   'C5_EFT1_induced_leading_log_per_N_gev_inverse':zero,
                   'C5_LLSS_total_one_loop_per_N_gev_inverse':zero}
        path=self.root/'wrong.json'
        path.write_text(json.dumps({'model':'D_+2','complete_one_loop_C5_certified':False,
                                    'response':responses}),encoding='utf-8')
        self.assertEqual(compare_archived_v6(r,path)['status'],
                         'SourceBackedLLSSKernelMismatch')

    def test_cli_end_to_end_and_resume(self):
        from scripts.RunT3OrdinaryHardFlavorBatch import main
        hard=self.root/'hard'
        running=self.root/'running'/'D_p2'
        running.mkdir(parents=True)
        stage(hard,'D_+2')
        (running/'running.json').write_text(json.dumps(example_running()),encoding='utf-8')
        argv=['--models','D_+2','--run-id','offline_check_01',
              '--hard-root',str(hard),'--running-root',str(running.parent),
              '--output-root',str(self.root/'output'),
              '--reference-root',str(self.root/'missing_refs')]
        self.assertEqual(main(argv+['--check-only']),0)
        self.assertEqual(main(argv),0)
        self.assertEqual(main(argv+['--resume']),0)
        result=self.root/'output'/'offline_check_01'/'D_p2'/'result.json'
        self.assertTrue(result.is_file())
        self.assertEqual(json.loads(result.read_text())['status'],
                         'HardThreeFlavorAndDirectLogComputed_NotFullC5')


if __name__=='__main__':
    unittest.main()
