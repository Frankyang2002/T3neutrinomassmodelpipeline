"""Batch7: exact source-Hessian algebra and complex scalar basis tests."""
from __future__ import annotations

from dataclasses import dataclass
import json
import tempfile
import sys
from pathlib import Path
import types
import unittest
from unittest.mock import patch

import numpy as np
import sympy as sp

ROOT=Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0,str(ROOT))
from Numerical.lo_certification.benchmark import native_scalar_basis_batch7 as m
from Numerical.lo_certification.benchmark import native_portal_batch6 as b6


@dataclass(frozen=True)
class Leg:
    name: str
    dummy: str|None
    rep: str|None
    conjugated: bool=False


@dataclass(frozen=True)
class Call:
    name: str


def synthetic_seed():
    rows=[]
    for idx, name in sorted(m.EXPECTED_HHSS.items()):
        rows.append({'Index':idx,'HiggsFieldCount':2,
                     'Scalar1FieldCount':2 if idx==2 else 1 if idx in (3,4) else 0,
                     'Scalar2FieldCount':0 if idx==2 else 1 if idx in (3,4) else 2,
                     'TermInputForm':f'EXACT_NATIVE_TERM_{idx}'})
    return {'IntegratedFields':['F'],'ActiveHeavyFields':['S1','S2'],
            'ScalarQuarticTermCount':4,'ScalarQuarticTerms':rows}


def fake_parser():
    H1=Leg('H','i','SU2',True)
    H2=Leg('H','i','SU2',False)
    S1=Leg('NewScalar1',None,None,False)
    B1=Leg('NewScalar1',None,None,True)
    S2a=Leg('NewScalar2','t','triplet',False)
    S2b=Leg('NewScalar2','t','triplet',True)
    legs={2:(H1,H2,S1,B1),
          3:(H1,H2,S1,S2b),
          4:(H1,H2,B1,S2a),
          5:(H1,H2,S2a,S2b)}
    def _find_scalar_legs(term):
        return legs[int(term.rsplit('_',1)[-1])]
    def _find_cg_calls(term):
        idx=int(term.rsplit('_',1)[-1]);return [Call('mix')] if idx in (3,4) else []
    def _cg_value(call,registry,assignment):
        return sp.Integer(1) if assignment['t','triplet']==1 else sp.Integer(0)
    def _find_coupling(term):
        idx=int(term.rsplit('_',1)[-1]);return m.EXPECTED_HHSS[idx],idx==4
    fake=types.SimpleNamespace(sp=sp,load_cg_registry=lambda _: {'mix':'checked'},
        _find_scalar_legs=_find_scalar_legs,_find_cg_calls=_find_cg_calls,
        _find_coupling=_find_coupling,_numeric_prefactor=lambda _:sp.Integer(1),
        _cg_value=_cg_value)
    package=types.ModuleType('RGE.running.intermediate')
    package.ScalarOnlyTensorAdapters=fake
    return package


def cg_matrix():
    z=np.sqrt(.5)
    return np.array([[-1j*z,0,1j*z],[0,1,0],[z,0,z]],complex)


def complex_scalar_hessian(a,b):
    x=np.block([[a.real+b.real,-a.imag-b.imag],
                [a.imag-b.imag,a.real-b.real]])
    ix=np.ravel(np.column_stack((np.arange(4),np.arange(4,8))))
    return x[np.ix_(ix,ix)]


class TestNativeScalarBasisBatch7(unittest.TestCase):
    def test_01_term_inventory_four_native_HHSS(self):
        self.assertEqual(m.EXPECTED_HHSS,{2:'lambdaH1Inv1',3:'lambdaT3',4:'lambdaT3',5:'lambdaH2Inv2'})

    def test_02_source_direct_hessian_separate_from_quartic_factorial(self):
        text=Path(m.__file__).read_text(encoding='utf8')
        self.assertIn('_numeric_prefactor',text)
        self.assertNotIn('build_eft1_quartic_tensor',text)

    def test_03_synthetic_exact_source_quadratic(self):
        fake=fake_parser()
        with patch.dict(sys.modules,{'RGE.running.intermediate':fake}):
            h,parts=m.neutral_vev_hessian(synthetic_seed(),
                {'lambdaH1Inv1':.02,'lambdaH2Inv2':.03,'lambdaT3':.1},246.22)
        v2=246.22**2
        self.assertEqual(set(parts),{'2','3','4','5'})
        self.assertAlmostEqual(h[0,0],-.02*v2/2,places=8)
        self.assertAlmostEqual(h[1,1],-.02*v2/2,places=8)
        for i in range(2,8):
            self.assertAlmostEqual(h[i,i],-.03*v2/2,places=8)
        self.assertAlmostEqual(h[0,2],-.1*v2/2,places=8)
        self.assertAlmostEqual(h[1,3],-.1*v2/2,places=8)
        self.assertLess(np.linalg.norm(h-h.T),1e-10)
        self.assertAlmostEqual(h[0,4],0.)

    def test_04_synthetic_complex_conjugate_mixing(self):
        fake=fake_parser()
        with patch.dict(sys.modules,{'RGE.running.intermediate':fake}):
            h,_=m.neutral_vev_hessian(synthetic_seed(),
                {'lambdaH1Inv1':0.,'lambdaH2Inv2':0.,
                 'lambdaT3':{'real':.12,'imag':.07}},200.)
        a,b=b6.hessian_to_complex(h+np.diag([1e6]*8))
        self.assertLess(np.linalg.norm(b),1e-8)
        self.assertAlmostEqual(abs(a[0,1]),abs(.12+.07j)*200**2/2,places=7)

    def test_05_missing_coupling_fails_closed(self):
        with self.assertRaisesRegex(m.Batch7Error,'incomplete/extra'):
            m.neutral_vev_hessian(synthetic_seed(),{'lambdaT3':.1},246.)

    def test_06_wrong_source_HHSS_inventory_rejected(self):
        seed=synthetic_seed();seed['ScalarQuarticTerms'][0]['Index']=6
        with self.assertRaisesRegex(m.Batch7Error,'inventory changed'):
            m.neutral_vev_hessian(seed,{'lambdaH1Inv1':.02,'lambdaH2Inv2':.03,'lambdaT3':.1},246.)

    def test_07_wrong_fermion_first_sequence_rejected(self):
        seed=synthetic_seed();seed['ActiveHeavyFields']=['S1']
        with self.assertRaisesRegex(m.Batch7Error,'all-fermion-first'):
            m.neutral_vev_hessian(seed,{'lambdaH1Inv1':.02,'lambdaH2Inv2':.03,'lambdaT3':.1},246.)

    def test_08_native_triplet_rotation_is_unitary(self):
        u=m.weight_rotation(cg_matrix())
        self.assertEqual(u.shape,(4,4))
        np.testing.assert_allclose(u.conj().T@u,np.eye(4),atol=1e-15)

    def test_09_nonunitary_rotation_fails_closed(self):
        with self.assertRaisesRegex(m.Batch7Error,'not unitary'):
            m.weight_rotation(np.ones((3,3)))

    def test_10_nondiagonal_native_A_rotates_exactly(self):
        u=m.weight_rotation(cg_matrix())
        a=np.diag([1.4e6,2e6,2e6,2e6]).astype(complex)
        a[0,1]=300+12j;a[1,0]=a[0,1].conjugate()
        aw,bw=m.rotate_native_complex_matrices(a,np.zeros((4,4),complex),u)
        np.testing.assert_allclose(aw,u.conj().T@a@u,rtol=1e-15,atol=1e-10)
        self.assertLess(np.linalg.norm(bw),1e-20)

    def test_11_complex_holomorphic_B_rotates_and_inverts(self):
        u=m.weight_rotation(cg_matrix())
        a=np.diag([1.4e6,2e6,2e6,2e6]).astype(complex)
        b=np.zeros((4,4),complex)
        b[0,2]=b[2,0]=40+120j
        aw,bw=m.rotate_native_complex_matrices(a,b,u)
        np.testing.assert_allclose(u.conj()@bw@u.conj().T,b,atol=1e-10)
        np.testing.assert_allclose(aw,u.conj().T@a@u,atol=1e-10)

    def test_12_hessian_to_weight_conversion_no_conjugation_error(self):
        u=m.weight_rotation(cg_matrix())
        a=np.diag([1.4e6,2e6,2e6,2e6]).astype(complex)
        a[0,3]=100+200j;a[3,0]=a[0,3].conjugate()
        b=np.zeros((4,4),complex);b[1,3]=b[3,1]=50j
        h=complex_scalar_hessian(a,b)
        an,bn=b6.hessian_to_complex(h)
        aw,bw=m.rotate_native_complex_matrices(an,bn,u)
        np.testing.assert_allclose(aw,u.conj().T@a@u,atol=1e-9)
        np.testing.assert_allclose(bw,u.T@b@u,atol=1e-9)

    def test_13_rotation_keeps_scalar_spectrum(self):
        u=m.weight_rotation(cg_matrix())
        a=np.diag([1.2e6,1.4e6,1.5e6,1.6e6]).astype(complex)
        a[0,2]=77+5j;a[2,0]=a[0,2].conjugate()
        aw,bw=m.rotate_native_complex_matrices(a,np.zeros_like(a),u)
        np.testing.assert_allclose(np.linalg.eigvalsh(aw),np.linalg.eigvalsh(a),rtol=1e-12)
        self.assertEqual(np.linalg.norm(bw),0)

    def test_14_legacy_complex_matrix_scale_not_sufficient_for_certification(self):
        self.assertIn('physical_C5_certified',Path(m.__file__).read_text(encoding='utf8'))
        self.assertIn("'gammaCC_to_master_Weyl_signed_factor':None",Path(m.__file__).read_text(encoding='utf8'))

    def test_15_genuine_difference_detected(self):
        a=np.eye(3);b=np.diag([1,2,1]);self.assertGreater(m.matrix_relative(a,b),0.2)

    def test_16_nonfinite_comparison_detected(self):
        with self.assertRaisesRegex(m.Batch7Error,'input mismatch'):
            m.matrix_relative(np.array([np.nan]),np.zeros(1))

    def test_17_bad_weight_shape_rejected(self):
        with self.assertRaisesRegex(m.Batch7Error,'Invalid native triplet'):
            m.weight_rotation(np.eye(4))

    def test_18_unphysical_complex_A_rejected(self):
        with self.assertRaisesRegex(m.Batch7Error,'not Hermitian'):
            m.rotate_native_complex_matrices(np.triu(np.ones((4,4))),np.zeros((4,4)),np.eye(4))

    def test_19_imag_B_symmetric_allowed(self):
        a=np.eye(4,dtype=complex);b=np.diag([1j,2j,0,0])
        aw,bw=m.rotate_native_complex_matrices(a,b,np.eye(4))
        np.testing.assert_allclose(bw,b)

    def test_20_scalar_hessian_only_same_eft1_threshold(self):
        code=Path(m.__file__).read_text(encoding='utf8')
        self.assertIn('integrate out BOTH scalars together at muS',code)
        self.assertNotIn('integrate S1 first',code)

    def test_21_output_owned_and_immutable(self):
        self.assertEqual(m.OWNED,Path('output/general_observable/lo_certification_benchmark'))
        self.assertIn('Refusing to overwrite existing Batch7 output',Path(m.__file__).read_text(encoding='utf8'))

    def test_22_runs_six_combinations(self):
        self.assertEqual(len(m.STAGES)*len(m.BRANCHES),6)

    def _read_batch6_matrix_fixture(self, payload):
        with tempfile.TemporaryDirectory() as folder:
            path=Path(folder)/'Hessian_real_interleaved_GeV2.json'
            path.write_text(json.dumps(payload,allow_nan=True),encoding='utf8')
            return m.read_real_hessian(path,'Batch6 test Hessian')

    def test_24_batch6_hessian_written_as_bare_json_array(self):
        matrix=np.diag([1.3e6,1.3e6]+[1e6]*6)
        result=self._read_batch6_matrix_fixture(matrix.tolist())
        np.testing.assert_array_equal(result,matrix)
        self.assertEqual(result.shape,(8,8))

    def test_25_generic_object_reader_still_rejects_matrix(self):
        with tempfile.TemporaryDirectory() as folder:
            path=Path(folder)/'bare_matrix.json'
            path.write_text(json.dumps(np.eye(8).tolist()),encoding='utf8')
            with self.assertRaisesRegex(m.Batch7Error,'Expected JSON object'):
                m.read(path)
            np.testing.assert_array_equal(m.read_real_hessian(path,'Hessian'),np.eye(8))

    def test_26_batch6_matrix_object_rejected(self):
        with self.assertRaisesRegex(m.Batch7Error,'expected bare 8x8'):
            self._read_batch6_matrix_fixture({'real':np.eye(8).tolist()})

    def test_27_wrong_hessian_shape_rejected(self):
        with self.assertRaisesRegex(m.Batch7Error,'expected bare 8x8'):
            self._read_batch6_matrix_fixture(np.eye(7).tolist())
        with self.assertRaisesRegex(m.Batch7Error,'expected bare 8x8'):
            self._read_batch6_matrix_fixture([[1]*8 for _ in range(7)]+[[1]*7])

    def test_28_nonfinite_or_nonreal_hessian_rejected(self):
        bad=np.eye(8).tolist();bad[0][0]=float('nan')
        with self.assertRaisesRegex(m.Batch7Error,'nonfinite or non-real'):
            self._read_batch6_matrix_fixture(bad)
        bad=np.eye(8).tolist();bad[0][0]='2+3j'
        with self.assertRaisesRegex(m.Batch7Error,'nonfinite or non-real'):
            self._read_batch6_matrix_fixture(bad)

    def test_29_nonsymmetric_batch6_hessian_rejected(self):
        bad=np.eye(8).tolist();bad[0][1]=1.0
        with self.assertRaisesRegex(m.Batch7Error,'nonsymmetric Hessian'):
            self._read_batch6_matrix_fixture(bad)

    def test_30_missing_batch6_hessian_rejected(self):
        with tempfile.TemporaryDirectory() as folder:
            with self.assertRaisesRegex(m.Batch7Error,'Missing Batch6'):
                m.read_real_hessian(Path(folder)/'missing.json','Batch6 Hessian')

    def test_23_optional_local_native_preflight(self):
        if not (ROOT/m.b4.BENCH).is_file():
            self.skipTest('Full local repository is not mounted in test environment')
        if not (ROOT/m.OWNED/'batch6_portal_source01'/'summary.json').is_file():
            self.skipTest('Completed Batch6 release is not mounted in test environment')
        status=m.preflight(ROOT,'source01','source01','source01')
        self.assertEqual(status['status'],'PREFLIGHT_OK_NO_OUTPUT_WRITTEN')


if __name__ == '__main__':
    unittest.main()
