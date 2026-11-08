"""Only independent numerical and fail-closed native contract tests; no Matchete mocked as certified."""
from __future__ import annotations

import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

from Numerical.lo_certification.spinor import a_minus4_native_phase_audit as probe


def valid_native():
    return {
        'format':probe.NATIVE_FORMAT, 'status':probe.NATIVE_STATUS,
        'model':'A_-4',
        'checks_action':{name:True for name in probe.ACTIONS},
        'checks_abstract_spinor':{name:True for name in probe.SPINOR},
        'native_heavy_fields':['F'], 'native_light_fields':['S1','S2'],
        'native_yukawa_1_terms_inputform':['Coupling[y1, {}] * DiracProduct[GammaCC]'],
        'native_yukawa_2_terms_inputform':['Coupling[y2, {}] * Field[NewScalar2]'],
        'native_heavy_free_inputform':'NewFermion*Coupling[MF] * NewFermion',
        'native_spinor_bilinear_explicit_inputform':'DiracProduct[GammaCC]',
        'native_concrete_4x4_cc_matrix':None,
        'native_4x4_CC_to_Weyl_epsilon_ratio':None,
        'native_connected_UV_amplitude_calculated':False,
        'native_source_matched_contact_amplitude_calculated':False,
        'full_F_matching_performed_by_this_probe':False,
        'physical_certification':False,
    }


class TestNativeAminus4Phase(unittest.TestCase):
    def test_exact_complex_spinor_signs_and_phase(self):
        w=probe.numeric_spinor_convention()
        self.assertEqual(w['C_PL_to_Weyl_ratio'],{'real':1.,'imag':0.})
        self.assertEqual(w['C_inverse_PL_to_Weyl_ratio'],{'real':-1.,'imag':0.})
        self.assertFalse(w['native_Matchete_sign_measured'])
        self.assertNotEqual(w['sample_C_amplitude'],w['sample_C_inverse_amplitude'])

    def test_exact_grassmann_Dirac_exchange(self):
        w=probe.exact_dirac_saddle()
        self.assertTrue(w['eom_equality'])
        self.assertFalse(w['native_Aminus4_UV_vertex_identification_established'])

    def test_abstract_native_report_keeps_physical_flag_false(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/'native.json'
            p.write_text(json.dumps(valid_native()),encoding='utf-8')
            data=probe.load_native(p)
            self.assertIsNone(data['native_concrete_4x4_cc_matrix'])

    def test_missing_GammaCC_native_action_rejected(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/'native.json'
            v=valid_native()
            v['native_yukawa_1_terms_inputform']=['Coupling[y1, {}]']
            p.write_text(json.dumps(v),encoding='utf-8')
            with self.assertRaises(probe.NativeBlocked):probe.load_native(p)

    def test_forged_native_numeric_sign_rejected(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/'native.json'
            v=valid_native()
            v['native_4x4_CC_to_Weyl_epsilon_ratio']=1
            p.write_text(json.dumps(v),encoding='utf-8')
            with self.assertRaises(probe.NativeBlocked):probe.load_native(p)

    def test_wrong_heavy_field_inventory_rejected(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/'native.json'
            v=valid_native()
            v['native_heavy_fields']=['F','S1','S2']
            p.write_text(json.dumps(v),encoding='utf-8')
            with self.assertRaises(probe.NativeBlocked):probe.load_native(p)

    def test_conditional_boundary_is_rebuilt_not_trusted_from_saved_status(self):
        with tempfile.TemporaryDirectory() as d:
            root=Path(d)
            path=root/'native.json'
            path.write_text(json.dumps(valid_native()),encoding='utf-8')
            pseudo_source={
                'physical_certification':False,
                'physical_C_ijab_at_MF':None,
                'sources':{'native_full_tree':{'sha256':'a'*64}},
                'numeric_flavour_checks':{
                    'mu_F_gev':100000.0,
                    'UV_evolved_3x3_flavour_recomputed':True},
                'independent_conditional_tensor_entries':[{'i':1}],
                'basis_contract':{'full_tensor_shape':[45,45,12,12]}}
            with patch.object(probe.boundary,'verify_from_sources',return_value=pseudo_source) as rebuilt:
                result=probe.evaluate(root,path)
            rebuilt.assert_called_once_with(root.resolve(),check_running=True)
            self.assertFalse(result['certified'])
            self.assertIsNone(result['physical_C_ijab_at_MF'])
            self.assertIsNone(result['native_GammaCC_to_Weyl_factor'])
            self.assertEqual(result['conditional_tensor_full_shape'],[45,45,12,12])

    def test_failed_native_cc_antisymmetry_rejected(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/'native.json'
            v=valid_native()
            v['checks_abstract_spinor']['native_CC_transpose_antisymmetric']=False
            p.write_text(json.dumps(v),encoding='utf-8')
            with self.assertRaises(probe.NativeBlocked):probe.load_native(p)


if __name__=='__main__':unittest.main()
