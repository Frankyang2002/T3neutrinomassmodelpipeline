"""Offline signed p=0 A_-4 UV Dirac-F exchange and input fail-closed tests.
Project path: output/general_observable/lo_fermion_boundary/test_native_zero_momentum_a_m4.py
"""
from __future__ import annotations
import hashlib
import unittest
from copy import deepcopy
import native_zero_momentum_a_m4 as check


def block(txt):
    return {
        'head': 'Association', 'inputform': txt, 'inputform_complete': True,
        'inputform_sha256': hashlib.sha256(txt.encode('utf8')).hexdigest(),
        'inputform_characters': len(txt),
    }


def fixture():
    source = 'NCM[Transp[Field[l, Fermion, {}]], DiracProduct[GammaCC, Proj[-1]], Field[l, Fermion, {}]])/Coupling[MF, {r}, 0]'
    return {
        'format': 't3-native-matchete-fr-evidence-v1',
        'model': 'A_-4',
        'status': 'NativeUV3pointAndEFT1ContactVerticesExtracted_AmplitudeOpen',
        'native_source_index': 2,
        'native_source_inputform': source,
        'uv_three_point_rules': block('<|{Conj[NewScalar1], NewFermion, l} -> I*Bar[Coupling[y1,{r}]], {NewScalar2, l, Conj[NewFermion]} -> (-I)*Bar[Coupling[y2,{r}]]|>'),
        'uv_yukawa_lagrangian': block('Bar[Coupling[y1,{r}]]*NCM[Transp[Field[l,{}]],DiracProduct[GammaCC, Proj[-1]], Field[NewFermion,{}]] + Bar[Coupling[y2,{r}]]*NCM[Bar[Field[NewFermion,{}]], DiracProduct[Proj[-1]], Field[l,{}]]'),
        'uv_mass_lagrangian_source': block('-(Coupling[MF,{r}]*NCM[Bar[Field[NewFermion,{}]],Field[NewFermion,{}]])'),
        'original_holomorphic_contact': block(source),
        'eft1_four_point_rules': block('{Conj[NewScalar1], NewScalar2, l, l} -> (I*Bar[Coupling[y1,{i}]]*Bar[Coupling[y2,{j}]] + I*Bar[Coupling[y1,{j}]]*Bar[Coupling[y2,{i}]])*NCM[Transp[Field[Ext[l, 1],{}]],DiracProduct[GammaCC, Proj[-1]], Field[Ext[l,2],{}]]|>'),
        'physical_master_wilson_normalization_certified': False,
        'independent_connected_uv_four_point_amplitude_computed': False,
        'matchete_gamma_cc_to_weyl_phase_certified': False,
    }


class TestSignedA(unittest.TestCase):
    def test_native_order_dirac_wick_exact(self):
        self.assertTrue(check.exact_check()['uv_ij_equals_native_contact_ij_exact'])

    def test_external_lepton_swap(self):
        a,b=check.connected_exchange()
        self.assertTrue((a+b).equal(
            check.spin_bilinear(check.LI,check.LJ).scale((check.Y1I*check.Y2J+check.Y1J*check.Y2I)/check.M)))

    def test_zero_momentum_minkowski_propagator_sign(self):
        r=check.exact_check()
        self.assertEqual(r['heavy_dirac_propagator_at_p_zero'],'-i/M times identity_4')

    def test_two_time_orderings_agree(self):
        v=check.yukawa_one(check.LI,check.Y1I)
        w=check.yukawa_two(check.LJ,check.Y2J)
        self.assertTrue(check.contract_free_dirac_at_zero(v*w).equal(check.contract_free_dirac_at_zero(w*v)))

    def test_charge_conjugation_phase_cannot_be_determined(self):
        plus,_=check.connected_exchange(check.C)
        minus,_=check.connected_exchange(-check.C)
        self.assertTrue(minus.equal(-plus))

    def test_free_berezin_propagator_sign(self):
        f=check.Exterior.generator(0)
        fb=check.Exterior.generator(4)
        weight=check.Exterior({():1})-(fb*f).scale(check.M)
        z=weight.terms[(0,4)]
        normalized=(f*fb*weight).terms[(0,4)]/z
        self.assertEqual(normalized,1/check.M)

    def test_input_native_sign_and_both_orderings(self):
        native=check.required_report(fixture())
        self.assertEqual(native['uv_vertex_signs']['y1_bar'],'+i')
        self.assertEqual(native['holomorphic_native_eft1_vertex_orderings'],2)

    def test_fails_wrong_uv_vertex_sign(self):
        f=fixture()
        txt=f['uv_three_point_rules']['inputform'].replace('-> (-I)*Bar[Coupling[y2','-> I*Bar[Coupling[y2')
        f['uv_three_point_rules']=block(txt)
        with self.assertRaises(ValueError):check.required_report(f)

    def test_fails_missing_fermion_ordering(self):
        f=fixture()
        txt=f['eft1_four_point_rules']['inputform'].replace(' + I*Bar[Coupling[y1,{j}]]*Bar[Coupling[y2,{i}]]','')
        f['eft1_four_point_rules']=block(txt)
        with self.assertRaises(ValueError):check.required_report(f)

    def test_fails_corrupt_digest(self):
        f=fixture();f['uv_mass_lagrangian_source']['inputform_sha256']='0'*64
        with self.assertRaises(ValueError):check.required_report(f)

    def test_fails_physical_claim(self):
        f=fixture();f['physical_master_wilson_normalization_certified']=True
        with self.assertRaises(ValueError):check.required_report(f)


if __name__=='__main__':
    unittest.main(verbosity=2)
