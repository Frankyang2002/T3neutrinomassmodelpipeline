"""B1 v2 action tests: synthetic strict negatives and native-shaped positives.

These are NOT a substitute for running the installed Matchete native probe.
The positive witness shapes mirror three source reports read from GitHub main.
"""
from __future__ import annotations

from copy import deepcopy
import json
import hashlib
import tempfile
from unittest.mock import patch
from pathlib import Path
import sys
import unittest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from Numerical.lo_certification.physical_normalisation.b1_native_action_contract import (
    ActionContractError, parse, multiplicative, verify_record, verify_term,
)
from scripts.RunPhysicalB1NativeCoefficientAudit import check_all16, check_cc, build
from Numerical.lo_certification.physical_normalisation.b1_matchete_weyl import model_specs


F = 'Field[NewFermion,Fermion,{Index[d$$1,NFlavor]},{}]'
L = 'Field[l,Fermion,{Index[d$$1,Flavor]},{}]'
S1 = 'Field[NewScalar1,Scalar,{},{}]'
S2 = 'Field[NewScalar2,Scalar,{},{}]'
M = 'Coupling[MF,{},0]'

def p(pos): return f'DiracProduct[Proj[{pos}]]'
def cp(pos): return f'DiracProduct[GammaCC,Proj[{pos}]]'
def cm(): return 'DiracProduct[GammaCC]'
def n(*args): return 'NCM['+','.join(args)+']'
def cc(y): return f'Coupling[{y},{{Index[d$$1,Flavor],Index[d$$1,NFlavor]}},0]'
def cg(y, hc=False):
    name = 'T3'+y.upper()+'CG'
    return f'CG[{f"Bar[{name}]" if hc else name},{{}}]'


def example(model):
    majorana = model != 'A_p0'
    mass = ('-1/2*('+M+'*'+n(f'Transp[{F}]',cm(),F)+')' if majorana else
            '-('+M+'*'+n(f'Bar[{F}]',F)+')')
    yuks = {}
    for v in ('y1','y2'):
        scalar = S1 if v == 'y1' else S2
        a_scalar = scalar if v == 'y1' else f'Bar[{scalar}]'
        h_scalar = f'Bar[{scalar}]' if v == 'y1' else scalar
        if not majorana and v == 'y1':
            op = n(f'Bar[{L}]',cp(1),f'Transp[Bar[{F}]]')
            hc = n(f'Transp[{L}]',cp(-1),F)
        else:
            op = n(f'Bar[{L}]',p(1),F)
            if majorana:
                hc = n(f'Transp[{F}]',cp(-1),L)
            else:
                hc = n(f'Bar[{F}]',p(-1),L)
        orig = '*'.join((cg(v),cc(v),a_scalar,op))
        conjugate = '*'.join((f'Bar[{cc(v)}]',h_scalar,cg(v,True),hc))
        yuks[v] = {'original_terms_inputform':[orig], 'hc_terms_inputform':[conjugate],
                   'original_term_count':1,'hc_term_count':1}
    terms = [s for v in ('y1','y2') for key in ('original_terms_inputform','hc_terms_inputform')
             for s in yuks[v][key]]
    return {'format':'t3-b1-native-restricted-action-v1',
            'status':'NativeActionTermsIsolated_PhysicalDictionaryPending',
            'model':model, 'tests':{'action_success':True,'four_flags_zero':True},
            'effective_flags':{'T3_EXPERIMENTAL_GENERAL_FERMIONS':'0',
                               'T3_EXPERIMENTAL_GENERAL_QUADRATICS':'0',
                               'T3_EXPERIMENTAL_GENERAL_QUARTICS':'0',
                               'T3_EXPERIMENTAL_FLAVOR_MASSES':'0'},
            'physical_source_certified':False,
            'source_gammaCC_concrete_matrix_measured':False,
            'cg_component_coefficients_measured_here':False,
            'mass_term_count':1,'mass_terms_inputform':[mass], 'yukawas':yuks,
            'source_action_LInt_inputform':' + '.join(terms),
            'source_action_LFree_inputform':mass}


def prior16():
    return {'status':'All16ConditionalCandidateWeylTensorsReconstructed_PhysicalDictionaryPending',
            'count_conditional_comparison_passed':16, 'native_action_observed':True,
            'physical_source_certified':False,
            'models':[{'model':s,'physically_verified':False,
                       'residuals':{'mass':0.0,'Y_after_global_minus':1e-16,'Z_after_global_minus':0.0}}
                      for s,_ in model_specs()]}


class NativeContractTest(unittest.TestCase):
    def test_A_dirac_original_and_hc(self):
        out = verify_record(example('A_p0'), 'A_p0')
        self.assertEqual(out['couplings_and_chiral_projectors']['mass']['coefficient'], '-1/1')
        self.assertTrue(out['couplings_and_chiral_projectors']['y1']['original']['gamma_cc'])
        self.assertFalse(out['couplings_and_chiral_projectors']['y2']['original']['gamma_cc'])

    def test_B_majorana_half_and_projectors(self):
        out = verify_record(example('B_m1'), 'B_m1')
        self.assertEqual(out['couplings_and_chiral_projectors']['mass']['coefficient'], '-1/2')
        self.assertTrue(out['couplings_and_chiral_projectors']['y2']['hc']['gamma_cc'])

    def test_C_triplet_same_chirality_contract(self):
        out = verify_record(example('C_m1'), 'C_m1')
        self.assertEqual(out['couplings_and_chiral_projectors']['mass']['coefficient'], '-1/2')

    def test_wrong_dirac_mass_half_rejected(self):
        src = example('A_p0')['mass_terms_inputform'][0].replace('-(', '-1/2*(')
        with self.assertRaisesRegex(ActionContractError, 'prefactor'):
            verify_term(src,'A_p0','mass')

    def test_wrong_majorana_mass_half_rejected(self):
        src = example('C_m1')['mass_terms_inputform'][0].replace('-1/2*','-1*')
        with self.assertRaisesRegex(ActionContractError, 'prefactor'):
            verify_term(src,'C_m1','mass')

    def test_wrong_yukawa_factor_two_rejected(self):
        src = example('A_p0')['yukawas']['y1']['original_terms_inputform'][0]
        with self.assertRaisesRegex(ActionContractError, 'prefactor'):
            verify_term('2*'+src,'A_p0','y1')

    def test_yukawa_wrong_complex_conjugation_rejected(self):
        src = example('A_p0')['yukawas']['y1']['hc_terms_inputform'][0]
        with self.assertRaisesRegex(ActionContractError, 'coupling/CG/scalar'):
            verify_term(src.replace(f'Bar[{cc("y1")}]',cc('y1'),1),
                        'A_p0','y1',True)

    def test_yukawa_wrong_projector_rejected(self):
        src=example('A_p0')['yukawas']['y2']['original_terms_inputform'][0]
        with self.assertRaisesRegex(ActionContractError, 'projector'):
            verify_term(src.replace('Proj[1]','Proj[-1]'), 'A_p0','y2')

    def test_dirac_y1_missing_charge_conjugation_rejected(self):
        src=example('A_p0')['yukawas']['y1']['original_terms_inputform'][0]
        with self.assertRaisesRegex(ActionContractError, 'projector'):
            verify_term(src.replace('GammaCC,',''), 'A_p0','y1')

    def test_majorana_hc_wrong_C_rejected(self):
        src=example('B_m1')['yukawas']['y2']['hc_terms_inputform'][0]
        with self.assertRaisesRegex(ActionContractError, 'projector'):
            verify_term(src.replace('GammaCC,',''), 'B_m1','y2',True)

    def test_field_swap_rejected(self):
        src=example('B_m1')['yukawas']['y2']['original_terms_inputform'][0]
        with self.assertRaisesRegex(ActionContractError, 'field ordering'):
            verify_term(src.replace('Bar['+L+']', 'Bar['+F+']'), 'B_m1','y2')

    def test_wrong_scalar_HC_orientation_rejected(self):
        src=example('C_m1')['yukawas']['y2']['hc_terms_inputform'][0]
        with self.assertRaisesRegex(ActionContractError, 'coupling/CG/scalar'):
            verify_term(src.replace('*'+S2+'*',f'*Bar[{S2}]*'), 'C_m1','y2',True)

    def test_false_native_certificate_rejected(self):
        rec=example('A_p0'); rec['physical_source_certified']=True
        with self.assertRaisesRegex(ActionContractError, 'certificate'):
            verify_record(rec,'A_p0')

    def test_fake_native_flags_rejected(self):
        rec=example('A_p0');rec['effective_flags']['T3_EXPERIMENTAL_GENERAL_QUADRATICS']='1'
        with self.assertRaisesRegex(ActionContractError, 'restricted'):
            verify_record(rec,'A_p0')

    def test_reject_removed_hc_or_duplicate(self):
        rec=example('A_p0');rec['yukawas']['y1']['hc_term_count']=2
        with self.assertRaisesRegex(ActionContractError, 'count'):
            verify_record(rec,'A_p0')

    def test_forged_term_not_in_full_action_rejected(self):
        rec=example('B_m1');rec['source_action_LInt_inputform']='no terms'
        with self.assertRaisesRegex(ActionContractError, 'full LInt'):
            verify_record(rec,'B_m1')

    def test_unexpected_operators_rejected(self):
        with self.assertRaisesRegex(ActionContractError, 'Unsupported'):
            multiplicative(parse('Coupling[y1,{},0]+NCM[x,y]'))
        with self.assertRaisesRegex(ActionContractError, 'Unsupported|denominator'):
            multiplicative(parse('Coupling[y1,{},0]/h'))

    def test_native_mathematica_index_nested_parser(self):
        s='NCM[Transp[Field[NewFermion,Fermion,{Index[d$$1,T3BSMd3],Index[d$$2,NFlavor]},{}]],DiracProduct[GammaCC],Field[NewFermion,Fermion,{Index[d$$1,T3BSMd3],Index[d$$2,NFlavor]},{}]]'
        r=parse(s)
        self.assertEqual(r.name,'NCM')
        self.assertEqual(len(r.args),3)

    def test_conditional_all16_success_and_residual_bound(self):
        r=check_all16(prior16())
        self.assertEqual(r['conditional_models'],16)
        self.assertLess(r['max_Y_relative_residual'],2e-9)

    def test_conditional_all16_wrong_physics_promotion_rejected(self):
        v=prior16();v['physical_source_certified']=True
        with self.assertRaisesRegex(ActionContractError, 'CONDITIONAL'):
            check_all16(v)

    def test_conditional_all16_wrong_residual_rejected(self):
        v=prior16();v['models'][0]['residuals']['mass']=0.1
        with self.assertRaisesRegex(ActionContractError, 'invalid residual'):
            check_all16(v)

    def test_cc_abstract_not_promoted(self):
        rec={'format':'t3-b1-native-cc-abstract-v2',
             'status':'NativeAbstractCCAlgebraVerified_NumericOrientationPending',
             'model':'A_p0','checks':{x:True for x in range(6)},
             'installed_symbol_manifest':{'symbol_full_name':'Matchete`GammaCC'},
             'chosen_numeric_gamma_matrix_test_performed_here':False,
             'installed_numeric_gamma_cc_sign_measured':False,
             'cross_package_spinor_phase_certified':False,
             'physical_source_certified':False}
        self.assertEqual(check_cc(rec)['abstract_native_cc_checks'],6)
        rec['installed_numeric_gamma_cc_sign_measured']=True
        with self.assertRaisesRegex(ActionContractError,'cannot claim'):
            check_cc(rec)


    def test_complete_runner_binds_three_native_hashes(self):
        with tempfile.TemporaryDirectory() as td:
            root = Path(td); native = root/'native'; native.mkdir()
            priorpath = root/'prior.json'
            refs = {}
            for m in ('A_p0','B_m1','C_m1'):
                p = native/(m+'.json')
                p.write_text(json.dumps(example(m)), encoding='utf-8')
                refs[m] = {'sha256':hashlib.sha256(p.read_bytes()).hexdigest()}
            prior = prior16(); prior['native'] = {'representatives':refs}
            priorpath.write_text(json.dumps(prior), encoding='utf-8')
            with patch('scripts.RunPhysicalB1NativeCoefficientAudit.guard_source', return_value={'stub':{'github_main_blob_match':True}}):
                r = build(root, native, priorpath, None)
                self.assertEqual(r['native_actions_structurally_verified_count'], 3)
                self.assertEqual(r['previous_conditional_16']['conditional_models'],16)
                self.assertFalse(r['physical_source_certified'])
                self.assertFalse(r['safe_to_integrate_into_production'])
                self.assertEqual(r['native_cc_abstract_probe']['status'],'NOT_EXECUTED')
                # An even syntactically valid, but different, native file must block.
                q = example('B_m1'); q['source_action_LFree_inputform'] += ' '
                (native/'B_m1.json').write_text(json.dumps(q), encoding='utf-8')
                with self.assertRaisesRegex(ActionContractError, 'changed since'):
                    build(root,native,priorpath,None)

    def test_source_term_with_extra_factor_rejected(self):
        src=example('A_p0')['yukawas']['y2']['hc_terms_inputform'][0]
        with self.assertRaisesRegex(ActionContractError,'extra factors'):
            verify_term(src+'*HoldForm[h]','A_p0','y2',True)

if __name__=='__main__':
    unittest.main()
