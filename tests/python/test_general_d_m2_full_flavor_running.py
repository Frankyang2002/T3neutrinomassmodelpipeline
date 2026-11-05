"""Isolated regression tests for the full-flavour D(-2) direct EFT1 insertion."""
from __future__ import annotations
import json
from pathlib import Path
import sys
import tempfile
import unittest

ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT))
from scripts.RunGeneralDMinus2FullFlavorRunning import (
    beta_pair_coefficients, pl_pair_inventory, wolfram_insertion,
    wolfram_kernel, format_kernel, calculate
)
import sympy as sp

BETA = ('4*sqrt(3)*T3S2Phase*(lambdaT3*conjugate(y2) + '
        'conjugate(lambdaT3ExtraQ020101I1)*conjugate(y2Conj))*'
        'conjugate(T3S2Phase)*conjugate(y1)/(3*MF)')


def fixture_seed():
    # Matchete's actual exported four PL structures may differ in field/CG
    # placement; our parser validates the exact indexed Yukawa contraction.
    terms=[]
    for k,(a,b) in enumerate((('y1','y2'),('y1','y2'),
                               ('y1','y2Conj'),('y1','y2Conj')),1):
        s=(f'Bar[Coupling[{a}, {{Index[p$, Flavor], Index[r$, NFlavor]}}, 0]]*'
           f'Bar[Coupling[{b}, {{Index[q$, Flavor], Index[r$, NFlavor]}}, 0]]/'
           'Coupling[MF, {}, 0]')
        terms.append({'Index':k,'Chirality':'PL','TermInputForm':s})
    for k in range(5,9):
        terms.append({'Index':k,'Chirality':'PR','TermInputForm':'placeholder'})
    return {'Model':'D_-2','Status':'Success',
            'PhysicalScalarNames':['NewScalar1','T3S2R','T3S2I'],
            'TreeWilsonTerms':terms,'TreeWilsonTermCount':8}


class TestD2FullFlavor(unittest.TestCase):
    def test_exact_beta_pair_extraction(self):
        d=beta_pair_coefficients(BETA)
        self.assertEqual(set(d),{('y1','y2'),('y1','y2Conj')})
        self.assertEqual(sp.simplify(d[('y1','y2')]-4*sp.sqrt(3)*sp.Symbol('lambdaT3')/3),0)
        self.assertEqual(sp.simplify(d[('y1','y2Conj')]-4*sp.sqrt(3)*sp.conjugate(sp.Symbol('lambdaT3ExtraQ020101I1'))/3),0)

    def test_neutral_phase_is_removed_using_unitarity(self):
        for v in beta_pair_coefficients(BETA).values():
            self.assertNotIn('T3S2Phase',str(v))

    def test_unexpected_yukawa_monomial_fails(self):
        with self.assertRaisesRegex(ValueError,'Unexpected Yukawa monomials'):
            beta_pair_coefficients('conjugate(y1)**2/MF')

    def test_tree_flavor_seed_has_both_mixed_pairs(self):
        x=pl_pair_inventory(fixture_seed())
        self.assertEqual(x['pl_terms'],4)
        self.assertEqual({tuple(p) for p in x['unique_yukawa_pairs']},
                         {('y1','y2'),('y1','y2Conj')})

    def test_missing_second_yukawa_pair_fails(self):
        s=fixture_seed()
        for term in s['TreeWilsonTerms'][:4]:
            term['TermInputForm']=term['TermInputForm'].replace('y2Conj','y2')
        with self.assertRaisesRegex(ValueError,'Missing full-flavour'):
            pl_pair_inventory(s)

    def test_unbarred_factor_fails(self):
        s=fixture_seed()
        s['TreeWilsonTerms'][0]['TermInputForm']=s['TreeWilsonTerms'][0]['TermInputForm'].replace('Bar[Coupling[y1, {Index[p$, Flavor], Index[r$, NFlavor]}, 0]]','Coupling[y1, {Index[p$, Flavor], Index[r$, NFlavor]}, 0]')
        with self.assertRaisesRegex(ValueError,'expected two barred'):
            pl_pair_inventory(s)

    def test_unshared_heavy_flavor_fails(self):
        s=fixture_seed()
        s['TreeWilsonTerms'][0]['TermInputForm']=s['TreeWilsonTerms'][0]['TermInputForm'].replace('Index[q$, Flavor], Index[r$, NFlavor]','Index[q$, Flavor], Index[s$, NFlavor]')
        with self.assertRaisesRegex(ValueError,'inconsistent heavy'):
            pl_pair_inventory(s)

    def test_kernel_external_exchange_symmetry(self):
        s=format_kernel('y1','y2')
        self.assertIn('conjugate(y1[p,r])*conjugate(y2[q,r])',s)
        self.assertIn('conjugate(y2[p,r])*conjugate(y1[q,r])',s)
        self.assertIn('/2',s)

    def test_wolfram_insertion_safe_and_independent(self):
        wl=wolfram_insertion(beta_pair_coefficients(BETA),'MF','MuScalar')
        self.assertIn('EFT1D2C12[p_, q_]',wl)
        self.assertIn('EFT1D2C1Conj[p_, q_]',wl)
        self.assertIn('Bar[Coupling[lambdaT3ExtraQ020101I1, {}, 0]]',wl)
        self.assertIn('Bar[Coupling[y2Conj, {q, Index[r, NFlavor]}, 0]]',wl)
        self.assertIn('Log[Coupling[MuScalar, {}, 0]/Coupling[MF, {}, 0]]',wl)
        self.assertIn('EFT1D2WilsonRunningHeavyInsertion = 0',wl)
        self.assertNotIn('hbar *',wl)

    def test_wolfram_invalid_scale_rejected(self):
        with self.assertRaisesRegex(ValueError,'Scale symbols'):
            wolfram_insertion(beta_pair_coefficients(BETA),'MF','foo/bar')

    def test_hard_threshold_not_claimed(self):
        with tempfile.TemporaryDirectory() as tmp:
            base=Path(tmp)
            (base/'data/general_d_m2_eft1_rge').mkdir(parents=True)
            (base/'data/general_d_m2_stage1_summary.json').write_text(json.dumps({'status':'GeneralDMinus2RealMassBasisFermionThresholdVerified'}))
            (base/'data/eft1_after_F_realbasis_wilson_tensor_seed.json').write_text(json.dumps(fixture_seed()))
            (base/'data/general_d_m2_eft1_rge/weinberg_mixing_beta.json').write_text(json.dumps({
                'status':'Success','checked_LLHH_components':30,
                'nonzero_LLHH_beta_components':10,
                'weinberg_subspace_validation':{
                    'matches_weinberg_subspace':True,'failure_count':0,
                    'beta_kappa_16pi2':BETA}}))
            (base/'scalar_mass_basis.json').write_text(json.dumps({
                'RealFields':['T3S2R','T3S2I'],
                'RealMassSquared':'MS2^2-2*Abs[kappa]',
                'ImagMassSquared':'MS2^2+2*Abs[kappa]'}))
            for stem in ('tree','loop','transition_tree','transition_full'):
                (base/f'data/eft1_after_F_realbasis_{stem}.wxf').write_bytes(b'fixture')
            result=calculate(base)
            self.assertEqual(result['status'],'GeneralDMinus2FullFlavorDirectRunningPrepared')
            self.assertTrue(result['full_flavor_direct_running_performed'])
            self.assertFalse(result['scalar_hard_threshold_performed'])
            self.assertFalse(result['full_wilson_self_rge_performed'])
            output=json.loads((base/'data/general_d_m2_full_flavor/full_flavor_direct_running.json').read_text())
            self.assertTrue(output['one_generation_reduction_matches'])
            self.assertIn('y2Conj', output['full_flavor_beta_tensor'])
            preflight=json.loads((base/'data/general_d_m2_full_flavor/scalar_threshold_preflight.json').read_text())
            self.assertEqual(len(preflight['input_eft1_continuation_files']),4)
            self.assertFalse(preflight['scalar_hard_threshold_performed'])

    def test_missing_continuation_fails_closed(self):
        with tempfile.TemporaryDirectory() as tmp:
            result=calculate(Path(tmp))
            self.assertEqual(result['status'],'Failed')
            self.assertIn('Missing required',result['error'])


if __name__=='__main__':unittest.main()
