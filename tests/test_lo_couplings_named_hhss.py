"""Standalone exact parser / algebraic source bridge regressions."""
from __future__ import annotations

import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace
import sympy as sp

from Numerical.lo_certification.couplings.named_hhss_projection import (
    generated_declaration, projection_report, named_hhss_sources,
    _literal_wl, _cg3, _tri, _bracketed, SourceInvariantError)
from Numerical.lo_certification.couplings.named_hhss_batch import (
    _matchete_name_evidence, run_named_hhss_batch, _model_slugs)


def decl(name,fields,expr):
    src=(f'AddQuartic[ {name}, {{{", ".join(fields)}}}, '
         f'GroupInvariant -> (({expr}) &), SelfConjugate -> False];')
    return generated_declaration(src,name)


class TestNamedHHSSSourceProofs(unittest.TestCase):
    def test_whole_wl_declaration_is_extracted_exactly(self):
        d=decl('Q',['H','Bar @ H','Bar @ S1','Bar @ S1'],
               'del[SU2L @ fund, #2, #1]')
        self.assertEqual(d['self_conjugate_option'],'False')
        self.assertEqual(d['group_invariant_expression'],'del[SU2L @ fund, #2, #1]')

    def test_simple_source_projection_not_certified_physical(self):
        d=decl('Q',['H','Bar @ H','Bar @ S1','Bar @ S1'],
               'del[SU2L @ fund, #2, #1]')
        r=projection_report('A_+0',(1,1,0,2,0,0),d)
        self.assertEqual(r['reference_dimension'],1)
        self.assertEqual(r['valid_canonical_invariant_orientation_count'],2)
        self.assertIn(['1'],[v['source_to_reference_polynomial_coefficients']
                                    for v in r['candidate_orientation_projection_maps']])
        self.assertIsNone(r['physical_Matchete_RGBeta_coupling_dictionary'])
        self.assertFalse(r['safe_physical_input'])

    def test_triplet_symmetric_cg_and_literal_triple(self):
        self.assertEqual(_cg3(0,0,0),1)
        self.assertEqual(sp.simplify(_cg3(1,0,1)-1/sp.sqrt(2)),0)
        self.assertEqual(_cg3(2,1,1),1)
        self.assertEqual(_cg3(1,0,0),0)
        self.assertEqual(sp.simplify(_tri(0,0,0)),0)

    def test_A_triplet_invariant_dual_map_ambiguity_retained(self):
        d=decl('Q',['Bar @ H','Bar @ H','Bar @ S1','Bar @ S2'],
               'delS2[SU2L, #4, #1, #2]')
        r=projection_report('A_-2',(0,2,0,1,0,1),d)
        self.assertEqual(r['valid_canonical_invariant_orientation_count'],2)
        self.assertIn(['BarS2'],[v['dual_map_field_types'] for v in r['candidate_orientation_projection_maps']])

    def test_E_triplet_compound_invariant_source_coefficient_exact(self):
        d=decl('Q',['Bar @ H','Bar @ H','Bar @ S1','Bar @ S2'],
               'delS2[SU2L, a, #1, #2] * T3RGBetaS2TripleInvariant[a, #3, #4]')
        r=projection_report('E_-2',(0,2,0,1,0,1),d)
        self.assertEqual(r['valid_canonical_invariant_orientation_count'],2)
        self.assertTrue(all(v['source_to_reference_polynomial_coefficients']==['-sqrt(2)/2']
                            for v in r['candidate_orientation_projection_maps']))

    def test_B_mixed_invariant_two_independent_exact_directions(self):
        d=decl('Q',['H','Bar @ H','Bar @ S1','Bar @ S2'],
               'del[SU2L @ fund, #4, #1] * eps[SU2L @ fund, #2, #3]')
        r=projection_report('B_-1',(1,1,0,1,0,1),d)
        self.assertEqual(r['reference_dimension'],2)
        self.assertGreater(r['valid_canonical_invariant_orientation_count'],0)
        self.assertTrue(all(len(v['source_to_reference_polynomial_coefficients'])==2
                            for v in r['candidate_orientation_projection_maps']))

    def test_B_repeated_doublet_eps_symmetrization(self):
        d=decl('Q',['Bar @ H','Bar @ H','Bar @ S2','Bar @ S2'],
           '(eps[SU2L @ fund, #1, #3]*eps[SU2L @ fund, #2, #4] + eps[SU2L @ fund, #1, #4]*eps[SU2L @ fund, #2, #3])/2')
        r=projection_report('B_-3',(0,2,0,0,0,2),d)
        self.assertEqual(r['reference_dimension'],1)
        self.assertEqual(r['valid_canonical_invariant_orientation_count'],2)

    def test_rejects_unsupported_wolfram_and_unsafe_syntax(self):
        for expr in ['tGen[SU2L @ S2, a, #1, #2]',
                     '__import__("os").system("whoami")',
                     'del[SU3, #1, #2]']:
            with self.subTest(expr=expr):
                with self.assertRaises((SourceInvariantError,ValueError)):
                    from Numerical.lo_certification.couplings.named_hhss_projection import _evaluate_raw
                    _evaluate_raw(expr,(0,0,0,0))

    def test_mismatched_field_multiset_fails(self):
        d=decl('Q',['H','Bar @ H','Bar @ S1','Bar @ S1'],'del[SU2L @ fund, #1, #2]')
        with self.assertRaisesRegex(SourceInvariantError,'field content'):
            projection_report('A_+0',(0,2,0,2,0,0),d)

    def test_zero_source_is_not_a_normalization(self):
        d=decl('Q',['H','Bar @ H','Bar @ S1','Bar @ S1'],'0')
        with self.assertRaisesRegex(SourceInvariantError,'identically zero'):
            projection_report('A_+0',(1,1,0,2,0,0),d)

    def test_real_wl_and_registry_must_concord(self):
        with tempfile.TemporaryDirectory() as td:
            root=Path(td);base=root/'RGE/general/generated';base.mkdir(parents=True)
            core=root/'RGE/running/rgbeta';core.mkdir(parents=True)
            (core/'T3RGBetaModel.wl').write_text('(* unused *)',encoding='utf8')
            (base/'A_+0_extra_quartics.wl').write_text(
                'AddQuartic[Q, {H, Bar @ H, Bar @ S1, Bar @ S1}, '
                'GroupInvariant -> (del[SU2L @ fund, #2, #1] &), SelfConjugate -> False];',encoding='utf8')
            reg={'new_rgbeta_symbols':[{'name':'Q','counts':[1,1,0,2,0,0],
                        'invariant':'del[SU2L @ fund, #1, #2]'}],
                 'preexisting_rgbeta_symbols':[]}
            with self.assertRaisesRegex(SourceInvariantError,'differs from registry'):
                named_hhss_sources(root,'A_+0',reg)
            reg['new_rgbeta_symbols'][0]['invariant']='del[SU2L @ fund, #2, #1]'
            r=named_hhss_sources(root,'A_+0',reg)
            self.assertEqual(len(r['sectors']),1)

    def test_native_seed_provenance_is_not_a_mapping(self):
        seed={'ScalarQuarticTerms':[{'Index':2,'TermInputForm':
             'Coupling[lambdaT3ExtraQ110200I1, {}, 0]', 'CGNames':['CGA']} ]}
        out=_matchete_name_evidence(seed,{'matchete_basis_name':'lambdaT3ExtraQ110200I#'},'new')
        self.assertEqual(len(out['matching_native_terms']),1)
        self.assertFalse(out['verified_named_CG_to_RGBeta_map'])

    def test_outside_owned_output_rejected_first(self):
        with self.assertRaisesRegex(ValueError,'ownership'):
            run_named_hhss_batch(Path('/not-a-repo'),Path('/tmp/outside'),('A_-4',))

    def test_slug_paths_match_prior_v3_runs(self):
        self.assertEqual(_model_slugs('A_-4'),'A_m4')
        self.assertEqual(_model_slugs('E_+2'),'E_p2')


if __name__=='__main__':unittest.main()
