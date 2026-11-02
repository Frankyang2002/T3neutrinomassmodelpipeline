"""Formal SU(2) tensor and diagnostic-only quartic bridge regression tests."""
from __future__ import annotations

from copy import deepcopy
from itertools import permutations
import json
from pathlib import Path
import tempfile
import unittest

import numpy as np

from RGE.running.intermediate.FullSU2SourceResponse import (
    QUARTIC_FORMAT, SourceKernelError, _pair_blocks, _contract_pair,
    project_weinberg, read_diagnostic_quartic, read_source_blocks,
    run_diagnostic, weinberg_template,
)
from RGE.running.intermediate.ProvisionalSourceKernel import SOURCE_FORMAT, _orbit, _parse_coordinates
from RGE.running.intermediate.GeneralEFT1C5Mixing import Context, calculate_master_beta
from scripts import RunFullSU2SourceResponse


def row(i,p,j,q,a,b,z):
    k=(i,p,j,q,a,b)
    assert k==min(_orbit(k))
    return {'coordinate':list(k),'orbit_size':len(_orbit(k)),
            'coefficient_real':z.real,'coefficient_imag':z.imag}


def boundary():
    r=row(1,1,2,2,5,6,3e-11+2e-11j)
    return {'format':SOURCE_FORMAT, 'model':'A_-4',
            'status':'NumericalLLSSBySourceVerified_SpinorNormalizationOpen',
            'mu_F_gev':100000.0,'units':'GeV^-1',
            'index_base':1,
            'index_order':['lepton_su2_i','lepton_flavour_p','lepton_su2_j',
                           'lepton_flavour_q','real_scalar_a','real_scalar_b'],
            'real_scalar_dimension':6,
            'scalar_field_offsets_one_based':{'X':5},
            'scalar_field_complex_dimensions':{'X':1},
            'spinor_normalization':{'certified':False,
                'gamma_cc_to_master_weyl_factor':None,
                'gamma_cc_to_master_weyl_sign':None,
                'do_not_use_as_physical_master_wilson_boundary':True},
            'completeness':{'F_tree_LLSS_included':True,
                'distinct_PL_sources_retained':True,
                'PR_Hermitian_conjugate_not_double_counted':True,
                'complete_physical_leading_one_loop_neutrino_mass_certified':False},
            'sources':[{'matchete_source_index':37,'native_chirality':'PL',
                        'zero_at_this_benchmark':False,
                        'independent_symmetric_coordinates':[r]}],
            'summed_independent_symmetric_coordinates':[deepcopy(r)]}


def quartic():
    return {'format':QUARTIC_FORMAT,'model':'A_-4','authority':'diagnostic_only',
            'mu_F_gev':100000.0,'source_boundary_sha256':'fixture-sha',
            'coverage':'partial_diagnostic',
            'quartic_components':[{'indices':[3,3,5,6],'value':0.375}]}


def decode(raw):
    return np.asarray(raw['real'])+1j*np.asarray(raw['imag'])


class SU2TemplateTests(unittest.TestCase):
    def test_normalization_with_half_c5_and_complex_higgs(self):
        t=weinberg_template()
        self.assertEqual(t.shape,(2,2,4,4))
        self.assertAlmostEqual(t[0,0,2,2].real, .25)
        self.assertAlmostEqual(t[0,0,2,3].imag, .25)
        self.assertAlmostEqual(t[0,0,3,3].real,-.25)
        self.assertAlmostEqual(t[1,1,0,0].real,.25)
        self.assertAlmostEqual(t[0,1,2,0].real,-.125)
        self.assertAlmostEqual(t[0,1,0,2].real,-.125)
        self.assertTrue(np.allclose(t,t.transpose(1,0,2,3)))
        self.assertTrue(np.allclose(t,t.transpose(0,1,3,2)))

    def test_full_weinberg_projection_and_pivot_agree_for_exact_tensor(self):
        m=np.array([[2+1j,0.1+2j,-1+0.2j],
                    [0.1+2j,3-2j,2.5+2j],
                    [-1+0.2j,2.5+2j,4-4j]])*1e-12
        beta=np.einsum('pq,ijAB->pqijAB',m,weinberg_template())
        d=project_weinberg(beta)
        self.assertTrue(d['SU2_Ward_tensor_algebra_pass'])
        self.assertTrue(d['C5_flavour_symmetric'])
        np.testing.assert_allclose(decode(d['projected_C5_beta_16pi2_per_NGamma']),m,rtol=1e-12,atol=1e-28)
        np.testing.assert_allclose(decode(d['neutrino_H0R_pivot_beta_16pi2_per_NGamma']),m,rtol=1e-12,atol=1e-28)

    def test_only_neutrino_pivot_can_miss_su2_noninvariance(self):
        m=np.zeros((3,3,2,2,4,4),complex)
        m[0,0,0,0,2,2]=1e-11
        d=project_weinberg(m)
        self.assertFalse(d['SU2_Ward_tensor_algebra_pass'])
        self.assertGreater(d['pivot_vs_projection_abs'],0)
        self.assertGreater(d['weinberg_residual_relative'],.1)


class FullSourceTests(unittest.TestCase):
    def test_complex_full_lepton_Su2_blocks_and_scalar_pair_symmetry(self):
        blocks,ns=read_source_blocks(boundary())
        self.assertEqual(ns,6)
        self.assertEqual(len(blocks),1)
        z=blocks[0][1][(5,6)]
        self.assertEqual(z.shape,(3,3,2,2))
        self.assertAlmostEqual(z[0,1,0,1].real,3e-11)
        self.assertAlmostEqual(z[1,0,1,0].imag,2e-11)
        self.assertEqual(z[0,1,1,0],0)

    def test_diagnostic_wo_quartic_is_factorised_not_beta(self):
        out=run_diagnostic(boundary(),source_sha256='fixture-sha')
        self.assertEqual(out['status'],'SourceFactorisedFullSU2_NoQuartic')
        self.assertNotIn('total',out)
        self.assertIsNone(out['N_GammaCC_to_Weyl'])
        self.assertFalse(out['physical_C5_computed'])
        raw=out['source_resolved_blocks'][0]['scalar_pairs'][0]
        self.assertEqual(raw['heavy_pair_one_based'],[5,6])
        self.assertEqual(raw['ordered_pair_multiplicity'],2)
        self.assertEqual(decode(raw['LLSS_C_pqij_per_NGamma_GeV_inv'])[0,1,0,1],3e-11+2e-11j)

    def test_quartic_matches_exact_master_RGE_component(self):
        b=boundary()
        src=_parse_coordinates(b['sources'][0]['independent_symmetric_coordinates'],6,'fixture')
        c=np.zeros((6,6,6,6),complex)
        for (i,p,j,q,a,s),z in src.items():
            c[2*(p-1)+(i-1),2*(q-1)+(j-1),a-1,s-1]=z
        lam=np.zeros((6,6,6,6),float)
        for k in set(permutations((2,2,4,5))):
            lam[k]=0.375
        ctx=Context('A_-4',6,6,((0,1),(2,3),(4,5)),(0,1,2,3),
                    100000.,25000.,np.zeros((6,6,6),complex),lam,(),
                    (('tree',c),),'synthetic',{})
        expected=calculate_master_beta(ctx,c)['scalar_pair']
        out=run_diagnostic(b,source_sha256='fixture-sha',quartic=quartic(),mu_s_gev=25000)
        self.assertFalse(out['physical_C5_computed'])
        self.assertAlmostEqual(out['log_muS_over_muF'],np.log(.25))
        # The new contraction is unprojected, full SU2, and must match beta master.
        blocks,_=read_source_blocks(b)
        lam_hh,_=read_diagnostic_quartic(quartic(),b,'fixture-sha')
        got=_contract_pair(blocks[0][1][(5,6)],lam_hh,5,6)
        for p in range(3):
            for q in range(3):
                for i in range(2):
                    for j in range(2):
                        for A in range(4):
                            for B in range(4):
                                self.assertAlmostEqual(
                                    abs(got[p,q,i,j,A,B]-
                                    expected[2*p+i,2*q+j,A,B]),0,delta=1e-25)
        self.assertFalse(out['total']['SU2_Ward_tensor_algebra_pass'])
        summed=decode(out['per_source'][0]['conditional_full_SU2_Weinberg']['projected_C5_beta_16pi2_per_NGamma'])
        np.testing.assert_allclose(summed,decode(out['total']['projected_C5_beta_16pi2_per_NGamma']),atol=1e-26)
        self.assertEqual(len(out['per_source'][0]['quartic_contributions'][0]['pair_contributions']),1)

    def test_independent_quartic_sources_same_component_add_without_overwrite(self):
        q=quartic()
        q.pop('quartic_components')
        q['quartic_terms']=[
            {'label':'operator_one', 'components':[{'indices':[3,3,5,6],'value':.125}]},
            {'label':'operator_two', 'components':[{'indices':[6,3,5,3],'value':.25}]},
        ]
        out=run_diagnostic(boundary(),source_sha256='fixture-sha',quartic=q)
        self.assertEqual(out['quartic_input']['quartic_channel_count'],2)
        self.assertEqual([x['quartic_label'] for x in out['per_quartic']],
                         ['operator_one','operator_two'])
        sums=[decode(z['conditional_full_SU2_Weinberg']['projected_C5_beta_16pi2_per_NGamma'])
              for z in out['per_quartic']]
        np.testing.assert_allclose(sum(sums),
                                  decode(out['total']['projected_C5_beta_16pi2_per_NGamma']),
                                  rtol=1e-13,atol=1e-27)
        self.assertEqual(len(out['per_source'][0]['quartic_contributions']),2)

    def test_missing_conflicting_wrong_origin_quartic_fail_closed(self):
        b=boundary()
        for action in ('sha','authority','complex','duplicate','wrong_fields'):
            q=quartic()
            if action=='sha':q['source_boundary_sha256']='wrong'
            elif action=='authority':q['authority']='certified'
            elif action=='complex':q['quartic_components'][0]['value']=1j
            elif action=='duplicate':q['quartic_components'].append({'indices':[6,5,3,3],'value':.5})
            elif action=='wrong_fields':q['quartic_components'][0]['indices']=[3,3,4,5]
            with self.subTest(action=action):
                with self.assertRaises(SourceKernelError):
                    run_diagnostic(b,source_sha256='fixture-sha',quartic=q)

    def test_changed_spinor_certification_rejected(self):
        b=boundary()
        b['spinor_normalization']['certified']=True
        with self.assertRaisesRegex(SourceKernelError,'normalization'):
            run_diagnostic(b,source_sha256='fixture-sha')

    def test_source_sum_rechecked_not_just_pivot(self):
        b=boundary()
        b['summed_independent_symmetric_coordinates'][0]['coefficient_real']+=.1
        with self.assertRaisesRegex(SourceKernelError,'sum mismatch'):
            run_diagnostic(b,source_sha256='fixture-sha')

    def test_cli_check_only_and_never_overwrite(self):
        with tempfile.TemporaryDirectory() as d:
            root=Path(d)
            release=root/'rel'; (release/'A_m4').mkdir(parents=True)
            (release/'A_m4/boundary.json').write_text(json.dumps(boundary()))
            summary={'absolute_spinor_normalization_certified':False,
                     'physical_eft1_workers_must_reject_this_as_certified_input':True,
                     'models':[{'model':'A_-4','status':'NumericalLLSSBySourceVerified_SpinorNormalizationOpen',
                                'PL_sources':1,'path':'A_m4/boundary.json'}]}
            (release/'summary.json').write_text(json.dumps(summary))
            output=root/'output/general_observable/lo_eft1_mixing/test_case'
            old=RunFullSU2SourceResponse.PROJECT
            try:
                RunFullSU2SourceResponse.PROJECT=root
                args=['--release-root',str(release),'--output-root',str(output),'--model','A_-4']
                self.assertEqual(RunFullSU2SourceResponse.main(args+['--check-only']),0)
                self.assertFalse(output.exists())
                self.assertEqual(RunFullSU2SourceResponse.main(args),0)
                self.assertTrue((output/'A_m4/full_su2_response.json').is_file())
                with self.assertRaises(FileExistsError):
                    RunFullSU2SourceResponse.main(args)
            finally:
                RunFullSU2SourceResponse.PROJECT=old


if __name__ == '__main__':
    unittest.main()
