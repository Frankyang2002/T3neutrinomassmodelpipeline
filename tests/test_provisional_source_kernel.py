"""Independent identities for the *conditional*, source-resolved scalar-pair kernel."""
from __future__ import annotations

import copy
import json
from pathlib import Path
import tempfile
from itertools import permutations
import unittest

import numpy as np

from RGE.running.intermediate.ProvisionalSourceKernel import (
    RESULT_FORMAT, SOURCE_FORMAT, SourceKernelError, _parse_coordinates,
    kernel_for_boundary, _orbit,
)
from RGE.running.intermediate.GeneralEFT1C5Mixing import Context, calculate_master_beta
from scripts import RunProvisionalSourceKernel


def row(i, p, j, q, a, b, z):
    k = (i,p,j,q,a,b)
    assert k == min(_orbit(k))
    return {"coordinate": list(k), "coefficient_real": z.real,
            "coefficient_imag": z.imag, "orbit_size": len(_orbit(k))}


def fixture():
    z = 1e-11 + 2e-12j
    # SU2 neutrino, flavour 1,2, heavy scalars 5,6.
    r = row(1,1,1,2,5,6,z)
    return {
        "format":SOURCE_FORMAT,"model":"A_-4",
        "status":"NumericalLLSSBySourceVerified_SpinorNormalizationOpen",
        "mu_F_gev":1e5,"units":"GeV^-1",
        "index_base":1,"index_order":["lepton_su2_i","lepton_flavour_p",
               "lepton_su2_j","lepton_flavour_q","real_scalar_a","real_scalar_b"],
        "real_scalar_dimension":6,"scalar_field_offsets_one_based":{"X":5},
        "scalar_field_complex_dimensions":{"X":1},
        "spinor_normalization":{
            "certified":False,"gamma_cc_to_master_weyl_factor":None,
            "gamma_cc_to_master_weyl_sign":None,
            "do_not_use_as_physical_master_wilson_boundary":True},
        "completeness":{"F_tree_LLSS_included":True,
            "distinct_PL_sources_retained":True,"PR_Hermitian_conjugate_not_double_counted":True,
            "complete_physical_leading_one_loop_neutrino_mass_certified":False},
        "sources":[{"matchete_source_index":17,"native_chirality":"PL",
                    "zero_at_this_benchmark":False,
                    "independent_symmetric_coordinates":[r]}],
        "summed_independent_symmetric_coordinates":[copy.deepcopy(r)],
    }


class SourceKernelTests(unittest.TestCase):
    def test_exact_distinct_orbit_and_complex_flavour(self):
        data=fixture()
        out=kernel_for_boundary(data)
        self.assertEqual(out['format'],RESULT_FORMAT)
        self.assertEqual(out['status'],'ProvisionalLinearQuarticKernel_NotPhysicalC5')
        self.assertEqual(out['source_count'],1)
        self.assertEqual(out['distinct_nonzero_quartic_pairs'],1)
        channel=out['summed_kernel_channels'][0]
        self.assertEqual(channel['heavy_scalar_pair_one_based'],[5,6])
        self.assertEqual(channel['multiplicity_in_ordered_ab_sum'],2)
        mat=np.array(channel['conditional_beta_C5_per_unit_quartic']['real'])+1j*np.array(channel['conditional_beta_C5_per_unit_quartic']['imag'])
        self.assertAlmostEqual(mat[0,1].real,8e-11)
        self.assertAlmostEqual(mat[0,1].imag,1.6e-11)
        self.assertEqual(mat[1,0],mat[0,1])
        self.assertEqual(mat[0,0],0)
        self.assertIsNone(out['physical_normalization']['N_GammaCC_to_Weyl'])

    def test_one_scalar_diagonal_pair_has_no_double_counting(self):
        data=fixture()
        r=row(1,1,1,1,5,5,3e-11+0j)
        data['sources'][0]['independent_symmetric_coordinates']=[r]
        data['summed_independent_symmetric_coordinates']=[r]
        out=kernel_for_boundary(data)
        item=out['summed_kernel_channels'][0]
        self.assertEqual(item['heavy_scalar_pair_one_based'],[5,5])
        self.assertEqual(item['multiplicity_in_ordered_ab_sum'],1)
        self.assertAlmostEqual(item['conditional_beta_C5_per_unit_quartic']['real'][0][0],1.2e-10)

    def test_zero_source_retained_and_recomposes(self):
        data=fixture()
        data['sources'].append({'matchete_source_index':18,'native_chirality':'PL',
                                'zero_at_this_benchmark':True,
                                'independent_symmetric_coordinates':[]})
        out=kernel_for_boundary(data)
        self.assertEqual(out['source_count'],2)
        self.assertTrue(out['per_source'][1]['zero_at_this_benchmark'])
        self.assertFalse(out['per_source'][1]['kernel_channels'])

    def test_source_sum_tamper_rejected(self):
        data=fixture()
        data['summed_independent_symmetric_coordinates'][0]['coefficient_real']+=1e-8
        with self.assertRaisesRegex(SourceKernelError,'sum mismatch'):
            kernel_for_boundary(data)

    def test_noncanonical_wrong_orbit_or_pr_rejected(self):
        for action in ('orbit','noncanonical','PR'):
            with self.subTest(action=action):
                data=fixture()
                if action=='orbit': data['sources'][0]['independent_symmetric_coordinates'][0]['orbit_size']=2
                elif action=='noncanonical':data['sources'][0]['independent_symmetric_coordinates'][0]['coordinate']=[1,2,1,1,5,6]
                else:data['sources'][0]['native_chirality']='PR'
                with self.assertRaises(SourceKernelError): kernel_for_boundary(data)

    def test_cannot_upgrade_certification_by_flag(self):
        data=fixture()
        data['spinor_normalization']['certified']=True
        with self.assertRaisesRegex(SourceKernelError,'normalization flags changed'):
            kernel_for_boundary(data)

    def test_no_non_quartic_direct_term_for_zero_y_and_block_scalar_pair(self):
        data=fixture()
        source=_parse_coordinates(data['sources'][0]['independent_symmetric_coordinates'],6,'fixture')
        c=np.zeros((6,6,6,6),complex)
        for (i,p,j,q,a,b),z in source.items():
            I=2*(p-1)+(i-1)
            J=2*(q-1)+(j-1)
            c[I,J,a-1,b-1]=z
        lam=np.zeros((6,6,6,6),float)
        for ind in set(permutations((2,2,4,5))):lam[ind]=0.3
        ctx=Context('A_-4',6,6,((0,1),(2,3),(4,5)),(0,1,2,3),
                    100000.0,2000.0,np.zeros((6,6,6),complex),lam,(),
                    (('tree',c),),'synthetic',{})
        terms=calculate_master_beta(ctx,c)
        for name, v in terms.items():
            if name!='scalar_pair':
                self.assertLess(np.max(np.abs(v)),1e-25,name)
        beta=4*terms['scalar_pair'][0,2,2,2]
        item=kernel_for_boundary(data)['summed_kernel_channels'][0]
        K=complex(item['conditional_beta_C5_per_unit_quartic']['real'][0][1],
                  item['conditional_beta_C5_per_unit_quartic']['imag'][0][1])
        np.testing.assert_allclose(beta,K*0.3,rtol=1e-12,atol=1e-25)

    def test_cli_one_model_and_no_overwrite(self):
        with tempfile.TemporaryDirectory() as temp:
            root=Path(temp)
            release=root/'release'
            source=release/'A_m4'/'boundary.json'
            source.parent.mkdir(parents=True)
            source.write_text(json.dumps(fixture()),encoding='utf-8')
            summary={"absolute_spinor_normalization_certified":False,
                     "physical_eft1_workers_must_reject_this_as_certified_input":True,
                     "models":[{"model":"A_-4",
                                "status":"NumericalLLSSBySourceVerified_SpinorNormalizationOpen",
                                "PL_sources":1,"path":"A_m4/boundary.json"}]}
            (release/'summary.json').write_text(json.dumps(summary),encoding='utf-8')
            out=root/'output/general_observable/lo_eft1_mixing/test_case'
            old=RunProvisionalSourceKernel.ROOT
            try:
                RunProvisionalSourceKernel.ROOT=root
                command=['--release-root',str(release),'--output-root',str(out),'--model','A_-4']
                self.assertEqual(RunProvisionalSourceKernel.main(command),0)
                self.assertTrue((out/'A_m4/quartic_pivot_kernel.json').is_file())
                self.assertEqual(json.loads((out/'summary.json').read_text())['model_count'],1)
                with self.assertRaises(FileExistsError):RunProvisionalSourceKernel.main(command)
                with self.assertRaises(SourceKernelError):
                    RunProvisionalSourceKernel.main(command[:3]+[str(root/'outside'),'--model','A_-4'])
            finally:
                RunProvisionalSourceKernel.ROOT=old


if __name__=='__main__':unittest.main()
