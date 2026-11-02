"""Independent source-backed seven-term batch checks (all synthetic in tests)."""
from __future__ import annotations
import copy
import json
from pathlib import Path
import tempfile
import unittest

import numpy as np

from RGE.running.intermediate.ProvisionalSourceKernel import (
    SOURCE_FORMAT, _orbit, MODELS, SourceKernelError, tag,
)
from RGE.running.intermediate.SourceBackedEFT1Batch import (
    COUPLING_FORMAT, CONTRACT_FORMAT, _operator_sources, _merge_quartics,
    _source_layout, _operator_fidelity, required_couplings, evaluate,
    _validate_couplings, stage22_inventory,
)
from RGE.running.intermediate.GeneralEFT1C5Mixing import (
    load_boundary, calculate_master_beta, encode_matrix, TERMS,
)
from scripts import RunSourceBackedEFT1Batch


def row(i,p,j,q,a,b,z):
    k=(i,p,j,q,a,b)
    assert k==min(_orbit(k))
    return {'coordinate':list(k),'coefficient_real':z.real,
            'coefficient_imag':z.imag,'orbit_size':len(_orbit(k))}


def source(model='A_-4', *, mixed=False, split=False):
    # 4 Higgs + 8 heavy real scalars (matching the A/B/C/D source dimension).
    if mixed: names={'NewScalarX':5,'NewScalarY':9}
    elif split: names={'NewScalar1':5,'T3Scalar2Minus':11,'T3Scalar2Plus':12}
    else: names={'NewScalar1':5,'NewScalar2':7}
    r=row(1,1,1,2,5,7, 1.8e-11+3.3e-12j)
    return {
        'format':SOURCE_FORMAT,'model':model,
        'status':'NumericalLLSSBySourceVerified_SpinorNormalizationOpen',
        'mu_F_gev':1e5,'units':'GeV^-1','index_base':1,
        'index_order':['lepton_su2_i','lepton_flavour_p','lepton_su2_j',
                       'lepton_flavour_q','real_scalar_a','real_scalar_b'],
        'real_scalar_dimension':12,
        'scalar_field_offsets_one_based':names,
        'scalar_field_complex_dimensions':{name:1 for name in names},
        'spinor_normalization':{'certified':False,
             'gamma_cc_to_master_weyl_factor':None,
             'gamma_cc_to_master_weyl_sign':None,
             'do_not_use_as_physical_master_wilson_boundary':True},
        'completeness':{'F_tree_LLSS_included':True,
           'distinct_PL_sources_retained':True,
           'PR_Hermitian_conjugate_not_double_counted':True,
           'complete_physical_leading_one_loop_neutrino_mass_certified':False},
        'sources':[
            {'matchete_source_index':17,'native_chirality':'PL',
             'zero_at_this_benchmark':False,
             'independent_symmetric_coordinates':[r]},
            {'matchete_source_index':19,'native_chirality':'PL',
             'zero_at_this_benchmark':True,
             'independent_symmetric_coordinates':[]},
        ],'summed_independent_symmetric_coordinates':[copy.deepcopy(r)],
    }


def diagnostic_couplings(b, *, gauge_zero=True):
    ns=b['real_scalar_dimension']
    # Two independent quartic vertices on the same component orbit must add.
    return {
        'format':COUPLING_FORMAT,'model':b['model'],
        'source_boundary_sha256':'abc', 'renormalization_scheme':'MSbar',
        'master_conventions':{
          'master_wilson':'C_ijab psi_i psi_j phi_a phi_b + hc',
          'weinberg':'(1/2) C5_ij (L_i.H)(L_j.H) + hc',
          'quartic':'V4=lambda_abcd phi_a phi_b phi_c phi_d/4!',
        },
        'n_real_scalars':ns,'n_fermions':6,'mu_F_gev':1e5,'mu_S_gev':1500.,
        'coverage_assertion':{'authority':'diagnostic_only',
           'all_active_weyl_fermion_components_enumerated':True,
           'all_yukawa_components_enumerated':True,
           'all_quartic_channels_enumerated':True,
           'all_gauge_generators_enumerated_or_explicitly_toy_zero':True,
           'toy_gauge_couplings_set_to_zero':gauge_zero,
        },
        'yukawa_components':[
           {'indices':[1,2,5],'value':{'real':0.05,'imag':0.02}},
           {'indices':[3,4,7],'value':{'real':0.03,'imag':-0.04}},
        ],
        'gauge_sectors':[],
        'quartic_channels':[
            {'label':'portal1','components':[{'indices':[3,3,5,7],'value':0.12}]},
            {'label':'portal2','components':[{'indices':[7,5,3,3],'value':0.08}]},
        ]
    }


class SourceConversionTests(unittest.TestCase):
    def test_mixed_XY_offsets_not_guessed(self):
        b=source('B_-1',mixed=True)
        self.assertEqual([len(x['indices']) for x in _source_layout(b)],[4,4,4])
        self.assertEqual([x['name'] for x in _source_layout(b)],['H','NewScalarX','NewScalarY'])

    def test_split_real_D_E_not_automatically_doubled(self):
        b=source('D_-2',split=True)
        self.assertEqual([len(x['indices']) for x in _source_layout(b)],[4,6,1,1])
        self.assertEqual(_source_layout(b)[2]['indices'],[11])

    def test_source_pl_zero_preserved_complex_and_separate(self):
        b=source()
        ops=_operator_sources(b)
        self.assertEqual([x['source_id'] for x in ops],['Matchete_PL_17','Matchete_PL_19'])
        self.assertEqual(ops[1]['components'][0]['value'],{'real':0.0,'imag':0.0})
        self.assertEqual(ops[0]['components'][0]['indices'],[1,3,5,7])
        self.assertEqual(ops[0]['components'][0]['value']['imag'],3.3e-12)

    def test_stage22_all16_46_candidate_blocks_not_physical(self):
        blocks={model:[] for model in MODELS}
        blocks['A_-2']=[{'block':'021010','status':'DirectGeometryCandidate',
                           'physical_mapping_certified':False}]
        report={'format':'stage22-aligned-geometry-diagnostics-v1',
            'status':'ExactConditionalGeometryNotPhysicalCouplingCertification',
            'native_rgbeta_component_basis_verified':False,
            'physical_prefactors_verified':False,
            'physical_coupling_map_certified':False,
            'models':blocks}
        audit=stage22_inventory(report)
        self.assertEqual(audit['A_-2']['extra_quartic_candidate_block_count'],1)
        self.assertFalse(audit['A_-2']['may_be_used_as_master_quartic_components'])
        report['physical_coupling_map_certified']=True
        with self.assertRaises(SourceKernelError):stage22_inventory(report)

    def test_source_contract_explicit_missing_inputs(self):
        c=required_couplings(source(),'abc')
        self.assertEqual(c['format'],CONTRACT_FORMAT)
        self.assertIsNone(c['physical_source_spinor_factor'])
        self.assertEqual(c['known_zero_at_benchmark'],[19])
        self.assertEqual(c['n_real_scalars'],12)
        self.assertFalse(c['physical_quartic_normalization_verified'])

    def test_quartic_duplicate_source_channels_add(self):
        b=source();q=diagnostic_couplings(b)
        vals, inventory=_merge_quartics(q['quartic_channels'],12)
        self.assertEqual(len(vals),1)
        self.assertAlmostEqual(vals[0]['value'],0.2)
        self.assertEqual([i['label'] for i in inventory],['portal1','portal2'])

    def test_missing_or_conflicting_quartic_fail_closed(self):
        b=source();q=diagnostic_couplings(b)
        del q['quartic_channels']
        with self.assertRaises(SourceKernelError): _validate_couplings(q,b,'abc')
        q=diagnostic_couplings(b)
        q['quartic_channels'][0]['components'].append({'indices':[3,3,7,5],'value':0.13})
        with self.assertRaisesRegex(SourceKernelError,'conflicting'):_validate_couplings(q,b,'abc')

    def test_changed_source_sha_or_norm_rejected(self):
        b=source();q=diagnostic_couplings(b)
        with self.assertRaises(SourceKernelError):evaluate(b,'different',q)
        b['spinor_normalization']['certified']=True
        with self.assertRaises(SourceKernelError):evaluate(b,'abc',q)

    def test_no_implicit_zero_gauge(self):
        b=source();q=diagnostic_couplings(b)
        q['coverage_assertion']['toy_gauge_couplings_set_to_zero']=False
        with self.assertRaisesRegex(SourceKernelError,'omitting gauge'):_validate_couplings(q,b,'abc')

    def test_wrong_model_or_scheme_or_muS_fails(self):
        b=source()
        for changed,val in [('model','B_-1'),('renormalization_scheme','on-shell'),('mu_S_gev',1e9)]:
            q=diagnostic_couplings(b);q[changed]=val
            with self.subTest(changed=changed):
                with self.assertRaises(SourceKernelError):evaluate(b,'abc',q)

    def test_source_recomposition_versus_master_engine(self):
        b=source();q=diagnostic_couplings(b)
        from RGE.running.intermediate.SourceBackedEFT1Batch import _validate_couplings
        payload,_=_validate_couplings(q,b,'abc')
        ctx=load_boundary(payload,allow_synthetic=True)
        audit=_operator_fidelity(ctx,b)
        self.assertEqual(audit['PL_sources_preserved'],2)
        self.assertEqual(audit['absolute_max_tensor_conversion_error'],0.)
        c=ctx.operators[0][1]
        self.assertEqual(c[0,2,4,6],1.8e-11+3.3e-12j)
        self.assertEqual(c[2,0,6,4],c[0,2,4,6])

    def test_all_seven_terms_and_weinberg_pivot_versus_direct(self):
        b=source();q=diagnostic_couplings(b)
        r=evaluate(b,'abc',q,trajectory_points=5)
        self.assertEqual(r['status'],'SourceBackedSevenTermFormalPerNGamma_NotPhysicalC5')
        self.assertFalse(r['physical_C5_computed'])
        self.assertIsNone(r['N_GammaCC_to_master_weyl'])
        self.assertEqual(len(r['trajectory_per_NGamma']),5)
        self.assertEqual(set(r['conditional_master_beta_by_rge_term']),set(TERMS))
        np.testing.assert_allclose(
           np.asarray(r['quartic_only_master_neutrino_pivot_C5_beta']['real'])[0,1],
           4*2*(1.8e-11)*.2, rtol=1e-13,atol=1e-25)
        self.assertFalse(r['quartic_only_SU2_projection']['SU2_Ward_tensor_algebra_pass'])
        # Incorrect full SU2 projection cannot be interpreted as physical C5.
        self.assertFalse(r['joint_scalar_matching_handoff']['ready_for_certified_joint_scalar_matching'])

    def test_source_7term_zero_yukawa_other_terms_except_quartic(self):
        b=source();q=diagnostic_couplings(b)
        q['yukawa_components']=[]
        r=evaluate(b,'abc',q)
        for term in TERMS:
            if term!='scalar_pair':
                c=r['conditional_master_beta_by_rge_term'][term]
                self.assertLess(np.max(abs(np.asarray(c['real'])+1j*np.asarray(c['imag']))),1e-25,term)


class BatchCLITests(unittest.TestCase):
    def test_one_model_contract_and_no_overwrite_and_all16_scan(self):
        with tempfile.TemporaryDirectory() as td:
            root=Path(td)
            release=root/'release'
            entries=[]
            for model in MODELS:
                b=source(model)
                path=release/tag(model)/'boundary.json'
                path.parent.mkdir(parents=True,exist_ok=True)
                path.write_text(json.dumps(b))
                entries.append({'model':model,'PL_sources':2,'path':str(Path(tag(model))/'boundary.json'),
                                'mu_F_gev':1e5})
            (release/'summary.json').write_text(json.dumps({
                'format':'t3-lo-fermion-boundary-source-resolved-provisional-v1',
                'n_models_verified':16,
                'absolute_spinor_normalization_certified':False,
                'physical_eft1_workers_must_reject_this_as_certified_input':True,
                'models':entries,
            }))
            out=root/'output/general_observable/lo_eft1_mixing/contracts'
            original_root=RunSourceBackedEFT1Batch.ROOT
            try:
                RunSourceBackedEFT1Batch.ROOT=root
                args=['--release-root',str(release),'--output-root',str(out)]
                self.assertEqual(RunSourceBackedEFT1Batch.main(args+['--check-only']),0)
                self.assertFalse(out.exists())
                self.assertEqual(RunSourceBackedEFT1Batch.main(args),0)
                summary=json.loads((out/'summary.json').read_text())
                self.assertEqual(summary['model_count'],16)
                self.assertFalse(summary['physical_C5_computed'])
                self.assertTrue((out/'C_m1/required_couplings_contract.json').exists())
                with self.assertRaises(FileExistsError):RunSourceBackedEFT1Batch.main(args)
                with self.assertRaises(SourceKernelError):
                    RunSourceBackedEFT1Batch.main(args+['--output-root',str(root/'illegal')])
            finally:
                RunSourceBackedEFT1Batch.ROOT=original_root


if __name__=='__main__':unittest.main()
