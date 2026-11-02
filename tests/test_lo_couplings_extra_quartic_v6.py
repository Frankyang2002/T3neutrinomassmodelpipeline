"""Standalone fail-closed V6 source/operator and live Matchete CG tests."""
from __future__ import annotations
import copy
import json
import tempfile
from pathlib import Path
import unittest

from Numerical.lo_certification.couplings.extra_quartic_source_v6 import (
    ALL_MODELS,V6SourceError,extra_name,canonical_counts,model_registry_quartics,
    generator_contract,parse_live_cg,parse_live_operator_export,
    proportional_candidate,seed_belongs_to_model,indexed_seed_evidence,
)
from Numerical.lo_certification.couplings.extra_quartic_batch_v6 import (
    _validate_output,_validate_registry,_load_v5,run_v6_batch,
)


COUNTS=(0,2,0,1,0,1)


def sample_entry(counts=COUNTS,idx=1):
    dims=[2,2,3] if counts==COUNTS else []
    values=['0']*12 if dims else ['1']
    if dims:values[0]='Sqrt[2]/2';values[2]='-I/2'
    return {'canonical_field_counts':list(counts),'invariant_index':idx,
        'coupling_symbol':extra_name(counts,idx),
        'CG_symbol':extra_name(counts,idx,cg=True),
        'CG_dimensions':dims,'CG_components_inputform':values,
        'implicit_singlet_CG':not bool(dims),
        'present_in_validated_Lagrangian':True}


def expected_entry(counts=COUNTS):
    return {'rgbeta_symbol':'lambdaT3GeneralQ020101I1',
      'canonical_field_counts':list(counts),'sector_HHSS':True}


def live_doc(entries):
    return {'format':'t3-ws2-live-Matchete-extra-quartic-CGs-v6',
            'model':'A_-2','status':'LiveExtraQuarticCGExportVerified',
            'full_matchete_lagrangian_validation_passed':True,
            'operators':entries}


class TestV6SourceNames(unittest.TestCase):
    def test_model_inventory_16(self):self.assertEqual(len(ALL_MODELS),16)
    def test_bosonic_canonical_conjugation(self):
        self.assertEqual(canonical_counts((2,0,1,0,0,1)),(0,2,0,1,1,0))
    def test_generator_name(self):
        self.assertEqual(extra_name(COUNTS,1),'lambdaT3ExtraQ020101I1')
        self.assertEqual(extra_name(COUNTS,1,cg=True),'T3ExtraCGQ020101I1')
    def test_invalid_counts_rejected(self):
        with self.assertRaises(V6SourceError):canonical_counts((1,1,1,1,0,1))
    def test_explicit_model_path_required(self):
        self.assertTrue(seed_belongs_to_model(Path('output/full16_symbolic/T3_A_alpha_m2/data/eft1_after_F_scalar_quartic_seed.json'),'A_-2'))
        self.assertFalse(seed_belongs_to_model(Path('output/full16_symbolic/T3_A_alpha_m4/data/eft1_after_F_scalar_quartic_seed.json'),'A_-2'))
    def test_two_rgbeta_orientations_share_canonical_sector(self):
        x=(1,0,1,0,1,1);y=(0,1,0,1,1,1)
        self.assertEqual(canonical_counts(x),canonical_counts(y))
    def test_registry_grouped_indices_remain_unverified(self):
        row={'models':{'A_-2':{'total_additional_quartics':2,
             'new_rgbeta_symbols':[{'name':'first','counts':list(COUNTS)},
              {'name':'second','counts':list(COUNTS)}],
             'preexisting_rgbeta_symbols':[]}}}
        rows=model_registry_quartics('A_-2',row)
        self.assertEqual(len(rows),2)
        self.assertEqual(rows[1]['Matchete_generator_candidate_symbol'],'lambdaT3ExtraQ020101I2')
        self.assertTrue(all(x['Matchete_candidate_index_unverified'] for x in rows))


class TestV6LiveSource(unittest.TestCase):
    def test_exact_radical_and_phase_survive(self):
        cg=parse_live_cg(sample_entry(),model='A_-2')
        self.assertIsNotNone(cg)
        self.assertTrue(cg.invariant_summary()['contains_explicit_complex_phase'])
        self.assertEqual(cg.invariant_summary()['nonzero_components'],2)
    def test_exact_implicit_singlet_not_missing_tensor(self):
        entry=sample_entry((0,0,0,4,0,0))
        self.assertIsNone(parse_live_cg(entry,model='A_-2'))
    def test_fake_number_fails(self):
        entry=sample_entry();entry['CG_components_inputform'][0]='1.234'
        with self.assertRaises(Exception):parse_live_cg(entry,model='A_-2')
    def test_bad_dimension_fails(self):
        entry=sample_entry();entry['CG_dimensions']=[2,2,2]
        with self.assertRaises(V6SourceError):parse_live_cg(entry,model='A_-2')
    def test_bad_operator_name_fails(self):
        entry=sample_entry();entry['coupling_symbol']='fake'
        with self.assertRaises(V6SourceError):parse_live_cg(entry,model='A_-2')
    def test_valid_live_batch_one(self):
        r=parse_live_operator_export(live_doc([sample_entry()]),'A_-2',[expected_entry()])
        self.assertEqual(r['live_operator_count'],1)
        self.assertEqual(r['live_CG_count'],1)
    def test_missing_live_operators_fail(self):
        with self.assertRaises(V6SourceError):
            parse_live_operator_export(live_doc([]),'A_-2',[expected_entry()])
    def test_invalid_live_validation_fail(self):
        d=live_doc([sample_entry()]);d['full_matchete_lagrangian_validation_passed']=False
        with self.assertRaises(V6SourceError):parse_live_operator_export(d,'A_-2',[expected_entry()])
    def test_duplicate_live_invariant_fail(self):
        with self.assertRaises(V6SourceError):
            parse_live_operator_export(live_doc([sample_entry(),sample_entry()]),'A_-2',
              [expected_entry(),{'canonical_field_counts':list(COUNTS),'sector_HHSS':True}])
    def test_source_ratio_exact_sqrt(self):
        self.assertEqual(proportional_candidate(['sqrt(2)','-I'],['1','-sqrt(2)*I/2']),'sqrt(2)')
    def test_bad_proportionality_does_not_fit(self):
        self.assertIsNone(proportional_candidate(['1','2'],['1','3']))
    def test_unknown_source_reference_fails_closed(self):
        with self.assertRaises(V6SourceError):proportional_candidate(['os.system(1)'],['1'])


class TestV6BatchContracts(unittest.TestCase):
    def test_wrong_output_owner(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp)
            with self.assertRaises(ValueError):_validate_output(root,root/'RGE/audit_v6')
    def test_existing_output_never_overwritten(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);p=root/'output/general_observable/lo_certification_couplings/repeat'
            p.mkdir(parents=True)
            with self.assertRaises(FileExistsError):_validate_output(root,p)
    def test_v5_missing_fails_before_output_creation(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);out=root/'output/general_observable/lo_certification_couplings/v6'
            with self.assertRaises(Exception):run_v6_batch(root,out)
            self.assertFalse(out.exists())
    def test_source_inventory_requires_50_and_24(self):
        with self.assertRaises(V6SourceError):_validate_registry({'status':'GeneratedCandidateBasisNotWolframVerified','models':{}})
    def test_complete_synthetic_all16_source_only_run_no_physical_values(self):
        from Numerical.lo_certification.couplings.extra_quartic_source_v6 import ALL_MODELS
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp)
            owner=root/'output/general_observable/lo_certification_couplings'
            v5=owner/'ws2_matchete_rgbeta_dictionary_v1'
            v5.mkdir(parents=True)
            (v5/'summary.json').write_text(json.dumps({
                'models_processed':16,'all16_source_and_CG_crosswalk_passed':True}))
            interaction=root/'Lagrangian/interactions';interaction.mkdir(parents=True)
            (interaction/'GeneralScalarOperators.wl').write_text("""
                DefineCoupling[coupling, SelfConjugate -> False]
                contribution = coupling[] product If[
                "Expression" -> PlusHc[contribution]
                T3ExtendedSymmetrizeGroup[tensor_, positions_List]
                T3ExtendedVerifiedStructureCount[model]
                T3ExtraCG lambdaT3Extra
                """)
            (interaction/'GeneralScalarVerifiedQuartics.wl').write_text('verified table')
            (interaction/'T3Topology.wl').write_text('T3ExtendedBuildExtraQuarticCandidates[model]')
            builder=root/'Lagrangian/model/LagrangianBuilder.wl'
            builder.parent.mkdir(parents=True)
            builder.write_text('BuildT3InteractionCandidates[model]')
            counts={'A_-4':(0,0),'A_-2':(5,2),'A_+0':(5,2),
              'B_-3':(1,1),'B_-1':(7,4),'B_+1':(1,1),
              'C_-3':(1,1),'C_-1':(7,4),'C_+1':(1,1),
              'D_-2':(5,2),'D_+0':(5,2),'D_+2':(0,0),
              'E_-4':(0,0),'E_-2':(6,2),'E_+0':(6,2),'E_+2':(0,0)}
            reg={'status':'GeneratedCandidateBasisNotWolframVerified','models':{}}
            for model,(total,hhss) in counts.items():
                entries=[{'name':f'lambdaTest{n}',
                   'counts':list(COUNTS if n<hhss else (0,0,0,4,0,0))}
                  for n in range(total)]
                reg['models'][model]={'total_additional_quartics':total,
                  'new_rgbeta_symbols':entries,'preexisting_rgbeta_symbols':[]}
            registry=root/'RGE/general/generated/general_rgbeta_quartic_registry.json'
            registry.parent.mkdir(parents=True)
            registry.write_text(json.dumps(reg))
            out=owner/'v6_synthetic'
            result=run_v6_batch(root,out)
            self.assertTrue(result['all16_source_inventory_passed'])
            self.assertEqual(result['additional_HHSS_operator_inventory'],24)
            self.assertEqual(result['additional_quartic_operator_inventory'],50)
            self.assertEqual(result['live_native_extra_Matchete_operators_verified'],0)
            self.assertEqual(result['physically_certified_HHSS_coupling_maps'],0)
            self.assertFalse(result['physical_C5_calculated'])
            self.assertEqual(len(list(out.glob('*/V6_source_evidence.json'))),16)

    def test_seed_exact_matching_does_not_claim_physics(self):
        with tempfile.TemporaryDirectory() as tmp:
            path=Path(tmp)/'eft1_after_F_scalar_quartic_seed.json'
            path.write_text(json.dumps({'IntegratedFields':['F'],'ActiveHeavyFields':['S1','S2'],
                'CGRegistry':[], 'ScalarQuarticTermCount':1,
                'ScalarQuarticTerms':[{'Index':1,'CGNames':[],
                 'HiggsFieldCount':2,'Scalar1FieldCount':1,'Scalar2FieldCount':1,
                 'TermInputForm':(
                   'Coupling[lambdaT3ExtraQ020101I1, {}, 0]*'
                   'Bar[Field[H, Scalar, {}, {}]]*Bar[Field[H, Scalar, {}, {}]]*'
                   'Bar[Field[NewScalar1, Scalar, {}, {}]]*'
                   'Bar[Field[NewScalar2, Scalar, {}, {}]]')}] }))
            e=expected_entry()
            r=indexed_seed_evidence(path,'A_-2',[{'rgbeta_symbol':e['rgbeta_symbol'],
              'canonical_field_counts':list(COUNTS),'sector_HHSS':True,
              'Matchete_generator_candidate_symbol':'lambdaT3ExtraQ020101I1'}])
            self.assertEqual(r['named_extra_HHSS_matches'],1)
            self.assertFalse(r['entries'][0]['identifies_physical_RGBeta_map'])

if __name__=='__main__':unittest.main()

class TestV6ExperimentalUVExports(unittest.TestCase):
    def test_source_coupling_presence_is_not_an_EFT1_certificate(self):
        from Numerical.lo_certification.couplings.extra_quartic_source_v6 import experimental_inputform_evidence
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp)
            p=root/'output/experimental/T3_A_alpha_m2/export.json.interaction.inputform.txt'
            p.parent.mkdir(parents=True)
            p.write_text('Coupling[lambdaT3ExtraQ020101I1, {}, 0]')
            p.with_name('export.json').write_text(json.dumps({
                'full_matchete_check_passed':True,'quartic_expected':1,'quartic_included':1}))
            e=experimental_inputform_evidence(p,'A_-2',[
                {'Matchete_generator_candidate_symbol':'lambdaT3ExtraQ020101I1'}])
            self.assertTrue(e['all_expected_extra_names_present'])
            self.assertTrue(e['experimental_full_Lagrangian_validation_recorded'])
            self.assertTrue(e['cannot_certify_source_CG_or_EFT1_boundary_from_inputform_alone'])
    def test_no_cross_model_experimental_source_misassociation(self):
        from Numerical.lo_certification.couplings.extra_quartic_source_v6 import experimental_inputform_evidence
        with tempfile.TemporaryDirectory() as tmp:
            p=Path(tmp)/'output/T3_A_alpha_m4/test.json.interaction.inputform.txt'
            p.parent.mkdir(parents=True);p.write_text('anything')
            with self.assertRaises(V6SourceError):experimental_inputform_evidence(p,'A_-2',[])
