"""v10 negative controls and optional live all-16 native evidence integration."""
from __future__ import annotations
from copy import deepcopy
import json
from pathlib import Path
import tempfile
import unittest

import numpy as np

from Numerical.parallel_t3.wilson.native_quartic_projectors_v10 import (
    FORMAT,NATIVE_STATUS,SPECIAL,V10Blocked,EXPECTED_MODELS,DIMENSIONS,
    expected_symbols,validate_native_evidence,compare_with_v9,audit,
    audit_from_paths,
)
from Numerical.parallel_t3.wilson.portal_real_action_v9 import STATUS as V9_STATUS


def fixture_native():
    rows={}
    for tag in sorted(EXPECTED_MODELS):
        d1,d2,df=DIMENSIONS[tag[0]]
        symbols=expected_symbols(tag)
        inventory=[{
            'symbol':name,'native_self_conjugate':real,
            'has_installed_RGBeta_projector':True,
            'native_unique_arrangement_count':24,
            'native_coupling_inputform':f'6*Matrix[{name}][##1] &',
            'native_invariant_inputform':'del[SU2L[fund],i,j]',
            'fields_inputform':'{Bar[H], H, Bar[S1], S1}',
        } for name,real in symbols.items()]
        diagonal=[{'symbol':name,'native_response':{'ok':True,'value':1.}}
                  for name,real in symbols.items() if real]
        blocks=[]
        for which,dim in ((1,d1),(2,d2)):
            ss=[f'lambdaH{which}']+([f'lambdaH{which}Adj'] if dim>1 else [])
            n=len(ss)
            arr=np.eye(n)
            if n==2:arr[0,1]=.12  # projectors need not be orthogonal
            blocks.append({'which':which,'dimension':dim,'symbols':ss,
               'native_response_matrix':[[{'ok':True,'value':float(arr[i,j])}
                  for j in range(n)] for i in range(n)]})
        rows[tag]={'ok':True,'dimensions':[d1,d2,df],
                   'alpha':int(tag.split('_')[1]),
                   'registered_quartics_count':len(symbols),
                   'registry_exact_source_inventory':True,
                   'inventory':inventory,
                   'native_real_self_projectors':diagonal,
                   'native_portal_response_blocks':blocks,
                   'native_real_self_projectors_all_unit':True,
                   'native_portal_response_blocks_evaluated':True,
                   'complex_quartics_are_inventory_only':True,
                   'native_numeric_S2_components_measured':False}
    return {'format':FORMAT,'status':NATIVE_STATUS,'registered_model_count':16,
      'native_Lam_evaluated':True,'native_quartic_projectors_evaluated':True,
      'native_coupling_response_is_symbolic_not_component_Lam':True,
      'native_numeric_S2_frame_measured':False,
      'native_Matchete_scalar_CG_evaluated':False,
      'full_scalar_potential_dictionary_certified':False,
      'physical_Wilson_evolution_authorized':False,
      'source_sha256':{s:'a'*64 for s in (
        'RGE/running/rgbeta/T3RGBetaModel.wl','scripts/RGBetaNativeAccess.wl',
        'scripts/RunT3ParallelANativeQuarticProjectorsV10.wl')},
      'models':rows}


def fixture_v9():
    rows={}
    for tag in sorted(EXPECTED_MODELS):
        d1,d2,_=DIMENSIONS[tag[0]]
        sectors=[]
        for which,dim in ((1,d1),(2,d2)):
            ss=[f'lambdaH{which}']+([f'lambdaH{which}Adj'] if dim>1 else [])
            n=len(ss)
            sectors.append({'which':which,'rgbeta_registered_couplings':ss,
                'native_matchete_couplings':[f'lambdaH{which}Inv{i+1}' for i in range(n)],
                'conditional_rgbeta_from_matchete_real_matrix':np.eye(n).tolist(),
                'max_relative_residual':1.e-15})
        rows[tag]=sectors
    return {'status':V9_STATUS,'models':rows,
      'native_Matchete_portal_coupling_slots':58,
      'triplet_sectors_spectral_normalization_verified':14,
      'anti_Hermitian_Matchete_phase_slots_checked':14,
      'certification':{'native_RGBeta_Lam_projector_numerically_executed':False}}


class TestNativeEvidenceSchema(unittest.TestCase):
    def test_inventory_178_150_28_and_16_classes(self):
        counts=[expected_symbols(t) for t in EXPECTED_MODELS]
        self.assertEqual(sum(map(len,counts)),178)
        self.assertEqual(sum(sum(v.values()) for v in counts),150)
        self.assertEqual(sum(sum(not x for x in v.values()) for v in counts),28)
        self.assertEqual(set(k for k,v in expected_symbols('B_-1').items() if not v),
                         set(SPECIAL)|{'lambdaT3'})
        self.assertEqual(set(k for k,v in expected_symbols('E_-2').items() if not v),
                         {'lambdaT3'})

    def test_native_symbolic_projector_coverage_and_v9_bridge(self):
        data=audit(fixture_native(),fixture_v9())
        self.assertEqual(data['registered_quartic_slots'],178)
        self.assertEqual(data['native_real_projector_diagonals_verified'],150)
        self.assertEqual(data['native_portal_response_matrix_elements'],110)
        self.assertEqual(data['native_portal_offdiagonal_responses'],52)
        self.assertEqual(data['complex_quartic_slots_inventory_only'],28)
        self.assertFalse(data['certification']['physical_Wilson_evolution_authorized'])
        row=data['models']['B_-1'][0]
        self.assertAlmostEqual(row['native_projector_observations_from_Matchete'][0][1],.12)

    def test_singular_native_projectors_rejected(self):
        native=fixture_native()
        native['models']['B_-1']['native_portal_response_blocks'][0][
           'native_response_matrix'][1][0]['value']=1.
        native['models']['B_-1']['native_portal_response_blocks'][0][
           'native_response_matrix'][0][1]['value']=1.
        with self.assertRaisesRegex(V10Blocked,'singular'):
            audit(native,fixture_v9())

    def test_wrong_native_6matrix_factor_rejected(self):
        native=fixture_native()
        native['models']['A_-4']['inventory'][0]['native_coupling_inputform'] = (
            '12*Matrix[lambdaH][##1] &')
        with self.assertRaisesRegex(V10Blocked,r'6\*Matrix'):
            validate_native_evidence(native)

    def test_wrong_native_scalar_diagonal_rejected(self):
        native=fixture_native()
        native['models']['A_-4']['native_real_self_projectors'][0][
            'native_response']['value']=2.
        with self.assertRaisesRegex(V10Blocked,'diagonal projector response'):
            validate_native_evidence(native)

    def test_native_coverage_fails_if_complex_channel_falsely_claimed(self):
        native=fixture_native()
        native['models']['C_-1']['inventory'][0]['native_self_conjugate']=False
        with self.assertRaisesRegex(V10Blocked,'self-conjugacy'):
            validate_native_evidence(native)

    def test_native_source_provenance_hash_guard(self):
        native=fixture_native()
        with tempfile.TemporaryDirectory() as td:
            root=Path(td)
            for path in native['source_sha256']:
                p=root/path;p.parent.mkdir(parents=True,exist_ok=True);p.write_text('real source')
            native['installed_RGBeta_main_file']='missing_RGBeta_file.m'
            with self.assertRaisesRegex(V10Blocked,'native source has changed'):
                validate_native_evidence(native,root)

    def test_no_unsupported_physical_certification(self):
        native=fixture_native()
        native['physical_Wilson_evolution_authorized']=True
        with self.assertRaisesRegex(V10Blocked,'unsupported native certification'):
            validate_native_evidence(native)

    def test_v9_tampered_transfer_matrix_fails(self):
        native=fixture_native();v9=fixture_v9()
        v9['models']['E_-2'][0]['conditional_rgbeta_from_matchete_real_matrix']=[
            [1.,0.],[0.,0.]]
        with self.assertRaisesRegex(V10Blocked,'singular'):
            audit(native,v9)

class TestCurrentRepositoryIntegration(unittest.TestCase):
    def test_native_full_178_150_110_from_local_wolfram_witness(self):
        root=Path(__file__).resolve().parents[1]
        home=root/'output/parallel_t3/wilson'
        n=home/'rgbeta_quartic_projectors_native_v10.json'
        v9=home/'portal_real_action_conditional_v9.json'
        if not n.is_file() or not v9.is_file():
            self.skipTest('Native v10 export not yet executed on installed Wolfram kernel')
        result=audit_from_paths(root,n,v9)
        self.assertEqual(result['native_real_projector_diagonals_verified'],150)
        self.assertEqual(result['native_portal_response_matrix_elements'],110)

if __name__=='__main__':unittest.main()
