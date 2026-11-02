"""Tests for all-model ordered tensor handoff; synthetic Gamma is NEVER physics."""
from __future__ import annotations

import json
import math
from pathlib import Path
import tempfile
import unittest

import numpy as np

from Numerical.running.GeneralT3OrderedTensor import scalar_basis
from Numerical.running.GeneralT3WilsonEvolution import (
    FORMAT, STATUS, MODELS, GAMMA_FORMAT, tag, boundary_coordinates, digest,
    save_json, verify_batch, extend_to_gamma_basis, aggregate_source_coordinates, evolve_given_gamma,
)
from scripts.RunT3GeneralWilsonEvolution import main


def fake_tensor(model, source_index=1):
    basis=scalar_basis(model)
    arr=np.zeros((3,3));arr[0,0]=1
    component={'lepton_su2':[1,2], 'scalar_real':[5,5],
               'ordered_flavour_coefficient_gev_inv':{'real':arr.tolist(),
                                                       'imag':np.zeros((3,3)).tolist()}}
    operator={'source_index':source_index,'basis_scalar_fields':dict(basis.dimensions),
              'real_scalar_dimension_including_higgs':basis.total_dimension,
              'ordered_not_canonical_symmetrized':True,'component_count':1,
              'components':[component]}
    return {'model':model,
            'format':FORMAT if model not in ('B_-1','C_-1') else 't3-mixed-native-matchete-ordered-tensor-v1',
            'status':STATUS,'mu_F_gev':1e5,'scalar_basis_explicit':True,
            'matchete_cg_contractions_retained':True,
            'external_three_family_flavour_order_not_symmetrized':True,
            'canonical_C_ijab_normalization_certified':False,
            'tree_operator_count_PL':1,'total_ordered_real_component_count':1,
            'operators':[operator]}


def fake_gamma(model, boundary, sha, closure=True):
    labels,_,_ = boundary_coordinates(boundary,model)
    coords, values = aggregate_source_coordinates(labels, np.zeros(len(labels),dtype=complex))
    coords=coords+['closure:Z']
    n=len(coords)
    G=np.zeros((n,n));G[-1,0]=1
    return {'format':GAMMA_FORMAT,'model':model,'ordered_tensor_sha256':sha,
            'coordinates':coords,'closure_includes_zero_boundary_components':closure,
            'basis':'ordered_Matchete_real_components_with_flavour',
            'equation':'16pi2*dC/dlnmu=Gamma(mu)@C',
            'mu_F_gev':1e5,'mu_S_gev':1e3,
            'samples':[{'mu_gev':m,'gamma':{'real':G.tolist(),
                                          'imag':np.zeros((n,n)).tolist()}}
                       for m in (1e3,1e5)]}


class WilsonEvolutionTests(unittest.TestCase):
    def test_all16_basis_inventory(self):
        self.assertEqual(len(MODELS),16)
        for model in MODELS:
            labels,vector,report=boundary_coordinates(fake_tensor(model),model)
            self.assertEqual(len(labels),9)
            self.assertEqual(len(set(labels)),9)
            self.assertEqual(vector[0],1)
            self.assertEqual(report['mu_F_gev'],1e5)

    def test_no_automatic_symmetrization(self):
        source=fake_tensor('A_-4')
        source['operators'].append(dict(source['operators'][0],source_index=2))
        source['tree_operator_count_PL']=2
        source['total_ordered_real_component_count']=2
        labels,values,_=boundary_coordinates(source,'A_-4')
        self.assertEqual(len(labels),18)
        self.assertEqual(sum(labels[0] in {l} for l in labels),1)
        self.assertEqual(np.count_nonzero(values),2)
        combined,combined_values=aggregate_source_coordinates(labels,values)
        self.assertEqual(len(combined),9)
        self.assertEqual(combined_values[0],2)

    def test_wrong_scalar_basis_rejected(self):
        src=fake_tensor('A_-4')
        src['operators'][0]['real_scalar_dimension_including_higgs']=99
        with self.assertRaisesRegex(ValueError,'basis mismatch'):
            boundary_coordinates(src,'A_-4')

    def test_duplicate_operator_rejected(self):
        src=fake_tensor('A_-4')
        src['operators'].append(src['operators'][0]);src['tree_operator_count_PL']=2
        src['total_ordered_real_component_count']=2
        with self.assertRaisesRegex(ValueError,'duplicate source'):
            boundary_coordinates(src,'A_-4')

    def test_gamma_can_generate_initially_zero_closure(self):
        model='B_-1';raw=fake_tensor(model)
        out=evolve_given_gamma(raw,model,'f'*64,
                               fake_gamma(model,raw,'f'*64),npoints=5)
        self.assertEqual(out['status'],'ConditionalExplicitGammaNumericalTransport')
        self.assertEqual(out['closure_zero_initial_coordinates'],1)
        self.assertEqual(out['trajectory_real'][0][-1],0)
        predicted=math.log(1e3/1e5)/(16*math.pi**2)
        self.assertAlmostEqual(out['trajectory_real'][-1][-1],predicted,places=9)
        self.assertFalse(out['physical_eft1_wilson_rge_certified'])

    def test_incomplete_gamma_rejected(self):
        model='C_-1';raw=fake_tensor(model)
        bad=fake_gamma(model,raw,'a'*64)
        bad['coordinates']=bad['coordinates'][1:]
        with self.assertRaisesRegex(ValueError,'omits source'):
            evolve_given_gamma(raw,model,'a'*64,bad)

    def test_missing_closure_declaration_rejected(self):
        raw=fake_tensor('E_+2')
        g=fake_gamma('E_+2',raw,'c'*64,False)
        with self.assertRaisesRegex(ValueError,'closure inventory'):
            evolve_given_gamma(raw,'E_+2','c'*64,g)

    def test_stale_gamma_rejected(self):
        raw=fake_tensor('A_-4')
        with self.assertRaisesRegex(ValueError,'SHA-256'):
            evolve_given_gamma(raw,'A_-4','0'*64,fake_gamma('A_-4',raw,'1'*64))

    def test_cli_preflight_reads_all16_without_gamma(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp)/'input';out=Path(tmp)/'result'
            rows=[]
            for model in MODELS:
                path=root/tag(model)/'ordered_tensor.json'
                save_json(path,fake_tensor(model))
                rows.append({'model':model,'status':STATUS,
                             'tree_operator_count_PL':1,'real_component_count':1,
                             'tensor_sha256':digest(path)})
            save_json(root/'summary.json',{'format':FORMAT,'results':rows,
                                       'ordered_real_tensor_count':16,'blocked_count':0})
            self.assertEqual(len(verify_batch(root)),16)
            self.assertEqual(main(['--all','--mode','preflight',
                                  '--boundary-root',str(root),
                                  '--output-root',str(out),
                                  '--emit-coordinates']),0)
            summary=json.loads((out/'summary.json').read_text())
            self.assertEqual(summary['selected_model_count'],16)
            self.assertEqual(summary['conditional_gamma_evolved_count'],0)
            self.assertFalse(summary['physical_eft1_wilson_running_certified'])
            self.assertTrue((out/'B_m1'/'gamma_basis_requirements.json').is_file())
            self.assertTrue((out/'B_m1'/'summed_ordered_boundary.json').is_file())

    def test_cli_one_model_explicit_gamma_with_closure(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp)/'input';out=Path(tmp)/'evolved';gamma_root=Path(tmp)/'gamma'
            rows=[]
            for model in MODELS:
                path=root/tag(model)/'ordered_tensor.json'
                save_json(path,fake_tensor(model))
                rows.append({'model':model,'status':STATUS,
                             'tree_operator_count_PL':1,'real_component_count':1,
                             'tensor_sha256':digest(path)})
            save_json(root/'summary.json',{'format':FORMAT,'results':rows,
                                           'ordered_real_tensor_count':16,'blocked_count':0})
            model='A_-4'
            tensor=root/tag(model)/'ordered_tensor.json'
            save_json(gamma_root/tag(model)/'gamma.json',
                      fake_gamma(model,fake_tensor(model),digest(tensor)))
            self.assertEqual(main(['--model',model,'--mode','evolve',
                                   '--boundary-root',str(root),
                                   '--gamma-root',str(gamma_root),
                                   '--output-root',str(out),
                                   '--points','5']),0)
            result=json.loads((out/tag(model)/'conditional_gamma_running.json').read_text())
            self.assertEqual(len(result['scales_gev']),5)
            self.assertEqual(result['closure_zero_initial_coordinates'],1)
            self.assertFalse(result['physical_eft1_wilson_rge_certified'])

    def test_summary_digest_mismatch(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp)
            rows=[]
            for model in MODELS:
                save_json(root/tag(model)/'ordered_tensor.json',fake_tensor(model))
                rows.append({'model':model,'status':STATUS,'tree_operator_count_PL':1,
                             'real_component_count':1,'tensor_sha256':'a'*64})
            save_json(root/'summary.json',{'format':FORMAT,'results':rows,
                                       'ordered_real_tensor_count':16,'blocked_count':0})
            with self.assertRaisesRegex(ValueError,'digest changed'):
                verify_batch(root)

if __name__=='__main__':unittest.main()
