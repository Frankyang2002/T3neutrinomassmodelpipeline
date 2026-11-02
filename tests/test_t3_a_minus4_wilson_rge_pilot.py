"""Regression tests for the A_-4 RGE structural audit and conditional solver."""
from __future__ import annotations

import math
import unittest
import tempfile
import json
from pathlib import Path
from unittest.mock import patch

import numpy as np
from scipy.sparse import csr_matrix

from Numerical.running.GeneralT3RealLinearWilsonFlow import (
    LOOP, FORMAT, EQUATION, BASIS, encode_sparse_matrix,
    evolve_real_linear, parse_sparse_matrix,
)
from Numerical.running.GeneralT3WilsonRGEPilot import (
    read_flavour_components, symmetry_audit, master_antilinearity_check, audit_a_minus_4,
)


def source(A=None, B=None, *, labels=('L11:S5,5:F11','L11:S5,5:F12'), n=None):
    n=len(labels) if n is None else n
    if A is None:A=np.zeros((n,n),dtype=complex)
    if B is None:B=np.zeros((n,n),dtype=complex)
    return {'format':FORMAT, 'equation':EQUATION, 'basis':BASIS,
            'model':'A_-4','ordered_tensor_sha256':'a'*64,
            'closure_includes_zero_boundary_components':True,
            'conjugate_term_included':True,
            'mu_F_gev':1e5,'mu_S_gev':1e3,
            'coordinates':list(labels),
            'samples':[{'mu_gev':mu,'A':encode_sparse_matrix(A),'B':encode_sparse_matrix(B)}
                       for mu in (1e3,1e5)]}


def solve(boundary, gamma, *, labels=('L11:S5,5:F11','L11:S5,5:F12')):
    return evolve_real_linear(list(labels),np.asarray(boundary,dtype=complex),gamma,
                              model='A_-4', tensor_sha256='a'*64,mu_f=1e5,points=7)


class RealLinearSolverTests(unittest.TestCase):
    def test_sparse_exact_round_trip(self):
        A=np.array([[0,1+2j],[3-4j,0]],dtype=complex)
        self.assertTrue(np.array_equal(parse_sparse_matrix(encode_sparse_matrix(A),2).toarray(),A))

    def test_analytical_antilinear_solution(self):
        b=2.7
        gamma=source(B=b*np.eye(2,dtype=complex))
        output=solve([2+3j, 1-2j],gamma)
        dt=math.log(1e3/1e5)*LOOP
        expected=np.array([2*math.exp(b*dt)+3j*math.exp(-b*dt),
                           math.exp(b*dt)-2j*math.exp(-b*dt)])
        self.assertTrue(np.allclose(np.asarray(output['trajectory_real'][-1])+1j*np.asarray(output['trajectory_imag'][-1]),
                                    expected,atol=1e-10,rtol=1e-10))
        self.assertFalse(output['physical_eft1_wilson_rge_certified'])

    def test_real_linear_mixing_generates_zero_boundary_mode(self):
        gamma=source(A=np.array([[0,0],[0.5,0]],complex),
                     B=np.array([[0,0],[1.0,0]],complex))
        output=solve([2+1j],gamma,labels=('L11:S5,5:F11',))
        self.assertEqual(output['zero_initial_closure_coordinates'],1)
        z=complex(output['trajectory_real'][-1][1],output['trajectory_imag'][-1][1])
        self.assertGreater(abs(z),1e-7)

    def test_zero_matrices_give_constant_boundary(self):
        out=solve([1j,3+2j],source())
        self.assertLess(out['relative_change'],1e-12)

    def test_omitted_B_matrix_rejected(self):
        g=source()
        del g['samples'][0]['B']
        with self.assertRaisesRegex(ValueError,'A, B'):
            solve([1,0],g)

    def test_duplicate_sparse_coordinate_rejected(self):
        bad={'format':'complex-coo-v1','shape':[2,2],
             'row':[0,0], 'col':[1,1], 'real':[1,2], 'imag':[0,0]}
        with self.assertRaisesRegex(ValueError,'Duplicate'):
            parse_sparse_matrix(bad,2)

    def test_wrong_provenance_rejected(self):
        g=source();g['ordered_tensor_sha256']='b'*64
        with self.assertRaisesRegex(ValueError,'SHA-256'):
            solve([1,0],g)

    def test_missing_boundary_coordinate_rejected(self):
        g=source(labels=('L11:S5,5:F11',))
        with self.assertRaisesRegex(ValueError,'omits'):
            solve([1,0],g)

    def test_explicit_closure_declaration_required(self):
        g=source();g['closure_includes_zero_boundary_components']=False
        with self.assertRaisesRegex(ValueError,'complete closure'):
            solve([1,0],g)


class BasisAuditTests(unittest.TestCase):
    def test_symmetric_tensor_has_zero_residual(self):
        # Scalar and fermion swaps generate exactly the other three components.
        seed=(1,1,2,2,5,7)
        inverse=(2,2,1,1,5,7)
        ss=(1,1,2,2,7,5)
        both=(2,2,1,1,7,5)
        info=symmetry_audit({seed:1+2j,inverse:1+2j,ss:1+2j,both:1+2j})
        self.assertLess(info['removed_by_separate_symmetry_projection_relative'],1e-14)

    def test_antisymmetric_scalar_tensor_detected(self):
        info=symmetry_audit({(1,1,1,1,5,7):1,(1,1,1,1,7,5):-1})
        self.assertGreater(info['scalar_exchange_residual_relative'],1)
        self.assertGreater(info['removed_by_separate_symmetry_projection_relative'],0.99)

    def test_real_basis_index_validation(self):
        vals=np.asarray([1+0j])
        components=read_flavour_components(['L12:S5,7:F13'], vals,'A_-4')
        self.assertEqual(components[(1,1,2,3,5,7)],1+0j)
        with self.assertRaisesRegex(ValueError,'exceeds'):
            read_flavour_components(['L12:S5,99:F13'],vals,'A_-4')

    def test_a_minus4_audit_fixture(self):
        from Numerical.running.GeneralT3OrderedTensor import scalar_basis
        sb=scalar_basis('A_-4')
        arr=np.zeros((3,3));arr[0,0]=1
        payload={'model':'A_-4',
            'format':'t3-ordered-source-flavour-and-real-scalar-tensor-v1',
            'status':'OrderedNumericalRealComponentTensorFromMatchete',
            'mu_F_gev':1e5,'scalar_basis_explicit':True,
            'matchete_cg_contractions_retained':True,
            'external_three_family_flavour_order_not_symmetrized':True,
            'canonical_C_ijab_normalization_certified':False,
            'tree_operator_count_PL':1,
            'total_ordered_real_component_count':1,
            'operators':[{'source_index':1,
                'basis_scalar_fields':dict(sb.dimensions),
                'real_scalar_dimension_including_higgs':sb.total_dimension,
                'ordered_not_canonical_symmetrized':True,
                'component_count':1,
                'components':[{'lepton_su2':[1,2], 'scalar_real':[5,6],
                    'ordered_flavour_coefficient_gev_inv':{
                        'real':arr.tolist(), 'imag':np.zeros((3,3)).tolist()}}]}]}
        result=audit_a_minus_4(payload,include_master_check=False)
        self.assertEqual(result['status'],'BoundaryAuditedPhysicalGammaNotConstructed')
        self.assertEqual(result['summed_ordered_flavour_coordinate_count'],9)
        self.assertEqual(result['hypothetical_symmetric_closure_sizes']['L_plus_eC_fermion_pair_by_scalar_pair'],3510)
        self.assertGreater(result['source_symmetry_diagnostic']['fermion_exchange_residual_relative'],0)
        self.assertFalse(result['physical_eft1_wilson_rge_certified'])

    def test_cli_audit_readonly_batch_fixture(self):
        from Numerical.running.GeneralT3WilsonEvolution import (
            MODELS, FORMAT as TFORMAT, STATUS as TSTATUS,
            save_json, digest, tag,
        )
        from Numerical.running.GeneralT3OrderedTensor import scalar_basis
        from scripts.RunT3Aminus4WilsonRGEPilot import main
        with tempfile.TemporaryDirectory() as td:
            root=Path(td)/'input'
            output=Path(td)/'audit'
            rows=[]
            for model in MODELS:
                sb=scalar_basis(model)
                arr=np.zeros((3,3));arr[0,0]=1.0
                raw={'model':model, 'format':TFORMAT if model not in ('B_-1','C_-1') else 't3-mixed-native-matchete-ordered-tensor-v1',
                     'status':TSTATUS, 'mu_F_gev':1e5,
                     'scalar_basis_explicit':True, 'matchete_cg_contractions_retained':True,
                     'external_three_family_flavour_order_not_symmetrized':True,
                     'canonical_C_ijab_normalization_certified':False,
                     'tree_operator_count_PL':1,'total_ordered_real_component_count':1,
                     'operators':[{'source_index':1, 'basis_scalar_fields':dict(sb.dimensions),
                                  'real_scalar_dimension_including_higgs':sb.total_dimension,
                                  'ordered_not_canonical_symmetrized':True,
                                  'component_count':1,'components':[
                                    {'lepton_su2':[1,2],'scalar_real':[5,5],
                                     'ordered_flavour_coefficient_gev_inv':{
                                         'real':arr.tolist(), 'imag':np.zeros((3,3)).tolist()}}]}]}
                file=root/tag(model)/'ordered_tensor.json'
                save_json(file,raw)
                rows.append({'model':model,'status':TSTATUS,'tree_operator_count_PL':1,
                             'real_component_count':1,'tensor_sha256':digest(file)})
            save_json(root/'summary.json',{'format':TFORMAT,'results':rows,
                      'ordered_real_tensor_count':16,'blocked_count':0})
            with patch('Numerical.running.GeneralT3WilsonRGEPilot.master_antilinearity_check',return_value={'complex_linearity_violated':True}):
                code=main(['--mode','audit','--boundary-root',str(root),
                           '--output-root',str(output)])
            self.assertEqual(code,0)
            report=json.loads((output/'A_m4'/'basis_and_master_rge_audit.json').read_text())
            self.assertFalse(report['physical_eft1_wilson_rge_certified'])
            self.assertTrue(report['actual_master_rge_antilinearity']['complex_linearity_violated'])
            self.assertTrue((root/'A_m4'/'ordered_tensor.json').exists())

    def test_actual_master_antilinarity(self):
        try:
            from RGE.general.WilsonTensorRGE import conjugate_coefficient_yukawa_term
        except ModuleNotFoundError:
            self.skipTest('Complete master RGE package unavailable in standalone test environment')
        evidence=master_antilinearity_check()
        self.assertTrue(evidence['complex_linearity_violated'])
        self.assertFalse(evidence['complete_A_B_constructed'])


if __name__=='__main__':unittest.main()
