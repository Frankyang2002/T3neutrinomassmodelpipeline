"""Synthetic exact negative controls; native Matchete is not mocked as executed."""
from __future__ import annotations
import json
from pathlib import Path
import tempfile
import unittest

import numpy as np
import sympy as sp

from Numerical.lo_certification.direct_uv.native_portal_dictionary import (
    PortalDictionaryError, FIVE, NATIVE_FORMAT, SOURCE_REL, FORMAT,
    _as_real, hermitian_portal_tensor, identity_tensor, parse_exact,
    portal_tensor, solve_singlet_basis, source_contract, validate_native_models,
    frozen_inputs, audit, hashfile,
)
from Numerical.lo_certification.direct_uv.native_portal_dictionary import all_models, CLASSES


def adjoint(d: int) -> sp.Matrix:
    # An exactly Hermitian, independent 4-index tensor built from a traceless
    # diagonal operator in the Higgs and scalar spaces.
    q=list(range(d))
    q=[sp.Integer(i)*2-(d-1) for i in q]
    return sp.Matrix(4,d*d,[int(i==j and a==b)*(sp.Integer(1) if i==0 else -1)*q[a]
        for i in range(2) for j in range(2) for a in range(d) for b in range(d)])


def fake_tensor(t:sp.Matrix,which:int,k:int,d:int)->dict:
    return {"name":f"HiggsPortal{which}Inv{k}","coupling":f"lambdaH{which}Inv{k}",
        "scalar":which,"invariant":k,"cg_symbol":f"T3HiggsPortal{which}CGInv{k}",
        "dimensions":[2,2]+([d,d] if d>1 else []),
        "components_inputform":[str(x) for x in list(t)],
        "tensor_accessor":"Matchete`GetCGTensor", "registered_in_assembled_interaction":True}


def fake_native():
    rows=[]
    for cls,(d1,d2,df,alpha) in FIVE.items():
        ports=[]
        for which,d in ((1,d1),(2,d2)):
            I=identity_tensor(d)
            if d==1:
                tensors=[3*I]
            else:
                X=adjoint(d)
                tensors=[2*I+3*X,4*X]
            ports += [fake_tensor(t,which,k,d) for k,t in enumerate(tensors,1)]
        rows.append({"model":f"{cls}_{alpha:+d}","class":cls,"alpha":alpha,
                     "dS1":d1,"dS2":d2,"dF":df,
                     "native_build_status":"Success", "native_full_lagrangian_validation":True,
                     "portals":ports})
    return {"format":NATIVE_FORMAT,
        "status":"five_classes_native_registered_portal_tensors_extracted",
        "restricted_environment":{"isolated_to_wolframscript_process":True,
            "effective":{x:"0" for x in (
                "T3_EXPERIMENTAL_GENERAL_FERMIONS",
                "T3_EXPERIMENTAL_GENERAL_QUADRATICS",
                "T3_EXPERIMENTAL_GENERAL_QUARTICS",
                "T3_EXPERIMENTAL_FLAVOR_MASSES")}},
        "native_rgbeta_portal_map_certified":False,
        "matching_or_running_performed":False,
        "models":rows}


class PortalDictionaryTests(unittest.TestCase):
    def test_mathematica_exact_parsing(self):
        self.assertEqual(parse_exact('Sqrt[2]/Sqrt[8]'),sp.Rational(1,2))
        self.assertEqual(parse_exact('Rational[7,8]'),sp.Rational(7,8))

    def test_expression_injection_rejected(self):
        with self.assertRaises(PortalDictionaryError):
            parse_exact('Sin[x]')

    def test_inexact_machine_precision_native_values_rejected(self):
        with self.assertRaises(PortalDictionaryError):
            parse_exact('0.3333333333333333')

    def test_singlet_source_normalization(self):
        r=solve_singlet_basis([3*identity_tensor(1)],1)
        self.assertEqual(r['coefficients_native_lambdaHiInvN_per_RGBeta_lambdaHi'],['1/3'])

    def test_two_invariant_basis_with_adjoint(self):
        I=identity_tensor(2);X=adjoint(2)
        r=solve_singlet_basis([2*I+3*X,4*X],2)
        self.assertEqual(r['coefficients_native_lambdaHiInvN_per_RGBeta_lambdaHi'],['1/2','-3/8'])

    def test_triplet_basis(self):
        I=identity_tensor(3);X=adjoint(3)
        r=solve_singlet_basis([2*I+3*X,4*X],3)
        self.assertEqual(r['coefficients_native_lambdaHiInvN_per_RGBeta_lambdaHi'],['1/2','-3/8'])

    def test_hermitian_part_retained_and_antihermitian_discarded(self):
        I=identity_tensor(2)
        raw=2*I+sp.I*adjoint(2)
        self.assertEqual(hermitian_portal_tensor(raw,2),2*I)
        self.assertEqual(solve_singlet_basis([raw],2)['coefficient_approx'],[0.5])

    def test_missing_singlet_span_rejected(self):
        with self.assertRaises(PortalDictionaryError):
            solve_singlet_basis([adjoint(3)],3)

    def test_duplicate_basis_rejected(self):
        with self.assertRaises(PortalDictionaryError):
            solve_singlet_basis([identity_tensor(3),2*identity_tensor(3)],3)

    def test_five_model_class_transfer(self):
        rows=validate_native_models(fake_native())
        self.assertEqual(set(rows),set(FIVE))
        self.assertEqual(sum(len(v) for v in rows.values()),10)

    def test_dependent_after_hermitian_projection_has_exact_nullspace(self):
        I=identity_tensor(3); X=adjoint(3)
        r=solve_singlet_basis([2*I,sp.I*X],3,label='E_+0/S1')
        self.assertEqual(r['raw_registered_CG_rank_over_complex'],2)
        self.assertEqual(r['real_Hermitian_operator_rank'],1)
        self.assertFalse(r['unique_coupling_dictionary'])
        self.assertEqual(r['coefficients_native_lambdaHiInvN_per_RGBeta_lambdaHi'],
                         ['1/2','0'])
        self.assertEqual(r['identically_zero_Hermitian_projection_indices_one_based'],[2])
        self.assertEqual(r['nullspace_basis_real_coupling_shifts'],[['0','1']])

    def test_proportional_hermitian_operators_are_nonunique(self):
        I=identity_tensor(2);X=adjoint(2)
        with self.assertRaisesRegex(PortalDictionaryError,
                                    'OUTSIDE native real-Hermitian'):
            solve_singlet_basis([2*I+3*X,4*I+6*X+sp.I*X],2)

    def test_degenerate_operator_span_with_singlet_present(self):
        I=identity_tensor(2);X=adjoint(2)
        r=solve_singlet_basis([2*I,4*I+sp.I*X],2)
        self.assertFalse(r['unique_coupling_dictionary'])
        self.assertEqual(r['real_Hermitian_operator_rank'],1)
        self.assertEqual(r['coefficients_native_lambdaHiInvN_per_RGBeta_lambdaHi'],
                         ['1/2','0'])
        self.assertEqual(r['nullspace_basis_real_coupling_shifts'],[['-2','1']])

    def test_no_singlet_in_dependent_span_blocks_with_diagnostics(self):
        I=identity_tensor(2);X=adjoint(2)
        with self.assertRaisesRegex(PortalDictionaryError,
                                    'OUTSIDE native real-Hermitian'):
            solve_singlet_basis([2*I+3*X,4*I+6*X+sp.I*X],2,label='B_-1/S1')

    def test_raw_dependent_registered_cgs_block(self):
        I=identity_tensor(2)
        with self.assertRaisesRegex(PortalDictionaryError,
                                    'RAW Matchete registered CG'):
            solve_singlet_basis([I,2*I],2)

    def test_nonunique_native_class_transfer(self):
        native=fake_native()
        r=native['models'][0]['portals'][1]
        r['components_inputform']=[str(z) for z in list(
            2*identity_tensor(3))]
        r2=native['models'][0]['portals'][2]
        r2['components_inputform']=[str(z) for z in list(
            sp.I*adjoint(3))]
        result=validate_native_models(native)
        self.assertFalse(result['A']['S2']['unique_coupling_dictionary'])
        self.assertTrue(result['B']['S2']['unique_coupling_dictionary'])

    def test_bad_source_portal_registration_rejected(self):
        native=fake_native()
        native['models'][0]['portals'][0]['registered_in_assembled_interaction']=False
        with self.assertRaises(PortalDictionaryError):
            validate_native_models(native)

    def test_nonrestricted_native_report_rejected(self):
        native=fake_native()
        native['restricted_environment']['effective']['T3_EXPERIMENTAL_GENERAL_FERMIONS']='1'
        with self.assertRaises(PortalDictionaryError):
            validate_native_models(native)

    def test_source_guards_and_native_hashes(self):
        with tempfile.TemporaryDirectory() as temp:
            root=Path(temp)
            src={
                SOURCE_REL[-1]: '1/2 PlusHc[ coupling[] fieldProduct If[cg === None, 1, CG[cg, labels]] ]',
                'RGE/running/rgbeta/T3RGBetaModel.wl': '''
T3RGBetaPairInvariant[1, i_, j_] := 1;
T3RGBetaPairInvariant[2, i_, j_] := del[SU2L @ fund, i, j];
T3RGBetaPairInvariant[3, i_, j_] := del[SU2L @ S2, i, j];
AddQuartic[lambdaH1,{Bar @ H,H,Bar @ S1,S1},GroupInvariant ->
(del[SU2L @ fund,#1,#2]*T3RGBetaPairInvariant[dS1,#3,#4] &)];
AddQuartic[lambdaH2,{Bar @ H,H,Bar @ S2,S2},GroupInvariant ->
(del[SU2L @ fund,#1,#2]*T3RGBetaPairInvariant[dS2,#3,#4] &)];
''',
                'scripts/BuildT3DirectOneLoopInputs.py': "-float(o['lambdaH1'])*v*v/2\n-float(o['lambdaH2'])*v*v/2",
            }
            for p in SOURCE_REL:
                src.setdefault(p,'native source stub')
            for path,content in src.items():
                p=root/path;p.parent.mkdir(parents=True,exist_ok=True)
                p.write_text(content,encoding='utf-8')
            native=fake_native()
            native['native_source_sha256']={p:hashfile(root/p) for p in SOURCE_REL}
            self.assertEqual(source_contract(root,native)['status'],
                             'native_source_bytes_match_and_portal_forms_guarded')
            (root/SOURCE_REL[0]).write_text('mutated',encoding='utf-8')
            with self.assertRaisesRegex(PortalDictionaryError,'stale'):
                source_contract(root,native)

    def test_all_16_synthetic_frozen_inputs(self):
        with tempfile.TemporaryDirectory() as temp:
            root=Path(temp)
            for model,key in all_models():
                ds1,ds2,df,_=CLASSES[model[0]]
                p=root/'Numerical/benchmarks/fitted_reference'/f'{key}.json'
                p.parent.mkdir(parents=True,exist_ok=True)
                p.write_text(json.dumps({'representation':{
                    'd_s1':ds1,'d_s2':ds2,'d_f':df,'alpha':int(model[2:])},
                    'base_state':{'ordinary':{'mS1Sq':1e6,'mS2Sq':2e6,
                    'lambdaH1':0.1,'lambdaH2':-0.2}},
                    'numerical':{'vev_gev':246.}}),encoding='utf-8')
                A=np.diag([1e6-0.1*246**2/2]*ds1+[2e6+0.2*246**2/2]*ds2)
                inp=root/'output/general_observable/lo_direct_one_loop/provisional_inputs'/f'{model}.json'
                inp.parent.mkdir(parents=True,exist_ok=True)
                inp.write_text(json.dumps({'model':model,'scalar_A_gev2':{'real':A.tolist(),'imag':np.zeros_like(A).tolist()}}),encoding='utf-8')
            report=frozen_inputs(root,validate_native_models(fake_native()))
            self.assertEqual(report['models_checked'],16)
            self.assertLess(report['max_relative_block_residual'],1e-15)
            first=root/'output/general_observable/lo_direct_one_loop/provisional_inputs/A_-4.json'
            obj=json.loads(first.read_text());obj['scalar_A_gev2']['real'][0][0]+=100
            first.write_text(json.dumps(obj),encoding='utf-8')
            with self.assertRaisesRegex(PortalDictionaryError,'frozen diagonal'):
                frozen_inputs(root,validate_native_models(fake_native()))


if __name__=='__main__':
    unittest.main()
