"""Stage 4 integration tests; synthetic states and conditional B/C reports only.

Requires native project modules. No missing matching coefficient is generated.
"""
from __future__ import annotations

from dataclasses import replace
import unittest

import numpy as np

from Numerical.core.State import SMNumericalState, T3Representation, T3UVState
from Numerical.running.FermionThresholdBoundary import build_intermediate_scalar_boundary
from Numerical.running.SMWeinbergEvolution import SMWeinbergInitialConditions
from Numerical.experimental.GeneralT3ScalarRGState import GeneralT3ScalarRGState
from Numerical.experimental.GeneralT3EFTHandoffController import (
    Callback, EFTHandoffController, ThresholdEvent, HandoffContractError,
)
from Numerical.experimental.GeneralT3NumericalStateAdapters import bridge_native_state
from Numerical.experimental.GeneralT3NumericalChecks import (
    EXPECTED_BC_MASSES_EV, bc_wilson_block, check_bc_conditional_report,
    check_neutrino_matrix, complex_matrix_payload, validate_bc_wilson_tensor,
    validate_c5,
)


def sm():
    diag = np.diag([1e-5,0.01,1.]).astype(complex)
    return SMNumericalState(gY=0.35,g2=0.65,g3=1.16,lambdaH=0.13,
                            yu=diag,yd=diag*.1,ye=diag*.01).validated()


def uv(alpha=-4):
    return T3UVState(mu_gev=1e7,
        representation=T3Representation(d_s1=1,d_s2=3,d_f=2,alpha=alpha),
        sm=sm(), y1=np.eye(3,dtype=complex)*0.02,
        y2=np.eye(3,dtype=complex)*0.03,
        MF=np.diag([1e5,1.1e5,1.2e5]).astype(complex),
        mS1Sq=1.2e6,mS2Sq=1.1e6,
        lambdaS1=0.03,lambdaS2=0.04,lambdaH1=0.01,lambdaH2=0.02,
        lambda12=0.05,lambdaT3=0.002+0.001j,
        lambdaH2Adj=0.006,lambdaS2Adj=0.007).validated()


def make_c5():
    return np.array([[1e-15, 0.3e-15j, 0],
                     [0.3e-15j, -2e-15, 0.2e-15],
                     [0,0.2e-15,3e-15]],dtype=complex)


def pair(z):
    a=np.asarray(z,dtype=complex)
    return {'real':a.real.tolist(),'imag':a.imag.tolist()}


def synthetic_report(model):
    # The synthetic test report reproduces the saved benchmark spectrum by
    # construction. It is NEVER used as a replacement for actual B/C reports.
    masses=np.asarray(EXPECTED_BC_MASSES_EV[model],dtype=float)
    matrix=np.diag(masses)
    c5=-2*matrix/(246.22**2*1e9)
    return {
        'format':'t3-bc-conditional-sequential-hx-to-sm-v1',
        'status':'EFT2HXRunningXCrossingSMWeinbergConditionalFiniteScalarMatchingMissing',
        'model':model,'mu_Y_gev':1190.,'mu_X_corrected_gev':1030.,
        'C5_at_X_tree_gev_inverse':pair(c5),
        'C5_at_100gev_gev_inverse':pair(c5),
        'mnu_at_100gev_ev':pair(matrix),
        'takagi_masses_ev_descending':masses.tolist(),
        'finite_Y_threshold_delta_C_computed':False,
        'finite_X_threshold_delta_C5_computed':False,
        'full_sequential_matching_performed':False,
        'physical_neutrino_mass_prediction_certified':False,
        'finite_threshold_corrections_status':'BOTH_MISSING_NOT_PHYSICAL_ZERO',
        'leading_log_sequential_candidate_completed':True,
    }


class TestStage4NumericalInterfaces(unittest.TestCase):
    def test_uv_roundtrip_preserves_native_complex_blocks(self):
        native=uv()
        bridge=bridge_native_state('A_-4',native,stage='UV_A',
              active_particles=('fermion:2','fermion:1','fermion:0','scalar:S1','scalar:S2'))
        self.assertIn('MF',bridge.state.parameters)
        self.assertEqual(bridge.state.parameters['MF'].shape,(3,3))
        self.assertEqual(bridge.contract.parameters[0].name,'gY')
        restored=bridge.restore(bridge.state).native
        np.testing.assert_array_equal(restored.MF,native.MF)
        np.testing.assert_array_equal(restored.sm.ye,native.sm.ye)
        self.assertEqual(restored.lambdaT3,native.lambdaT3)

    def test_native_model_mismatch_rejected(self):
        with self.assertRaises(HandoffContractError):
            bridge_native_state('D_+2',uv(),stage='UV_A',active_particles=('fermion:2',))

    def test_uv_mutation_roundtrips_at_new_scale_without_matching(self):
        b=bridge_native_state('A_-4',uv(),stage='UV_A',active_particles=('fermion:2',))
        x=dict(b.state.parameters)
        x['MF']=np.diag([111.,222.,333.]).astype(complex)
        restored=b.restore(replace(b.state,mu_gev=9e6,parameters=x))
        self.assertEqual(restored.native.mu_gev,9e6)
        np.testing.assert_array_equal(np.diag(restored.native.MF),[111,222,333])

    def test_explicit_mu2_sidecar_roundtrip(self):
        z=930.+17.j
        b=bridge_native_state('A_-2',uv(alpha=-2),stage='UV_A_m2',
                              active_particles=('fermion:2',),
                              quadratic_kind='self_s2',quadratic_mu2_gev2=z)
        self.assertEqual(b.state.parameters['mu2'],z)
        self.assertEqual(b.restore(b.state).quadratic_mu2_gev2,z)
        self.assertEqual(b.contract.parameters[-1].units,'GeV^2')

    def test_wrong_mu2_contract_rejected(self):
        with self.assertRaises(HandoffContractError):
            bridge_native_state('A_-2',uv(alpha=-2),stage='UV',active_particles=(),
                                quadratic_kind='self_s2')
        with self.assertRaises(HandoffContractError):
            bridge_native_state('A_-4',uv(),stage='UV',active_particles=(),
                                quadratic_mu2_gev2=1j)

    def test_eft1_state_roundtrips_original_source(self):
        native=build_intermediate_scalar_boundary(uv())
        b=bridge_native_state('A_-4',native,stage='EFT1_A',
                             active_particles=('scalar:S1','scalar:S2'))
        self.assertNotIn('MF',b.state.parameters)
        self.assertEqual(b.source_kind,'eft1')
        restored=b.restore(b.state).native
        np.testing.assert_array_equal(restored.sm.yu,native.sm.yu)
        self.assertEqual(restored.lambdaT3,native.lambdaT3)

    def test_general_scalar_registry_roundtrip(self):
        # Use the actual registry when present. It is source material, not a
        # substitute hand-written inventory. Missing registry skips this test.
        from pathlib import Path
        from Numerical.experimental.GeneralT3ScalarRGState import load_registry
        registry_path=Path(__file__).resolve().parents[1]/'RGE/general/generated/general_rgbeta_quartic_registry.json'
        if not registry_path.is_file():
            self.skipTest('Candidate RGBeta quartic registry not installed')
        registry=load_registry(registry_path)
        extras={entry['name']:0.002+0.001j
                for entry in registry['models']['A_+0']['new_rgbeta_symbols']}
        g=GeneralT3ScalarRGState(uv(alpha=0),extras).validated(registry)
        b=bridge_native_state('A_+0',g,stage='UV_A0',active_particles=('fermion:2',),
                              registry=registry,quadratic_kind='self_s1',
                              quadratic_mu2_gev2=300.+9j)
        restored=b.restore(b.state)
        self.assertEqual(restored.quadratic_mu2_gev2,300+9j)
        self.assertEqual(restored.native.extra,g.extra)

    def test_final_sm_weinberg_native_k_to_c5_roundtrip(self):
        orig=SMWeinbergInitialConditions(gY=.35,g2=.65,g3=1.16,lambdaH=.13,
               ye=sm().ye,yu=sm().yu,yd=sm().yd,K=make_c5()).validated()
        b=bridge_native_state('B_-1',orig,stage='SM_W',active_particles=(),mu_gev=100.,
                    c5_operator_basis='saved conditional project Weinberg',
                    c5_normalization='mnu=-(v^2/2)C5')
        np.testing.assert_array_equal(b.state.wilson['C5'],orig.K)
        np.testing.assert_array_equal(b.restore(b.state).native.K,orig.K)
        self.assertEqual(b.source_kind,'sm_weinberg')

    def test_final_sm_c5_missing_provenance_fails_closed(self):
        orig=SMWeinbergInitialConditions(gY=.35,g2=.65,g3=1.16,lambdaH=.13,
               ye=sm().ye,yu=sm().yu,yd=sm().yd,K=make_c5())
        with self.assertRaises(HandoffContractError):
            bridge_native_state('B_-1',orig,stage='SM',active_particles=(),mu_gev=100.)

    def test_complex_symmetric_c5_and_native_takagi(self):
        x=check_neutrino_matrix(make_c5())
        self.assertEqual(x.masses_ev_descending.shape,(3,))
        self.assertLess(x.takagi_residual,1e-7)
        self.assertLess(x.svd_max_relative_error,1e-7)
        self.assertFalse(np.allclose(x.neutrino_mass_gev,x.neutrino_mass_gev.conj().T,
                                     rtol=1e-12,atol=0))

    def test_c5_small_scale_antisymmetric_part_rejected(self):
        matrix=make_c5()
        matrix[1,0]+=1e-16
        with self.assertRaises(HandoffContractError):
            validate_c5(matrix)

    def test_wilson_tensor_shape_and_pair_exchange(self):
        z=np.zeros((3,3,8,8,6),dtype=complex)
        z[1,2,0,4,2]=3+2j
        z[2,1,0,4,2]=3+2j
        z[1,2,4,0,2]=3+2j
        z[2,1,4,0,2]=3+2j
        np.testing.assert_array_equal(validate_bc_wilson_tensor(z),z)
        spec=bc_wilson_block('C_BC',8,operator_basis='one-generation L/eC conditional',
                             normalization='source tensor entry conventions',mass_dimension=-1)
        self.assertEqual(spec.shape,z.shape)
        z[1,2,0,4,2]+=0.1
        with self.assertRaises(HandoffContractError):
            validate_bc_wilson_tensor(z)

    def test_wilson_tensor_wrong_shapes_rejected(self):
        for x in (np.zeros((3,3,8,8)),np.zeros((3,3,10,10,6)),
                  np.zeros((3,3,8,8,5))):
            with self.subTest(shape=x.shape),self.assertRaises(HandoffContractError):
                validate_bc_wilson_tensor(x)
        with self.assertRaises(HandoffContractError):
            bc_wilson_block('C',8,operator_basis='',normalization='',mass_dimension=-1)

    def test_synthetic_report_checks_both_benchmarks(self):
        for model in EXPECTED_BC_MASSES_EV:
            with self.subTest(model=model):
                out=check_bc_conditional_report(synthetic_report(model),expected_model=model)
                self.assertFalse(out['physical_neutrino_mass_prediction_certified'])

    def test_c5_report_tampering_or_false_certification_rejected(self):
        record=synthetic_report('B_-1')
        record['finite_X_threshold_delta_C5_computed']=True
        with self.assertRaises(HandoffContractError):
            check_bc_conditional_report(record)
        record=synthetic_report('B_-1')
        record['C5_at_100gev_gev_inverse']['real'][0][1]=1e-15
        with self.assertRaises(HandoffContractError):
            check_bc_conditional_report(record)

    def test_no_physical_matching_registered_handoff_stops(self):
        b=bridge_native_state('A_-4',uv(),stage='UV_only',
                              active_particles=('fermion:2','fermion:1','fermion:0',
                                                'scalar:S1','scalar:S2'))
        ctl=EFTHandoffController()
        ctl.add_stage(b.contract)
        evt=(ThresholdEvent('fermion:2','fermion',125000.,'next_unknown'),)
        outcome=ctl.run(b.state,evt,final_scale_gev=1000.)
        self.assertEqual(outcome.status,'Blocked')
        self.assertIn('Missing RG running',outcome.reason)
        self.assertEqual(outcome.state.mu_gev,1e7)
        self.assertFalse(outcome.physical_matching_certified)

    def test_report_real_imag_encoding_rejected_if_malformed(self):
        with self.assertRaises(HandoffContractError):
            complex_matrix_payload({'real':[[1]],'imag':[[2]]},shape=(3,3))


if __name__ == '__main__':
    unittest.main()
