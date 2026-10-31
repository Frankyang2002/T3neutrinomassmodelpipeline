"""Stage 3 self-contained controller tests. All callbacks are SYNTHETIC, not physics."""
from __future__ import annotations

from dataclasses import replace
import copy
import json
import unittest

import numpy as np

from Numerical.experimental.GeneralT3EFTHandoffController import (
    BlockSpec, StageContract, EFTState, ThresholdEvent, Callback,
    EFTHandoffController, HandoffContractError, IncompleteCandidateInventory,
    candidate_inventory_from_stage2,
)
from scripts.AuditGeneralT3EFTHandoffReadiness import build_summary


C5 = BlockSpec('C5', (3, 3), 'complex', True, 'GeV^-1',
               operator_basis='synthetic-Weinberg-test-only',
               normalization='synthetic-test-only', mass_dimension=-1)
G = BlockSpec('g')


def contracts():
    return (
        StageContract('UV', ('fermion:2', 'scalar:S1'), (G, BlockSpec('MF', (1, 1), 'complex'))),
        StageContract('EFT1', ('scalar:S1',), (G,), (C5,)),
        StageContract('SM', (), (G,), (C5,)),
    )


def initial_state():
    return EFTState('SYNTHETIC', 'UV', 200.0, ('fermion:2', 'scalar:S1'),
                    {'g': 0.5, 'MF': [[100+1j]]})


def evolve(state, goal):
    # A deliberate test callback, NOT any T3 beta function.
    new = dict(state.parameters)
    new['g'] = float(np.asarray(new['g'])) + 0.001 * np.log(state.mu_gev/goal)
    return replace(state, mu_gev=goal, parameters=new)


def match_fermion(state, event):
    return EFTState(state.model, event.next_stage, state.mu_gev,
                    ('scalar:S1',), {'g': state.parameters['g']},
                    {'C5': np.eye(3, dtype=complex) * 1e-12})


def match_scalar(state, event):
    return EFTState(state.model, event.next_stage, state.mu_gev, (),
                    {'g': state.parameters['g']}, {'C5': state.wilson['C5']})


def make_controller(*, all_runners=False, fermion_match=False, scalar_match=False):
    obj = EFTHandoffController()
    for spec in contracts():
        obj.add_stage(spec)
    if all_runners:
        for stage in ('UV', 'EFT1', 'SM'):
            obj.add_runner(stage, Callback(stage+'-synthetic-running', evolve, 'unit-test'))
    if fermion_match:
        obj.add_matcher('UV', 'fermion:2', 'EFT1',
                        Callback('synthetic-F2-projection', match_fermion, 'unit-test'))
    if scalar_match:
        obj.add_matcher('EFT1', 'scalar:S1', 'SM',
                        Callback('synthetic-scalar-projection', match_scalar, 'unit-test'))
    return obj


def event_sequence():
    return (ThresholdEvent('fermion:2', 'fermion', 150., 'EFT1'),
            ThresholdEvent('scalar:S1', 'scalar', 110., 'SM'))


def audit_row(kind='mixed', missing=()):
    branches = {
        'mixed': ('scalar:Y', 'scalar:X'),
        'none': ('scalar:S1', 'scalar:S2'),
        'self_s1': ('scalar:S1_plus', 'scalar:S1_minus', 'scalar:S2'),
        'self_s2': ('scalar:S1', 'scalar:S2_plus', 'scalar:S2_minus'),
    }[kind]
    found_fermions = [f'fermion:{k}' for k in (2, 1, 0) if f'fermion:{k}' not in missing]
    found_scalar = [s for s in branches if s not in missing]
    events = ([{'particle': x, 'kind': 'fermion', 'scale_gev': 130000.-i*10000,
                'status': 'candidate_crossing_matching_required_not_applied'}
               for i, x in enumerate(('fermion:2','fermion:1','fermion:0')) if x in found_fermions]
              + [{'particle': x, 'kind': 'scalar', 'scale_gev': 1400.-j*100,
                  'status': 'matching_required_not_applied'} for j,x in enumerate(branches) if x in found_scalar])
    ordered = sorted(events, key=lambda x: -x['scale_gev'])
    return {'model':'TEST', 'status':'CandidateRunningFermionCrossingAuditNoMatching',
            'quadratic_kind':kind, 'uv_scalar_crossings':[],
            'uv_fermion_crossings':[x for x in events if x['kind']=='fermion'],
            'eft1_scalar_crossings':[x for x in events if x['kind']=='scalar'],
            'unresolved_fermion_crossings':[{'particle':x} for x in missing if x.startswith('fermion:')],
            'candidate_event_order':[{'particle':x['particle'],'scale_gev':x['scale_gev']}
                                     for x in ordered]}


class TestEFTHandoffController(unittest.TestCase):
    def test_complete_with_explicit_synthetic_callbacks(self):
        ctl = make_controller(all_runners=True, fermion_match=True, scalar_match=True)
        out = ctl.run(initial_state(), event_sequence(), final_scale_gev=50.)
        self.assertEqual(out.status, 'Complete')
        self.assertEqual(out.state.stage, 'SM')
        self.assertEqual(out.state.mu_gev, 50.)
        self.assertEqual([x.operation for x in out.history], ['run','match','run','match','run'])
        self.assertFalse(out.physical_matching_certified)
        self.assertEqual(tuple(out.state.wilson['C5'].shape), (3,3))

    def test_missing_uv_running_stops_before_any_progress(self):
        out = make_controller().run(initial_state(), event_sequence(), final_scale_gev=50.)
        self.assertEqual(out.status, 'Blocked')
        self.assertIn('Missing RG running', out.reason)
        self.assertEqual(out.state.mu_gev, 200.)
        self.assertFalse(out.history)

    def test_missing_matching_stops_at_threshold_and_preserves_uv(self):
        out = make_controller(all_runners=True).run(initial_state(), event_sequence(), final_scale_gev=50.)
        self.assertEqual(out.status, 'Blocked')
        self.assertIn('Missing matching', out.reason)
        self.assertEqual(out.state.stage, 'UV')
        self.assertEqual(out.state.mu_gev, 150.)
        self.assertEqual(len(out.history), 1)
        self.assertEqual(out.next_event.particle, 'fermion:2')

    def test_missing_later_matching_stops_after_only_applied_prior_callback(self):
        out = make_controller(all_runners=True, fermion_match=True).run(initial_state(), event_sequence(), final_scale_gev=50.)
        self.assertEqual(out.status, 'Blocked')
        self.assertEqual(out.state.stage, 'EFT1')
        self.assertEqual(out.state.mu_gev, 110.)
        self.assertEqual(len([x for x in out.history if x.operation=='match']), 1)

    def test_missing_intermediate_running_stops_after_validated_match(self):
        ctl = make_controller(fermion_match=True)
        ctl.add_runner('UV', Callback('uv-test-run',evolve,'test'))
        out = ctl.run(initial_state(),event_sequence(),final_scale_gev=50.)
        self.assertEqual(out.state.stage,'EFT1')
        self.assertEqual(out.state.mu_gev, 150.)
        self.assertIn('Missing RG running', out.reason)

    def test_wrong_active_particle_transition_cannot_register(self):
        ctl = make_controller()
        with self.assertRaises(HandoffContractError):
            ctl.add_matcher('UV', 'scalar:S1', 'EFT1', Callback('wrong',match_scalar,'test'))

    def test_bad_event_order_and_duplicate_fail(self):
        ctl = make_controller()
        for events in (
            (event_sequence()[1],event_sequence()[0]),
            (event_sequence()[0],event_sequence()[0]),
            (ThresholdEvent('fermion:2','fermion',210.,'EFT1'),),
            (ThresholdEvent('fermion:2','fermion',40.,'EFT1'),),
        ):
            with self.subTest(events=events), self.assertRaises(HandoffContractError):
                ctl.run(initial_state(), events, final_scale_gev=50.)

    def test_nearly_simultaneous_mass_roots_fail_closed(self):
        events=(ThresholdEvent('fermion:2','fermion',150.,'EFT1'),
                ThresholdEvent('scalar:S1','scalar',150.-1e-9,'SM'))
        with self.assertRaises(HandoffContractError):
            make_controller().run(initial_state(),events,final_scale_gev=50.)

    def test_bad_callback_endpoints_and_models_rejected(self):
        for callback in (
            lambda state, goal: replace(evolve(state, goal), mu_gev=goal*1.01),
            lambda state, goal: replace(evolve(state, goal), model='WRONG'),
            lambda state, goal: replace(evolve(state, goal), parameters={'g':1.0}),
        ):
            ctl=make_controller()
            ctl.add_runner('UV',Callback('bad',callback,'test'))
            with self.subTest(callback=callback),self.assertRaises(HandoffContractError):
                ctl.run(initial_state(),event_sequence(),final_scale_gev=50.)

    def test_bad_match_output_rejected(self):
        ctl=make_controller(all_runners=True)
        def bad_match(st,ev):
            after=match_fermion(st,ev)
            return replace(after, wilson={'C5':np.array([[1+1j,2],[3,4]])})
        ctl.add_matcher('UV','fermion:2','EFT1',Callback('bad',bad_match,'test'))
        with self.assertRaises(HandoffContractError):
            ctl.run(initial_state(),event_sequence(),final_scale_gev=50.)

    def test_wilson_c5_requires_explicit_basis_and_normalization(self):
        with self.assertRaises(HandoffContractError):
            StageContract('bad',(),(),(BlockSpec('C5',(3,3),'complex',True),))
        with self.assertRaises(HandoffContractError):
            StageContract('bad',(),(),(BlockSpec('C5',(3,3),'complex',False,
                                                   operator_basis='b',normalization='n',mass_dimension=-1),))

    def test_valid_c5_must_be_symmetric_not_hermitian(self):
        spec=StageContract('withC5',(),(),(C5,))
        sym=np.array([[1.,1j,0],[1j,2.,0],[0,0,3.]],dtype=complex)
        state=EFTState('TEST','withC5',100.,(),{}, {'C5':sym})
        validated=spec.validate(state)
        np.testing.assert_allclose(validated.wilson['C5'],sym)
        sym[0,1]=17
        self.assertNotEqual(validated.wilson['C5'][0,1],17)  # deep copy
        with self.assertRaises(HandoffContractError):
            spec.validate(EFTState('TEST','withC5',100.,(),{}, {'C5':np.array([[1.,1j,0],[-1j,2.,0],[0,0,3.]])}))

    def test_four_scalar_structures_and_binding_require_explicit_next_stages(self):
        for kind in ('none','mixed','self_s1','self_s2'):
            inv=candidate_inventory_from_stage2(audit_row(kind))
            with self.subTest(kind=kind):
                self.assertTrue(inv.complete)
                self.assertEqual(len(inv.events), 5 if kind in ('none','mixed') else 6)
                self.assertEqual(inv.events[0].particle,'fermion:2')
                with self.assertRaises(HandoffContractError):
                    inv.bind_stages({})
                mapping={e.particle:'STAGE_AFTER_'+e.particle for e in inv.events}
                self.assertEqual(len(inv.bind_stages(mapping)),len(inv.events))

    def test_unresolved_fermion_or_scalar_blocks_full_stage_binding(self):
        for kind,missing in (('mixed',('fermion:0',)),('self_s2',('scalar:S2_minus',))):
            inv=candidate_inventory_from_stage2(audit_row(kind, missing))
            with self.subTest(kind=kind):
                self.assertFalse(inv.complete)
                self.assertEqual(inv.unresolved,missing)
                with self.assertRaises(IncompleteCandidateInventory):
                    inv.bind_stages({e.particle:'AFTER' for e in inv.events})

    def test_candidate_cached_event_order_tampering_detected(self):
        row=audit_row('mixed')
        row['candidate_event_order'][0]['scale_gev']+=1.
        with self.assertRaises(HandoffContractError):
            candidate_inventory_from_stage2(row)

    def test_source_readiness_summary_contains_no_fake_matching(self):
        source={'format':'general-t3-running-heavy-mass-crossing-audit-v1',
                'status':'CandidateAuditComplete','physical_matching_certified':False,
                'models':[audit_row('mixed')]}
        out=build_summary(source,source_hash='abc')
        self.assertEqual(out['models_checked'],1)
        self.assertEqual(out['candidate_mass_roots_found'],5)
        self.assertEqual(out['matching_callbacks_registered'],0)
        self.assertFalse(out['physical_matching_certified'])
        self.assertEqual(out['models'][0]['status'],'MissingMatchingAndIntermediateRGEs')

    def test_fail_on_wrong_audit_status(self):
        row=audit_row('mixed')
        row['status']='Blocked'
        with self.assertRaises(HandoffContractError):
            candidate_inventory_from_stage2(row)


if __name__=='__main__':
    unittest.main()
