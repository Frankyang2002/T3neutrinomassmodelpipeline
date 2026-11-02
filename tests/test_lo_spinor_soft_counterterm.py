"""Fail-closed tests for native F-only pole/soft/CT preparation."""
from __future__ import annotations

import copy
from hashlib import sha256
import json
import math
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
from unittest.mock import patch

from Numerical.lo_certification.spinor.boundary import MODELS, BoundaryError, tag
from Numerical.lo_certification.spinor.soft_counterterm import (
    FORMAT, EXPECTED_TOTALS, SOURCE_SEMANTICS, POLE_CONVENTION,
    bubble_b0_finite, scalar_master_table, _validate_mass_basis,
    validate_bridge_snapshot, build_readiness,
)

ROOT = Path(__file__).resolve().parents[1]
COUNTS = {'A_-4':2, 'A_-2':4, 'A_+0':4,
          'B_-3':2, 'B_-1':68, 'B_+1':2,
          'C_-3':2, 'C_-1':68, 'C_+1':2,
          'D_-2':4, 'D_+0':4, 'D_+2':2,
          'E_-4':2, 'E_-2':4, 'E_+0':4, 'E_+2':2}


def make_all(root: Path):
    """Synthetic complete 16-model bridge records for output-contract tests."""
    plans = []
    model_results = {}
    for m in MODELS:
        pair = []
        for i in range(COUNTS[m]):
            chirality = 'PL' if i < COUNTS[m]//2 else 'PR'
            text = f'+ (hbar*Coupling[y1, {{}}, 0]*Coupling[MF, {{}}, 0]^-1*GammaCC*{tag(m)}*{i})'
            pair.append({'chirality': chirality,
                         'kernel_inputform': text,
                         'kernel_sha256': sha256(text.encode()).hexdigest(),
                         'mass_log_inputform': 'Log[μbar2/Coupling[MF, {}, 0]^2]',
                         'source_indices': {'constant':3*i+1,'pole':3*i+2,'scale_log':3*i+3}})
        model_results[m] = {'records': pair, 'PL_groups': COUNTS[m]//2, 'PR_groups': COUNTS[m]//2, 'source_inventory_terms_sha256':'a'*64}
        plan = {'model': m, 'mass_basis_input': None, 'mass_basis_input_sha256': None}
        if m in ('B_-1', 'C_-1'):
            path = root/'inputs'/f'{tag(m)}.json'
            path.parent.mkdir(parents=True, exist_ok=True)
            data = {'model':m, 'format':'t3-experimental-mass-eigen-matchete-parameters-v1',
                    'status':'NumericMassEigenstateBoundaryPreparedMatcheteNotRun',
                    'heavy_mass_squared_gev2':[1.0e6,4.0e6]}
            path.write_text(json.dumps(data), encoding='utf8')
            plan['mass_basis_input'] = str(path.relative_to(root)).replace('/','\\')
            plan['mass_basis_input_sha256'] = sha256(path.read_bytes()).hexdigest()
        plans.append(plan)
    bridge = {'source_release_summary_sha256':'e515d82d118284ba'+'0'*48,
              'mu_F_gev':1.0e5}
    native = {'model_plans': plans}
    return bridge, model_results, native


class ScalarB0Tests(unittest.TestCase):
    def test_equal_masses(self):
        self.assertAlmostEqual(bubble_b0_finite(4e6,4e6,1e10),
                               -math.log(4e6/1e10), places=13)

    def test_mass_permutation(self):
        self.assertAlmostEqual(bubble_b0_finite(1.2e6,4e6,1e10),
                               bubble_b0_finite(4e6,1.2e6,1e10))

    def test_zero_mass_single_leg(self):
        self.assertAlmostEqual(bubble_b0_finite(0,1.2e6,1e10),
                               1-math.log(1.2e6/1e10), places=13)

    def test_zero_zero_rejected_as_ir_singular(self):
        with self.assertRaisesRegex(BoundaryError, 'infrared'):
            bubble_b0_finite(0,0,1e10)

    def test_near_degenerate_stable(self):
        a=4.0e6
        b=a*(1-1e-12)
        ref=-math.log(a/1e10)+(a-b)/(2*a)
        self.assertAlmostEqual(bubble_b0_finite(a,b,1e10), ref, delta=2e-12)

    def test_nontrivial_quadrature(self):
        a,b,mu2=1.919e6,2.557e6,1e10
        n=4096
        # Midpoint quadrature of -int_0^1 log[(xa+(1-x)b)/mu2] dx.
        quadrature=sum(-math.log((b+(a-b)*(i+.5)/n)/mu2) for i in range(n))/n
        self.assertAlmostEqual(bubble_b0_finite(a,b,mu2), quadrature, delta=1e-9)

    def test_scale_derivative_wrt_ln_mubar_squared(self):
        a,b,mu2=2.0e6,3.0e6,1e10
        delta=0.075
        self.assertAlmostEqual(bubble_b0_finite(a,b,mu2*math.exp(delta))-
                               bubble_b0_finite(a,b,mu2), delta, places=12)

    def test_finite_pole_normalization_not_operator_matching(self):
        t=scalar_master_table([2e6,3e6],1e5,'B_-1','a'*64)
        self.assertEqual([x['scalar_pair'] for x in t['scalar_pairs']], ['XX','XY','YY'])
        self.assertTrue(all(x['B0_1_over_epsilon_bar_pole']==1 for x in t['scalar_pairs']))
        self.assertFalse(t['EFT1_soft_diagram_calculated'])
        self.assertFalse(t['loop_normalization_and_operator_symmetry_factors_applied'])

    def test_bad_scale_or_masses_fail(self):
        for x,y,mu in [(-1,4,1e10),(4,4,0),(4,4,-1),(float('nan'),4,1e10),
                       (4,float('inf'),1e10),(True,4,1e10),(4,4,True)]:
            with self.subTest(x=x,y=y,mu=mu),self.assertRaises(BoundaryError):
                bubble_b0_finite(x,y,mu)

    def test_no_sneakily_false_schema(self):
        with self.assertRaisesRegex(BoundaryError,'Two scalar'):
            scalar_master_table([1.2e6],1e5,'C_-1','f'*64)


class ReadinessTests(unittest.TestCase):
    def test_all16_source_resolved_and_physical_blocked(self):
        with tempfile.TemporaryDirectory() as td:
            root=Path(td); b,models,native=make_all(root)
            record,details=build_readiness(root,b,models,native,'b'*64,'c'*64)
            self.assertEqual(record['totals'], {'kernels':176,'PL':88,'PR':88})
            self.assertEqual(len(details['ledgers']),16)
            self.assertEqual(len(details['scalar_masters']),16)
            self.assertEqual(record['physical_C5'],None)
            self.assertFalse(record['EFT1_soft_loop_matching_amplitude_computed'])
            self.assertFalse(record['physical_MSbar_hard_matching_certified'])
            self.assertIsNone(record['native_one_loop_source_already_soft_subtracted'])
            self.assertIsNone(record['absolute_matchete_gamma_cc_to_master_weyl_factor'])
            self.assertEqual(record['native_loop_source_interpretation'],SOURCE_SEMANTICS)
            self.assertEqual(record['conditional_pole_balance_equation'],POLE_CONVENTION)
            self.assertEqual(sum(len(r['records']) for r in details['ledgers'].values()),176)
            for m in MODELS:
                for item in details['ledgers'][m]['records']:
                    self.assertEqual(item['raw_pole_residue_ratio_in_native_kernel_basis'],'1')
                    self.assertIsNone(item['external_EFT1_soft_pole_residue_inputform'])
                    self.assertIsNone(item['native_counterterm_pole_inputform'])
                    self.assertIsNone(item['physical_hard_one_loop_wilson'])
            self.assertEqual(len(details['scalar_masters']['B_-1']['scalar_pairs']),3)
            self.assertIsNone(details['scalar_masters']['A_-4']['scalar_masses_squared_gev2'])

    def test_mass_basis_hash_tamper_blocks(self):
        with tempfile.TemporaryDirectory() as td:
            root=Path(td); b,models,native=make_all(root)
            plan=next(x for x in native['model_plans'] if x['model']=='B_-1')
            path=root/plan['mass_basis_input'].replace('\\','/')
            path.write_text(path.read_text()+' ',encoding='utf8')
            with self.assertRaisesRegex(BoundaryError, 'hash mismatch'):
                build_readiness(root,b,models,native,'b'*64,'c'*64)

    def test_mass_basis_wrong_model_blocks(self):
        with tempfile.TemporaryDirectory() as td:
            root=Path(td); b,models,native=make_all(root)
            p=next(x for x in native['model_plans'] if x['model']=='B_-1')
            path=root/p['mass_basis_input'].replace('\\','/')
            data=json.loads(path.read_text());data['model']='C_-1'
            path.write_text(json.dumps(data),encoding='utf8')
            p['mass_basis_input_sha256']=sha256(path.read_bytes()).hexdigest()
            with self.assertRaisesRegex(BoundaryError,'schema'):
                build_readiness(root,b,models,native,'b'*64,'c'*64)

    def test_outside_project_mass_basis_path_blocks(self):
        with tempfile.TemporaryDirectory() as td:
            root=Path(td)/'repo';root.mkdir()
            p=Path(td)/'out.json';p.write_text('{}',encoding='utf8')
            plan={'mass_basis_input':'../out.json','mass_basis_input_sha256':sha256(p.read_bytes()).hexdigest()}
            with self.assertRaisesRegex(BoundaryError,'outside project'):
                _validate_mass_basis(root,plan,'B_-1',1e5)

    def test_wrong_number_of_native_plans(self):
        with tempfile.TemporaryDirectory() as td:
            root=Path(td); b,models,native=make_all(root)
            native['model_plans'].pop()
            with self.assertRaisesRegex(BoundaryError,'16 model'):
                build_readiness(root,b,models,native,'b'*64,'c'*64)

    def test_changed_kernel_detected(self):
        with tempfile.TemporaryDirectory() as td:
            root=Path(td); b,models,native=make_all(root)
            models['A_-4']['records'][0]['kernel_inputform']='different'
            with self.assertRaisesRegex(BoundaryError,'Modified native'):
                build_readiness(root,b,models,native,'b'*64,'c'*64)

    def test_wrong_mu_f_blocks(self):
        with tempfile.TemporaryDirectory() as td:
            root=Path(td); b,models,native=make_all(root)
            b['mu_F_gev']=0
            with self.assertRaisesRegex(BoundaryError,'threshold'):
                build_readiness(root,b,models,native,'b'*64,'c'*64)

    def test_no_duplicate_model_status(self):
        with tempfile.TemporaryDirectory() as td:
            root=Path(td); b,models,native=make_all(root)
            native['model_plans'][0]['model']='B_-1'
            with self.assertRaisesRegex(BoundaryError,'16 model'):
                build_readiness(root,b,models,native,'b'*64,'c'*64)

    def test_missing_bridge_blocks(self):
        with tempfile.TemporaryDirectory() as td:
            root=Path(td); (root/'output/general_observable/lo_certification_spinor').mkdir(parents=True)
            owned=root/'output/general_observable/lo_certification_spinor'
            with self.assertRaisesRegex(BoundaryError,'existing owned'):
                validate_bridge_snapshot(root,owned/'n',owned/'i',owned/'b')

    def test_cli_rejects_unsafe_runid_before_writes(self):
        cp=subprocess.run([sys.executable,str(ROOT/'scripts/BuildLOSpinorSoftSubtractionLedger.py'),
                           '--repo-root',str(ROOT),'--run-id','../evil', '--check-only'],
                           capture_output=True, text=True, encoding='utf8')
        self.assertEqual(cp.returncode,2)
        self.assertIn('Unsafe run-id',cp.stderr)

    def test_format_labels_nonphysical(self):
        self.assertIn('readiness',FORMAT)
        self.assertEqual(EXPECTED_TOTALS['kernels'],176)
        self.assertIn('CONDITIONAL',POLE_CONVENTION)
        self.assertIn('not independently verified',SOURCE_SEMANTICS)


if __name__=='__main__':
    unittest.main()

# Intentional standalone fixtures for byte-for-byte bridge provenance checks.
def make_mock_bridge_tree(root: Path):
    from Numerical.lo_certification.spinor.pole_log_bridge import expressions_for
    base = root/'output/general_observable/lo_certification_spinor'
    native = base/'native'; inventory = base/'inventory'; bridge = base/'bridge'
    for d in (native, inventory, bridge): d.mkdir(parents=True)
    _, models, _ = make_all(root)
    summary={'format':'t3-f-first-source-resolved-pole-log-bridge-v1',
             'status':'ALL16_EXACT_NATIVE_LAURENT_TRIPLETS_AUDITED_PHYSICAL_BLOCKED',
             'totals': dict(EXPECTED_TOTALS), 'source_release_summary_sha256':'e515d82d118284ba'+'0'*48,
             'mu_F_gev':100000.0}
    (bridge/'summary.json').write_text(json.dumps(summary),encoding='utf8')
    (native/'summary.json').write_text(json.dumps({
        'format':'t3-f-only-raw-one-loop-native-evidence-v1',
        'n_models_with_raw_native_one_loop_evidence':16,
        'hard_one_loop_MSbar_certified':False,
        'native_GammaCC_to_master_Weyl_phase':None
    }), encoding='utf8')
    for m in MODELS:
        d=bridge/tag(m);d.mkdir()
        (d/'trios.json').write_text(json.dumps(models[m]),encoding='utf8')
        for ch in ('PL','PR'):
            for kind, e in expressions_for(models[m],ch).items():
                (d/f'{ch}_{kind}.inputform.txt').write_text(e+'\n',encoding='utf8')
    return native, inventory, bridge, summary, models


class IndependentBridgeProvenanceTests(unittest.TestCase):
    def test_all16_bridge_exports_exactly_rechecked(self):
        with tempfile.TemporaryDirectory() as td:
            root=Path(td)
            n,i,b,summary,models=make_mock_bridge_tree(root)
            with patch('Numerical.lo_certification.spinor.soft_counterterm.audit_all',
                       return_value=(summary,models)):
                got, checked, native=validate_bridge_snapshot(root,n,i,b)
            self.assertEqual(got,summary)
            self.assertEqual(len(checked),16)
            self.assertEqual(native['n_models_with_raw_native_one_loop_evidence'],16)

    def test_expression_tampering_blocks(self):
        with tempfile.TemporaryDirectory() as td:
            root=Path(td)
            n,i,b,summary,models=make_mock_bridge_tree(root)
            path=b/'A_m4/PL_pole_residue_raw.inputform.txt'
            path.write_text(path.read_text() + '+ 123\n',encoding='utf8')
            with patch('Numerical.lo_certification.spinor.soft_counterterm.audit_all',
                       return_value=(summary,models)):
                with self.assertRaisesRegex(BoundaryError,'expression mismatch'):
                    validate_bridge_snapshot(root,n,i,b)

    def test_trios_tampering_blocks(self):
        with tempfile.TemporaryDirectory() as td:
            root=Path(td)
            n,i,b,summary,models=make_mock_bridge_tree(root)
            path=b/'C_m1/trios.json'
            obj=json.loads(path.read_text());obj['records'][0]['kernel_sha256']='0'*64
            path.write_text(json.dumps(obj),encoding='utf8')
            with patch('Numerical.lo_certification.spinor.soft_counterterm.audit_all',
                       return_value=(summary,models)):
                with self.assertRaisesRegex(BoundaryError,'trios mismatch'):
                    validate_bridge_snapshot(root,n,i,b)

    def test_native_uncertified_status_enforced(self):
        with tempfile.TemporaryDirectory() as td:
            root=Path(td)
            n,i,b,summary,models=make_mock_bridge_tree(root)
            path=n/'summary.json';obj=json.loads(path.read_text())
            obj['hard_one_loop_MSbar_certified']=True
            path.write_text(json.dumps(obj),encoding='utf8')
            with patch('Numerical.lo_certification.spinor.soft_counterterm.audit_all',
                       return_value=(summary,models)):
                with self.assertRaisesRegex(BoundaryError,'Native status unexpectedly'):
                    validate_bridge_snapshot(root,n,i,b)
