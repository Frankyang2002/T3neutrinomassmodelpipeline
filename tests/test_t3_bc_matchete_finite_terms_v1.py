"""Regression tests for the fail-closed Matchete local-amplitude partition."""
from __future__ import annotations

import json
from pathlib import Path
from tempfile import TemporaryDirectory
import unittest

from Numerical.experimental.T3BCMatcheteFiniteTerms import (
    AuditError, additive_terms, split_inputform,
)
from scripts.ExtractT3BCFiniteScalarTerms import process, run, MODELS, Y_DIR, X_DIR

H = 'Field[H, Scalar, {i}, {}]'
X = 'Field[NewScalarX, Scalar, {i}, {}]'
L = 'NCM[Bar[Field[l, Fermion, {i}, {}]], GammaCC]'
DENOM = 'Coupling[MF, {}, 0]'


def term(stage='Y', sector='LLHH', component='pole', multiplier='1'):
    scalar = {'LLHH': H+'*'+H, 'LLHX': H+'*'+X, 'LLXX': X+'*'+X}[sector]
    num = f'hbar*{multiplier}*{scalar}*{L}'
    if component == 'pole':
        return f'({num})/(ϵ*{DENOM})'
    if component == 'log':
        return f'({num}*Log[μbar2/Coupling[MS{stage}, {{}}, 0]^2])/{DENOM}'
    return f'({num})/{DENOM}'


class TestParser(unittest.TestCase):
    def test_top_level_split_preserves_nested_signs(self):
        raw = '-' + term() + ' + ' + term(component='log') + ' - ' + term(component='constant')
        out = additive_terms(raw)
        self.assertEqual(len(out), 3)
        self.assertTrue(out[0].startswith('-'))
        self.assertTrue(out[1].startswith('+'))
        self.assertTrue(out[2].startswith('-'))

    def test_equal_mass_formal_split(self):
        raw = '-' + term() + ' + ' + term(component='log') + ' + ' + term(component='constant')
        s = split_inputform(raw, stage='Y', sector='LLHH')
        self.assertEqual((len(s.pole_terms), len(s.log_terms), len(s.constant_terms)), (1,1,1))
        self.assertNotIn('ϵ', s.pole_residue)
        self.assertEqual(s.scale_equal_mass.replace(' ', ''), ('+'+term(component='constant')).replace(' ', ''))
        self.assertIn('Log[', s.pole_free)
        self.assertNotIn('Log[', s.scale_equal_mass)
        self.assertFalse(s.summary()['physical_matching_certified'])

    def test_x_no_constant(self):
        raw = term(stage='X') + ' + ' + term(stage='X', component='log')
        s = split_inputform(raw, stage='X', sector='LLHH')
        self.assertEqual(s.scale_equal_mass, '0')
        self.assertEqual(len(s.pole_terms), 1)

    def test_zero_llhx(self):
        s=split_inputform('0', stage='Y', sector='LLHX')
        self.assertEqual(s.pole_residue, '0')
        self.assertEqual(s.pole_free, '0')
        self.assertEqual(s.summary()['terms_total'], 0)

    def test_zero_llhh_rejected(self):
        with self.assertRaises(AuditError):
            split_inputform('0', stage='Y', sector='LLHH')

    def test_higher_pole_rejected(self):
        with self.assertRaises(AuditError):
            split_inputform(term().replace('ϵ*','ϵ^2*'), stage='Y', sector='LLHH')

    def test_epsilon_numerator_rejected(self):
        with self.assertRaises(AuditError):
            split_inputform(term().replace('/(ϵ*', '*ϵ*/('), stage='Y', sector='LLHH')

    def test_mixed_log_pole_rejected(self):
        with self.assertRaises(AuditError):
            split_inputform(term()+'*Log[μbar2/Coupling[MSY, {}, 0]^2]', stage='Y', sector='LLHH')

    def test_other_log_rejected(self):
        with self.assertRaises(AuditError):
            split_inputform(term(component='log').replace('MSY','MSX'), stage='Y', sector='LLHH')

    def test_wrong_scalar_structure_rejected(self):
        with self.assertRaises(AuditError):
            split_inputform(term(sector='LLXX'), stage='Y', sector='LLHH')

    def test_unsupported_stage(self):
        with self.assertRaises(AuditError):
            split_inputform(term(), stage='F', sector='LLHH')

    def test_unbalanced(self):
        with self.assertRaises(AuditError):
            split_inputform(term()+'(', stage='Y', sector='LLHH')

    def test_sqrt_denominator(self):
        raw=term().replace('/(ϵ*Coupling[MF', '/(Sqrt[3]*ϵ*Coupling[MF')
        s=split_inputform(raw, stage='Y', sector='LLHH')
        self.assertNotIn('ϵ',s.pole_residue)
        self.assertIn('Sqrt[3]',s.pole_residue)


class TestSourceAudit(unittest.TestCase):
    def _fixture(self, root: Path, *, tamper=False):
        for model in MODELS:
            short=model.replace('_-', '_m')
            yd=root/Y_DIR; xd=root/X_DIR
            yd.mkdir(parents=True, exist_ok=True)
            xd.mkdir(parents=True, exist_ok=True)
            ys=f'{short}_final_y_only_canonical_match'
            xs=f'{short}_final_x_from_saved_y_match'
            (yd / (ys+'.json')).write_text(json.dumps({
                'model':model,'status':'YOnlyIndependentCanonicalRawEFT2NotCertified',
                'tree_matching_completed': True, 'one_loop_matching_completed': True,
                'raw_EFT2_operator_sector_counts': {'LLHH':2,'LLHX':0,'LLXX':2},
            }),encoding='utf-8')
            (xd / (xs+'_x_followup.json')).write_text(json.dumps({
                'model': model, 'x_tree_matching_completed': True,
                'rebased_check_passed': True,'source_is_y_tree_EOMEFT':True,
            }),encoding='utf-8')
            (xd / (xs+'.log')).write_text('Matching pipeline completed successfully.\n'*2+
                'x-followup-one-loop-matching-returned\nx-followup-matching-exports-returned',
                encoding='utf-8')
            for sector,source in {
                'llhh':term()+' + '+term(component='log'),
                'llhx':'0',
                'llxx':term(sector='LLXX')+' + '+term(sector='LLXX',component='log'),
            }.items():
                (yd / (ys+f'_raw_{sector}.inputform.txt')).write_text(source,encoding='utf-8')
            (xd / (xs+'_x_followup_x_raw_llhh.inputform.txt')).write_text(
                ' + '.join(term(stage='X',component='pole',multiplier=str(i))
                           for i in range(1,37)),encoding='utf-8')
            if tamper:
                (yd / (ys+'_raw_llxx.inputform.txt')).write_text('0',encoding='utf-8')

    def test_end_to_end_both_models(self):
        with TemporaryDirectory() as tmp:
            root=Path(tmp);self._fixture(root)
            rows=run(root,MODELS,Y_DIR,X_DIR,root/'out')
            self.assertEqual(len(rows), 8)
            self.assertTrue((root/'out/summary.json').exists())
            self.assertFalse(json.loads((root/'out/summary.json').read_text())[
                'physical_sequential_matching_certified'])
            self.assertEqual(rows[0]['terms_total'],2)
            self.assertTrue((root/'out/B_m1/Y_LLHH_pole_residue.inputform.txt').is_file())

    def test_bad_source_inventory_fails_closed_before_outputs(self):
        with TemporaryDirectory() as tmp:
            root=Path(tmp);self._fixture(root,tamper=True)
            with self.assertRaises(AuditError):
                process('B_-1',root,Y_DIR,X_DIR,root/'out')
            self.assertFalse((root/'out/B_m1').exists())

    def test_missing_log_fails_closed(self):
        with TemporaryDirectory() as tmp:
            root=Path(tmp);self._fixture(root)
            (root/X_DIR/'B_m1_final_x_from_saved_y_match.log').unlink()
            with self.assertRaises(AuditError):
                process('B_-1',root,Y_DIR,X_DIR,root/'out')

if __name__ == '__main__':
    unittest.main()
