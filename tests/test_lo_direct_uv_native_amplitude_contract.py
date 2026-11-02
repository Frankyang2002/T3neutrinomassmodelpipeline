"""Source guard, SU2 epsilon and Majorana 1/2 negative controls."""
from pathlib import Path
import tempfile
import unittest

from Numerical.lo_certification.direct_uv.native_amplitude_contract import (
    CLASSES, GUARDS, SOURCES, _wl_without_comments, source_audit,
    epsilon_doublet_mass_insertion, majorana_singlet_mass_insertion,
    model_dictionary, report, ContractError,
)


class NativeAmplitudeContractTests(unittest.TestCase):
    def test_model_count_and_restrictions(self):
        rows = model_dictionary()
        self.assertEqual(len(rows), 16)
        self.assertEqual(sum(x['fermion_type'] == 'Majorana' for x in rows), 2)
        self.assertEqual(sum(x['extra_symmetry_allowed_yukawa'] is not None
                             for x in rows), 6)
        b = next(x for x in rows if x['model'] == 'B_-1')
        a = next(x for x in rows if x['model'] == 'A_+0')
        self.assertEqual(b['RGBeta_mass'], '{F,F}')
        self.assertEqual(a['RGBeta_SU2_mass'], 'epsilon')
        self.assertFalse(a['physical_component_normalisation_certified'])

    def test_epsilon_source_mass_amplitude(self):
        w = epsilon_doublet_mass_insertion()
        self.assertEqual(w['inverse_mixed_component_plus'], '1/7')
        self.assertEqual(w['inverse_mixed_component_minus'], '-1/7')
        self.assertTrue(w['fermion_mass_square_MF2_identity'])

    def test_majorana_source_mass_no_factor_two(self):
        w = majorana_singlet_mass_insertion()
        self.assertEqual(w['mixed_LLSS_coefficient_in_units_y1y2'], '1/7')
        self.assertNotEqual(w['mixed_LLSS_coefficient_in_units_y1y2'], '2/7')

    def test_invalid_source_mass_rejected(self):
        with self.assertRaises(ContractError):
            epsilon_doublet_mass_insertion(0)
        with self.assertRaises(ContractError):
            majorana_singlet_mass_insertion(0)

    def test_nested_mathematica_comments_not_evidence(self):
        sample = 'x (* example (* nested *) valid-looking-form *) z'
        self.assertEqual(_wl_without_comments(sample), 'xz')

    def test_source_integrity_guard_and_mutation_controls(self):
        # Minimal structural fixture, not a claim of native package execution.
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            for name, rel in SOURCES.items():
                p = root / rel
                p.parent.mkdir(parents=True, exist_ok=True)
                content = '\n'.join(GUARDS[name])
                if p.suffix == '.py':
                    # Separate every Python guard in a syntactically valid
                    # string literal so the tokenizer sees its contents.
                    content = '\n'.join(repr(x) for x in GUARDS[name])
                p.write_text(content + '\n', encoding='utf-8')
            self.assertEqual(source_audit(root)['status'], 'source_forms_verified')
            self.assertEqual(report(root)['physical_source_normalisation_certified_count'], 0)
            # Wrong scalar orientation must fail, not quietly alter the map.
            p = root / SOURCES['topology']
            p.write_text(p.read_text().replace('scalarBar=which===2;',
                                                'scalarBar=which===1;'))
            self.assertEqual(source_audit(root)['status'], 'blocked_source_forms')
            with self.assertRaises(ContractError):
                report(root)
            # Correct scalar but wrong Dirac mass epsilon must fail closed.
            p.write_text(p.read_text().replace('scalarBar=which===1;',
                                                'scalarBar=which===2;'))
            p = root / SOURCES['rgbeta']
            p.write_text(p.read_text().replace('dF===2,(eps[', 'dF===2,(del['))
            with self.assertRaises(ContractError):
                report(root)


if __name__ == '__main__':
    unittest.main()

class ExistingFrozenInputShapeTests(unittest.TestCase):
    def test_a0_and_bminus1_readonly_frozen_tensor_audit(self):
        # All 16 are generated ONLY as synthetic test shapes. The actual local
        # frozen benchmark data are never claimed to have been run here.
        import json
        import numpy as np
        from Numerical.lo_certification.direct_uv.native_amplitude_contract import audit_frozen_inputs
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp)
            d=root/'output/general_observable/lo_direct_one_loop/provisional_inputs'
            d.mkdir(parents=True)
            for row in model_dictionary():
                ds1,ds2,df=(row['SU2_dimensions'][j] for j in ('S1','S2','F'))
                maj=row['fermion_type']=='Majorana'
                half=3*df
                nf=half if maj else 2*half
                ns=ds1+ds2
                M=np.eye(nf,dtype=complex)*7 if maj else np.block([
                    [np.zeros((half,half)),np.eye(half)*7],
                    [np.eye(half)*7,np.zeros((half,half))]])
                Y1=np.zeros((3,nf,ns),complex)
                Y2=np.zeros((3,nf,ns),complex)
                Y1[:,0 if maj else half,0]=0.1
                Y2[:,0,ds1]=0.2
                def enc(x):return {'real':x.real.tolist(),'imag':x.imag.tolist()}
                (d/(row['model']+'.json')).write_text(json.dumps({
                    'model':row['model'], 'fermion_mass_gev':enc(M),
                    'scalar_A_gev2':enc(np.eye(ns,dtype=complex)),
                    'scalar_B_gev2':enc(np.zeros((ns,ns),complex)),
                    'vertices':[{'name':'y1_HC','conjugate':enc(Y1)},
                                {'name':'y2_HC','field':enc(Y2)}]}))
            result=audit_frozen_inputs(root)
            self.assertEqual(result['models_checked'],16)
            self.assertFalse(result['physical_source_mapping_certified'])
            a=d/'A_+0.json'
            obj=json.loads(a.read_text())
            obj['vertices'][1]['name']='y1_HC'  # cannot accidentally pass wrong pair
            a.write_text(json.dumps(obj))
            with self.assertRaises(ContractError):
                audit_frozen_inputs(root)
