import json
import tempfile
import unittest
from pathlib import Path
from Numerical.lo_certification.completeness.inventory import ALL_MODELS,build_ledger,STATUSES,ENTRIES
from scripts.RunLOCompletenessAudit import main,compute_conditional

class InventoryTests(unittest.TestCase):
    def test_all_16_and_explicit_statuses(self):
        x=build_ledger()
        self.assertEqual(len(ALL_MODELS),16)
        self.assertEqual(x['n_rows'],16*len(ENTRIES))
        self.assertTrue({r['status'] for r in x['rows']}<=STATUSES)
        self.assertFalse(x['fully_physics_complete'])
        self.assertEqual(len({(r['model'],r['contribution']) for r in x['rows']}),len(x['rows']))
        self.assertIn('B_-1',ALL_MODELS)
        self.assertIn('E_+2',ALL_MODELS)

    def test_no_claim_of_complete_hard_matching(self):
        rows=build_ledger()['rows']
        hard=[r for r in rows if r['contribution']=='F_THRESHOLD_HARD_LLHH']
        self.assertTrue(all(r['status']=='MISSING' for r in hard))
        self.assertEqual(sum(r['raw_hard_pole_log_factorization_available'] for r in hard),2)
        self.assertTrue(all(not r['physical_certified'] for r in rows))
        extra=[r for r in rows if r['contribution']=='ONE_LOOP_HIGHER_DIMENSION_POWER_TERMS']
        self.assertTrue(all(r['status']!='HIGHER_ORDER' for r in extra))

    def test_local_inspection_only_and_no_overwrite(self):
        with tempfile.TemporaryDirectory() as tmp:
            d=Path(tmp)
            self.assertFalse((d/'dummy').exists())
            ledger=build_ledger(d)
            self.assertFalse(ledger['sources']['tree']['present_locally'])
            self.assertFalse((d/'dummy').exists())

    def test_six_opt_in_neutral_conjugate_yukawas(self):
        x=build_ledger()
        r=[i for i in x['rows'] if i['contribution']=='EXTRA_NEUTRAL_SCALAR_CONJUGATE_YUKAWA']
        self.assertEqual(sum(v['status']=='UNRESOLVED' for v in r),6)
        self.assertEqual(sum(v['status']=='ZERO_BY_SYMMETRY' for v in r),10)
        self.assertEqual({v['model'] for v in r if v['status']=='UNRESOLVED'},
                         {'A_-2','A_+0','D_-2','D_+0','E_-2','E_+0'})

    def test_missing_hard_input_is_not_zero(self):
        with tempfile.TemporaryDirectory() as tmp:
            v=compute_conditional('A_-4',Path(tmp))
            self.assertEqual(v['status'],'MISSING_INPUT')
            self.assertFalse(v['physical_C5_certified'])

if __name__=='__main__': unittest.main()
