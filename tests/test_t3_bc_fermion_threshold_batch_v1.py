import unittest
from scripts import RunT3BCFermionThresholdBatch as m

class TestFThreshold(unittest.TestCase):
    def test_native_certificates(self):
        tree={'model':'B_-1','status':'FermionFirstTreeLLSSExtractedNotSequentialEFT1',
              'fermion_field_absent':True,'x_scalar_present':True,'y_scalar_present':True,
              'physical_matching_performed':False,'matchete_rgbeta_quartic_normalization_certified':False,
              'tree_matching_completed':True,'llss_term_count':12}
        loop=dict(tree, **{'status':'FermionFirstOneLoopRawSectorsExtractedNotRenormalized',
                  'one_loop_matching_completed':True,'raw_hard_llhh_scheme_renormalized':False,
                  'raw_hard_llhh_pole_rge_certified':False,
                  'raw_hard_llhh_candidate_term_count':204})
        from pathlib import Path
        import tempfile
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/'a';p.write_text('x')
            out=m.validate_one('B_m1',tree,loop,(p,p))
            self.assertFalse(out['physical_C5_certified'])
            self.assertTrue(out['F_integrated_out'])
            loop['raw_hard_llhh_scheme_renormalized']=True
            with self.assertRaises(ValueError):m.validate_one('B_m1',tree,loop,(p,p))

    def test_formal_sectors(self):
        def row(model):
            x={'conventional_C5_2x2_without_hbar_MF':[['a','b'],['b','c']]}
            return {'model':model,'raw_sector_counts':{'pole':68,'log':68,'constant':68},
              'pole_equals_MF_log_coefficient':True,
              'F_to_EFT1_independent_counterterms_calculated':False,
              'F_hard_finite_physical_matching_certified':False,'source_sha256':'xyz',
              'sectors':{k:x for k in ('pole','log','constant')}}
        r={'format':'t3-bc-F-hard-formal-boundary-plus-Y-conditional-ledger-v1',
           'full_physical_matching_certified':False,'physical_C5':None,
           'F_hard_constant_and_Y_conditional_are_SEPARATE_boundaries':True,
           'models':{'B_m1':row('B_-1'),'C_m1':row('C_-1')}}
        self.assertFalse(m.check_formal(r)['B_-1']['renormalized_finite_boundary'])
        r['models']['B_m1']['raw_sector_counts']['constant']=67
        with self.assertRaises(ValueError):m.check_formal(r)

if __name__=='__main__':unittest.main()
