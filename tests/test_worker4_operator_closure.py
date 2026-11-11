"""Production-source contract checks, run from repo root."""
import importlib.util
from pathlib import Path
import unittest

SCRIPT=Path('scripts/Worker4OperatorClosureAudit.py')
spec=importlib.util.spec_from_file_location('worker4_operator_audit',SCRIPT)
mod=importlib.util.module_from_spec(spec)
spec.loader.exec_module(mod)

class Worker4Tests(unittest.TestCase):
    def test_live_seven_term_source(self):
        data=mod.source_inventory()
        self.assertEqual(tuple(data['source_term_names']),mod.TERMS)
        self.assertTrue(data['quark_gamma_added_to_scalar_anomalous_dimension'])
        self.assertTrue(data['full_45_weyl_dictionary_exists'])
        self.assertTrue(data['simultaneous_solver_unchanged'])
    def test_scalar_pair_feedback_witness(self):
        data=mod._witnesses()
        self.assertGreater(data['omitted_scalar_pair_into_retained']['scalar_pair'],0)
        self.assertEqual(data['omitted_scalar_pair_into_retained']['yukawa_wavefunction'],0)
    def test_fermion_selection_witness(self):
        data=mod._witnesses()
        self.assertGreater(sum(data['omitted_fermion_index_into_retained'].values()),0)
    def test_reference_missing_blocks_flagged(self):
        self.assertEqual({x[1] for x in mod.OTHER_CLASSES},
                         {'4.84','4.86','4.87','4.88','4.89','4.90','4.91','4.92','4.93'})
if __name__=='__main__':unittest.main()
