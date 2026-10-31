import sys
import unittest
from pathlib import Path
from unittest.mock import patch

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
from scripts import RunT3BCSameTrajectoryWilson as subject

class TestSameTrajectory(unittest.TestCase):
    def test_audit_accepts_pass(self):
        subject.check_audit({'status':'ConditionalUVToEFT1SourceBoundaryAudited',
                             'physical_C5_certified':False})
    def test_audit_rejects_missing(self):
        with self.assertRaises(ValueError):subject.check_audit({'status':'Wrong'})
    def test_audit_rejects_physical_claim(self):
        with self.assertRaises(ValueError):subject.check_audit({
            'status':'ConditionalUVToEFT1SourceBoundaryAudited','physical_C5_certified':True})
    def test_interval_gate(self):
        with self.assertRaises(ValueError):subject.require(False,'expected')
    def test_only_one_eft(self):
        self.assertEqual(subject.MODELS,('B_-1','C_-1'))
        self.assertIn('SameTrajectory',subject.STATUS)

if __name__=='__main__':unittest.main()
