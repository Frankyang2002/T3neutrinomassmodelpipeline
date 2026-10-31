"""Offline source-adaptation tests; these do NOT run Matchete/Wolfram."""
from __future__ import annotations
import tempfile
import unittest
from pathlib import Path

from Numerical.experimental.T3BCYOnlyMatcheteBuilder import (
    BASE_NAME, GENERATED_NAME, adapt, build, replace_once,
)


def fixture():
    # The guarded anchors below mirror the inspected native probe. Filler only
    # meets the explicit large-native-program precondition; NOT a Matchete file.
    return ('(* source fixture *)\n' + 'x'*81000 + '\n'
        + 'T3MELeadingC5OneLoopProjection;\n'
        + 'T3MEVerifiedRenormalizableReconstruction;\n'
        + 'RunT3Matching[scalarEFT,5,1];\n'
        + 'checkHeavyX=Quiet@Check[\n'
        + '   DefineField[NewScalarX,Scalar,SelfConjugate->False,\n'
        + '     Indices->{T3BSMd2},Charges->{U1Y[-1/2]},Mass->{Heavy,MSX}];True,$Failed];\n'
        + 'checkHeavyY=Quiet@Check[\n'
        + '    DefineField[NewScalarY,Scalar,SelfConjugate->False,\n'
        + '     Indices->{T3BSMd2},Charges->{U1Y[-1/2]},Mass->{Heavy,MSY}];True,$Failed];\n'
        + 'StringContainsQ[completeText,"Field[NewScalarX"] ||\n'
        + '    StringContainsQ[completeText,"Field[NewScalarY"] ||\n'
        + '{"Field[NewScalarX","Field[NewScalarY","Field[NewFermion"}\n'
        + '"heavy_scalar_fields"->{"NewScalarX","NewScalarY"}\n'
        + '"no_scalar_fields_in_matched_eft"->True,\n'
        + '"excluded_four_heavy_scalar_quartic_terms"->Lookup[leadingC5Partition,\n'
        + ' T3MECheckpoint["scalar-threshold-EFT1-lagrangian-check-started"];\n'
        + ' <|\n   "format"->"t3-mass-eigen-eft1-scalar-threshold-probe-v1",\n'
        + '"raw_scalar_llhh_term_count"->Length[scalarCandidates],\n'
        + 'ScalarThresholdOneLoopRawSectorsExtractedNoEFT1Running\n'
        + 'ScalarThresholdOneLoopRawSectorsExtractedNoEFT1Running\n')


class TestYOnlyMatcheteGenerator(unittest.TestCase):
    def test_generates_full_source(self):
        result=adapt(fixture())
        self.assertIn('Mass->{Light,MSX}',result)
        self.assertIn('Mass->{Heavy,MSY}',result)
        self.assertIn('"heavy_scalar_fields"->{"NewScalarY"}',result)
        self.assertIn('"light_scalar_fields"->{"H","NewScalarX"}',result)
        self.assertIn('LLXX',result)
        self.assertIn('LLHX',result)
        self.assertIn('LLHH',result)
        self.assertIn('YOnlyOneLoopRawEFT2SectorsExtractedNotPhysicalMatching',result)
        self.assertNotIn('ScalarThresholdOneLoopRawSectorsExtractedNoEFT1Running',result)
        self.assertIn('"physical_matching_performed"',result) if '"physical_matching_performed"' in fixture() else None

    def test_only_y_is_heavy(self):
        s=adapt(fixture())
        self.assertNotIn('Mass->{Heavy,MSX}];True,$Failed]',s)
        self.assertIn('!StringContainsQ[completeText,"Field[NewScalarX"]',s)
        self.assertIn('StringContainsQ[completeText,"Field[NewScalarY"]',s)

    def test_preserves_native_matchete_internals(self):
        s=adapt(fixture())
        self.assertIn('RunT3Matching[scalarEFT,5,1]',s)
        self.assertIn('T3MEVerifiedRenormalizableReconstruction',s)
        self.assertIn('T3MELeadingC5OneLoopProjection',s)

    def test_drift_guard_missing_anchor(self):
        with self.assertRaisesRegex(ValueError,'retain-X-as-light'):
            adapt(fixture().replace('Mass->{Heavy,MSX}','Mass->{Light,MSX}'))

    def test_drift_guard_duplicate(self):
        p='"no_scalar_fields_in_matched_eft"->True,'
        with self.assertRaisesRegex(ValueError,'report-output-field-inventory'):
            adapt(fixture()+p)

    def test_drift_guard_status(self):
        with self.assertRaisesRegex(ValueError,'stage-status count'):
            adapt(fixture().replace('ScalarThresholdOneLoopRawSectorsExtractedNoEFT1Running','WRONG',1))

    def test_short_source_rejected(self):
        with self.assertRaisesRegex(ValueError,'short'):
            adapt(fixture()[:1000])

    def test_safe_replacement_helper(self):
        self.assertEqual(replace_once('abc','b','B','test'),'aBc')
        with self.assertRaises(ValueError):replace_once('bbb','b','B','test')

    def test_file_builder_preserves_input(self):
        with tempfile.TemporaryDirectory() as td:
            src=Path(td)/BASE_NAME;target=Path(td)/GENERATED_NAME
            original=fixture();src.write_text(original,encoding='utf8')
            result=build(src,target)
            self.assertEqual(src.read_text(),original)
            self.assertTrue(target.is_file())
            self.assertFalse(result['physical_matching_certified'])
            self.assertEqual(build(src,target)['generated_source_sha256'],result['generated_source_sha256'])

    def test_refuse_to_overwrite_modified_generated_file(self):
        with tempfile.TemporaryDirectory() as td:
            src=Path(td)/BASE_NAME;target=Path(td)/GENERATED_NAME
            src.write_text(fixture(),encoding='utf8');build(src,target)
            target.write_text('accidental overwrite')
            with self.assertRaisesRegex(ValueError,'differs'):
                build(src,target)

    def test_wrong_generated_name(self):
        with tempfile.TemporaryDirectory() as td:
            src=Path(td)/BASE_NAME;src.write_text(fixture())
            with self.assertRaisesRegex(ValueError,'basename'):
                build(src,Path(td)/'other.wl')

    def test_wrong_directory(self):
        with tempfile.TemporaryDirectory() as td:
            src=Path(td)/BASE_NAME;src.write_text(fixture())
            with self.assertRaisesRegex(ValueError,'must be next'):
                build(src,Path(td)/'sub'/GENERATED_NAME)


if __name__=='__main__':unittest.main()
