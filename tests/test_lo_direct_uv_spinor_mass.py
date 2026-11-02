"""Fail-closed source-MF Ward and conditional y1/y2/MF/lambda phase tests."""
from __future__ import annotations
from pathlib import Path
import json
import subprocess
import sys
import tempfile
import unittest

from Numerical.lo_certification.direct_uv.spinor_mass_dictionary import (
    DictionaryError, DIRAC, MAJORANA, audit, guard_sources, mass_ward,
    optional_evidence, phase_sum, spurion_witness, without_wl_comments,
)

RGB_TEXT='''
T3RGBetaRep[1] := None;
T3RGBetaRep[2] := SU2L[fund];
T3RGBetaRep[3] := SU2L[S2];
T3RGBetaAddHeavyFermion[dF_Integer, yF_] := Module[{},
 AddFermion[F, GaugeRep -> T3RGBetaGaugeRep[dF, yF]];
 AddFermion[Fc, GaugeRep -> T3RGBetaGaugeRep[dF, -yF]]];
secondHeavy = If[yF === 0, F, Fc];
AddYukawa[y1, {Bar @ S1, l, F}, Chirality -> Right];
AddYukawa[y2, {S2, l, secondHeavy}, Chirality -> Right];
T3RGBetaAddMasses[
 dS1_Integer, dS2_Integer, dF_Integer, yF_
] := Module[{fermionInvariant},
 If[yF === 0,
  fermionInvariant = If[dF === 1, (1 &), (del[T3RGBetaRep[dF], #1, #2] &)];
  AddFermionMass[MF, {F, F}, GroupInvariant -> fermionInvariant,
   MassIndices -> ({heavy[#1], heavy[#2]} &), Chirality -> Right],
  fermionInvariant = Which[
   dF === 1, (1 &),
   dF === 2, (eps[T3RGBetaRep[dF], #1, #2] &),
   dF === 3, (del[T3RGBetaRep[dF], #1, #2] &),
   True, Return[$Failed]
  ];
  AddFermionMass[MF, {F, Fc}, GroupInvariant -> fermionInvariant,
   MassIndices -> ({heavy[#1], heavy[#2]} &), Chirality -> Right]
 ]
];
'''
TOP_TEXT='''
 useCConj = which === 1;
 scalarBar = which === 2;
 f = If[dF === 1, If[useCConj, CConj[NewFermion[r]], NewFermion[r]]];
 "Expression" -> PlusHc[y[p, r] NCM[Bar[l[i, p]], f] scalar CG[cg, labels]];
 conjugated = which === 1;
 "Expression" -> PlusHc[coupling[p, r] NCM[Bar[l[i, p]], f] scalar CG[cg, labels]];
 scalar = If[dS === 1, ScalarFieldValue[which, conjugated], ScalarFieldValue[which, conjugated, k]];
 s1 = If[d1 === 1, NewScalar1[], NewScalar1[a]];
 s2 = If[d2 === 1, Bar[NewScalar2[]], Bar[NewScalar2[b]]];
 lambdaT3[] H[i] H[j] s1 s2 CG[T3MixCG, labels];
'''
FIELDS_TEXT='''
 DefineCoupling[y1, Indices -> {Flavor, NFlavor}, SelfConjugate -> False];
 DefineCoupling[y2, Indices -> {Flavor, NFlavor}, SelfConjugate -> False];
 DefineCoupling[lambdaT3, SelfConjugate -> False];
 selfConj = TrueQ[PossibleZeroQ[f["Y"]]] && OddQ[f["SU2"]];
 SelfConjugate -> True;
 SelfConjugate -> False;
'''


def fixture(root:Path,*,old_doublet:bool=False):
    files={
        'RGE/running/rgbeta/T3RGBetaModel.wl':(
            RGB_TEXT.replace('dF === 2, (eps[T3RGBetaRep[dF], #1, #2] &)',
                             'dF === 2, (del[T3RGBetaRep[dF], #1, #2] &)')
            if old_doublet else RGB_TEXT),
        'Lagrangian/interactions/T3Topology.wl':TOP_TEXT,
        'Lagrangian/model/T3Fields.wl':FIELDS_TEXT,
    }
    for name,content in files.items():
        p=root/name;p.parent.mkdir(parents=True,exist_ok=True);p.write_text(content,encoding='utf-8')


class NativeMassAndSpurions(unittest.TestCase):
    def test_all_three_mass_ward_invariants(self):
        for d in (1,2,3):
            self.assertTrue(mass_ward(d)['all_su2_generators_ward_zero'])
        self.assertEqual(mass_ward(2)['tensor'],'epsilon')

    def test_doublet_delta_rejected_by_gauge_generators(self):
        wrong=mass_ward(2,force_delta_doublet=True)
        self.assertFalse(wrong['all_su2_generators_ward_zero'])
        self.assertEqual(wrong['nonzero_generator_numbers'],[1,3])

    def test_phase_invariant_dirac_and_majorana(self):
        for majorana in (False,True):
            info=spurion_witness(majorana)
            self.assertEqual(info['loop_spurion_rephasing_charge'],[0]*4)
            self.assertNotEqual(info['incorrect_nonconjugated_lambda_charge'],[0]*4)
            self.assertFalse(info['absolute_Matchete_RGBeta_coefficient_certified'])
            if not majorana:
                self.assertNotEqual(info['incorrect_same_chirality_dirac_assignment_charge'],[0]*4)

    def test_independent_extra_yukawa_operator_phases(self):
        for majorana in (False,True):
            a=spurion_witness(majorana,'y1Conj')
            b=spurion_witness(majorana,'y2Conj')
            self.assertNotEqual(a['additional_source_allowed_independent_yukawa']['operator_charge'],
                                (MAJORANA if majorana else DIRAC)['y1'])
            self.assertNotEqual(b['additional_source_allowed_independent_yukawa']['operator_charge'],
                                (MAJORANA if majorana else DIRAC)['y2'])

    def test_wrong_spurion_conjugation_negative_control(self):
        for majorana in (False,True):
            self.assertEqual(phase_sum({'y1':1,'y2':1,'MF':-1,'lambdaT3':-1},majorana),(0,0,0,0))
            self.assertNotEqual(phase_sum({'y1':1,'y2':1,'MF':1,'lambdaT3':-1},majorana),(0,0,0,0))

    def test_strip_nested_comments_cannot_hide_source_drift(self):
        self.assertEqual(without_wl_comments('x (* one (* nested *) end *) y'),'x  y')
        with self.assertRaises(DictionaryError):without_wl_comments('x (* unterminated')

    def test_frozen_source_contract_all_16(self):
        with tempfile.TemporaryDirectory() as t:
            root=Path(t);fixture(root)
            r=audit(root)
            self.assertEqual(r['model_count'],16)
            self.assertEqual(r['model_counts'],{'Dirac':14,'Majorana':2,'extra_yukawa_allowed':6})
            self.assertEqual(r['physical_source_certified_models'],0)
            self.assertEqual({x['model'] for x in r['models']},
                {'A_-4','A_-2','A_+0','B_-3','B_-1','B_+1','C_-3','C_-1','C_+1',
                 'D_-2','D_+0','D_+2','E_-4','E_-2','E_+0','E_+2'})

    def test_native_old_doublet_delta_blocks(self):
        with tempfile.TemporaryDirectory() as t:
            root=Path(t);fixture(root,old_doublet=True)
            with self.assertRaisesRegex(DictionaryError,'RGBeta fermion mass source contract changed'):
                audit(root)

    def test_source_y2_or_mf_drift_blocks(self):
        for filename,old,new in [
            ('Lagrangian/interactions/T3Topology.wl','scalarBar = which === 2','scalarBar = which === 1'),
            ('RGE/running/rgbeta/T3RGBetaModel.wl','{F, Fc}', '{Fc, Fc}'),
        ]:
            with tempfile.TemporaryDirectory() as t:
                root=Path(t);fixture(root)
                p=root/filename;p.write_text(p.read_text().replace(old,new))
                with self.assertRaises(DictionaryError):guard_sources(root)

    def test_stale_native_A0_report_is_not_a_certificate(self):
        with tempfile.TemporaryDirectory() as t:
            root=Path(t);fixture(root)
            p=root/'Reports/output/workstream_a/a0_dirac_mass_ward/a0_dirac_mass_ward.json'
            p.parent.mkdir(parents=True)
            p.write_text(json.dumps({'status':'NativeA0DiracMassSU2WardViolationConfirmed',
                                     'registered_MF_invariant':'del[SU2L[fund], dm1, dm2]',
                                     'actual_su2_mass_ward_passed':False,
                                     'local_rgbeta_source_sha256':'old_sha'}))
            report=optional_evidence(root)
            self.assertFalse(report['a0_native_mass_report']['current_source_mass_ward_certificate'])
            self.assertTrue(report['a0_native_mass_report']['stale_or_incompatible_native_report'])
            self.assertEqual(audit(root)['physical_source_certified_models'],0)

    def test_full_cli_honours_output_guard(self):
        with tempfile.TemporaryDirectory() as t:
            root=Path(t);fixture(root)
            script=Path(__file__).resolve().parents[1]/'scripts/RunLODirectUVSpinorMassDictionary.py'
            out=root/'output/general_observable/lo_certification_direct_uv/test_01'
            cmd=[sys.executable,str(script),'--project-root',str(root),'--output-dir',str(out)]
            proc=subprocess.run(cmd,capture_output=True,text=True)
            self.assertEqual(proc.returncode,0,proc.stderr)
            result=json.loads((out/'all_16_spinor_mass_dictionary.json').read_text())
            self.assertEqual(result['model_count'],16)
            self.assertEqual(subprocess.run(cmd,capture_output=True,text=True).returncode,2)
            self.assertEqual(subprocess.run([sys.executable,str(script),'--project-root',str(root),
                '--output-dir',str(root/'bad')],capture_output=True,text=True).returncode,2)

if __name__=='__main__':unittest.main()
