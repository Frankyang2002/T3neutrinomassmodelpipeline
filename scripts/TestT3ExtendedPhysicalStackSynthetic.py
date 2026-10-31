#!/usr/bin/env python3
"""Synthetic end-to-end orchestrator controls. No Wolfram matching is simulated as physical evidence."""
from __future__ import annotations

import json
from pathlib import Path
import sys
from tempfile import TemporaryDirectory
from types import SimpleNamespace
from unittest.mock import patch

ROOT=Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:sys.path.insert(0,str(ROOT))
from scripts import RunT3ExtendedMatchingBatch as batch
from Numerical.experimental import T3ExtendedC5KernelProxy as proxy
from Numerical.experimental import T3ExtendedC5Observables as physical


def fixture(root:Path):
    # Only minimal source guards are constructed; does NOT claim GitHub sources.
    sources={
       'Lagrangian/model/LagrangianBuilder.wl':'BuildT3Lagrangian[model_Association]',
       'Lagrangian/interactions/T3Topology.wl':(
           'T3GeneralExtraYukawaCount[model_Association]\nBuildT3ExtraYukawaCandidates'),
       'Lagrangian/RunMatching.wl':'RunT3Matching',
       'Lagrangian/RunT3ExtendedC5Native.wl':(
           'RunT3Matching[built["LUV"],5,1]\nExtractWeinbergCoefficient\n'
           'T3_EXPERIMENTAL_GENERAL_QUADRATICS\n'
           'native_coefficient_extra_difference_path'),
       'RGE/running/rgbeta/T3RGBetaModel.wl':'T3RGBetaBuild',
    }
    for path,text in sources.items():
        p=root/path;p.parent.mkdir(parents=True,exist_ok=True)
        p.write_text(text,encoding='utf-8')
    summary=root/'Reports/output/workstream_a/full_coupled_general_yukawa_uv/summary.json'
    batch.atomic(summary,{'status':'AllRequestedFullCoupledUVRGBetaOnlyPassed',
       'rows':[{'model':name,'status':'CoupledGeneralYukawaRGBetaOnlyUVPassed'}
               for name in (*batch.SPECS, *(f'control_{i}' for i in range(10)))]})
    batch.atomic(root/'Reports/output/workstream_a/integrated_general_yukawa_uv/summary.json',
                 {'all_requested_models_succeeded':True})
    for model,(_,_,extra) in batch.SPECS.items():
        batch.atomic(root/'Reports/output/workstream_a/full_coupled_general_yukawa_uv/models'/f'{model}.json',
            {'status':'CoupledGeneralYukawaRGBetaOnlyUVPassed','extra_yukawa':extra,
             'endpoint_above_evolved_heavy_mass':True})


def fake_native_run(cmd,cwd,env,stdout,stderr,text,errors,timeout,check):
    assert cmd[1]=='-file' and cmd[4] in ('A','D','E') and cmd[5] in ('p0','m2')
    assert env['T3_EXPERIMENTAL_GENERAL_FERMIONS']=='1'
    assert env['T3_EXPERIMENTAL_GENERAL_QUADRATICS']=='0'
    out=Path(cmd[3]);out.parent.mkdir(parents=True,exist_ok=True)
    model=cmd[4]+'_'+cmd[5]
    data={'format':'t3-extended-c5-native-matchete-one-step-v1',
          'model':model,'heavy_flavours':int(cmd[6]),
          'status':batch.SUCCESS,'additional_yukawa':batch.SPECS[model][2],
          'one_loop_matching_performed':True,
          'native_matchete_lagrangian_check_passed':True,
          'matched_c5_in_matchete_native_coordinates':True,
          'absolute_matchete_rgbeta_dictionary_certified':False,
          'hierarchical_matching_certified':False,
          'cross_package_numeric_substitution_authorized':False,
          'physical_neutrino_mass_prediction':False}
    for key in ('native_coefficient_full_path','native_coefficient_zero_extra_path',
                'native_coefficient_extra_difference_path'):
        file=out.parent/(out.stem+key+'.txt');file.write_text('c5 + yConj',encoding='utf-8')
        data[key]=str(file)
    batch.atomic(out,data)
    return SimpleNamespace(stdout='synthetic native result\n',returncode=0)


def main():
    batch.self_test();proxy.self_test();physical.self_test()
    with TemporaryDirectory() as td:
        r=Path(td)
        fixture(r)
        assert len(batch.fingerprints(r))==len(batch.SOURCES)
        assert batch.previous_runs(r)['status']=='Previous16RGESucceededSixNativeTargetsAvailable'
        with patch.object(batch.subprocess,'run',side_effect=fake_native_run):
            for model in batch.SPECS:
                rec=batch.run_one(model,heavy_flavours=3,wl='fake-wolfram',timeout=5,
                                  root=r,output=r/'native_outputs')
                assert rec['status']==batch.SUCCESS and rec['exit_code']==0,rec
        with patch.object(batch.subprocess,'run',return_value=SimpleNamespace(
                stdout='NATIVE_EXT_C5_STAGE: match\nMatchete native error\n',returncode=6)):
            failure=batch.run_one('A_p0',heavy_flavours=3,wl='fake-wolfram',timeout=5,
                                  root=r,output=r/'native_outputs')
            assert failure['status']=='Blocked' and 'Matchete native error' in failure['log_tail']
        # Negative provenance control: a single missing model must block.
        missing=r/'Reports/output/workstream_a/full_coupled_general_yukawa_uv/models/A_p0.json'
        missing.unlink()
        try:batch.previous_runs(r)
        except FileNotFoundError:pass
        else:raise AssertionError('Missing UV evidence passed')
    print('PASS: six synthetic native export validations, per-model log failure, prior evidence gate')
    print('SYNTHETIC: actual Matchete matching and physical coupling dictionary NOT demonstrated')


if __name__=='__main__':main()
