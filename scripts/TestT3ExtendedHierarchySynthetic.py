#!/usr/bin/env python3
"""Synthetic-only source, contract and orchestration tests for the native hierarchy pilot."""
from __future__ import annotations

import importlib.util
import json
from pathlib import Path
from types import SimpleNamespace
import tempfile
from unittest.mock import patch

HERE = Path(__file__).resolve().parent
spec = importlib.util.spec_from_file_location('t3_hierarchy_runner', HERE/'RunT3ExtendedHierarchyBatch.py')
assert spec and spec.loader
mod = importlib.util.module_from_spec(spec)
spec.loader.exec_module(mod)


def expect_block(fn, label):
    try:
        fn()
    except mod.GateError:
        return
    raise AssertionError('Bad witness accepted: '+label)


def make_native(model:str, mode:str, variant:str)->dict:
    entry = {'format':'t3-extra-yukawa-hierarchical-native-v1',
       'model':model,'mode':mode,'variant':variant,'heavy_flavours':3,
       'status':mod.SUCCESS[mode],
       'fermion_threshold_integrated_out':True,
       'scalar_threshold_integrated_out':mode=='hierarchy',
       'stage_count':1 if mode=='stage1' else 2,
       'eft1_Wilson_RGE_evolved':False,
       'absolute_matchete_rgbeta_dictionary_certified':False,
       'physical_low_energy_prediction_certified':False,
       'full_scalar_quadratic_basis_included':False}
    if mode=='hierarchy':
        entry.update({'stage2_c5_is_finite_renormalized':False,
                      'EFT1_Wilson_running_insertion_added':False})
    return entry


def make_diff(model:str,mode:str)->dict:
    req=['stage1_tree','stage1_one_loop','stage1_full']
    if mode=='hierarchy':req.extend(['stage2_full','stage2_c5_ordered'])
    return {'status':'NativeEFTDifferenceExported','model':model,'mode':mode,
            'seeded_vs_uv_zero_expressions_compared':req,
            'difference_not_proven_zero_by_PossibleZeroQ':{x:True for x in req},
            'absolute_coupling_dictionary_certified':False,
            'full_renormalized_hierarchical_C5_certified':False,
            'EFT1_running_inserted':False,
            'cross_kernel_dummy_index_canonicalisation_certified':False}


def test_fixture_run():
    with tempfile.TemporaryDirectory(prefix='t3_hierarchy_mock_') as tmp:
        out=Path(tmp)
        calls=[]
        def mock_invoke(command,log,timeout,env):
            log.parent.mkdir(parents=True,exist_ok=True)
            log.write_text('MOCK ONLY: never a live Matchete result\n',encoding='utf-8')
            if 'CompareT3ExtendedHierarchy.wl' in command[2]:
                *_, model = log.parent.parts
                mode=command[-1]
                (log.parent/'difference.json').write_text(json.dumps(make_diff(model,mode)),encoding='utf-8')
            else:
                dest=Path(command[3]);model=command[4]+'_'+command[5]
                variant,mode=command[6],command[7]
                (dest/'native.json').write_text(json.dumps(make_native(model,mode,variant)),encoding='utf-8')
            calls.append(command)
            return 0
        args=SimpleNamespace(wolframscript='wolframscript',timeout_seconds=3,use_existing=False)
        with patch.object(mod,'invoke',side_effect=mock_invoke):
            for mode in ('stage1','hierarchy'):
                for model in mod.MODELS:
                    row=mod.run_model(model,mode,out,args)
                    assert row['status'].startswith('NATIVE_')
                    assert row['comparison']['difference_not_proven_zero_by_PossibleZeroQ']['stage1_full'] is True
        assert len(calls)==6*2*3,'Each model/mode must have two independent native runs and a comparator'
        # A blocked native status must prevent the comparison stage.
        a=out/'stage1/A_p0/seeded/native.json'
        damaged=json.loads(a.read_text(encoding='utf-8'))
        damaged['physical_low_energy_prediction_certified']=True
        a.write_text(json.dumps(damaged),encoding='utf-8')
        args.use_existing=True
        row=mod.run_model('A_p0','stage1',out,args)
        assert row['status']=='Blocked'
        assert 'physical certification' in row['variants']['seeded']['reason']
    print('PASS: 6x2 synthetic native seeded/zero matching and expression comparison; tampered native report blocks')


def test_wolfram_sources():
    a=mod.NATIVE.read_text(encoding='utf-8')
    b=mod.COMPARE.read_text(encoding='utf-8')
    assert a.index('Needs["Matchete`"]') < a.index('RunSequentialT3Matching[')
    assert 'SetT3HeavyFields[{"F"}]' in a
    assert '{{"F"},{"S1","S2"}}' in a
    assert 'HoldPattern[Bar[Coupling[y1Conj,___]]]' in a
    assert 'HoldPattern[Bar[Coupling[y2Conj,___]]]' in a
    assert '"stage1_one_loop"' in a and '"stage2_c5_ordered"' in a
    assert '"stage2_c5_is_finite_renormalized"->False' in a
    assert '"EFT1_running_inserted"->False' in b
    assert 'Expand[s-z]' in b
    # Repeated bracket balance check, ignoring Mathematica comments/strings.
    for path in (mod.NATIVE,mod.COMPARE):
        src=path.read_text(encoding='utf-8')
        i=0;filtered=[];depth=0;string=False
        while i<len(src):
            if not string and src.startswith('(*',i):depth+=1;i+=2;continue
            if depth and src.startswith('*)',i):depth-=1;i+=2;continue
            if depth:i+=1;continue
            if src[i]=='"':
                string=not string;i+=1;continue
            if string and src[i]=='\\':i+=2;continue
            if not string:filtered.append(src[i])
            i+=1
        assert depth==0 and not string, 'Unclosed string/comment in '+str(path)
        stack=[]
        for ch in filtered:
            if ch in '([{':stack.append(ch)
            elif ch in ')]}':
                assert stack and {'(':')','[':']','{':'}'}[stack.pop()]==ch, 'Wolfram bracket error in '+str(path)
        assert not stack,'Unclosed Wolfram bracket in '+str(path)
    print('PASS: Wolfram full sources, separate load, UV-zero substitution, fixed-order provenance, balanced tokens')


def main():
    mod.self_test()
    for mode in ('stage1','hierarchy'):
        for model in mod.MODELS:
            for var in ('seeded','zero'):
                mod.check_native(make_native(model,mode,var),model,mode,var)
            mod.check_difference(make_diff(model,mode),model,mode)
    for bad in ({'EFT1_running_inserted':True},
                {'absolute_coupling_dictionary_certified':True}):
        expect_block(lambda:mod.check_difference({**make_diff('A_p0','hierarchy'),**bad},'A_p0','hierarchy'),str(bad))
    test_wolfram_sources()
    test_fixture_run()
    print('SYNTHETIC ONLY: actual Matchete sequential thresholds and EFT1 Wilson running not demonstrated')


if __name__=='__main__':main()
