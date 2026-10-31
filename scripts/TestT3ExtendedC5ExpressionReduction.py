#!/usr/bin/env python3
"""Test real uploaded native expressions plus tampered-input negative controls."""
import argparse
import json
import sys
import tempfile
from pathlib import Path
from zipfile import ZipFile,ZIP_DEFLATED

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
from Numerical.experimental.T3ExtendedC5ExpressionReduction import MODELS,SECTORS,C5SourceError,verify_model,self_test


def must_reject(fn, name):
    try: fn()
    except (C5SourceError, ValueError, KeyError): return
    raise AssertionError('Tampered source accepted: '+name)


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--zip',type=Path,help='Native expressions ZIP, otherwise use local matching output')
    args=ap.parse_args()
    self_test()
    original=args.zip or ROOT/'Reports/output/workstream_a/extended_c5_native/native'
    actual={m:verify_model(original,m) for m in MODELS}
    assert len(actual)==6
    assert {v['common_group_prefactor_exact'] for v in actual.values()}=={'-2*sqrt(3)/3','2*sqrt(6)/3'}
    with tempfile.TemporaryDirectory() as temp:
        out=Path(temp)
        name='A_p0_h3.extra_difference_C5.inputform.txt'
        for k in ('sign','quartic','index','loop'):
            # Copy only A_p0 to avoid growing the temporary ZIP.
            entries={}
            for f in ['A_p0_h3.json']+[f'A_p0_h3.{s}.inputform.txt' for s in SECTORS]:
                if original.suffix.lower()=='.zip':
                    with ZipFile(original) as z:entries[f]=z.read(f)
                else:entries[f]=(original/f).read_bytes()
            old=entries[name].decode('utf-8')
            if k=='sign':new=old.replace('(-2*hbar*','(2*hbar*',1)
            elif k=='quartic':new=old.replace('lambdaT3ExtraQ021010I1','lambdaT3ExtraQ020101I1',1)
            elif k=='index':new=old.replace('d$$2, Flavor','d$$1, Flavor',1)
            else:new=old.replace('Coupling[MF, {}, 0]^2, Coupling[MS1, {}, 0]^2, Coupling[MS2, {}, 0]^2]',
                                 'Coupling[MS2, {}, 0]^2, Coupling[MS1, {}, 0]^2, Coupling[MF, {}, 0]^2]',1)
            assert old != new
            entries[name]=new.encode('utf-8')
            file=out/(k+'.zip')
            with ZipFile(file,'w',compression=ZIP_DEFLATED) as z:
                for fname,v in entries.items():z.writestr(fname,v)
            must_reject(lambda:verify_model(file,'A_p0'),k)
        entries={}
        for f in ['A_p0_h3.json']+[f'A_p0_h3.{s}.inputform.txt' for s in SECTORS]:
            if original.suffix.lower()=='.zip':
                with ZipFile(original) as z:entries[f]=z.read(f)
            else:entries[f]=(original/f).read_bytes()
        meta=json.loads(entries['A_p0_h3.json'])
        meta['cross_package_numeric_substitution_authorized']=True
        entries['A_p0_h3.json']=json.dumps(meta).encode()
        f=out/'wrong_meta.zip'
        with ZipFile(f,'w',compression=ZIP_DEFLATED) as z:
            for n,v in entries.items():z.writestr(n,v)
        must_reject(lambda:verify_model(f,'A_p0'),'wrong physical certification')
    print('PASS: real native 6/6; all three exact sectors; tamper checks (sign, quartic, flavour index, loop, metadata)')

if __name__=='__main__':main()
