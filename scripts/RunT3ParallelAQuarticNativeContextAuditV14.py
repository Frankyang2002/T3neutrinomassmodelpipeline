"""Read-only native RGBeta PackageScope quartic audit (v14 with complex v15 probes).

--inspect reports blocked or passed evidence without promoting a failure;
--audit requires all 16 models and fresh, source-rechecked v9 portal evidence.
Never modifies an existing JSON.
"""
from __future__ import annotations
import argparse
import json
from pathlib import Path
import sys
from Numerical.parallel_t3.wilson.native_context_quartic_v14 import (
    V14Blocked, audit_witness, audit_from_paths, file_hash, read,
)


def write_immutable(out: Path, report: dict) -> None:
    encoded = json.dumps(report, indent=2, sort_keys=True, allow_nan=False) + '\n'
    if out.exists():
        if out.read_text(encoding='utf-8') != encoded:
            raise V14Blocked('refusing to overwrite different existing report: '+str(out))
        print('v14 existing report matches; preserved:',out)
        return
    out.parent.mkdir(parents=True,exist_ok=True)
    tmp=out.with_suffix(out.suffix+'.partial')
    if tmp.exists():raise V14Blocked('partial output already exists: '+str(tmp))
    tmp.write_text(encoded,encoding='utf-8')
    tmp.replace(out)


def main(argv=None) -> int:
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--project-root',type=Path,default=Path('.'))
    ap.add_argument('--native',type=Path,required=True)
    group=ap.add_mutually_exclusive_group(required=True)
    group.add_argument('--inspect',action='store_true')
    group.add_argument('--audit',action='store_true')
    ap.add_argument('--v9',type=Path)
    ap.add_argument('--output',type=Path)
    a=ap.parse_args(argv)
    try:
        root=a.project_root.resolve()
        native=a.native if a.native.is_absolute() else root/a.native
        if a.inspect:
            report=audit_witness(read(native),root,allow_blocked=True)
            print(f"v14: {report['status']}"
                  f"; real diagonals {report['source_census'].get('real_diagonal_unit',0)}/"
                  f"{report['source_census']['real_diagonal_attempted']}; portal responses "
                  f"{report['source_census'].get('portal_response_numeric',0)}/"
                  f"{report['source_census']['portal_response_attempted']}; complex direct+conj "
                  f"{report['complex_probe_census'].get('complex_direct_unit_conj_zero',0)}/"
                  f"{report['complex_probe_census']['complex_attempted']}")
            if report['failures']:
                print(f'v14: {len(report["failures"])} diagnostic records (first 10):')
                for row in report['failures'][:10]:
                    r=row.get('native_response')
                    print('  ',row['model'],row['type'],row.get('symbol',row.get('target','')),
                          r.get('reason','non-unit response') if isinstance(r,dict) else '')
            if a.output:
                out=a.output if a.output.is_absolute() else root/a.output
                write_immutable(out,report)
            if not report['real_projector_conditionally_verified']:
                print('v14 BLOCKED. Inspection succeeded, not projector certification.')
            return 0
        require_file=a.v9 if a.v9 else Path('output/parallel_t3/wilson/portal_real_action_conditional_v9.json')
        v9=require_file if require_file.is_absolute() else root/require_file
        report=audit_from_paths(root,native,v9)
        out=a.output if a.output else Path('output/parallel_t3/wilson/native_quartic_projectors_conditional_v14.json')
        out=out if out.is_absolute() else root/out
        write_immutable(out,report)
        print('v14: all 16 native real-projector sectors passed; 150/150 diagonals, 110/110 portal response entries')
        print('v14: complex source/conjugate exploratory channels',report['complex_probe_census'])
        print(report['status'],'->',out)
        print('No numeric S2 frame / complete physical potential / Wilson evolution authorized.')
        return 0
    except (V14Blocked,ValueError,TypeError,KeyError,IndexError,ArithmeticError,OSError) as exc:
        print('v14 BLOCKED:',exc,file=sys.stderr)
        return 7

if __name__=='__main__':raise SystemExit(main())
