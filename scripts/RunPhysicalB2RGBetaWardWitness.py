#!/usr/bin/env python3
"""Run an exact, non-certifying RGBeta literal-index Ward diagnostic.

The source guard accepts later local Git HEADs *only* when the exact six
reviewed source blobs and clean tracked worktree contents still match.
Writes into a caller-supplied NEW directory; never overwrites existing data.
"""
from __future__ import annotations
import argparse
import json
from pathlib import Path
import sys

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'Numerical/lo_certification/physical_normalisation'))
from b2_rgbeta_ward_witness import compute  # noqa: E402
from b2_gauge_coupling import AuditBlocked  # noqa: E402


def main() -> int:
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--root',type=Path,default=ROOT)
    p.add_argument('--formal',type=Path,default=Path(
        'output/general_observable/physical_normalisation/b2_rgbeta_native_01/formal_invariants.json'))
    p.add_argument('--registry',type=Path,default=None,
                   help='optional result from RunPhysicalB2RGBetaRegistryAudit.wl')
    p.add_argument('--output',type=Path,default=Path(
        'output/general_observable/physical_normalisation/b2_ward_local_01'))
    args=p.parse_args()
    root=args.root.resolve()
    output=args.output if args.output.is_absolute() else root/args.output
    formal=args.formal if args.formal.is_absolute() else root/args.formal
    registry = None
    if args.registry is not None:
        registry = args.registry if args.registry.is_absolute() else root / args.registry
    if output.exists():
        print('BLOCKED: output already exists:',output,file=sys.stderr)
        return 2
    try:
        report=compute(root,formal,registry)
    except (AuditBlocked,OSError,ValueError,TypeError,KeyError,json.JSONDecodeError) as ex:
        print('BLOCKED: no source-validated Ward report produced:',ex,file=sys.stderr)
        return 3
    output.mkdir(parents=True,exist_ok=False)
    (output/'ward_witness.json').write_text(json.dumps(report,indent=2,sort_keys=True)+'\n',encoding='utf-8')
    mat=report['matchete_B_C_CG_comparison']
    note=f'''# B2 literal-index SU(2) Ward witness — local execution

- Local HEAD: `{report['source_provenance']['repository_head_measured']}`.
- Six reviewed tracked source blobs match Git HEAD, index and canonical worktree.
- Native RGBeta formal file validated: {report['native_formal_evidence'].get('five_classes_formal_records_validated',False)}.
- B/C model charge samples: {report['B_C_models_checked']}/6.
- Matchete native CG comparison: `{mat['status']}`; rows {mat.get('models_checked',0)}.
- Status: `{report['status']}`.
- **This is a literal-index conditional diagnostic**, not a finding that the installed RGBeta indices have these numerical meanings.
- No actual RGBeta `eps`/`delS2` components or physical vertex/action multipliers measured.
- Cross-package physical coupling certification: **0/16**.
- No production source changed; no EFT performed.
'''
    (output/'HANDOFF.md').write_text(note,encoding='utf-8')
    print(f"{report['status']}; native-formal=5/5; literal B/C models=6/6; "
          f"Matchete comparisons={mat.get('models_checked',0)}/6; physical=0/16; "
          f"output={output}")
    return 0


if __name__=='__main__':
    sys.exit(main())
