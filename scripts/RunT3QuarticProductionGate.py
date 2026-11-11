"""Read-only readiness check for production full-quartic EFT1 tensor responses."""
from __future__ import annotations
import argparse
import json
from pathlib import Path
from Numerical.running.T3QuarticProductionGate import audit_gate
from Numerical.running.T3RunningQuarticTensor import QuarticEvidenceError

ROOT = Path(__file__).resolve().parents[1]
DEFAULT = ROOT / 'output/general_t3_stage22'
REGISTRY = ROOT / 'RGE/general/generated/general_rgbeta_quartic_registry.json'


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--model', required=True, help='e.g. A_-4 or B_-1')
    p.add_argument('--baseline-symbols', nargs='+', required=True,
                   help='Exact baseline scalar quartic coupling names from EFT1 state')
    p.add_argument('--registry', type=Path, default=REGISTRY)
    p.add_argument('--stage22', type=Path, default=DEFAULT/'stage22_summary.json')
    p.add_argument('--response', type=Path, required=True)
    p.add_argument('--binding', type=Path, required=True)
    p.add_argument('--output', type=Path, default=None)
    args = p.parse_args(argv)
    try:
        result = audit_gate(args.registry,args.stage22,args.response,args.binding,
                            args.model,set(args.baseline_symbols))
        result.pop('basis')
        code=0
    except (QuarticEvidenceError,KeyError,ValueError,OSError) as exc:
        result = {'model':args.model,'status':'BLOCKED', 'reason':str(exc)}
        code=2
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(result,indent=2))
    return code


if __name__ == '__main__':
    raise SystemExit(main())
