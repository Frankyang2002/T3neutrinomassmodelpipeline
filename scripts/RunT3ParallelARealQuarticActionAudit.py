"""Conditional real-scalar RGBeta quartic-action algebra check (50 extra slots)."""
from __future__ import annotations
import argparse,json,sys
from pathlib import Path
from Numerical.parallel_t3.wilson.scalar_action_audit import audit_from_paths
from Numerical.parallel_t3.wilson.matchete_rgbeta_dictionary import DictionaryEvidenceError

def main(argv=None):
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--matchete',type=Path,required=True)
    p.add_argument('--rgbeta',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True)
    p.add_argument('--samples',type=int,default=3)
    p.add_argument('--seed',type=int,default=3241)
    a=p.parse_args(argv)
    if a.output.exists():p.error(f'Refusing to overwrite {a.output}')
    try:
        result=audit_from_paths(a.matchete,a.rgbeta,samples=a.samples,seed=a.seed)
    except (ValueError,DictionaryEvidenceError,KeyError,OSError,TypeError) as exc:
        print(f'Blocked_RealQuarticActionAudit: {exc}',file=sys.stderr);return 2
    a.output.parent.mkdir(parents=True,exist_ok=True)
    a.output.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n',encoding='utf-8')
    print(f"{result['extra_quartic_slots']}/50 complex extra-quartic actions; residual={result['max_real_action_relative_residual']:.3e}")
    print(f"{result['baseline_norm_audit']['baseline_norm_coupling_slots']}/96 basic norm-quartic model-coupling actions; residual={result['baseline_norm_audit']['max_baseline_norm_action_residual']:.3e}")
    print(result['status'], '->',a.output)
    print('NOT full scalar normalization; native RGBeta coordinates and baseline portals remain uncertified.')
    return 0
if __name__=='__main__':raise SystemExit(main())
