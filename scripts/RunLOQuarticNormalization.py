"""All-16 source-audited RGBeta/Matchete portal normalization and HHSS inventory.

Examples:
  python -m scripts.RunLOQuarticNormalization --repo-root . --run-name ws2_portals_v1
  python -m unittest discover -s tests -p 'test_lo_couplings_normalization*.py' -v
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path


def main(argv=None)->int:
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--repo-root',type=Path,default=Path('.'))
    p.add_argument('--run-name',default='ws2_portal_normalization_v1')
    p.add_argument('--model',action='append',default=None,
                   help='Diagnostic subset; omit for all 16 models')
    args=p.parse_args(argv)
    name=args.run_name
    if not name or name in ('.','..') or '/' in name or '\\' in name:
        p.error('--run-name must be a fresh single folder name')
    root=args.repo_root.resolve()
    from Numerical.lo_certification.couplings.quartic_audit import run_portal_batch
    output=root/'output/general_observable/lo_certification_couplings'/name
    doc=run_portal_batch(root,output,tuple(args.model) if args.model else None)
    print(json.dumps({'output_directory':str(output),
                      'models_processed':doc['models_processed'],
                      'all16_audit_passed':doc['all16_audit_passed'],
                      'verified_full_rgbeta_portal_span_models':doc['verified_full_rgbeta_portal_span_models'],
                      'newly_discovered_baseline_portal_deficiency':doc['newly_discovered_baseline_portal_deficiency'],
                      'extra_HHSS_registry_entries':doc['extra_HHSS_registry_entries'],
                      'physical_C5_calculated':False},indent=2))
    return 0 if doc['all16_audit_passed'] else 2

if __name__=='__main__':
    raise SystemExit(main())
