"""Run all-16 V6 live/archival Matchete extra-quartic source audit."""
from __future__ import annotations
import argparse
import json
from pathlib import Path
from Numerical.lo_certification.couplings.extra_quartic_batch_v6 import run_v6_batch


def main()->int:
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--repo-root',default='.')
    p.add_argument('--run-name',default='ws2_extra_quartic_source_v1')
    p.add_argument('--previous-v5',default='ws2_matchete_rgbeta_dictionary_v1')
    p.add_argument('--v4-run-name',default='ws2_named_hhss_v1')
    p.add_argument('--run-live',action='store_true',help='Run live Matchete GetCGTensor for all models')
    p.add_argument('--wolframscript',default='wolframscript')
    p.add_argument('--live-export-root',type=Path)
    p.add_argument('--probe-orientations',action='store_true',help='Slow conditional exact CG tensor probes; requires live export')
    p.add_argument('--wolfram-timeout-seconds',type=int,default=900)
    p.add_argument('--model',action='append')
    a=p.parse_args()
    if not a.run_name or a.run_name in ('.','..') or '/' in a.run_name or '\\' in a.run_name:
        p.error('--run-name must be one bare directory name')
    root=Path(a.repo_root).resolve()
    out=root/'output/general_observable/lo_certification_couplings'/a.run_name
    result=run_v6_batch(root,out,models=tuple(a.model) if a.model else None,
       previous_v5=a.previous_v5,v4_run_name=a.v4_run_name,
       live=a.run_live,wolframscript=a.wolframscript,
       live_export_root=a.live_export_root,probe_orientations=a.probe_orientations,
       wolfram_timeout_seconds=a.wolfram_timeout_seconds)
    print(json.dumps({k:result[k] for k in ('models_processed',
       'all16_source_inventory_passed','all16_live_Matchete_extra_operator_exports_passed',
       'additional_quartic_operator_inventory',
       'additional_HHSS_operator_inventory',
       'live_native_extra_Matchete_operators_verified',
       'live_native_extra_Matchete_CGs_verified',
       'native_named_extra_HHSS_in_EFT1_seeds',
       'physically_certified_HHSS_coupling_maps','physical_C5_calculated')},indent=2))
    complete=result['all16_source_inventory_passed']
    if a.run_live or a.live_export_root is not None:
        complete=complete and result['all16_live_Matchete_extra_operator_exports_passed']
    return 0 if complete else 1

if __name__=='__main__':raise SystemExit(main())
