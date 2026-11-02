"""Independent T3 all-16 direct UV source certification/benchmark launcher.

PowerShell: python scripts/RunLODirectUVCertification.py --project-root .
This never modifies source, frozen benchmarks, EFT results or old UV reports.
"""
from __future__ import annotations
import argparse
from pathlib import Path
import sys

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
from Numerical.lo_certification.direct_uv.runner import run


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--project-root',type=Path,default=ROOT)
    p.add_argument('--input-dir',type=Path,default=None,
                   help='Existing source-derived component inputs; defaults to previous provisional_inputs')
    p.add_argument('--output-dir',type=Path,default=None,
                   help='Must be below owned lo_certification_direct_uv directory')
    p.add_argument('--eft-reference-dir',type=Path,default=None,
                   help='Strict same-input same-scale physical mass matrix JSONs')
    p.add_argument('--suppressed-fraction',type=float,default=.01)
    p.add_argument('--convention-scan',action='store_true',
                   help='Native portal geometry, exact field-rephasing, conditional signed-coupling mass sensitivity')
    a=p.parse_args()
    if not (0<a.suppressed_fraction<1):p.error('--suppressed-fraction must be between 0 and 1')
    out=a.output_dir or a.project_root/'output/general_observable/lo_certification_direct_uv'
    summary=run(a.project_root,out,a.input_dir,a.eft_reference_dir,a.suppressed_fraction,
                convention_scan_enabled=a.convention_scan)
    if a.convention_scan:
        print('Conditional portal/phase scans: '
              f"{summary['counts']['convention_scan_completed']}/16; "
              f"scan blocked: {summary['counts']['convention_scan_blocked']}/16")
    print('Exact provisional: {exact_uv_provisional}/16; physical: '
          '{physical_source_certified}/16; blocked: {blocked}/16'.format(**summary['counts']))
    valid_base=summary['counts']['exact_uv_provisional']==16
    valid_scan=(not a.convention_scan or summary['counts']['convention_scan_completed']==16)
    return 0 if valid_base and valid_scan else 2

if __name__=='__main__':raise SystemExit(main())
