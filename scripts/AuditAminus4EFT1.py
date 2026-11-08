"""Audit the *existing* A(-4) seven-scale EFT1 tensor export without rewriting it."""
from __future__ import annotations
import argparse
import json
import sys
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:sys.path.insert(0,str(ROOT))
from Numerical.running.Aminus4EFT1Audit import audit_export


def main(argv=None):
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input-dir',type=Path,
        default=ROOT/'output/a_minus4_eft1_batch2')
    parser.add_argument('--beta-file',type=Path,
        default=ROOT/'output/full16_symbolic/T3_A_alpha_m4/data/eft1_rgbeta_rge.json')
    parser.add_argument('--output',type=Path,
        default=ROOT/'output/a_minus4_eft1_batch2/A_m4_eft1_audit.json')
    parser.add_argument('--no-ode-check',action='store_true',
        help='Skip the independent re-integration of each saved RGE interval')
    args=parser.parse_args(argv)
    report=audit_export(args.input_dir,args.beta_file,root=ROOT,
        check_sources=True,verify_ode=not args.no_ode_check)
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8')
    print('A(-4) tensor/gauge reconstruction:',
          'PASS' if report['tensor_reconstruction_and_gauge_ward_passed'] else 'FAIL')
    print('Saved scales verified:',report['sample_count'])
    print('Independent ODE interval checks:',len(report['ode_interval_reintegration']))
    print('Native/reference scalar beta agreement:',
          'YES' if report['all_native_reference_scalar_beta_coefficients_agree'] else 'NO')
    for item in report['scalar_beta_diagnostics']:
        print(f"  mu={item['mu_gev']:g}: {', '.join(item['mismatched_coupling_betas']) or '(none)'}")
    print('Physical Matchete coupling map certified: NO')
    print('Report:',args.output)
    return 0

if __name__=='__main__':raise SystemExit(main())
