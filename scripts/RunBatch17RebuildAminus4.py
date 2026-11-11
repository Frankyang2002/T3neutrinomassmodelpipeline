"""Rebuild A(-4) matching tensors from new Batch14 muS endpoint, then run Batch17.

Project path: scripts/RunBatch17RebuildAminus4.py

All source files and previous runs are read-only. The new output directory
must not contain an existing matching background. Physical C5 remains blocked
by spinor normalisation, the candidate quartic mapping, and complete hard
matching / EFT1 Wilson operator coverage.
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:sys.path.insert(0,str(ROOT))

from Numerical.experimental.Batch17RebuildBackground import rebuild_a_m4_background


def main(argv=None) -> int:
    from Numerical.joint_scalar_matching.fermion_boundary_bridge import load_release_manifest
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--source-release',type=Path,default=ROOT/'output/general_observable/lo_fermion_boundary/releases/source_021dcc1e8fe8ba0f')
    ap.add_argument('--batch14',type=Path,default=ROOT/'output/batch14/eft1/A_-4_EFT1_running.json')
    ap.add_argument('--batch15',type=Path,default=ROOT/'output/batch15/scalar_threshold/A_-4_scalar_threshold_handoff.json')
    ap.add_argument('--batch16',type=Path,default=ROOT/'output/batch16/scalar_mass_audit/A_-4_mass_hessian.json')
    ap.add_argument('--quartic-cache',type=Path,default=ROOT/'output/general_observable/lo_joint_scalar_matching/source_quartic_v4/native_quartic_channels/A_m4.json')
    ap.add_argument('--quartic-seed',type=Path,default=ROOT/'output/full16_symbolic/T3_A_alpha_m4/data/eft1_after_F_scalar_quartic_seed.json')
    ap.add_argument('--registry',type=Path,default=ROOT/'RGE/general/generated/general_rgbeta_quartic_registry.json')
    ap.add_argument('--out-dir',type=Path,default=ROOT/'output/batch17/rebuilt_A_m4')
    ap.add_argument('--background-only',action='store_true',help='Validate/rebuild HHSS background without running C5')
    args=ap.parse_args(argv)
    manifest=load_release_manifest(args.source_release)
    source_path=manifest['A_-4']
    needed=(args.batch14,args.batch15,args.batch16,args.quartic_cache,args.quartic_seed,args.registry,source_path)
    missing=[str(p) for p in needed if not p.is_file()]
    if missing:
        raise FileNotFoundError('Required source files missing:\n'+'\n'.join(missing))
    bg=args.out_dir/'A_m4.json'
    build=rebuild_a_m4_background(
        endpoint_path=args.batch14,handoff_path=args.batch15,mass_path=args.batch16,
        source_path=source_path,quartic_cache_path=args.quartic_cache,
        quartic_seed_path=args.quartic_seed,registry_path=args.registry,
        destination=bg)
    print('[A_-4] REBUILT | channels={} | new Batch14 scalar mass parity=True'.format(build['n_quartic_channels']),flush=True)
    if args.background_only:
        print('Background:',bg)
        return 0
    from Numerical.experimental.Batch17ConditionalWilson import run_conditional
    result=run_conditional('A_-4',handoff_path=args.batch15,mass_path=args.batch16,
        endpoint_path=args.batch14,source_path=source_path,background_path=bg)
    output=args.out_dir/'A_-4_conditional_c5.json'
    if output.exists():
        raise FileExistsError('Refusing to overwrite previous candidate C5 result')
    output.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n',encoding='utf-8')
    summary={
        'model':'A_-4',
        'status':'ConditionalLLSSMatchingAndSMEFTIntegrated_NotPhysicalC5',
        'background':build,
        'conditional_c5_result':str(output),
        'physical_C5_certified':False,
        'historical_fermion_boundary_matches_new_uv_benchmark_certified':False,
        'quartic_normalisation_certified':False,
        'smeft_nfev':result['smeft_nfev'],
        'candidate_neutrino_masses_per_unit_N_eV':result['candidate_neutrino_masses_per_unit_N_eV'],
    }
    path=args.out_dir/'summary.json'
    path.write_text(json.dumps(summary,indent=2,allow_nan=False)+'\n',encoding='utf-8')
    print('[A_-4] CONDITIONAL | LLSS joint scalar + SMEFT; physical C5=UNCERTIFIED',flush=True)
    print('Summary:',path)
    return 0


if __name__=='__main__':
    raise SystemExit(main())
