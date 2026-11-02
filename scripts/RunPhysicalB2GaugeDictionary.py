r"""B2 conditional all-16 audit. Read-only, fresh output directory, no overwrites.

PowerShell: .\.venv\Scripts\python.exe scripts\RunPhysicalB2GaugeDictionary.py
           --output output\general_observable\physical_normalisation\b2_local_01
"""
from __future__ import annotations
import argparse
from pathlib import Path
import json
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
from Numerical.lo_certification.physical_normalisation.b2_gauge_coupling import ( # noqa:E402
    compute, HEAD, AuditBlocked, spurion_phase_ward,
)


def atomic_new_json(path:Path, value:dict) -> None:
    if path.exists():
        raise AuditBlocked('Refusing to overwrite pre-existing evidence: ' + str(path))
    path.write_text(json.dumps(value,indent=2,sort_keys=True) + '\n',encoding='utf-8')


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--project-root',type=Path,default=ROOT)
    parser.add_argument('--output',type=Path,default=Path(
        'output/general_observable/physical_normalisation/b2_local_01'))
    args=parser.parse_args()
    root=args.project_root.resolve()
    out=(root/args.output).resolve() if not args.output.is_absolute() else args.output.resolve()
    if out.exists():
        print('BLOCKED: output already exists; choose a NEW directory: '+str(out))
        return 2
    if not root.is_dir() or not (root/'Lagrangian').is_dir():
        print('BLOCKED: project root is not available: '+str(root));return 2
    try:
        dictionary, channels=compute(root)
        # Use the project-root-validated HEAD measured by source_inventory.
        # An independent `git -C root rev-parse HEAD` could accept an ancestor
        # repository and incorrectly overwrite the validated provenance.
        revision = dictionary['source_inventory']['current_git_head']
        dictionary['spurion_dirac']=spurion_phase_ward(False)
        dictionary['spurion_majorana']=spurion_phase_ward(True)
        status=('CONDITIONAL_NATIVE_CG_SNAPSHOT_LOADED' if
                dictionary['native_cg_snapshot_loaded'] else
                'BLOCKED_NATIVE_TENSOR_AUDIT_OR_SOURCE_DRIFT')
        dictionary['status']=status
        out.mkdir(parents=True)
        atomic_new_json(out/'cross_package_dictionary.json',dictionary)
        atomic_new_json(out/'all16_channel_coefficients.json',channels)
        handoff=(
          '# B2 handoff (conditional; NOT physically source-certified)\n\n'
          f'- GitHub inspected main: `{HEAD}`; local git HEAD: `{revision or "unavailable"}`.\n'
          f'- Status: `{status}`.\n'
          '- Source contract fields recorded in `cross_package_dictionary.json`.\n'
          '- Common weight basis comes from the recorded Matchete SU(2) generator audit; '
          'the RGBeta SU(2) component axes are still not measured.\n'
          '- Dirac mass doublet invariant must be antisymmetric epsilon; triplet '
          'mass metric is `U.T*U` in the chosen Matchete weight basis.\n'
          '- Loop spurion **conditional on B1**: `y1*y2*conj(MF)*conj(lambdaT3)`.\n'
          '- RGBeta Yukawa internal sqrt(2) and MF internal factors are NOT physical action maps.\n'
          '- Basic singlet portals previously established for 10 sectors; preserve '
          'nonunique A/S2,D/S1,E/S1,E/S2 Hermitian nullspaces; no physical factor claimed.\n'
          '- `conditional_CG_product` is a gauge-component diagnostic, not a '
          'source-to-source physical neutrino-mass coefficient.\n'
          '- All 16 charge skeletons are algebraic; frozen RGE/matching not run.\n'
          '- Remaining required witness: native numeric RGBeta epsilon/del/delS2 '
          'component values in its actual registered indices + physical action '
          'prefactors and B1 spinor phase, measured against current local sources.\n'
          '- **Safe to integrate**: charge skeletons and guarded source metadata only.\n'
          '- **NOT safe to integrate**: absolute y1/y2/MF/lambdaT3 ratios, complex '
          'relative phases or purported physical mass normalisation.\n'
          '- `physical_crosspackage_complete=false`; `physical_source_certified=false`.\n')
        (out/'HANDOFF.md').write_text(handoff,encoding='utf-8')
        print(status, '; models=16; physical certificates=0; output=',out)
        return 0 if dictionary['native_cg_snapshot_loaded'] else 3
    except (AuditBlocked,ValueError,KeyError,TypeError,OSError) as exc:
        print('BLOCKED: '+str(exc));return 4


if __name__=='__main__':
    raise SystemExit(main())
