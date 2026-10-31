"""Stage 4 independent read-only numerical validation for Workstream D.

Checks Stage 2 -> Stage 3 candidate inventory provenance and both pre-existing
B/C conditional C5/mnu benchmark reports. Does not integrate new RGEs, derive
missing matching, or modify prior numerical results.

  python -m scripts.AuditGeneralT3NumericalCompletion
  python -m scripts.AuditGeneralT3NumericalCompletion --bc-reports PATH --output NEW_PATH
"""
from __future__ import annotations

import argparse
from hashlib import sha256
import json
from pathlib import Path
from typing import Any

from Numerical.experimental.GeneralT3EFTHandoffController import (
    candidate_inventory_from_stage2,
)
from Numerical.experimental.GeneralT3NumericalChecks import (
    EXPECTED_BC_MASSES_EV, check_bc_conditional_report,
)

ROOT = Path(__file__).resolve().parents[1]
STUDY = (ROOT / 'output/studies/general_scalar_extended/continuous_full_scalar/continuous_v1')
SOURCE_STAGE2 = STUDY / 'workstream_d_fermion_audit_v2/audit_summary.json'
SOURCE_STAGE3 = STUDY / 'workstream_d_handoff_readiness_v1/summary.json'
SOURCE_BC = STUDY / 'bc_eft2_d_regression_v1'
DEFAULT_OUT = STUDY / 'workstream_d_numerical_validation_v1/summary.json'
FORMAT = 'general-t3-workstream-d-stage4-contract-validation-v1'


def _load(path: Path) -> tuple[dict[str, Any], str]:
    if not path.is_file():
        raise ValueError('Missing existing source report: '+str(path))
    raw = path.read_bytes()
    data = json.loads(raw.decode('utf-8-sig'))
    if not isinstance(data, dict):
        raise ValueError('Source report must be a JSON object: '+str(path))
    return data, sha256(raw).hexdigest()


def verify_inventory_pair(stage2: dict, stage3: dict, *, stage2_sha256: str) -> dict:
    """Check completeness, root order, status and exact Stage 3 source hash."""
    if stage2.get('format') != 'general-t3-running-heavy-mass-crossing-audit-v1' or \
            stage2.get('status') != 'CandidateAuditComplete' or \
            stage2.get('physical_matching_certified') is not False:
        raise ValueError('Stage 2 source is not the complete conditional candidate audit')
    if stage3.get('format') != 'general-t3-eft-handoff-readiness-v1' or \
            stage3.get('status') != 'CandidateHandoffDependenciesAuditedNotExecuted' or \
            stage3.get('physical_matching_certified') is not False or \
            stage3.get('source_sha256') != stage2_sha256 or \
            stage3.get('matching_callbacks_registered') != 0 or \
            stage3.get('evolution_callbacks_registered') != 0:
        raise ValueError('Stage 3 handoff audit provenance/status/callback inventory changed')
    source = stage2.get('models')
    readiness = stage3.get('models')
    if not isinstance(source,list) or not isinstance(readiness,list) or len(source)!=16 or len(readiness)!=16:
        raise ValueError('Stage 4 requires all 16 trial models')
    if len({r['model'] for r in source}) != 16 or len({r['model'] for r in readiness}) != 16:
        raise ValueError('Repeated candidate model identifier')
    stage3_by_model = {r['model']:r for r in readiness}
    total = 0
    missing: dict[str, list[str]] = {}
    for row in source:
        inv = candidate_inventory_from_stage2(row)
        rec = stage3_by_model.get(inv.model)
        if not isinstance(rec,dict) or rec.get('status') != 'MissingMatchingAndIntermediateRGEs' or \
                rec.get('candidate_mass_inventory_complete') is not inv.complete or \
                rec.get('candidate_crossing_count') != len(inv.events) or \
                rec.get('unresolved_mass_branches') != list(inv.unresolved) or \
                rec.get('runtime_handoffs_executed') != 0 or \
                rec.get('intermediate_particle_removed_rge_callbacks_registered') != 0 or \
                rec.get('threshold_matching_callbacks_registered') != 0 or \
                rec.get('physical_matching_certified') is not False:
            raise ValueError(f'{inv.model}: readiness record disagrees with Stage 2')
        published = rec.get('ordered_candidate_events')
        if not isinstance(published,list) or len(published)!=len(inv.events):
            raise ValueError(f'{inv.model}: readiness events missing')
        for a,b in zip(published,inv.events):
            if (a.get('particle'),a.get('kind'),a.get('scale_gev'),a.get('status')) != \
                    (b.particle,b.kind,b.scale_gev,b.source_status):
                raise ValueError(f'{inv.model}: Stage 3 events differ from Stage 2')
        total += len(inv.events)
        if inv.unresolved:
            missing[inv.model] = list(inv.unresolved)
    if stage3.get('models_checked') != 16 or stage3.get('candidate_mass_roots_found') != total or \
            stage3.get('unresolved_candidate_branches') != sum(map(len,missing.values())) or \
            stage3.get('models_with_all_candidate_branches') != 16-len(missing):
        raise ValueError('Stage 3 aggregate counts disagree with validated event inventory')
    return dict(models_checked=16,candidate_mass_roots_found=total,
                models_with_complete_candidate_inventories=16-len(missing),
                unresolved_candidate_branches=missing,
                sequential_particle_matching_executed=False)


def build_summary(stage2: dict, stage3: dict, reports: dict[str,dict], *,
                  stage2_sha256: str) -> dict:
    source = verify_inventory_pair(stage2,stage3,stage2_sha256=stage2_sha256)
    if set(reports) != set(EXPECTED_BC_MASSES_EV):
        raise ValueError('Both B/C original or freshly reproduced conditional reports are required')
    bc = {model:check_bc_conditional_report(reports[model],expected_model=model)
          for model in sorted(EXPECTED_BC_MASSES_EV)}
    return {
        'format':FORMAT,
        'status':'NumericalInterfacesValidatedPhysicalMatchingMissing',
        'source_stage2_sha256':stage2_sha256,
        'physical_neutrino_predictions_certified':False,
        'full_sequential_particle_decoupling_performed':False,
        'running_and_matching_callbacks_supplied_by_B_C':False,
        'numerical_inventory':source,
        'bc_conditional_c5_takagi_regressions':bc,
        'remaining_physics_dependencies':[
            'B: verified EFT RGEs below individual heavy-fermion thresholds and scalar thresholds',
            'C: certified per-particle matching including B/C finite Y and X corrections',
            'A/B/C: verified operator-field basis maps and Wilson coefficient normalizations',
        ],
    }


def main() -> None:
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--stage2',type=Path,default=SOURCE_STAGE2)
    p.add_argument('--stage3',type=Path,default=SOURCE_STAGE3)
    p.add_argument('--bc-reports',type=Path,default=SOURCE_BC)
    p.add_argument('--output',type=Path,default=DEFAULT_OUT)
    args=p.parse_args()
    if args.output.exists():
        p.error('Refusing to replace existing report: '+str(args.output))
    stage2,h2=_load(args.stage2)
    stage3,_=_load(args.stage3)
    reports={}
    for model in EXPECTED_BC_MASSES_EV:
        tag=model.replace('-','m').replace('+','p')
        reports[model],_=_load(args.bc_reports/f'{tag}_eft2_sequential.json')
    summary=build_summary(stage2,stage3,reports,stage2_sha256=h2)
    args.output.parent.mkdir(parents=True,exist_ok=True)
    with args.output.open('x',encoding='utf-8') as f:
        json.dump(summary,f,indent=2,ensure_ascii=False,allow_nan=False)
        f.write('\n')
    print('Stage 4: models=',summary['numerical_inventory']['models_checked'],
          'candidate roots=',summary['numerical_inventory']['candidate_mass_roots_found'],
          'missing=',summary['numerical_inventory']['unresolved_candidate_branches'])
    for model,rec in summary['bc_conditional_c5_takagi_regressions'].items():
        print(model,'conditional masses(eV)=',rec['masses_ev_descending'],
              'finite matching=MISSING')
    print('Stage 4 numerical validation:',args.output)


if __name__=='__main__':
    main()
