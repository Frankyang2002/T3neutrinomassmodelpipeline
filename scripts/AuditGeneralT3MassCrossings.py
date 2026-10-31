"""Read-only audit of existing 16-model candidate scalar RG trajectory reports.

The report is a numerical threshold-scheduling diagnostic, NOT an EFT runner.
It does not use or guess finite matching coefficients or modify source reports.

Run from repo root:
  python -m scripts.AuditGeneralT3MassCrossings --all
  python -m scripts.AuditGeneralT3MassCrossings --model B_-1
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any, Mapping

from Numerical.experimental.GeneralT3ScalarQuadraticMasses import load_manifest
from Numerical.experimental.GeneralT3SpectrumThresholdAdapter import (
    DEFAULT_MANIFEST, RunningMassSample, scalar_spectrum,
    ordered_mass_events, find_mass_crossings,
)

ROOT=Path(__file__).resolve().parents[1]
SOURCE_ROOT=(ROOT/'output/studies/general_scalar_extended/'
             'continuous_full_scalar/continuous_v1')
STATUS='CandidateMassThresholdAuditNoMatching'


def complex_pair(x: Any) -> complex:
    if not isinstance(x,dict) or set(x)!= {'real','imag'}:
        raise ValueError('Expected explicit {real,imag} complex coefficient')
    return complex(float(x['real']),float(x['imag']))


def scalar_trajectory(rows: Any) -> list[RunningMassSample]:
    if not isinstance(rows,list) or len(rows)<2:
        raise ValueError('Missing sampled scalar RG trajectory')
    return [RunningMassSample(
        mu_gev=float(r['mu_gev']),
        m_s1_sq_gev2=float(r['mS1Sq_gev2']),
        m_s2_sq_gev2=float(r['mS2Sq_gev2']),
        mu2_gev2=complex_pair(r['mu2_gev2']),
    ) for r in rows]


def audit_record(data: Mapping, manifest: Mapping) -> dict:
    if data.get('format')!='t3-continuous-general-scalar-candidate-v1' or \
       data.get('status')!='CandidateContinuousUVFermionThresholdEFT1Integrated':
        raise ValueError('Input not an existing continuous-v1 candidate report')
    model=data['model']
    if model not in manifest['models']:
        raise ValueError('Input model missing from 16-model manifest')
    uv=scalar_trajectory(data['uv_trajectory'])
    eft=scalar_trajectory(data['eft1_trajectory'])
    # This saved file contains the heavy F masses at the grouped F boundary,
    # but not their continuously evolving UV Takagi masses. No UV fermion
    # crossing is inferred from these endpoint masses alone.
    state=data['uv_fermion_threshold_state']
    fermion_masses=data['heavy_fermion_threshold_masses_gev']
    spec=scalar_spectrum(model,float(state['mS1Sq']),float(state['mS2Sq']),
                         complex_pair(data['quadratic_mu2_at_fermion_threshold_gev2']),
                         manifest=manifest)
    try:
        events=[dict(particle=e.particle,kind=e.kind,mass_gev=e.mass_gev,
                     real_dof=e.real_dof,status=e.status)
                for e in ordered_mass_events(spec,fermion_masses)]
        at_f='positive'
    except ValueError as err:
        events=[]
        at_f='blocked_nonpositive_spectrum: '+str(err)
    result={
       'model':model,'status':STATUS,
       'quadratic_kind':spec.kind,
       'at_grouped_F_mass_order':events,
       'grouped_F_spectrum_status':at_f,
       'uv_scalar_crossings':[vars(x) for x in find_mass_crossings(model,uv,manifest=manifest)],
       'eft1_scalar_crossings':[vars(x) for x in find_mass_crossings(model,eft,manifest=manifest)],
       'uv_fermion_crossings_available':False,
       'physical_mass_basis_certified':False,
       'finite_matching_applied':False,
       'note':'Scalar crossings are candidate unbroken-phase quadratic mu=M(mu) diagnostics within saved stages. They do not change the legacy F-first EFT sequence or integrate any matching.',
    }
    return result


def main() -> int:
    p=argparse.ArgumentParser(description=__doc__)
    choice=p.add_mutually_exclusive_group(required=True)
    choice.add_argument('--all',action='store_true')
    choice.add_argument('--model')
    p.add_argument('--source-root',type=Path,default=SOURCE_ROOT)
    p.add_argument('--output',type=Path,default=None,
                   help='New report JSON path; existing output is never overwritten')
    a=p.parse_args()
    manifest=load_manifest(DEFAULT_MANIFEST)
    selected=tuple(manifest['models']) if a.all else (a.model,)
    if any(m not in manifest['models'] for m in selected):
        p.error('Unknown ordinary T3 model')
    results=[]
    for model in selected:
        path=a.source_root/model.replace('+','p').replace('-','m')/'continuous_trajectory.json'
        if not path.is_file():
            results.append({'model':model,'status':'MissingSavedContinuousTrajectory',
                            'source':str(path)})
            continue
        try:
            results.append(audit_record(json.loads(path.read_text(encoding='utf-8-sig')),
                                        manifest))
        except (ValueError,KeyError,TypeError,RuntimeError) as error:
            results.append({'model':model,'status':'Blocked','reason':str(error),
                            'source':str(path)})
    output=a.output or (a.source_root/'workstream_d_threshold_audit_v1/audit_summary.json')
    if output.exists():
        p.error('Refusing to overwrite existing output: '+str(output))
    output.parent.mkdir(parents=True,exist_ok=True)
    payload={'format':'general-t3-mass-crossing-audit-v1',
             'physical_matching_certified':False,'source_root':str(a.source_root),
             'models':results}
    temporary=output.with_name(output.name+'.tmp')
    temporary.write_text(json.dumps(payload,indent=2,allow_nan=False)+'\n',encoding='utf-8')
    temporary.replace(output)
    for item in results:
        print(item['model'],item['status'])
    print('New candidate threshold audit:',output)
    return 2 if any(x['status']!=STATUS for x in results) else 0


if __name__=='__main__':
    raise SystemExit(main())
