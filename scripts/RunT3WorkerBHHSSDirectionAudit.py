"""Read-only A_-4/D_-2 source-vs-Batch3B HHSS *direction* audit.

Compares derivatives with respect to each running quartic coordinate. Unlike
fitting one global scalar, this identifies unmatched/misnormalised couplings.
Both inputs remain candidate conventions: no physical normalization certified.
"""
from __future__ import annotations

import argparse
from hashlib import sha256
import json
from pathlib import Path
from typing import Mapping

import numpy as np

MODELS = ('A_-4', 'D_-2')
SCHEMA = 't3-worker-b-hhss-direction-audit-v1'
SHAPE = (4, 4, 8, 8)


def number(raw):
    if isinstance(raw, Mapping):
        if not {'real', 'imag'} <= set(raw):
            raise ValueError('complex scalar must contain real/imag')
        z = complex(float(raw['real']), float(raw['imag']))
    else:
        z = complex(raw)
    if not np.isfinite([z.real, z.imag]).all():
        raise ValueError('nonfinite quartic coordinate')
    return z


def put(state, name, z, complex_coordinate):
    out = dict(state)
    if complex_coordinate:
        out[name] = {'real': float(z.real), 'imag': float(z.imag)}
    else:
        if abs(z.imag) > 1e-12:
            raise ValueError(f'{name}: real coordinate received imaginary value')
        out[name] = float(z.real)
    return out


def relative(a, b):
    den = float(np.linalg.norm(b))
    return float(np.linalg.norm(a-b)/den) if den else None


def metrics(candidate, source):
    c=np.asarray(candidate,float); s=np.asarray(source,float)
    if c.shape != SHAPE or s.shape != SHAPE or not np.isfinite(c).all() or not np.isfinite(s).all():
        raise ValueError('expected two finite 4x4x8x8 HHSS tensors')
    nc = float(np.linalg.norm(c)); ns = float(np.linalg.norm(s))
    if ns:
        scale = float(np.vdot(s, c).real/(ns*ns))
        residual = float(np.linalg.norm(c-scale*s)/ns)
        cosine = float(np.vdot(s,c).real/(nc*ns)) if nc else None
    else:
        scale=residual=cosine=None
    return {'candidate_norm':nc,'source_norm':ns,
            'candidate_over_source_scale':scale,
            'shape_residual_relative_to_source':residual,
            'candidate_vs_source_relative':relative(c,s),
            'cosine':cosine,
            'source_direction_absent':bool(ns<=1e-12),
            'candidate_direction_absent':bool(nc<=1e-12)}


def check_symmetry(q, label):
    q=np.asarray(q,float)
    if q.shape != SHAPE or not np.isfinite(q).all():
        raise ValueError(f'{label}: invalid HHSS shape/values')
    norm = max(float(np.linalg.norm(q)), 1e-100)
    err=max(float(np.linalg.norm(q-q.swapaxes(0,1))/norm),
            float(np.linalg.norm(q-q.swapaxes(2,3))/norm))
    if err > 1e-8: raise ValueError(f'{label}: Bose residual {err:.3g}')
    return err


def audit_responses(state, candidate, candidate_directions, candidate_columns,
                    direct_source, *, tolerance=2e-8):
    """Source evaluator is independent of candidate columns; both use same fixed basis."""
    c=np.asarray(candidate,float)
    columns=np.asarray(candidate_columns,float)
    directions=tuple(tuple(x) for x in candidate_directions)
    if c.shape!=SHAPE or columns.shape!=(len(directions),*SHAPE):
        raise ValueError('candidate HHSS direction shapes mismatch')
    if len(set(directions))!=len(directions):
        raise ValueError('duplicate candidate directions')
    for q in [c, *columns]: check_symmetry(q,'candidate')
    base=np.asarray(direct_source(state),float)
    check_symmetry(base,'source reference')
    keys=sorted(k for k,v in state.items() if k.startswith('lambda') and v is not None)
    if not keys: raise ValueError('no running scalar quartic inventory')
    for name,axis in directions:
        if name not in keys or axis not in ('real','imag'):
            raise ValueError(f'bad candidate direction {name}/{axis}')
    union=[]
    direction_set=set(directions)
    for name in keys:
        z=number(state[name])
        is_complex=(isinstance(state[name],Mapping) or abs(z.imag)>0 or
                    (name,'imag') in direction_set)
        union.append((name,'real'))
        if is_complex: union.append((name,'imag'))
    lookup={direction:columns[j] for j,direction in enumerate(directions)}
    source_cols={}; records={}; zero=np.zeros(SHAPE,dtype=float)
    for name,axis in union:
        z=number(state[name]); is_complex=(name,'imag') in union
        step=1j if axis=='imag' else 1.
        probe=put(state,name,z+step,is_complex)
        source_col=np.asarray(direct_source(probe),float)-base
        check_symmetry(source_col,f'source {name}.{axis}')
        source_cols[name,axis]=source_col
        col=lookup.get((name,axis),zero)
        record=metrics(col,source_col)
        record['coupling']=name; record['axis']=axis
        record['source_coefficient_at_endpoint']=(z.imag if axis=='imag' else z.real)
        record['source_contribution_norm']=abs(record['source_coefficient_at_endpoint'])*record['source_norm']
        record['candidate_contribution_norm']=abs(record['source_coefficient_at_endpoint'])*record['candidate_norm']
        record['candidate_cached_direction_present']=(name,axis) in lookup
        records[f'{name}.{axis}']=record
    # Independent multi-axis affine interpolation; source evaluator evaluated
    # at a new point instead of merely resumming its one-hot probes.
    trial=dict(state); predicted=base.copy()
    for j,(name,axis) in enumerate(union):
        d=0.07*np.sin(j+1)
        z=number(trial[name]); is_complex=(name,'imag') in union
        trial=put(trial,name,z+(1j*d if axis=='imag' else d),is_complex)
        predicted+=d*source_cols[name,axis]
    actual=np.asarray(direct_source(trial),float)
    check_symmetry(actual,'source multi-axis probe')
    interpolation=relative(predicted,actual)
    if interpolation is None or interpolation>tolerance:
        raise ValueError(f'source HHSS not affine at multi-axis probe: {interpolation}')
    # Cached candidate consistency at its endpoint already follows from the
    # loaded Batch3B hash; here compare independent source endpoint values.
    records_sorted=sorted(records.values(),
       key=lambda x: (-(x['source_contribution_norm']+x['candidate_contribution_norm']),x['coupling'],x['axis']))
    return {'endpoint_comparison':metrics(c,base),
            'source_affine_interpolation_relative':interpolation,
            'directions':records_sorted,
            'candidate_direction_count':len(directions),
            'source_evaluated_direction_count':len(union),
            'candidate_only_directions':[list(x) for x in directions if x not in union],
            'source_only_directions':[list(x) for x in union if x not in lookup],
            'source_hhss_reconstructed_independently':True,
            'physical_quartic_normalization_certified':False}


class SourceEvaluator:
    """Prepare source inputs ONCE; all probes reuse identical CG/registry geometry."""
    def __init__(self, model, release_dir, v4_dir, seed_dir, registry_path, reference_state):
        from Numerical.joint_scalar_matching.fermion_boundary_bridge import (load_release_manifest,load_source)
        from Numerical.lo_certification.couplings.extra_hhss_geometry import model_extra_geometry
        from scripts import RunT3GeneralSequentialFitAblation as original
        from Numerical.joint_scalar_matching.special_mass_v6 import build_special_basis
        self.model=model
        self.registry_path=Path(registry_path)
        self.registry=json.loads(self.registry_path.read_text(encoding='utf-8-sig'))
        manifest=load_release_manifest(release_dir)
        self.boundary=load_source(manifest[model],expected_model=model)
        self.channels=original._load_channels(model,self.boundary,v4=v4_dir,seeds=seed_dir)
        self.geometry=model_extra_geometry(model,self.registry['models'][model],include_components=True)
        self.rotation=None
        if model=='D_-2':
            reference_snapshot=original.snapshot(model,dict(reference_state),self.registry_path)
            self.rotation=np.asarray(build_special_basis(model,reference_snapshot,self.boundary)['R_original_from_boundary'],float)
            if self.rotation.shape!=(8,8) or not np.allclose(self.rotation.T@self.rotation,np.eye(8),atol=1e-10):
                raise ValueError('invalid special scalar basis rotation')
        self.evaluations=0

    def __call__(self, state):
        from Numerical.running.GeneralT3SequentialLeadingLog import scalar_tensors
        from scripts import RunT3GeneralSequentialFitAblation as original
        from scripts.RunT3WorkerBSourceHHSSAudit import fixed_native_hhss
        snap=original.snapshot(self.model,dict(state),self.registry_path)
        bundle=scalar_tensors(self.model,self.boundary,self.channels,snap,self.registry,geometry=self.geometry)
        h=fixed_native_hhss(bundle.quartic,self.rotation)
        self.evaluations+=1
        return h


def _pairs(items,models):
    result={}
    for x in items:
        model,sep,path=x.partition('=')
        if not sep or model not in models or not path or model in result:
            raise ValueError('expected exactly one --state-file MODEL=PATH per model')
        result[model]=Path(path)
    return result


def main(argv=None):
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--project',type=Path,default=Path(__file__).resolve().parents[1])
    p.add_argument('--models',default='A_-4,D_-2')
    p.add_argument('--response-root',type=Path,required=True)
    p.add_argument('--state-file',action='append',default=[],metavar='MODEL=PATH')
    p.add_argument('--release-dir',type=Path,default=Path('output/general_observable/lo_fermion_boundary/releases/source_021dcc1e8fe8ba0f'))
    p.add_argument('--v4-dir',type=Path,default=Path('output/general_observable/lo_joint_scalar_matching/source_quartic_v4/native_quartic_channels'))
    p.add_argument('--seed-dir',type=Path,default=Path('output/full16_symbolic'))
    p.add_argument('--registry',type=Path,default=Path('RGE/general/generated/general_rgbeta_quartic_registry.json'))
    p.add_argument('--out',type=Path,default=Path('output/audits/workerb_hhss_direction_v1/report.json'))
    args=p.parse_args(argv)
    models=tuple(x.strip() for x in args.models.split(',') if x.strip())
    if not models or len(models)!=len(set(models)) or not set(models)<=set(MODELS):p.error('unique subset of A_-4,D_-2 required')
    project=args.project.resolve()
    def absolute(path):
        x=Path(path)
        return (x if x.is_absolute() else project/x).resolve()
    state_paths=_pairs(args.state_file,models)
    if set(state_paths)!=set(models):p.error('explicit --state-file required per model')
    out=absolute(args.out)
    if out.exists():p.error(f'refusing to overwrite historical output {out}')
    response_root=absolute(args.response_root)
    srcpaths=[absolute(args.release_dir),absolute(args.v4_dir),absolute(args.seed_dir),absolute(args.registry)]
    from Numerical.running.T3CoupledQuarticContext import CandidateQuarticContext
    from scripts.RunT3WorkerBEndpointAudit import _extract_endpoint
    results={}
    for model in models:
        path=absolute(state_paths[model]);mu,state,h=_extract_endpoint(path,model)
        ctx=CandidateQuarticContext.load(project,response_root,model)
        reference=ctx.full_reference(state)[:4,:4,4:,4:]
        columns=np.array(ctx.columns[:,:4,:4,4:,4:],float)
        direct=SourceEvaluator(model,*srcpaths,state)
        diagnostics=audit_responses(state,reference,ctx.directions,columns,direct)
        diagnostics.update({'model':model,'mu_s_gev':mu,'endpoint_path':str(path),
                            'endpoint_sha256':h,'source_evaluator_calls':direct.evaluations,
                            'basis_convention':'native fixed real coordinates, D_-2 source mass-basis back rotation',
                            'source_release_path':str(direct.boundary.path),
                            'source_quartic_seed_path':str(direct.channels.seed_path)})
        results[model]=diagnostics
        worst=sorted(diagnostics['directions'],key=lambda row:-(row['shape_residual_relative_to_source'] or 0))[:3]
        print(f'{model}: source affine={diagnostics["source_affine_interpolation_relative"]:.3e}; '
              f'candidate/source endpoint={diagnostics["endpoint_comparison"]["candidate_vs_source_relative"]:.6g}; '
              f'directions={len(diagnostics["directions"])}')
        for row in worst:
            print(f'  {row["coupling"]}.{row["axis"]}: scale={row["candidate_over_source_scale"]}, '
                  f'shape={row["shape_residual_relative_to_source"]}')
    doc={'schema':SCHEMA,'status':'ComparedCandidateDirectionGeometry_NotCertified',
         'response_root':str(response_root),'results':results,
         'native_rgbeta_mapping_certified':False,'physical_quartic_certified':False,
         'fermion_or_joint_scalar_matching_performed':False}
    out.parent.mkdir(parents=True,exist_ok=True)
    out.write_text(json.dumps(doc,indent=2,allow_nan=False)+'\n',encoding='utf-8')
    print('Report:',out)
    return 0


if __name__=='__main__':
    raise SystemExit(main())
