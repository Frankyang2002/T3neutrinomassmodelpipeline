"""Worker 5: read-only finite joint-scalar matching completeness audit.

Reuses Worker C's audited kernel and saved profiles. No production imports are
monkey-patched; no Matchete Laurent or hard finite piece is added to C5.
"""
from __future__ import annotations
import argparse
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import re
import numpy as np

ROOT = Path(__file__).resolve().parents[2]
OUTROOT = Path('output/audits/worker5_scalar_matching')
SOURCE_TERMS = ('D_-2', 'y2Conj', 'FullFlavorHardC5Verified', 'FiniteHardMinusBubbleRemainderCertified')


def digest(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()


def read(p):
    return json.loads(p.read_text(encoding='utf-8-sig'))


def mat(x):
    if not isinstance(x, dict) or not {'real', 'imag'} <= x.keys():
        raise ValueError('requires explicit real/imag matrix')
    z = np.asarray(x['real'], float) + 1j * np.asarray(x['imag'], float)
    if z.shape != (3, 3) or not np.isfinite(z).all():
        raise ValueError('expected finite 3x3 C5')
    return z


def matrix_metrics(a, b):
    """A numerical *shape* comparison, never proof of physical independence."""
    a, b = mat(a), mat(b)
    na, nb = np.linalg.norm(a), np.linalg.norm(b)
    ip = np.vdot(a, b)
    overlap = abs(ip)/(na*nb) if na*nb else None
    phase = ip/abs(ip) if abs(ip) else 1.0
    aligned_residual = float(np.linalg.norm(b - (ip/(na*na))*a)/nb) if na*nb else None
    # Complex symmetric flavour tensors only. Orthogonal residual is to a complex line.
    return {'norm_a_gev_inverse':float(na),'norm_b_gev_inverse':float(nb),
            'norm_ratio_b_over_a':float(nb/na) if na else None,
            'absolute_cosine':float(overlap) if overlap is not None else None,
            'least_squares_complex_line_residual':aligned_residual,
            'a_flavour_symmetry_residual':float(np.linalg.norm(a-a.T)/max(na,1e-100)),
            'b_flavour_symmetry_residual':float(np.linalg.norm(b-b.T)/max(nb,1e-100)),
            'inner_product_phase_real':float(phase.real),
            'inner_product_phase_imag':float(phase.imag)}


def archive_probe(root):
    parent = root/'output/audits'
    runs = sorted(parent.glob('worker_c_numeric_overlap_D_m2_*/result.json'))
    results=[]
    for f in runs:
        try:
            r=read(f); c=r['components']; h=c['hard']; t=c['total']
            metrics=matrix_metrics(h['finite_c5_gev_inv'],t['finite_c5_gev_inv'])
            bench=f.parent/'restricted_benchmark.json'
            b=read(bench) if bench.exists() else {}
            equal_scale=float(b.get('MuScalar_GeV',float('nan')))==1000.0
            results.append({'file':str(f.relative_to(root)),'sha256':digest(f),
              'status':r.get('status'), 'model':r.get('model'), 'nf':r.get('NF'),
              'benchmark_mu_s_gev':b.get('MuScalar_GeV'),
              'production_joint_mu_s_gev':1000.0,
              'same_scale_as_production':equal_scale,
              'benchmarks_source_sha256':digest(bench) if bench.exists() else None,
              'hard_total_flavour_alignment_ONLY':metrics,
              'hard_pole_residue_max_abs':h.get('pole_residue_max_abs_gev_inv'),
              'total_pole_residue_max_abs':t.get('pole_residue_max_abs_gev_inv'),
              'normalization_certified':r.get('operator_and_CG_normalization_certified') is True,
              'nonbubble_remainder_certified':r.get('nonbubble_remainder_certified') is True,
              'source_overlap_status':r.get('overlap',{}).get('status'),
              'valid_additional_physical_coefficient':False,
              'interpretation':'Finite-part hard and total source expressions: total includes inherited/running; NOT hard minus bubble.'})
        except (KeyError, ValueError, OSError, TypeError) as exc:
            results.append({'file':str(f.relative_to(root)), 'status':'UNREADABLE','reason':str(exc)})
    return results


def source_profile_audit(root):
    """Run unchanged Worker C gate against each existing complete saved profile."""
    from Numerical.joint_scalar_matching.worker_c_matchete_audit import FILES, audit_matchete_scalar_profile
    candidates=[]
    for f in (root/'data').rglob(FILES['match_result']) if (root/'data').exists() else []:
        if len(candidates)>=30: break
        candidates.append(f.parent)
    records=[]
    for p in candidates:
        try:
            v=audit_matchete_scalar_profile(p)
            records.append({'path':str(p.relative_to(root)), 'status':v['status'],
              'gates':v['gates'], 'candidate_inventory':v['candidate_inventory'],
              'blocking_reasons':v['blocking_reasons'],
              'safe_to_add':v['safe_to_add_as_independent_c5_correction'],
              'sources':v['source_files']})
        except (OSError,ValueError,KeyError) as exc:
            records.append({'path':str(p.relative_to(root)), 'status':'ERROR','reason':str(exc)})
    return records


def saved_restricted(root):
    files=sorted((root/'output/audits').glob('worker_c_restricted_D_m2_*/result.json'))
    rows=[]
    for p in files:
        try:
            r=read(p)
            rows.append({'path':str(p.relative_to(root)),'sha256':digest(p),
              'status':r.get('status'), 'nf':r.get('NF'),
              'restricted_y2conj_zeroed':r.get('independent_y2Conj_set_zero_in_saved_outputs'),
              'total_pole_cancelled_nf1':r.get('restricted_total_pole_zero_exact_NF1'),
              'finite_renormalization_certified':r.get('finite_renormalization_certified'),
              'operator_normalization_certified':r.get('weinberg_operator_normalization_certified'),
              'physical_mass_relations_imposed':r.get('physical_mass_relations_imposed'),
              'nonbubble_remainder_certified':r.get('independent_nonbubble_remainder_certified'),
              'safe_to_add':False})
        except (OSError, ValueError) as exc:
            rows.append({'path':str(p.relative_to(root)), 'error':str(exc)})
    return rows


def audit_source_files(root):
    """Indexed local references, not a new matching calculation."""
    targets=['Numerical/joint_scalar_matching/core.py',
             'Numerical/joint_scalar_matching/worker_c_joint.py',
             'Numerical/joint_scalar_matching/worker_c_matchete_audit.py',
             'Numerical/joint_scalar_matching/worker_c_stream_v2.py',
             'Lagrangian/RunGeneralDMinus2C5PoleCancellationAudit.wl',
             'Lagrangian/RunGeneralDMinus2C5FlavourStructureAudit.wl',
             'Lagrangian/RunGeneralDMinus2C5NormalizationReference.wl']
    rows=[]
    for name in targets:
        p=root/name
        if not p.is_file():
            rows.append({'path':name,'exists':False});continue
        s=p.read_text(encoding='utf-8-sig',errors='replace').splitlines()
        excerpts=[]
        for i,line in enumerate(s,1):
            if any(token.lower() in line.lower() for token in SOURCE_TERMS) or any(
                token in line for token in ('c5_total', 'c5_inherited', 'c5_finite_scalar','match_c5(','match_joint_scalars(')):
                if len(excerpts)<35: excerpts.append({'line':i,'text':line[:240]})
        rows.append({'path':name,'exists':True,'sha256':digest(p), 'line_count':len(s),
                     'matching_lines':excerpts})
    return rows


def endpoint_inventory(root):
    """Only compare when *actual* canonical endpoint JSON/NPZ exists."""
    from Numerical.joint_scalar_matching.worker_c_stream_v2 import validate_endpoint_contract, run_endpoint
    candidates=[]
    for p in (root/'output').rglob('*.json') if (root/'output').exists() else []:
        if len(candidates)>10000: break
        if p.stat().st_size>100000: continue
        try:
            obj=read(p)
            if isinstance(obj,dict) and obj.get('tensor_convention')=='T3_REAL_LLSS_HHSS_WEINBERG_v1' and 'tensor_npz' in obj:
                candidates.append(p)
        except (OSError,ValueError,UnicodeError): pass
    rows=[]
    for p in candidates[:40]:
        try:
            a=validate_endpoint_contract(p)
            result=run_endpoint(p)
            rows.append({'path':str(p.relative_to(root)), 'status':'CANONICAL_ENDPOINT_EVALUATED',
             'model':result['model'], 'mu_s_gev':result['mu_s_gev'],
             'norms':result['c5_norms_gev_inv'],
             'projection_residual':result['weinberg_projection_residual'],
             'source_flags':a['verification_flags'],
             'additional_finite_correction_compared':False})
        except (ValueError,OSError,KeyError,ImportError) as exc:
            rows.append({'path':str(p.relative_to(root)),'status':'ENDPOINT_REJECTED','reason':str(exc)})
    return rows


def run(root, output):
    root=Path(root).resolve()
    if not (root/'Numerical/joint_scalar_matching/core.py').is_file():
        raise FileNotFoundError('Not a checkout containing Numerical/joint_scalar_matching/core.py')
    output=Path(output).resolve()
    if output.exists(): raise FileExistsError(f'Refusing to overwrite previous audit: {output}')
    report={'format':'worker5-joint-scalar-source-audit-v1',
      'timestamp_utc':datetime.now(timezone.utc).isoformat(),
      'repo_root':str(root), 'physics_chain':'UV -> all heavy fermions -> EFT1(SM+S1+S2) -> both scalars jointly at 1000 GeV -> SMEFT',
      'implemented':'C5_SMEFT(muS) = C5_EFT1(muS) + finite LLSS*HHSS zero-momentum joint bubble',
      'production_modified':False,
      'potential_missing':'Other EFT1 D5 insertions, scalar wavefunction/counterterm terms and other scalar topologies: unverified, NOT assumed nonzero or independent',
      'flavour_note':'Hard vs total alignment tests are NOT bubble-subtracted candidate-vs-bubble tests; a new tensor can align with inherited C5 even if diagrams differ.',
      'status':'NO_CERTIFIED_INDEPENDENT_FINITE_REMAINDER',
      'archive_overlap':archive_probe(root),
      'restricted_source':saved_restricted(root),
      'source_profiles':source_profile_audit(root),
      'source_file_inventory':audit_source_files(root),
      'canonical_endpoints':endpoint_inventory(root),
      'gates_for_physical_comparison':[
          'restricted Y1/Y2 theory (no independent y2Conj) and exact scalar mass relations',
          'renormalized and pole-subtracted UV/EFT1 minus SMEFT amplitude at common muS=1000 GeV',
          'consistent Lorentz, Weyl, CG, Weinberg-operator and finite-scheme normalization',
          'explicit subtraction of inherited EFT1 and already implemented LLSS*HHSS bubble',
          'same-coupling canonical evolved C5 and scalar background at muS',
          'complex-symmetric non-overlapping finite C5 matrix with independently checked flavour tensor'],
      'no_raw_laurent_promoted':True, 'no_physical_missing_matrix_invented':True}
    output.mkdir(parents=True,exist_ok=False)
    (output/'audit.json').write_text(json.dumps(report,indent=2,allow_nan=False)+'\n',encoding='utf-8')
    summary=[f'# Worker 5 scalar threshold audit', '', f'Status: {report["status"]}',
      'Production: unchanged. Both scalar multiplets jointly integrated out at 1000 GeV.',
      f'Archived provisional D(-2) comparisons: {len(report["archive_overlap"])}',
      f'Restricted source snapshots: {len(report["restricted_source"])}',
      f'Worker C saved profile gates evaluated: {len(report["source_profiles"])}',
      f'Canonical endpoint candidates: {len(report["canonical_endpoints"])}',
      '', 'No new physically normalized finite threshold correction established.',
      'See audit.json for hashes, pole/scale/flavour diagnostics, and explicit blockers.']
    (output/'SUMMARY.md').write_text('\n'.join(summary)+'\n',encoding='utf-8')
    return report


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root',type=Path,default=ROOT)
    parser.add_argument('--output',type=Path,default=None)
    a=parser.parse_args()
    dest=a.output or (a.root/OUTROOT/datetime.now(timezone.utc).strftime('%Y%m%d_%H%M%S_%f'))
    v=run(a.root,dest)
    print(json.dumps({'output':str(dest),'status':v['status'],
      'archives':len(v['archive_overlap']),'source_profiles':len(v['source_profiles']),
      'canonical_endpoints':len(v['canonical_endpoints'])},indent=2))

if __name__=='__main__': main()
