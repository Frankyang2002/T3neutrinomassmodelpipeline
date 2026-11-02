"""Check whether native Matchete Higgs-portal CG tensors span RGBeta's identity.

Independent OFF-SHELL four-index test; does not change any benchmark, UV input,
matching, running or neutrino mass result. A successful geometric identity is
NOT a physical cross-package coupling sign/normalization certificate.

Project path: scripts/ResolveT3DirectUVPortalIdentity.py
Usage (repository root):
 python scripts/ResolveT3DirectUVPortalIdentity.py \
  --project-root . \
  --portal-root output/general_observable/lo_direct_one_loop/source_portal_cg_all16 \
  --audit-dir output/general_observable/lo_direct_one_loop/source_audits
"""
from __future__ import annotations

import argparse
from hashlib import sha256
import json
from pathlib import Path
import re
import sys

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
for _root in (str(ROOT), str(ROOT / "scripts")):
    if _root not in sys.path:
        sys.path.insert(0, _root)

from scripts.RunT3DirectOneLoop import models, CLASSES  # noqa: E402
from scripts.BuildT3DirectOneLoopInputs import _rotations, _tensor  # noqa: E402
from scripts.AuditT3DirectUVPortalBasis import FLAGS, source_contract  # noqa: E402


FORMAT = 't3-direct-uv-portal-offshell-identity-span-v1'
TOL = 2e-10


def require(ok: bool, reason: str) -> None:
    if not ok:
        raise ValueError(reason)


def read_json(path: Path) -> dict:
    require(path.is_file(), f'Missing required source: {path}')
    obj = json.loads(path.read_text(encoding='utf-8-sig'))
    require(isinstance(obj, dict), f'Expected JSON object: {path}')
    return obj


def hash_file(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def rgbeta_identity_source_contract(project: Path) -> dict:
    """Check that the compared RGBeta structure is actually the pair identity."""
    path = project / 'RGE/running/rgbeta/T3RGBetaModel.wl'
    text = path.read_text(encoding='utf-8-sig')
    required = (
        'T3RGBetaPairInvariant[1, i_, j_] := 1;',
        'T3RGBetaPairInvariant[2, i_, j_] := del[SU2L @ fund, i, j];',
        'T3RGBetaPairInvariant[3, i_, j_] := del[SU2L @ S2, i, j];',
        'T3RGBetaPairInvariant[dS1, #3, #4] &',
        'T3RGBetaPairInvariant[dS2, #3, #4] &',
        'del[SU2L @ fund, #1, #2] *',
    )
    require(all(s in text for s in required), 'RGBeta pair-invariant source changed')
    for coupling in ('lambdaH1', 'lambdaH2'):
        require(bool(re.search(r'AddQuartic\s*\[\s*'+coupling+r'\s*,\s*\{\s*Bar\s*@\s*H\s*,\s*H\s*,\s*Bar\s*@\s*S[12]\s*,\s*S[12]', text)),
                f'RGBeta {coupling} field orientation changed')
    return {'path': str(path.relative_to(project)), 'sha256': hash_file(path)}


def hermitian_portal(t: np.ndarray) -> np.ndarray:
    """Real-coupling 1/2 PlusHc action in BarH,H,BarS,S order."""
    t = np.asarray(t, dtype=np.complex128)
    require(t.ndim == 4 and t.shape[0:2] == (2, 2)
            and t.shape[2] == t.shape[3], 'Invalid native portal tensor shape')
    require(np.isfinite(t).all(), 'Nonfinite native portal tensor')
    return (t + t.transpose(1, 0, 3, 2).conj())/2


def identity_tensor(dim: int) -> np.ndarray:
    return np.einsum('ij,ab->ijab', np.eye(2), np.eye(dim)).astype(complex)


def fit_native_identity(native_tensors: list[np.ndarray], tol: float = TOL) -> dict:
    """Find REAL native Lagrangian coefficients reproducing delta(H)delta(S).

    Native source terms have 1/2 PlusHc with real lambdaHInvk; hence the
    physical four-index tensors, and fitted coefficients, are real-linear.
    """
    require(1 <= len(native_tensors) <= 2, 'Expect 1 or 2 native portal invariants')
    dim = native_tensors[0].shape[2]
    require(all(t.shape == (2, 2, dim, dim) for t in native_tensors),
            'Inconsistent native portal tensor shapes')
    tensors = [hermitian_portal(t) for t in native_tensors]
    target = identity_tensor(dim)
    a = np.column_stack([t.ravel() for t in tensors])
    b = target.ravel()
    areal = np.vstack([a.real, a.imag])
    breal = np.concatenate([b.real, b.imag])
    coeff, _, rank, sing = np.linalg.lstsq(areal, breal, rcond=1e-12)
    residual = a @ coeff - b
    max_residual = float(np.max(abs(residual)))
    frob_residual = float(np.linalg.norm(residual)/np.linalg.norm(b))
    weight_projection = sum((c*t[1, 1, :, :] for c, t in zip(coeff, tensors)),
                            np.zeros((dim, dim), dtype=complex))
    vev_residual = float(np.max(abs(weight_projection-np.eye(dim))))
    success = max_residual <= tol and frob_residual <= tol and vev_residual <= tol
    return {
        'full_tensor_identity_in_real_native_span': bool(success),
        'native_invariant_linear_rank': int(rank),
        'native_invariant_count': len(tensors),
        'native_coupling_coefficients_for_unit_identity': [float(v) for v in coeff],
        'max_abs_four_index_residual': max_residual,
        'relative_four_index_frobenius_residual': frob_residual,
        'max_abs_H0_projected_residual': vev_residual,
        'condition_ratio_smallest_over_largest': float(sing[-1]/sing[0]) if sing[0] else 0.,
        'native_vev_projector_diagonals': [
            [float(x.real) for x in np.diag(t[1, 1, :, :])]
            for t in tensors
        ],
        'does_not_certify_coupling_sign_or_normalization': True,
    }


def source_native_tensors(payload: dict, dim: int, which: int, us) -> list[np.ndarray]:
    from Numerical.diagnostics.AuditT3ComponentMultiplicity import expand_tensor
    from RGE.running.intermediate.IntermediateMatcheteParsing import load_cg_registry
    reg = load_cg_registry(payload)
    records = [r for r in payload['CGRegistry'] if int(r['WhichScalar']) == which]
    records.sort(key=lambda row: int(row['Invariant']))
    require([int(r['Invariant']) for r in records] == list(range(1,len(records)+1)),
            f'Invalid portal invariant numbering S{which}')
    require(1 <= len(records) <= 2, f'Unexpected invariant inventory for S{which}')
    result = []
    for k, row in enumerate(records, start=1):
        name = f'T3HiggsPortal{which}CGInv{k}'
        require(row['Name'] == name and
                row['Coupling'] == f'lambdaH{which}Inv{k}' and
                row['SourceInteraction'] == f'HiggsPortal{which}Inv{k}' and
                row['Dimensions'] == [2,2,dim,dim] and name in reg,
                f'Native portal metadata mismatch: {name}')
        components = [{'indices': list(idx), 'cg': str(value)}
                      for idx, value in expand_tensor(reg[name].tensor, (2,2,dim,dim))]
        result.append(_tensor({'dimensions': [2,2,dim,dim],
                               'components': components}, FLAGS, us))
    return result


def run(project: Path, portal_root: Path, audit_dir: Path, output: Path) -> dict:
    source_hashes = source_contract(project)
    rgbeta = rgbeta_identity_source_contract(project)
    covariance = read_json(audit_dir/'cg_covariance_audit.json')
    require(covariance.get('status') == 'cg_generator_covariance_verified'
            and len(covariance.get('models', [])) == 16
            and all(x.get('all_covariant') is True for x in covariance['models']),
            'Previous all-16 CG covariance validation missing')
    us = _rotations(read_json(audit_dir/'matchete_generator_analysis.json'))
    previous = read_json(audit_dir/'native_higgs_portal_basis_audit.json')
    require(previous.get('models_checked') == 16
            and len(previous.get('models', [])) == 16,
            'Prior all-16 portal projection report missing')
    prev_by_model = {x['model']: x for x in previous['models']}
    require(len(prev_by_model) == 16, 'Duplicate prior portal model record')
    results = []
    for short, frozen_name in models():
        frozen_path = project/'Numerical/benchmarks/fitted_reference'/f'{frozen_name}.json'
        frozen = read_json(frozen_path)
        rep = frozen['representation']
        dims = tuple(int(rep[k]) for k in ('d_s1','d_s2','d_f','alpha'))
        expected = (*CLASSES[short[0]][:3], int(short[2:]))
        require(dims == expected, f'Frozen metadata mismatch {short}')
        native_path = portal_root/frozen_name/'native_higgs_portal_registry.json'
        native = read_json(native_path)
        require(native.get('status') == 'Success' and
                native.get('format') == 't3-native-direct-uv-portal-cg-v1' and
                tuple(int(native['metadata'][k]) for k in ('dS1','dS2','dF','alpha'))==dims,
                f'Invalid native portal export {short}')
        old = prev_by_model[short]
        rows=[]
        for which, dim in ((1,dims[0]),(2,dims[1])):
            tensors=source_native_tensors(native,dim,which,us)
            sol=fit_native_identity(tensors)
            old_rows=next(x for x in old['scalars'] if x['scalar'] == f'S{which}')
            old_diag=[x['H0bar_H0_scalar_diagonal'] for x in old_rows['native_matchete_invariants']]
            require(len(old_diag) == len(sol['native_vev_projector_diagonals']) and
                    all(np.allclose(a,b,rtol=0,atol=1e-11)
                        for a,b in zip(old_diag,sol['native_vev_projector_diagonals'])),
                    f'Native tensor data disagree with prior diagonal projection: {short}/S{which}')
            sol['scalar']=f'S{which}'
            sol['dimension']=dim
            sol['native_coupling_names']=[f'lambdaH{which}Inv{k}' for k in range(1,len(tensors)+1)]
            sol['baseline_rgbeta_name']=f'lambdaH{which}'
            rows.append(sol)
            print(f'{short:6s} S{which} d={dim} '
                  f"span={'PASS' if sol['full_tensor_identity_in_real_native_span'] else 'BLOCKED'} "
                  f"c={sol['native_coupling_coefficients_for_unit_identity']} "
                  f"max_res={sol['max_abs_four_index_residual']:.3e}",flush=True)
        results.append({'model':short,'frozen_sha256':hash_file(frozen_path),
                        'native_export_sha256':hash_file(native_path), 'scalars':rows})
    require(len(results)==16, 'Expected all 16 models')
    success=sum(r['full_tensor_identity_in_real_native_span'] for m in results for r in m['scalars'])
    payload={'format':FORMAT,
             'status':'all_32_geometric_identities_representable_physical_dictionary_pending'
                      if success==32 else 'some_full_tensor_identity_decompositions_blocked',
             'models_checked':16,
             'scalar_portal_sectors_checked':32,
             'full_tensor_identities_representable':success,
             'full_tensor_identities_blocked':32-success,
             'matchete_source_sha256':source_hashes,
             'rgbeta_source_sha256':rgbeta,
             'prior_native_projector_report_sha256':hash_file(audit_dir/'native_higgs_portal_basis_audit.json'),
             'method':'real-linear least squares on full Hermitian 2 x 2 x d x d source tensors',
             'target':'delta_ij(H) delta_ab(S): geometric RGBeta lambdaH{1,2} pair invariant',
             'important':'A geometric tensor identity is NOT a physical coupling dictionary; '
                         'relative Lagrangian sign, interaction prefactor, fermion phases, '
                         'and scalar bilinears remain uncertified. No mass matrix is changed.',
             'models':results}
    output.parent.mkdir(parents=True,exist_ok=True)
    value=json.dumps(payload,indent=2,allow_nan=False)+'\n'
    if output.is_file():
        require(output.read_text(encoding='utf-8-sig')==value,
                f'Existing non-identical report not overwritten: {output}')
    else:
        output.write_text(value,encoding='utf-8')
    print(f'Full-index identity sectors: {success}/32; geometric blocks: {32-success}')
    print('Physical portal coupling sign, normalization, and all-16 masses: NOT CERTIFIED')
    print('Report:',output)
    return payload


def main() -> int:
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--project-root',type=Path,default=ROOT)
    parser.add_argument('--portal-root',type=Path,required=True)
    parser.add_argument('--audit-dir',type=Path,required=True)
    parser.add_argument('--output',type=Path,default=None)
    args=parser.parse_args()
    project=args.project_root.resolve()
    audit=args.audit_dir.resolve()
    out=args.output.resolve() if args.output else audit/'native_portal_identity_decomposition.json'
    try:
        run(project,args.portal_root.resolve(),audit,out)
    except (ValueError,KeyError,TypeError,OSError,IndexError,RuntimeError) as e:
        parser.exit(2,'BLOCKED: '+str(e)+'\n')
    return 0


if __name__=='__main__':
    raise SystemExit(main())
