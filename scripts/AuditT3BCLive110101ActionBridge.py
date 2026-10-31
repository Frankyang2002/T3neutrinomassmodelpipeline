"""B/C(-1) 110101: exact cross-tool quartic action identity, gated by live evidence.

Inputs are (i) verified live Matchete CheckLagrangian expressions,
(ii) verified live RGBeta Lam 64/model coefficients, and (iii) independent
native Stage22 Matchete/RGBeta component tensors.  Proves the signed
polynomial identity IN THE EXPORTED SU(2) COMPONENT BASIS.  Does not evaluate
Matchete's live CG numerically here, check off-diagonal heavy masses, or Match[].

Run from project root:
 python -u -X utf8 -m scripts.AuditT3BCLive110101ActionBridge \\
  --run-id bc_live_action_bridge_v1
"""
from __future__ import annotations

import argparse
from hashlib import sha256
import json
from pathlib import Path
import re

import sympy as sp
import numpy as np

from scripts.AuditT3BC110101MassBasisMap import (
    MODELS, MAT, RGB, NATIVE, RGB_NATIVE, RGB_LIVE, ROOT,
    check_live, exact_map, exact_number, load, original_matchete_records,
)
from scripts.VerifyT3BCNativeLam110101Coefficients import verify as verify_lam
from scripts.VerifyT3BCMatchete110101Action import (
    verify as verify_matchete, verify_source as verify_matchete_source,
)

MAT_LIVE = ROOT / 'output/audits/t3_bc_matchete_action/matchete_110101_v2.json'
RGB_LAM_LIVE = ROOT / 'output/audits/t3_bc_native_lam/native_lam_coeff_v4.json'
MAT_SOURCE = ROOT / 'Lagrangian/interactions/GeneralScalarOperators.wl'
OUT = ROOT / 'output/audits/t3_bc_110101_live_action_bridge'
EXPECTED_A = sp.Matrix([[sp.Rational(1, 2), 1], [-sp.sqrt(3)/2, 0]])


def require(cond: bool, explanation: str) -> None:
    if not cond:
        raise ValueError('B/C 110101 live-evidence action bridge: ' + explanation)


def digest(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def exact_action_identity(m: sp.Matrix, r: sp.Matrix, A: sp.Matrix) -> dict:
    """Exact direct and h.c. coefficient identity; reject wrong sign/factor.

    Columns of m,r are the two independent native contractions; each has
    16 ordered SU(2) components in H,Hdag,S1dag,S2dag order.
    R = M*A; L_M = -V_R means h=-A*g, including h.c.
    """
    require(m.shape == (16, 2) and r.shape == (16, 2),
            'expected 16 components, two contractions in both bases')
    require(A == EXPECTED_A and sp.simplify(A.det()) != 0,
            'unexpected exact two-invariant basis map')
    require(m.rank() == 2 and r.rank() == 2,
            'both sets of 110101 tensors must be linearly independent')
    signed = -A
    identity_direct = m*signed + r
    identity_conjugate = m.conjugate()*signed.conjugate() + r.conjugate()
    direct_residuals = [sp.simplify(x) for x in identity_direct]
    hc_residuals = [sp.simplify(x) for x in identity_conjugate]
    require(all(x == 0 for x in direct_residuals + hc_residuals),
            'exact signed action identity fails for direct or h.c. components')
    # Nontrivial negative controls: changing sign, retaining spurious /4,
    # or swapping the two Matchete contractions must fail the full action.
    controls = {
        'wrong_potential_sign': m*A + r,
        'spurious_extra_quarter': m*(signed/4) + r,
        'swapped_matchete_contractions': m*sp.Matrix([[0, 1], [1, 0]])*signed + r,
    }
    for name, matrix in controls.items():
        require(any(sp.simplify(x) != 0 for x in matrix),
                'negative control unexpectedly passed: ' + name)
    # A formal complex interaction has independent holomorphic coefficients
    # and their conjugates; the equalities above cover both orientations.
    return {
        'signed_matchete_lagrangian_from_rgbeta_potential': [
            [str(signed[i, j]) for j in range(2)] for i in range(2)],
        'rgbeta_in_matchete_invariant_basis': [
            [str(A[i, j]) for j in range(2)] for i in range(2)],
        'direct_component_coefficients_compared': 32,
        'hermitian_conjugate_component_coefficients_compared': 32,
        'exact_nonzero_component_residuals': 0,
        'exact_action_identity_in_exported_su2_component_basis': True,
        'incorrect_sign_rejected': True,
        'spurious_extra_quarter_rejected': True,
        'swapped_contractions_rejected': True,
    }


def run_one(model: str, native: dict, invariants: dict, covariance: dict) -> dict:
    require(model in MODELS, 'not B/C(-1)')
    rg_tensors = check_live(covariance, invariants, model)
    A = exact_map(native, model, rg_tensors)
    mats = original_matchete_records(native, model)
    m = sp.Matrix.hstack(*(sp.Matrix([exact_number(s) for s in record]) for record in mats))
    r = sp.Matrix.hstack(*(sp.Matrix([sp.Integer(int(round(float(np.real(v))))) for v in t.flat])
                           for t in rg_tensors))
    checks = exact_action_identity(m, r, A)
    return {
        'model': model,
        'status': 'LiveEvidenceGatedExactExportedBasisActionIdentity',
        'matchete_couplings': list(MAT),
        'rgbeta_couplings': list(RGB),
        'matchete_live_individual_checks_passed': 2,
        'matchete_live_combined_check_passed': True,
        'rgbeta_live_lam_component_tests_passed': 64,
        'matchete_live_cg_numerical_components_independently_extracted': False,
        'direct_live_cross_tool_action_subtraction_evaluated': False,
        'physical_coupling_mapping_certified': False,
        'mass_eigenstate_matchete_matching_certified': False,
        'sequential_eft_matching_certified': False,
        **checks,
    }


def run(mat_live_path: Path, rgb_live_path: Path, native_path: Path,
        invariants_path: Path, covariance_path: Path, source_path: Path) -> dict:
    sources = [mat_live_path, rgb_live_path, native_path, invariants_path,
               covariance_path, source_path]
    for path in sources:
        require(path.is_file(), 'missing source ' + str(path))
    verify_matchete_source(source_path)
    verified_m = verify_matchete(load(mat_live_path))
    verified_r = verify_lam(load(rgb_live_path))
    require(verified_m['matchete_native_terms_validated'] == 4 and
            verified_m['matchete_combined_checks_passed'] == 2 and
            verified_r['total_components_checked'] == 128,
            'live prerequisite reports incomplete')
    native, invariants, covariance = (load(x) for x in sources[2:5])
    models = [run_one(model, native, invariants, covariance) for model in MODELS]
    return {
        'format': 't3-bc-live-evidence-gated-110101-action-identity-v1',
        'status': 'ExactExportedComponentActionIdentityWithSeparateLiveGates',
        'source_sha256': {str(path.resolve()): digest(path) for path in sources},
        'models': models,
        'total_direct_and_hc_coefficients_checked': 128,
        'live_matchete_interactions_checked': 4,
        'live_rgbeta_lam_coefficients_checked': 128,
        'direct_live_matchete_cg_component_equality_verified': False,
        'direct_cross_tool_action_subtraction_in_live_matchete_verified': False,
        'physical_coupling_mapping_certified': False,
        'full_mass_basis_matchete_matching_certified': False,
        'sequential_eft_matching_certified': False,
        'production_modified': False,
    }


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--matchete-live', type=Path, default=MAT_LIVE)
    p.add_argument('--rgbeta-live', type=Path, default=RGB_LAM_LIVE)
    p.add_argument('--matchete-native', type=Path, default=NATIVE)
    p.add_argument('--rgbeta-invariants', type=Path, default=RGB_NATIVE)
    p.add_argument('--rgbeta-covariance', type=Path, default=RGB_LIVE)
    p.add_argument('--matchete-source', type=Path, default=MAT_SOURCE)
    p.add_argument('--run-id', required=True)
    p.add_argument('--output-root', type=Path, default=OUT)
    args = p.parse_args(argv)
    require(re.fullmatch(r'[A-Za-z0-9][A-Za-z0-9_-]{0,63}',args.run_id) is not None,
            'invalid run ID')
    output=args.output_root/args.run_id/'summary.json'
    require(not output.exists(), 'refusing to overwrite ' + str(output))
    result=run(args.matchete_live,args.rgbeta_live,args.matchete_native,
               args.rgbeta_invariants,args.rgbeta_covariance,args.matchete_source)
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n',encoding='utf-8')
    for row in result['models']:
        print(f"{row['model']}: exact exported-basis action coefficients 64/64; "
              'sign and factor controls rejected',flush=True)
    print('Separate live Matchete and RGBeta gates passed; direct live CG/action equality NOT certified')
    print('Summary:',output)
    return 0

if __name__ == '__main__':
    raise SystemExit(main())
