"""Check the independent A_-4 native Yukawa-CG product against Matchete's LLSS CG.

Read-only input: native_cg_a_m4_*.json from the installed native_cg_a_m4_probe.wl.
Output: a new, non-overwritten JSON in normalization_probes/.

IMPORTANT: This is a *representation-level* test. Native UV Feynman-rule spinor
flow, which was separately tested, is not re-evaluated; the absolute phase of
the master two-component Weyl operator is not established here.
"""
from __future__ import annotations
import argparse
import cmath
from datetime import datetime, timezone
import hashlib
import json
import math
from pathlib import Path
import sys

HERE = Path(__file__).resolve().parent
PROBES = HERE / 'normalization_probes'
EXPECTED_FIELDS = (
    ('native_yukawa1_cg', 'T3Y1CG', [2, 2]),
    ('native_yukawa2_cg', 'T3Y2CG', [2, 2, 3]),
    ('native_matched_mix_cg', 'T3MixCG', [2, 2, 3]),
)


def digest(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def parse_cg(cg: dict, name: str, dims: list[int]) -> list[complex]:
    if not isinstance(cg, dict) or cg.get('name') != name or cg.get('dimensions') != dims:
        raise ValueError(f'{name}: native CG name or dimension mismatch')
    if not isinstance(cg.get('representation_inputform'), str) or not cg['representation_inputform']:
        raise ValueError(f'{name}: missing representation orientation')
    if not isinstance(cg.get('tensor_inputform'), str) or 'SparseArray[' not in cg['tensor_inputform']:
        raise ValueError(f'{name}: missing Matchete native CG tensor source')
    n = math.prod(dims)
    re, im = cg.get('real_flat'), cg.get('imag_flat')
    if not isinstance(re, list) or not isinstance(im, list) or len(re) != n or len(im) != n:
        raise ValueError(f'{name}: invalid flattened CG values')
    if any(not isinstance(x, (int, float)) or isinstance(x, bool) or not math.isfinite(x)
           for x in re + im):
        raise ValueError(f'{name}: non-finite component')
    z = [complex(a, b) for a, b in zip(re, im)]
    if max(abs(x) for x in z) <= 0:
        raise ValueError(f'{name}: identically zero CG tensor')
    return z


def product_yukawa_cgs(y1: list[complex], y2: list[complex], *,
                        bar_y1=True, bar_y2=True) -> list[complex]:
    """Contract the Dirac internal doublet index k.

    Native A_-4 *holomorphic* vertices, from the saved UV interaction:
      Bar[T3Y1CG][Bar[L_i],Bar[F_k]] (L_i^T C PL F_k) S1bar;
      Bar[T3Y2CG][Bar[L_j], F_k, S2_u] (Fbar_k PL L_j) S2.
    Both 'Bar[CG]' conjugate the CG components. Heavy Dirac mass contracts
    fundamental and conjugate doublet with delta[k,k'], not epsilon[k,k'].
    """
    if len(y1) != 4 or len(y2) != 12:
        raise ValueError('Unexpected Yukawa CG dimensions')
    b1 = (lambda z: z.conjugate()) if bar_y1 else (lambda z: z)
    b2 = (lambda z: z.conjugate()) if bar_y2 else (lambda z: z)
    return [sum(b1(y1[2*i + k]) * b2(y2[(2*j + k)*3 + u])
                for k in range(2))
            for i in range(2) for j in range(2) for u in range(3)]


def proportionality(candidate: list[complex], native: list[complex]) -> dict:
    if len(candidate) != len(native) or not candidate:
        raise ValueError('Tensor length mismatch')
    denom = sum(abs(x)**2 for x in native)
    if denom == 0:
        raise ValueError('Zero target CG tensor')
    scale = sum(t.conjugate() * x for x, t in zip(candidate, native)) / denom
    diff = [x - scale * t for x, t in zip(candidate, native)]
    cnorm = math.sqrt(sum(abs(x)**2 for x in candidate))
    tnorm = math.sqrt(denom)
    residual = math.sqrt(sum(abs(x)**2 for x in diff)) / max(cnorm, 1e-300)
    return dict(coefficient_real=scale.real, coefficient_imag=scale.imag,
                residual_relative=residual, candidate_norm=cnorm,
                target_norm=tnorm, proportional=bool(residual <= 1e-10),
                unity_ratio=bool(residual <= 1e-10 and abs(scale-1) <= 1e-10),
                minus_unity_ratio=bool(residual <= 1e-10 and abs(scale+1) <= 1e-10))


def analyze(data: dict) -> dict:
    if data.get('format') != 't3-native-a-m4-yukawa-mixing-cg-evidence-v1' or data.get('model') != 'A_-4':
        raise ValueError('Wrong source format/model')
    if data.get('status') != 'NativeCGTensorsExtracted_NotPhaseCertified':
        raise ValueError('Source CG extraction not successful')
    if data.get('previous_release_changed') is not False or data.get('cg_product_identity_computed') is not False:
        raise ValueError('Unexpected native source certification flags')
    y1, y2, mix = [parse_cg(data[k], name, dims)
                   for k, name, dims in EXPECTED_FIELDS]
    oriented = product_yukawa_cgs(y1, y2)
    direct = proportionality(oriented, mix)
    swapped = proportionality(oriented, [mix[(j*2+i)*3+u]
                                          for i in range(2) for j in range(2) for u in range(3)])
    # Additional diagnostic only, never substitute a more favourable convention
    # for the native field orientation without source justification.
    unbarred = proportionality(product_yukawa_cgs(y1,y2,bar_y1=False,bar_y2=False),mix)
    return {
        'format': 't3-a-m4-native-cg-product-audit-v1',
        'model': 'A_-4',
        'status': ('NativeOrientedCGProductMatchesMixCG' if direct['unity_ratio']
                   else 'OrientedCGMismatchOrNonunitFactor_ReviewBeforeCertification'),
        'orientation_assumption': 'Bar[T3Y1CG](i,k)*Bar[T3Y2CG](j,k,u), summed over heavy Dirac SU2 doublet k; compare T3MixCG(i,j,u)',
        'direct_oriented': direct,
        'swapped_mix_order_diagnostic_only': swapped,
        'unbarred_cg_diagnostic_only': unbarred,
        'native_cg_identity_unit_ratio_verified': direct['unity_ratio'],
        'native_cg_identity_proportional_verified': direct['proportional'],
        'full_uv_feynman_amplitude_and_master_weyl_phase_certified': False,
        'master_weyl_phase_certified': False,
        'sixteen_model_boundary_changed': False,
        'native_cg_representations': {k:data[k]['representation_inputform'] for k,_,_ in EXPECTED_FIELDS},
        'native_cg_sources_sha256': {k: digest(data[k]['tensor_inputform'].encode())
                                      for k,_,_ in EXPECTED_FIELDS},
    }


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--report', type=Path,
                    help='Exact path to native_cg_a_m4_<uuid>.json')
    ap.add_argument('--check-only', action='store_true')
    args = ap.parse_args(argv)
    try:
        path = args.report
        if path is None:
            matches = list(PROBES.glob('native_cg_a_m4_*.json'))
            if len(matches) != 1:
                raise ValueError(f'Expected exactly one native CG source, found {len(matches)}; specify --report')
            path = matches[0]
        raw = path.read_bytes()
        data = json.loads(raw.decode('utf-8-sig'))
        result = analyze(data)
        result['native_cg_report_path'] = str(path)
        result['native_cg_report_sha256'] = digest(raw)
        print('A_-4 native CG product:', result['status'])
        for name in ('direct_oriented','swapped_mix_order_diagnostic_only','unbarred_cg_diagnostic_only'):
            v = result[name]
            print(f"  {name}: ratio=({v['coefficient_real']:.16g},{v['coefficient_imag']:.16g}) "
                  f"relative residual={v['residual_relative']:.3e}")
        if args.check_only:
            print('CHECK ONLY: no file written')
        else:
            result_dir = PROBES
            result_dir.mkdir(parents=True, exist_ok=True)
            short = digest(raw)[:16]
            out = result_dir / f'native_cg_product_a_m4_{short}.json'
            # Do not overwrite certified or previous probe results.
            with out.open('x', encoding='utf8') as f:
                json.dump(result, f, indent=2, ensure_ascii=False)
                f.write('\n')
            print('CREATED:', out)
        if not result['native_cg_identity_unit_ratio_verified']:
            print('NOT CERTIFIED: oriented native CG not identically +1 times the matched CG; inspect evidence')
            return 2
        print('PASS: oriented native CG product reproduces matched mixing CG')
        print('OPEN: absolute Matchete GammaCC -> master Weyl phase')
        return 0
    except (OSError, ValueError, KeyError, TypeError, AssertionError) as exc:
        print(f'BLOCKED: {type(exc).__name__}: {exc}', file=sys.stderr)
        return 2


if __name__ == '__main__':
    raise SystemExit(main())
