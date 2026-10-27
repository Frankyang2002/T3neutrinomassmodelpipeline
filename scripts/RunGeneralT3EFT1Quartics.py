"""Evolve all Stage-29 candidate scalar quartics in the *renormalisable* T3 EFT1.

Uses the EXISTING Stage-33 threshold boundary, experimental model-dependent
GeneralT3ScalarRGState, and full RGBeta quartic-class projector. Wolfram betas
are regenerated from local source for every invocation: older Stage-20 exports
may predate the Stage-24 B/C SU(2)-covariance corrections.

This is a convention-labelled experimental scalar-RGE calculation, NOT a
Matchete/RGBeta physical normalization certificate, Wilson matching, or a
complete generalized neutrino-mass prediction. Does not modify Stage 34--36.

Run from repository root:
    python -m scripts.RunGeneralT3EFT1Quartics --model B_+1
    python -m scripts.RunGeneralT3EFT1Quartics --all
"""
from __future__ import annotations

import argparse
from hashlib import sha256
import json
import math
from pathlib import Path
import subprocess
from typing import Any

from scripts import Stage33FermionThresholdHandoff as s33
from scripts import Stage34IntermediateEFTRunning as s34
from Numerical.experimental.GeneralT3ScalarRGState import (
    GeneralT3ScalarRGState,
    evaluate_general_scalar_derivative,
    load_registry,
    pack_general_scalar_state,
    registry_model,
    run_general_scalar_segment,
)

ROOT = s34.ROOT
REGISTRY = ROOT / 'RGE' / 'general' / 'generated' / 'general_rgbeta_quartic_registry.json'
RGBETA_RUNNER = ROOT / 'RGE' / 'running' / 'rgbeta' / 'RunGeneralT3RGBeta.wl'
DEFAULT_OUTPUT = ROOT / 'output' / 'general_t3_eft1_quartics'
FORMAT = 'general-t3-renormalisable-eft1-all-quartics-v1'


def fail(message: str) -> None:
    raise ValueError('General T3 EFT1 quartic check: ' + message)


def sha(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def decode_complex(raw: Any, name: str) -> complex:
    if not isinstance(raw, dict) or set(raw) != {'real', 'imag'}:
        fail(f'{name}: expected explicit real/imag scalar')
    try:
        value = complex(float(raw['real']), float(raw['imag']))
    except (ValueError, TypeError) as exc:
        fail(f'{name}: cannot decode complex scalar: {exc}')
    if not (math.isfinite(value.real) and math.isfinite(value.imag)):
        fail(f'{name}: non-finite complex scalar')
    return value


def export_fresh_eft1_betas(model: str, rep: Any, path: Path,
                            *, wolframscript: str) -> None:
    """Require *live* full projected betas; never trust pre-Stage-24 caches."""
    if not RGBETA_RUNNER.is_file():
        fail(f'Cannot find RGBeta general runner: {RGBETA_RUNNER}')
    path.parent.mkdir(parents=True, exist_ok=True)
    path.unlink(missing_ok=True)
    # WolframScript treats an argument beginning with '-' as a CLI option
    # on Windows, so a negative alpha can disappear from $ScriptCommandLine.
    # The repository's T3RGBetaParseIntegerToken explicitly supports m2 -> -2.
    alpha_token = f'm{-rep.alpha}' if rep.alpha < 0 else str(rep.alpha)
    command = [
        wolframscript, '-file', str(RGBETA_RUNNER), 'EFT1',
        str(rep.d_s1), str(rep.d_s2), str(rep.d_f), alpha_token,
        model, str(path),
    ]
    try:
        process = subprocess.run(command, cwd=ROOT, capture_output=True,
                                 text=True, check=False)
    except OSError as exc:
        fail(f'Wolfram process could not start: {exc}')
    if process.returncode != 0 or not path.is_file():
        fail(f'{model}: live RGBeta EFT1 export failed (exit={process.returncode})\n'
             + f'Invocation: {command!r}\n'
             + process.stdout[-3000:] + '\n' + process.stderr[-3000:])


def checked_payload(model: str, state: GeneralT3ScalarRGState,
                    registry: dict[str, Any], path: Path) -> tuple[dict, tuple[str, ...]]:
    """Verify complete projected beta inventory and *no* silent zero betas."""
    payload = s34.read_json(path)
    _, layout = pack_general_scalar_state(state, registry)
    rep = state.baseline.representation
    info = registry_model(registry, rep)
    meta = payload.get('metadata')
    if payload.get('status') != 'Success' or not isinstance(meta, dict):
        fail(f'{model}: unsuccessful RGBeta export')
    expected = {
        'GeneralScalarModel': model,
        'GeneralScalarStage': 'EFT1',
        'GeneralScalarBasisStatus': 'CandidateNativeRGBetaUnmappedToMatchete',
        'GeneralScalarQuarticBetaBasis': 'RGBetaFullQuarticProjection',
        'GeneralScalarQuadraticMassMixing': 'NotIncluded',
        'GeneralScalarAdditionalQuartics': info['total_additional_quartics'],
        'GeneralScalarExistingStructuralCount': len(info['preexisting_rgbeta_symbols']),
        'dS1': rep.d_s1, 'dS2': rep.d_s2,
        'dF': rep.d_f, 'alpha': rep.alpha,
    }
    for key, value in expected.items():
        if meta.get(key) != value:
            fail(f'{model}: RGBeta metadata {key!r} is not {value!r}')
    extras = set(state.extra)
    if set(meta.get('GeneralScalarNewBetaKeys', [])) != extras:
        fail(f'{model}: projected extra-quartic symbols differ from registry')
    betas = payload.get('report_betas')
    if not isinstance(betas, dict) or set(betas) != set(layout.names):
        missing = set(layout.names) - set(betas or {})
        excess = set(betas or {}) - set(layout.names)
        fail(f'{model}: incomplete projected EFT1 betas: missing={sorted(missing)}, '
             f'excess={sorted(excess)}')
    for name, expr in betas.items():
        if not isinstance(expr, str) or not expr.strip() or any(
            bad in expr for bad in ('BetaTerm[', 'Missing[', '$Failed', '$Aborted')
        ):
            fail(f'{model}: unresolved beta expression for {name}')
    # Also checks all numeric expressions and the exact ODE-layout contract.
    derivative, checked_layout = evaluate_general_scalar_derivative(payload, state, registry)
    if checked_layout != layout or len(derivative) != layout.size:
        fail(f'{model}: beta-vector dimension mismatch')
    return payload, tuple(layout.names)


def run_model(model: str, *, folder: Path, outdir: Path,
              wolframscript: str) -> dict[str, Any]:
    registry = load_registry(REGISTRY)
    if model not in s34.MODELS:
        fail(f'Unsupported ordinary T3 model: {model}')
    # Reuse the existing Stage-33 mass/Yukawa reconstruction and hash checks.
    single, baseline, sidecar, reused, input_hashes = s34.verify_boundary(model, folder)
    info = registry_model(registry, baseline.representation)
    expected_extra = {x['name'] for x in info['new_rgbeta_symbols']}
    expected_reused = {x['existing_name'] for x in info['preexisting_rgbeta_symbols']}
    if set(sidecar) != expected_extra or set(reused) != expected_reused:
        fail(f'{model}: Stage-33 sidecar/reused inventory differs from Stage-29 registry')

    initial = GeneralT3ScalarRGState(
        baseline, {name: decode_complex(value, name) for name, value in sidecar.items()}
    ).validated(registry)
    outdir.mkdir(parents=True, exist_ok=True)
    tag = s34.safe_model(model)
    betafile = outdir / f'{tag}_eft1_live_full_rgbeta.json'
    export_fresh_eft1_betas(model, baseline.representation, betafile,
                            wolframscript=wolframscript)
    payload, names = checked_payload(model, initial, registry, betafile)
    evolution = run_general_scalar_segment(
        initial, payload, registry, s34.STOP, rtol=1e-8, atol=1e-11,
    )
    final = evolution['state'].validated(registry)
    if not math.isclose(final.baseline.mu_gev, s34.STOP, rel_tol=1e-12):
        fail(f'{model}: scalar EFT1 endpoint is wrong')
    if set(final.extra) != expected_extra:
        fail(f'{model}: lost an extra quartic during integration')

    output = {
        'format': FORMAT,
        'status': 'CandidateCompleteRenormalisableQuarticEFT1',
        'model': model,
        'mu_initial_gev': s34.START,
        'mu_final_gev': s34.STOP,
        'fermion_threshold_boundary': 'Stage33 tree-level renormalisable projection',
        'full_quartic_class_projected': True,
        'additional_quartics_evolved_in_eft1': True,
        'new_quartic_count': len(expected_extra),
        'reused_quartic_count': len(expected_reused),
        'quadratic_mass_mixings_excluded': info['quadratics_not_implemented'],
        'coupling_names': list(names),
        'solver': 'existing GeneralT3ScalarRGState.run_general_scalar_segment / DOP853',
        'solver_nfev': evolution['nfev'],
        'initial_baseline': s33.json_state(initial.baseline),
        'final_baseline': s33.json_state(final.baseline),
        'extra_quartics_at_mu_f': {
            name: s33.json_state(value) for name, value in initial.extra.items()},
        'extra_quartics_at_mu_s': {
            name: s33.json_state(value) for name, value in final.extra.items()},
        'physical_matchete_rgbeta_coupling_map_certified': False,
        'quadratic_mass_basis_complete': info['quadratics_not_implemented'] == 0,
        'fermion_wilson_matching_performed': False,
        'wilson_eft1_running_performed': False,
        'scalar_threshold_matching_performed': False,
        'general_neutrino_mass_calculation_performed': False,
        'source_sha256': {
            **input_hashes,
            str(REGISTRY): sha(REGISTRY),
            str(betafile): sha(betafile),
            str(RGBETA_RUNNER): sha(RGBETA_RUNNER),
        },
        'limitations': [
            'General scalar quartics are evolved in the candidate RGBeta basis; '
            'cross-package physical normalization is uncertified.',
            'Eight per-model quadratic mass-mixing structures are not included.',
            'No threshold-induced higher-dimensional Wilson matching/running.',
            'No scalar-threshold rematching or new neutrino-mass prediction.',
        ],
    }
    destination = outdir / f'{tag}_eft1_full_quartic_running.json'
    destination.write_text(json.dumps(output, indent=2, allow_nan=False) + '\n',
                           encoding='utf-8')
    return output


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('--model', choices=s34.MODELS)
    group.add_argument('--all', action='store_true')
    parser.add_argument('--stage33-dir', type=Path, default=s34.DATA)
    parser.add_argument('--output-dir', type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument('--wolframscript', default='wolframscript')
    args = parser.parse_args()
    rows = []
    for model in (s34.MODELS if args.all else (args.model,)):
        report = run_model(model, folder=args.stage33_dir, outdir=args.output_dir,
                           wolframscript=args.wolframscript)
        rows.append(report)
        print(f'{model}: {report["new_quartic_count"]} new + '
              f'{report["reused_quartic_count"]} reused quartics in full EFT1; '
              f'nfev={report["solver_nfev"]}; '
              'physical map UNCERTIFIED', flush=True)
    if args.all:
        actual = tuple(sum(row[key] for row in rows) for key in (
            'new_quartic_count', 'reused_quartic_count', 'quadratic_mass_mixings_excluded'))
        if actual != (38, 12, 8):
            fail(f'16-model coverage differs from (38,12,8): {actual}')
        print('16/16 candidate full-renormalisable-quartic EFT1 evolutions; '
              '38 new, 12 reused, 8 quadratic entries excluded.')


if __name__ == '__main__':
    main()
