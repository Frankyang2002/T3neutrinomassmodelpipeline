"""Experimental B_-1/C_-1 Matchete heavy-eigenfield validation.

Consumes the existing *final or initial* EFT1 mixed-doublet rotation JSON
and supplies a numerical U and positive physical masses to a new Matchete
UV interaction-basis validation. The original 16-model production Lagrangian
and threshold matchers remain untouched.

python -m scripts.RunGeneralT3MassEigenMatchete --all
python -m scripts.RunGeneralT3MassEigenMatchete --all --prepare-only
"""
from __future__ import annotations

import argparse
from hashlib import sha256
import json
import math
import os
from pathlib import Path
import shutil
import subprocess
from typing import Any, Mapping

import numpy as np

from Numerical.experimental.GeneralT3MixedScalarMassBasis import (
    EPS, MIXED_MODELS, MixedDoubletMassBasis,
)

ROOT = Path(__file__).resolve().parents[1]
INPUT = ROOT / 'output/general_t3_mass_basis'
OUTPUT = ROOT / 'output/general_t3_mass_eigen_matchete'
WOLFRAM = ROOT / 'Lagrangian/RunGeneralT3MassEigenMatchete.wl'
FORMAT = 't3-experimental-mass-eigen-matchete-parameters-v1'
STATUS_PREPARED = 'NumericMassEigenstateBoundaryPreparedMatcheteNotRun'
STATUS_PASSED = 'MatcheteMassEigenHeavyBasisPassedPhysicalMatchingNotPerformed'
MODEL_INFO = {'B_-1': (2, 2, 1, -1), 'C_-1': (2, 2, 3, -1)}


def fail(message: str) -> None:
    raise ValueError('T3 Matchete mass-eigenstate fail-closed: ' + message)


def _read_json(path: Path) -> dict:
    if not path.is_file():
        fail(f'Missing {path}; first run scripts.RunGeneralT3MixedScalarMassBasis --all')
    obj = json.loads(path.read_text(encoding='utf-8'))
    if not isinstance(obj, dict):
        fail(f'Not a JSON object: {path}')
    return obj


def _complex_matrix(value: Any, label: str, n: int) -> np.ndarray:
    if not isinstance(value, dict) or set(value) != {'real', 'imag'}:
        fail(label + ': expected real/imag matrix')
    try:
        a = np.asarray(value['real'], dtype=float)
        b = np.asarray(value['imag'], dtype=float)
    except (ValueError, TypeError) as exc:
        fail(label + ': malformed matrix: ' + str(exc))
    if a.shape != (n, n) or b.shape != (n, n) or not (np.isfinite(a).all() and np.isfinite(b).all()):
        fail(label + f': expected finite {n}x{n} complex matrix')
    return a + 1j*b


def _encode_matrix(a: np.ndarray) -> dict:
    return {'real': np.asarray(a.real).tolist(), 'imag': np.asarray(a.imag).tolist()}


def _require_close(a: Any, b: Any, msg: str, *, rtol: float=3e-10, atol: float=1e-8) -> None:
    if not np.allclose(a, b, rtol=rtol, atol=atol):
        fail(msg)


def _normalised_source_hashes(value: Any) -> dict[str, str]:
    """Read provenance created on Windows or POSIX without weakening SHA checks.

    The numerical rotation originally serialized ``str(Path.relative_to(ROOT))``;
    Windows therefore produced ``RGE\\general\\...`` keys, whereas the
    consumer used POSIX-style literal keys. Reject collisions so different
    spellings cannot silently shadow distinct source records.
    """
    if not isinstance(value, dict):
        fail('Mass-basis source hashes missing')
    result: dict[str, str] = {}
    for raw_key, digest in value.items():
        if not isinstance(raw_key, str):
            fail('Non-string mass-basis provenance path')
        key = raw_key.replace('\\', '/')
        parts = key.split('/')
        if (not key or key.startswith('/') or ':' in key or
                any(part in ('', '.', '..') for part in parts)):
            fail('Malformed mass-basis provenance path: ' + repr(raw_key))
        if key in result:
            fail('Duplicate cross-platform provenance path: ' + key)
        if (not isinstance(digest, str) or len(digest) != 64 or
                any(ch not in '0123456789abcdef' for ch in digest)):
            fail('Malformed SHA-256 provenance digest: ' + key)
        result[key] = digest
    return result


def _model_file(model: str, which: str) -> str:
    if model not in MODEL_INFO or which not in ('final', 'initial'):
        fail(f'Unsupported model/endpoint {model}, {which}')
    return model.replace('+','p').replace('-','m') + '_' + which + '_mass_basis.json'


def prepare_model(model: str, *, which: str='final', folder: Path=INPUT,
                  outdir: Path=OUTPUT) -> dict:
    source = folder / _model_file(model, which)
    original = _read_json(source)
    if (original.get('format') != 't3-mixed-doublet-mass-basis-transport-v1' or
        original.get('status') != 'ExperimentalMassBasisTensorTransportNoPhysicalMatching' or
        original.get('model') != model or original.get('endpoint') != which or
        original.get('source_quartic_count') != 17 or
        original.get('source_new_quartics') != 1 or
        original.get('source_reused_quartics') != 6 or
        original.get('MatcheteRGBetaPhysicalCouplingMapVerified') is not False or
        original.get('physical_mass_matching_performed') is not False):
        fail('Unexpected numerical transport schema/certification: ' + model)
    quartics = original.get('quartics')
    expected_quartics = {
        'lambdaH','lambdaS1','lambdaS2','lambdaH1','lambdaH2','lambda12',
        'lambdaT3','lambdaH1Adj','lambdaH2Adj','lambda12Adj',
        'lambdaHHdagS2S2','lambdaHHdagS1barS1bar','lambdaS1bar2S2bar2',
        'lambdaS1barS2S2bar2','lambdaS1S1bar2S2bar',
        'lambdaHHdagS1barS2barCross','lambdaT3GeneralQ110101I1',
    }
    if not isinstance(quartics, dict) or set(quartics) != expected_quartics:
        fail('Mass-basis source quartic dictionary incomplete or changed')
    yukawa = original.get('yukawa_scalar_leg_transforms')
    if not isinstance(yukawa, dict) or set(yukawa) != {'y1','y2'}:
        fail('Both Yukawa scalar-leg transports are required')
    b = original.get('basis')
    if not isinstance(b, dict) or b.get('matchete_field_redefinition_installed') is not False:
        fail('Input mass basis certification unexpectedly changed')
    U = _complex_matrix(b.get('unitary_flavour_rotation_U'), 'unitary U', 2)
    A = _complex_matrix(b.get('mass_squared_flavour_matrix_gev2'), 'mass matrix A', 2)
    values = np.asarray(b.get('masses_squared_gev2'),dtype=float)
    masses = np.asarray(b.get('masses_gev'),dtype=float)
    if values.shape != (2,) or masses.shape != (2,) or not np.isfinite(values).all() or not np.isfinite(masses).all():
        fail('Missing finite two-eigenstate scalar mass data')
    if np.any(values<=0) or np.any(masses<=0):
        fail('Cannot introduce heavy scalar fields with nonpositive benchmark squared masses')
    if abs(A[0,0].imag)>1e-8 or abs(A[1,1].imag)>1e-8:
        fail('Original squared-mass diagonals not real')
    computed = MixedDoubletMassBasis.build(model,float(A[0,0].real),float(A[1,1].real),complex(A[0,1]))
    _require_close(A,A.conj().T,'Squared-mass matrix non-Hermitian')
    _require_close(U.conj().T@U, np.eye(2), 'U is nonunitary',atol=4e-12,rtol=4e-12)
    _require_close(U.conj().T@A@U, np.diag(values), 'U does not diagonalize saved mass matrix')
    _require_close(values, computed.eigenvalues_sq, 'Eigenvalues do not match independent diagonalization')
    _require_close(masses**2, values, 'Masses not square roots of eigenvalues')
    err=float(original.get('max_source_tensor_transport_relative_error',float('nan')))
    if not math.isfinite(err) or err>3e-12:
        fail('Input source tensor transport identity not verified')
    if original.get('quadratic_hessian_agreement') is not True:
        fail('Input eight-real-mode quadratic Hessian not verified')
    scale=float(original.get('scale_gev',float('nan')))
    if not math.isfinite(scale) or scale<=0:
        fail('Invalid EFT1 scale')
    # The 2 heavy complex doublets carry Y=-1/2 after the epsilon map.
    # The fermion is the UV T3 fermion. This is a mathematical UV
    # interaction-basis check, not a claim that EFT1 contains active Yukawas.
    source_sha = _normalised_source_hashes(original.get('source_sha256'))
    for relative in ('RGE/general/generated/general_rgbeta_quartic_registry.json',
                     'docs/general_scalar_trial_operator_manifest.json'):
        path = ROOT / relative
        if not path.is_file() or source_sha.get(relative) != sha256(path.read_bytes()).hexdigest():
            fail('Upstream source manifest changed since numerical rotation: ' + relative)
    payload = {
        'format': FORMAT, 'status': STATUS_PREPARED, 'model':model,
        'dimensions_alpha': list(MODEL_INFO[model]),
        'endpoint':which, 'source_eft1_scale_gev':scale,
        'unitary_U':_encode_matrix(U), 'original_mass_squared_A_gev2':_encode_matrix(A),
        'heavy_mass_squared_gev2':values.tolist(),
        'heavy_masses_gev':masses.tolist(),
        'threshold_plan': {
            'decreasing_scale_order':['NewScalarY','NewScalarX'],
            'mass_scales_gev':[float(masses[1]),float(masses[0])],
            'scale_ratio':float(masses[1]/masses[0]),
            'log_separation':float(np.log(masses[1]/masses[0])),
            'existing_eft1_endpoint_below_both_mass_thresholds':bool(scale < masses[0]),
            'numerical_eft1_state_must_be_evaluated_at_each_threshold':True,
            'sequential_matching_performed':False,
        },
        'epsilon':EPS.real.tolist(),
        'heavy_scalar_fields': ['NewScalarX','NewScalarY'],
        'heavy_scalar_hypercharge':-0.5,
        'quartic_source_count':17,'yukawa_source_count':2,
        'verified_eight_real_mode_hessian':True,
        'verified_17_tensor_rotation':True,
        'source_sha256':sha256(source.read_bytes()).hexdigest(),
        'source_path':source.relative_to(ROOT).as_posix() if source.is_relative_to(ROOT) else str(source),
        'matchete_mass_basis_verified':False,
        'physical_matching_performed':False,
        'matchete_rgbeta_component_map_certified':False,
        'coefficient_map_certified':False,
        'matching_next_actions': [
            'Certify native Matchete epsilon/CG tensor phases against numerical spherical basis',
            'Transform exact source Matchete quartic tensor coefficients to eigenbasis',
            'Match both non-degenerate scalar eigenfields and Wilson coefficients',
        ],
    }
    dest = outdir / (model.replace('+','p').replace('-','m') + '_' + which + '_parameters.json')
    dest.parent.mkdir(parents=True,exist_ok=True)
    dest.write_text(json.dumps(payload, indent=2, allow_nan=False)+'\n',encoding='utf-8')
    return {'model':model,'parameters':dest, 'input':source,
            'eigenvalues':values,'masses':masses,'payload':payload}


def validate_wolfram_output(payload: Mapping[str, Any], expected_model: str,
                            expected_unitary: Any = None) -> None:
    if not isinstance(payload, Mapping) or payload.get('status') != STATUS_PASSED:
        fail(f'Matchete did not pass heavy eigenbasis for {expected_model}: '
             + str(payload.get('status') if isinstance(payload,Mapping) else payload))
    if (payload.get('model') != expected_model or
        payload.get('mass_basis_ordinary_check_passed') is not True or
        payload.get('heavy_mass_basis_detailed_check_passed') is not True):
        fail('Matchete result lacks actual heavy mass-basis validation')
    if payload.get('all_old_scalar_fields_removed') is not True or payload.get('all_17_quartics_present') is not True:
        fail('Matchete result lacks full interaction inventory')
    if payload.get('both_yukawa_interactions_present') is not True:
        fail('Matchete result omits Yukawa source vertices')
    native = payload.get('native_topology_epsilon_tensor')
    if (not isinstance(native, Mapping) or
        native.get('OriginalTimesEpsilonTensorExact') is not True or
        native.get('SymmetricMassScalarCG') is not True or
        native.get('OriginalHHScalarContractions') != 1 or
        native.get('NativeMassScalarContractions') != 1 or
        native.get('MappedRGBetaQuarticNormalizationCertified') is not False or
        native.get('InvariantCG') != 'T3MENativeHHMassCG'):
        fail('Matchete topology tensor rotation lacks an exact native invariant certificate')
    portals = payload.get('native_H2_portal_epsilon_tensors')
    if (not isinstance(portals, Mapping) or
        portals.get('OldInvariantCount') != 2 or
        portals.get('NewInvariantCount') != 2 or
        portals.get('TransportedRank') != 2 or
        portals.get('CombinedBasisRank') != 2 or
        portals.get('ExactTensorProjectionVerified') is not True or
        portals.get('PhysicalCouplingMapCertified') is not False):
        fail('Matchete H2 portal tensor transport lacks an exact native invariant certificate')
    for label, count in (('native_cross_scalar_epsilon_tensors', 2),
                         ('native_S2_self_epsilon_tensor', 1)):
        cert = payload.get(label)
        if (not isinstance(cert, Mapping) or
            cert.get('OldInvariantCount') != count or
            cert.get('NewInvariantCount') != count or
            cert.get('TransportedRank') != count or
            cert.get('CombinedBasisRank') != count or
            cert.get('ExactTensorProjectionVerified') is not True or
            cert.get('PhysicalCouplingMapCertified') is not False):
            fail('Matchete S2 scalar-quartic tensor transport certificate invalid: ' + label)
    native_y2 = payload.get('native_yukawa2_epsilon_tensor')
    if (not isinstance(native_y2, Mapping) or
        native_y2.get('OldInvariantCount') != 1 or
        native_y2.get('NewInvariantCount') != 1 or
        native_y2.get('ExactTensorProjectionVerified') is not True or
        native_y2.get('FermionSU2Dimension') != MODEL_INFO[expected_model][2] or
        native_y2.get('NativeCG') != 'T3MENativeY2CG' or
        native_y2.get('PhysicalCouplingMapCertified') is not False or
        not isinstance(native_y2.get('OldTimesEpsilonCGFactor'), str) or
        not native_y2['OldTimesEpsilonCGFactor']):
        fail('Matchete Yukawa2 tensor transport certificate invalid')
    native_y1 = payload.get('native_yukawa1_tensor')
    if (not isinstance(native_y1, Mapping) or
        native_y1.get('OldInvariantCount') != 1 or
        native_y1.get('NewInvariantCount') != 1 or
        native_y1.get('ExactTensorProjectionVerified') is not True or
        native_y1.get('FermionSU2Dimension') != MODEL_INFO[expected_model][2] or
        native_y1.get('NativeCG') != 'T3MENativeY1CG' or
        native_y1.get('PhysicalCouplingMapCertified') is not False or
        not isinstance(native_y1.get('OldToNewCGFactor'), str) or
        not native_y1['OldToNewCGFactor']):
        fail('Matchete Yukawa1 tensor certificate invalid')
    extra = payload.get('native_additional_quartic_tensors')
    if (not isinstance(extra, Mapping) or
        extra.get('SourceOperatorCount') != 7 or
        extra.get('RecoupledOperatorCount') != 7 or
        extra.get('AllProjectedTensorsInvariant') is not True or
        extra.get('PhysicalCouplingMapCertified') is not False):
        fail('Matchete native additional-quartic inventory certificate invalid')
    certs = extra.get('Certificates')
    expected_extra = {
        'T3ExtraQ_0_0_0_1_1_2_1': (0,0,0,1,1,2,1),
        'T3ExtraQ_0_0_0_2_0_2_1': (0,0,0,2,0,2,1),
        'T3ExtraQ_0_0_1_2_0_1_1': (0,0,1,2,0,1,1),
        'T3ExtraQ_0_2_0_0_2_0_1': (0,2,0,0,2,0,1),
        'T3ExtraQ_0_2_0_2_0_0_1': (0,2,0,2,0,0,1),
        'T3ExtraQ_1_1_0_1_0_1_1': (1,1,0,1,0,1,1),
        'T3ExtraQ_1_1_0_1_0_1_2': (1,1,0,1,0,1,2),
    }
    if not isinstance(certs, Mapping) or set(certs) != set(expected_extra):
        fail('Matchete native additional-quartic operator list incomplete')
    for label, info in expected_extra.items():
        cert = certs[label]
        counts = info[:6]
        tensor_index = info[6]
        digits = ''.join(map(str, counts))
        expected_cg = f'T3MENativeExtraQ{digits}I{tensor_index}CG'
        expected_coupling = f'lambdaT3ExtraQ{digits}I{tensor_index}'
        if (not isinstance(cert, Mapping) or
            cert.get('SourceCoupling') != expected_coupling or
            cert.get('SourceCounts') != list(counts) or
            cert.get('SourceTensorIndex') != tensor_index or
            cert.get('NativeCG') != expected_cg or
            cert.get('S2EpsilonProjectedLegs') != counts[4] + counts[5] or
            cert.get('ExactTensorProjectionVerified') is not True or
            cert.get('NativeInvariantSpanRank') != 2 or
            cert.get('CombinedBasisRank') != 2):
            fail('Matchete additional-quartic tensor certificate invalid: ' + label)
    if payload.get('matchete_rotation_uses_exact_symbolic_spurions') is not True:
        fail('Matchete mass-basis result did not use symbolic mixing couplings')
    numerical_map = payload.get('rotation_spurion_values_from_numerical_U')
    expected_labels = ('T3MEU11', 'T3MEU12', 'T3MEU21', 'T3MEU22')
    if not isinstance(numerical_map, Mapping) or set(numerical_map) != set(expected_labels):
        fail('Matchete symbolic mixing couplings lack a complete numerical U map')
    try:
        mapped = np.asarray([
            [complex(float(numerical_map[f'T3MEU{i}{j}']['real']),
                     float(numerical_map[f'T3MEU{i}{j}']['imag']))
             for j in (1, 2)]
            for i in (1, 2)
        ], dtype=np.complex128)
    except (KeyError, TypeError, ValueError, OverflowError) as exc:
        fail('Malformed symbolic mixing U map: ' + str(exc))
    if not (np.isfinite(mapped.real).all() and np.isfinite(mapped.imag).all()):
        fail('Nonfinite numerical map for symbolic mixing couplings')
    if expected_unitary is not None:
        expected = _complex_matrix(expected_unitary, 'prepared numeric unitary U', 2)
        if not np.allclose(mapped, expected, atol=0.0, rtol=0.0):
            fail('Wolfram symbolic mixing map does not match prepared numerical U')
    for key in ('epsilon_bridge_BB','epsilon_bridge_UU'):
        bridge=payload.get(key)
        if not isinstance(bridge,Mapping) or bridge.get('phase_relative_to_matchete_component_epsilon_verified') is not True:
            fail('Missing Matchete-native epsilon phase alignment: '+key)
    if payload.get('physical_matching_performed') is not False or payload.get('component_basis_map_certified') is not False:
        fail('Matchete output claims an unsupported matching or normalization certification')


def run_model(model: str, *, which: str='final', folder:Path=INPUT,
              outdir:Path=OUTPUT, wolframscript:str='wolframscript',
              prepare_only:bool=False, timeout_sec:int=1800) -> dict:
    result = prepare_model(model,which=which,folder=folder,outdir=outdir)
    if prepare_only:
        return dict(result,status=STATUS_PREPARED)
    if not shutil.which(wolframscript) and not Path(wolframscript).is_file():
        fail('wolframscript unavailable; run --prepare-only or install Wolfram command line')
    if not WOLFRAM.is_file():
        fail('Missing project file ' + str(WOLFRAM))
    name=model.replace('+','p').replace('-','m')+'_'+which
    target=outdir/(name+'_matchete.json')
    log=outdir/(name+'_matchete.log')
    progress = target.with_name(target.stem + '.progress.json')
    # Never report a checkpoint left by a previous run as the current result.
    for stale in (target, progress):
        if stale.exists():
            stale.unlink()
    # Do not inherit experimental general-quartic opt-in: duplicates the
    # explicitly registered seven source quartics in the fresh builder.
    env=dict(os.environ)
    env.pop('T3_EXPERIMENTAL_GENERAL_QUARTICS', None)
    command=[wolframscript,'-file',str(WOLFRAM),str(result['parameters']),str(target)]
    proc=subprocess.run(command,cwd=ROOT,env=env,capture_output=True,text=True,
                        timeout=timeout_sec,check=False)
    log.write_text('command: '+' '.join(command)+'\nreturncode: '+str(proc.returncode)
                   +'\nstdout:\n'+proc.stdout+'\nstderr:\n'+proc.stderr,encoding='utf-8')
    progress_report = None
    if progress.is_file():
        try:
            progress_report = _read_json(progress)
        except (ValueError, OSError, json.JSONDecodeError):
            progress_report = {'last_stage': 'unreadable-progress-json'}
    last_stage = (progress_report or {}).get('last_stage', 'no-stage-checkpoint')
    log_tail = (proc.stdout + '\n' + proc.stderr)[-4200:]
    if not target.is_file():
        fail(f'{model}: Wolfram produced no JSON result (exit={proc.returncode}); '
             f'last stage={last_stage}; see {log} and {progress}\n'
             + log_tail)
    report=_read_json(target)
    if proc.returncode:
        fail(f'{model}: Matchete check failed (exit={proc.returncode}); '
             f'last stage={report.get("last_stage",last_stage)}; '
             f'see {log} and {progress}\n'
             + str(report.get('status')) + ': ' + str(report.get('diagnostic',''))[-1100:]
             + '\nLive Wolfram output tail:\n' + log_tail)
    validate_wolfram_output(report,model,result['payload']['unitary_U'])
    return dict(result,status=STATUS_PASSED,wolfram_report=report,log=log)


def main() -> None:
    parser=argparse.ArgumentParser(description=__doc__)
    selector=parser.add_mutually_exclusive_group(required=True)
    selector.add_argument('--model',choices=MIXED_MODELS)
    selector.add_argument('--all',action='store_true')
    parser.add_argument('--which',choices=['initial','final'],default='final')
    parser.add_argument('--input-dir',type=Path,default=INPUT)
    parser.add_argument('--output-dir',type=Path,default=OUTPUT)
    parser.add_argument('--wolframscript',default='wolframscript')
    parser.add_argument('--prepare-only',action='store_true')
    parser.add_argument('--timeout-sec',type=int,default=1800)
    args=parser.parse_args()
    if args.timeout_sec<=0:parser.error('timeout-sec must be positive')
    chosen=MIXED_MODELS if args.all else (args.model,)
    successes=[]
    for name in chosen:
        r=run_model(name,which=args.which,folder=args.input_dir,
                    outdir=args.output_dir,wolframscript=args.wolframscript,
                    prepare_only=args.prepare_only,timeout_sec=args.timeout_sec)
        masses=r['masses']
        print(f'{name}: heavy eigen-doublet masses {masses[0]:.6f}, {masses[1]:.6f} GeV; '
              + ('Matchete not run' if args.prepare_only else 'Matchete HeavyMassBasis passed')
              + '; physical matching NOT performed',flush=True)
        successes.append(r)
    summary={'format':'t3-mass-eigen-matchete-batch-v1',
             'status':('ParametersPreparedOnly' if args.prepare_only else 'MassEigenHeavyBasisPassed'),
             'models':[r['model'] for r in successes],
             'count':len(successes),
             'physical_matching_performed':False,
             'component_map_certified':False}
    path=args.output_dir/('summary_'+args.which+'.json')
    path.parent.mkdir(parents=True,exist_ok=True)
    path.write_text(json.dumps(summary,indent=2)+'\n',encoding='utf-8')
    print('Summary:',path)


if __name__=='__main__':
    main()
