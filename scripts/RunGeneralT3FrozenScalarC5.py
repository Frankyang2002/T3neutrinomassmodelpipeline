"""Direct one-scale C5 with frozen scalar diagonals and an explicit quadratic-mixing scenario.

By default, mu2=0 is an ASSUMPTION, not a value from fitted_reference.
Reuses the already extracted Matchete symbolic coefficient; does NOT run new
Matchete matching, EFT1 running, or reproduce the fitted oscillation results.

    python -m scripts.RunGeneralT3FrozenScalarC5 --all
    python -m scripts.RunGeneralT3FrozenScalarC5 --model B_-1 --mu2-real-gev2 10000
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

from Numerical.experimental.GeneralT3MixedScalarMassBasis import MixedDoubletMassBasis
from scripts.RunGeneralT3NumericalC5 import (
    ROOT, DEFAULT_SOURCE, DEFAULT_PARAMETERS, MODELS, ALIGNED_FORMAT,
    ALIGNED_STATUS, fail, read_json, numerical_model,
)
from scripts.RunGeneralT3ReferenceC5 import (
    MANIFEST, REFERENCE_ID, make_transfer,
)

OUTPUT = ROOT / 'output/general_t3_numerical_c5_frozen_scalar_scenario'


def build_scenario(model: str, manifest: dict, frozen_path: Path, params_path: Path,
                   mu2_gev2: complex, *, mixing_explicit_on_cli: bool = False):
    """Return immutable-reference Yukawas, actual frozen diagonals and a verified U.

    Frozen inputs do NOT contain a general quadratic mixing value. Consequently
    the zero-mixing default is a scenario, not a physical inference about mu2.
    """
    if model not in REFERENCE_ID:
        fail('unsupported frozen reference model')
    old=make_transfer(model,manifest,frozen_path,params_path)
    ordinary=read_json(frozen_path).get('base_state',{}).get('ordinary',{})
    # Refuse to override any quadratic mixing input that appears in future files.
    if any(k in ordinary for k in ('mu2', 'muSq', 'mS12Sq', 'muS12Sq', 'scalarMixingSq')):
        fail('frozen reference now supplies a quadratic mixing; do not overwrite it with a scenario')
    original=old['reference_provenance']
    diagonals=[float(original['original_reference_mS1Sq_gev2']),
               float(original['original_reference_mS2Sq_gev2'])]
    try:mu2=complex(mu2_gev2)
    except (TypeError,ValueError) as exc:raise ValueError('invalid mu2') from exc
    if not np.isfinite([mu2.real,mu2.imag]).all():fail('nonfinite quadratic mixing')
    mass_basis=MixedDoubletMassBasis.build(model,*diagonals,mu2)
    if mass_basis.tachyonic or mass_basis.degenerate or float(np.min(mass_basis.eigenvalues_sq))<=0:
        fail('non-positive or degenerate physical eigenstate mass scenario')
    if (not np.allclose(np.diag(mass_basis.A).real,diagonals,rtol=0,atol=1e-8)
        or not np.allclose(mass_basis.U.conj().T@mass_basis.A@mass_basis.U,
                           np.diag(mass_basis.eigenvalues_sq),rtol=1e-11,atol=1e-8)):
        fail('frozen scalar-mass matrix consistency test failed')
    old.update({
        'format':ALIGNED_FORMAT,
        'benchmark_label':'frozen_reference_yukawa_and_scalar_diagonals_explicit_quadratic_mixing',
        'reference_transfer_only':False,
        'frozen_reference_scalar_diagonals':True,
        'quadratic_mixing_is_explicit_scenario':True,
        'mass_basis_compatible_with_frozen_reference':True,
        'physical_input_certified':False,
        'sequential_EFT1':False,
        'reference_provenance':{
            'frozen_config_path':original['frozen_config_path'],
            'frozen_sha256':original['frozen_sha256'],
            'frozen_scalar_diagonal_masses_squared_gev2':diagonals,
            'quadratic_mixing_mu2_gev2':{'real':float(mu2.real),'imag':float(mu2.imag)},
            'quadratic_mixing_source':('user_supplied_command_line_not_a_fitted_parameter'
                if mixing_explicit_on_cli else 'zero_by_explicit_default_scenario_not_in_fitted_reference'),
            'mass_basis_generated_by':'MixedDoubletMassBasis.build',
            'frozen_yukawas_are_saved_base_state_not_threshold_evolved':True,
            'reproduces_original_fit':False,
            'physical_EFT1_matching_performed':False,
            'nonzero_matchete_rgbeta_quartic_map_certified':False,
        },
    })
    return old,mass_basis


def run_model(model: str, manifest_path: Path, benchmark_dir: Path,
              source_dir: Path, parameters_dir: Path, output_dir: Path,
              mu2_gev2: complex, *, mixing_explicit_on_cli: bool=False):
    manifest=read_json(manifest_path)
    frozen_id=REFERENCE_ID[model][0]
    frozen_path=benchmark_dir/(frozen_id+'.json')
    params_path=parameters_dir/(MODELS[model]+'_final_parameters.json')
    config,basis=build_scenario(model,manifest,frozen_path,params_path,mu2_gev2,
                                mixing_explicit_on_cli=mixing_explicit_on_cli)
    result,destination=numerical_model(model,config,source_dir,parameters_dir,
                                       output_dir,scalar_basis_override=basis)
    if (result.get('status')!=ALIGNED_STATUS or
        result.get('sequential_EFT1_performed') is not False or
        result.get('physical_c5_certified') is not False or
        result.get('frozen_scalar_aligned_one_scale_scenario') is not True or
        result.get('source_mass_parameters_used_as_numeric_scalar_input') is not False):
        fail('scalar-aligned result incorrectly claims certified matching')
    return result,destination


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    select=parser.add_mutually_exclusive_group(required=True)
    select.add_argument('--all',action='store_true')
    select.add_argument('--model',choices=tuple(REFERENCE_ID))
    parser.add_argument('--mu2-real-gev2',type=float,default=None)
    parser.add_argument('--mu2-imag-gev2',type=float,default=None)
    parser.add_argument('--manifest',type=Path,default=MANIFEST)
    parser.add_argument('--benchmark-dir',type=Path,default=MANIFEST.parent)
    parser.add_argument('--source-dir',type=Path,default=DEFAULT_SOURCE)
    parser.add_argument('--parameters-dir',type=Path,default=DEFAULT_PARAMETERS)
    parser.add_argument('--output-dir',type=Path,default=OUTPUT)
    args=parser.parse_args()
    cli_mix=args.mu2_real_gev2 is not None or args.mu2_imag_gev2 is not None
    mu2=complex(args.mu2_real_gev2 or 0.0,args.mu2_imag_gev2 or 0.0)
    records=[]
    for model in (tuple(REFERENCE_ID) if args.all else (args.model,)):
        result,path=run_model(model,args.manifest,args.benchmark_dir,
            args.source_dir,args.parameters_dir,args.output_dir,mu2,
            mixing_explicit_on_cli=cli_mix)
        records.append({'model':model,'output':str(path),
            'takagi_masses_eV_descending':result['takagi_singular_values_eV_descending'],
            'scalar_masses_gev':result['scalar_masses_gev']})
        print(f'{model}: frozen scalar diagonals, mu2={mu2} GeV^2; '
              f'eigenstate masses {result["scalar_masses_gev"]} GeV; '
              f'largest one-scale Takagi mass {result["largest_takagi_mass_eV"]:.8g} eV; '
              'NOT sequential EFT1 or a reproduced oscillation fit',flush=True)
    args.output_dir.mkdir(parents=True,exist_ok=True)
    summary=args.output_dir/'summary_final.json'
    summary.write_text(json.dumps({
        'status':ALIGNED_STATUS,
        'quadratic_mixing_mu2_gev2':{'real':mu2.real,'imag':mu2.imag},
        'quadratic_mixing_input_kind':('explicit_command_line_scenario' if cli_mix else 'explicit_default_zero_scenario'),
        'frozen_scalar_diagonal_masses_used':True,
        'physical_c5_certified':False,
        'sequential_EFT1_performed':False,
        'reproduces_original_fitted_result':False,
        'models':records,
    },indent=2)+'\n',encoding='utf-8')
    print('Summary:',summary)

if __name__=='__main__':main()
