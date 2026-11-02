"""B_-1/C_-1 exact Matchete-CG numerical F-threshold Wilson boundaries.

Works with synchronized mu_F inputs; obtains missing native CG arrays by
instrumenting a temporary copy of the *current local* Matchete matching .wl.
The original .wl and all existing 14-model results remain unchanged.

Example PowerShell:
 python -X utf8 -m scripts.RunT3GeneralMixedWilsonTensor --all --source-only
 python -X utf8 -m scripts.RunT3GeneralMixedWilsonTensor --all
 python -X utf8 -m scripts.RunT3GeneralMixedWilsonTensor --all --reuse-cg
"""
from __future__ import annotations

import argparse
from hashlib import sha256
import json
from pathlib import Path
import shutil
import subprocess
import tempfile
import time
from typing import Mapping

import numpy as np

from Numerical.running.GeneralT3MixedWilsonTensor import (
    FORMAT, MIXED, build_mixed_source_terms, project_mixed_terms, read_json,
    save_json, sha,
)
from scripts.RunT3GeneralWilsonBoundary import _running_basis

ROOT = Path(__file__).resolve().parents[1]
SYNC = Path('output/general_observable/mixed_f_synchronised_v1')
PRIOR = Path('output/general_observable/ordered_wilson_tensor_v2')
OUTPUT = Path('output/general_observable/ordered_wilson_tensor_v3')
WOLFRAM_SOURCE = Path('Lagrangian/RunGeneralT3FermionFirstEFT1Tree.wl')
EXPECTED = {"B_-1":20, "C_-1":24}
CG_FORMAT = 't3-mixed-matchete-native-cg-registry-v1'
CG_STAGE_ANCHOR = ' T3MECheckpoint["fermion-first-LLSS-sector-extracted"];'

# Inserted into a fresh temporary copy of the user's local Wolfram program
# immediately after it has matched F and selected the LLSS terms, while
# Matchete still retains the exact CG definitions (including eps[SU2L]).
CG_INSERTION = r'''
 (* Added by RunT3GeneralMixedWilsonTensor to an EPHEMERAL script only.
    Same extraction algorithm as the project's RunMatching.wl. *)
 $T3MEWilsonCGNames = DeleteDuplicates @ Cases[
   HoldComplete @@ llssTerms,
   CG[name_, ___] :> (name /. Bar[z_] :> z), Infinity];
 $T3MEWilsonCGData = DeleteCases[Map[
   Function[name, Module[{reps,tensor},
     reps = Quiet@Check[Matchete`PackageScope`$CGproperties[name, Indices], $Failed];
     tensor = Quiet@Check[Matchete`CGManipulations`PackagePrivate`$CGtensors[name], $Failed];
     If[reps === $Failed || tensor === $Failed || MissingQ[reps] || MissingQ[tensor],
       Nothing,
       <|"Name" -> ToString[Unevaluated[name], InputForm],
         "RepsInputForm" -> ToString[InputForm[reps]],
         "TensorInputForm" -> ToString[InputForm[tensor]]|>]
   ]], $T3MEWilsonCGNames], Nothing];
 If[Length[$T3MEWilsonCGNames] < 1 ||
    Length[$T3MEWilsonCGNames] =!= Length[$T3MEWilsonCGData],
    T3MEFail["MixedWilsonCGRegistryIncomplete",
       ToString[$T3MEWilsonCGNames, InputForm]]];
 $T3MEWilsonCGPath = FileNameJoin[{DirectoryName[outPath],
      FileBaseName[outPath] <> "_llss_cg_registry.json"}];
 If[Quiet@Check[Export[$T3MEWilsonCGPath,
      <|"format" -> "t3-mixed-matchete-native-cg-registry-v1",
        "model" -> modelName,
        "source" -> "Native CG tensors in same F-only Matchete matching run",
        "CGRegistry" -> $T3MEWilsonCGData|>, "RawJSON"], $Failed] === $Failed,
    T3MEFail["MixedWilsonCGRegistryExportFailed"]];
 T3MECheckpoint["fermion-first-LLSS-native-CG-registry-exported"];
'''


def tag(model: str) -> str:
    return model.replace('-', 'm').replace('+', 'p')


def absolute(p: Path | str) -> Path:
    q = Path(p)
    return q.resolve() if q.is_absolute() else (ROOT / q).resolve()


def validate_sync(model: str, sync: Path) -> tuple[dict, Path, Path, Path, dict, dict]:
    s_path = sync / tag(model) / 'status.json'
    status = read_json(s_path)
    if (status.get('format') != 't3-mixed-F-boundary-synchronised-to-running-v1' or
        status.get('status') != 'CurrentRunningMassBasisFThresholdTreeLLSSExported' or
        status.get('model') != model or
        status.get('source_operator_matching_at_running_mu_F') is not True or
        status.get('historical_mu_S_basis_reused') is not False or
        status.get('llss_term_count') != EXPECTED[model]):
        raise ValueError(f'{model}: source is not the verified synchronized F-threshold export')
    llss = absolute(status['matched_llss_file'])
    mass = absolute(status['mass_basis_file'])
    running = absolute(status['running_file'])
    for p in (llss, mass, running):
        if not p.is_file():
            raise FileNotFoundError(f'Missing local synchronized source: {p}')
    basis = read_json(mass)
    provenance = basis.get('running_boundary_provenance', {})
    if (basis.get('model') != model or
        provenance.get('running_sha256') != sha(running) or
        float(provenance.get('boundary_scale_gev', -1)) != float(status['mu_F_gev'])):
        raise ValueError(f'{model}: rotated mass basis is stale relative to numerical mu_F state')
    running_data = read_json(running)
    return status, llss, mass, running, basis, running_data


def instrument_source(wl_source: str) -> str:
    if wl_source.count(CG_STAGE_ANCHOR) != 1:
        raise ValueError('Local Wolfram exporter changed; missing/duplicate injection anchor')
    if 'T3MEPerform[]' not in wl_source or 'llssTerms=Select' not in wl_source:
        raise ValueError('Unexpected Matchete source structure; refusing instrumentation')
    return wl_source.replace(CG_STAGE_ANCHOR, CG_INSERTION + CG_STAGE_ANCHOR, 1)


def new_wolfram_output(workdir: Path) -> tuple[Path, Path, Path]:
    target = workdir / 'mixed_native_tree.json'
    return (target, workdir / 'mixed_native_tree_llss_tree.inputform.txt',
            workdir / 'mixed_native_tree_llss_cg_registry.json')


def verify_native_parameters(model: str, status: Mapping) -> Path:
    """Fail if the Matchete input U differs from the synchronized mu_F U."""
    param = absolute(status['native_parameters_file'])
    mass_path = absolute(status['mass_basis_file'])
    if not param.is_file() or not mass_path.is_file():
        raise FileNotFoundError(f'{model}: missing synchronized native Matchete parameters or basis')
    obj = read_json(param)
    mass = read_json(mass_path)
    if obj.get('model') != model or obj.get('endpoint') != 'final':
        raise ValueError(f'{model}: native parameter model/endpoint mismatch')
    a = obj.get('unitary_U')
    b = mass['basis']['unitary_flavour_rotation_U']
    if not isinstance(a,dict) or not isinstance(b,dict):
        raise ValueError('Uncertified native mass-basis U')
    au = np.asarray(a['real'],float)+1j*np.asarray(a['imag'],float)
    bu = np.asarray(b['real'],float)+1j*np.asarray(b['imag'],float)
    if au.shape != (2,2) or bu.shape != (2,2) or not np.allclose(au,bu,atol=1e-12,rtol=1e-12):
        raise ValueError(f'{model}: native Matchete U is not current running-mass U')
    return param


def run_native_registry(model: str, status: Mapping, workdir: Path, *,
                        source: Path, executable: str, timeout: int) -> tuple[Path,Path]:
    param = verify_native_parameters(model, status)
    if not source.is_file():
        raise FileNotFoundError(f'Current local Wolfram source absent: {source}')
    content = instrument_source(source.read_text(encoding='utf-8-sig'))
    target, llss, registry = new_wolfram_output(workdir)
    workdir.mkdir(parents=True, exist_ok=True)
    if not shutil.which(executable) and not Path(executable).is_file():
        raise FileNotFoundError(f'Cannot execute {executable}; this export requires local Wolfram/Matchete')
    # The temporary file must live BESIDE the original WL program because
    # that program resolves its model/RunMatching.wl dependencies from $InputFileName.
    fd, temp_name = tempfile.mkstemp(prefix='__t3_mixed_wilson_cg_', suffix='.wl', dir=str(source.parent))
    try:
        import os
        with os.fdopen(fd,'w',encoding='utf-8') as f:
            f.write(content)
        cmd = [executable, '-file', temp_name, str(param), str(target)]
        proc = subprocess.run(cmd, cwd=str(ROOT), capture_output=True, text=True,
                              encoding='utf-8', errors='replace', timeout=timeout, check=False)
        (workdir/'mixed_matchete_cg_export.log').write_text(
            'returncode: '+str(proc.returncode)+'\nstdout:\n'+proc.stdout+'\nstderr:\n'+proc.stderr,
            encoding='utf-8')
    finally:
        Path(temp_name).unlink(missing_ok=True)
    if proc.returncode != 0 or not target.is_file():
        raise RuntimeError(f'{model}: instrumented Matchete matching failed, exit={proc.returncode}; '
                           f'see {workdir / "mixed_matchete_cg_export.log"}')
    result = read_json(target)
    if (result.get('status') != 'FermionFirstTreeLLSSExtractedNotSequentialEFT1' or
        result.get('model') != model or
        result.get('tree_matching_completed') is not True or
        result.get('llss_term_count') != EXPECTED[model] or
        not llss.is_file() or not registry.is_file()):
        raise ValueError(f'{model}: native tree/CG registry output incomplete')
    registry_data = read_json(registry)
    if registry_data.get('format') != CG_FORMAT or registry_data.get('model') != model:
        raise ValueError(f'{model}: bad CG registry schema')
    save_json(workdir / 'mixed_source_cg_certificate.json', {
        'model': model,
        'native_parameters_sha256': sha(param),
        'native_wolfram_source_sha256': sha(source),
        'new_llss_sha256': sha(llss),
        'new_native_cg_sha256': sha(registry),
        'matching_report_sha256': sha(target),
    })
    return llss, registry


def reuse_native_registry(model: str, status: Mapping, workdir: Path, source: Path) -> tuple[Path,Path]:
    target, llss, registry = new_wolfram_output(workdir)
    cert = read_json(workdir/'mixed_source_cg_certificate.json')
    params = verify_native_parameters(model,status)
    expected = {'model':model, 'native_parameters_sha256':sha(params),
                'native_wolfram_source_sha256':sha(source),
                'new_llss_sha256':sha(llss), 'new_native_cg_sha256':sha(registry),
                'matching_report_sha256':sha(target)}
    if cert != expected:
        raise ValueError(f'{model}: cached Matchete CG registry provenance does not match local sources')
    data = read_json(target)
    if data.get('status') != 'FermionFirstTreeLLSSExtractedNotSequentialEFT1':
        raise ValueError('Cached Matchete matching not successful')
    return llss, registry


def _coefficient_sources(model: str, status: dict, basis: dict, running: Path,
                         llss: Path) -> tuple[list[dict],dict]:
    run, masses, yukawas, diagnostics = _running_basis(model, running, {'y1','y2'})
    if float(run['scales_gev']['F']) != float(status['mu_F_gev']):
        raise ValueError(f'{model}: running mu_F disagrees with matched basis')
    operators = build_mixed_source_terms(
        model, llss.read_text(encoding='utf-8-sig'), basis, masses,
        yukawas, expected_count=EXPECTED[model])
    info = {'mu_F_gev': float(status['mu_F_gev']),
            'heavy_fermion_masses_gev': np.asarray(masses,float).tolist(),
            'heavy_flavour_basis_diagnostics': diagnostics,
            'llss_sha256': sha(llss), 'running_sha256':sha(running)}
    return operators, info


def carry_previous_14(prior: Path, output: Path, selected: set[str]) -> list[dict]:
    if not prior.is_file():
        raise FileNotFoundError(f'Previous 14-model summary missing: {prior}')
    report = read_json(prior)
    rows = report.get('results')
    if (report.get('ordered_real_tensor_count') != 14 or
        report.get('blocked_count') != 0 or not isinstance(rows,list)):
        raise ValueError('Prior output is not the verified 14/16 ordered tensor batch')
    copied=[]
    for entry in rows:
        model = entry['model']
        if model in MIXED or model in selected:
            continue
        if entry.get('status') != 'OrderedNumericalRealComponentTensorFromMatchete':
            raise ValueError('Prior ordinary-model tensor record not successful')
        source=prior.parent / tag(model) / 'ordered_tensor.json'
        payload=read_json(source)
        if payload.get('status') != entry['status'] or payload.get('model') != model:
            raise ValueError('Previous 14-model tensor/source mismatch')
        target=output/tag(model)/'ordered_tensor.json'
        target.parent.mkdir(parents=True,exist_ok=True)
        shutil.copyfile(source,target)
        copied.append({**entry, 'tensor_file':str(target), 'carried_previous_batch':True,
                       'tensor_sha256':sha(target)})
    if len(copied)!=14:
        raise ValueError(f'Expected 14 previous tensors, found {len(copied)}')
    return copied


def main(argv: list[str] | None = None) -> int:
    ap=argparse.ArgumentParser(description=__doc__)
    selection=ap.add_mutually_exclusive_group(required=True)
    selection.add_argument('--all',action='store_true')
    selection.add_argument('--model',choices=MIXED)
    ap.add_argument('--source-only',action='store_true',help='Numerical Yukawa factorization only; no CG certification')
    ap.add_argument('--reuse-cg',action='store_true',help='Reuse registry from this runner only if all SHA checks pass')
    ap.add_argument('--synchronised-root',type=Path,default=SYNC)
    ap.add_argument('--prior-ordered-root',type=Path,default=PRIOR)
    ap.add_argument('--output-root',type=Path,default=OUTPUT)
    ap.add_argument('--wolfram-source',type=Path,default=WOLFRAM_SOURCE)
    ap.add_argument('--wolframscript',default='wolframscript')
    ap.add_argument('--timeout-sec',type=int,default=3600)
    args=ap.parse_args(argv)
    if args.source_only and args.reuse_cg:
        ap.error('--source-only and --reuse-cg are mutually exclusive')
    if args.timeout_sec < 1:
        ap.error('--timeout-sec must be positive')
    selected=MIXED if args.all else (args.model,)
    out=absolute(args.output_root)
    sync=absolute(args.synchronised_root)
    wolfram_source=absolute(args.wolfram_source)
    rows=[]
    if args.all and not args.source_only:
        try:
            rows=carry_previous_14(absolute(args.prior_ordered_root)/'summary.json',out,set(selected))
        except Exception as err:
            ap.error(f'Cannot safely carry existing 14 tensors: {err}')
    for model in selected:
        dest=out/tag(model)
        try:
            status, old_llss, mass, running, basis, _ = validate_sync(model,sync)
            if args.source_only:
                llss, registry = old_llss, None
            elif args.reuse_cg:
                llss, registry=reuse_native_registry(model,status,dest,wolfram_source)
            else:
                llss, registry=run_native_registry(model,status,dest,source=wolfram_source,
                                                 executable=args.wolframscript,timeout=args.timeout_sec)
            terms,info=_coefficient_sources(model,status,basis,running,llss)
            source_payload={'format':'t3-mixed-factorised-LLSS-source-v1','model':model,
                 'status':'NumericalOrderedSourceFlavourNoCG' if registry is None else 'NumericalOrderedSourceFlavourCGExported',
                 'source_status_sha256':sha(sync/tag(model)/'status.json'),
                 'source_mass_basis_sha256':sha(mass), **info,
                 'operator_terms':terms,
                 'canonical_C_ijab_normalization_certified':False,
                 'eft1_wilson_rge_computed':False}
            save_json(dest/'source_flavour.json',source_payload)
            if registry is None:
                row={'model':model,'status':'MixedOrderedFlavourSourcesOnly_CGNotExported',
                     'tree_operator_count_PL':sum(t['chirality']=='PL' for t in terms),
                     'raw_llss_term_count':len(terms), 'real_component_count':None,
                     'flavour_file':str(dest/'source_flavour.json')}
            else:
                cg_data=read_json(registry)
                operators=project_mixed_terms(model,terms,cg_data)
                payload={'format':FORMAT,'model':model,
                         'status':'OrderedNumericalRealComponentTensorFromMatchete',
                         'source_synced_matchete_llss_sha256':sha(old_llss),
                         'source_native_fresh_llss_sha256':sha(llss),
                         'source_native_cg_sha256':sha(registry),
                         'source_status_sha256':sha(sync/tag(model)/'status.json'),
                         'source_mass_basis_sha256':sha(mass),
                         'mu_F_gev':info['mu_F_gev'],
                         'scalar_basis_explicit':True,
                         'matchete_cg_contractions_retained':True,
                         'external_three_family_flavour_order_not_symmetrized':True,
                         'tree_operator_count_PL':len(operators),
                         'total_ordered_real_component_count':sum(x['component_count'] for x in operators),
                         'operators':operators,
                         'canonical_C_ijab_normalization_certified':False,
                         'physical_matchete_rgbeta_coupling_dictionary_certified':False,
                         'full_eft1_wilson_anomalous_dimension_implemented':False,
                         'eft1_wilson_rge_computed':False,
                         'scalar_threshold_c5_computed':False}
                save_json(dest/'ordered_tensor.json',payload)
                row={'model':model, 'status':payload['status'],
                     'tree_operator_count_PL':len(operators),
                     'real_component_count':payload['total_ordered_real_component_count'],
                     'tensor_file':str(dest/'ordered_tensor.json'),
                     'source':str(dest/'source_flavour.json')}
        except Exception as exc:
            row={'model':model,'status':'Blocked','error':f'{type(exc).__name__}: {exc}'}
            save_json(dest/'error.json',row)
        rows.append(row)
        print(f"{model}: {row['status']}" + ('; '+row['error'] if 'error' in row else ''),flush=True)
        summary={'format':'t3-ordered-source-flavour-and-real-scalar-tensor-v1','results':rows,
                 'ordered_real_tensor_count':sum(x.get('status')=='OrderedNumericalRealComponentTensorFromMatchete' for x in rows),
                 'synchronized_mixed_source_pending_count':sum(x.get('status')=='MixedOrderedFlavourSourcesOnly_CGNotExported' for x in rows),
                 'blocked_count':sum(x.get('status')=='Blocked' for x in rows),
                 'physical_sequential_matching_certified':False,
                 'full_eft1_wilson_rge_implemented':False}
        save_json(out/'summary.json',summary)
    print(f"Ordered real tensors: {summary['ordered_real_tensor_count']}; "
          f"pending mixed CG: {summary['synchronized_mixed_source_pending_count']}; "
          f"blocked: {summary['blocked_count']}; EFT1 running NOT CERTIFIED")
    return 1 if summary['blocked_count'] else 0


if __name__=='__main__':
    raise SystemExit(main())
