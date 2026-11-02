"""All-16 audit of current formal fermion-tree release and optional hard inputs.

The default always exports source-backed FORMAL (per unknown signed N_Gamma)
three-heavy-flavour LLSS matrices. Producing numerical conditional hard
inputs requires an independently supplied all-channel scalar NPZ/JSON bundle;
without it no model-specific hard matching is asserted.

Outputs ONLY under output/general_observable/lo_certification_completeness/.
"""
from __future__ import annotations
import argparse
import json
from pathlib import Path
import sys
import numpy as np

ROOT=Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path: sys.path.insert(0,str(ROOT))

from Numerical.joint_scalar_matching.fermion_boundary_bridge import load_release_manifest, load_source
from Numerical.lo_certification.completeness.source_bridge import (
    digest,read_json,split_formal_by_flavour,prepare_conditional_input,SourceBridgeError,
)
from Numerical.lo_certification.completeness.inventory import ALL_MODELS

DEFAULT_RELEASE=Path('output/general_observable/lo_fermion_boundary/releases/source_021dcc1e8fe8ba0f')
DEFAULT_RUNNING=Path('output/general_observable/all16_v1')
OWNED=Path('output/general_observable/lo_certification_completeness')


def slug(model):return model.replace('-','m').replace('+','p')


def main(argv=None):
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--release-dir',type=Path,default=DEFAULT_RELEASE)
    p.add_argument('--running-dir',type=Path,default=DEFAULT_RUNNING)
    p.add_argument('--scalar-input-dir',type=Path,default=None,
                   help='Optional explicitly normalized and source-basis aligned scalar JSON+NPZ bundles')
    p.add_argument('--output-dir',type=Path,required=True)
    p.add_argument('--check-only',action='store_true')
    args=p.parse_args(argv)
    dest=(ROOT/args.output_dir).resolve()
    own=(ROOT/OWNED).resolve()
    if dest==own or own not in dest.parents:
        p.error(f'--output-dir must be new descendant of {own}')
    if dest.exists() and not args.check_only:
        raise FileExistsError('Use a NEW output directory: '+str(dest))
    release=load_release_manifest(ROOT/args.release_dir)
    results=[];objects=[]
    for model in ALL_MODELS:
        source=load_source(release[model],expected_model=model)
        running_path=ROOT/args.running_dir/slug(model)/'running.json'
        running=read_json(running_path)
        if running.get('model')!=model:
            raise SourceBridgeError(f'{model}: wrong saved UV trajectory model')
        formal=split_formal_by_flavour(source,running,digest(running_path))
        status='FORMAL_PER_K_READY_UNSIGNED'
        cinput=None
        if args.scalar_input_dir is not None:
            base=ROOT/args.scalar_input_dir
            metadata_path=base/(slug(model)+'.json')
            raw=read_json(metadata_path)
            npz_name=raw.get('tensor_npz')
            if not isinstance(npz_name,str) or Path(npz_name).name!=npz_name:
                raise SourceBridgeError(f'{model}: scalar NPZ must be adjacent to metadata')
            data_path=metadata_path.parent/npz_name
            with np.load(data_path,allow_pickle=False) as z:
                arrays={key:np.array(z[key]) for key in z.files}
            cinput=prepare_conditional_input(formal,raw,arrays,digest(data_path))
            status='CONDITIONAL_PER_UNKNOWN_SIGNED_NORMALIZATION_NUMERICAL_READY'
        results.append({'model':model,'status':status,'source_sha256':source.hash,
                        'running_sha256':formal.running_hash,
                        'n_source_PL_operators':len(source.source_indices),
                        'n_fermion_flavours':len(formal.fermion_masses_gev),
                        'n_real_heavy_scalars':len(formal.labels),
                        'max_abs_fermion_sum_residual':float(np.max(abs(formal.c_by_fermion.sum(axis=0)-formal.c_total))),
                        'gamma_cc_to_weyl_sign_and_factor':'UNRESOLVED',
                        'complete_HHSS_and_mass_input':cinput is not None,
                        'physical_C5_certified':False})
        objects.append((formal,cinput))
    summary={'format':'t3-lo-hard-source-input-builder-v1','models':results,
             'model_count':len(results),'conditional_numerical_input_count':sum(ci is not None for _,ci in objects),
             'all_fermions_matched_at_one_muF':True,'both_scalar_multiplets_kept_for_joint_muS':True,
             'full_complex_flavour_general_case_supported_by_this_builder':False,
             'source_signed_normalization_certified':False,'physical_full_one_loop_certified':False,
             'missing_for_physical_matching':['signed GammaCC-to-Weyl normalization',
                'complete certified HHSS scalar tensor and mass eigenbasis at common muS',
                'general non-diagonal fermion Takagi vertex decomposition',
                'complete hard counterterm and dim-5 operator closure']}
    if args.check_only:
        print(json.dumps({'model_count':len(results),'statuses':{x['model']:x['status'] for x in results}},indent=2))
        return 0
    dest.mkdir(parents=True,exist_ok=False)
    (dest/'formal_fermion_by_flavour').mkdir()
    if args.scalar_input_dir is not None: (dest/'conditional_inputs').mkdir()
    for formal,cinput in objects:
        name=slug(formal.model)
        fdir=dest/'formal_fermion_by_flavour'
        np.savez_compressed(fdir/(name+'.npz'),
            C_LLSS_by_fermion_formal_per_NGamma=formal.c_by_fermion,
            C_LLSS_total_formal_per_NGamma=formal.c_total,
            fermion_masses_gev=formal.fermion_masses_gev)
        (fdir/(name+'.json')).write_text(json.dumps({
            'model':formal.model,'status':'FORMAL_DIAGONAL_BENCHMARK_ONLY',
            'mu_f_gev':formal.mu_f_gev,'source_sha256':formal.source_hash,
            'running_sha256':formal.running_hash,'source_indices':list(formal.source_indices),
            'scalar_basis':list(formal.labels),'spinor_normalization_certified':False,
            'do_not_treat_as_physical_wilson':True,'tensor_npz':name+'.npz',
            'independent_heavy_fermion_thresholds_used':False},indent=2)+'\n',encoding='utf-8')
        if cinput is not None:
            doc,arrays=cinput
            target=dest/'conditional_inputs'
            np.savez_compressed(target/(name+'.npz'),**arrays)
            (target/(name+'.json')).write_text(json.dumps(doc,indent=2)+'\n',encoding='utf-8')
    (dest/'summary.json').write_text(json.dumps(summary,indent=2)+'\n',encoding='utf-8')
    print(json.dumps({'output':str(dest),'model_count':len(results),
         'conditional_numerical_input_count':summary['conditional_numerical_input_count'],
         'physical_full_one_loop_certified':False},indent=2))
    return 0

if __name__=='__main__': raise SystemExit(main())
