"""Numerically evaluate the EXACT uploaded Matchete ordered Weinberg expressions.

Experimental simultaneous all-heavy one-scale result only.  No sequential
EFT1 matching, parameter fit, or Matchete/RGBeta quartic map certification.

From the repository root:
  python -m scripts.RunGeneralT3NumericalC5 --all
  python -m scripts.RunGeneralT3NumericalC5 --model B_-1 --config configs/general_t3/c5_illustrative.json
"""
from __future__ import annotations

import argparse
import ast
import csv
from hashlib import sha256
import json
from pathlib import Path
from typing import Any

import mpmath as mp
import numpy as np

from Numerical.experimental.GeneralT3MixedScalarMassBasis import (
    MixedDoubletMassBasis, matrix_json,
)

ROOT = Path(__file__).resolve().parents[1]
DEFAULT_SOURCE = ROOT / 'output/general_t3_mass_eigen_matchete/one_loop_probe'
DEFAULT_PARAMETERS = ROOT / 'output/general_t3_mass_eigen_matchete'
DEFAULT_OUTPUT = ROOT / 'output/general_t3_numerical_c5_exploratory'
DEFAULT_CONFIG = ROOT / 'configs/general_t3/c5_illustrative.json'
MODELS = {'B_-1': 'B_m1', 'C_-1': 'C_m1'}
STATUS = 'IllustrativeDirectOneScaleC5NumericallyEvaluatedNotPhysicalEFT1'
REFERENCE_FORMAT = 't3-frozen-reference-yukawa-transfer-one-scale-v1'
REFERENCE_STATUS = 'FrozenReferenceYukawaTransferToExperimentalScalarBasisNotPhysicalPrediction'
ALIGNED_FORMAT = 't3-frozen-scalar-diagonal-direct-C5-input-v1'
ALIGNED_STATUS = 'FrozenScalarDiagonalDirectOneScaleC5NotPhysicalEFT1'
ALLOWED_COUPLINGS = {'MF','MSX','MSY','T3MEU11','T3MEU12','T3MEU21','T3MEU22',
                     'lambdaT3','lambdaT3ExtraQ020020I1','lambdaT3ExtraQ020200I1','y1','y2'}


def fail(msg: str) -> None:
    raise ValueError('T3 numerical Matchete C5 fail-closed: '+msg)


def read_json(path: Path) -> dict:
    if not path.is_file():
        fail(f'missing input file {path}')
    data=json.loads(path.read_text(encoding='utf-8'))
    if not isinstance(data,dict):fail(f'not a JSON object: {path}')
    return data


def complex_value(obj: Any, name: str) -> mp.mpc:
    if not isinstance(obj,dict) or set(obj) != {'real','imag'}:fail(f'{name}: complex value requires real and imag')
    try:
        z=mp.mpc(str(obj['real']),str(obj['imag']))
    except (TypeError,ValueError) as exc:raise ValueError(f'{name}: invalid complex number') from exc
    if not(mp.isfinite(z.real) and mp.isfinite(z.imag)):fail(f'{name}: nonfinite')
    return z


def matrix_values(obj: Any, name: str, *, rows: int | None = None):
    if not isinstance(obj,list) or not obj or not all(isinstance(row,list) and row for row in obj):fail(f'{name}: invalid matrix')
    width=len(obj[0]);
    if any(len(row)!=width for row in obj) or (rows is not None and len(obj)!=rows):fail(f'{name}: invalid matrix shape')
    return [[complex_value(x,f'{name}[{i},{j}]') for j,x in enumerate(row)] for i,row in enumerate(obj)]


def _interpolate_couplings(source: str) -> str:
    """Translate a tiny, checked subset of Mathematica InputForm to Python.

    No eval of the source string or user code. Unknown heads, indices or
    couplings fail closed.  d$$1/d$$2 are external lepton labels; the one
    internal NFlavor index is summed explicitly in evaluate_ordered().
    """
    s=source.strip()
    if not s or len(s)>2_000_000:fail('empty or oversize ordered coefficient')
    for field in ('y1','y2'):
        for external,tag in [('d$$1','a'),('d$$2','b')]:
            key=f'Coupling[{field}, {{Index[{external}, Flavor], Index[d$$1, NFlavor]}}, 0]'
            s=s.replace(f'Bar[{key}]',f'conj({field}{tag})').replace(key,f'{field}{tag}')
    for key in sorted(ALLOWED_COUPLINGS-{'y1','y2'},key=len,reverse=True):
        literal=f'Coupling[{key}, {{}}, 0]'
        s=s.replace('Bar['+literal+']','conj('+key+')').replace(literal,key)
    for head in ('Sqrt','Log'):
        s=s.replace(head+'[',head.lower()+'(')
    s=s.replace(']',')')
    if '[' in s or ']' in s or 'Coupling' in s or 'Bar' in s or 'Index' in s or '$' in s:
        fail('unrecognized Matchete expression/index/coupling; refusing evaluation')
    if any(ch not in '0123456789.+-*/()^ _abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ\n\r\t' for ch in s):
        fail('unsupported characters in Matchete expression')
    s=s.replace('^','**')
    try:
        tree=ast.parse(s,mode='eval')
    except (SyntaxError,MemoryError,RecursionError) as exc:
        raise ValueError('T3 numerical Matchete C5 fail-closed: expression parse failed') from exc
    allowed_vars=ALLOWED_COUPLINGS|{'y1a','y1b','y2a','y2b','hbar'}
    # Explicit safety gate: prohibit source from executing arbitrary Python.
    for node in ast.walk(tree):
        if isinstance(node,(ast.Expression,ast.Load,ast.Add,ast.Sub,ast.Mult,
                            ast.Div,ast.Pow,ast.UAdd,ast.USub,ast.BinOp,
                            ast.UnaryOp,ast.Constant,ast.Call,ast.Name)):
            pass
        else:fail(f'unsupported AST node {type(node).__name__}')
        if isinstance(node,ast.Constant) and (not isinstance(node.value,(int,float)) or isinstance(node.value,bool)):
            fail('unsupported literal')
        if isinstance(node,ast.Name) and node.id not in allowed_vars|{'conj','sqrt','log'}:
            fail('unknown Matchete name: '+node.id)
        if isinstance(node,ast.Call) and (not isinstance(node.func,ast.Name) or
                                          node.func.id not in {'conj','sqrt','log'} or
                                          len(node.args)!=1 or node.keywords):
            fail('unrecognized Mathematica function')
    return s


def compile_coefficient(expr: str):
    checked=_interpolate_couplings(expr)
    return compile(checked,'<checked-Matchete-C5>', 'eval')


def _as_array(values):
    return np.array([[complex(z) for z in row] for row in values],dtype=np.complex128)


def _spec_json(z: complex):
    return {'real':float(z.real),'imag':float(z.imag)}


def numerical_model(model: str,config: dict,source_dir: Path,parameters_dir: Path,output_dir: Path,
                    *, scalar_basis_override: MixedDoubletMassBasis | None = None):
    if model not in MODELS:fail('unsupported model')
    stem=MODELS[model]+'_final'
    probe_path=source_dir/(stem+'_one_loop_probe.json')
    coef_path=source_dir/(stem+'_one_loop_probe_weinberg_ordered_coefficient.inputform.txt')
    params_path=parameters_dir/(stem+'_parameters.json')
    report=read_json(probe_path);params=read_json(params_path)
    if (report.get('status')!='OneLoopDirectSymbolicWeinbergProbePassedNumericalC5Pending' or
        report.get('model')!=model or report.get('one_loop_match_stage_completed') is not True or
        report.get('weinberg_extraction_status')!='Success' or
        report.get('weinberg_holomorphic_term_count')!=118 or
        report.get('ordered_coefficient_is_canonical_flavour_symmetric_c5') is not False or
        report.get('physical_matching_performed') is not False or
        report.get('production_sequential_matching_performed') is not False):
        fail(model+': missing/unexpected one-loop source certification')
    if (params.get('model')!=model or params.get('verified_eight_real_mode_hessian') is not True or
        params.get('verified_17_tensor_rotation') is not True):
        fail(model+': mass input certification missing')
    if not coef_path.is_file():fail('missing ordered coefficient '+str(coef_path))
    text=coef_path.read_text(encoding='utf-8')
    used=set(__import__('re').findall(r'Coupling\[([A-Za-z0-9]+)',text))
    if used != ALLOWED_COUPLINGS:fail(f'{model}: unexpected or incomplete symbol set: {sorted(used)}')
    code=compile_coefficient(text)
    if (config.get('format') not in {'t3-illustrative-direct-C5-input-v1', REFERENCE_FORMAT, ALIGNED_FORMAT} or
        config.get('physical_input_certified') is not False or
        config.get('sequential_EFT1') is not False):
        fail('direct one-scale input cannot claim physical certification')
    if config['format']==REFERENCE_FORMAT and (config.get('reference_transfer_only') is not True or
                                              config.get('mass_basis_compatible_with_frozen_reference') is not False):
        fail('reference transfer must explicitly record mass mismatch')
    aligned=config.get('format')==ALIGNED_FORMAT
    if aligned != (scalar_basis_override is not None):
        fail('aligned direct C5 format requires explicit validated scalar mass-basis override')
    if aligned:
        if (config.get('frozen_reference_scalar_diagonals') is not True or
            config.get('quadratic_mixing_is_explicit_scenario') is not True or
            config.get('mass_basis_compatible_with_frozen_reference') is not True or
            config.get('reference_transfer_only') is not False):
            fail('frozen scalar-aligned benchmark missing explicit scenario guard')
        if not isinstance(scalar_basis_override,MixedDoubletMassBasis) or scalar_basis_override.model != model:
            fail('scalar-basis override is not a validated mass basis for the selected model')
        scalar_basis_override.verify()
        if scalar_basis_override.degenerate or scalar_basis_override.tachyonic or np.min(scalar_basis_override.eigenvalues_sq)<=0:
            fail('singular or tachyonic scalar mass basis in frozen benchmark')
        reference= config.get('reference_provenance',{})
        diagonals=reference.get('frozen_scalar_diagonal_masses_squared_gev2')
        if (not isinstance(diagonals,list) or len(diagonals)!=2 or
            not np.allclose(np.diag(scalar_basis_override.A).real,np.asarray(diagonals,dtype=float),rtol=0,atol=1e-8)):
            fail('scalar basis does not reproduce the frozen reference diagonal squared masses')
        mu=reference.get('quadratic_mixing_mu2_gev2')
        if (not isinstance(mu,dict) or set(mu)!={'real','imag'} or
            abs(complex(float(mu['real']),float(mu['imag']))-scalar_basis_override.mu_sq)>1e-8):
            fail('quadratic mixing not equal to explicitly recorded scenario')
    mp.mp.dps=60
    if aligned:
        masses=[float(np.sqrt(x)) for x in scalar_basis_override.eigenvalues_sq]
        u=matrix_json(scalar_basis_override.U)
    else:
        masses=params.get('heavy_masses_gev',[])
        u=params['unitary_U']
    if not isinstance(masses,list) or len(masses)!=2 or any(float(x)<=0 for x in masses):fail('invalid scalar masses')
    re=u['real'];im=u['imag']
    if len(re)!=2 or any(len(x)!=2 for x in re):fail('invalid mixing U shape')
    u_vars={f'T3MEU{i+1}{j+1}':mp.mpc(str(re[i][j]),str(im[i][j])) for i in range(2) for j in range(2)}
    y1=matrix_values(config.get('y1'), 'y1',rows=3)
    y2=matrix_values(config.get('y2'), 'y2',rows=3)
    if len(y1[0])!=len(y2[0]) or len(y1[0])<1:fail('Yukawa flavour dimensions incompatible')
    if 'fermion_masses_gev' in config:
        raw_masses=config['fermion_masses_gev']
        if not isinstance(raw_masses,list) or len(raw_masses)!=len(y1[0]):
            fail('fermion mass eigenvalue count does not match Yukawa columns')
        try:
            fermion_masses=[mp.mpf(str(v)) for v in raw_masses]
        except (ValueError, TypeError):
            fail('invalid fermion mass eigenvalues')
    else:
        try: common_mass=mp.mpf(str(config['fermion_mass_gev']))
        except (KeyError, ValueError, TypeError): fail('invalid common fermion mass')
        fermion_masses=[common_mass]*len(y1[0])
    if any(not mp.isfinite(m) or m<=0 for m in fermion_masses):
        fail('fermion masses must be finite and strictly positive')
    couplings={x:complex_value(config['quartics'][x],x)
               for x in ('lambdaT3','lambdaT3ExtraQ020020I1','lambdaT3ExtraQ020200I1')}
    env={'MSX':mp.mpf(str(masses[0])),'MSY':mp.mpf(str(masses[1])),
         'hbar':mp.mpf(1)/(16*mp.pi*mp.pi),
         'conj':mp.conj,'sqrt':mp.sqrt,'log':mp.log,**u_vars,**couplings}
    ordered=[[mp.mpc(0) for _ in range(3)] for _ in range(3)]
    for i in range(3):
        for j in range(3):
            for k in range(len(y1[0])):
                local=dict(env,MF=fermion_masses[k],y1a=y1[i][k],y1b=y1[j][k],y2a=y2[i][k],y2b=y2[j][k])
                try:
                    val=eval(code,{'__builtins__':{}},local)
                except (ZeroDivisionError,OverflowError,ValueError) as exc:
                    raise ValueError(model+': numerical coefficient evaluation failed') from exc
                if not (mp.isfinite(mp.re(val)) and mp.isfinite(mp.im(val))):fail('nonfinite coefficient')
                ordered[i][j]+=val
    # Production RunMatching.wl: C5_ij = A5_ij + A5_ji; do not double-count loop factor.
    canonical=[[ordered[i][j]+ordered[j][i] for j in range(3)] for i in range(3)]
    v=mp.mpf(str(config.get('vev_gev',246.0)))
    if v<=0:fail('vev must be positive')
    mnu=[[-v*v*canonical[i][j]/2 for j in range(3)] for i in range(3)]
    mnu_np=_as_array(mnu)
    takagi=np.linalg.svd(mnu_np,compute_uv=False) * 1.e9
    if any(not np.isfinite(x) for x in takagi):fail('nonfinite Takagi masses')
    max_abs=max(abs(complex(z)) for row in mnu for z in row)
    sym_err=float(max(abs(complex(canonical[i][j]-canonical[j][i])) for i in range(3) for j in range(3)))
    payload={
        'format':'t3-direct-one-scale-numerical-c5-v1',
        'status':(ALIGNED_STATUS if aligned else REFERENCE_STATUS if config['format']==REFERENCE_FORMAT else STATUS),
        'model':model,
        'reference_yukawa_transfer_only':config['format']==REFERENCE_FORMAT,
        'frozen_scalar_aligned_one_scale_scenario':aligned,
        'reproduces_frozen_fitted_result':False if aligned or config['format']==REFERENCE_FORMAT else None,
        'mass_basis_compatible_with_frozen_reference':(True if aligned else False if config['format']==REFERENCE_FORMAT else None),
        'input_semantics':('Frozen fitted-reference Yukawas and scalar mass DIAGONALS; additional quadratic mixing is explicit scenario input; not EFT1 evolved or fitted'
                           if aligned else 'Frozen fitted-reference Yukawas transferred to a DIFFERENT experimental scalar mass benchmark; '
                           'not original frozen fit or EFT1 evolved'
                           if config['format']==REFERENCE_FORMAT else
                           'Illustrative couplings supplied by user or labelled default; not EFT1 evolved'),
        'physical_c5_certified':False,'sequential_EFT1_performed':False,
        'matching_approach':'SimultaneousAllHeavyOneScaleExploratory',
        'weinberg_conversion':'C5_ij=A5_ij+A5_ji; mnu=-v^2*C5/2',
        'loop_factor':'hbar=1/(16*pi^2)',
        'source_coefficient_sha256':sha256(coef_path.read_bytes()).hexdigest(),
        'source_probe_report_sha256':sha256(probe_path.read_bytes()).hexdigest(),
        'source_mass_parameters_sha256':sha256(params_path.read_bytes()).hexdigest(),
        'source_mass_parameters_used_as_numeric_scalar_input':not aligned,
        'effective_scalar_basis_A_gev2':matrix_json(scalar_basis_override.A) if aligned else params.get('original_mass_squared_A_gev2'),
        'effective_scalar_basis_U':u,
        'effective_scalar_basis_source':('Frozen reference scalar diagonals with explicit quadratic mixing scenario'
                                         if aligned else 'Original experimentally prepared scalar mass parameters'),
        'source_coupling_names':sorted(used),'benchmark_config':config,
        'fermion_mass_gev':(float(fermion_masses[0]) if len(set(fermion_masses))==1 else None),
        'fermion_masses_gev':[float(m) for m in fermion_masses],
        'scalar_masses_gev':masses,'n_fermion_flavours':len(y1[0]),
        'ordered_A5_gev_inv':[[_spec_json(complex(z)) for z in row] for row in ordered],
        'symmetrized_C5_gev_inv':[[_spec_json(complex(z)) for z in row] for row in canonical],
        'mnu_gev':[[_spec_json(complex(z)) for z in row] for row in mnu],
        'takagi_singular_values_eV_descending':takagi.tolist(),
        'max_abs_mnu_eV':float(max_abs*1.e9),
        'largest_takagi_mass_eV':float(takagi[0]),
        'mass_matrix_entry_max_vs_takagi_different_quantities':True,
        'symmetry_max_abs_C5_gev_inv':sym_err,
        'note':'Numbers evaluate the original one-loop Matchete InputForm, not a postulated neutrino-mass loop formula.'
    }
    output_dir.mkdir(parents=True,exist_ok=True)
    dest=output_dir/(stem+'_numerical_c5.json')
    dest.write_text(json.dumps(payload,indent=2)+'\n',encoding='utf-8')
    mass_csv=output_dir/(stem+'_mnu_eV.csv')
    with mass_csv.open('w',newline='',encoding='utf-8') as handle:
        writer=csv.writer(handle)
        writer.writerow(['flavour_i','flavour_j','mnu_real_eV','mnu_imag_eV'])
        for i in range(3):
            for j in range(3):
                z=complex(mnu[i][j])*1e9
                writer.writerow([i+1,j+1,format(z.real,'.15g'),format(z.imag,'.15g')])
    return payload,dest


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    sel=parser.add_mutually_exclusive_group(required=True)
    sel.add_argument('--model',choices=tuple(MODELS));sel.add_argument('--all',action='store_true')
    parser.add_argument('--config',type=Path,default=DEFAULT_CONFIG)
    parser.add_argument('--source-dir',type=Path,default=DEFAULT_SOURCE)
    parser.add_argument('--parameters-dir',type=Path,default=DEFAULT_PARAMETERS)
    parser.add_argument('--output-dir',type=Path,default=DEFAULT_OUTPUT)
    args=parser.parse_args()
    config=read_json(args.config)
    summary=[]
    for model in (tuple(MODELS) if args.all else (args.model,)):
        result,path=numerical_model(model,config,args.source_dir,args.parameters_dir,args.output_dir)
        summary.append({'model':model,'file':(str(path.relative_to(ROOT)) if path.is_relative_to(ROOT) else str(path)),'max_abs_mnu_eV':result['max_abs_mnu_eV'],
                        'masses_eV':result['takagi_singular_values_eV_descending']})
        print(f'{model}: illustrative direct one-scale max|mnu_ij|={result["max_abs_mnu_eV"]:.6g} eV; '
              f'largest Takagi mass={result["largest_takagi_mass_eV"]:.6g} eV; '
              f'singular masses={result["takagi_singular_values_eV_descending"]} eV; '
              'physical EFT1 NOT performed',flush=True)
    args.output_dir.mkdir(parents=True,exist_ok=True)
    path=args.output_dir/'summary_final.json'
    path.write_text(json.dumps({'status':STATUS,'models':summary,'physical_c5_certified':False,
                    'sequential_EFT1_performed':False},indent=2)+'\n',encoding='utf-8')
    with (args.output_dir/'summary_final.csv').open('w',newline='',encoding='utf-8') as handle:
        writer=csv.writer(handle)
        writer.writerow(['model','max_abs_mnu_eV','m1_eV_desc','m2_eV_desc','m3_eV_desc','physical_EFT1'])
        for item in summary:
            writer.writerow([item['model'],item['max_abs_mnu_eV'],*item['masses_eV'],False])
    print('Summary:',path)

if __name__=='__main__':main()
