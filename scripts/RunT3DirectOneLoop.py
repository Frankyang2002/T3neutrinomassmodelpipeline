"""All-16 direct UV one-loop runner (strict by default; provisional opt-in).

No EFT imports. Does not reconstruct or guess Matchete component tensors.
Uncertified inputs can only run with --allow-uncertified and are reported as such.
"""
from __future__ import annotations
import argparse
import hashlib
import json
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from Numerical.direct_one_loop.core import evaluate  # noqa: E402
from Numerical.direct_one_loop.comparison import compare_same_point # noqa: E402

CLASSES = {'A':(1,3,2,(-4,-2,0)), 'B':(2,2,1,(-3,-1,1)),
           'C':(2,2,3,(-3,-1,1)), 'D':(3,1,2,(-2,0,2)),
           'E':(3,3,2,(-4,-2,0,2))}


def models():
    for label,(s1,s2,f,alphas) in CLASSES.items():
        for alpha in alphas:
            sign = 'm' if alpha<0 else 'p'
            yield f'{label}_{alpha:+d}', f'T3_dS1_{s1}_dS2_{s2}_dF_{f}_alpha_{sign}{abs(alpha)}'


def atomic_json(p: Path, obj: object) -> None:
    p.parent.mkdir(parents=True,exist_ok=True)
    temporary=p.with_suffix(p.suffix+'.tmp')
    temporary.write_text(json.dumps(obj,indent=2,allow_nan=False)+'\n',encoding='utf-8')
    temporary.replace(p)


def source_sha256(p: Path) -> str:
    return hashlib.sha256(p.read_bytes()).hexdigest()


def run(project: Path, physical: Path, output: Path, *, allow_uncertified: bool=False,
        retain_individual_terms: bool=True, eft_reference: Path | None=None) -> dict:
    if project.resolve()==output.resolve():
        raise ValueError('Output must be a dedicated subdirectory, not repository root')
    if not (project/'Numerical'/'benchmarks'/'fitted_reference').is_dir():
        raise FileNotFoundError('Missing source frozen benchmark directory under '+str(project))
    output.mkdir(parents=True,exist_ok=True)
    records=[]
    for short,key in models():
        reference=project/'Numerical'/'benchmarks'/'fitted_reference'/f'{key}.json'
        source=physical/f'{short}.json'
        status={'model':short,'source_benchmark':str(reference.relative_to(project)),
                'source_benchmark_sha256':None,'component_input':str(source),
                'status':'input_missing','reason':None,'output_file':None}
        if not reference.is_file():
            status.update(status='benchmark_missing',reason='No frozen reference benchmark')
        else:
            status['source_benchmark_sha256']=source_sha256(reference)
            benchmark=json.loads(reference.read_text(encoding='utf-8-sig'))
            expected=CLASSES[short[0]]
            rep=benchmark['representation']
            if (rep['d_s1'],rep['d_s2'],rep['d_f'],rep['alpha']) != (*expected[:3],int(short[2:])):
                status.update(status='source_mismatch',reason='Frozen representation metadata mismatch')
            elif not source.is_file():
                status['reason']='Source-derived, convention-verified component tensor input not supplied'
            else:
                try:
                    uv=json.loads(source.read_text(encoding='utf-8-sig'))
                    provenance=uv.get('provenance') or {}
                    if provenance.get('benchmark_sha256')!=status['source_benchmark_sha256']:
                        raise ValueError('Physical component input does not match frozen benchmark bytes')
                    if uv.get('model')!=short:
                        raise ValueError('Wrong model id in physical component input')
                    if provenance.get('validated_uv_component_conventions') is not True and not allow_uncertified:
                        raise ValueError('Missing explicit validation of CG, spinor and scalar-potential conventions')
                    answer=evaluate(uv,retain_individual_terms=retain_individual_terms)
                    if not answer['diagnostics']['uv_finite']:
                        raise ValueError('Uncancelled UV pole: incomplete or inconsistent source vertices')
                    answer['status']=('uv_mass_evaluated_source_certified' if
                        provenance.get('validated_uv_component_conventions') is True
                        else 'uv_mass_evaluated_uncertified_source_mapping')
                    if eft_reference is not None:
                        eft_path=eft_reference/(short+'.json')
                        if eft_path.is_file():
                            eft=json.loads(eft_path.read_text(encoding='utf-8-sig'))
                            answer['eft_comparison']=compare_same_point(answer,eft)
                        else:
                            answer['eft_comparison']={'status':'no_same_point_eft_reference'}
                    name=f'{short}_direct_one_loop.json'
                    atomic_json(output/name,answer)
                    status.update(status=answer['status'],output_file=name,reason=None,
                                  masses_ev=answer['total']['neutrino_masses_ev'],
                                  cancellation_ratio=answer['diagnostics']['cancellation_ratio_total_over_sum'])
                except (ValueError,KeyError,RuntimeError,TypeError,ArithmeticError) as err:
                    status.update(status='not_evaluable',reason=str(err))
        records.append(status)
        print(f"{short:6s} {status['status']:45s} {status['reason'] or ''}",flush=True)
    report={'purpose':'independent_UV_direct_one_loop_benchmark',
            'count':16,'evaluated':sum('uv_mass_evaluated' in v['status'] for v in records),
            'no_parameters_fabricated':True,
            'models':records,
            'limitations':('Benchmark numerical inputs and source-generated component CG tensors '
               'are distinct prerequisites. No result is marked certified without an '
               'explicit source-specific convention audit. This calculation does not '
               'include EFT matching or RGE running.')}
    atomic_json(output/'all_16_status.json',report)
    return report


def main() -> int:
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--project-root',type=Path,default=ROOT)
    ap.add_argument('--physical-input',type=Path,default=ROOT/'input'/'lo_direct_one_loop')
    ap.add_argument('--output',type=Path,default=ROOT/'output'/'general_observable'/'lo_direct_one_loop')
    ap.add_argument('--allow-uncertified',action='store_true',help='For explicitly provisional diagnostics only')
    ap.add_argument('--no-individual-terms',action='store_true')
    ap.add_argument('--eft-reference-dir',type=Path,help='Optional same-scale, same-benchmark 3x3 EFT mass-matrix reference directory')
    a=ap.parse_args()
    p=run(a.project_root,a.physical_input,a.output,allow_uncertified=a.allow_uncertified,
          retain_individual_terms=not a.no_individual_terms,eft_reference=a.eft_reference_dir)
    certified = sum(row["status"] == "uv_mass_evaluated_source_certified"
                    for row in p["models"])
    print(f"Direct UV numerical results: {p['evaluated']}/16; "
          f"physically source-certified: {certified}/16; "
          f"status: {a.output/'all_16_status.json'}")
    return 0 if p['evaluated']==16 else 2

if __name__=='__main__':raise SystemExit(main())
