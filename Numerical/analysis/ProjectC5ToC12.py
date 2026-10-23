"""Exact representation-dependent Weinberg -> mixed-LLSS operator projection.

Uses repository Matchete CG seed, EFT1 quartic seed and RGBeta metadata.
No fitting, numerical benchmarking or RGEs generated from scratch.

The one-loop reverse mixing comes from the quartic scalar-pair term.
Gauge / Yukawa / anomalous-dimension terms cannot replace both Higgs
external legs of the input Weinberg tensor by BSM scalar legs because
SM Yukawa and gauge generators conserve scalar blocks. The script only
reports the quartic-induced reverse mixing and never labels it 'full' C12 RGE.
"""
from __future__ import annotations
import argparse
import csv
import json
from pathlib import Path
import sympy as sp

from RGE.general.WilsonTensorRGE import WilsonRGEInputs, scalar_pair_term
from RGE.matching.WeinbergTensorAdapter import (
    build_weinberg_wilson_tensor, SparseWilsonLookup,
)
from RGE.running.intermediate.ScalarOnlyWilsonTensorRGE import (
    build_scalar_only_rge_context,
)

ROOT = Path(__file__).resolve().parents[2]
DATA_ROOT = ROOT / 'output/full/comparison/smallY_smallL'
OUT = ROOT / 'Reports/output/analytical/GroupFactors/C12_reverse'
PARTS = ('scalar_pair',)


def simplify(x):
    return sp.simplify(sp.expand(x))


def as_text(x):
    return sp.sstr(sp.factor(simplify(x)))


def canonical(i, j, a, b):
    return min(i, j), max(i, j), min(a, b), max(a, b)


def weight(k):
    i, j, a, b = k
    return (1 if i == j else 2) * (1 if a == b else 2)


def scalar_components(context):
    f = context.fermion_basis
    leptons = (f.global_index('L', 1), f.global_index('L', 2))
    return [(i, j, a, b)
            for ni, i in enumerate(leptons)
            for j in leptons[ni:]
            for a in range(5, context.model.total_real_scalar_dimension + 1)
            for b in range(a, context.model.total_real_scalar_dimension + 1)]


def seed_basis(context, components):
    """Use matched CG tensor, with its common Yukawa/MF prefactor removed."""
    values = {k: sp.sympify(context.coefficient[k]) for k in components}
    symbols = sorted(set().union(*(x.free_symbols for x in values.values())), key=str)
    required = {'y1', 'y2', 'MF'}
    present = {str(x) for x in symbols}
    if not required.issubset(present):
        raise ValueError(f'Cannot identify matched seed common factor: {present}')
    sub = {x: sp.S.One for x in symbols if str(x) in required}
    values = {k: simplify(x.subs(sub)) for k, x in values.items()}
    if not any(x != 0 for x in values.values()):
        raise ValueError('Matched C12 basis vanished after prefactor stripping')
    remaining = set().union(*(x.free_symbols for x in values.values()))
    if remaining:
        raise ValueError(f'C12 basis unexpectedly retains symbols: {remaining}')
    return values


def project(values, basis, components):
    """Exact Hermitian projection using independent symmetric tensor entries."""
    den = simplify(sum(weight(k) * sp.conjugate(basis[k]) * basis[k]
                       for k in components))
    if den == 0:
        raise ValueError('Zero C12 tensor norm')
    num = simplify(sum(weight(k) * sp.conjugate(basis[k]) * values[k]
                       for k in components))
    alpha = simplify(num / den)
    residual = {','.join(map(str, k)): as_text(values[k] - alpha * basis[k])
                for k in components
                if simplify(values[k] - alpha * basis[k]) != 0}
    return alpha, residual


def find_models(root, requested):
    names = sorted(p.name for p in root.iterdir() if p.is_dir())
    if requested:
        names = [n for n in names if requested in n]
    if not names:
        raise FileNotFoundError('No matching model artifacts under '+str(root))
    return names


def locate(folder):
    choices = {
        'wilson': list(folder.rglob('eft1_after_F_wilson_seed.json')),
        'quartic': list(folder.rglob('eft1_after_F_scalar_quartic_seed.json')),
        'rgbeta': list(folder.rglob('eft1_rgbeta_rge.json')),
    }
    for kind, hits in choices.items():
        if len(hits) != 1:
            raise FileNotFoundError(f'{folder}: expected one {kind}, got {len(hits)}')
    return {k: v[0] for k, v in choices.items()}


def evaluate(model, base):
    paths = locate(base / model)
    context = build_scalar_only_rge_context(paths['wilson'], paths['quartic'], paths['rgbeta'])
    components = scalar_components(context)
    basis = seed_basis(context, components)
    w5 = build_weinberg_wilson_tensor(context.model, context.fermion_basis, sp.S.One)
    coefficient = SparseWilsonLookup(w5)
    # The quartic is a linear functional of its scalar couplings.
    inputs = WilsonRGEInputs(
        fermion_dimension=context.fermion_basis.dimension,
        yukawa=context.inputs.yukawa,
        quartic=context.quartic,
        gauge_sectors=context.inputs.gauge_sectors,
    )
    values = {k: simplify(scalar_pair_term(context.model, inputs, k, coefficient))
              for k in components}
    alpha, residual = project(values, basis, components)
    # Exact separate coefficients in quartic variables. We require linearity.
    couplings = sorted(alpha.free_symbols, key=str)
    coefficients = {}
    for coupling in couplings:
        derivative = simplify(sp.diff(alpha, coupling))
        if derivative != 0:
            coefficients[str(coupling)] = as_text(derivative)
    leftover = simplify(alpha - sum(sp.Symbol(k)*sp.sympify(v) for k,v in coefficients.items()))
    if leftover != 0:
        raise ValueError(f'Unaccounted nonlinear/free part of reverse mixing: {leftover}')
    return {'model': model, 'representation': context.metadata,
            'coefficient_convention': 'C5=1; C12 matched CG with y1=y2=MF=1',
            'beta_C12_from_C5': as_text(alpha),
            'coefficients': coefficients,
            'residual_component_count': len(residual),
            'residuals': residual,
            'c12_basis_nonzero': sum(x != 0 for x in basis.values()),
            'c5_input_nonzero': len(w5),
            'paths': {k: str(v.relative_to(ROOT)) for k,v in paths.items()}}


def write_outputs(results, output):
    output.mkdir(parents=True, exist_ok=True)
    (output/'c5_to_c12.json').write_text(json.dumps(results, indent=2)+'\n',encoding='utf-8')
    with (output/'c5_to_c12.csv').open('w', newline='', encoding='utf-8') as file:
        writer=csv.DictWriter(file,fieldnames=['model','coupling','factor','residual_component_count'])
        writer.writeheader()
        for result in results:
            for coupling, factor in result['coefficients'].items():
                writer.writerow({'model':result['model'],'coupling':coupling,'factor':factor,
                                 'residual_component_count':result['residual_component_count']})
    tex = [r'\documentclass[11pt]{article}',r'\usepackage[margin=18mm]{geometry}',
           r'\usepackage{amsmath,amssymb,longtable,booktabs}',
           r'\begin{document}',r'\section*{Representation-specific $C_5\to C_{12}$ mixing}',
           r'Computed from the exact scalar quartic tensor contraction in the one-loop $\psi^2\phi^2$ RGE.',
           r'\[16\pi^2\beta_{C_{12}}\supset\gamma_{12,5}C_5.\]',
           r'Nonzero projection residuals mean that the matched $C_{12}$ direction is not closed under this mixing.',
           r'\begin{longtable}{p{0.33\linewidth}p{0.42\linewidth}p{0.15\linewidth}}',
           r'\toprule Model & $\gamma_{12,5}$ (quartic) & Residual components\\\midrule']
    for result in results:
        name=result['model'].replace('_',r'\_')
        expr=sp.latex(sp.sympify(result['beta_C12_from_C5']))
        tex.append(r'\texttt{\small '+name+r'} & $\displaystyle '+expr+r'$ & '+str(result['residual_component_count'])+r'\\')
    tex += [r'\bottomrule\end{longtable}',r'\end{document}']
    (output/'C5_to_C12.tex').write_text('\n'.join(tex)+'\n',encoding='utf-8')


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--model',default='',help='Substring of representation directory')
    p.add_argument('--input-root',type=Path,default=DATA_ROOT)
    p.add_argument('--output',type=Path,default=OUT)
    args=p.parse_args()
    results=[]
    for model in find_models(args.input_root,args.model):
        result=evaluate(model,args.input_root)
        results.append(result)
        print(model, 'gamma12,5 =',result['beta_C12_from_C5'],
              'residuals =',result['residual_component_count'],flush=True)
    write_outputs(results,args.output)
    print('Saved',len(results),'representations to',args.output)

if __name__=='__main__':
    main()
