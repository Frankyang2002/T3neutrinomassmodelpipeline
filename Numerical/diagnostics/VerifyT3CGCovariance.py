"""Verify infinitesimal SU(2) covariance of exported T3 CG arrays.

Inputs: Reports/output/full/multiplicity/component_support_audit.json
        Reports/output/full/multiplicity/matchete_generator_probe.json
Output: Reports/output/full/multiplicity/cg_covariance_audit.json

No loop multiplicity or coupling compensation is inferred.
"""
from __future__ import annotations
import itertools
import json
from pathlib import Path
import sympy as sp

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / 'Reports/output/full/multiplicity'
FLAGS = {
    'T3Y1CG': (True, True, False),  # Bar[L] F^c S1
    'T3Y2CG': (True, False, True),  # Bar[L] F Bar[S2]
    'T3MixCG': (False, False, False, True),  # H H S1 Bar[S2]
}

def generators(probe: dict) -> dict[int, list[sp.Matrix]]:
    out = {1: [sp.zeros(1) for _ in range(3)]}
    for rep in probe['representations']:
        d = int(rep['dimension'])
        raw = [sp.Matrix([[sp.sympify(entry) for entry in row] for row in matrix])
               for matrix in rep['generators_input_form']]
        if len(raw) != 3 or any(m.shape != (d, d) for m in raw):
            raise ValueError(f'Invalid generator matrices for dimension {d}')
        # Matchete original order: G1, G2(Cartan), G3. Positive orientation: G1,G3,G2.
        gs = [raw[i] for i in (0, 2, 1)]
        for a, b, c in ((0, 1, 2), (1, 2, 0), (2, 0, 1)):
            if sp.simplify(gs[a]*gs[b] - gs[b]*gs[a] - sp.I*gs[c]) != sp.zeros(d):
                raise ValueError(f'SU(2) algebra failed for dimension {d}')
        out[d] = gs
    for d in (2, 3):
        if d not in out:
            raise ValueError(f'Missing SU(2) dimension {d}')
    return out

def tensor_residual(vertex: dict, flags: tuple[bool, ...], gs: dict[int, list[sp.Matrix]]) -> list[dict]:
    dimensions = vertex['dimensions']
    if len(flags) != len(dimensions):
        raise ValueError('Mismatched tensor dimensions and flags')
    values = {tuple(i-1 for i in component['indices']): sp.sympify(component['cg'])
              for component in vertex['components']}
    failures = []
    for a in range(3):
        for idx in itertools.product(*(range(d) for d in dimensions)):
            residual = sp.S.Zero
            for axis, d in enumerate(dimensions):
                T = gs[d][a]
                # For a field in the conjugate representation: Tbar = -conjugate(T).
                if flags[axis]:
                    T = -T.conjugate()
                for j in range(d):
                    other = idx[:axis] + (j,) + idx[axis+1:]
                    residual += T[j, idx[axis]] * values.get(other, sp.S.Zero)
            residual = sp.simplify(residual)
            if residual != 0:
                failures.append({'generator': a+1, 'component_1_based': [i+1 for i in idx],
                                 'residual': str(residual)})
    return failures

def audit(components: dict, probe: dict) -> dict:
    gs = generators(probe)
    cases = []
    for case in components['cases']:
        vertices = {}
        for name, tensor in case['vertices'].items():
            failures = tensor_residual(tensor, FLAGS[name], gs)
            vertices[name] = {'covariant': not failures, 'nonzero_residual_count': len(failures),
                              'residual_examples': failures[:5]}
        cases.append({'model': case['model'], 'vertices': vertices,
                      'all_covariant': all(item['covariant'] for item in vertices.values()),
                      'verified_loop_multiplicity': None,
                      'compensation_factor': None})
    return {'status': 'cg_generator_covariance_verified' if all(x['all_covariant'] for x in cases)
            else 'cg_generator_covariance_failed', 'count': len(cases), 'models': cases,
            'convention': 'For each vertex: sum_axis (T_axis)^T C=0; conjugated fields use -conjugate(T).',
            'warning': 'CG covariance does not determine the full propagator contraction or physical loop multiplicity.'}

def main() -> None:
    components = json.loads((OUT/'component_support_audit.json').read_text(encoding='utf-8'))
    probe = json.loads((OUT/'matchete_generator_probe.json').read_text(encoding='utf-8'))
    result = audit(components, probe)
    destination = OUT/'cg_covariance_audit.json'
    destination.parent.mkdir(parents=True, exist_ok=True)
    destination.write_text(json.dumps(result, indent=2)+'\n', encoding='utf-8')
    print(f"Models checked: {result['count']}; covariant: {sum(x['all_covariant'] for x in result['models'])}")
    print(f'Output: {destination}')
    if result['status'] != 'cg_generator_covariance_verified':
        raise SystemExit(1)
if __name__ == '__main__':
    main()
