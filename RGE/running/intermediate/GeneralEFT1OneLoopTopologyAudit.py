"""Exact graph-topology enumeration for one-loop psi^2 phi^2 counterterms.

Sources motivating this *structural* (not loop-integral) check:
* Aebischer--Bresciani--Selimovic, arXiv:2512.16890v1, Eq. (4.85).
* Fonseca--Olgoso--Santiago, JHEP 07 (2026) 094, Eq. (22).

Two-component *massless* Weyl rules:
    <psi psi-dagger> and <phi phi> may be contracted;
    <psi psi> and <psi-dagger psi-dagger> vanish.
The vertices considered are C psi psi phi phi, C* psi-dagger psi-dagger
phi phi and y psi psi phi, y* psi-dagger psi-dagger phi.

We enumerate connected, non-self-contracted diagrams with exactly three
internal propagators and external legs psi psi phi phi. A one-loop graph is
1PI iff removing any single internal propagator does not disconnect it.

IMPORTANT: This is NOT a diagrammatic loop integral or a proof that no
EOM/redundant-operator contribution could induce a crossed structure. It
checks the 1PI graph classification independently of Wilson beta formulas.
"""
from __future__ import annotations

from collections import Counter, defaultdict
from dataclasses import dataclass
from itertools import combinations, permutations
from typing import Iterable

VERTICES = {
    'C': ('psi', 'psi', 'phi', 'phi'),
    'Cstar': ('psidag', 'psidag', 'phi', 'phi'),
    'y': ('psi', 'psi', 'phi'),
    'ystar': ('psidag', 'psidag', 'phi'),
}
TARGET_EXTERNAL = Counter({'psi': 2, 'phi': 2})


@dataclass(frozen=True)
class Diagram:
    """One representative contraction; field-identical permutations count separately."""
    species: tuple[str, str, str]
    edges: tuple[tuple[int, int, str], ...]
    external: tuple[tuple[str, ...], ...]
    one_particle_irreducible: bool

    @property
    def external_pattern(self) -> str:
        return '|'.join(f'{self.species[i]}:{",".join(self.external[i]) or "none"}'
                        for i in range(3))

    @property
    def edge_pattern(self) -> str:
        return ','.join(f'{self.species[a]}-{self.species[b]}:{kind}'
                        for a, b, kind in self.edges)


def _connected(edges: Iterable[tuple[int, int, str]]) -> bool:
    visited = {0}
    edges = tuple(edges)
    active = True
    while active:
        active = False
        for u, v, _ in edges:
            if u in visited and v not in visited:
                visited.add(v)
                active = True
            if v in visited and u not in visited:
                visited.add(u)
                active = True
    return len(visited) == 3


def enumerate_three_vertex_diagrams(species: tuple[str, str, str]) -> tuple[Diagram, ...]:
    """Enumerate all chirality-compatible, 1-loop 4-point contractions.

    The three vertices are labelled: a Wilson insertion and two Yukawas.
    No self-contractions are allowed (massless tadpole graphs excluded).
    Distinct identical-field slots are counted, so multiplicities are not
    physical symmetry factors; use topology/external signatures only.
    """
    if len(species) != 3 or species[0] not in ('C', 'Cstar') or any(
            name not in ('y','ystar') for name in species[1:]):
        raise ValueError('Expected (C or Cstar, y or ystar, y or ystar)')
    slots = tuple((i, field, k) for i, name in enumerate(species)
                  for k, field in enumerate(VERTICES[name]))
    by_kind = {field: [slot for slot in slots if slot[1] == field]
               for field in ('psi','psidag','phi')}
    total = Counter(field for _, field, _ in slots)
    nferm = total['psi']-TARGET_EXTERNAL['psi']
    nscal = total['phi']-TARGET_EXTERNAL['phi']
    if nferm < 0 or nferm != total['psidag'] or nscal < 0 or nscal % 2:
        return ()
    nscal //= 2
    if nferm + nscal != 3:  # exactly one loop: E=3, V=3
        return ()

    results = []
    for left in combinations(by_kind['psi'], nferm):
        for right in permutations(by_kind['psidag'], nferm):
            ferm_pairs = tuple(zip(left,right))
            if any(a[0] == b[0] for a,b in ferm_pairs):
                continue
            for scalar_set in combinations(by_kind['phi'], 2*nscal):
                for scalar_pairs in _perfect_scalar_matchings(scalar_set):
                    if any(a[0] == b[0] for a,b in scalar_pairs):
                        continue
                    all_pairs = ferm_pairs + scalar_pairs
                    used = set(h for pair in all_pairs for h in pair)
                    external = tuple(tuple(sorted(slot[1] for slot in slots
                                                  if slot[0] == idx and slot not in used))
                                     for idx in range(3))
                    out = Counter(field for group in external for field in group)
                    if out != TARGET_EXTERNAL:
                        continue
                    edges = tuple(sorted((min(a[0],b[0]), max(a[0],b[0]),
                                          'fermion' if a[1]!='phi' else 'scalar')
                                         for a,b in all_pairs))
                    if not _connected(edges):
                        continue
                    one_pi = all(_connected(edges[:k]+edges[k+1:])
                                 for k in range(len(edges)))
                    results.append(Diagram(species,edges,external,one_pi))
    # Keep only genuinely distinct labelled Wick contractions, not repeated
    # scalar-pair orderings; output topology counts need not be weighted.
    return tuple(results)


def _perfect_scalar_matchings(slots):
    slots = tuple(slots)
    if not slots:
        yield ()
        return
    first, *tail = slots
    for i, second in enumerate(tail):
        remainder = tail[:i]+tail[i+1:]
        for rest in _perfect_scalar_matchings(remainder):
            yield ((first, second),) + rest


def classify() -> dict:
    """Classify distinct external/edge patterns, not symmetry factors."""
    sectors = {}
    for species in (('C','y','y'), ('Cstar','y','y'), ('C','y','ystar')):
        diagrams = enumerate_three_vertex_diagrams(species)
        tag = '_'.join(species)
        patterns = defaultdict(lambda: {'wick_contractions':0, '1PI':None})
        for diagram in diagrams:
            key = (diagram.external_pattern,diagram.edge_pattern)
            patterns[key]['wick_contractions'] += 1
            patterns[key]['1PI'] = diagram.one_particle_irreducible
        sectors[tag] = {
            'species': species,
            'number_of_wick_contractions': len(diagrams),
            'number_1pi_wick_contractions': sum(d.one_particle_irreducible for d in diagrams),
            'patterns': [
                {'external':key[0],'edges':key[1],
                 'is_1PI':value['1PI'], 'wick_contractions':value['wick_contractions']}
                for key,value in sorted(patterns.items())
            ],
        }
    cross = sectors['Cstar_y_y']
    crossed = [p for p in cross['patterns'] if 'Cstar:phi|' in p['external']]
    noncrossed = [p for p in cross['patterns'] if 'Cstar:phi,phi|' in p['external']]
    if not crossed or any(p['is_1PI'] for p in crossed):
        raise RuntimeError('1PR crossed-topology assertion failed')
    if not noncrossed or not any(p['is_1PI'] for p in noncrossed):
        raise RuntimeError('1PI non-crossed-topology assertion failed')
    if sectors['C_y_y']['number_of_wick_contractions']:
        raise RuntimeError('Chiral <psi psi> obstruction failed')
    if not sectors['C_y_ystar']['number_1pi_wick_contractions']:
        raise RuntimeError('Expected mixed-Yukawa 1PI control')
    return {
        'format': 't3-batch34-one-loop-weyl-graph-audit-v1',
        'status': '1PI_TOPOLOGY_FAVOURS_SIX_TERM_REFERENCE_UNVERIFIED_FULL_RGE',
        'assumptions': ['massless Weyl/real scalar propagators',
                        'exactly one C/Cstar insertion and two Yukawa vertices',
                        'external psi psi phi phi',
                        'one loop, connected; self-tadpoles excluded'],
        'sectors': sectors,
        'crossed_external_scalar_attachment_only_1PR': True,
        'noncrossed_Cstar_y_y_has_1PI_triangle': True,
        'unstarred_C_y_y_chirality_forbids_diagram': True,
        'is_complete_counterterm_derivation': False,
        'is_physical_Aminus4_certification': False,
        'caveat': ('Topology alone cannot exclude on-shell/EOM-induced mixing or ' 
                   'settle the September 2026 version-of-record equation.'),
    }
