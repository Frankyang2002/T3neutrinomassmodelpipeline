"""Source-exact B(-1)/C(-1) Matchete F-tree LLSS, at fitted flavour/mixing.

UV -> integrate out all F -> EFT1(S1,S2) -> integrate out all scalars -> SMEFT.

The existing *native* B/C X/Y Matchete source and CG registry are used, not a
single flavour proxy. Each PL term is separately re-evaluated with current
heavy-flavour Yukawas and complex scalar mixing U. PR conjugates are excluded.
Each exact ordered SU(2)/scalar term is mapped to the established symmetric
formal commuting-spinor polynomial. No spinor convention is inferred.

Caveat: the Matchete X/Y F-tree source is evaluated with U(muS); identifying
its source scalar basis with the fitted muS basis is an exploratory fixed-basis
assumption. No complete EFT1 anomalous dimension/hard matching is asserted.
"""
from __future__ import annotations

from hashlib import sha256
import json
import math
from pathlib import Path
import re
from typing import Any, Mapping

import numpy as np

from Numerical.running.GeneralT3SequentialLeadingLog import require

MIXED = ('B_-1', 'C_-1')
EXPECTED_PL = {'B_-1': 10, 'C_-1': 12}
SOURCE_ROOT = 'output/general_observable/ordered_wilson_tensor_v4'
_U = re.compile(
    r'(?P<bar>Bar\[)?Coupling\[(?P<name>T3MEU[12][12]),\s*\{\},\s*0\]'
    r'(?P<close>\])?'
)
# Numeric denominators may contain products such as 2*Sqrt[3]*MF (C_-1).
_DENOM = re.compile(
    r'/\(\s*(?:(?:\d+(?:/\d+)?|Sqrt\[\d+\])\s*\*\s*)+$'
)


def source_tag(model: str) -> str:
    require(model in MIXED, f'not a mixed Majorana model: {model}')
    return model.replace('-', 'm')


def read_native(model: str, root: Path) -> tuple[dict, list[dict], dict]:
    """Read exact saved PL terms and native CGs; never use reference flavour values."""
    directory = Path(root) / SOURCE_ROOT / source_tag(model)
    def read(name):
        path = directory / name
        obj = json.loads(path.read_text(encoding='utf-8-sig'))
        require(isinstance(obj, dict), f'{path}: expected JSON object')
        return obj
    source = read('source_flavour.json')
    cg = read('mixed_native_tree_llss_cg_registry.json')
    ordered = read('ordered_tensor.json')
    require(source.get('model') == cg.get('model') == ordered.get('model') == model,
            f'{model}: mixed source model mismatch')
    require(source.get('format') == 't3-mixed-factorised-LLSS-source-v1' and
            ordered.get('format') == 't3-mixed-native-matchete-ordered-tensor-v1' and
            cg.get('format') == 't3-mixed-matchete-native-cg-registry-v1',
            f'{model}: unexpected native mixed source format')
    require(source.get('source_status_sha256') == ordered.get('source_status_sha256')
            and source.get('source_mass_basis_sha256') == ordered.get('source_mass_basis_sha256')
            and source.get('llss_sha256') == ordered.get('source_synced_matchete_llss_sha256'),
            f'{model}: mixed source and CG-expanded operator provenance differs')
    cg_file = directory / 'mixed_native_tree_llss_cg_registry.json'
    native_bytes = cg_file.read_bytes()
    # Git checkouts on Windows can materialise a text JSON with CRLF bytes.
    digests = {sha256(native_bytes).hexdigest(),
               sha256(native_bytes.replace(b'\r\n', b'\n')).hexdigest()}
    require(ordered.get('source_native_cg_sha256') in digests,
            f'{model}: native CG registry checksum differs from saved operator source')
    require(source.get('mu_F_gev') == ordered.get('mu_F_gev') == 1e5,
            f'{model}: F threshold mismatch')
    terms = source.get('operator_terms')
    require(isinstance(terms, list), f'{model}: no native mixed term inventory')
    pl = [row for row in terms if row.get('chirality') == 'PL']
    require(len(pl) == EXPECTED_PL[model] and
            len(terms) == 2 * EXPECTED_PL[model] and
            len({row['source_index'] for row in pl}) == len(pl),
            f'{model}: incomplete PL/PR operator inventory')
    observed = {int(o['source_index']) for o in ordered.get('operators', [])}
    require(observed == {int(row['source_index']) for row in pl},
            f'{model}: saved ordered source differs from PL terms')
    registry = cg.get('CGRegistry')
    require(isinstance(registry, list) and len(registry) == (2 if model == 'B_-1' else 3),
            f'{model}: native CG registry incomplete')
    return source, pl, cg, ordered


def normalise_mixed_term(raw: str) -> dict:
    """Retain the exact symbolic scalar/spinor/CG source and ordered Yukawas.

    In particular recognise C(-1)'s 2*Sqrt[3]*MF denominator: the existing
    nonmixed normaliser only recognises a single numeric denominator factor.
    """
    from scripts.RunT3GeneralWilsonBoundary import _factors, _MF
    require(isinstance(raw, str) and 'GammaCC' in raw and 'Proj[-1]' in raw
            and 'Proj[1]' not in raw and 'Field[NewFermion' not in raw,
            'expected a PL F-tree LLSS term')
    factors = _factors(raw)
    require(len(factors) == 2 and len({f.heavy_index for f in factors}) == 1
            and len({f.external_index for f in factors}) == 2,
            'not a two-Yukawa Majorana source with one internal heavy index')
    mf = list(_MF.finditer(raw))
    require(len(mf) == 1 and (mf[0].group('bar') is None) ==
            (mf[0].group('close') is None), 'expected one balanced MF factor')
    mass = mf[0]
    require(mass.group('internal') in (None, factors[0].heavy_index),
            'heavy mass/Yukawa flavour indices disagree')
    pre = raw[:mass.start()].rstrip()
    require(pre.endswith('/') or (_DENOM.search(pre) is not None and
            raw[mass.end():].lstrip().startswith(')')),
            'unrecognised inverse heavy-fermion mass denominator')
    spans = sorted([f.span for f in factors] + [mass.span()])
    require(not any(x[1] > y[0] for x,y in zip(spans,spans[1:])),
            'overlapping flavour factors')
    replaced = raw
    for a,b in reversed(spans):
        replaced = replaced[:a] + '1' + replaced[b:]
    require(not re.search(r'Coupling\[(?:MF|y1|y2|y1Conj|y2Conj)\b', replaced),
            'unreplaced Yukawa or MF source factor')
    return {
        'operator_template_unit_flavor_inputform': replaced,
        'factors': [dict(name=f.name,conjugated=f.conjugated,
                        external_index=f.external_index,
                        heavy_index=f.heavy_index) for f in factors],
        'flavour_matrix_axes_external_indices': [f.external_index for f in factors],
        'symmetric_flavour_combination_assumed': False,
    }


def insert_current_mixing(template: str, U: np.ndarray) -> tuple[str,tuple[str,...]]:
    """Replace only native T3MEUij vertices by current fitted U, preserving Bar."""
    U = np.asarray(U, complex)
    require(U.shape == (2,2) and np.isfinite(U).all(), 'expected finite 2x2 U')
    seen=[]
    def substitute(match):
        require((match.group('bar') is None) == (match.group('close') is None),
                'unbalanced Matchete Bar[U]')
        key=match.group('name')
        z=complex(U[int(key[-2])-1,int(key[-1])-1])
        if match.group('bar') is not None:
            z=z.conjugate()
        seen.append(key)
        return f'({z.real:.17g}+({z.imag:.17g})*I)'
    result=_U.sub(substitute,template)
    require(bool(seen) and 'Coupling[T3MEU' not in result and
            'Bar[Coupling[T3MEU' not in result and 'Coupling[' not in result,
            'missing/unreplaced current scalar-mixing spurion')
    return result,tuple(seen)


def archived_reference_reproduction(model: str, pl: list[dict], cg: dict,
                                   ordered: dict) -> float:
    """Regenerate the archived *ordered* source from its own native CGs.

    This is a regression and provenance check, not the new fitted flavour
    evaluation. Dynamic matching below starts from original symbolic terms.
    """
    from Numerical.running.GeneralT3OrderedTensor import project_ordered_term
    from Numerical.running.GeneralT3SequentialLeadingLog import decode
    saved = {int(x['source_index']): x for x in ordered['operators']}
    worst = 0.0
    for term in pl:
        idx = int(term['source_index'])
        observed = project_ordered_term(term, cg['CGRegistry'], model)
        reference = saved[idx]
        def components(obj):
            return {(tuple(x['lepton_su2']), tuple(x['scalar_real'])):
                    decode(x['ordered_flavour_coefficient_gev_inv'])
                    for x in obj['components']}
        a,b = components(observed),components(reference)
        require(set(a)==set(b), f'{model} PL {idx}: archived CG component keys changed')
        scale=max((float(np.linalg.norm(x)) for x in b.values()),default=0.0)
        difference=max((float(np.linalg.norm(a[k]-b[k])) for k in a),default=0.0)
        residual=difference/max(scale,1e-35)
        require(residual < 5e-10,
                f'{model} PL {idx}: cannot reproduce saved native CG components '
                f'(relative discrepancy {residual:.3g})')
        worst=max(worst,residual)
    return worst


def project_source(model: str, source: dict, pl: list[dict], cg: dict,
                   masses: np.ndarray, yukawas: Mapping[str,np.ndarray],
                   U: np.ndarray):
    from Numerical.running.GeneralT3OrderedTensor import project_ordered_term
    from Numerical.running.GeneralT3DynamicFermionWilson import _dense_operator
    from scripts.RunT3GeneralWilsonBoundary import flavour_matrix, encode_matrix
    operator_rows=[]; rows=[]; source_metadata=[]
    for saved in pl:
        row=normalise_mixed_term(saved['source_term_inputform'])
        row['source_index']=int(saved['source_index'])
        exact,seen=insert_current_mixing(row['operator_template_unit_flavor_inputform'],U)
        row['operator_template_unit_flavor_inputform']=exact
        # Validate that every distinct source U spurion survives dynamically.
        require(set(seen)==set(saved['spurion_names']),
                f'{model} PL {row["source_index"]}: source mixing inventory differs')
        assert_factors = saved['factors']
        require(row['factors']==assert_factors and
                row['flavour_matrix_axes_external_indices'] ==
                saved['flavour_matrix_axes_external_indices'],
                f'{model} PL {row["source_index"]}: source flavour inventory differs')
        W=flavour_matrix(row,masses,yukawas)
        row['ordered_flavour_matrix_gev_inv']=encode_matrix(W)
        operator=project_ordered_term(row,cg['CGRegistry'],model)
        operator_rows.append(operator)
        rows.append(row)
        source_metadata.append({'source_index':row['source_index'],
                 'yukawa_pair':[x['name'] for x in row['factors']],
                 'conjugated':[x['conjugated'] for x in row['factors']],
                 'scalar_mixings':list(seen),
                 'ordered_bilinear_frobenius':float(np.linalg.norm(W))})
    by, indices, polynomial_audit = _dense_operator(operator_rows,rows,model)
    return by,indices,source_metadata,polynomial_audit


def mixed_exact_boundary(model: str, uv_state: Any, scalar_state: Mapping,
                         *, repo_root: Path, mu_f: float=1e5):
    """Actual fitted Majorana flavour, all source CGs, and current U(muS)."""
    from Numerical.matching.GeneralT3FlavourC5 import prepare_heavy_flavour
    from Numerical.experimental.GeneralT3MixedScalarMassBasis import MixedDoubletMassBasis
    from Numerical.joint_scalar_matching.fermion_boundary_bridge import SourceBoundary
    from Numerical.joint_scalar_matching.special_mass_v6 import expected_special_labels

    require(model in MIXED and abs(float(uv_state.mu_gev)-mu_f) <= 1e-6,
            'invalid mixed model or F matching scale')
    src,pl,cg,archived=read_native(model,repo_root)
    archived_residual=archived_reference_reproduction(model,pl,cg,archived)
    m1=float(scalar_state['mS1Sq']);m2=float(scalar_state['mS2Sq'])
    mu=scalar_state['muT3GeneralS1S2']
    if isinstance(mu,Mapping) and {'real','imag'} <= set(mu):
        mu=complex(float(mu['real']),float(mu['imag']))
    else:
        mu=complex(mu)
    basis=MixedDoubletMassBasis.build(model,m1,m2,mu)
    require(not basis.tachyonic and not basis.degenerate and
            float(min(basis.eigenvalues_sq))>0,
            'mixed scalar X/Y mass basis invalid')
    U=np.asarray(basis.U,complex)
    masses,ys,heavy_audit=prepare_heavy_flavour(uv_state,{'y1','y2'})
    by,indices,terms,audit=project_source(model,src,pl,cg,masses,ys,U)
    C=np.sum(by,axis=0)
    require(C.shape==(3,3,2,2,8,8) and np.isfinite(C).all(),
            'invalid source-exact mixed LLSS dimension')
    release_path=Path(repo_root)/SOURCE_ROOT/source_tag(model)/'source_flavour.json'
    signature=sha256(release_path.read_bytes()+C.tobytes()).hexdigest()
    source=SourceBoundary(model,release_path,signature,float(mu_f),C,by,
               indices,expected_special_labels(model),False,
               'MixedNativeCG_ExactMajoranaFlavour_MuSMixingCandidate')
    info={'status':'source_exact_mixed_native_CG_candidate',
          'source_path':str(release_path),
          'source_sha256':sha256(release_path.read_bytes()).hexdigest(),
          'source_count':len(pl),'source_indices':list(indices),
          'archived_native_CG_reproduction_relative_residual':archived_residual,
          'distinct_yukawa_pairs':sorted({tuple(x['yukawa_pair']) for x in terms}),
          'mixing_U_muS':{'real':U.real.tolist(),'imag':U.imag.tolist()},
          'heavy_mass_eigenvalues_gev':np.asarray(masses,float).tolist(),
          'heavy_basis_diagnostics':str(heavy_audit),
          'operator_details':terms, **audit,
          'muF_to_muS_scalar_basis_identification_assumed':True,
          'spinor_normalization_certified':False,
          'independent_hard_matching_computed':False}
    return source, info
