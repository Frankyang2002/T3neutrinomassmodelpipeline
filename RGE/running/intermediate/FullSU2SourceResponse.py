"""Source-resolved *conditional* SU(2) LLSS -> LLHH quartic response.

This independent downstream bridge reads Chat 1's provisional 16-model F-tree
LLSS release. It does NOT guess a GammaCC/Weyl normalization or a physical
Matchete/RGBeta scalar-quartic dictionary. It contracts the actual independent
3-family complex source coordinates against an explicitly supplied diagnostic
fully symmetric real-scalar hhss tensor and checks the *entire* Weinberg
SU(2) tensor, not just a neutrino/H0R pivot.

Master convention:
  L5 = C_ijab psi_i psi_j phi_a phi_b + h.c.
  V4 = lambda_abcd phi_a phi_b phi_c phi_d/4!
  (16 pi^2) beta(C_ijAB)|lambda = sum_{a,b heavy} lambda_ABab C_ijab.

For P5=(1/2) sum_pq C5_pq (nu_p H0 - e_p H+) (nu_q H0 - e_q H+),
  C_(i,p)(j,q)AB = C5_pq T_ijAB,
  T_ijAB = (z_iA z_jB + z_iB z_jA)/4,
  z[nu,H0R]=1/sqrt(2), z[nu,H0I]=i/sqrt(2),
  z[e,H+R]=-1/sqrt(2), z[e,H+I]=-i/sqrt(2).

The output C5 values are mathematical projections of an *uncertified* source
boundary with artificial/diagnostic quartics. They are not physical Wilson
coefficients or evidence that other one-loop mixing classes are absent.
"""
from __future__ import annotations

from hashlib import sha256
from math import isfinite, log, pi, sqrt
from pathlib import Path
from typing import Any, Mapping
import json

import numpy as np

from RGE.running.intermediate.ProvisionalSourceKernel import (
    SOURCE_FORMAT, SourceKernelError, _parse_coordinates, kernel_for_boundary,
)

FORMAT = 't3-eft1-full-su2-source-response-v1'
QUARTIC_FORMAT = 't3-eft1-diagnostic-hhss-quartic-v1'
LOOP = 1.0 / (16.0 * pi**2)


def digest(path: Path) -> str:
    return sha256(Path(path).read_bytes()).hexdigest()


def _require(ok: bool, detail: str) -> None:
    if not ok:
        raise SourceKernelError(detail)


def _encoded(a: np.ndarray) -> dict[str, list]:
    z = np.asarray(a, dtype=complex)
    return {'real': z.real.tolist(), 'imag': z.imag.tolist()}


def _pair_blocks(coordinates: Mapping[tuple[int, ...], complex]
                 ) -> dict[tuple[int, int], np.ndarray]:
    """One heavy scalar unordered pair -> C[p,q,i,j], no orbit double counting."""
    out: dict[tuple[int, int], np.ndarray] = {}
    for (i,p,j,q,a,b), value in coordinates.items():
        if a > b:
            continue
        key = (a,b)
        block = out.setdefault(key, np.zeros((3,3,2,2), dtype=complex))
        block[p-1,q-1,i-1,j-1] = value
    for pair, block in out.items():
        _require(np.array_equal(block, block.transpose(1,0,3,2)),
                 f'Lepton/flavour simultaneous exchange fails for pair {pair}')
    return out


def read_source_blocks(boundary: Mapping[str, Any]) -> tuple[list[tuple[int,dict]], int]:
    """Validate independently against the source-summed polynomial first."""
    validated = kernel_for_boundary(boundary)
    ns = int(validated['real_scalar_dimension'])
    result = []
    sums={}
    absolute_scale={}
    for term in boundary['sources']:
        source_id = term['matchete_source_index']
        entries = _parse_coordinates(term['independent_symmetric_coordinates'],
                                     ns, f'source_{source_id}')
        for key,z in entries.items():
            sums[key]=sums.get(key,0j)+z
            absolute_scale[key]=absolute_scale.get(key,0.0)+abs(z)
        result.append((source_id, _pair_blocks(entries)))
    reference=_parse_coordinates(boundary['summed_independent_symmetric_coordinates'],
                                 ns,'full_su2_source_sum')
    for key in set(sums)|set(reference):
        lhs=sums.get(key,0j)
        rhs=reference.get(key,0j)
        # No fixed absolute floor: do not erase numerically suppressed sources.
        scale=max(absolute_scale.get(key,0.0),abs(rhs),1e-300)
        _require(abs(lhs-rhs)<=1e-10*scale,
                 f'Precision-sensitive source recomposition failed at {key}')
    return result, ns


def weinberg_template() -> np.ndarray:
    """Shape (2,2,4,4); C5=1 gives C[nu nu H0R H0R]=1/4."""
    z = np.zeros((2,4), complex)
    z[0,2], z[0,3] = 1/sqrt(2), 1j/sqrt(2)
    z[1,0], z[1,1] = -1/sqrt(2), -1j/sqrt(2)
    t = (np.einsum('ia,jb->ijab', z, z) +
         np.einsum('ib,ja->ijab', z, z))/4.0
    _require(abs(t[0,0,2,2]-.25)<1e-14, 'Internal Weinberg normalization error')
    _require(np.allclose(t, t.transpose(1,0,2,3)), 'Lepton Ward template asymmetry')
    _require(np.allclose(t, t.transpose(0,1,3,2)), 'Higgs Ward template asymmetry')
    return t


def project_weinberg(beta: np.ndarray, *, atol: float = 1e-300,
                     rtol: float = 1e-8) -> dict:
    """Return projection and residual; only the full tensor can certify a Ward test.

    beta shape: (3,3,2,2,4,4), already symmetric in (p,i)<->(q,j)
    and A<->B. The coefficient of T_ijAB is the Weinberg C5 in *this*
    master polynomial convention (GammaCC still external/unknown).
    """
    beta = np.asarray(beta, complex)
    _require(beta.shape == (3,3,2,2,4,4) and np.isfinite(beta).all(),
             'Expected finite complex beta[p,q,i,j,H_A,H_B]')
    t = weinberg_template()
    norm = float(np.vdot(t,t).real)
    projected = np.einsum('ijab,pqijab->pq', t.conj(), beta)/norm
    residual = beta - np.einsum('pq,ijab->pqijab', projected, t)
    pivot = 4*beta[:,:,0,0,2,2]
    resid_abs = float(np.linalg.norm(residual))
    beta_norm = float(np.linalg.norm(beta))
    resid_rel = resid_abs/max(beta_norm,1e-300) if beta_norm > 0 else 0.0
    pivot_error = float(np.linalg.norm(pivot-projected))
    # No physical validation follows just from this diagnostic check.
    return {
        'projected_C5_beta_16pi2_per_NGamma': _encoded(projected),
        'neutrino_H0R_pivot_beta_16pi2_per_NGamma': _encoded(pivot),
        'weinberg_residual_abs': resid_abs,
        'weinberg_residual_relative': resid_rel,
        'pivot_vs_projection_abs': pivot_error,
        'SU2_Ward_tensor_algebra_pass': bool(resid_abs <= atol + rtol*beta_norm),
        'C5_flavour_symmetric': bool(np.linalg.norm(projected-projected.T) <=
            1e-11*max(float(np.linalg.norm(projected)),1e-300)),
    }


def read_diagnostic_quartic_channels(payload: Mapping[str, Any], boundary: Mapping,
                                     source_sha: str
                                     ) -> tuple[list[tuple[str,np.ndarray]],dict]:
    """Validate labelled independent quartic channels and build symmetric tensors.

    A model can have *distinct* quartic vertices contributing to the same hhss
    real-component entry; those channels are summed, never overwritten. Within
    a single channel, duplicate symmetry permutations must agree exactly.
    """
    ns = int(boundary['real_scalar_dimension'])
    _require(payload.get('format') == QUARTIC_FORMAT and
             payload.get('model') == boundary['model'] and
             payload.get('authority') == 'diagnostic_only',
             'Expected model-matching diagnostic-only hhss quartic input')
    _require(payload.get('source_boundary_sha256') == source_sha,
             'Quartic input must bind to exact fermion-boundary SHA256')
    try:
        scale = float(payload['mu_F_gev'])
    except (KeyError, TypeError, ValueError) as e:
        raise SourceKernelError('Quartic input missing mu_F_gev') from e
    _require(isfinite(scale) and scale == float(boundary['mu_F_gev']),
             'Quartic input scale disagrees with source boundary')
    _require(payload.get('coverage') in ('partial_diagnostic', 'all_hhss_uncertified'),
             'Explicit partial/full-uncertified quartic coverage required')
    if 'quartic_terms' in payload:
        _require('quartic_components' not in payload,
                 'Supply either quartic_terms or quartic_components, never both')
        terms=payload['quartic_terms']
    else:
        terms=[{'label':'unlabelled_diagnostic',
                'components':payload.get('quartic_components')}]
    _require(isinstance(terms,list) and len(terms)>0,
             'Must provide at least one explicit quartic term')
    channels=[]
    labels=set()
    component_count=0
    for t,entry in enumerate(terms):
        _require(isinstance(entry,dict),f'Quartic term {t} must be a mapping')
        label=entry.get('label')
        _require(isinstance(label,str) and label.strip() and label not in labels,
                 f'Quartic term {t} has empty/duplicate label')
        labels.add(label)
        rows=entry.get('components')
        _require(isinstance(rows,list), f'{label}: components must be a list')
        four=np.zeros((4,4,ns-4,ns-4),float)
        used: dict[tuple[int,int,int,int],float]={}
        for n,row in enumerate(rows):
            _require(isinstance(row,dict),f'{label}: quartic row {n} malformed')
            inds=row.get('indices')
            _require(isinstance(inds,list) and len(inds)==4 and all(type(x) is int for x in inds),
                     f'{label}: quartic row {n} must have four 1-based integer indices')
            _require(sum(1 <= x <= 4 for x in inds)==2 and
                     sum(5 <= x <= ns for x in inds)==2,
                     f'{label}: quartic row {n} must have exactly two H and two heavy legs')
            h=sorted(x for x in inds if x<=4)
            heavy=sorted(x for x in inds if x>4)
            key=(*h,*heavy)
            raw=row.get('value')
            _require(type(raw) in (float,int) and isfinite(raw),
                     f'{label}: quartic row {n} must be finite real (real-scalar potential)')
            value=float(raw)
            if key in used:
                _require(value==used[key],f'{label}: conflicting quartic orbit at {key}')
            used[key]=value
            for A,B in ((h[0],h[1]),(h[1],h[0])):
                for a,b in ((heavy[0],heavy[1]),(heavy[1],heavy[0])):
                    four[A-1,B-1,a-5,b-5]=value
        component_count+=len(used)
        channels.append((label,four))
    return channels,{'entry_count':component_count,
                     'quartic_channel_count':len(channels),
                     'quartic_channel_labels':[label for label,_ in channels],
                     'coverage':payload['coverage'],'certified':False,
                     'quartic_real_potential_convention':'V4=lambda_abcd phi_a phi_b phi_c phi_d/4!'}


def read_diagnostic_quartic(payload: Mapping[str, Any], boundary: Mapping,
                            source_sha: str) -> tuple[np.ndarray, dict]:
    """Summed tensor convenience adapter; labelled channels remain available."""
    channels,metadata=read_diagnostic_quartic_channels(payload,boundary,source_sha)
    total=sum((tensor for _,tensor in channels),
              np.zeros_like(channels[0][1]))
    return total,metadata

def _contract_pair(block: np.ndarray, lam: np.ndarray, a: int, b: int
                   ) -> np.ndarray:
    """Contract one source + one unordered pair -> [p,q,i,j,H_A,H_B]."""
    multiplicity = 1 if a==b else 2
    return multiplicity*np.einsum('pqij,AB->pqijAB',block,lam[:,:,a-5,b-5])


def run_diagnostic(boundary: Mapping[str, Any], *, source_sha256: str,
                   quartic: Mapping[str, Any] | None = None,
                   mu_s_gev: float | None = None,
                   include_pair_matrices: bool = True) -> dict:
    """Source-resolved conditional response; never promote a formal input to certified.

    When quartic=None, write factorised SU(2) source blocks, not a numerical beta.
    With a diagnostic quartic, give every source/pair beta and its full SU(2)
    Ward residual as well as the sum. N_Gamma remains unknown and external.
    """
    if boundary.get('format') != SOURCE_FORMAT:
        raise SourceKernelError('Not the accepted provisional fermion-first source')
    blocks, ns = read_source_blocks(boundary)
    mu_f = float(boundary['mu_F_gev'])
    if mu_s_gev is not None:
        try:
            mu_s_gev=float(mu_s_gev)
        except (TypeError,ValueError,OverflowError) as e:
            raise SourceKernelError('mu_S is not a finite positive scale') from e
        _require(isfinite(mu_s_gev) and 0 < mu_s_gev < mu_f,
                 'Require 0 < mu_S < mu_F')
    output = {
        'format': FORMAT,
        'model': boundary['model'],
        'status': 'SourceFactorisedFullSU2_NoQuartic' if quartic is None
                  else 'DiagnosticFullSU2Quartic_NotPhysicalC5',
        'source_boundary_sha256':source_sha256,
        'mu_F_gev':mu_f, 'mu_S_gev':mu_s_gev,
        'source_count':len(blocks),
        'N_GammaCC_to_Weyl':None,
        'absolute_spinor_normalization_certified':False,
        'physical_C5_computed':False,
        'complete_one_loop_operator_inventory_certified':False,
        'one_loop_fermion_hard_matching_included':False,
        'finite_joint_scalar_matching_included':False,
        'full_su2_contraction_convention':'beta[p,q,i,j,H_A,H_B] = sum_{a<=b}(2-delta_ab)*lambda[H_A,H_B,a,b]*C[p,q,i,j,a,b]',
        'weinberg_tensor_convention':'(1/2) C5[p,q] (nu_p H0 - e_p H+) (nu_q H0 - e_q H+)',
    }
    if quartic is None:
        output['source_resolved_blocks']=[{
            'matchete_source_index':sid,
            'scalar_pairs':[{
                'heavy_pair_one_based':list(pair),
                'ordered_pair_multiplicity':1 if pair[0]==pair[1] else 2,
                'LLSS_C_pqij_per_NGamma_GeV_inv':_encoded(block),
            } for pair,block in sorted(src.items())]
        } for sid,src in blocks]
        return output
    channels,info=read_diagnostic_quartic_channels(quartic,boundary,source_sha256)
    output['quartic_input']=info
    base=np.zeros((3,3,2,2,4,4),complex)
    source_outputs=[]
    quartic_totals={label:np.zeros_like(base) for label,_ in channels}
    for sid,pairs in blocks:
        by_source=np.zeros_like(base)
        per_quartic=[]
        for label,lam in channels:
            by_quartic=np.zeros_like(base)
            individual=[]
            for (a,b),block in sorted(pairs.items()):
                term=_contract_pair(block,lam,a,b)
                by_quartic+=term
                if include_pair_matrices and np.any(term):
                    individual.append({'heavy_pair_one_based':[a,b],
                                       'multiplicity':1 if a==b else 2,
                                       **project_weinberg(term)})
            by_source+=by_quartic
            quartic_totals[label]+=by_quartic
            per_quartic.append({'quartic_label':label,
                                'conditional_full_SU2_Weinberg':project_weinberg(by_quartic),
                                'pair_contributions':individual})
        base+=by_source
        source_outputs.append({'matchete_source_index':sid,
                               'zero_at_this_benchmark':not bool(pairs),
                               'conditional_full_SU2_Weinberg':project_weinberg(by_source),
                               'quartic_contributions':per_quartic})
    output['per_quartic']=[{'quartic_label':label,
                           'conditional_full_SU2_Weinberg':project_weinberg(quartic_totals[label])}
                          for label,_ in channels]
    grand=project_weinberg(base)
    _require(grand['C5_flavour_symmetric'],
             'Projected aggregate is not symmetric in flavour')
    output['per_source']=source_outputs
    output['total']=grand
    output['su2_ward_test_interpretation']=('Tensor algebra for supplied quartic only; '
        'a passing test neither certifies source GammaCC sign nor the completeness '
        'of the quartic input. A failing test forbids interpreting the projected '
        'matrix as the Weinberg coefficient.')
    if mu_s_gev is not None:
        lratio=log(mu_s_gev/mu_f)
        output['log_muS_over_muF']=lratio
        delta = (np.asarray(grand['projected_C5_beta_16pi2_per_NGamma']['real']) +
                 1j*np.asarray(grand['projected_C5_beta_16pi2_per_NGamma']['imag']))*LOOP*lratio
        output['conditional_deltaC5_log_per_NGamma_GeV_inv'] = _encoded(delta)
        output['scale_note']=('Fixed-order log only, does not include one-loop hard '
             'F matching, heavy-scalar finite threshold or downstream SMEFT running.')
    return output


def load_json(path: Path) -> dict:
    return json.loads(Path(path).read_text(encoding='utf-8-sig'))
