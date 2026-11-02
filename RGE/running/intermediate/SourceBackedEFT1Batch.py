"""Read-only provisional Chat-1 LLSS -> full seven-term EFT1 master diagnostic.

This module does NOT certify GammaCC, the scalar potential, hard matching,
additional dimension-5 operators or physical C5. All computed matrices are
*per unit unknown source-to-master Weyl factor* N_Gamma.

Written as a new downstream bridge; no modifications to source modules.
"""
from __future__ import annotations

from hashlib import sha256
from math import isfinite, log
from pathlib import Path
from typing import Any, Mapping
import json
import numpy as np

from RGE.running.intermediate.ProvisionalSourceKernel import (
    MODELS, SOURCE_FORMAT, SourceKernelError, kernel_for_boundary,
    _parse_coordinates, tag,
)
from RGE.running.intermediate.GeneralEFT1C5Mixing import (
    FORMAT as MASTER_FORMAT, CERTIFICATES, TERMS, LOOP,
    BoundaryError, encode_matrix, load_boundary, calculate_master_beta,
    solve_boundary, project_c5,
)
from RGE.running.intermediate.FullSU2SourceResponse import (
    weinberg_template, project_weinberg,
)

FORMAT = 't3-eft1-source-backed-seven-term-diagnostic-v1'
COUPLING_FORMAT = 't3-eft1-explicit-master-couplings-diagnostic-v1'
CONTRACT_FORMAT = 't3-eft1-source-backed-required-couplings-contract-v1'
RESULT_FORMAT = 't3-eft1-source-backed-seven-term-result-v1'
SOURCE_TO_MASTER = 'FORMAL_UNNORMALIZED_N_GAMMA_EXTERNALLY_UNKNOWN'


def _need(condition: bool, reason: str) -> None:
    if not condition:
        raise SourceKernelError(reason)


def digest(path: Path) -> str:
    return sha256(Path(path).read_bytes()).hexdigest()


def read_json(path: Path) -> dict:
    obj = json.loads(Path(path).read_text(encoding='utf-8-sig'))
    _need(isinstance(obj, dict), f'JSON object required: {path}')
    return obj


def _source_layout(boundary: Mapping[str, Any]) -> list[dict]:
    """Use source offsets, NOT representation guesses (D/E have split fields)."""
    ns = boundary['real_scalar_dimension']
    offsets = boundary['scalar_field_offsets_one_based']
    _need(isinstance(offsets, dict) and len(offsets) >= 2,
          'Both scalar species must be present in the source layout')
    ordered = sorted((int(start), name) for name, start in offsets.items())
    _need(len({name for _, name in ordered}) == len(ordered), 'Repeated scalar name')
    _need(ordered[0][0] == 5 and len({p for p, _ in ordered}) == len(ordered),
          'Scalar offsets must begin after the four Higgs real fields')
    blocks = [{'name': 'H', 'role': 'Higgs', 'indices': [1, 2, 3, 4]}]
    for k, (start, name) in enumerate(ordered):
        stop = ordered[k+1][0] if k+1 < len(ordered) else ns+1
        _need(5 <= start < stop <= ns+1, f'Bad real scalar field offsets: {name}')
        blocks.append({'name': name, 'role': 'heavy', 'indices': list(range(start,stop))})
    _need(sum(len(x['indices']) for x in blocks) == ns,
          'Scalar field blocks leave real scalar components unaccounted for')
    return blocks


def _operator_sources(boundary: Mapping[str, Any]) -> list[dict]:
    ns = boundary['real_scalar_dimension']
    outputs = []
    for src in boundary['sources']:
        sid = src['matchete_source_index']
        # _parse_coordinates checks the full distinct exchange orbit, its
        # representative and relative factors before we pass to the master.
        _parse_coordinates(src['independent_symmetric_coordinates'], ns, f'source PL {sid}')
        entries = []
        for row in src['independent_symmetric_coordinates']:
            i,p,j,q,a,b = row['coordinate']
            entries.append({'indices':[2*(p-1)+i,2*(q-1)+j,a,b],
                            'value':{'real':row['coefficient_real'],
                                     'imag':row['coefficient_imag']}})
        if not entries:
            # The master input schema requires one representative per source;
            # represent a *known zero* source by an explicitly zero component.
            entries=[{'indices':[1,1,5,5], 'value': {'real':0.0,'imag':0.0}}]
        outputs.append({'source_id':f'Matchete_PL_{sid}', 'components':entries})
    return outputs


def required_couplings(boundary: Mapping[str, Any], source_sha: str) -> dict:
    """Ready-now 16-model metadata, not invented physical coupling values."""
    checked = kernel_for_boundary(boundary, source_sha256=source_sha)
    blocks = _source_layout(boundary)
    return {
        'format':CONTRACT_FORMAT, 'model':boundary['model'],
        'status':'SourceBoundaryVerified_CouplingsAndSpinorNormalizationMissing',
        'source_boundary_sha256':source_sha,
        'mu_F_gev':float(boundary['mu_F_gev']),
        'n_real_scalars':boundary['real_scalar_dimension'],
        'n_fermions_minimum':6,
        'lepton_components_master_one_based':[[1,2],[3,4],[5,6]],
        'higgs_real_components_one_based':[1,2,3,4],
        'source_scalar_blocks':blocks,
        'PL_matchete_source_indices':[s['matchete_source_index'] for s in boundary['sources']],
        'known_zero_at_benchmark':[s['matchete_source_index'] for s in boundary['sources']
                                   if s['zero_at_this_benchmark']],
        'nonzero_neutral_quartic_pivot_pairs':checked['distinct_nonzero_quartic_pairs'],
        'required_external_input_format':COUPLING_FORMAT,
        'required_external_keys':['format','model','source_boundary_sha256','mu_F_gev',
          'mu_S_gev','n_fermions','n_real_scalars','renormalization_scheme',
          'master_conventions','coverage_assertion','yukawa_components',
          'quartic_channels','gauge_sectors'],
        'full_input_must_contain':['all EFT1 two-Weyl-one-real-scalar Yukawas',
          'all relevant independent real-scalar quartic channels in V4=lambda phi4/4!',
          'gauge couplings and complete SU2L/U1Y generators for every Weyl/scalar component',
          'all other active fermions required for anomalous dimensions',
          'matched mu_F and mu_S in same MSbar scheme'],
        'physical_source_spinor_factor':None,
        'physical_quartic_normalization_verified':False,
        'full_dim5_operator_inventory_verified':False,
        'hard_fermion_one_loop_matching_included':False,
        'joint_scalar_matching_performed':False,
    }


def _merge_quartics(channels: list[dict], ns: int) -> tuple[list[dict], list[dict]]:
    """Add distinct source invariants; never overwrite two channels on one orbit."""
    _need(isinstance(channels, list) and channels, 'Nonempty quartic_channels required')
    total: dict[tuple[int, ...], float] = {}
    source_inventory: list[dict] = []
    labels: set[str] = set()
    for channel in channels:
        _need(isinstance(channel,dict), 'Quartic channels must be objects')
        name=channel.get('label')
        _need(isinstance(name,str) and name and name not in labels,
              'Quartic channel labels must be distinct and nonempty')
        labels.add(name)
        coordinates:dict[tuple[int,...],float]={}
        rows=channel.get('components')
        _need(isinstance(rows,list),f'{name}: components must be a list')
        for row in rows:
            inds=row.get('indices')
            value=row.get('value')
            _need(isinstance(inds,list) and len(inds)==4 and
                  all(type(x) is int and 1<=x<=ns for x in inds),
                  f'{name}: expected 4 one-based scalar indices')
            _need(type(value) in (int,float) and isfinite(float(value)),
                  f'{name}: real finite lambda_abcd required')
            key=tuple(sorted(inds))
            if key in coordinates:
                _need(np.isclose(coordinates[key],value,rtol=0,atol=0),
                      f'{name}: conflicting equivalent quartic permutations {key}')
            coordinates[key]=float(value)
        for key,value in coordinates.items():
            total[key]=total.get(key,0.0)+value
        source_inventory.append({'label':name,'independent_nonzero_components':
                                 sum(bool(v) for v in coordinates.values()),
                                 'entry_count':len(coordinates)})
    merged=[{'indices':list(key),'value':v} for key,v in sorted(total.items()) if v!=0.]
    return merged,source_inventory


def _validate_couplings(couplings: Mapping[str, Any], boundary: Mapping[str,Any],
                        source_sha: str) -> tuple[dict, list[dict]]:
    model=boundary['model']; ns=boundary['real_scalar_dimension']
    _need(isinstance(couplings,Mapping) and couplings.get('format')==COUPLING_FORMAT,
          f'{model}: explicit diagnostic coupling schema required')
    _need(couplings.get('model')==model and
          couplings.get('source_boundary_sha256')==source_sha,
          f'{model}: coupling model/source hash mismatch')
    _need(couplings.get('renormalization_scheme')=='MSbar',
          f'{model}: MSbar scheme required')
    _need(couplings.get('master_conventions') == {
       'master_wilson':'C_ijab psi_i psi_j phi_a phi_b + hc',
       'weinberg':'(1/2) C5_ij (L_i.H)(L_j.H) + hc',
       'quartic':'V4=lambda_abcd phi_a phi_b phi_c phi_d/4!',
    }, f'{model}: master operator definitions missing or changed')
    _need(couplings.get('n_real_scalars')==ns,
          f'{model}: scalar dimension mismatch')
    nf=couplings.get('n_fermions')
    _need(type(nf) is int and 6<=nf<=100, f'{model}: invalid full Weyl basis dimension')
    _need(couplings.get('mu_F_gev')==boundary['mu_F_gev'],
          f'{model}: mu_F disagrees with source')
    mu_s=couplings.get('mu_S_gev')
    _need(type(mu_s) in (int,float) and isfinite(float(mu_s)) and
          0<mu_s<boundary['mu_F_gev'],f'{model}: invalid mu_S')
    coverage=couplings.get('coverage_assertion')
    _need(isinstance(coverage,dict) and coverage.get('authority')=='diagnostic_only',
          f'{model}: mandatory diagnostic-only authority')
    _need(all(coverage.get(k) is True for k in (
        'all_active_weyl_fermion_components_enumerated',
        'all_yukawa_components_enumerated',
        'all_quartic_channels_enumerated',
        'all_gauge_generators_enumerated_or_explicitly_toy_zero',
    )), f'{model}: explicit coupling coverage assertions required (not independently certified)')
    _need(isinstance(couplings.get('yukawa_components'),list),
          f'{model}: Yukawa components must be explicit (empty is a diagnostic zero)')
    gauge=couplings.get('gauge_sectors')
    _need(isinstance(gauge,list),f'{model}: gauge generators must be explicit')
    if not gauge:
        _need(coverage.get('toy_gauge_couplings_set_to_zero') is True,
              f'{model}: omitting gauge generators is only permitted for toy diagnostics')
    else:
        _need(coverage.get('toy_gauge_couplings_set_to_zero') is False,
              f'{model}: gauge status contradictory')
        _need(len(gauge)==2 and {x.get('name') for x in gauge}=={'SU2L','U1Y'},
              f'{model}: need named SU2L and U1Y sectors')
        for sector in gauge:
            _need(len(sector.get('generators',[]))==(3 if sector['name']=='SU2L' else 1),
                  f'{model}: missing gauge generators')
    _need('quartic_components' not in couplings,
          f'{model}: use labelled quartic_channels; do not double count a summed quartic')
    merged,channels=_merge_quartics(couplings.get('quartic_channels'),ns)
    return {
        'format':MASTER_FORMAT,'model':model,'authority':'synthetic',
        'synthetic_notice':'ENGINE-GATE ONLY: real source data, provisional N_Gamma; ' 
          'coupling dictionary supplied diagnostically; do not treat as physical',
        'scheme':'MSbar','mu_F_gev':float(boundary['mu_F_gev']),
        'mu_S_gev':float(mu_s),'n_fermions':nf,'n_real_scalars':ns,
        'lepton_components':[[1,2],[3,4],[5,6]],
        'higgs_real_components':[1,2,3,4],
        'scalar_field_blocks':_source_layout(boundary),
        'master_lagrangian_convention':'C_ijab psi_i psi_j phi_a phi_b + hc',
        'weinberg_lagrangian_convention':'(1/2) C5_ij (L_i.H)(L_j.H) + hc',
        'certificates':{key:False for key in CERTIFICATES},
        'yukawa_components':couplings['yukawa_components'],
        'quartic_components':merged,
        'gauge_sectors':gauge,
        'operators':_operator_sources(boundary),
    },channels


def _operator_fidelity(context, boundary: Mapping[str,Any]) -> dict:
    """Independent actual source tensor -> mapped master tensor check."""
    ns=boundary['real_scalar_dimension']; nf=context.nf
    _need(len(context.operators)==len(boundary['sources']),
          'Master operator count differs from the PL source inventory')
    model=boundary['model']
    max_errors=[]
    for (id_,observed),src in zip(context.operators,boundary['sources']):
        _need(id_==f"Matchete_PL_{src['matchete_source_index']}",
              f'{model}: source ordering or ID changed')
        expected=np.zeros((nf,nf,ns,ns),complex)
        expanded=_parse_coordinates(src['independent_symmetric_coordinates'],ns,id_)
        for (i,p,j,q,a,b),z in expanded.items():
            expected[2*(p-1)+i-1,2*(q-1)+j-1,a-1,b-1]=z
        error=float(np.max(abs(expected-observed)))
        _need(error<=1e-11*max(float(np.max(abs(expected))),1e-300),
              f'{model}/{id_}: master source coordinates not preserved')
        max_errors.append(error)
    return {'PL_sources_preserved':len(context.operators),
            'absolute_max_tensor_conversion_error':max(max_errors,default=0.),
            'fermion_index_map':'source (i,p) 1-based -> master 2*(p-1)+i',
            'scalar_index_map':'source real a -> master same 1-based a',
            'spinor_normalization_factor_applied':None}


def evaluate(boundary: Mapping[str,Any], source_sha: str,
             couplings: Mapping[str,Any], *, trajectory_points: int=21) -> dict:
    """All seven one-loop anomalous-dimension terms: conditional, per N_Gamma."""
    kernel_for_boundary(boundary,source_sha256=source_sha)
    payload,channels=_validate_couplings(couplings,boundary,source_sha)
    try:
        context=load_boundary(payload,allow_synthetic=True)
        engine=solve_boundary(payload,allow_synthetic=True,
                              trajectory_points=trajectory_points)
    except BoundaryError as exc:
        raise SourceKernelError(f"Master RGE input rejected: {exc}") from exc
    fidelity=_operator_fidelity(context,boundary)
    # Recompute the quartic-only response *independently* from the explicit
    # full master scalar_pair channel (which also contains gauge exchange).
    # Store the full SU2 projection and its independent neutrino pivot.
    qtensor=np.zeros((context.nf,context.nf,context.ns,context.ns),complex)
    term_sums={name:np.zeros((3,3),complex) for name in TERMS}
    for _,op in context.operators:
        qtensor+=np.einsum('abcd,ijcd->ijab',context.quartic,op,optimize=True)
        master=calculate_master_beta(context,op)
        for name in TERMS:
            q,_ward,_sym=project_c5(master[name],context)
            term_sums[name]+=q
    pivot, ward_error, flavour_error=project_c5(qtensor,context)
    leptons=np.array(context.leptons)
    f_indices=leptons.flatten()
    h_indices=np.array(context.higgs)
    sub=qtensor[np.ix_(f_indices,f_indices,h_indices,h_indices)]
    full=np.zeros((3,3,2,2,4,4),complex)
    for p in range(3):
        for q in range(3):
            full[p,q]=sub[2*p:2*p+2,2*q:2*q+2]
    qreport=project_weinberg(full)
    projection=np.asarray(qreport['projected_C5_beta_16pi2_per_NGamma']['real']) +\
               1j*np.asarray(qreport['projected_C5_beta_16pi2_per_NGamma']['imag'])
    _need(np.allclose(pivot,4*full[:,:,0,0,2,2],rtol=1e-12,atol=1e-28),
          'Master-neutral pivot disagrees with source SU2 block')
    scale=float(np.max(abs(pivot)))
    if qreport['SU2_Ward_tensor_algebra_pass']:
        _need(np.linalg.norm(pivot-projection)<=1e-8*max(np.linalg.norm(pivot),1e-300),
              'Passing SU2 Ward test nevertheless disagrees with neutrino pivot')
    for name,mat in term_sums.items():
        got=engine['summed_by_rge_term_beta_C5_16pi2'][name]
        a=np.asarray(got['real'])+1j*np.asarray(got['imag'])
        _need(np.allclose(a,mat,rtol=1e-11,atol=1e-27),
              f'Independent per-operator recomposition failed: {name}')
    return {
        'format':RESULT_FORMAT, 'model':boundary['model'],
        'status':'SourceBackedSevenTermFormalPerNGamma_NotPhysicalC5',
        'source_boundary_sha256':source_sha,
        'source_to_master_fidelity':fidelity,
        'quartic_channel_inventory':channels,
        'mu_F_gev':float(boundary['mu_F_gev']),
        'mu_S_gev':float(context.mu_s),
        'conditional_master_beta_by_rge_term':engine['summed_by_rge_term_beta_C5_16pi2'],
        'conditional_total_beta_C5_16pi2_per_NGamma':engine['total_beta_C5_16pi2'],
        'conditional_total_delta_C5_log_per_NGamma_gev_inv':engine['total_delta_C5_log_gev_inverse'],
        'source_resolved_seven_term_results':engine['operator_results'],
        'trajectory_per_NGamma':engine['trajectory'],
        'quartic_only_SU2_projection':qreport,
        'quartic_only_master_neutrino_pivot_C5_beta':encode_matrix(pivot),
        'quartic_only_master_weinberg_residual_max':ward_error,
        'unphysical_diagnostic_coverage':couplings['coverage_assertion'],
        'physical_C5_computed':False,
        'N_GammaCC_to_master_weyl':None,
        'absolute_spinor_normalization_certified':False,
        'complete_EFT1_coupling_dictionary_independently_certified':False,
        'hard_fermion_one_loop_matching_included':False,
        'other_dim5_operator_classes_audited':False,
        'finite_joint_scalar_matching_included':False,
        'joint_scalar_matching_handoff':{
            'C5_EFT1_log_per_NGamma_gev_inv':engine['total_delta_C5_log_gev_inverse'],
            'full_physical_C5_EFT1_muS':None,
            'ready_for_certified_joint_scalar_matching':False,
            'never_add_this_log_twice':True,
            'note':'Scalar matcher must consume same-muS C5 once, after normalisation.'
        },
    }


def stage22_inventory(stage22: Mapping[str,Any]) -> dict[str,dict]:
    """Read-only geometry inventory. Never interprets symbolic candidates as couplings.

    Stage22 covers extra quartic candidates, not the complete ordinary potential.
    """
    _need(stage22.get('format')=='stage22-aligned-geometry-diagnostics-v1' and
          stage22.get('status')=='ExactConditionalGeometryNotPhysicalCouplingCertification',
          'Unexpected Stage22 diagnostic schema/status')
    for key in ('native_rgbeta_component_basis_verified','physical_prefactors_verified',
                'physical_coupling_map_certified'):
        _need(stage22.get(key) is False,
              'Stage22 is now apparently certified: re-audit the new version, do not assume')
    models=stage22.get('models')
    _need(isinstance(models,dict) and set(models)==set(MODELS),
          'Stage22 diagnostics must describe exactly the sixteen T3 models')
    inventory={}
    for model,blocks in models.items():
        _need(isinstance(blocks,list),'Stage22 model candidate blocks must be a list')
        _need(all(isinstance(b,dict) and b.get('physical_mapping_certified') is False
                  and isinstance(b.get('block'),str) and
                  isinstance(b.get('status'),str) for b in blocks),
              f'{model}: a new or certified quartic geometry record needs inspection')
        inventory[model]={
            'status':'ConditionalGeometryOnly_NotANumericalPotential',
            'extra_quartic_candidate_block_count':len(blocks),
            'candidate_statuses':sorted(set(b['status'] for b in blocks)),
            'physical_prefactors_and_coupling_map_certified':False,
            'ordinary_scalar_quartics_covered_by_this_inventory':False,
            'may_be_used_as_master_quartic_components':False,
        }
    return inventory
