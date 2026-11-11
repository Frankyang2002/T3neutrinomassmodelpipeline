"""Opt-in 14-model exploratory nonmixed full-quartic simultaneous EFT1 runner.

Preserves existing fitted baseline UV trajectory and matching, but sets the
NEW (nonbaseline) quartics to zero at muF. Their EFT1 evolution uses fresh
full-quartic projected RGBeta betas and is simultaneous with the existing
seven-term Wilson master tensor. Matchete invariant order is provisionally
identified with registry RGBeta coupling order (NOT certified).

This is NOT full UV extra-quartic running: extras are initialized to zero at
muF even in the U=1 configurations. It produces an explicitly exploratory
16-case U/E/W/S ablation using existing scalar matching and SMEFT code.
"""
from __future__ import annotations

from dataclasses import replace
import hashlib
import json
from pathlib import Path
import sys

import numpy as np

from Numerical.experimental.GeneralT3ScalarRGState import (
    GeneralT3ScalarRGState, pack_general_scalar_state,
    unpack_general_scalar_state, load_registry,
)
from Numerical.orchestration.GeneralRunningBenchmark import _state_values
from Numerical.running.T3ExploratoryFullQuartic import build_overlay, run_with_adapter
from scripts.T3ExtendedSimultaneousAdapters import ordinary_adapter
from scripts import RunGeneralT3EFT1Quartics as full_beta
from scripts import RunT3SimultaneousEFT1All16 as general

ROOT = Path(__file__).resolve().parents[1]
REGISTRY = ROOT/'RGE/general/generated/general_rgbeta_quartic_registry.json'
CANDIDATES = ROOT/'output/general_t3_stage22/matchete_cartesian_candidates.json'
OUTPUT = ROOT/'output/studies/general_scalar_extended/experimental_simultaneous_eft1_nonmixed_fullq'
LABEL = 'Exploratory Matchete I_k = RGBeta block-order; extra quartics zero at muF'
# These independent holomorphic scalar bilinears are not in the currently
# generated full-quartic RGBeta beta inventory. Holding them constant is an
# explicit exploratory approximation, NOT a calculated beta function.
FROZEN_QUADRATICS = frozenset({'muT3GeneralS1Sq', 'muT3GeneralS2Sq',
                                'muT3GeneralS1S2'})


def require_sources():
    if not REGISTRY.is_file() or not CANDIDATES.is_file():
        raise FileNotFoundError('Missing general quartic registry or Batch 6 Cartesian candidates')
    reg = load_registry(REGISTRY)
    candidates = json.loads(CANDIDATES.read_text(encoding='utf-8-sig'))
    if candidates.get('format') != 't3-matchete-candidate-cartesian-quartics-v1':
        raise ValueError('Incorrect Cartesian candidate evidence format')
    if set(reg['models']) != set(candidates['models']):
        raise ValueError('Cartesian candidates do not cover all 16 models')
    return reg


def extra_boundary(model, eft_start, registry):
    """Explicit muF-zero approximation for all NEW RGBeta scalar quartics."""
    entry = registry['models'][model]
    extra = {row['name']: 0j for row in entry['new_rgbeta_symbols']}
    return GeneralT3ScalarRGState(eft_start, extra).validated(registry)


def checked_eft1_betas_with_frozen_quadratics(model, state, registry,
                                               original_path, effective_path):
    """Complete only absent *quadratic* beta keys via explicit frozen-RG ansatz.

    Preserve native RGBeta quartic and every other beta without modification.
    Write the approximation to a SEPARATE file; leave native evidence intact.
    """
    _, layout = pack_general_scalar_state(state, registry)
    native = json.loads(Path(original_path).read_text(encoding='utf-8-sig'))
    raw = native.get('report_betas')
    if native.get('status') != 'Success' or not isinstance(raw, dict):
        raise ValueError(f'{model}: malformed native RGBeta beta payload')
    missing = set(layout.names) - set(raw)
    excess = set(raw) - set(layout.names)
    if excess or not missing <= FROZEN_QUADRATICS:
        raise ValueError(f'{model}: unexpected missing/extra RGBeta beta keys: '
                         f'missing={sorted(missing)}, excess={sorted(excess)}')
    if not missing:
        # Delegate to unchanged repository validator for native complete payload.
        beta, names = full_beta.checked_payload(model, state, registry,
                                               Path(original_path))
        return beta, names, ()
    effective = dict(native)
    effective['report_betas'] = dict(raw)
    effective['report_betas'].update({name: '0' for name in sorted(missing)})
    effective['exploratory_frozen_quadratic_beta_names'] = sorted(missing)
    effective['frozen_quadratic_approximation'] = (
        'd(mu_complex)/dln(mu)=0 for listed absent independent quadratic betas; '
        'not supplied by native RGBeta; all native quartic betas unchanged')
    effective['native_beta_source'] = str(original_path)
    target = Path(effective_path)
    target.parent.mkdir(parents=True, exist_ok=True)
    target.write_text(json.dumps(effective, indent=2, allow_nan=False) + '\n',
                      encoding='utf-8')
    # The existing full RGBeta validator still checks every dimension, metadata
    # item and numerical beta expression in the augmented vector.
    beta, names = full_beta.checked_payload(model, state, registry, target)
    return beta, names, tuple(sorted(missing))


class HHSSOverlay:
    """Thin adapter for Batch 7 overlay with instantaneous HHSS input."""
    def __init__(self, inner, slot):
        self.inner = inner
        self._hhss_slot = slot
        self.model = inner.model
        self.symbols = inner.symbols
        self.source = inner.source
        self.mapping_assumption = inner.mapping_assumption
    def evaluate(self, state):
        return self.inner.evaluate(state)
    def validate_hhss(self, state, expected, *, rtol=1e-9):
        self._hhss_slot['value'] = expected
        return self.inner.validate_hhss(state, expected, rtol=rtol)


def run_joint(model, eft_start, source_c, q_basis, mu_s, registry, beta,
              *, eft_on):
    state0 = extra_boundary(model, eft_start, registry)
    entry = registry['models'][model]
    required = {row['name'] for row in entry['new_rgbeta_symbols']}
    required |= {row['existing_name'] for row in entry['preexisting_rgbeta_symbols']}
    slot = {'value': None}
    raw = build_overlay(model,CANDIDATES,REGISTRY,general.general_quartic,
                        lambda _state: slot['value'], required_quartics=required)
    overlay = HHSSOverlay(raw,slot)
    adapter, vec0 = ordinary_adapter(model,state0,beta,registry,source_c,q_basis,mu_s)
    # This exploratory overlay's explicit directions cover the additional
    # non-HHSS tensors; the ordinary baseline and HHSS are provided separately.
    adapter = replace(adapter,required_quartics=overlay.symbols)
    old_unpack = adapter.unpack
    def unpack(vec,mu):
        data,q = old_unpack(vec,mu)
        slot['value'] = q
        return data,q
    adapter = replace(adapter,unpack=unpack)
    master0 = general.general_embed(model,source_c)
    combined = run_with_adapter(adapter,vec0,master0,overlay,eft_on=eft_on)
    final = unpack_general_scalar_state(combined['end_ren_vector'],
                    pack_general_scalar_state(state0,registry)[1],registry,
                    uv=False,mu_gev=float(mu_s))
    combined['end_state'] = final.baseline
    combined['extra_endpoint'] = dict(final.extra)
    q0 = q_basis.tensor(_state_values(eft_start) | dict(state0.extra))
    from Numerical.joint_scalar_matching.core import project_weinberg
    first = general.seven.project_leptonic_weinberg(combined['initial_scalar_pair'])[0]
    direct = np.einsum('fgijab,ABab->fgijAB',source_c,q0,optimize=True)
    projected = project_weinberg(direct,strict=False)[0]
    parity = float(np.linalg.norm(first-projected)/max(np.linalg.norm(projected),1e-300))
    if parity > 1e-8:
        raise ValueError(f'{model}: initial scalar-pair parity failed: {parity:.3e}')
    combined['initial_scalar_quartic_source_parity_rel'] = parity
    combined['ren_layout'] = pack_general_scalar_state(state0,registry)[1]
    combined['extended_state'] = final
    combined['quartic_mapping'] = LABEL
    return combined, overlay


def _nonmixed_full(model,saved,canonical_run_id,quartic_source,registry,
                   *,spinor_factor,wolframscript,vev):
    from scripts import RunT3GeneralAll16FromReference as g
    from scripts import RunT3GeneralSequentialFitAblation as original
    from Numerical.orchestration.GeneralRunningBenchmark import initial_state
    from Numerical.running.UVRunner import run_uv_segment
    from Numerical.running.FermionThresholdBoundary import build_intermediate_scalar_boundary
    from Numerical.fitting.OscillationFit import OscillationFitTarget
    from Numerical.orchestration.RunSmallLambdaRefit import TARGET
    from Numerical.lo_certification.couplings.extra_hhss_geometry import model_extra_geometry

    if model in ('B_-1','C_-1'):
        raise ValueError('Use Batch 9 mixed runner for mixed models')
    mu_f,mu_s = original.SCALES[1:3]
    beta_uv,_beta_old,_,kernel,matchete_path,kind = g.nonmixed_inputs(
        model,canonical_run_id,generate_general_c5=False,
        wolframscript=wolframscript)
    overrides = original.load(saved/'general_uv_overrides.json')
    _,uv_initial=initial_state(beta_uv,mu_uv=original.SCALES[0],
                     mu_f=mu_f,mu_s=mu_s,overrides=overrides)
    _,channels = quartic_source
    target = OscillationFitTarget.from_json(TARGET)
    geometry = model_extra_geometry(model,registry['models'][model],include_components=True)
    cases,source_info,endpoint_comparisons = {},{},{}
    beta_files = {}
    frozen_quad_by_uv = {}
    for uv_on in (False,True):
        if uv_on:
            uv_sol=run_uv_segment(uv_initial,beta_uv,mu_f,rtol=1e-7,atol=1e-10)
            uv_end,uv_nfev = uv_sol.final_state,int(uv_sol.nfev)
        else:
            uv_end,uv_nfev=original.frozen_state(uv_initial,mu_f),0
        tree,tree_diag=original.nonmixed_boundary(model,uv_end,
                                            repo_root=original.ROOT,mu_f=mu_f)
        source_info[str(int(uv_on))]=tree_diag
        eft_start=build_intermediate_scalar_boundary(uv_end,matching_scale_gev=mu_f)
        start_data=_state_values(eft_start)
        start_data.update(extra_boundary(model,eft_start,registry).extra)
        q_basis=general.HHSSLinearBasis(model,tree,channels,registry,geometry,
                                       saved/'fitted_config.json',start_data)
        eft_ext=extra_boundary(model,eft_start,registry)
        beta_path=OUTPUT/'live_beta'/f'{model.replace("-","m").replace("+","p")}_eft1.json'
        # One fresh full quartic RGBeta export per distinct representation. The
        # U switch changes boundary couplings, not symbolic beta expressions.
        if not uv_on:
            full_beta.export_fresh_eft1_betas(model,eft_start.representation,
                                              beta_path,wolframscript=wolframscript)
        elif not beta_path.is_file():
            raise FileNotFoundError('Expected freshly generated U=0 EFT1 beta export')
        effective_path = beta_path.with_name(beta_path.stem + '_frozen_quadratic.json')
        beta, _, frozen_quadratics = checked_eft1_betas_with_frozen_quadratics(
            model, eft_ext, registry, beta_path, effective_path)
        beta_files[str(int(uv_on))]=str(beta_path)
        frozen_quad_by_uv[str(int(uv_on))] = list(frozen_quadratics)
        for eft_on in (False,True):
            source_c=np.asarray(tree.raw,complex)*complex(spinor_factor)
            combined,overlay=run_joint(model,eft_start,source_c,q_basis,
                                         mu_s,registry,beta,eft_on=eft_on)
            eft_end=combined['end_state']
            end_data=_state_values(eft_end)
            end_data.update(combined['extra_endpoint'])
            parity=q_basis.assert_reproduces(end_data,'muS')
            data_for_source=general._source_scalar_state(end_data)
            snap=original.snapshot(model,data_for_source,saved/'fitted_config.json')
            tensor=original.scalar_tensors(model,tree,channels,snap,registry,
                                           geometry=geometry)
            scalar=original.ScalarEndpoint(tensor,data_for_source)
            one_scale=kernel.evaluate(uv_end,eft_end)
            # The existing diagnostic frozen-seven-term solver needs the
            # extras at muS as well. Keep the replacement local to this call.
            old=general.general_quartic
            slot=overlay._hhss_slot
            def end_quartic(inner_model,data,hhss):
                slot['value']=hhss
                supplied=dict(data)
                supplied.update(combined['extra_endpoint'])
                return overlay.evaluate(supplied)
            general.general_quartic=end_quartic
            try:
                block=general._scalar_endpoint_prediction(
                    model,one_scale,tree,scalar,target,
                    uv_on=uv_on,eft_on=eft_on,
                    uv_nfev=uv_nfev,eft_nfev=combined['nfev'],
                    spinor_factor=spinor_factor,vev=vev,simultaneous=combined)
            finally:
                general.general_quartic=old
            for rec in block.values():
                rec.update(exploratory_full_non_HHSS_quartics_in_Wilson_RHS=True,
                           physical_rgbeta_mapping_certified=False,
                           quartic_mapping=LABEL,
                           extra_quartics_zero_at_muF=True,
                           frozen_quadratic_beta_names=list(frozen_quadratics),
                           frozen_quadratic_approximation_used=bool(frozen_quadratics),
                           extra_quartics_at_muS={k:{'real':v.real,'imag':v.imag}
                                                  for k,v in combined['extra_endpoint'].items()},
                           hhss_linear_basis_parity_rel=parity)
            cases.update(block)
            endpoint_comparisons[f'{int(uv_on)}{int(eft_on)}']={
                'hhss_basis_parity_rel':parity,
                'combined_nfev':combined['nfev'],
                'extra_eft1_running_enabled':bool(eft_on),
                'extra_quartic_slots':len(eft_ext.extra),
            }
    return cases,{'matching_kind':kind,'fermion_tree':source_info,
                  'full_general_C5_matchete_path':str(matchete_path),
                  'fitted_uv_sha256':original.hash_file(saved/'general_uv_overrides.json'),
                  'endpoint_comparisons':endpoint_comparisons,
                  'live_full_eft1_beta_paths':beta_files,
                  'exploratory_frozen_quadratic_beta_names_by_uv':frozen_quad_by_uv,
                  'simultaneous_rhs':True,
                  'exploratory_full_non_HHSS_quartics_in_Wilson_RHS':True,
                  'extra_quartics_zero_at_muF':True,
                  'quartic_mapping':LABEL,
                  'physical_rgbeta_mapping_certified':False}


def main(argv=None):
    require_sources()
    args=list(sys.argv[1:] if argv is None else argv)
    if '--resume' in args:
        raise ValueError('Use a fresh --run-id; no resume in exploratory full-quartic runner')
    from scripts import RunT3GeneralSequentialFitAblation as original
    old_nonmixed=original._nonmixed
    old_mixed=original._mixed
    old_out=original.OUT
    old_save=original.save
    def save_with_scope(path, data):
        if isinstance(data,dict):
            data=dict(data)
            if data.get('format') == 't3-full-general-16_U_E_W_S_ablation-v1':
                data.update(
                    format='t3-exploratory-nonmixed-joint-eft1-fullq-mufzero-v1',
                    one_scale_C5_used_as_W_on_reference=False,
                    running_background='extended_full_RGBeta_EFT1_plus_seven_term_Wilson',
                    additional_quartics_initialized_zero_at_muF=True,
                    additional_quartics_uv_evolution_performed=False,
                    missing_quadratic_betas_held_constant=True,
                    matchete_rgbeta_normalization_certified=False,
                    model_scope='14_nonmixed; B_-1_C_-1_separate',
                    quartic_mapping=LABEL,
                    complete_physical_C5=False)
            elif data.get('status')=='Computed' and isinstance(data.get('cases'),dict):
                data['additional_quartics_initialized_zero_at_muF']=True
                data['physical_quartic_normalization_certified']=False
        old_save(path,data)
    original.save=save_with_scope
    original._nonmixed=_nonmixed_full
    original._mixed=lambda model,*a,**k: (_ for _ in ()).throw(
        ValueError(f'{model}: use Batch 9 mixed-quartic CLI'))
    original.OUT=OUTPUT
    try:
        return original.main(args)
    finally:
        original._nonmixed=old_nonmixed
        original._mixed=old_mixed
        original.OUT=old_out
        original.save=old_save


if __name__=='__main__':
    raise SystemExit(main())
