"""Selected extra-quartic LLHH source -> projected Weinberg -> SM running.

Not a full extended two-loop neutrino-mass calculation. Uses the same
renormalisable B/C(-1) backgrounds as Stage16/18 and the existing Weinberg ODE.
"""
from __future__ import annotations
from pathlib import Path
from typing import Any, Mapping

import numpy as np

from Numerical.experimental.GeneralT3WeinbergProjection import (
    encode_matrix, project_stage18_integral,
    verify_full_seed_weinberg_projection,
    verify_production_c5_normalization,
)


def transport_stage19_model(model:str, config:Mapping[str,Any],
                            uv_payload:Mapping[str,Any],eft_payload:Mapping[str,Any],
                            stage17:Mapping[str,Any],stage18:Mapping[str,Any],
                            seed:Mapping[str,Any],wilson_rge:Mapping[str,Any],
                            final_weinberg:Mapping[str,Any]) -> dict[str,Any]:
    """Certify SU2/physical normalization, then evolve only selected delta K.

    No full-Wilson-system claims are made; this propagates a gauge-invariant
    source which is already LLHH through the *tree-level* scalar boundary.
    """
    if model not in ('B_m1','C_m1'):
        raise ValueError('Only B/C(-1) supported')
    tensor_certificate=verify_full_seed_weinberg_projection(model,seed,stage17)
    normalization=verify_production_c5_normalization(wilson_rge,final_weinberg)
    k_s, projection=project_stage18_integral(stage18,model,tensor_certificate,normalization)

    from Numerical.experimental.GeneralT3CoupledUV import construct_uv_state
    from Numerical.experimental.GeneralT3HierarchicalTrajectory import (
        _threshold_scales,_guard_running_hierarchy,_validate_payload,
    )
    from Numerical.running.T3Trajectory import (
        run_renormalisable_t3_trajectory,continue_with_weinberg_running,
    )
    from physics.NeutrinoMass import numerical_neutrino_mass_matrix

    state=construct_uv_state(model,config)
    scales=_threshold_scales(config,state)
    for key,expected in (('mu_F',scales.mu_f),('mu_S',scales.mu_s)):
        got=float(stage18['scales_gev'][key])
        if not np.isfinite(got) or not np.isclose(got,expected,rtol=1e-10,atol=1e-10):
            raise ValueError(f'Stage18 {key} does not match current benchmark')
    _validate_payload(uv_payload,'UV',model)
    _validate_payload(eft_payload,'EFT1',model)
    trajectory=run_renormalisable_t3_trajectory(
        state,uv_payload,eft_payload,
        mu_fermion_threshold_gev=scales.mu_f,
        mu_scalar_threshold_gev=scales.mu_s,
        uv_save_scales_gev=np.geomspace(state.mu_gev,scales.mu_f,9),
        eft1_save_scales_gev=np.geomspace(scales.mu_f,scales.mu_s,9),
        rtol=3e-8,atol=3e-12,
    )
    masses=_guard_running_hierarchy(trajectory,scales)
    # Re-run consistency: it must be the same Stage18 scalar background.
    recorded=stage18.get('eft1_source_quartics',[])
    if len(recorded)<5:
        raise ValueError('Stage18 scalar trajectory audit is missing')
    for ix,state_at in ((0,trajectory.eft1_initial_state),(-1,trajectory.eft1_threshold_state)):
        row=recorded[ix]
        if not np.isclose(float(row['mu_gev']),state_at.mu_gev,rtol=1e-10,atol=1e-8):
            raise ValueError('Stage18 EFT1 endpoint scale differs from rerun')
        for name in ('lambdaHHdagS2S2','lambdaHHdagS1barS1bar'):
            ref=complex(*row[name])
            current=complex(getattr(state_at,name))
            if not np.isclose(ref,current,rtol=1e-8,atol=1e-22):
                raise ValueError(f'Stage18 EFT1 {name} differs from rerun')
    final=continue_with_weinberg_running(trajectory,k_s,scales.mu_low)
    k_low=np.asarray(final.weinberg.K,dtype=complex)
    mass_s=numerical_neutrino_mass_matrix(k_s)*1e9
    mass_low=numerical_neutrino_mass_matrix(k_low)*1e9
    return {
        'status':'SelectedExtraWeinbergSourceTransported',
        'model':model,'projection_certificate':tensor_certificate,
        'normalization_certificate':normalization,
        'integral_projection':projection,
        'scales_gev':{'mu_uv':float(state.mu_gev),'mu_F':scales.mu_f,
                      'mu_S':scales.mu_s,'mu_low':scales.mu_low},
        'running_masses_gev':masses,
        'uv_nfev':int(trajectory.uv.nfev),'eft1_nfev':int(trajectory.eft1.nfev),
        'sm_weinberg_nfev':int(final.weinberg.nfev),
        'selected_extra_kappa_at_mu_s':encode_matrix(k_s,unit='GeV^-1'),
        'selected_extra_kappa_at_mu_low':encode_matrix(k_low,unit='GeV^-1'),
        'selected_extra_mass_at_mu_s_ev':encode_matrix(mass_s,unit='eV'),
        'selected_extra_mass_at_mu_low_ev':encode_matrix(mass_low,unit='eV'),
        'tree_level_scalar_projection_used':True,
        'finite_scalar_matching_implemented':False,
        'full_coupled_eft1_wilson_evolution_implemented':False,
        'complete_two_loop_result':False,
        'physical_full_neutrino_mass_prediction':False,
        'physics_scope':(
            'Source-selected O(hbar^2) if extra UV quartics begin at zero; '
            'template-calibrated tree-level LLHH scalar projection and final SM+Weinberg running. '
            'Not a complete O(hbar^2) matching, no finite scalar threshold correction '
            'or full EFT1 Wilson self-running.'
        ),
    }
