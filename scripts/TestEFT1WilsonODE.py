"""Quick synthetic numerical regression for the master-basis Wilson solver.

No T3 physical claim: the tiny test has two dynamical spectator scalar blocks,
vanishing gauge/Yukawa couplings, and explicitly chosen scalar quartics.
Run from repository root:
    python -X utf8 -m scripts.TestEFT1WilsonODE
"""
from __future__ import annotations

from copy import deepcopy
import numpy as np

from Numerical.running.GeneralEFT1WilsonODE import solve_coupled_wilson


def boundary() -> dict:
    return {
        "format": "t3-eft1-certified-master-weyl-v1",
        "model": "D_-2", "authority": "synthetic", "certificates": {},
        "scheme": "MSbar",
        "master_lagrangian_convention": "C_ijab psi_i psi_j phi_a phi_b + hc",
        "weinberg_lagrangian_convention": "(1/2) C5_ij (L_i.H)(L_j.H) + hc",
        "n_fermions": 6, "n_real_scalars": 6,
        "lepton_components": [[1,2],[3,4],[5,6]],
        "higgs_real_components": [1,2,3,4],
        "scalar_field_blocks": [
            {"name":"H","role":"Higgs","indices":[1,2,3,4]},
            {"name":"S1","role":"heavy","indices":[5]},
            {"name":"S2","role":"heavy","indices":[6]},
        ],
        "mu_F_gev": 1e5, "mu_S_gev": 1e3,
        "yukawa_components":[], "gauge_sectors":[],
        "quartic_components": [
            {"indices":[3,3,5,5],"value":0.10},
            {"indices":[5,5,5,5],"value":0.20},
        ],
        "operators": [{"source_id":"LLS1S1", "components":[
            {"indices":[1,1,5,5],"value":1e-5}]}],
    }


def main() -> int:
    b = boundary()
    a = deepcopy(b); a["mu_gev"] = 1e5
    z = deepcopy(b); z["mu_gev"] = 1e3
    res = solve_coupled_wilson(b, [a,z], allow_synthetic=True,
                               rtol=1e-9, atol=1e-19, max_step_log=0.1)
    assert res["solver"]["success"]
    assert res["operator_shape"] == [6,6,6,6]
    assert res["wilson_pair_symmetry_residual"] < 1e-17
    assert res["flavour_symmetry_residual"] < 1e-17
    c5 = res["C5_at_muS_component_gev_inverse"]
    num = np.array(c5["real"]) + 1j*np.array(c5["imag"])
    ll = res["C5_fixed_couplings_leading_log_gev_inverse"]
    lead = np.array(ll["real"]) + 1j*np.array(ll["imag"])
    assert np.isfinite(num).all() and np.isfinite(lead).all()
    assert abs(num[0,0]) > 0 and np.linalg.norm(num[1:,:]) < 1e-17
    # The integrated equation reduces to the fixed-order LL in small-log
    # expansion; self running creates higher-order logarithmic corrections.
    assert np.linalg.norm(num-lead)/np.linalg.norm(lead) < 0.05
    print("Synthetic master-basis solver PASS")
    print(f"C5[0,0] numerical = {num[0,0]:.12e} GeV^-1")
    print(f"C5[0,0] leading log = {lead[0,0]:.12e} GeV^-1")
    print(f"nfev = {res['solver']['nfev']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
