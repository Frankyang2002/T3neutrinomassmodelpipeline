"""Independent SYNTHETIC A(-4)-sized (6 Weyl, 12 real scalar) contracts.

Not gauge-invariant, not a real general-T3 coupling set, not Matchete-normalized.
Only the field *dimensions* and the canonical test parser are reused.
"""
from __future__ import annotations

from copy import deepcopy
import numpy as np

CERTIFICATES = (
    "matchete_gamma_cc_to_weyl_factor_and_sign",
    "canonical_C_ijab_operator_normalization",
    "complete_tree_dim5_operator_inventory",
    "other_dim5_sources_audited_at_one_loop",
    "complete_EFT1_yukawa_and_gauge_dictionary",
    "complete_scalar_quartic_counterterms",
    "all_heavy_scalar_components_and_mixed_masses",
    "weinberg_C5_symmetric_convention",
    "MSbar_subtraction_and_matching_log_partition",
)


def z(x):
    return {"real":float(complex(x).real),"imag":float(complex(x).imag)}


def matrix(x):
    x=np.asarray(x,complex)
    return {"real":x.real.tolist(),"imag":x.imag.tolist()}


def fixtures():
    # Source C: symmetric LLSS in the full real-component basis, with CP phases.
    op = []
    for p,q,v in ((1,1,2.0e-11+1.0e-11j),(1,3,1.1e-11-0.5e-11j),
                  (3,3,1.5e-11+0.7e-11j),(3,5,0.8e-11+0.3e-11j)):
        op.append({"indices":[p,q,5,7],"value":z(v)})
    op2=[{"indices":[p,q,6,8],"value":z(complex(0.27,-0.19)*complex(v["real"],v["imag"]))}
         for p,q,v in ((a["indices"][0], a["indices"][1], a["value"]) for a in op)]
    boundary = {
        "format":"t3-eft1-certified-master-weyl-v1", "model":"A_-4",
        "authority":"synthetic", "synthetic_notice":"UNPHYSICAL closure test only",
        "scheme":"MSbar", "mu_F_gev":100000.0, "mu_S_gev":2000.0,
        "n_fermions":6, "n_real_scalars":12,
        "lepton_components":[[1,2],[3,4],[5,6]],
        "higgs_real_components":[1,2,3,4],
        "scalar_field_blocks":[{"name":"H","role":"Higgs","indices":[1,2,3,4]},
          {"name":"S1","role":"heavy","indices":[5,6]},
          {"name":"S2","role":"heavy","indices":list(range(7,13))}],
        "master_lagrangian_convention":"C_ijab psi_i psi_j phi_a phi_b + hc",
        "weinberg_lagrangian_convention":"(1/2) C5_ij (L_i.H)(L_j.H) + hc",
        "certificates":{k:False for k in CERTIFICATES},
        "yukawa_components":[],"quartic_components":[],"gauge_sectors":[],
        "operators":[{"source_id":"SYN_LLSS_1","components":op},
                     {"source_id":"SYN_LLSS_2","components":op2}],
    }
    # NOTE: these arbitrary synthetic Yukawas are not gauge-invariant.
    y = [([1,1,5],0.10+0.08j),([1,3,7],-0.07+0.03j),
         ([3,3,6],0.06-0.01j),([5,1,8],0.04+0.02j),
         ([1,2,3],0.01+0.01j),([3,4,4],0.02-0.03j)]
    lam = [([3,3,5,7],0.09),([4,4,5,7],-0.09),
           ([1,1,5,7],0.09),([2,2,5,7],-0.09),
           ([3,3,6,8],-0.04),([4,4,6,8],0.04),
           ([5,5,7,7],0.13),([5,6,7,8],0.06)]
    theta = np.zeros((12,12),complex)
    # U(1)-like block-pair real-scalar Hermitian generators.
    for k in range(0,12,2):
        theta[k,k+1]=0.5j
        theta[k+1,k]=-0.5j
    tf = np.diag([-0.5,-0.5,-0.5,-0.5,-0.5,-0.5]).astype(complex)
    generator = {"name":"SYN_U1", "generators":[{"fermion":matrix(tf),"scalar":matrix(theta)}]}
    rows = []
    for mu,factor in ((2000.0,0.83),(12000.0,0.92),(100000.0,1.0)):
        rows.append({"mu_gev":mu,"yukawa_components":[
            {"indices":k,"value":z(v*factor)} for k,v in y],
            "quartic_components":[{"indices":k,"value":v*(0.75+0.25*factor)} for k,v in lam],
            "gauge_couplings":{"SYN_U1":0.22+0.08*factor}})
    couplings={"format":"t3-eft1-seven-term-coupling-samples-v1",
               "model":"A_-4", "scheme":"MSbar",
               "basis":"canonical_left_Weyl_real_scalar",
               "n_fermions":6,"n_real_scalars":12,
               "authority":"synthetic", "couplings_scale_dependent":True,
               "complete_general_EFT1_RGE_verified":False,
               "gauge_generators":[generator], "samples":rows,
               "fermion_field_blocks":[
                   {"name":f"L{q}","indices":[2*q-1,2*q]} for q in (1,2,3)],
               "scalar_field_blocks":deepcopy(boundary["scalar_field_blocks"]),
               "notice":"UNPHYSICAL: deliberately arbitrary Yukawas and quartics"}
    return boundary,couplings
