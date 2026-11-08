"""Check crossed-Yukawa Ward-covariance candidate, no physical certification.

Run from project root:
  python scripts/CheckBatch34CrossedCandidate.py --output output/batch34_candidate/ward_candidate.json
The output is intentionally never labelled a verified A(-4) calculation.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys
import numpy as np

ROOT=Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:sys.path.insert(0,str(ROOT))

from scripts.CheckBatch34GaugeCovariance import example, ward_c, ward_y, ward_l, diagnose
from RGE.running.intermediate.GeneralEFT1SevenTermFlow import (
    Couplings, TERMS, master_beta, total_beta,
)
from RGE.running.intermediate.GeneralEFT1CrossedWardCandidate import (
    crossed_candidate, candidate_beta, CANDIDATE_STATUS,
)
TOL=2e-11


def u1_nontrivial_fixture():
    """C_{22} ~ XX and y_{02} ~ X, with q=(0,0,-1), q(X)=+1.

    A crossed C* contribution induces C_{00} ~ X*X. Both are gauge invariant.
    The channel is nonzero, making this more than a null Ward check.
    """
    v=np.array([1.0,1.0j])/np.sqrt(2)
    t=np.diag([0,0,-1]).astype(complex)
    theta=np.array([[0,1j],[-1j,0]],complex)
    y=np.zeros((3,3,2),complex)
    y[0,2]=y[2,0]=.35*v
    C=np.zeros((3,3,2,2),complex)
    C[2,2]=(.2+.1j)*np.outer(v,v)
    inp=Couplings(y,np.zeros((2,)*4),((.30,((t,theta),)),))
    return C,inp


def _ward_abs(coeff, inputs):
    return max(float(np.max(abs(ward_c(coeff,t,s))))
               for _,pairs in inputs.gauge for t,s in pairs)


def diagnose_candidate():
    input_reports=[]
    for name, builder in [
            ("SM_L_eC_H_with_Weinberg",example),
            ("nonzero_U1_XX_to_XstarX",u1_nontrivial_fixture)]:
        c, x=builder()
        y_ward=max(float(np.max(abs(ward_y(x.yukawa,t,s)))) for _,ps in x.gauge for t,s in ps)
        c_ward=_ward_abs(c,x)
        lam_ward=max(float(np.max(abs(ward_l(x.quartic,s)))) for _,ps in x.gauge for t,s in ps)
        if max(y_ward,c_ward,lam_ward)>TOL:
            raise RuntimeError("Diagnostic input Ward identity failed")
        terms=candidate_beta(c,x)
        beta=total_beta(terms)
        crossed=terms["crossed_yukawa"]
        original=master_beta(c,x,corrected_crossed=True)["crossed_yukawa"]
        cross_norm=float(np.linalg.norm(crossed))
        full_norm=float(np.max(abs(beta)))
        ward_beta=_ward_abs(beta,x)
        input_reports.append({
            "name":name,"y_ward_abs":y_ward,"C_ward_abs":c_ward,
            "lambda_ward_abs":lam_ward,
            "crossed_candidate_norm":cross_norm,
            "crossed_original_norm":float(np.linalg.norm(original)),
            "crossed_candidate_ward_abs":_ward_abs(crossed,x),
            "beta_ward_abs":ward_beta,
            "beta_ward_relative":ward_beta/max(full_norm,1e-300),
            "term_ward_abs":{key:_ward_abs(terms[key],x) for key in TERMS},
        })
    if input_reports[1]["crossed_candidate_norm"]<1e-8:
        raise RuntimeError("Candidate nontrivial test produced a vanishing crossed term")
    passed=all(row["beta_ward_relative"]<TOL for row in input_reports)
    old=diagnose()
    return {
        "format":"t3-batch34-crossed-candidate-diagnostic-v1",
        "source_head":"12128de602ac6bf109dbad50dbac8446ed629434",
        "status":"WARD_DIAGNOSTICS_PASSED_UNVERIFIED_CANDIDATE" if passed else "BLOCKED_WARD_DIAGNOSTICS",
        "candidate_equation":"crossed = sum_{sigma({i,j} x {a,b})} y_jla y_ikc conjugate(C_lkbc)",
        "original_v1_contraction":"crossed = sum_{sigma({i,j} x {a,b})} y_jla y_ikc C_lkbc",
        "anti_linear_blocks":"both crossed_yukawa and conjugate_coefficient_yukawa under candidate",
        "reference_conflict":"Uploaded arXiv:2512.16890v1 Eq.(4.85) explicitly shows unstarred C in the crossed term; this candidate is NOT an asserted correction to the paper",
        "original_ward_status":old["status"],
        "candidate_status":CANDIDATE_STATUS,
        "tests":input_reports,
        "physical_A_minus4_rge_certified":False,
        "chatA_canonical_boundary_verified":False,
        "chatB_running_general_EFT1_couplings_verified":False,
        "notes":"Necessary gauge-covariance checks pass, but an independent RGE derivation is still required. No physical ODE trajectory is modified here.",
    }


def main(argv=None):
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output",type=Path)
    a=parser.parse_args(argv)
    r=diagnose_candidate()
    if a.output:
        a.output.parent.mkdir(parents=True,exist_ok=True)
        with a.output.open("x",encoding="utf8") as f:
            json.dump(r,f,indent=2,allow_nan=False)
            f.write("\n")
    print(r["status"])
    for row in r["tests"]:
        print(f'{row["name"]}: beta Ward={row["beta_ward_relative"]:.3g}, crossed norm={row["crossed_candidate_norm"]:.3g}')
    return 0 if r["status"].startswith("WARD_DIAGNOSTICS_PASSED") else 2

if __name__=="__main__":raise SystemExit(main())
