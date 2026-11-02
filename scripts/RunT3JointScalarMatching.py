"""CLI: single-threshold scalar bubble C5 for 16 T3 models, with strict provenance."""
from __future__ import annotations
import argparse
import json
from pathlib import Path
import sys

import numpy as np

ROOT=Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path: sys.path.insert(0,str(ROOT))

from Numerical.joint_scalar_matching.adapter import (ALL_MODELS,load_input,save_input,
                                                        normalize_model,scalar_groups)
from Numerical.joint_scalar_matching.fixtures import make_fixture
from Numerical.joint_scalar_matching.core import match_c5, project_weinberg, scalar_bubble_matching

DEFAULT_OUTPUT=Path("output/general_observable/lo_joint_scalar_matching")


def matrix_json(a: np.ndarray) -> dict:
    arr=np.asarray(a,complex)
    return {"real":arr.real.tolist(),"imag":arr.imag.tolist()}


def compute(inp, *, strict_projection=True) -> dict:
    r=match_c5(inp.c_llss,inp.lambda_hhss,inp.mass_squared,inp.mu_s_gev,
               c5_eft1=inp.c5_eft1,reference_scale=inp.reference_scale_gev,
               group_names=scalar_groups(inp.model),strict_projection=strict_projection)
    contributions={"C5_inherited_EFT1_at_muS":matrix_json(r["c5_eft1"]),
                   "scalar_bubble_all_eigenstate_pairs":matrix_json(r["c5_scalar_bubble"]),
                   "scalar_bubble_finite_at_reference":matrix_json(r["c5_reference_finite"]),
                   "scalar_bubble_explicit_matching_scale_log":matrix_json(r["c5_explicit_mu_log"])}
    pair_output=[]
    for row in r["pair_data"]:
        pair_residual=project_weinberg(r["pair_amplitudes"][row["pair"]],strict=False)[1]
        pair_output.append({**row,"delta_C5":matrix_json(r["pair_c5"][row["pair"]]),
                            "individual_pair_Weinberg_projection_residual":pair_residual})
    channel_outputs=[]
    channels=(inp.quartic_channels if inp.quartic_channels is not None
              else inp.lambda_hhss[np.newaxis,...])
    channel_names=(inp.quartic_channel_names if inp.quartic_channels is not None
                   else ("aggregated_quartics_independence_not_certified",))
    channel_total=np.zeros((3,3),complex)
    for name,qchannel in zip(channel_names,channels):
        chunk=scalar_bubble_matching(inp.c_llss,qchannel,inp.mass_squared,inp.mu_s_gev,
                                     reference_scale=inp.reference_scale_gev)
        c5, residual=project_weinberg(chunk["amplitude"],strict=False)
        channel_total+=c5
        channel_outputs.append({"name":name,"C5_gev_inverse":matrix_json(c5),
                                "individual_channel_projection_residual":residual})
    if not np.allclose(channel_total,r["c5_scalar_bubble"],rtol=1e-10,atol=1e-28):
        raise RuntimeError("Quartic-channel decomposition does not reproduce full C5")
    return {"schema_version":1,"model":inp.model,"matching":"both scalars together at mu_S",
            "sequence":"UV -> F all at mu_F -> EFT1(S1,S2) -> S1+S2 all at mu_S -> SMEFT",
            "mu_f_gev":inp.mu_f_gev,"mu_s_gev":inp.mu_s_gev,
            "reference_scale_gev":inp.reference_scale_gev,
            "scheme":"MSbar","tensors":"T3_REAL_LLSS_HHSS_WEINBERG_v1",
            "truncation":"one scalar-pair bubble with one dimension-5 LLSS insertion and one scalar quartic; zeroth momentum, leading 1/M_F",
            "scalar_mass_origin":inp.scalar_mass_origin,
            "synthetic_fixture":inp.is_fixture,"provenance":inp.provenance,
            "physical_full_one_loop_certified":False,
            "note":"Finite bubble normalization follows source EFT1 quartic RGE; independent full native amplitude and counterterm certification not supplied.",
            "weinberg_projection_residual":r["projection_residual"],
            "weinberg_projection_passed":r["projection_passed"],
            "C5_total_gev_inverse":matrix_json(r["c5_total"]),
            "contributions_gev_inverse":contributions,
            "scalar_eigenstate_masses_gev":np.sqrt(r["eigenvalues_sq_gev2"]).tolist(),
            "real_orthogonal_rotation_old_equals_U_eigen":r["mass_rotation_U"].tolist(),
            "individual_eigenstate_pair_contributions":pair_output,
            "independent_quartic_channel_contributions":channel_outputs,
            "quartic_channels_explicitly_supplied":inp.quartic_channels is not None,
            "RG_accounting":"C5_EFT1(mu_S) is input once; no ln(mu_F/mu_S) is added by this module. Reference/explicit-mu pieces sum to the one threshold bubble and are NOT extra running corrections.",
            "mu_log_derivative_identity":r["rge_log_identity"]}


def main(argv=None) -> int:
    p=argparse.ArgumentParser(description=__doc__)
    mode=p.add_mutually_exclusive_group(required=True)
    mode.add_argument("--fixtures",action="store_true",help="Generate + match 16 explicitly synthetic fixtures")
    mode.add_argument("--input-dir",type=Path,help="Folder containing 16 model metadata JSON and NPZ files")
    p.add_argument("--model",help="One model (default all 16)")
    p.add_argument("--output-dir",type=Path,default=None)
    p.add_argument("--allow-nonweinberg-projection",action="store_true",
                   help="Diagnostic only: save projected C5 even if residual is nonzero")
    args=p.parse_args(argv)
    models=[normalize_model(args.model)] if args.model else list(ALL_MODELS)
    output=(args.output_dir if args.output_dir is not None else
            (DEFAULT_OUTPUT/"fixtures_demo" if args.fixtures else DEFAULT_OUTPUT)).resolve()
    # Only the agreed output subtree can be modified, including when using --fixtures.
    safe=(ROOT/DEFAULT_OUTPUT).resolve()
    if output != safe and safe not in output.parents:
        p.error(f"output directory must be within {safe}")
    output.mkdir(parents=True,exist_ok=True)
    summary=[]
    for model in models:
        if args.fixtures:
            inp=make_fixture(model)
            meta_path=save_input(output/"fixtures"/(model.replace("+","p").replace("-","m")+".json"),inp)
        else:
            assert args.input_dir is not None
            inputs=[path for path in args.input_dir.glob("*.json")]
            hits=[]
            for path in inputs:
                try:
                    content=json.loads(path.read_text(encoding="utf-8"))
                    if normalize_model(content.get("model",""))==model:
                        hits.append(path)
                except (ValueError,KeyError,TypeError):
                    continue
            if len(hits)!=1:
                raise RuntimeError(f"Expected exactly one metadata JSON for {model}; found {len(hits)}")
            meta_path=hits[0]
            inp=load_input(meta_path)
        result=compute(inp,strict_projection=not args.allow_nonweinberg_projection)
        name=model.replace("+","p").replace("-","m")
        target=output/(name+"_c5.json")
        target.write_text(json.dumps(result,indent=2,allow_nan=False)+"\n",encoding="utf-8")
        summary.append({"model":model,"output":str(target.relative_to(ROOT)) if ROOT in target.parents else str(target),
                        "synthetic_fixture":inp.is_fixture,
                        "weinberg_projection_residual":result["weinberg_projection_residual"],
                        "physical_full_one_loop_certified":False})
        print(f"{model}: {target} (projection residual {result['weinberg_projection_residual']:.2e})")
    (output/"summary.json").write_text(json.dumps({"model_count":len(models),
        "models":summary,"scope":"joint scalar pair, no scalar EFT2",
        "certified_full_physical_matching":False},indent=2)+"\n",encoding="utf-8")
    return 0


if __name__=="__main__":
    raise SystemExit(main())
