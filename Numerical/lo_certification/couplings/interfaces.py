"""Optional read-only adapters to EXISTING seven-term and joint scalar solvers.

No alternative RGE, beta function or scalar bubble is implemented here.
Solver execution only consumes files independently supplied at the matching
conventions of the current repository; candidate previews are NOT substituted.
"""
from __future__ import annotations
from pathlib import Path
from .batch import write_json, read_json, digest


def check_existing_solver_inputs(model: str, boundary: Path,
                                 master_couplings: Path, joint_background: Path,
                                 output_path: Path) -> dict:
    """If both complete existing-format inputs exist, run the existing bridge.

    Both source-basis/muF and matching-basis/muS are validated by that solver.
    This diagnostic remains per unit uncertified N_Gamma, even on success.
    """
    from RGE.running.intermediate.IntegratedEFT1JointScalarBridge import evaluate
    if not master_couplings.is_file() or not joint_background.is_file():
        missing=[str(p) for p in (master_couplings, joint_background) if not p.is_file()]
        record={'model':model,'status':'blocked','reason':'missing external fully enumerated inputs',
                'missing':missing, 'physical_C5_certified':False}
        write_json(output_path,record)
        return record
    try:
        result=evaluate(boundary,master_couplings,joint_background)
        record={'model':model,'status':'conditional_existing_solver_passed',
                'master_input_sha256':digest(master_couplings),
                'joint_background_sha256':digest(joint_background),
                'existing_seven_term_joint_solver_result':result,
                'physical_C5_certified':False,
                'source_to_master_GammaCC_factor':None}
    except (ValueError,TypeError,KeyError,ImportError,RuntimeError) as exc:
        record={'model':model,'status':'blocked',
                'reason':'existing solver rejected physical interface: '+str(exc),
                'physical_C5_certified':False}
    write_json(output_path, record)
    return record
