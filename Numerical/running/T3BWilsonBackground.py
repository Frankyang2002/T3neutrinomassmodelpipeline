"""Export the exact one-generation T3-B background from an intermediate scalar run.

Use export_background(result, path) on an IntermediateScalarRunningResult.
For three-generation states, diagonal first-generation Yukawas are NOT a
valid replacement for the flavour tensor RGE; this adapter refuses them.
"""
from __future__ import annotations
import json
from pathlib import Path
import numpy as np


def _complex_pair(value):
    z = complex(value)
    return [float(z.real), float(z.imag)]


def _one_generation_yukawa(matrix, name):
    array = np.asarray(matrix, dtype=complex)
    if array.shape != (1, 1):
        raise ValueError(f"{name} has shape {array.shape}; one-generation diagnostic requires (1,1). "
                         "Do not replace 3x3 Yukawa matrices by their (0,0) entries.")
    return _complex_pair(array[0, 0])


def background_from_trajectory(result):
    """Convert the existing IntermediateScalarRunningResult to Wilson samples."""
    if not result.solver_success or result.n_points < 2:
        raise ValueError("A successful scalar trajectory with two or more points is required")
    samples = []
    for index in range(result.n_points):
        state = result.state_at_index(index)
        rep = state.representation
        if rep.shared_scalar or (rep.d_s1, rep.d_s2, rep.d_f, rep.alpha) != (2, 2, 1, -1):
            raise ValueError("Only ordinary T3-B alpha=-1 is supported")
        if state.lambda12Adj is None:
            raise ValueError("The T3-B state has no lambda12Adj coupling")
        samples.append({
            "mu_gev": float(result.mu_gev[index]),
            "g2": float(state.sm.g2),
            "lambdaH": float(state.sm.lambdaH),
            "lambda12": float(state.lambda12),
            "lambda12Adj": float(state.lambda12Adj),
            "lambdaT3": _complex_pair(state.lambdaT3),
            "yu": _one_generation_yukawa(state.sm.yu, "yu"),
            "yd": _one_generation_yukawa(state.sm.yd, "yd"),
            "ye": _one_generation_yukawa(state.sm.ye, "ye"),
        })
    return {"status":"Success","model":"T3-B alpha=-1",
            "scope":"one-generation diagnostic background", "samples":samples}


def export_background(result, output_path):
    output = background_from_trajectory(result)
    path = Path(output_path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(output, indent=2), encoding="utf-8")
    return output
