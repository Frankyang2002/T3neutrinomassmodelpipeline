"""Fast hermetic tests for candidate tensor insertion and Wilson pair mixing.

Run: python -X utf8 scripts/TestT3CoupledScalarCached.py
No fitted model files or Wolfram installation required.
"""
from __future__ import annotations

import numpy as np
from scipy.integrate import solve_ivp

from Numerical.running.T3CoupledQuarticContext import (
    CandidateQuarticContext, _put_hhss, symmetry_residual,
)


def run():
    h = np.zeros((4, 4, 8, 8))
    h[2, 2, 0, 0] = 0.40
    full = np.zeros((12,) * 4)
    full[4, 4, 4, 4] = 0.7
    full[5, 5, 5, 5] = -0.3
    column = np.zeros_like(full)
    column[4, 4, 4, 4] = 0.6
    c = CandidateQuarticContext(
        'A_-4', {'g': 0.15}, (('g','real'),), full,
        np.stack([column]), {}, 'synthetic',
    )
    s0 = {'g':0.15}
    s1 = {'g':0.83}
    q0 = c.tensor(s0, h)
    q1 = c.tensor(s1, h)
    assert abs(q1[4,4,4,4] - q0[4,4,4,4] - .6 * (.83-.15)) < 1e-12
    assert q1[2,2,4,4] == h[2,2,0,0]
    assert q1[4,4,2,2] == h[2,2,0,0]
    assert abs(q1[5,5,5,5] + .3) < 1e-12
    assert symmetry_residual(q1) < 1e-15
    assert c.number_of_calls == 2

    # Real symmetric quartic mixes scalar Wilson sectors via
    # beta_ijab[scalar_pair] = lambda_abcd C_ijcd.  A full tensor changes
    # C_SS self-running even when the LLSS->LLHH portal is unchanged.
    # Reduce the source to (C_HH, C_SS) = (H, S), and evolve downward.
    t0, t1 = np.log(1e5), np.log(1e3)
    loop = 1 / (16 * np.pi**2)
    def solve(self_quartic):
        def rhs(_, vec):
            hcoef, scoef = vec
            return [loop*.4*scoef, loop*(.4*hcoef+self_quartic*scoef)]
        sol=solve_ivp(rhs,(t0,t1),[0.,1.],method='DOP853',rtol=1e-12,atol=1e-14)
        assert sol.success
        return sol.y[:,-1]
    old = solve(q0[4,4,4,4]); new=solve(q1[4,4,4,4])
    assert abs(new[1]-old[1]) > 1e-6
    assert abs(new[0]-old[0]) > 1e-8
    assert np.isfinite(new).all()
    print('PASS: cached 12^4 coupling dependence and exact HHSS override')
    print('PASS: non-HHSS scalar self-coupling changes coupled Wilson IVP at second order')


if __name__ == '__main__':
    run()
