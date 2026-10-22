"""Evaluate optimised and shared A(-2) 5e-10 benchmarks using saved artefacts.

Run from project root after fitting all 16 representations:
    python -m Numerical.orchestration.RunLambda5e10Comparisons --mode optimised --resume
    python -m Numerical.orchestration.RunLambda5e10Comparisons --mode shared --resume

Reuses the verified 5e-9 comparison numerical evaluator, ODE continuation,
NuFIT cross-check and per-model figure generation. No symbolic regeneration.
"""
from __future__ import annotations

from Numerical.orchestration import RunLambda5e9Refit as study
from Numerical.orchestration import RunLambda5e9Comparisons as runner

study.OUT = study.base.ROOT / 'output' / 'studies' / 'lambda_5e-10_refit'
study.LAMBDA = {'normal': 5e-10, 'E': 2.5e-10}


def main() -> None:
    runner.main()


if __name__ == '__main__':
    main()
