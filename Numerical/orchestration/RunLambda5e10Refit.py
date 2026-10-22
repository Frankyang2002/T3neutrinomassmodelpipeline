"""5e-10 Sobol/refinement study using the verified 5e-9 optimiser.

Run from the repository root:
    python -m Numerical.orchestration.RunLambda5e10Refit prepare
    python -m Numerical.orchestration.RunLambda5e10Refit fit --resume

Starts from the individually fitted 5e-9 Yukawas. Writes only to
output/studies/lambda_5e-10_refit. Existing runs remain unchanged.
"""
from __future__ import annotations

from Numerical.orchestration import RunLambda5e9Refit as runner

# The original module uses these globals in prepare/fit/cross/analyse and its
# argparse entry point. Patch only this Python process, not source files.
PREVIOUS = runner.OUT
runner.OLD = PREVIOUS
runner.OUT = runner.base.ROOT / 'output' / 'studies' / 'lambda_5e-10_refit'
runner.LAMBDA = {'normal': 5e-10, 'E': 2.5e-10}


def main() -> None:
    runner.main()


if __name__ == '__main__':
    main()
