"""T3 first-fermion-threshold spinor certification (fail-closed).

This package does not assert a physical GammaCC/Weyl normalization without
independent matching evidence, and never constructs an intermediate EFT with
only a subset of the heavy BSM fields removed.
"""
from .boundary import (MODELS, EXPECTED, BoundaryError, certify_source_release,
                       validate_model, weave_source_handoff)
from .spinors import convention_probe, exchange_probe

__all__ = ['MODELS', 'EXPECTED', 'BoundaryError', 'certify_source_release',
           'validate_model', 'weave_source_handoff', 'convention_probe', 'exchange_probe']
