"""Compatibility import for the canonical Numerical.core module."""

from Numerical.core.RGBetaEvaluator import *  # noqa: F401,F403
from Numerical.core.RGBetaEvaluator import __dict__ as _impl_dict

globals().update(
    {
        name: value
        for name, value in _impl_dict.items()
        if not name.startswith("__")
    }
)

del _impl_dict
