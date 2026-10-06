"""Compatibility import for Numerical.running.T3Trajectory."""

from Numerical.running.T3Trajectory import *  # noqa: F401,F403
from Numerical.running.T3Trajectory import __dict__ as _impl_dict

globals().update(
    {
        name: value
        for name, value in _impl_dict.items()
        if not name.startswith("__")
    }
)

del _impl_dict
