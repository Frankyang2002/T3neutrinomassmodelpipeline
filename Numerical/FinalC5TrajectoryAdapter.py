"""Compatibility import for Numerical.running.FinalC5TrajectoryAdapter."""

from Numerical.running.FinalC5TrajectoryAdapter import *  # noqa: F401,F403
from Numerical.running.FinalC5TrajectoryAdapter import __dict__ as _impl_dict

globals().update(
    {
        name: value
        for name, value in _impl_dict.items()
        if not name.startswith("__")
    }
)

del _impl_dict
