"""Compatibility import for Numerical.diagnostics.IntermediateWeinbergDiagnostics."""

from Numerical.diagnostics.IntermediateWeinbergDiagnostics import *  # noqa: F401,F403
from Numerical.diagnostics.IntermediateWeinbergDiagnostics import __dict__ as _impl_dict

globals().update(
    {
        name: value
        for name, value in _impl_dict.items()
        if not name.startswith("__")
    }
)

del _impl_dict
