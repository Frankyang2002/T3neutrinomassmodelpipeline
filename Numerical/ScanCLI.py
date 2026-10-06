"""Compatibility import for Numerical.fitting.ScanCLI."""

from Numerical.fitting.ScanCLI import *  # noqa: F401,F403
from Numerical.fitting.ScanCLI import __dict__ as _impl_dict

globals().update(
    {
        name: value
        for name, value in _impl_dict.items()
        if not name.startswith("__")
    }
)

del _impl_dict
