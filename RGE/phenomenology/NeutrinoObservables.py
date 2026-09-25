"""Compatibility imports for the historical neutrino-observables path.

Physical neutrino observables now live in ``physics.NeutrinoObservables``.
This module preserves existing imports used by numerical scans and trajectory
code while containing no phenomenology implementation of its own.
"""

from physics.NeutrinoObservables import (
    GEV_TO_EV,
    NeutrinoObservables,
    calculate_neutrino_observables,
    run_neutrino_observables_stage,
    takagi_factorization,
)

__all__ = [
    "GEV_TO_EV",
    "NeutrinoObservables",
    "calculate_neutrino_observables",
    "run_neutrino_observables_stage",
    "takagi_factorization",
]
