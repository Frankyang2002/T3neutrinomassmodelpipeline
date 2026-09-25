"""Compatibility surface for the historical full-flavor Weinberg stage path.

The canonical RGE-stage implementation now lives in
``RGE.running.weinberg.FullFlavorWeinbergStage``. Historical function and
neutrino-mass imports remain available here for existing callers.
"""

from physics.NeutrinoMass import (
    build_neutrino_mass_matrix,
    run_symbolic_neutrino_mass_stage,
)
from RGE.running.weinberg.FullFlavorWeinbergStage import (
    run_full_flavor_weinberg_stage,
)

# Historical public function name.
run_flavor_matched_weinberg_rge = run_full_flavor_weinberg_stage

__all__ = [
    "build_neutrino_mass_matrix",
    "run_flavor_matched_weinberg_rge",
    "run_full_flavor_weinberg_stage",
    "run_symbolic_neutrino_mass_stage",
]
