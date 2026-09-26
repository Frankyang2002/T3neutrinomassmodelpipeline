from __future__ import annotations

from pathlib import Path


PROJECT_ROOT = Path(__file__).resolve().parents[1]
FIGURES = PROJECT_ROOT / "Numerical" / "RunningResultFigures.py"


def test_default_generation_uses_focused_result_set() -> None:
    source = FIGURES.read_text(encoding="utf-8")
    generate_source = source.split(
        "def generate_running_result_figures",
        1,
    )[1]

    generated = {
        "uv_coupling_running.png",
        "intermediate_direct_weinberg_running.png",
        "c5_threshold_contributions.png",
        "c5_running.png",
        "neutrino_mass_splitting_running.png",
        "neutrino_mixing_running.png",
        "oscillation_observable_sensitivity.png",
    }
    diagnostic_only = {
        "intermediate_scalar_coupling_running.png",
        "neutrino_mass_running.png",
        "mnu_matrix_element_running.png",
        "mnu_matrix_abs_ev.png",
        "pmns_abs.png",
    }

    for filename in generated:
        assert filename in generate_source

    for filename in diagnostic_only:
        assert filename not in generate_source
