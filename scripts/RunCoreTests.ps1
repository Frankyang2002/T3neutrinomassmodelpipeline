# Fast project-critical regression selection. Use -Full for all remaining tests.
# This is not a proof of complete general-T3 physics matching.
param([switch]$Full)
$ErrorActionPreference = 'Stop'
$root = (Resolve-Path (Join-Path $PSScriptRoot '..')).Path
Push-Location $root
try {
    $env:PYTHONUTF8 = '1'
    if ($Full) {
        python -X utf8 -m pytest -q --tb=line
    }
    else {
        $tests = @(
            'tests/test_t3_field_content.py',
            'tests/test_t3_model_names.py',
            'tests/test_t3_model_selection.py',
            'tests/test_pipeline_plan.py',
            'tests/test_scalar_first_safety.py',
            'tests/test_eft_transitions.py',
            'tests/test_pipeline_cli.py',
            'tests/test_pipeline_backbone.py',
            'tests/test_numerical_state.py',
            'tests/test_state_vector.py',
            'tests/test_neutrino_mass_convention.py',
            'tests/test_eft_foundations.py',
            'tests/test_full_flavor_weinberg.py',
            'tests/test_weinberg_trajectory.py',
            'tests/python/test_general_coupling_inventory.py',
            'tests/python/test_full16_general_symbolic_physics_audit.py',
            'tests/python/test_general_t3_fermion_first_eft1_tree.py',
            'tests/python/test_general_t3_fermion_first_eft1_one_loop.py'
        )
        foreach ($test in $tests) {
            if (!(Test-Path -LiteralPath $test)) {
                throw "Required core test is missing: $test"
            }
        }
        python -X utf8 -m pytest -q --tb=line @tests
    }
    if ($LASTEXITCODE -ne 0) {
        throw "Pytest failed (exit code $LASTEXITCODE)."
    }
}
finally {
    Pop-Location
}
