# Batch 1 execution gate

No source files modified.

| Component | Present |
|---|---|
| `general_uv` | yes |
| `extended_c5` | yes |
| `low_energy` | yes |
| `c5_observables` | yes |
| `production` | yes |
| `coupled_b` | yes |
| `mass_observables` | yes |
| `charged_basis` | yes |
| `eft1_wilson` | yes |
| `general_state` | yes |

## Focused test result

- `tests/python/test_general_t3_eft1_hard_pole_audit.py`: return code 1
- `tests/python/test_general_t3_quartic_action_identity.py`: return code 1
- `tests/python/test_general_t3_triplet_intertwiner.py`: return code 1
- `tests/python/test_general_t3_triplet_native_structure.py`: return code 1
- `tests/python/test_general_t3_coupled_uv_running.py`: return code 0
- `tests/python/test_general_t3_eft1_scalar_pair_mixing.py`: return code missing/error

## Scientific certification

Not certified by this gate. Inspect JSON for explicit gaps.
