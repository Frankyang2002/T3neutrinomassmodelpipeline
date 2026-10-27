# Matchete mass-eigenstate: typed-field transport and baseline isolation

**Scope:** `B_-1`, `C_-1`, experimental `RunGeneralT3MassEigenMatchete` only.

## Observed Windows live failure

After the symbolic mixing-spurion and bridge-index batches, the live run passes original model construction, all seven additional quartic source registrations, SU(2) bridge setup, heavy eigenfield registration, and eigenfield Lagrangian assembly. The *combined* `CheckLagrangian` aborts with `DefineCG::notInvariant`; the previous isolator only identifies `Baseline11Combined`, which is not an individual operator. Neither `HeavyMassBasis` nor matching is certified.

## Source-level correction

The original source's already-evaluated scalar object is `Field[NewScalar1, Scalar, {Index[i,T3BSMd2]}, {}]`. The previous rule sent the **already typed** `Index[i,T3BSMd2]` back through the public accessor `NewScalarX[...]`, risking an `Index[Index[i,T3BSMd2],T3BSMd2]` nested index. The new rule builds the native `Field[NewScalarX,Scalar,{index},{}]` using the pre-existing typed index. It does not alter the SU(2) epsilon bridge or the numerical mass rotation. The S2-side new dummy label remains passed to the public field accessor exactly once, and the bridge receives its independently typed `Index[j,T3BSMd2]`.

The Wolfram script checks the shapes of new scalar `Field` indices immediately after the interaction substitution and again after final `RelabelIndices` on the complete Lagrangian. The preflight fails on any non-`T3BSMd2` or doubly wrapped field index, and still rejects untyped bridge CG slots.

## Improved live failure isolation

If `CheckLagrangian` still fails, the runner now checks `LSM + FreeLag[F,X,Y]` by itself, then partitions the **already validated** original `built["LInt"]` by coupling symbol and CG label, and checks these disjoint baseline partitions before any of the seven additional quartics. This avoids rebuilding/duplicating `DefineCG` registrations. The partitioner verifies its parts sum back to the exact expanded original interaction Lagrangian. It prints each tested signature, and returns the first failing one in `isolated_source_failure`.

If the baseline partition check is inconclusive because a partition is not separately Hermitian, its failure must be interpreted with the message/log, **not** as a proof that the physical interaction is gauge non-invariant. `DefineCG::notInvariant` specifically indicates a tensor issue; `CheckLagrangian == False` alone may have other causes.

## Installation paths

Full replacements:

- `Lagrangian/RunGeneralT3MassEigenMatchete.wl`
- `tests/python/test_general_t3_mass_eigen_matchete.py`

New documentation:

- `docs/T3_MATCHETE_EIGENFIELD_TYPED_INDEX_AND_BASELINE_ISOLATION.md`

No change to the Python Matchete runner, numerical mass rotation, RGEs, production Lagrangian, or source hash logic. No GitHub changes.

## Commands

```powershell
cd C:\thesistemp\T3neutrinomassmodelpipeline
Expand-Archive -LiteralPath 'C:\path\to\T3_Matchete_Eigenfield_Index_and_Isolation_Fix.zip' -DestinationPath . -Force
python -m pytest -q tests\python\test_general_t3_mass_eigen_matchete.py tests\python\test_general_t3_mixed_scalar_mass_basis.py
python -m scripts.RunGeneralT3MassEigenMatchete --model B_-1
python -m scripts.RunGeneralT3MassEigenMatchete --model C_-1
```

Expected diagnostic destinations:

- `output/general_t3_mass_eigen_matchete/B_m1_final_matchete.log`
- `output/general_t3_mass_eigen_matchete/B_m1_final_matchete.json`
- `output/general_t3_mass_eigen_matchete/C_m1_final_matchete.log`
- `output/general_t3_mass_eigen_matchete/C_m1_final_matchete.json`

Offline verification: Python test suite and standalone representative Wolfram-language index-shape and expression-partition checks. **Live Matchete has not been run on the user's environment.** A real `HeavyMassBasis` pass is not claimed. External numeric U phases, Matchete–RGBeta normalization, and two-threshold matching are still uncertified.
