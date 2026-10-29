# General T3: exact RGBeta invariant geometry diagnostic

## Scope

A new **non-certifying** source-to-tensor comparison is implemented by
`scripts/RGBetaExplicitBasis.py`. It consumes existing native outputs:

- `output/general_t3_stage22/matchete_tensors.json`
- `output/general_t3_stage22/rgbeta_native_declarations.json`
- `output/general_t3_stage22/rgbeta_native_invariants.json`

It writes `output/general_t3_stage22/stage22_explicit_geometry_candidates.json`.
**It NEVER writes `rgbeta_native_tensors.json` and NEVER changes the production
RGEs, matching, or scalar-sector conventions.** The original Stage 22
certification runner remains fail-closed.

## What is actually checked

- All 16 model and 50 slot inventories and consistency of two independent RGBeta
  metadata exports.
- Field species and ordering against the Matchete tensor exporter (no guessed
  symbol pairing for reused RGBeta couplings).
- Exact, component-by-component evaluation of the *literal* exported RGBeta
  invariant expressions built from `del[SU2L[fund],...]`,
  `del[SU2L[S2],...]`, and `eps[SU2L[fund],...]`.
- Exact rank, span, Bose symmetries, and change-of-basis matrices via the
  existing RGBetaTensorAlgebra solver. No floating point fitting.
- An explicit diagnostic convention `eps[0,1]=+1`, ordinary Kronecker delta
  components and matching native index labels. The comparison does **not**
  independently establish which complex/daggered SU(2) representation basis
  the packages use. The native RGBeta quartic function `6*Matrix[lambda]`
  is **not** interpreted as a verified scalar-potential coefficient.
- Triplet Clebsch-Gordan `delS2` structures are blocked pending a common
  representation basis: even when algebraic expressions look similar, this
  script does not insert a guessed `3 x 3` unitary basis transformation.

The exact symbolic comparisons for the uploaded Fix 2 diagnostic evidence yield:

| Category | Coupling slots | Details |
|---|---:|---|
| `ExactGeometryCandidate` | 32 | Zero component residual in the stated provisional common component labels |
| `ExplicitGeometryMismatch` | 12 | B(-1)/C(-1), five field-content blocks per model, including a 2D mixing block (two slots) |
| `TripletBasisPending` | 6 | One H/H/triplet invariant per A(+0), A(-2), D(+0), D(-2), E(+0), E(-2) |
| Total | 50 | 46 field-content blocks |

The two **E** blocks have a nontrivial formal basis change. Define `R = M A`,
where `M` and `R` are matrices whose columns are the Matchete and RGBeta
invariant components in `I1, I2` order. Then for E(+0)/000211 and E(-2)/001102:

```
A = [[1/2,               1],
     [sqrt(33)/6,         0]]
```

This is **not** a map between physical Matchete and RGBeta scalar couplings.
For that, independent signed scalar-potential prefactors and basis orientation
must be established (as well as triplet basis alignment in the remaining blocks).

### Significant warning for the B/C 2x2 block

The RGBeta `110101` invariant pair is linearly independent but falls outside
Matchete's 2D tensor span under the simple index-identification convention.
This can indicate a mismatch in raising/lowering conjugate SU(2) indices or a
problem with one generated invariant; the diagnostic cannot distinguish the two.
**Do not silently coerce or discard this mismatch.** Before certification,
verify the actual gauge invariance of both registered quartics under the
package's native representation and conjugation conventions.

## Run from repository root (PowerShell)

```
python -m pytest -q tests\python\test_general_t3_tensor_explicit_basis.py
python -m scripts.RGBetaExplicitBasis
```

Do not change `RunGeneralT3ScalarTensorCertification.py` to accept this provisional report as
native RGBeta component evidence. The old Stage 22 summary remains
`Uncertified`, and the next dependency is **independent native RGBeta tensor
basis and physical Lagrangian normalisation certification**.

## Source evidence

- Native Matchete tensor components: uploaded `matchete_tensors.json` (50 slots).
- Native RGBeta fields, invariant functions and `6*Matrix[...]` declaration
  functions: uploaded `rgbeta_native_declarations.json` and
  `rgbeta_native_invariants.json` (50 slots).
- RGBeta `Package/FieldsAndCouplings.m` `AddQuartic` and `Lam` source, read-only,
  defines the internal quartic normalisation factor and permutation machinery.
- RGBeta `Package/GroupsAndIndices.m` source, read-only, describes `delS2`
  contraction reductions. These definitions do not directly supply a certified
  matrix relating RGBeta triplet basis to Matchete's triplet component basis.
