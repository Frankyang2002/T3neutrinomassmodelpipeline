# Worker A handover (2026-11-12)

## Scope and verified upstream evidence

Read-only upstream repository: https://github.com/Frankyang2002/T3neutrinomassmodelpipeline

- `Numerical/running/GeneralT3DynamicFermionWilson.py`: current tree LLSS source tensor axes `(flavour,flavour,SU2,SU2,scalar,scalar)`; 14 nonmixed source-factored boundaries; mixed B/C use provisional reference lift. Explicitly states `y1Conj/y2Conj` retained if in source.
- `Numerical/parallel_t3/complete/hard_llhh_source.py`: A_-4 six and D_-2 twelve candidate LLHH raw Matchete terms. The file states counterterms, EFT1 soft subtraction and GammaCC/Weyl conversion unresolved. D_-2 extra quartic channel proportional to `y2Conj`, inconsistent with restricted sector unless set to zero and consistently projected.
- `Numerical/parallel_t3/wilson/general_batch2_laurent.py`: formal Laurent separation, explicit status `FullOneLoopFormalLaurentAudited_TreeMassRemainderPreserved_NotFiniteWilsons`.
- `Numerical/parallel_t3/complete/full_general_completion.py`: requires independent renormalised matching source rather than raw Laurent term.

No accessible local working repository or output/Matchete artifacts in this sandbox. Public upstream can differ from local copy. No claim made that a final local finite match does not exist.

## New files

- `Numerical/running/T3FermionOneLoopBoundary.py`: independent safe tree + finite matching interface with certification guards; it does *not* produce physical finite loop arrays.
- `scripts/AuditT3FermionOneLoopBoundary.py`: read-only evidence discovery; synthetic smoke tests for A_-4/D_-2; optional exclusive-new JSON report.

## PowerShell (from local repository root)

```powershell
cd C:\thesistemp\T3neutrinomassmodelpipeline
Expand-Archive -LiteralPath "C:\path\to\worker_a_fermion_boundary_guarded.zip" -DestinationPath . -Force
python -m py_compile Numerical\running\T3FermionOneLoopBoundary.py scripts\AuditT3FermionOneLoopBoundary.py
$runId = Get-Date -Format 'yyyyMMdd_HHmmss'
python scripts\AuditT3FermionOneLoopBoundary.py --repo . --output "output\audits\worker_a_fermion_$runId.json"
```

Do not connect this module to the coupled solver before checking the local current code, finite soft-subtracted Matchete EFT1 export, source-to-Weyl normalisation and Yukawa restriction. Synthetic smoke test is not an A_-4/D_-2 numerical benchmark. Existing outputs remain untouched.
