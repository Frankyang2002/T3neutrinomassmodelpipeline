# Workstream A: symbolic SU(2) reference and conditional physical normalization boundary

## Intent

Run one fail-closed, opt-in evidence audit that combines the *already implemented* Stage23–26 full signed quartic action polynomial code and the Stage27–28 native Matchete/RGBeta representation probes. This does not change the RGBeta or Matchete production code and never calls `Match[]`.

## Mathematics

For the defining doublet generators `sigma_a/2`, the normalized symmetric-pair embedding `E` has columns `|00>`, `(|01>+|10>)/sqrt(2)`, and `|11>`. Define

\[
 J_a=E^\dagger(t_a\otimes I+I\otimes t_a)E,
 \quad K=\begin{pmatrix}0&0&1\\0&-1&0\\1&0&0\end{pmatrix}.
\]

The runner verifies exactly `E E†=(I+Swap)/2`, `[J_a,J_b]=i ε_{abc}J_c`, `sum J_a²=2I`, `tr(J_a J_b)=2δ_ab`, `J_a^T K+K J_a=0`, and `K(-conj J_a)K^-1=J_a`. This is a **declared abstract spin-one reference**; it does not claim that the installed RGBeta uses these numerical component coordinates.

The stage then re-executes `scripts/QuarticActionIdentity.py`'s existing `audit()` on the saved native Matchete tensors, native RGBeta *symbolic* invariant records, Stage23 candidate coupling conversion matrices, and Stage25 actual RGBeta permutation/field-metric evidence. Each 46-block polynomial is compared to the saved Stage26 action report; all 50 quartic slots must be covered. The action identity includes formal Hermitian conjugation, correct Lagrangian-versus-potential sign and a negative `1/4` control. Those action polynomials use a **provisional common component basis** and are not independent native RGBeta numerical tensors.

## Cross-package limitations

1. RGBeta's live `delS2` and `tGen` symbolic structure is supported by the saved Stage28 evidence; native triplet numeric components were unavailable in the preceding coordinate gate.
2. Source inspection of `RGE/running/rgbeta/T3RGBetaModel.wl` shows only the baseline `y1` and `y2` Yukawa registrations. The additional six conjugate-neutral-scalar Yukawa matrix types built by generalized Matchete models are not registered in that RGBeta model source. This is a coverage gap between tools, **not a failure of the independent UV operator census**.
3. No complete native one-session Matchete/RGBeta action subtraction, full fermion normalization identification, explicit Matchete Z2 registration, or actual EFT matching is claimed.

## Paths and commands

New project files:

- `Lagrangian/RunT3UVSymbolicClosure.py`
- `docs/WORKSTREAM_A_SYMBOLIC_CLOSURE.md`

Run from the project root in the current Python virtual environment:

```powershell
python Lagrangian\RunT3UVSymbolicClosure.py --self-test
python Lagrangian\RunT3UVSymbolicClosure.py --all
```

Produces, without overwriting existing output:

- `Reports/output/workstream_a/symbolic_uv_closure/symbolic_uv_closure.json`
- `Reports/output/workstream_a/symbolic_uv_closure/symbolic_uv_closure.md`

The expected status is `ConditionalSymbolicUVClosurePassedPhysicalRGBetaDictionaryBlocked`. Its successful exit code means all supported *conditional* tests passed—not physical normalization certification.
