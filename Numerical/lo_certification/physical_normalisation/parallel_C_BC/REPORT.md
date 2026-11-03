# Track C — C4 final four-model B/C Dirac batch report

**Scope:** only B_-3, B_+1, C_-3, C_+1. No B/C Majorana, A/D/E, B2, RGE, EFT or GitHub changes. One shared native action extractor and one shared independent numerical reviewer.

**Execution classification:** The four numerical matrices below were **executed** using a source-shaped transcription of the four GitHub fitted-reference Yukawa matrices, masses and archived Matchete CG/portal tensors. The real local checkout and installed Matchete were **not** run here. The report does not certify any charged Dirac model physically. New local Matchete witness extraction remains required. Frozen comparisons performed only after the independent calculation, without coefficient fitting.

## C0 — source and benchmark

- GitHub source tree `d2163e6bb95109c0fac419f7e41f335f0d347b62` inspected read-only. File blob IDs in `GITHUB_SOURCE_PROVENANCE.json`. All four model definitions and six physics source files inspected; prior B1 / BM1 / CM1 witnesses inspected but not imported as numerical engine.
- `Y(S1)=alpha/2`, `Y(S2)=(alpha+2)/2`, `Y(F)=(alpha+1)/2`; SU(2) dimensions B=(2,2,1), C=(2,2,3). Both scalar doublets mixed; all F charged Dirac. Three benchmark diagonal masses `(100000,110000,120000) GeV`; `v=246.22 GeV`, `mS1^2=1306122.448979592 GeV²`, `mS2^2=1000000 GeV²`, `lambdaT3=0.1`, `lambdaH1=0.02`, `lambdaH2=0.03`.

## C1 — native contraction and quadratic normalisation

- Verified from existing read-only Matchete CG JSON: `C1_B=+1`, `C2_B=+1`, `C1_C=+3^(-1/4)`, `C2_C=-3^(-1/4)` on charge-selected `L[1]` with `S1[1]` and `Bar[S2[2]]`. For C the contributing Cartesian fermion axis is `F[2]`, which has `T3=0` and electric charge `Q(F)=±1`. There is **one**, not three, mixed mass loop channels per C model.
- Native Yukawa source `y1: Bar[L] CConj[F] S1`, `y2: Bar[L] F Bar[S2]`, each `+ H.c.`; fermion `-MF[r] Bar[F_r] F_r` (source contract). In triplet Cartesian representation the Dirac mass metric is `delta_ab`; in weight basis `U†U=I3`. `U^T U` is a possible real-representation Majorana pairing, **not** the charged-Dirac mass contraction.
- `H[2]=v/sqrt(2)` and `Cmix[2,2,1,2]=-2/sqrt(3)` give `A_(S2[2],S1[1])=+lambdaT3*v²/sqrt(3)=3500.14492271695 GeV²`. All other H0 mixing entries vanish.
- Full native *four-index* portal identity is `I_H⊗I_S = (1/2) Inv1 + (sqrt(3)/2) Inv2` for **both** scalars, obtained by a two-column least-squares solve over all 16 components; residual ~1e-16. After H0, `Inv1=diag(0,1)`, `Inv2=diag(2,1)/sqrt(3)`. Identity-portal assignment: for i=1,2, `lambdaHiInv1=lambdaHi/2`, `lambdaHiInv2=sqrt(3)*lambdaHi/2`. This is a **chosen benchmark direction**, not an automatic native name identification.
- Native zero-portal point: set all four portal invariant couplings zero. The full 4-complex-scalar Hermitian Hessian has one mixed 2×2 block (Q=-1 or +1) plus two spectator 1×1 blocks (Q=-2,0 or Q=0,+2), all strictly positive in both points.

## C2 — independent direct one-loop formula

With conditional Weyl convention `g1_ar=-conj(C1*y1_ar)`, `g2_ar=-conj(C2*y2_ar)`, `A21=+lambdaT3*v²/sqrt(3)` and mixed scalar eigenvalues `s_-<s_+`,

\[
\mathcal D_r=-\frac{1}{s_+-s_-}\int_0^1 dx\,\log\left(1+\frac{(1-x)(s_+-s_-)}{xM_r^2+(1-x)s_-}\right),
\]

\[
(m_\nu)_{ab}=\sum_{r=1}^3\frac{M_rA_{21}\mathcal D_r}{16\pi^2}\left(g_{1,ar}g_{2,br}+g_{2,ar}g_{1,br}\right).
\]

The quadrature uses `mpmath` with 55 decimal digits and independently matches the difference of finite analytic `B0(0,M²,s)` functions. The relative factor comes from actual source CG tensors, *not* the provisional mass. Yukawas are the fitted-reference arrays, not the comparison-point schematic textures.

## C3 — four-model numerical validation (conditional identity portal, three diagonal MF)

| Model | Native evidence | Scalar normalisation | Fermion normalisation | Takagi masses, descending (eV) | Provisional Frobenius residual | Remaining physical blockers |
|---|---|---|---|---|---:|---|
| `B_-3` | archived Matchete tensor; fresh native pending | full 4-component Hessian positive, identity-map conditional | charged Dirac singlet, diagonal MF source-convention conditional | 0.04833239942, 0.01098144853, 0.007218230529 | 5.978e-12 | live four-model Matchete; Track-A GammaCC; identity portal; numeric MF assignment |
| `B_+1` | archived Matchete tensor; fresh native pending | full 4-component Hessian positive, identity-map conditional | charged Dirac singlet, diagonal MF source-convention conditional | 0.04807431857, 0.0109461042, 0.00721582712 | 7.278e-12 | live four-model Matchete; Track-A GammaCC; identity portal; numeric MF assignment |
| `C_-3` | archived Matchete tensor; fresh native pending | full 4-component Hessian positive, identity-map conditional | charged Dirac triplet, diagonal MF source-convention conditional | 0.0550407444, 0.02819228569, 0.02694687434 | 6.161e-12 | live four-model Matchete; Track-A GammaCC; identity portal; numeric MF assignment |
| `C_+1` | archived Matchete tensor; fresh native pending | full 4-component Hessian positive, identity-map conditional | charged Dirac triplet, diagonal MF source-convention conditional | 0.04778610821, 0.008268705382, 0.0003654877874 | 1.268e-11 | live four-model Matchete; Track-A GammaCC; identity portal; numeric MF assignment |

### Explicit independent mass matrices

Units: eV, conditional identity portals and three fitted Dirac masses. All imaginary parts are zero at the original real fitted-reference point.

**B_-3:**

\[
m_\nu=\begin{pmatrix}
+0.002641601890 & -0.001311924273 & +0.010733181921 \\
-0.001311924273 & +0.020649561432 & +0.026350851187 \\
+0.010733181921 & +0.026350851187 & +0.021278018096
\end{pmatrix}\;\mathrm{eV},
\]

**B_+1:**

\[
m_\nu=\begin{pmatrix}
+0.002645111785 & -0.001326686229 & +0.010692133969 \\
-0.001326686229 & +0.020540856815 & +0.026213301494 \\
+0.010692133969 & +0.026213301494 & +0.021158072888
\end{pmatrix}\;\mathrm{eV},
\]

**C_-3:**

\[
m_\nu=\begin{pmatrix}
+0.010925307259 & +0.013755376417 & -0.022132463966 \\
+0.013755376417 & +0.023303758646 & +0.032764642179 \\
-0.022132463966 & +0.032764642179 & +0.019566267141
\end{pmatrix}\;\mathrm{eV},
\]

**C_+1:**

\[
m_\nu=\begin{pmatrix}
-0.003327026578 & +0.007501734156 & -0.002206454075 \\
+0.007501734156 & -0.024376141599 & +0.020508996734 \\
-0.002206454075 & +0.020508996734 & -0.027986157627
\end{pmatrix}\;\mathrm{eV},
\]

### Validation ledger

- **Pole cancellation:** for the mixed pair, `sum_k U_(S2,k) U*_(S1,k)=0`, so the common B0 UV pole cancels before taking the finite difference; measured off-diagonal pole residual zero to numerical precision in all four.
- **Scalar positivity:** smallest positive eigenvalue over all four scalar components at selected identity portals is `999050.6604968358 GeV²` in each model. Spectator charges are separate blocks and explicitly checked.
- **Symmetry and Takagi:** `||m-m^T||/||m||` < 1e-13; full numerical `U^T m U` Takagi reconstruction < 3e-8 (full exact errors in preview JSON).
- **Complex rephasing:** nonuniform source Yukawa phases + complex portal phase, followed by global shifts `(theta_y1,theta_y2,theta_lambda)=(0.21,-0.17,0.43)` satisfy `m_new=exp(i*0.39)*m_old`, with relative matrix error ≤ 3.5e-16. Takagi masses invariant.
- **Negative controls:** `lambdaT3=0` gives exactly zero mass; the wrong scalar-charge mixing component vanishes; wrong complex Yukawa conjugation changes the result by 0.98%–35%; changing only a relative Yukawa vertex sign flips the matrix but not Takagi masses, illustrating why masses alone do not certify the spinor convention.
- **Zero-portal versus identity-portal:** distinct calculations for both, and for common native MF versus three benchmark diagonal masses; see `C2.points` in JSON. The frozen comparison uses only the selected identity-portal / three-diagonal-MF point.
- **Provisional comparison:** performed last, using `FROZEN_REFERENCE_READONLY_SNAPSHOT.json`, relative residual definition `||m_ind-m_provisional||_F/||m_provisional||_F`; no fitted scale, phase or sign.

## C4 — single ZIP / exact PowerShell commands

All generated project-relative files:

```text
scripts/RunPhysicalCBCDNative.wl
scripts/ReviewPhysicalCBCDIndependent.py
Numerical/lo_certification/physical_normalisation/parallel_C_BC/ROADMAP.md
Numerical/lo_certification/physical_normalisation/parallel_C_BC/REPORT.md
Numerical/lo_certification/physical_normalisation/parallel_C_BC/GITHUB_SOURCE_PROVENANCE.json
Numerical/lo_certification/physical_normalisation/parallel_C_BC/FROZEN_REFERENCE_READONLY_SNAPSHOT.json
output/general_observable/physical_normalisation/parallel_C_BC_20261103_132040/track_c_four_models.json
output/general_observable/physical_normalisation/parallel_C_BC_20261103_132040/track_c_four_rows.csv
output/general_observable/physical_normalisation/parallel_C_BC_20261103_132040/EXECUTION_CONTEXT.md
```

From PowerShell, first inspect the target paths and refuse overwrites:

```powershell
cd C:\thesistemp\T3neutrinomassmodelpipeline
$zip = Join-Path $HOME "Downloads\TrackC_BC_20261103.zip"  # use the actual downloaded ZIP filename
Add-Type -AssemblyName System.IO.Compression.FileSystem
$archive = [System.IO.Compression.ZipFile]::OpenRead($zip)
try {
  foreach ($entry in $archive.Entries) {
    if ($entry.FullName.EndsWith("/")) { continue }
    $target = Join-Path (Get-Location).Path ($entry.FullName.Replace("/","\"))
    if (Test-Path -LiteralPath $target) { throw "Refuse overwrite: $target" }
  }
} finally { $archive.Dispose() }
Expand-Archive -LiteralPath $zip -DestinationPath .
python -m pip install numpy sympy mpmath
```

New timestamped local installed-Matchete extraction followed by a four-model review:

```powershell
$stamp = Get-Date -Format "yyyyMMdd_HHmmss"
$base = "output\general_observable\physical_normalisation\parallel_C_BC_$stamp"
if (Test-Path -LiteralPath $base) { throw "Existing batch output: $base" }
New-Item -ItemType Directory -Path "$base\native" -Force | Out-Null
$failed = @()
foreach ($m in @("B_-3","B_+1","C_-3","C_+1")) {
  wolframscript -file scripts\RunPhysicalCBCDNative.wl $m "$base\native\${m}_native_action.json"
  if ($LASTEXITCODE -ne 0) { $failed += $m; Write-Warning "Blocked native: $m" }
}
python -X utf8 scripts\ReviewPhysicalCBCDIndependent.py `
  --repo-root . --native-dir "$base\native" --source-only `
  --frozen-reference "Numerical\lo_certification\physical_normalisation\parallel_C_BC\FROZEN_REFERENCE_READONLY_SNAPSHOT.json" `
  --output-dir "$base\review"
if ($LASTEXITCODE -ne 0) { throw "Track C numerical batch blocked; inspect $base\review\track_c_four_models.json" }
Import-Csv "$base\review\track_c_four_rows.csv" | Format-Table model,native_evidence,scalar_normalisation,fermion_normalisation,Takagi_masses_eV,comparison_residual -AutoSize
if ($failed.Count -gt 0) { Write-Warning "Native failures: $($failed -join ", "); archived-only results are not native-certified" }
```

For a strict all-native acceptance (no archived-only fallback), omit `--source-only`. The reviewer never opens old production engine or frozen mass before calculating each independent matrix. `wolframscript` must resolve to the installed Matchete-capable kernel.

### Remaining blockers and stop

1. **Not executed here:** newly compiled charged-Dirac `CheckLagrangian` and live native CG/portal extraction for the four models; users must run the exact PowerShell commands.
2. **Native source convention:** one diagonal `MF[r]` action is tested, but three numeric MF eigenvalues are not substituted inside Matchete; off-diagonal mass matrices are not certified.
3. **Spinor normalisation:** `GammaCC` absolute convention and shared four-spinor↔Weyl phase remains Track A-owned, despite complex rephasing tests. Never declare full physical certification on agreement alone.
4. **Scalar portal:** benchmark identity is a *selected* exact Inv1/Inv2 combination; it is not automatically equivalent to a single native `lambdaH` term. Other symmetry-allowed interactions are conditionally set to zero.

**Batch status:** four independent *conditional* matrices and spectra produced; four native action extractions pending on installed Matchete; **0 newly physically certified**. **STOP at C4.**
