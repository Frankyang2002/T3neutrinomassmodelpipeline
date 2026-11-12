# Thesis-consistent EFT1 LLSS isolation — A_-4

Same frozen denominator for both comparisons: original EFT1-off (1001).

Running = 100 × (O₁₁₁₁ − O₁₀₀₁) / O₁₀₀₁.
Mixing off = 100 × (O_LLSS-off − O₁₀₀₁) / O₁₀₀₁.

The Running numbers are copied from the source thesis stage-freeze CSV and audited against its raw observables. The LLSS-off numbers are freshly calculated.

| Model | Δm²₂₁ Running / Mixing off (%) | Δm²₃ℓ Running / Mixing off (%) | sin²θ₁₂ Running / Mixing off (%) | sin²θ₁₃ Running / Mixing off (%) | sin²θ₂₃ Running / Mixing off (%) |
|---|---:|---:|---:|---:|---:|
| A_-4 | +101013 / +174.193 | +101013 / +174.192 | -5.94753e-05 / +7.48323e-05 | +0.0001297 / -0.00016319 | -0.000158258 / +0.000199121 |
| A_-2 | +149753 / +94.7294 | +149753 / +94.7291 | -6.03367e-05 / +7.48221e-05 | +0.000131578 / -0.000163168 | -0.000160549 / +0.000199094 |
| A_+0 | +149614 / +105.537 | +149614 / +105.537 | -6.03409e-05 / +7.48296e-05 | +0.000131588 / -0.000163185 | -0.000160561 / +0.000199114 |
| B_-3 | +138008 / +99.7728 | +138009 / +99.7725 | -6.06396e-05 / +7.47809e-05 | +0.000132239 / -0.000163078 | -0.000161356 / +0.000198985 |
| B_-1 | +300.681 / +300.681 | +300.682 / +300.682 | +2.40707e-05 / +2.40707e-05 | +5.57968e-06 / +5.57968e-06 | +0.000104368 / +0.000104368 |
| B_+1 | +144971 / +97.0594 | +144971 / +97.0591 | -6.0363e-05 / +7.47809e-05 | +0.000131636 / -0.000163078 | -0.00016062 / +0.000198985 |
| C_-3 | +141221 / +79.1509 | +141221 / +79.1507 | -5.94954e-05 / +7.49146e-05 | +0.000129744 / -0.00016337 | -0.000158311 / +0.00019934 |
| C_-1 | +289.327 / +289.327 | +289.328 / +289.328 | +2.27935e-05 / +2.27935e-05 | +5.84055e-06 / +5.84055e-06 | +9.6006e-05 / +9.6006e-05 |
| C_+1 | +148182 / +76.4301 | +148183 / +76.4298 | -5.92466e-05 / +7.49146e-05 | +0.000129201 / -0.00016337 | -0.000157649 / +0.00019934 |
| D_-2 | +153087 / +201.079 | +153087 / +201.079 | -5.71916e-05 / +7.45664e-05 | +0.00012472 / -0.000162611 | -0.000152185 / +0.000198414 |
| D_+0 | +155642 / +166.205 | +155642 / +166.205 | -5.76479e-05 / +7.46007e-05 | +0.000125714 / -0.000162685 | -0.000153388 / +0.000198505 |
| D_+2 | +113004 / +287.308 | +113004 / +287.308 | -5.62399e-05 / +7.44835e-05 | +0.000122645 / -0.00016243 | -0.000149652 / +0.000198194 |
| E_-4 | +67300.7 / +509.01 | +67300.8 / +509.009 | -5.83479e-05 / +7.45168e-05 | +0.000127241 / -0.000162503 | -0.000155258 / +0.000198282 |
| E_-2 | +92109.4 / +310.201 | +92109.5 / +310.2 | -5.84171e-05 / +7.46412e-05 | +0.000127392 / -0.000162774 | -0.000155443 / +0.000198613 |
| E_+0 | +96979.8 / +335.333 | +96979.9 / +335.332 | -5.7652e-05 / +7.45534e-05 | +0.000125724 / -0.000162582 | -0.000153405 / +0.000198379 |
| E_+2 | +70400.8 / +489.836 | +70400.8 / +489.835 | -5.81387e-05 / +7.45168e-05 | +0.000126785 / -0.000162502 | -0.000154701 / +0.000198282 |

## Scope and interpretation

Frozen: UV on; both renormalisable and Wilson EFT1 running off; scalar matching and SMEFT on.
Running: UV, renormalisable EFT1, Wilson EFT1, scalar matching and SMEFT all on.
LLSS off: identical to Running except the quartic LLSS→LLHH Wilson beta source is removed at each RHS evaluation.
Both scalar multiplets are integrated out together at the scalar threshold.
Exploratory only: spinor normalization, full basis map, ADM and independent finite hard matching are not certified.
Thesis-vs-new source mismatches invalidate claims of exact reproducibility; input provenance is logged.
