# Thesis-consistent EFT1 LLSS isolation — D_-2

Same frozen denominator for both comparisons: original EFT1-off (1001).

Running = 100 × (O₁₁₁₁ − O₁₀₀₁) / O₁₀₀₁.
Mixing off = 100 × (O_LLSS-off − O₁₀₀₁) / O₁₀₀₁.

The Running numbers are copied from the source thesis stage-freeze CSV and audited against its raw observables. The LLSS-off numbers are freshly calculated.

| Model | Δm²₂₁ Running / Mixing off (%) | Δm²₃ℓ Running / Mixing off (%) | sin²θ₁₂ Running / Mixing off (%) | sin²θ₁₃ Running / Mixing off (%) | sin²θ₂₃ Running / Mixing off (%) |
|---|---:|---:|---:|---:|---:|
| A_-4 | +101064 / +174.248 | +101064 / +174.248 | +2.12636e-05 / -2.67548e-05 | -0.000414522 / +0.000521566 | -0.000174882 / +0.000220042 |
| A_-2 | +149841 / +94.7731 | +149841 / +94.7732 | +2.15716e-05 / -2.67512e-05 | -0.000420526 / +0.000521495 | -0.000177414 / +0.000220012 |
| A_+0 | +149701 / +105.59 | +149701 / +105.59 | +2.15732e-05 / -2.67539e-05 | -0.000420557 / +0.000521547 | -0.000177428 / +0.000220034 |
| B_-3 | +138007 / +99.7683 | +138007 / +99.7684 | +2.16802e-05 / -2.67364e-05 | -0.000422644 / +0.000521208 | -0.000178308 / +0.000219891 |
| B_-1 | +300.678 / +300.678 | +300.678 / +300.678 | -2.53728e-05 / -2.53728e-05 | +0.000321389 / +0.000321389 | +0.000710827 / +0.000710827 |
| B_+1 | +144967 / +97.0609 | +144967 / +97.0609 | +2.15814e-05 / -2.67364e-05 | -0.000420717 / +0.000521208 | -0.000177495 / +0.000219891 |
| C_-3 | +141217 / +79.144 | +141217 / +79.1441 | +2.12712e-05 / -2.67842e-05 | -0.000414669 / +0.00052214 | -0.000174944 / +0.000220284 |
| C_-1 | +289.322 / +289.322 | +289.322 / +289.322 | -2.41145e-05 / -2.41145e-05 | +0.000328739 / +0.000328739 | +0.000662324 / +0.000662324 |
| C_+1 | +148176 / +76.4329 | +148176 / +76.4329 | +2.11823e-05 / -2.67842e-05 | -0.000412936 / +0.00052214 | -0.000174213 / +0.000220284 |
| D_-2 | +152996 / +200.997 | +152996 / +200.997 | +2.04478e-05 / -2.66597e-05 | -0.000398618 / +0.000519713 | -0.000168176 / +0.00021926 |
| D_+0 | +155545 / +166.139 | +155545 / +166.14 | +2.0611e-05 / -2.6672e-05 | -0.000401799 / +0.000519952 | -0.000169507 / +0.00021936 |
| D_+2 | +112943 / +287.22 | +112943 / +287.221 | +2.01076e-05 / -2.66301e-05 | -0.000391987 / +0.000519135 | -0.000165378 / +0.000219017 |
| E_-4 | +67300.5 / +508.994 | +67300.5 / +508.994 | +2.08609e-05 / -2.6642e-05 | -0.00040667 / +0.000519367 | -0.000171569 / +0.000219114 |
| E_-2 | +92108.5 / +310.195 | +92108.5 / +310.195 | +2.08857e-05 / -2.66865e-05 | -0.000407154 / +0.000520234 | -0.000171774 / +0.00021948 |
| E_+0 | +96978.1 / +335.331 | +96978.1 / +335.331 | +2.06122e-05 / -2.66551e-05 | -0.000401822 / +0.000519622 | -0.000169522 / +0.000219222 |
| E_+2 | +70399.3 / +489.841 | +70399.3 / +489.841 | +2.07862e-05 / -2.6642e-05 | -0.000405214 / +0.000519367 | -0.000170955 / +0.000219114 |

## Scope and interpretation

Frozen: UV on; both renormalisable and Wilson EFT1 running off; scalar matching and SMEFT on.
Running: UV, renormalisable EFT1, Wilson EFT1, scalar matching and SMEFT all on.
LLSS off: identical to Running except the quartic LLSS→LLHH Wilson beta source is removed at each RHS evaluation.
Both scalar multiplets are integrated out together at the scalar threshold.
Exploratory only: spinor normalization, full basis map, ADM and independent finite hard matching are not certified.
Thesis-vs-new source mismatches invalidate claims of exact reproducibility; input provenance is logged.
