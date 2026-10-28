# Batch 16 — hierarchy-safe neutrino normalization probes

**Not a fitted or certified full-general-scalar neutrino calculation.**
Candidate native RGBeta UV/EFT1, existing conditional one-scale Matchete C5,
and SM+Weinberg running to 100 GeV. No UV threshold matching is changed.

Mass screen: rhoF <= 0.5; scalar masses at 1 TeV in [0.5, 5] TeV.

| Model | Yukawa seed | Y factor | Lambda multiplier | Lambda UV | rhoF | Scalars at 1 TeV (TeV) | chi2 | Status |
|---|---|---:|---:|---:|---:|---|---:|---|
| B_-1 | saved_poor_fit | 0.005 | 1 | 0.0002 | 0.014787 | 1.33–1.51 | 1.732266e+09 | CandidateOscillationEvaluated |
| B_-1 | saved_poor_fit | 0.01 | 1 | 5e-05 | 0.023762 | 1.67–2.42 | 1.136486e+09 | CandidateOscillationEvaluated |
| B_-1 | saved_poor_fit | 0.015 | 1 | 2.222e-05 | 0.033783 | 2.13–3.44 | 7.746794e+08 | CandidateOscillationEvaluated |
| B_-1 | saved_poor_fit | 0.02 | 1 | 1.25e-05 | 0.044141 | 2.64–4.5 | 5.551106e+08 | CandidateOscillationEvaluated |
| B_-1 | historical_fit | 0.005 | 1 | 0.0002 | 0.012945 | 1.11–1.32 | 20129.18 | CandidateOscillationEvaluated |
| B_-1 | historical_fit | 0.01 | 1 | 5e-05 | 0.016246 | 1.31–1.66 | 13611.89 | CandidateOscillationEvaluated |
| B_-1 | historical_fit | 0.015 | 1 | 2.222e-05 | 0.020603 | 1.59–2.1 | 8530.29 | CandidateOscillationEvaluated |
| B_-1 | historical_fit | 0.02 | 1 | 1.25e-05 | 0.025481 | 1.92–2.6 | 5436.687 | CandidateOscillationEvaluated |
| C_-1 | saved_poor_fit | 0.005 | 1 | 0.0002 | 0.015812 | 1.44–1.82 | 3.754668e+08 | CandidateOscillationEvaluated |
| C_-1 | saved_poor_fit | 0.01 | 1 | 5e-05 | 0.026126 | 2.25–3.01 | 2.197392e+08 | CandidateOscillationEvaluated |
| C_-1 | saved_poor_fit | 0.015 | 1 | 2.222e-05 | 0.037461 | 3.17–4.31 | 1.407243e+08 | CandidateOscillationEvaluated |
| C_-1 | historical_fit | 0.005 | 1 | 0.0002 | 0.015346 | 1.3–1.77 | 1.568063e+07 | CandidateOscillationEvaluated |
| C_-1 | historical_fit | 0.01 | 1 | 5e-05 | 0.02499 | 1.88–2.88 | 9175325 | CandidateOscillationEvaluated |
| C_-1 | historical_fit | 0.015 | 1 | 2.222e-05 | 0.035676 | 2.57–4.11 | 5783477 | CandidateOscillationEvaluated |

## Interpretation

- The 1/f² lambda scaling is a starting hypothesis, NOT an RG or C5 identity.
- Any chi2 <=5.86 is a seed-level candidate fit conditional on the existing one-scale C5.
- No Yukawa or lambda optimisation is performed; numeric trials are distinct from fits.
- The beta functions and physical threshold ordering have not been newly certified.
- Inputs are source-hash-checked and historical outputs remain untouched.

### B_-1: lowest computed chi2

Yukawa source `historical_fit`, factor 0.02, lambda multiplier 1, chi2 = 5436.6872.
rhoF = 0.025481, scalar masses = 1.919–2.6 TeV.

### C_-1: lowest computed chi2

Yukawa source `historical_fit`, factor 0.015, lambda multiplier 1, chi2 = 5783477.2.
rhoF = 0.035676, scalar masses = 2.571–4.106 TeV.

