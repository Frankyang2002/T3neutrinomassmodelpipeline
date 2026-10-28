# Batch 15 — UV scalar mass feasibility (B/C)

Candidate UV/EFT1 RGBeta only; this is **not a refit** or certified matching.
Yukawa factors multiply all 18 saved fitted real Yukawas, and scalar mass factors multiply
the original positive UV scalar masses (so m² is multiplied by the factor squared).
MF and general-quartic boundary convention are held fixed.

Screen: mSmax(muF)/MFmin(muF) <= 0.5, and all scalar running masses at muS in [0.5, 5] TeV.
A screen pass is NOT an accepted oscillation fit; chi2 is NOT recomputed.

| Model | Y factor | UV scalar mass factor | rho(muF) | Scalar masses at muS (TeV) | Screen |
|---|---:|---:|---:|---|---|
| B_-1 | 1 | 1 | 1.722 | not run / failed | CandidateScreenFail |
| B_-1 | 0.005 | 1 | 0.01479 | 1.33–1.51 | CandidateScreenPass |
| B_-1 | 0.01 | 1 | 0.02376 | 1.67–2.42 | CandidateScreenPass |
| B_-1 | 0.015 | 1 | 0.03378 | 2.13–3.44 | CandidateScreenPass |
| B_-1 | 0.02 | 1 | 0.04414 | 2.64–4.5 | CandidateScreenPass |
| B_-1 | 0.025 | 1 | 0.05464 | 3.18–5.57 | CandidateScreenFail |
| B_-1 | 0.03 | 1 | 0.06522 | 3.73–6.64 | CandidateScreenFail |
| B_-1 | 0.05 | 1 | 0.1078 | 6.01–11 | CandidateScreenFail |
| C_-1 | 1 | 1 | 1.732 | not run / failed | CandidateScreenFail |
| C_-1 | 0.005 | 1 | 0.01581 | 1.44–1.82 | CandidateScreenPass |
| C_-1 | 0.01 | 1 | 0.02613 | 2.25–3.01 | CandidateScreenPass |
| C_-1 | 0.015 | 1 | 0.03746 | 3.17–4.31 | CandidateScreenPass |
| C_-1 | 0.02 | 1 | 0.04911 | 4.12–5.65 | CandidateScreenFail |
| C_-1 | 0.025 | 1 | 0.0609 | 5.09–7.01 | CandidateScreenFail |
| C_-1 | 0.03 | 1 | 0.07275 | 6.07–8.37 | CandidateScreenFail |
| C_-1 | 0.05 | 1 | 0.1204 | 10–13.9 | CandidateScreenFail |

## UV mass-squared beta dependence at fixed state

Values are d(mS_i²)/dln(mu), in GeV²; signed terms may cancel.
The inclusion/exclusion varies y1,y2 and MF *only during expression evaluation*,
without changing any solved trajectory. This is an algebraic diagnostic, not a
unique diagrammatic decomposition.

### B_-1: uv_start
| Mass | Full | No Y/MF | Y only (MF=0) | MF only (Y=0) | Joint Y–MF |
|---|---:|---:|---:|---:|---:|
| mS1Sq | -2.7709e+09 | -9774.3 | 38036 | 0 | -2.7709e+09 |
| mS2Sq | -9.5897e+09 | -6586.7 | 1.0235e+05 | 0 | -9.5898e+09 |

### B_-1: uv_at_muF
| Mass | Full | No Y/MF | Y only (MF=0) | MF only (Y=0) | Joint Y–MF |
|---|---:|---:|---:|---:|---:|
| mS1Sq | -1.1339e+09 | 2.0503e+08 | 1.6924e+08 | 0 | -1.5082e+09 |
| mS2Sq | -1.7888e+09 | 7.7117e+08 | 1.1974e+09 | 0 | -3.7574e+09 |

### C_-1: uv_start
| Mass | Full | No Y/MF | Y only (MF=0) | MF only (Y=0) | Joint Y–MF |
|---|---:|---:|---:|---:|---:|
| mS1Sq | -1.2661e+10 | -9774.3 | 1.6671e+05 | 0 | -1.2662e+10 |
| mS2Sq | -6.6378e+09 | -6586.7 | 71390 | 0 | -6.6379e+09 |

### C_-1: uv_at_muF
| Mass | Full | No Y/MF | Y only (MF=0) | MF only (Y=0) | Joint Y–MF |
|---|---:|---:|---:|---:|---:|
| mS1Sq | -3.8601e+09 | 1.3124e+09 | 2.3658e+09 | 0 | -7.5383e+09 |
| mS2Sq | -3.0892e+09 | 6.784e+08 | 9.9721e+08 | 0 | -4.7648e+09 |

## Interpretation

- UV masses, heavy-fermion masses and threshold ordering are RUNNING parameters, not pole masses.
- EFT1 is deliberately skipped for inverted-ordering UV endpoints. No scalar-first matching is performed.
- Scanning Yukawa magnitudes changes neutrino observables; candidate points MUST be refitted with the five-observable chi2.
- Baseline UV evolution is replayed with the source-fingerprint-validated fitted Yukawas.
- If the prior Batch 14 diagnostic exists, its 100 TeV scalar endpoint is independently reproduced.
- The additional quartic mapping and sequential general-potential C5 matching remain uncertified.
- If no candidate passes, the result is limited to this finite positive-UV-mass/Yukawa grid; no no-go theorem follows.
- No historical or previous general-scalar output is overwritten.

