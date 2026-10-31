## Experimental B/C Y-only native Matchete matching attempt (2026-10-31)

- The B(−1) and C(−1) `continuous_v1` sequential numerical EFT2/SM runs have completed with corrected muX thresholds of 1035.77775471 GeV and 1032.85523509 GeV, respectively. These remain **conditional**: finite Y/X matching absent.
- A scalar-bubble diagnostic passed 13 standalone tests. Its full-mass candidate neutrino spectra are (0.0393113530, 0.0157895632, 0.0005317403) eV and (0.0470776247, 0.0467726114, 0.0223335971) eV, respectively, but **physical matching remains uncertified**.
- The experimental `T3BCYOnlyMatcheteBuilder.py` and `RunT3BCYOnlyMatchete.py` now create and run a distinct Matchete calculation with Y heavy and X light, using the existing exact mass-eigenstate CG, EOMEFT, and `RunMatching` infrastructure. Source transformations are fail-closed and preserve the installed baseline.
- Its output is a raw one-loop EFT2 operator diagnostic (`LLHH`, `LLHX`, `LLXX`) rather than a renormalized numerical threshold coefficient. The saved native parameter point has not yet been reconciled with the running `muY` benchmark. `physical_matching_certified=False` must remain unchanged.
- 12 **offline** source-adaptation tests passed; live Matchete execution and numeric coupling handoff require a local run.

Further work: native full-EFT Y-only Matchete result -> EFT1/EFT2 infrared subtraction/counterterms -> correct boundary values at running muY -> X-only one-loop matching -> complete flavor/Wilson operator closure and numerical result.
