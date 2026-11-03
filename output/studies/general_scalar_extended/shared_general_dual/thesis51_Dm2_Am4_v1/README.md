# General T3 dual benchmark — thesis Section 5.1

Benchmarks: D(-2) and A(-4) independently fitted full-general Yukawa matrices.
For each benchmark the same 3x3 real y1,y2 is used in every target model.
Other fixed UV inputs come from target reference configurations (not necessarily identical across models).
Added quartic boundary = 0.1 * |lambdaT3|; scalar bilinear boundary = 0.01 * min(mS1Sq,mS2Sq).
Additional Yukawas: exact source y2Conj copied for D(-2),
otherwise 0.05 * shared y1/y2 where that extra channel exists.
The exact-source and assumed-fraction cases are not equivalent controls.

Completed conditional evaluations: 32/32.
Blocked: 0; preflight-ready: 0.

summary.csv: C5 norm, source-relative residual, masses and oscillation observables.
six_pair_comparison.csv: A(-2)/D(-2), B(-3)/C(-3), B(-1)/C(-1) under both Yukawa benchmarks.
same_model_texture_comparison.csv: D(-2) vs A(-4) Yukawas in each fixed target model.
See per-model result.json and status.json for the full 3x3 C5 matrices and provenance.

STATUS: conditional one-scale general matching, not physically certified sequential EFT.
The physical Matchete-RGBeta coupling dictionary is uncertified.
B/C(-1) uses a separate matching evaluator; cross-adapter normalisation is uncertified.
D(0) can be blocked by an unsupported negative real self-bilinear rotation.
A large fixed-Yukawa chi2 does not exclude a model after independent Yukawa refitting.
