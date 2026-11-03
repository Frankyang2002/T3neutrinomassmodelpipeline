# B2 internal RGBeta tensor-normalisation witness — pre-run handoff

Source project commit inspected: `e2a5def5326cc452f668ef593bd81ae13fdff30b`.
Relevant unchanged tracked blobs are rechecked by the live Python runner;
upstream RGBeta implementation source observed at
`FieldsAndCouplings.m` blob `d524fa948fdbcf70ea5a61358b13a78b16a38faa`.

Derived at the level of **RGBeta registry tensors**, not physical action:
`y1=y2: Sqrt[2]`; `MF: 2`; `lambdaT3=lambdaH1=lambdaH2: 6`.
`Coupling` and `CouplingBar` signs/orientations are to be confirmed by the
fresh-kernel native script. For restricted five-class representatives the
probe intends 5 × 6 = 30 checks. Python tests are synthetic and do not
substitute for native execution. The conditional B/C Ward mismatch and
portal rank-deficiency results are preserved unmodified.

**Physical source certificates: 0/16**. No signed physical couplings, no
SMEFT matching, no source code updates. Not safe to integrate as a physical
coupling dictionary. See the full derivation and run the native script first.
