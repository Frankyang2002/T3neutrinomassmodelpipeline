# LLSS-off initial parity-check correction

Replace `scripts/RunT3OptimisedYukawaLLSSOff.py` with the complete file included here.

The simultaneous integrators perform their mu_F scalar-pair source parity check after entering the integration wrapper but before the first RHS evaluation. The previous wrapper subtracted the LLSS source during that check, causing a relative parity mismatch of 1.0. This version preserves the first beta evaluation, then subtracts the source on subsequent ODE evaluations. It does not change saved donor inputs, benchmark definitions, or GitHub.

Because the previous attempt saved Blocked entries, use a fresh counterfactual run ID. See the accompanying chat response for PowerShell commands.

Syntax check passed; full end-to-end physics run not verified in this environment.
