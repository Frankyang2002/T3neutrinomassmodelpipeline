T3 gamma filter for ORIGINAL fitted-model envelope HTML

Replace:
Numerical/plotting/BuildFittedModelEnvelopes.py

Command from project root:
python -m Numerical.plotting.BuildFittedModelEnvelopes --scenario smallY_smallL --output Reports/output/display/fitted_model_envelopes_all.html
Start-Process "Reports/output/display/fitted_model_envelopes_all.html"

The two gamma inputs filter saved evaluated points immediately (no re-evaluation).
Original NuFIT 1sigma/3sigma, JUNO, trajectories and model sensitivity ovals retained.
The fitted centre is always displayed if its checkbox is enabled.
Extreme scan values are only visible when successfully evaluated and saved.
Preview: example_original_data_with_gamma_controls.html (demonstration only).
