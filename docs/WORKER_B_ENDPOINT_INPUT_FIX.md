# Worker B endpoint audit v2: input discovery

Replaces only scripts/RunT3WorkerBEndpointAudit.py. The old script hardcoded
output/general_observable/all16_v1/<MODEL>/running.json; that file may not
exist in the user's current local project.

The v2 script searches output/**/running.json and accepts only JSON with
matching model, scales_gev.S and explicit eft1_at_scalar_threshold entries
mS1Sq/mS2Sq (plus muT3GeneralS2Sq for D_-2). This does not imply all
saved runs contain that schema. It refuses ambiguous multiple matches.
Use --state-file A_-4=PATH (repeat per model) for an explicit authoritative
input. It reports input paths and hashes, and never overwrites output.

PowerShell from repository root:

python -X utf8 -m unittest discover -s tests -p 'test_t3_worker_b_input_discovery.py' -v
python -X utf8 -m scripts.RunT3WorkerBEndpointAudit --find-inputs --models A_-4,D_-2
python -X utf8 -m scripts.RunT3WorkerBEndpointAudit --models A_-4,D_-2 --out 'output/audits/workerb_mass_endpoint_v2/report.json'

If discovery reports zero candidates, DO NOT invent a running.json. The newer
rough16 consolidated files may have a different schema and require a separately
validated adapter or exported endpoint. If more than one candidate, supply
--state-file MODEL=PATH for each model.

The audit checks quadratic endpoint consistency only, not the full scalar
potential or the source HHSS tensor.
