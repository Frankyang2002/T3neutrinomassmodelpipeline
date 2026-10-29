# B/C coupled Wilson parser API correction

This replaces the whole `Numerical/experimental/T3BCCoupledWilsonTensor.py` from the original coupled-Wilson batch. The public parser API in `RGE/running/intermediate/IntermediateMatcheteParsing.py` is `parse_mathematica_scalar`, not `_parse_mathematica_scalar`. The `_matchete_complex` helper now imports and calls the public function. No Wilson RGE, normalization, matching formula, or output path is changed.

Extract the ZIP into the repository root and overwrite the existing `Numerical/experimental/T3BCCoupledWilsonTensor.py`.

```powershell
python -X utf8 -m pytest -q tests/test_t3_bc_coupled_parser_api.py tests/test_t3_bc_coupled_wilson.py
python -m scripts.RunT3BCCoupledWilsonSequential --model B_-1 --run-id continuous_v1 --resume
python -m scripts.RunT3BCCoupledWilsonSequential --all --run-id continuous_v1 --resume
```

The failure occurred before the numerical calculation; rerunning with `--resume` preserves existing outputs. Further runtime/interface problems may still arise and must be resolved from their actual tracebacks.
