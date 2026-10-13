Replace Numerical/orchestration/RunFourBenchmarkCorrelations.py in the project root.

Rebuild display from existing CSV:
python -m Numerical.orchestration.RunFourBenchmarkCorrelations --display-only

Retry failed LH/HH and other failed cases, then rebuild display:
python -m Numerical.orchestration.RunFourBenchmarkCorrelations --run --retry-failed

If LH/HH still fail, the HTML status line displays the actual Python exception. Check lambdaT3=5 against model validation and numerical integration limits before changing the physics.
