Compact physical C5 comparison

Install ZIP into repository root.
Generated PDF/TeX: Reports/output/analytical/Lagrangian/C5.pdf and C5.tex
Exact source: C5_original_matchete.tex; common kernel: C5_common_kernel.txt
Generator: scripts/BuildCompactC5Report.py
Test: python -m pytest -q tests/test_compact_c5.py
Rebuild: python scripts/BuildCompactC5Report.py --compile

The generator extracts all 16 model coefficients from the saved Matchete report, verifies an identical class-normalised mass kernel, and displays only physical symmetric C5. It does not run matching or numerical physics. The mass kernel retains the regulator epsilon exactly as the source does. Do not run the old C5 visual polisher on this new report.
