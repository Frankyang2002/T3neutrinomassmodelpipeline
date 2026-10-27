T3 Matchete DetailedOutput key-parser correction (27 October 2026)

Full replacement files:
  Lagrangian/RunGeneralQuadraticMatcheteAudit.wl
  scripts/RunGeneralT3QuadraticMatcheteAudit.py
  tests/python/test_general_t3_matchete_mass_basis.py
  docs/T3_MATCHETE_DETAILED_MASS_BASIS_AUDIT.md

Fix: Matchete DetailedOutput contains string-keyed Association entries.
The previous Wolfram key parser applied ToString[..., InputForm] to strings,
introducing quote characters and causing BaselineDetailedCheckFailed even
for A_+0. The new parser preserves strings, adds raw-key diagnostics, and
keeps explicit HeavyMassBasis failures separate for B_-1 and C_-1.

Install from repository root (PowerShell):
  Expand-Archive -LiteralPath 'C:\path\to\T3_Matchete_Diagnostic_Key_Fix.zip' -DestinationPath . -Force
  python -m pytest -q tests\python\test_general_t3_matchete_mass_basis.py
  python -m scripts.RunGeneralT3QuadraticMatcheteAudit --model A_+0
  python -m scripts.RunGeneralT3QuadraticMatcheteAudit --model B_-1
  python -m scripts.RunGeneralT3QuadraticMatcheteAudit --all

No GitHub changes. No production physics changes. Live Matchete check required.
