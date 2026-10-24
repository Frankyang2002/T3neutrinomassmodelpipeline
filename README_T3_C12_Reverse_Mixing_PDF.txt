T3 C12 / C5 GroupFactors PDF: reverse mixing audit and duplicate removal
=======================================================================

INSTALL (PowerShell)
--------------------
Extract the ZIP in C:\thesistemp\T3neutrinomassmodelpipeline, preserving
its folder structure. This replaces only the report-source writer and its
previous regression test; the physics master tensor RGE is unchanged.

FILES
-----
Full replacement:
  Reports/C12GroupFactorPDF.py
  tests/test_c12_c5_mixing_pdf.py
New:
  RGE/group_factors/validation/C12ReverseMixingGroupFactors.py
  tests/test_c12_reverse_mixing_group_factors.py
  README_T3_C12_Reverse_Mixing_PDF.txt

WHAT CHANGED
------------
- GF_EFT_1_after_F.pdf now shows forward C12->C5 once, in the existing R_W
  table. The duplicate forward row is removed from the C12 block.
- The C12 block retains only lambda5* C5 -> beta_C12.
- The existing T3-B alpha=-1 factor 8*sqrt(3)/3 is preserved.
- Other reverse coefficients can be populated from exact, C5-seeded tensor
  projections saved beside each representation's model RGE outputs.
- The audit tests collinearity over all independent LLSS tensor components,
  including components outside the matched C12 template. Nonzero residuals
  are not silently converted to a group factor.

TESTS
-----
python -m pytest -q tests/test_c12_c5_mixing_pdf.py tests/test_c12_reverse_mixing_group_factors.py tests/test_c12_group_factor_pdf.py tests/test_c12_group_factor_block_tables.py

FIRST VERIFY THE REFERENCE T3-B ALPHA=-1 CASE
--------------------------------------------
python -m RGE.group_factors.validation.C12ReverseMixingGroupFactors `
    --data-dir "output\T3_B_alpha_m1\data"

Expected, if the local tensor source has the same normalization as the
benchmark: status Verified, factor 8*sqrt(3)/3. If the status differs, read
the emitted c12_reverse_c5_group_factor.json and do not treat the extracted
factor as verified. The audit does not overwrite the established T3-B result.

EXTRACT ALL REPRESENTATIONS FROM THE SAVED ANALYTICAL STUDY
----------------------------------------------------------
python -m RGE.group_factors.validation.C12ReverseMixingGroupFactors `
    --source "output\full\comparison\smallY_smallL"

The audit writes one cache file in each physical run's data directory:
    data/c12_reverse_c5_group_factor.json
It uses existing matched Wilson, quartic and RGBeta seed files, but DOES
recompute a symbolic pure-C5 master-beta insertion. It does not rerun
Matchete, RGBeta, fitting, or the numerical RGE pipeline.

REGENERATE EXISTING PDF
-----------------------
python -m Reports.C12GroupFactorPDF `
    --tex "Reports\output\analytical\GroupFactors\GF_EFT_1_after_F.tex" `
    --source "output\full\comparison\smallY_smallL" `
    --compile

Start-Process "Reports\output\analytical\GroupFactors\GF_EFT_1_after_F.pdf"

NOTES
-----
- Verified means that the one-generation lambda5-star-generated LLSS tensor
  matches the Matchete C12 direction exactly, including zero-seed components.
- Not verified means no matching successful saved reverse audit, or a failed
  proportionality/normalization check. It does NOT mean the factor is zero.
- Other generated operator channels and full 3-generation operator closure
  remain separate questions.
- Other T3 representations are never assigned the T3-B factor by assumption.
