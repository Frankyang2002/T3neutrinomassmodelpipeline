# Workstream C — v15 fixed-order Y->X inheritance

**Purpose:** One source-anchored physics-level resolution from the successful v14 evidence: the archived X-only calculation used the Y **tree** EFT2 as input. Does that omission invalidate the complete one-loop LLHH coefficient? At strict local dimension five and overall order `hbar`, Y's one-loop **LLHH** term is a light-field spectator through the X threshold, whereas a Y one-loop `LLXX` term inside an additional X loop is order `hbar^2`.

Input facts verified from the current GitHub evidence for B(-1)/C(-1): saved Y-tree EFT2 has no LLHH operators and every X occurrence is even, while the full saved Y EFT2 has one `hbar` LLHH term for each chirality. The local run rechecks these facts on your actual files and exports the unaltered PL/PR InputForms.

## Files and outputs

- `scripts/BuildT3BCOneLoopInheritance.experimental.py` — **new**, complete source file.
- `scripts/RunT3BCOneLoopInheritance.ps1` — **new**, single-run wrapper.
- `docs/WorkstreamC_OneLoopInheritance_v15_README.md` — this file.
- `output/audits/t3_bc_one_loop_inheritance_v15/summary.json` — summary containing explicit provenance and status flags.
- `output/audits/t3_bc_one_loop_inheritance_v15/{B_m1,C_m1}_Y_one_loop_LLHH_{PL,PR}_inherited_through_X.UNRENORMALIZED.inputform.txt` — exact unexpanded Y `LLHH` spectator coefficients. These are not renormalized Wilson coefficients.

The script also projects the native F/Y/X `LLHH` pole sectors onto the same neutral SU2 component and checks **both** formal algebraic residue relations, without assuming cancellation. A result here alone is not UV/IR subtraction. This comparison is recorded separately and cannot fail the proven spectator-inheritance calculation; any unsupported source is recorded as `NotEstablished`.

## Run (PowerShell)

```powershell
cd C:\thesistemp\T3neutrinomassmodelpipeline
Expand-Archive -LiteralPath "$HOME\Downloads\T3_WorkstreamC_OneLoopInheritance_v15.zip" -DestinationPath . -Force
& .\scripts\RunT3BCOneLoopInheritance.ps1
```

Alternative:

```powershell
python scripts\BuildT3BCOneLoopInheritance.experimental.py --self-test
python scripts\BuildT3BCOneLoopInheritance.experimental.py --output output\audits\t3_bc_one_loop_inheritance_v15\summary.json
```

## Power counting

Write the EFT2 Lagrangian as `L2 = L2^(0) + hbar L2^(1)`. For an X-odd-free (Z2-preserving) tree theory, integrating out X at tree level sets X to its zero classical solution for the relevant LLHH projection. The existing tree-source X-loop contribution is `hbar Δ_X^(1)[L2^(0)]`. The one-loop Y `LLHH` term, containing no X, is transmitted directly to EFT3 at tree level. Applying an additional X loop to the one-loop Y `LLXX` term produces order `hbar^2` and must not be counted at one loop. This does not assert that the **physical finite** X correction vanishes.

**Still blocking physical C5:** the original fermion-first EFT1 versus canonical scalar EFT algebraic equivalence, the finite UV-minus-EFT one-loop subtractions, relevant counterterms/EOM/field/parameter matching, consistent threshold running and scheme. No production files, UV inputs or RG conventions are modified.
