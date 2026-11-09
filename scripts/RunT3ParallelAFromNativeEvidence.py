r"""Fit ALL 50 extra-quartic Matchete/RGBeta raw invariant polynomials.

The resulting dictionary is CONDITIONAL: the RGBeta input is the v2
abstract-component export and the SU(2) reference triplet CG is aligned to
the Matchete anchor.  Do not insert the result as a certified native coupling
map into the Wilson ODE until native scalar and interaction normalization has
been established, including baseline quartics and the A(-4) portal.

PowerShell, from repository root:
  .\.venv\Scripts\python.exe -m scripts.RunT3ParallelAFromNativeEvidence `
    --matchete output\parallel_t3\wilson\matchete_native_tensors.json `
    --rgbeta output\parallel_t3\wilson\rgbeta_abstract_components_conditional_v2.json `
    --output output\parallel_t3\wilson\additional_quartic_dictionary_conditional_v3.json
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys
import zipfile
from hashlib import sha256

from Numerical.parallel_t3.wilson.matchete_rgbeta_dictionary import (
    DictionaryEvidenceError, solve_from_paths, solve_from_uploaded_evidence,
)


def _from_archive(path: Path, *, tolerance: float, seed: int) -> dict:
    if not path.is_file():
        raise DictionaryEvidenceError(f"Archive missing: {path}")
    with zipfile.ZipFile(path) as z:
        expected = {"matchete_native_tensors.json",
                    "rgbeta_abstract_components_conditional_v2.json"}
        if set(z.namelist()) != expected:
            raise DictionaryEvidenceError("Expected exactly two flat JSON evidence files in archive")
        raw_m = z.read("matchete_native_tensors.json")
        raw_r = z.read("rgbeta_abstract_components_conditional_v2.json")
    m = json.loads(raw_m.decode("utf-8-sig"))
    r = json.loads(raw_r.decode("utf-8-sig"))
    return solve_from_uploaded_evidence(m,r,msha=sha256(raw_m).hexdigest(),
                                         rsha=sha256(raw_r).hexdigest(),
                                         tol=tolerance,seed=seed)


def main(argv: list[str] | None = None) -> int:
    ap=argparse.ArgumentParser(description=__doc__)
    source=ap.add_mutually_exclusive_group(required=True)
    source.add_argument("--archive", type=Path,
                        help="Archive with the two Stage22/v2 JSON files")
    source.add_argument("--matchete", type=Path,
                        help="Exact Matchete 50-slot tensor JSON")
    ap.add_argument("--rgbeta",type=Path,
                    help="Conditional v2 RGBeta reconstruction; required with --matchete")
    ap.add_argument("--output",type=Path,required=True)
    ap.add_argument("--tolerance",type=float,default=1e-9)
    ap.add_argument("--seed",type=int,default=1129)
    args=ap.parse_args(argv)
    if args.output.exists():
        ap.error(f"Refusing to overwrite existing numerical output: {args.output}")
    if args.archive is None and args.rgbeta is None:
        ap.error("--rgbeta is required with --matchete")
    if args.archive is not None and args.rgbeta is not None:
        ap.error("--rgbeta is not used with --archive")
    try:
        if args.archive is not None:
            result = _from_archive(args.archive,tolerance=args.tolerance,seed=args.seed)
        else:
            result = solve_from_paths(args.matchete,args.rgbeta,
                                      tol=args.tolerance,seed=args.seed)
    except (DictionaryEvidenceError, ValueError, OSError, json.JSONDecodeError,
            zipfile.BadZipFile) as exc:
        print(f"Blocked_NativeEvidenceDictionary: {exc}",file=sys.stderr)
        return 2
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(result,indent=2,allow_nan=False)+"\n",
                           encoding="utf-8")
    print(f"{result['additional_quartic_slots']}/50 raw invariant slots, "
          f"{result['model_count']}/16 models; fit={result['max_component_fit_residual']:.3e}; "
          f"holdout={result['max_independent_potential_holdout_residual']:.3e}")
    print(f"{result['status']} -> {args.output}")
    print("NOT a certified native RGBeta coupling normalization or full potential map.")
    return 0 if result['n_failed_models']==0 else 1


if __name__=="__main__":
    raise SystemExit(main())
