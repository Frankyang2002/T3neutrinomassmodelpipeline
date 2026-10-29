"""Collect accepted frozen T3 fits into a single CSV/JSON summary.

No physics is rerun. This is a convenience step for display/reporting.
"""
from __future__ import annotations

import argparse
import csv
import json

from Numerical.plotting import BuildFrozenGeneralScalarLambdaReports as reports


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--run-ids', default='frozen_massnorm_a0_v1,frozen_massnorm_remaining_v1')
    ap.add_argument('--output-dir', default='output/studies/general_scalar_extended/frozen_display/accepted_summary_v1')
    args = ap.parse_args(argv)

    out = reports.ROOT / args.output_dir
    out.mkdir(parents=True, exist_ok=True)
    run_ids = [x.strip() for x in args.run_ids.split(',') if x.strip()]
    rows, provenance = reports.summary_rows(reports.ROOT / reports.OUT_BASE, run_ids)

    with (out / 'summary.json').open('w', encoding='utf-8') as stream:
        json.dump({'rows': rows, 'provenance': provenance}, stream, indent=2)
    with (out / 'summary.csv').open('w', newline='', encoding='utf-8') as stream:
        writer = csv.DictWriter(stream, fieldnames=['model', 'case', 'chi2', *reports.OBS, 'run_id', 'status', 'source_path'])
        writer.writeheader()
        writer.writerows(rows)
    print(f'Wrote {out} | {len(rows)} rows')
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
