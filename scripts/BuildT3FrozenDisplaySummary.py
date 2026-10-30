from __future__ import annotations

import argparse
import csv
import json

from Numerical.plotting.BuildFrozenGeneralScalarLambdaReports import summary_rows, ROOT


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description='Collect accepted frozen fits into one summary.')
    ap.add_argument('--run-ids', default='frozen_massnorm_a0_v1,frozen_massnorm_remaining_v1')
    ap.add_argument('--output-dir', default='output/studies/general_scalar_extended/frozen_display/accepted_summary_v1')
    args = ap.parse_args(argv)
    out = ROOT / args.output_dir
    out.mkdir(parents=True, exist_ok=True)
    rows, provenance = summary_rows([x.strip() for x in args.run_ids.split(',') if x.strip()])
    with (out / 'summary.json').open('w', encoding='utf-8') as f:
        json.dump({'rows': rows, 'provenance': provenance}, f, indent=2)
    with (out / 'summary.csv').open('w', newline='', encoding='utf-8') as f:
        writer = csv.DictWriter(f, fieldnames=['model', 'case', 'chi2', 'mass_hi_ev', 'mass_mid_ev', 'mass_lo_ev', 'delta_m21_sq_ev2', 'delta_m3l_sq_ev2', 'sin2_theta12', 'sin2_theta13', 'sin2_theta23', 'run_id', 'status', 'source_path'])
        writer.writeheader(); writer.writerows(rows)
    print(f'Wrote {out} | {len(rows)} rows')
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
