"""Label the exact T1/T2/T3 UV deformations in saved T3 sensitivity HTML.

Project path: tools/LabelSensitivityDirections.py

From the project root:
    python -m tools.LabelSensitivityDirections
    python -m tools.LabelSensitivityDirections --file Reports/output/display/fitted_model_envelopes_all.html

Presentation only: edits generated HTML in place; no numerical calculation.
Definitions match Numerical/fitting/BenchmarkSensitivity.py.
This script can also be invoked by tools/RebuildCurrentReports.py after rebuilding
fitted and common-input envelope dashboards.
"""
from __future__ import annotations

import argparse
from pathlib import Path
import re

ROOT = Path(__file__).resolve().parents[1]
DISPLAY = ROOT / "Reports" / "output" / "display"
TARGET_FILENAMES = (
    "fitted_model_envelopes_all.html",
    "fitted_model_envelopes.html",
    "common_input_model_envelopes_all.html",
)
MARKER = "t3-uv-direction-labels-v1"

# The saved sensitivity definition is in Numerical/fitting/BenchmarkSensitivity.py:
# T1: lambdaT3 -> gamma*lambdaT3; y1,y2 -> gamma**(-1/2)*(y1,y2)
# T2: lambdaT3 -> gamma*lambdaT3; y2 -> gamma**(-1)*y2
# T3: y1 -> gamma*y1; y2 -> gamma**(-1)*y2
DIRECTION_LABEL_SCRIPT = r"""
<!-- t3-uv-direction-labels-v1: labels only; no effect on numeric data -->
<style>
#dirs label.choice { display: flex; align-items: flex-start; gap: 7px; line-height: 1.4; }
#dirs label.choice > input[type="checkbox"] { flex: none; margin-top: 3px; }
#dirs .t3-direction-definition { flex: 1; min-width: 0; overflow-wrap: anywhere; }
#dirs .t3-direction-note { font-size: 11px; line-height: 1.4; margin: 9px 0 12px; opacity: 0.85; }
</style>
<script>
(function () {
  'use strict';
  const formulas = {
    T1: '<strong>T1</strong> — &lambda;<sub>T3</sub> &rarr; &gamma;&lambda;<sub>T3</sub>; '
        + 'y<sub>1</sub>, y<sub>2</sub> &rarr; &gamma;<sup>&minus;1/2</sup>y<sub>1</sub>, '
        + '&gamma;<sup>&minus;1/2</sup>y<sub>2</sub>',
    T2: '<strong>T2</strong> — &lambda;<sub>T3</sub> &rarr; &gamma;&lambda;<sub>T3</sub>; '
        + 'y<sub>2</sub> &rarr; &gamma;<sup>&minus;1</sup>y<sub>2</sub> '
        + '(y<sub>1</sub> unchanged)',
    T3: '<strong>T3</strong> — y<sub>1</sub> &rarr; &gamma;y<sub>1</sub>; '
        + 'y<sub>2</sub> &rarr; &gamma;<sup>&minus;1</sup>y<sub>2</sub> '
        + '(&lambda;<sub>T3</sub> unchanged)'
  };
  const container = document.getElementById('dirs');
  if (!container) return;
  for (const label of container.querySelectorAll('label')) {
    const key = label.textContent.trim();
    if (!Object.prototype.hasOwnProperty.call(formulas, key)) continue;
    const checkbox = label.querySelector('input[type="checkbox"]');
    if (!checkbox) continue;
    const description = document.createElement('span');
    description.className = 't3-direction-definition';
    description.innerHTML = formulas[key]; // Static equations, no user inputs.
    label.replaceChildren(checkbox, description); // Keeps the existing change handler.
  }
  if (!container.querySelector('.t3-direction-note')) {
    const note = document.createElement('p');
    note.className = 't3-direction-note';
    note.innerHTML = '&gamma; = 1 is the original UV input. Each direction preserves '
        + '&lambda;<sub>T3</sub>(y<sub>1</sub>y<sub>2</sub> + y<sub>2</sub>y<sub>1</sub>) '
        + 'at the input scale; running and matching can still change the low-energy result.';
    container.appendChild(note);
  }
})();
</script>
"""


def label_html(path: Path) -> bool:
    """Return True for a new label injection, False if already present."""
    path = Path(path)
    if not path.is_file():
        raise FileNotFoundError(path)
    html = path.read_text(encoding="utf-8-sig")
    if MARKER in html:
        return False
    if not re.search(r'<div\s+id=["\']dirs["\']', html, re.IGNORECASE):
        raise ValueError(f"No T1/T2/T3 direction checkbox container: {path}")
    if not all(re.search(r'\b' + d + r'\b', html) for d in ("T1", "T2", "T3")):
        raise ValueError(f"Expected T1/T2/T3 controls are missing: {path}")
    match = list(re.finditer(r"</body\s*>", html, flags=re.IGNORECASE))
    if len(match) != 1:
        raise ValueError(f"Expected exactly one HTML body closing tag: {path}")
    insertion = match[0].start()
    updated = html[:insertion] + DIRECTION_LABEL_SCRIPT + "\n" + html[insertion:]
    path.write_text(updated, encoding="utf-8")
    return True


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--file", type=Path, action="append",
                        help="Explicit generated HTML file; may be supplied multiple times")
    parser.add_argument("--display-root", type=Path, default=DISPLAY,
                        help="Display root for the retained fitted/common-input sensitivity HTML names")
    args = parser.parse_args(argv)

    paths = (args.file if args.file else
             [args.display_root / name for name in TARGET_FILENAMES])
    changed = 0
    for path in paths:
        path = Path(path)
        if not path.is_file() and not args.file:
            print(f"[SKIP] Missing display: {path}")
            continue
        try:
            if label_html(path):
                changed += 1
                print(f"[LABELED] {path}")
            else:
                print(f"[ALREADY LABELED] {path}")
        except (OSError, ValueError) as exc:
            parser.error(str(exc))
    print(f"Direction-label injection complete: {changed} file(s) updated")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
