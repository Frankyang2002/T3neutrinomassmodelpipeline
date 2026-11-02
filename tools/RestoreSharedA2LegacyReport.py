"""Restore the original T3 RGE dashboard template for the shared-A2 general16 study.

Presentation-only. Loads the *actual* InteractiveModelComparison.py HTML template
from this checkout, never maintains a separate chart implementation. Reads only
saved trajectories and endpoint results. Does not integrate, match, or fit.
"""
from __future__ import annotations

import argparse
import ast
import hashlib
import html
import json
import math
from pathlib import Path
import shutil

ROOT = Path(__file__).resolve().parents[1]
RUN = (ROOT / 'output/studies/general_scalar_extended/shared_a2_general16_reports'
       / 'shared_A2_general16_gamma08_12_v1')
LEGACY = ROOT / 'Numerical/plotting/InteractiveModelComparison.py'
OBSOLETE = ROOT / 'tools/RebuildSharedA2ClassicRGE.py'
SCALES = dict(mu_uv=1e7, mu_fermion_threshold=1e5,
              mu_scalar_threshold=1e3, mu_low=100)


def read(path: Path):
    return json.loads(path.read_text(encoding='utf-8-sig'))


def legacy_template() -> str:
    """Extract literal HTML/JS from the installed original generator verbatim."""
    tree = ast.parse(LEGACY.read_text(encoding='utf-8'))
    for node in tree.body:
        if isinstance(node, ast.FunctionDef) and node.name == 'write_dashboard':
            for item in ast.walk(node):
                if (isinstance(item, ast.Assign) and any(isinstance(t, ast.Name) and
                    t.id == 'template' for t in item.targets) and
                    isinstance(item.value, ast.Constant) and isinstance(item.value.value, str)):
                    template = item.value.value
                    for token in ('__TITLE__', '__OPTIONS__', '__DATA__', '__META__',
                                  'const COLORS=', 'const regionColours='):
                        if token not in template:
                            raise ValueError(f'Original template missing required marker {token}')
                    return template
    raise ValueError('Could not locate original write_dashboard() literal HTML template')


def finite(values, expected):
    if not isinstance(values, list) or len(values) != expected:
        return False
    return all(isinstance(v, (int, float)) and math.isfinite(v) for v in values)


def format_matrix(row):
    if not isinstance(row, dict) or 'real' not in row:
        return None
    re = row['real']
    im = row.get('imag', [[0]*3 for _ in range(3)])
    if len(re)!=3 or len(im)!=3 or any(len(r)!=3 for r in re+im):
        return None
    m = [[complex(re[i][j], im[i][j]) for j in range(3)] for i in range(3)]
    if not all(math.isfinite(z.real) and math.isfinite(z.imag) for row in m for z in row):
        return None
    return math.sqrt(sum(abs(z)**2 for row in m for z in row))


def dataset(folder: Path):
    status = read(folder/'rge_status.json')
    if status.get('no_weinberg_trajectory') is not True:
        raise ValueError('Unexpected new Weinberg data: re-evaluate report assumptions')
    models = list(status['models_available'])
    if len(models) != 14 or len(set(models)) != 14:
        raise ValueError('Unexpected RGE model coverage; inspect saved data')
    quantities = {}
    metadata = {}
    selector_groups = {'UV Yukawa and scalar running':[],
                       'EFT1 scalar running':[]}
    labels = {
        'gY':'g₁ hypercharge coupling', 'g2':'g₂ weak coupling', 'g3':'g₃ strong coupling',
        'lambdaH':'λH Higgs quartic', 'lambdaT3_abs':'|λT3|',
        'lambdaT3_real':'Re(λT3)', 'y1_norm':'‖Y₁‖F', 'y2_norm':'‖Y₂‖F',
        'y1Conj_norm':'‖Y₁conj‖F', 'y2Conj_norm':'‖Y₂conj‖F',
        'mS1_gev':'mS1 [GeV]', 'mS2_gev':'mS2 [GeV]',
        'mS1_sq_gev2':'mS1² [GeV²]', 'mS2_sq_gev2':'mS2² [GeV²]',
        'muS2_real_gev2':'Re(μS2²) [GeV²]', 'muS2_imag_gev2':'Im(μS2²) [GeV²]',
        'muS1_real_gev2':'Re(μS1²) [GeV²]', 'muS1_imag_gev2':'Im(μS1²) [GeV²]',
    }
    endpoints=[]
    for model in models:
        path=folder/model.replace('+','p').replace('-','m')/'rge_trajectories.json'
        # Actual saved directory names: A_m2, B_p1, etc.
        if not path.exists():
            raise FileNotFoundError(path)
        record=read(path)
        if record.get('status')!='Success' or record.get('model')!=model:
            raise ValueError(f'Bad trajectory provenance: {path}')
        for stage in ('uv','eft1'):
            block=record[stage]
            x=block['mu_gev']
            if not finite(x,len(x)) or len(x)<2 or any(v<=0 for v in x):
                raise ValueError(f'Invalid scale samples: {model} {stage}')
            for name,y in block['series'].items():
                if not finite(y,len(x)):
                    continue
                key=f'{stage}_{name}'
                quantities.setdefault(key,{})[model]={'x':x,'y':y}
                if key not in metadata:
                    title=f'{labels.get(name,name)} — {"UV" if stage=="uv" else "EFT1"}'
                    metadata[key]={'label':title,'unit':'','weinberg':False,'block':stage}
                    selector_groups['UV Yukawa and scalar running' if stage=='uv' else 'EFT1 scalar running'].append(key)
    for model in read(folder/'settings.json')['models']:
        centre_file=folder/model.replace('+','p').replace('-','m')/'centre.json'
        if not centre_file.is_file():
            continue
        centre=read(centre_file)
        pred=centre.get('prediction',{})
        num=centre.get('numerical_details',{})
        c5_scalar=format_matrix(num.get('c5_at_scalar_threshold'))
        c5_low=format_matrix(num.get('c5_at_low_scale'))
        endpoints.append({'model':model,'status':centre.get('status','unknown'),
                          'c5_scalar':c5_scalar,'c5_low':c5_low,
                          'prediction':pred})
    metadata['llss_kernel']={'label':'LLSS analytic one-loop kernel η (not CLLSS running)',
                             'unit':'dimensionless','weinberg':False,'block':'analytic_kernel'}
    return ({'models':models,'scales':SCALES,'quantities':quantities,
             'matched':{},'components':{},'overlays':{},
             'continuation_methods':{},'continuation_issues':{}},
            metadata,selector_groups,endpoints,status)


def endpoint_html(endpoints,status):
    rows=[]
    def fmt(v):
        return '—' if v is None or not isinstance(v,(float,int)) or not math.isfinite(v) else f'{v:.5g}'
    for item in endpoints:
        p=item['prediction']
        cells=[item['model'],fmt(item['c5_scalar']),fmt(item['c5_low']),
               fmt(p.get('delta_m21_sq_ev2')),fmt(p.get('delta_m3l_sq_ev2')),
               fmt(p.get('sin2_theta12')),fmt(p.get('sin2_theta13')),
               fmt(p.get('sin2_theta23'))]
        rows.append('<tr>'+''.join(f'<td>{html.escape(c)}</td>' for c in cells)+'</tr>')
    missing=', '.join(html.escape(x['model']) for x in status.get('unavailable',[]))
    return '''<section class="legacy-endpoints" style="max-width:1450px;margin:12px auto 30px;padding:16px 22px;">
<h3 style="margin-top:0">C₅ in EFT1 and SM — saved endpoints only</h3>
<p class="small">The available source contains C₅ at the scalar threshold and low-energy endpoint. There is <strong>no stored continuous C₅ running curve</strong>; the table does not interpolate between endpoints. The EFT1 C₅ trajectory is likewise not saved.</p>
<h3>C<sub>LLSS</sub> in EFT1</h3>
<p class="small">Select “LLSS analytic one-loop kernel η” in the original plotting control. This is the universal η(μ)=ln(μ/MF)/(16π²) kernel, <strong>not</strong> a numerical full-flavour C<sub>LLSS</sub> evolution.</p>
<h3>Neutrino observables and C₅ endpoints</h3>
<div style="overflow-x:auto"><table style="width:100%;border-collapse:collapse;font-size:12px;white-space:nowrap">
<thead><tr><th>Model</th><th>‖C₅(10³ GeV)‖F [GeV⁻¹]</th><th>‖C₅(low)‖F [GeV⁻¹]</th><th>Δm²₂₁ [eV²]</th><th>Δm²₃ℓ [eV²]</th><th>sin²θ₁₂</th><th>sin²θ₁₃</th><th>sin²θ₂₃</th></tr></thead><tbody>''' + ''.join(rows) + '''</tbody></table></div>
<p class="small">UV/EFT1 plot coverage: 14/16 models. No dense UV/EFT1 trajectories for '''+missing+'''. C₅ endpoint predictions are conditional, not certified sequential matching; these are saved forward results, not refits.</p>
</section><style>.legacy-endpoints th,.legacy-endpoints td{border-bottom:1px solid #dce3ee;padding:7px;text-align:right}.legacy-endpoints th:first-child,.legacy-endpoints td:first-child{text-align:left}.legacy-endpoints h3{font-size:16px;margin-top:22px}</style>'''


def build(folder:Path, output:Path, *, remove_obsolete=True):
    template=legacy_template()
    data,meta,groups,endpoints,status=dataset(folder)
    def group(label,keys):
        return '<optgroup label="'+html.escape(label)+'">'+''.join(
            '<option value="'+html.escape(k)+'">'+html.escape(meta[k]['label'])+'</option>'
            for k in keys)+'</optgroup>'
    # The original selector and canvas code remain unchanged: every series is
    # represented as a valid one-stage curve, with no artificial threshold join.
    options=(group('C₅ / SM — endpoints listed below',[])
             +group('C_LLSS in EFT1 — analytic kernel only',['llss_kernel'])
             +''.join(group(label,keys) for label,keys in groups.items())
             +group('Observables — saved endpoints listed below',[]))
    page=(template.replace('__TITLE__','Shared A(−2) general16 — original interactive RGE report')
          .replace('__OPTIONS__',options)
          .replace('__DATA__',json.dumps(data,separators=(',',':'),allow_nan=False))
          .replace('__META__',json.dumps(meta,separators=(',',':'),allow_nan=False)))
    # Append endpoint material outside the original renderer. Original chart
    # styling, controls, colours and JS plotting logic are not copied or edited.
    page=page.replace('</body></html>', '</body></html>',1)
    insert=endpoint_html(endpoints,status)
    if '</body>' not in page:
        raise ValueError('Invalid legacy HTML template')
    page=page.replace('</body>',insert+'</body>',1)
    # Default to real saved UV data rather than the model-independent LLSS kernel.
    first=groups['UV Yukawa and scalar running'][0]
    page=page.replace('resetView();\n</script>',
                      'select.value='+json.dumps(first)+';\nresetView();\n</script>')
    output.parent.mkdir(parents=True,exist_ok=True)
    backup=output.with_name('interactive_RGE.before_legacy_restore.html')
    if output.is_file() and not backup.exists():
        shutil.copy2(output,backup)
    tmp=output.with_suffix('.html.tmp')
    tmp.write_text(page,encoding='utf-8')
    tmp.replace(output)
    if remove_obsolete and OBSOLETE.is_file():
        OBSOLETE.unlink()  # only exact previous generated script
    return {'output':str(output),'backup':str(backup),
            'legacy_source':str(LEGACY),'sha256':hashlib.sha256(page.encode()).hexdigest(),
            'models':len(data['models']),'endpoint_rows':len(endpoints),
            'obsolete_script_deleted':remove_obsolete and not OBSOLETE.exists()}


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--run-dir',type=Path,default=RUN)
    ap.add_argument('--output',type=Path)
    ap.add_argument('--keep-obsolete-script',action='store_true')
    a=ap.parse_args()
    folder=a.run_dir.resolve()
    output=(a.output or folder/'interactive_RGE.html').resolve()
    print(json.dumps(build(folder,output,remove_obsolete=not a.keep_obsolete_script),indent=2))

if __name__=='__main__':
    main()
