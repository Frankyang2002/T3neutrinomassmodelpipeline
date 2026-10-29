"""Presentation-only synthetic fixture tests; these are NOT physics validation."""
import csv
import importlib.util
import json
from pathlib import Path

import pytest

MODULE = Path(__file__).resolve().parents[1] / 'Numerical/plotting/BuildGeneralScalarLambdaReports.py'
spec = importlib.util.spec_from_file_location('gs_lambda_reports', MODULE)
mod = importlib.util.module_from_spec(spec)
spec.loader.exec_module(mod)


def fixture(root: Path):
    run=root/'output/studies/general_scalar_extended/continuous_full_scalar/continuous_v1'
    run.mkdir(parents=True)
    for model in mod.MODELS:
        samples=[]
        for a, b in ((1e7,1e5),(1e5,1e3)):
            samples.append([{'mu_gev':mu,'mS1Sq_gev2':1.25e6,
                             'mS2Sq_gev2':1e6,'extra_quartics':{'testQuartic':{'real':0.1,'imag':0}},
                             'mu2_gev2':{'real':10,'imag':0}}
                            for mu in (a,(a*b)**0.5,b)])
        payload={'format':mod.FORMAT,'status':mod.STATUS,'model':model,
                 'quadratic_kind':'none', 'scalar_endpoint_physically_admissible':True,
                 'scalar_endpoint_spectrum':{'minimum_mass_squared_gev2':10},
                 'settings':{'uv_start_gev':1e7,'fermion_threshold_gev':1e5,'scalar_endpoint_gev':1e3},
                 'uv_trajectory':samples[0],'eft1_trajectory':samples[1]}
        folder=run/mod.slug(model);folder.mkdir()
        (folder/'continuous_trajectory.json').write_text(json.dumps(payload),encoding='utf-8')
    baseline=run/'neutrino_v1';baseline.mkdir()
    with (baseline/'summary.csv').open('w',newline='') as stream:
        writer=csv.DictWriter(stream,fieldnames=['model','status','largest_neutrino_mass_ev',
            'middle_neutrino_mass_ev','smallest_neutrino_mass_ev','sin2_theta12'])
        writer.writeheader()
        for m in mod.MODELS:
            writer.writerow({'model':m,'status':'Fixture','largest_neutrino_mass_ev':0.05,
                'middle_neutrino_mass_ev':0.009,'smallest_neutrino_mass_ev':0.001,
                'sin2_theta12':0.3})
    coupled=run/'bc_coupled_wilson_v1';coupled.mkdir()
    with (coupled/'summary.csv').open('w',newline='') as stream:
        writer=csv.DictWriter(stream,fieldnames=['model','status','new_largest_mnu_ev',
            'new_middle_mnu_ev','new_smallest_mnu_ev']);writer.writeheader()
        writer.writerow({'model':'B_-1','status':'FixtureCoupled','new_largest_mnu_ev':0.048,
            'new_middle_mnu_ev':0.009,'new_smallest_mnu_ev':0.001})
    return run


RGE='''<html><script>const DATA=__DATA__;const META=__META__;</script>
__TITLE__ __OPTIONS__
Numerical integration direction: UV → IR. Saved benchmark trajectories; no matching or RGE recalculation.
<button id="zoomUV">Full T3</button>
<button id="zoomInt">Intermediate EFT</button>
<span><i class="dot" style="background:#eaf2ff"></i>Full T3</span>
<span><i class="dot" style="background:#fff1cf"></i>Intermediate scalar EFT</span>
rectRegion(s.mu_uv,s.mu_fermion_threshold,regionColours.uv,'Full T3');
rectRegion(s.mu_fermion_threshold,s.mu_scalar_threshold,regionColours.int,'Intermediate scalar EFT');
<button id="zoomFinal">SM + Weinberg</button>
<span><i class="dot" style="background:#e5f5e9"></i>SM + Weinberg</span></html>'''
SENS='''<html><div class="layout"><aside></aside></div><script>const D=__DATA__;
let x='sin2_theta12',y='delta_m21_sq_ev2';</script>
Individually fitted T3 model sensitivity
Each model uses its own fitted λT3, Yukawa matrices and other fitted parameters. Data from saved scans only — no RGE recomputation.
<input type="checkbox" id="experimental" checked>NuFIT Gaussian 1σ/3σ
<input type="checkbox" id="future" checked>JUNO 6-year projection (solar only)
Fitted centres</html>'''


def test_full_build_with_lambda_template_fixtures(tmp_path,monkeypatch):
    run=fixture(tmp_path)
    src=tmp_path/'Numerical/plotting';src.mkdir(parents=True)
    (src/'InteractiveModelComparison.py').write_text('def write_dashboard():\n    template = '+repr(RGE)+'\n')
    (src/'BuildFittedModelEnvelopes.py').write_text('HTML = '+repr(SENS)+'\n')
    monkeypatch.setattr(mod,'ROOT',tmp_path)
    out=tmp_path/'report'
    manifest=mod.build(run,out)
    assert manifest['model_count']==16
    assert manifest['preferred_model_count']==16
    assert manifest['variant_count']==17
    assert manifest['sensitivity_model_count']==0
    for name in ('index.html','interactive_rge.html','sensitivity_report.html','all_results.html',
                 'all_results.csv','matching_variants.csv','provenance.json'):
        assert (out/name).is_file(),name
    assert 'No sensitivity ovals have been manufactured' in (out/'sensitivity_report.html').read_text()
    assert 'mS1²' in (out/'interactive_rge.html').read_text()
    assert 'testQuartic' in (out/'interactive_rge.html').read_text()
    with (out/'all_results.csv').open(newline='',encoding='utf-8') as f:
        selected={r['model']:r for r in csv.DictReader(f)}
    assert len(selected)==16
    assert selected['B_-1']['source']=='coupled'
    assert selected['B_-1']['sin2_theta12']==''  # Do not inherit baseline mixing
    assert selected['A_+0']['sin2_theta12']=='0.3'
    assert selected['A_+0']['scalar_endpoint_admissible']=='True'


def test_old_sensitivity_rejected(tmp_path):
    run=fixture(tmp_path)
    old=tmp_path/'output/full/fitted_sensitivity/old.csv';old.parent.mkdir(parents=True)
    old.write_text('model,direction,gamma\nA_+0,T1,1\n')
    with pytest.raises(ValueError,match='historical reduced-potential'):
        mod.read_sensitivity(run,old,{})


def test_saved_gamma_points_are_required_to_render_ovals(tmp_path):
    run=fixture(tmp_path)
    scan=run/'general_scalar_gamma_points.csv'
    with scan.open('w',newline='') as stream:
        wr=csv.DictWriter(stream,fieldnames=('model','direction','gamma','sin2_theta12','delta_m21_sq_ev2'))
        wr.writeheader()
        for d in ('T1','T2','T3'):
            for gamma in (0.8,1.0,1.2):
                # All centres identical; direction-dependent response is intentionally a fixture.
                wr.writerow({'model':'B_-1','direction':d,'gamma':gamma,
                             'sin2_theta12':0.3+(gamma-1)*(0.1 if d=='T1' else 0.2),
                             'delta_m21_sq_ev2':7.5e-5+(gamma-1)*(1e-5 if d=='T3' else 5e-6)})
    models,skipped,state=mod.read_sensitivity(run,scan,{})
    assert state=='SAVED_GENERAL_SCALAR_GAMMA_CSV'
    assert len(models)==1
    assert models[0]['model']=='B_-1'
    rendered=mod.build_sensitivity(models,skipped,state,SENS,scan)
    assert 'Saved general-scalar gamma scans only' in rendered
    assert '"gamma":1.2' in rendered


def test_missing_or_invalid_curve_fails_before_writing(tmp_path,monkeypatch):
    run=fixture(tmp_path)
    p=run/'A_p0/continuous_trajectory.json'
    j=json.loads(p.read_text());j['uv_trajectory'][1]['mu_gev']=1e9
    p.write_text(json.dumps(j))
    monkeypatch.setattr(mod,'ROOT',tmp_path)
    dest=tmp_path/'report'
    with pytest.raises(ValueError,match='Non-descending'):
        mod.build(run,dest)
    assert not dest.exists()
