"""Read-only post-processing tests. No Wolfram, Matchete, or RGEs required."""
import csv
import json
import math
from pathlib import Path

import numpy as np
import pytest

from scripts import PlotT3GeneralZeeOvals as plot


def _fake_point(model, direction, gamma):
    k = {"T1": 1., "T2": 3., "T3": -2.}[direction]
    return {
        "model": model, "direction": direction, "gamma": gamma, "status": "Success",
        "chi2": 0.01 + 2 * (gamma - 1)**2,
        "delta_m21_sq_ev2": 7.537e-5 + k * (gamma-1) * 1e-5,
        "delta_m3l_sq_ev2": 2.511e-3 + k*(gamma-1)*1e-4,
        "sin2_theta12": 0.3088 + k*(gamma-1)*0.1,
        "sin2_theta13": 0.02248 + k*(gamma-1)*0.01,
        "sin2_theta23": 0.47 + k*(gamma-1)*0.05,
        "error": "",
    }


def _fixture(tmp_path, *, status="Complete", bad_direction=False):
    src = tmp_path / "source"
    src.mkdir()
    rows = [_fake_point("A_-2", d, g) for d in plot.DIRECTIONS
            for g in (0.99, 1., 1.01)]
    if bad_direction:
        rows.pop()
    with (src / "points.csv").open("w", newline="", encoding="utf-8") as file:
        writer = csv.DictWriter(file, fieldnames=rows[0].keys())
        writer.writeheader()
        writer.writerows(rows)
    (src / "summary.json").write_text(json.dumps({
        "schema": "t3-general-zee-compensated-sensitivity-v1",
        "models": [{"model": "A_-2", "status": status, "failed_points": 0,
                    "baseline_chi2_saved": 0.01}],
        "gamma_grid": [0.99, 1., 1.01]
    }), encoding="utf-8")
    nu = tmp_path / "nufit.json"
    nu.write_text(json.dumps({
        "ordering": "NO", "dataset": "NuFIT 6.1 (2025), IC24 NO",
        "central_values": {k: _fake_point("A_-2", "T1", 1.)[k] for k in plot.OBS},
        "one_sigma_errors": {
            "delta_m21_sq_ev2": 9.7e-7, "delta_m3l_sq_ev2": 2.05e-5,
            "sin2_theta12": 0.00665, "sin2_theta13": 0.00057,
            "sin2_theta23": 0.0155,
        }
    }), encoding="utf-8")
    return src, nu


def test_juno_latest_measured_not_projection():
    assert plot.JUNO_2026["central_values"]["sin2_theta12"] == 0.3036
    assert plot.JUNO_2026["one_sigma_errors"]["sin2_theta12"] == 0.0064
    assert plot.JUNO_2026["central_values"]["delta_m21_sq_ev2"] == 7.388e-5
    assert plot.JUNO_2026["one_sigma_errors"]["delta_m21_sq_ev2"] == 7.8e-7
    assert "2609.37895" in plot.JUNO_2026["url"]


def test_joint_gaussian_radii():
    assert math.isclose(plot.JOINT_RADII["68% joint"]**2, 2.30, rel_tol=.003)
    assert math.isclose(plot.JOINT_RADII["99.73% joint"]**2, 11.83, rel_tol=.003)


def test_read_complete_scan(tmp_path):
    src, _ = _fixture(tmp_path)
    data, meta = plot.load_points(src / "points.csv", src / "summary.json", models=["A_-2"])
    assert len(data["A_-2"]["T1"]) == 3
    assert plot.centre(data["A_-2"])["gamma"] == 1.
    assert meta["A_-2"]["baseline_chi2_saved"] == 0.01


def test_centered_slope(tmp_path):
    src, _ = _fixture(tmp_path)
    data, _ = plot.load_points(src / "points.csv", src / "summary.json", models=["A_-2"])
    deriv, delta_g = plot.central_slope(data["A_-2"]["T2"], "sin2_theta12")
    expected = 0.3*(0.02)/(math.log(1.01) - math.log(0.99))
    assert deriv == pytest.approx(expected, rel=1e-12)
    assert delta_g == pytest.approx(0.02)


def test_model_covariance_is_about_unmoved_centre(tmp_path):
    src, _ = _fixture(tmp_path)
    data, _ = plot.load_points(src / "points.csv", src / "summary.json", models=["A_-2"])
    cov = plot.centered_covariance(data["A_-2"], "sin2_theta12", "delta_m21_sq_ev2")
    assert cov.shape == (2,2)
    assert np.allclose(cov, cov.T)
    assert np.linalg.eigvalsh(cov)[-1] > 0
    assert np.linalg.eigvalsh(cov)[0] >= -1e-17


def test_incomplete_status_rejected(tmp_path):
    src, _ = _fixture(tmp_path, status="PartialSuccess")
    with pytest.raises(ValueError, match="incomplete"):
        plot.load_points(src / "points.csv", src / "summary.json", models=["A_-2"])


def test_bad_direction_rejected(tmp_path):
    src, _ = _fixture(tmp_path, bad_direction=True)
    with pytest.raises(ValueError, match="Bad gamma grid|Inconsistent gamma grids"):
        plot.load_points(src / "points.csv", src / "summary.json", models=["A_-2"])


def test_full_end_to_end_one_model(tmp_path):
    src, nu = _fixture(tmp_path)
    out = tmp_path / "out"
    assert plot.main(["--source", str(src), "--nufit", str(nu),
                      "--output-dir", str(out), "--models", "A_-2",
                      "--no-per-model"]) == 0
    assert (out / "model_grid_zee_zoom_out.png").stat().st_size > 1000
    assert (out / "model_grid_zee_zoom_out.pdf").stat().st_size > 1000
    assert (out / "solar_sensitivity_ranking.csv").is_file()
    assert (out / "baseline_nufit_juno_comparison.csv").is_file()
    prov = json.loads((out / "provenance.json").read_text())
    assert prov["plots_use_saved_points_only"] is True
    assert prov["models"] == ["A_-2"]
    assert "NOT" in prov["physics_scope"]


def test_per_model_files_and_no_grid(tmp_path):
    src, nu = _fixture(tmp_path)
    out = tmp_path / "pictures"
    assert plot.main(["--source", str(src), "--nufit", str(nu),
                      "--output-dir", str(out), "--models", "A_-2",
                      "--no-grid"]) == 0
    assert (out / "models" / "A_m2_zee_zoom_out.png").stat().st_size > 1000
    assert (out / "models" / "A_m2_zee_zoom_out.pdf").stat().st_size > 1000
    assert not (out / "model_grid_zee_zoom_out.png").exists()


def test_save_figure_avoids_pillow_path_open(tmp_path,monkeypatch):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.figure import Figure
    fig, ax = plt.subplots()
    ax.plot([0,1],[2,3])
    real = Figure.savefig
    seen = []
    def wrapper(self, target, *args, **kwargs):
        assert hasattr(target, "write"), "Savefig must target in-memory buffer"
        seen.append(kwargs.get("format"))
        return real(self,target,*args,**kwargs)
    monkeypatch.setattr(Figure,"savefig",wrapper)
    plot._save_figure(fig, tmp_path / "output_image.png",dpi=120)
    plt.close(fig)
    assert seen == ["png", "pdf"]
    assert (tmp_path / "output_image.png").is_file()
    assert (tmp_path / "output_image.pdf").is_file()


def test_model_sensitivity_survives_zee_style(tmp_path):
    src, nu_path = _fixture(tmp_path)
    data, _ = plot.load_points(src / "points.csv",src / "summary.json",models=["A_-2"])
    nu = plot.load_nufit(nu_path)
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots()
    handles=plot._plot_panel(ax,"A_-2",data["A_-2"],nu,
                             "sin2_theta12","delta_m21_sq_ev2")
    labels=[h.get_label() for h in handles]
    assert "NuFIT 6.1 ~1σ (joint Gaussian)" in labels
    assert "JUNO 2026 1σ (measured) (joint Gaussian)" in labels
    assert "T1" in labels and "T3" in labels
    assert "Model response oval (descriptive)" in labels
    plt.close(fig)


def test_figure9_palette_and_current_juno_mode():
    assert plot.PALETTE["nufit_3"] != plot.PALETTE["nufit_1"]
    assert plot.PALETTE["juno_1"] not in {plot.PALETTE["nufit_3"],plot.PALETTE["nufit_1"]}
    assert [r[1] for r in plot.juno_regions("measured")] == ["JUNO 2026 1σ (measured)"]
    assert [r[1] for r in plot.juno_regions("future")] == ["JUNO 6-y 1σ (projection)"]
    assert len(plot.juno_regions("both")) == 2
    with pytest.raises(ValueError,match="Unsupported JUNO mode"):
        plot.juno_regions("unrecognised")


def test_future_is_not_mislabelled_as_measured(tmp_path):
    src,nu_path=_fixture(tmp_path)
    data,_=plot.load_points(src/"points.csv",src/"summary.json",models=["A_-2"])
    nufit=plot.load_nufit(nu_path)
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig,ax=plt.subplots()
    labels=[h.get_label() for h in plot._plot_panel(ax,"A_-2",data["A_-2"],nufit,
                              "sin2_theta12","delta_m21_sq_ev2",juno_mode="future")]
    assert "JUNO 6-y 1σ (projection) (joint Gaussian)" in labels
    assert "JUNO 2026 1σ (measured) (joint Gaussian)" not in labels
    plt.close(fig)


def test_detail_plot_generated_and_uses_same_saved_points(tmp_path):
    src,nu=_fixture(tmp_path)
    out=tmp_path/"figures"
    assert plot.main(["--source",str(src),"--nufit",str(nu),
                      "--output-dir",str(out),"--models","A_-2","--no-grid"])==0
    assert (out/"models"/"A_m2_zee_zoom_in.pdf").stat().st_size>1000
    assert (out/"models"/"A_m2_zee_zoom_in.png").stat().st_size>1000
    provenance=json.loads((out/"provenance.json").read_text())
    assert provenance["juno_plot_mode"]=="measured"
    assert provenance["individual_close_up_drawn"] is True


def test_experimental_and_model_ovals_not_strips(tmp_path):
    src,nu=_fixture(tmp_path)
    data,_=plot.load_points(src/"points.csv",src/"summary.json",models=["A_-2"])
    nufit=plot.load_nufit(nu)
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig,ax=plt.subplots()
    plot._plot_panel(ax,"A_-2",data["A_-2"],nufit,
                     "sin2_theta12","delta_m21_sq_ev2")
    from matplotlib.patches import Ellipse
    # NuFIT 3sigma + 1sigma and JUNO measured 1sigma, plus model response.
    assert len(ax.patches) == 4
    assert all(isinstance(p, Ellipse) for p in ax.patches)
    assert not any(line.get_marker()=="*" for line in ax.lines)
    assert ax.get_xlim()[0]<nufit["central_values"]["sin2_theta12"]<ax.get_xlim()[1]
    plt.close(fig)


def test_two_4x4_grid_views_created(tmp_path):
    src,nu = _fixture(tmp_path)
    out = tmp_path/"dual"
    assert plot.main(["--source",str(src),"--nufit",str(nu),
                      "--output-dir",str(out),"--models","A_-2",
                      "--no-per-model"]) == 0
    for view in ("zoom_out","zoom_in"):
        assert (out/f"model_grid_zee_{view}.png").stat().st_size>1000
        assert (out/f"model_grid_zee_{view}.pdf").stat().st_size>1000
    prov=json.loads((out/"provenance.json").read_text())
    assert prov["baseline_star_drawn"] is False
    assert prov["grid_views"] == ["zoom_out","zoom_in"]
    assert prov["model_response_oval_drawn"] is True
    assert "independent-error Gaussian" in prov["experimental_oval_definition"]


def test_zoom_frames_all_gamma_points_and_response_oval(tmp_path):
    src,nu=_fixture(tmp_path)
    data,_=plot.load_points(src/"points.csv",src/"summary.json",models=["A_-2"])
    nufit=plot.load_nufit(nu)
    model=data["A_-2"]
    xlim,ylim=plot._axis_limits(model,nufit,"measured","sin2_theta12",
                                 "delta_m21_sq_ev2",zoom="trajectory")
    for d in plot.DIRECTIONS:
        for pt in model[d]:
            assert xlim[0] < pt["sin2_theta12"] < xlim[1]
            assert ylim[0] < pt["delta_m21_sq_ev2"]*1e5 < ylim[1]
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.patches import Ellipse
    fig,ax=plt.subplots()
    handles=plot._plot_panel(ax,"A_-2",model,nufit,"sin2_theta12",
                              "delta_m21_sq_ev2",zoom="trajectory")
    assert len(ax.patches)==1  # Only descriptive model oval on zoom.
    assert isinstance(ax.patches[0],Ellipse)
    assert "Model response oval (descriptive)" in [h.get_label() for h in handles]
    assert not any(line.get_marker()=="*" for line in ax.lines)
    plt.close(fig)


def test_no_model_oval_switch_keeps_experimental_ellipses(tmp_path):
    src,nu=_fixture(tmp_path)
    data,_=plot.load_points(src/"points.csv",src/"summary.json",models=["A_-2"])
    nufit=plot.load_nufit(nu)
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig,ax=plt.subplots()
    plot._plot_panel(ax,"A_-2",data["A_-2"],nufit,
                     "sin2_theta12","delta_m21_sq_ev2",show_model_oval=False)
    assert len(ax.patches)==3
    plt.close(fig)


def test_zoom_labels_are_baseline_subtracted_without_scientific_offset(tmp_path):
    src, nu_path = _fixture(tmp_path)
    data, _ = plot.load_points(src / "points.csv", src / "summary.json",
                                models=["A_-2"])
    nu = plot.load_nufit(nu_path)
    model = data["A_-2"]
    frame = plot._detail_frame(model,nu,"sin2_theta12","delta_m21_sq_ev2")
    assert frame["powers"]["sin2_theta12"] < 0
    assert frame["powers"]["delta_m21_sq_ev2"] < 0
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig,ax=plt.subplots()
    plot._plot_panel(ax,"A_-2",model,nu,"sin2_theta12","delta_m21_sq_ev2",
                     zoom="trajectory")
    fig.canvas.draw()
    assert ax.get_title() == "A(-2)"   # no 'zoom to trajectories' suffix
    assert r"\delta(\sin^2\theta_{12})" in ax.get_xlabel()
    assert r"\delta(\Delta m^2_{21})" in ax.get_ylabel()
    assert r"\mathrm{eV}^2" in ax.get_ylabel()
    assert "1e-" not in ax.xaxis.get_offset_text().get_text()
    assert "1e-" not in ax.yaxis.get_offset_text().get_text()
    assert ax.get_xlim()[0] < 0 < ax.get_xlim()[1]
    assert ax.get_ylim()[0] < 0 < ax.get_ylim()[1]
    from matplotlib.patches import FancyArrowPatch
    assert not any(isinstance(p,FancyArrowPatch) for p in ax.patches)
    assert len(ax.lines) == 3
    assert all(line.get_markersize() >= 4.8 for line in ax.lines)
    for d,line in zip(plot.DIRECTIONS,ax.lines):
        expected_x=[(p["sin2_theta12"]-frame["origin"]["sin2_theta12"])
                    /frame["units"]["sin2_theta12"] for p in model[d]]
        expected_y=[(p["delta_m21_sq_ev2"]-frame["origin"]["delta_m21_sq_ev2"])
                    /frame["units"]["delta_m21_sq_ev2"] for p in model[d]]
        assert np.allclose(line.get_xdata(), expected_x)
        assert np.allclose(line.get_ydata(), expected_y)
    plt.close(fig)


def test_closeup_independent_powers_preserve_very_small_response(tmp_path):
    src,nu_path=_fixture(tmp_path)
    data,_=plot.load_points(src/"points.csv",src/"summary.json",models=["A_-2"])
    nu=plot.load_nufit(nu_path)
    d=data["A_-2"]
    # Almost-flat responses like A(-2), D(-2), E(-2) in the real study.
    for direction in plot.DIRECTIONS:
        for row in d[direction]:
            slope={"T1": 1., "T2": -1.5, "T3": .3}[direction]
            row["sin2_theta12"]=.3088 + (row["gamma"]-1)*slope*1e-7
            row["delta_m21_sq_ev2"]=7.537e-5+(row["gamma"]-1)*slope*1e-10
    frame=plot._detail_frame(d,nu,"sin2_theta12","delta_m21_sq_ev2")
    assert frame["powers"]["sin2_theta12"] <= -9
    assert frame["powers"]["delta_m21_sq_ev2"] <= -12
    for direction in plot.DIRECTIONS:
        for row in d[direction]:
            for obs in ("sin2_theta12","delta_m21_sq_ev2"):
                val=(row[obs]-frame["origin"][obs])/frame["units"][obs]
                a,b=frame["limits"][obs]
                assert a<val<b


def test_grid_titles_no_bottom_footer_and_readable_detail(monkeypatch,tmp_path):
    src,nu=_fixture(tmp_path)
    rows,_=plot.load_points(src/"points.csv",src/"summary.json",models=["A_-2"])
    ref=plot.load_nufit(nu)
    from matplotlib.figure import Figure
    saved=[]
    def fake_save(fig,path,*,dpi):
        if str(path).endswith("model_grid_zee_zoom_in.png") or str(path).endswith("model_grid_zee_zoom_out.png"):
            saved.append((str(path),[text.get_text() for text in fig.texts]))
    monkeypatch.setattr(plot,"_save_figure",fake_save)
    plot.render(rows,ref,tmp_path/"view","sin2_theta12","delta_m21_sq_ev2",
                per_model=False)
    assert len(saved)==2
    assert all(len(texts)==2 for _,texts in saved)
    assert any("Detail: compensated-coupling" in texts[0] for _,texts in saved)
    assert not any("Zoom in: T1/T2/T3" in t for _,ts in saved for t in ts)
