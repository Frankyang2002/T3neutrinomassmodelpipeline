"""Regression for the direct-Weinberg transport summary schema."""
from pathlib import Path
from types import SimpleNamespace

from RGE.running.backends.ScalarOnlyAfterFermion import _record_full_flavor_success


def test_transport_summary_uses_nested_scales_and_real_term_count(tmp_path: Path):
    record = SimpleNamespace(summary={}, output_dir=tmp_path)
    paths = {
        name: tmp_path / "data" / filename
        for name, filename in {
            "flavor_seed": "flavor_seed.json",
            "flavor_transport": "flavor_transport.json",
            "insertion": "insertion.wl",
            "final_c5": "final_c5.json",
        }.items()
    }
    transport = {
        "status": "Success",
        "scales": {"mu_high": "MF", "mu_low": "MS", "log_ratio": "log(MS/MF)"},
        "running_corrections": {"Weinberg": {"terms": [{"source": "C12"}]}},
    }
    _record_full_flavor_success(
        record, paths,
        flavor_seed={"status": "Success"},
        flavor_transport=transport,
        insertion={"status": "Success"},
        resume={"status": "Success"},
        final_c5={
            "status": "Success",
            "combined": {"ready_for_physical_majorana_numerics": True},
            "direct_running": {"entries": [{}]},
        },
    )
    s = record.summary
    assert s["EFT1FullFlavorBridgeStatus"] == "Success"
    assert s["EFT1WilsonTransportMuHigh"] == "MF"
    assert s["EFT1WilsonTransportMuLow"] == "MS"
    assert s["EFT1WilsonTransportLogRatio"] == "log(MS/MF)"
    assert s["EFT1WilsonTransportCorrectionCount"] == 1
    assert s["FinalWeinbergDirectRunningTermCount"] == 1
