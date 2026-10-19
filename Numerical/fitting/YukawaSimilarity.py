"""Compare frozen complex Yukawa pairs and form complete-linkage clusters.

This is an entrywise, basis-dependent comparison. Cluster membership is not a
claim of physically equivalent flavour models or acceptable cross-model fits.
"""
from __future__ import annotations

import hashlib
import json
import math
from copy import deepcopy
from pathlib import Path
from typing import Any

import numpy as np


def read_object(path: Path) -> dict[str, Any]:
    value = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(value, dict):
        raise ValueError(f"Expected JSON object: {path}")
    return value


def yukawa_pair(payload: dict[str, Any]) -> tuple[np.ndarray, np.ndarray]:
    ordinary = payload["base_state"]["ordinary"]
    matrices: list[np.ndarray] = []
    for name in ("y1", "y2"):
        raw = ordinary[name]
        real = np.asarray(raw["real"], dtype=float)
        imag = np.asarray(raw["imag"], dtype=float)
        if real.shape != (3, 3) or imag.shape != (3, 3):
            raise ValueError(f"{name} must have 3x3 real and imaginary components")
        matrix = real + 1j * imag
        if not np.all(np.isfinite(matrix)):
            raise ValueError(f"Nonfinite values in {name}")
        matrices.append(matrix)
    return matrices[0], matrices[1]


def pair_distance(left: tuple[np.ndarray, np.ndarray], right: tuple[np.ndarray, np.ndarray]) -> float:
    """Dimensionless squared-norm-normalised Frobenius pair distance."""
    numerator = sum(float(np.linalg.norm(a - b, "fro") ** 2) for a, b in zip(left, right))
    denominator = sum(float(np.linalg.norm(x, "fro") ** 2) for pair in (left, right) for x in pair)
    if denominator == 0:
        return 0.0
    return math.sqrt(numerator / denominator)


def load_collection(manifest_path: Path, root: Path) -> list[dict[str, Any]]:
    """Load manifest-listed original configs, validating hashes and identities."""
    manifest = read_object(manifest_path)
    entries = manifest.get("models")
    if not isinstance(entries, list) or not entries:
        raise ValueError("Manifest must contain a nonempty models list")
    collection = []
    seen: set[str] = set()
    for entry in entries:
        name = str(entry["model"])
        if name in seen:
            raise ValueError(f"Duplicate model: {name}")
        seen.add(name)
        relative = Path(str(entry["config_path"]))
        if relative.is_absolute() or ".." in relative.parts:
            raise ValueError(f"Unsafe manifest path: {relative}")
        path = root / relative
        digest = hashlib.sha256(path.read_bytes()).hexdigest()
        if entry.get("sha256") and entry["sha256"] != digest:
            raise ValueError(f"SHA-256 mismatch: {path}")
        payload = read_object(path)
        if payload.get("representation") != entry["representation"]:
            raise ValueError(f"Representation mismatch: {name}")
        collection.append({"model": name, "path": path, "sha256": digest,
                           "payload": payload, "pair": yukawa_pair(payload)})
    return collection


def group_collection(collection: list[dict[str, Any]], threshold: float = 0.1) -> dict[str, Any]:
    """Complete-linkage agglomerative clustering with deterministic tie-breaking."""
    if not math.isfinite(threshold) or threshold < 0:
        raise ValueError("Threshold must be finite and nonnegative")
    items = sorted(collection, key=lambda x: x["model"])
    n = len(items)
    distances = [[pair_distance(a["pair"], b["pair"]) for b in items] for a in items]
    groups = [[i] for i in range(n)]
    while len(groups) > 1:
        choices = []
        for i in range(len(groups)):
            for j in range(i + 1, len(groups)):
                maximum = max(distances[a][b] for a in groups[i] for b in groups[j])
                choices.append((maximum, i, j))
        value, i, j = min(choices)
        if value > threshold:
            break
        groups[i] = sorted(groups[i] + groups[j])
        del groups[j]
    results = []
    for group in groups:
        ranked = sorted(group, key=lambda i: (sum(distances[i][j] for j in group), items[i]["model"]))
        medoid = ranked[0]
        results.append({"size": len(group), "reference_model": items[medoid]["model"],
                        "members": [items[i]["model"] for i in group],
                        "maximum_pairwise_distance": max(distances[a][b] for a in group for b in group),
                        "reference_distances": {items[i]["model"]: distances[medoid][i] for i in group}})
    results.sort(key=lambda g: (-g["size"], g["reference_model"]))
    return {"method": "complete_linkage", "distance": "normalised_complex_pair_frobenius",
            "threshold": threshold, "model_count": n, "groups": results,
            "models": [x["model"] for x in items], "distance_matrix": distances,
            "warning": "Entrywise basis-dependent Yukawa similarity; no cross-model chi2 validation performed."}


def make_candidate(member: dict[str, Any], reference: dict[str, Any]) -> dict[str, Any]:
    """Copy only y1/y2; preserve target model representation and remaining couplings."""
    payload = deepcopy(member["payload"])
    for name in ("y1", "y2"):
        payload["base_state"]["ordinary"][name] = deepcopy(reference["payload"]["base_state"]["ordinary"][name])
    payload["benchmark_search"] = {"enabled": False}
    payload["sensitivity"] = {"enabled": False}
    payload.pop("comparison_point", None)  # Legacy metadata may describe obsolete fixed textures.
    payload["yukawa_grouping"] = {
        "status": "UNVALIDATED_REFERENCE_REPLACEMENT",
        "reference_model": reference["model"],
        "reference_sha256": reference["sha256"],
        "target_model": member["model"],
        "target_original_sha256": member["sha256"],
        "changed_parameters": ["base_state.ordinary.y1", "base_state.ordinary.y2"],
    }
    return payload
