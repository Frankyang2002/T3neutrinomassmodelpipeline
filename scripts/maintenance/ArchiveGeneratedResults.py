"""Archive generated results outside the repository before a fresh rebuild.

Preview is default. --apply moves output/ and Reports/output/ into a dated
external archive, creates empty replacement directories, and records a manifest.
--restore ARCHIVE_DIR reverses the operation but refuses to overwrite new output.
No source, configs, generated beta files outside those directories, or GitHub
objects are touched.
"""
from __future__ import annotations
import argparse
from datetime import datetime, timezone
import json
import os
from pathlib import Path
import shutil
import sys

TARGETS = ("output", "Reports/output")
MARKER = "t3-generated-output-archive-v1"

def inventory(folder: Path):
    count = 0
    size = 0
    for parent, _, files in os.walk(folder):
        for name in files:
            q = Path(parent) / name
            if q.is_symlink():
                raise RuntimeError(f"Refusing symlink in generated outputs: {q}")
            st = q.stat()
            count += 1
            size += st.st_size
    return {"files": count, "bytes": size}

def within(child: Path, parent: Path) -> bool:
    return child == parent or parent in child.parents

def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--root", type=Path, default=Path(__file__).resolve().parents[2])
    ap.add_argument("--archive-root", type=Path)
    ap.add_argument("--apply", action="store_true")
    ap.add_argument("--restore", type=Path, metavar="ARCHIVE_DIRECTORY")
    args = ap.parse_args()
    root = args.root.resolve()
    if not (root / "pipeline.py").is_file() or not (root / "Numerical").is_dir():
        ap.error("Not the expected T3 repository root")
    base = (args.archive_root or root.parent / "T3_repository_archives").resolve()
    if within(base, root) or within(root, base):
        ap.error("Archive root must be outside the repository")
    if args.restore:
        arc = args.restore.resolve()
        if not within(arc, base):
            ap.error("Restore directory must be inside the configured external archive root")
        manifest_path = arc / "manifest.json"
        manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
        if manifest.get("format") != MARKER or manifest.get("repo_root") != str(root):
            raise RuntimeError("Archive does not belong to this repository")
        restored = []
        for key in TARGETS:
            entry = manifest["entries"].get(key)
            if not entry or not entry["present"]:
                continue
            src, dst = arc / "files" / key, root / key
            if not src.is_dir():
                raise RuntimeError(f"Archived directory missing: {src}")
            if dst.exists() and any(dst.iterdir()):
                raise RuntimeError(f"Refusing overwrite: {dst} is not empty. Move new outputs aside first.")
            if dst.exists():
                dst.rmdir()
            dst.parent.mkdir(parents=True, exist_ok=True)
            shutil.move(str(src), str(dst))
            restored.append(key)
            print("RESTORED", key, flush=True)
        print("RESTORE_COMPLETE", restored)
        return

    entries = {}
    for key in TARGETS:
        folder = root / key
        if folder.is_symlink():
            raise RuntimeError(f"Refusing to archive symlink: {folder}")
        entries[key] = {"present": folder.is_dir(),
                        **(inventory(folder) if folder.is_dir() else {"files": 0, "bytes": 0})}
    preview = {"format": MARKER, "repo_root": str(root), "entries": entries,
               "note": "Only output/ and Reports/output/ are moved. Source and configs remain unchanged."}
    preview_path = root / "generated_archive_preview.json"
    preview_path.write_text(json.dumps(preview, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(entries, indent=2), flush=True)
    if not args.apply:
        print("PREVIEW_ONLY:", preview_path)
        return
    stamp = datetime.now(timezone.utc).strftime("%Y%m%d_%H%M%S_UTC")
    arc = base / ("generated_outputs_" + stamp)
    if arc.exists():
        raise FileExistsError(arc)
    (arc / "files").mkdir(parents=True)
    manifest = dict(preview, created_utc=stamp, status="in_progress", moved=[])
    mp = arc / "manifest.json"
    mp.write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    for key in TARGETS:
        if not entries[key]["present"]:
            continue
        src, dst = root / key, arc / "files" / key
        dst.parent.mkdir(parents=True, exist_ok=True)
        # On same drive this is normally a fast directory rename.
        shutil.move(str(src), str(dst))
        src.mkdir(parents=True, exist_ok=False)
        manifest["moved"].append(key)
        mp.write_text(json.dumps(manifest, indent=2), encoding="utf-8")
        print("ARCHIVED:", key, "->", dst, flush=True)
    manifest["status"] = "complete"
    mp.write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    print("ARCHIVE_COMPLETE:", arc)
    print("New generated outputs may now be created in the empty source directories.")
    print("Old generated matching/RGE inputs will NOT be available until restored or regenerated.")

if __name__ == "__main__":
    try:
        main()
    except Exception as e:
        print("ARCHIVE_ERROR:", str(e), file=sys.stderr)
        raise
