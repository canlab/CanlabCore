"""Build every codemap spec under docs/_codemap_tools/specs and render to png.

Usage:
    PYTHONPATH=. python3 build_all.py [name_pattern1 name_pattern2 ...]

Without arguments, builds every spec. With names (without .py), builds only
the matching specs. Each spec module must expose a build() function returning
a Slide.
"""

from __future__ import annotations

import glob
import importlib.util
import os
import sys
import time

from codemap_lib import render_pptx_to_png


def main() -> int:
    here = os.path.dirname(os.path.abspath(__file__))
    spec_dir = os.path.join(here, "specs")
    pptx_dir = os.path.normpath(os.path.join(here, "..", "code_maps_pptx"))
    png_dir = os.path.normpath(os.path.join(here, "..", "code_maps_png"))

    patterns = sys.argv[1:]
    paths = sorted(glob.glob(os.path.join(spec_dir, "*.py")))
    if patterns:
        paths = [p for p in paths if any(pat in os.path.basename(p) for pat in patterns)]
    if not paths:
        print("No specs matched.", file=sys.stderr)
        return 1

    failures = []
    t0 = time.time()
    for path in paths:
        name = os.path.splitext(os.path.basename(path))[0]
        spec = importlib.util.spec_from_file_location(name, path)
        mod = importlib.util.module_from_spec(spec)
        try:
            spec.loader.exec_module(mod)
            slide = mod.build()
            pptx_path = os.path.join(pptx_dir, f"{name}_codemap.pptx")
            slide.save(pptx_path)
            png_path = render_pptx_to_png(pptx_path, png_dir)
            print(f"[OK] {name} -> {os.path.basename(png_path)}")
        except Exception as exc:  # noqa: BLE001
            failures.append((name, exc))
            print(f"[FAIL] {name}: {exc}", file=sys.stderr)
    dt = time.time() - t0
    print(f"\nBuilt {len(paths) - len(failures)}/{len(paths)} specs in {dt:.1f}s")
    if failures:
        print("Failures:")
        for n, e in failures:
            print(f"  {n}: {e}")
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
