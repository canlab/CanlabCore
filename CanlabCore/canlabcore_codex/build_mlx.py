#!/usr/bin/env python3
"""Build a MATLAB .mlx by replacing document.xml and output.xml in a template .mlx.

Usage:
  python3 canlabcore_codex/build_mlx.py \
    --template "HRF_Est_Toolbox2/EstHRF_inAtlas/EstHRF_inAtlas Tutorial.mlx" \
    --document "canlabcore_codex/document.xml" \
    --output "canlabcore_codex/output.xml" \
    --out "canlabcore_codex/My_Live_Script.mlx"
"""

import argparse
import zipfile
from pathlib import Path


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--template", required=True, help="Path to template .mlx")
    parser.add_argument("--document", required=True, help="Path to document.xml")
    parser.add_argument("--output", required=True, help="Path to output.xml")
    parser.add_argument("--out", required=True, help="Path to output .mlx (no spaces)")
    return parser.parse_args()


def main():
    args = parse_args()

    template_path = Path(args.template)
    document_path = Path(args.document)
    output_path = Path(args.output)
    out_path = Path(args.out)

    if " " in out_path.name:
        raise SystemExit("Output .mlx filename must not contain spaces.")

    if not template_path.exists():
        raise SystemExit(f"Template not found: {template_path}")

    if not document_path.exists():
        raise SystemExit(f"document.xml not found: {document_path}")

    if not output_path.exists():
        raise SystemExit(f"output.xml not found: {output_path}")

    document_xml = document_path.read_bytes()
    output_xml = output_path.read_bytes()

    with zipfile.ZipFile(template_path, "r") as zin:
        with zipfile.ZipFile(out_path, "w", compression=zipfile.ZIP_DEFLATED) as zout:
            for item in zin.infolist():
                if item.filename == "matlab/document.xml":
                    zout.writestr(item, document_xml)
                elif item.filename == "matlab/output.xml":
                    zout.writestr(item, output_xml)
                else:
                    zout.writestr(item, zin.read(item.filename))

    print(f"Wrote {out_path}")


if __name__ == "__main__":
    main()
