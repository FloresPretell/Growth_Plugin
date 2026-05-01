#!/usr/bin/env python3
"""
Create a manifest of PNG/u_t directories for SWC export.

Handles both:
  - Old paths: /path/to/Cases/D*_V*/Simulation2D/merged
  - New paths: /path/to/Cases/D*_V*/Simulation2D
"""

import argparse
from pathlib import Path


def parse_args():
    parser = argparse.ArgumentParser(
        description="Create a manifest of PNG/u_t directories for SWC export."
    )
    parser.add_argument("--path", help="Existing PNG root containing per-case folders")
    parser.add_argument("--source-dirs-file", help="Stage 1 manifest (VTKHDF or merged)")
    parser.add_argument("--png-root", help="PNG root to build case paths")
    parser.add_argument("--output", required=True, help="Output manifest file path")
    return parser.parse_args()


def extract_case_code(path_str: str) -> str:
    """
    Extract case code (D*_V*) from either:
      - /path/to/Cases/D2.5_V2p5/Simulation2D/merged
      - /path/to/Cases/D2.5_V2p5/Simulation2D
    
    Returns: "D2.5_V2p5"
    """
    path = Path(path_str.strip()).expanduser().resolve()
    
    # Handle both patterns
    if path.name == "merged":
        # Old pattern: .../Simulation2D/merged
        case_code = path.parent.parent.name
    elif path.name == "Simulation2D":
        # New pattern: .../Simulation2D
        case_code = path.parent.name
    else:
        raise ValueError(
            f"Expected path ending in 'Simulation2D' or 'merged', got: {path}"
        )
    
    return case_code


def build_from_manifest(source_file: Path, png_root: Path):
    """Build u_t paths from Stage 1 manifest"""
    entries = []
    with source_file.open("r", encoding="utf-8") as handle:
        for line in handle:
            stripped = line.strip()
            if not stripped:
                continue
            try:
                case_code = extract_case_code(stripped)
                u_t_path = (png_root / case_code / "u_t").resolve()
                entries.append(u_t_path)
            except ValueError as e:
                print(f"[warning] Skipping line: {e}")
    return entries


def build_from_png_root(png_root: Path):
    """Build u_t paths by scanning PNG root"""
    entries = []
    for case_dir in sorted(png_root.iterdir()):
        if not case_dir.is_dir():
            continue
        u_t_dir = (case_dir / "u_t").resolve()
        if u_t_dir.is_dir():
            entries.append(u_t_dir)
    return entries


def main():
    args = parse_args()
    output_path = Path(args.output).expanduser().resolve()
    output_path.parent.mkdir(parents=True, exist_ok=True)

    if args.source_dirs_file:
        source_file = Path(args.source_dirs_file).expanduser().resolve()
        png_root = Path(args.png_root).expanduser().resolve()
        
        if not source_file.is_file():
            raise SystemExit(f"ERROR: Source manifest not found: {source_file}")
        
        print(f"[info] Reading manifest: {source_file}")
        print(f"[info] PNG root: {png_root}")
        
        entries = build_from_manifest(source_file, png_root)
        
    elif args.path:
        png_root = Path(args.path).expanduser().resolve()
        
        if not png_root.is_dir():
            raise SystemExit(f"ERROR: PNG root not found: {png_root}")
        
        print(f"[info] Scanning PNG root: {png_root}")
        entries = build_from_png_root(png_root)
        
    else:
        raise SystemExit(
            "ERROR: Provide either --path or --source-dirs-file with --png-root"
        )

    if not entries:
        raise SystemExit("ERROR: No u_t directories found")

    with output_path.open("w", encoding="utf-8") as handle:
        for entry in entries:
            handle.write(f"{entry}\n")

    print(f"\n[success] Wrote {len(entries)} u_t paths to {output_path}")


if __name__ == "__main__":
    main()