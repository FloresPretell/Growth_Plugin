#!/usr/bin/env python3

import argparse
from pathlib import Path


def parse_args():
    parser = argparse.ArgumentParser(
        description="Create a manifest of PNG/u_t directories for SWC export."
    )
    parser.add_argument("--path", help="Existing PNG root containing per-case folders")
    parser.add_argument("--source-dirs-file", help="Stage 1 manifest used to derive case order")
    parser.add_argument("--png-root", help="Run-specific PNG root used with --source-dirs-file")
    parser.add_argument("--output", required=True, help="Output manifest file path")
    return parser.parse_args()


def case_code_from_stage1_path(path_str: str) -> str:
    merged_dir = Path(path_str.strip()).expanduser().resolve()
    if not merged_dir.name == "merged":
        raise ValueError(f"Expected a merged directory path, got: {merged_dir}")
    return merged_dir.parent.parent.name


def build_from_manifest(source_file: Path, png_root: Path):
    entries = []
    with source_file.open("r", encoding="utf-8") as handle:
        for line in handle:
            stripped = line.strip()
            if not stripped:
                continue
            case_code = case_code_from_stage1_path(stripped)
            entries.append((png_root / case_code / "u_t").resolve())
    return entries


def build_from_png_root(png_root: Path):
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
            raise SystemExit(f"ERROR: Stage 1 manifest does not exist: {source_file}")
        entries = build_from_manifest(source_file, png_root)
    elif args.path:
        png_root = Path(args.path).expanduser().resolve()
        if not png_root.is_dir():
            raise SystemExit(f"ERROR: PNG root does not exist: {png_root}")
        entries = build_from_png_root(png_root)
    else:
        raise SystemExit("ERROR: Provide either --path or --source-dirs-file with --png-root")

    with output_path.open("w", encoding="utf-8") as handle:
        for entry in entries:
            handle.write(f"{entry}\n")

    print(f"Wrote {len(entries)} paths to {output_path}")


if __name__ == "__main__":
    main()
