#!/usr/bin/env python3
"""
Create a manifest of case Simulation2D directories with Simulation.vtkhdf.

Usage:
    python3 detect_vtkhdf_cases.py --path /data/cases --output manifest.txt
"""

import argparse
from pathlib import Path


def parse_args():
    parser = argparse.ArgumentParser(
        description="Create a manifest of cases with Simulation2D/Simulation.vtkhdf"
    )
    parser.add_argument("--path", required=True, help="Root directory containing D*_V* cases")
    parser.add_argument("--output", required=True, help="Output manifest file path")
    parser.add_argument(
        "--max-cases", default="", help="Optional limit on detected cases"
    )
    parser.add_argument(
        "--verbose", action="store_true", help="Print file sizes"
    )
    return parser.parse_args()


def detect_case_dirs(dataset_root: Path, verbose: bool = False):
    """Find all D*_V* directories containing Simulation.vtkhdf"""
    case_dirs = []

    for candidate in sorted(dataset_root.iterdir()):
        if not candidate.is_dir():
            continue

        # Check naming pattern D*_V*
        if not candidate.name.startswith("D") or "_V" not in candidate.name:
            continue

        sim_dir = candidate / "Simulation2D"

        if not sim_dir.is_dir():
            print(f"[skip] {candidate.name}: no Simulation2D/")
            continue

        # Check for Simulation.vtkhdf
        vtkhdf_path = sim_dir / "Simulation.vtkhdf"

        if vtkhdf_path.exists() and vtkhdf_path.is_file():
            size_gb = vtkhdf_path.stat().st_size / 1e9
            if verbose:
                print(f"[found] {sim_dir} ({size_gb:.3f} GB)")
            else:
                print(f"[found] {sim_dir}")
            case_dirs.append(sim_dir.resolve())
            continue

        print(f"[skip] {candidate.name}: no Simulation.vtkhdf")

    return case_dirs


def main():
    args = parse_args()
    dataset_root = Path(args.path).expanduser().resolve()
    output_path = Path(args.output).expanduser().resolve()

    if not dataset_root.is_dir():
        raise SystemExit(f"ERROR: Dataset root does not exist: {dataset_root}")

    print(f"[info] Scanning: {dataset_root}")
    cases = detect_case_dirs(dataset_root, verbose=args.verbose)

    if not cases:
        raise SystemExit("ERROR: No cases with Simulation.vtkhdf found")

    # Apply max limit
    if args.max_cases:
        try:
            max_cases = int(args.max_cases)
        except ValueError as exc:
            raise SystemExit(f"ERROR: Invalid --max-cases: {args.max_cases}") from exc
        if max_cases < 1:
            raise SystemExit("ERROR: --max-cases must be > 0")
        cases = cases[:max_cases]

    # Write manifest
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", encoding="utf-8") as f:
        for case in cases:
            f.write(f"{case}\n")

    print(f"\n[success] Wrote {len(cases)} cases to {output_path}")


if __name__ == "__main__":
    main()