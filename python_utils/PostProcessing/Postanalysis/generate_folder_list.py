#!/usr/bin/env python3

import argparse
from pathlib import Path


def parse_args():
    parser = argparse.ArgumentParser(
        description="Create a manifest of case Simulation2D/merged directories."
    )
    parser.add_argument("--path", required=True, help="Root directory containing D*_V* cases")
    parser.add_argument("--output", required=True, help="Output manifest file path")
    parser.add_argument("--max-cases", default="", help="Optional limit on detected cases")
    return parser.parse_args()


def detect_case_dirs(dataset_root: Path):
    case_dirs = []

    for candidate in sorted(dataset_root.iterdir()):
        if not candidate.is_dir():
            continue

        if not candidate.name.startswith("D") or "_V" not in candidate.name:
            continue

        sim_dir = candidate / "Simulation2D"
        merged_dir = sim_dir / "merged"

        if not sim_dir.is_dir():
            print(f"Skipping {candidate}: missing Simulation2D")
            continue

        has_vtu_in_merged = False
        has_vtu_in_simdir = False

        if merged_dir.is_dir():
            has_vtu_in_merged = (
                any(merged_dir.glob("Output_t*.vtu"))
                or any(merged_dir.glob("Output_t*.pvtu"))
            )

        has_vtu_in_simdir = (
            any(sim_dir.glob("Output_t*.vtu"))
            or any(sim_dir.glob("Output_t*.pvtu"))
        )

        if has_vtu_in_merged:
            print(f"Detected VTU/PVTU in merged: {merged_dir}")
            case_dirs.append(sim_dir.resolve())
            continue

        if has_vtu_in_simdir:
            print(f"Detected VTU/PVTU in Simulation2D: {sim_dir}")
            case_dirs.append(sim_dir.resolve())
            continue

        print(
            f"Skipping {candidate}: no Output_t*.vtu or Output_t*.pvtu files found "
            f"in Simulation2D/merged or Simulation2D"
        )

    return case_dirs


def main():
    args = parse_args()
    dataset_root = Path(args.path).expanduser().resolve()
    output_path = Path(args.output).expanduser().resolve()

    if not dataset_root.is_dir():
        raise SystemExit(f"ERROR: Dataset root does not exist: {dataset_root}")

    cases = detect_case_dirs(dataset_root)

    if args.max_cases:
        try:
            max_cases = int(args.max_cases)
        except ValueError as exc:
            raise SystemExit(f"ERROR: Invalid --max-cases value: {args.max_cases}") from exc
        if max_cases < 1:
            raise SystemExit("ERROR: --max-cases must be greater than zero")
        cases = cases[:max_cases]

    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", encoding="utf-8") as handle:
        for case in cases:
            handle.write(f"{case}\n")

    print(f"Wrote {len(cases)} paths to {output_path}")


if __name__ == "__main__":
    main()
