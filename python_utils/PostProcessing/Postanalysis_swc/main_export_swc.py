#!/usr/bin/env python3

import argparse
from pathlib import Path

from util_postanalisis_export_swc import export_one


def infer_case_code(simdir: Path) -> str:
    if simdir.name == "u_t":
        return simdir.parent.name
    parts = simdir.parts
    if "Simulation2D" in parts:
        sim_idx = parts.index("Simulation2D")
        if sim_idx > 0:
            return parts[sim_idx - 1]
    if "Cases" in parts:
        cases_idx = parts.index("Cases")
        if cases_idx + 1 < len(parts):
            return parts[cases_idx + 1]
    raise ValueError(f"Could not infer case code from simdir: {simdir}")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--simdir", required=True, help="Path to one case PNG directory")
    parser.add_argument("--outdir", required=True, help="Root directory for exported SWC results")
    args = parser.parse_args()

    simdir = Path(args.simdir).expanduser().resolve()
    outroot = Path(args.outdir).expanduser().resolve()

    if not simdir.is_dir():
        raise SystemExit(f"ERROR: Input PNG directory does not exist: {simdir}")

    case_code = infer_case_code(simdir)
    final_outdir = outroot / case_code
    final_outdir.mkdir(parents=True, exist_ok=True)

    print(f"simdir      : {simdir}")
    print(f"case_code   : {case_code}")
    print(f"final_outdir: {final_outdir}")

    export_one(str(simdir), str(final_outdir))


if __name__ == "__main__":
    main()
