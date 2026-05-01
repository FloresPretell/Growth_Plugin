#!/usr/bin/env pvpython
"""
Main PNG export script for ParaView VTKHDF/VTU files.
Supports both VTKHDF (preferred) and VTU fallback.
"""

import argparse
from pathlib import Path
from util_postanalisis_export_png import export_one, export_one_local_safe


def main():
    ap = argparse.ArgumentParser(
        description="Export PNGs from VTKHDF or VTU files"
    )
    ap.add_argument(
        "--simdir", required=True,help="Path to Simulation2D directory (containing Simulation.vtkhdf or Output_t*.vtu)")
    ap.add_argument(
        "--outdir", required=True,help="Root folder to store exported PNGs (case subdirectory will be created)")
    ap.add_argument(
        "--export-mode", choices=["render", "local_safe"], default="render",help="Export mode: 'render' (ParaView native, no numpy), 'local_safe' (requires numpy)")
    ap.add_argument(
        "--max-timesteps", type=int, default=None,help="Optional limit on timesteps for testing")
    ap.add_argument(
        "--fields_to_export", nargs="+",default=["u_u", "u_p", "u_t", "u_b", "u_ca_cyt"],help="Fields to export")

    args = ap.parse_args()

    simdir = Path(args.simdir).resolve()
    outroot = Path(args.outdir).resolve()

    # Extract case code from parent directory
    # Expected: /path/Cases/D2.5_V2p5/Simulation2D
    case_code = simdir.parent.name
    final_outdir = outroot / case_code
    final_outdir.mkdir(parents=True, exist_ok=True)

    print(f"simdir        : {simdir}")
    print(f"case_code     : {case_code}")
    print(f"final_outdir  : {final_outdir}")
    print(f"export_mode   : {args.export_mode}")
    print(f"fields        : {args.fields_to_export}")
    print(f"max_timesteps : {args.max_timesteps if args.max_timesteps else 'all'}")

    # === Detect input source (VTKHDF or VTU) ===
    vtkhdf_path = simdir / "Simulation.vtkhdf"

    if vtkhdf_path.exists() and vtkhdf_path.stat().st_size > 0:
        input_source = str(vtkhdf_path)
        print(f"[info] Using VTKHDF: {input_source}")
    else:
        # Fallback to VTU
        input_source = str(simdir)
        print(f"[info] Using VTU files from: {input_source}")

    # === Export ===
    if args.export_mode == "local_safe":
        export_one_local_safe(
            input_source,
            str(final_outdir),
            fields_to_export=args.fields_to_export,
            max_timesteps=args.max_timesteps,
        )
    else:  # render mode
        export_one(
            input_source,
            str(final_outdir),
            fields_to_export=args.fields_to_export,
        )

    print(f"[success] PNG export complete: {final_outdir}")


if __name__ == "__main__":
    main()