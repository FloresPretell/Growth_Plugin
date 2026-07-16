#!/usr/bin/env pvpython

import argparse
import re
import sys
from pathlib import Path

from paraview.simple import XMLUnstructuredGridReader, SaveData, Delete


def parse_args():
    parser = argparse.ArgumentParser(
        description="Rebuild Simulation.vtkhdf from Output_t*.vtu using ParaView only."
    )
    parser.add_argument("--root", required=True, help="Simulation2D directory")
    parser.add_argument("--output-name", default="Simulation.vtkhdf")
    parser.add_argument("--replace", action="store_true")
    return parser.parse_args()


def timestep_from_vtu(path: Path):
    match = re.search(r"Output_t(\d+)\.vtu$", path.name)
    if match is None:
        return None
    return int(match.group(1))


def collect_vtu_files(simdir: Path):
    merged_dir = simdir / "merged"

    files = []

    if merged_dir.is_dir():
        files = list(merged_dir.glob("Output_t*.vtu"))
        source = "merged"
    else:
        source = "Simulation2D"

    if not files:
        files = list(simdir.glob("Output_t*.vtu"))
        source = "Simulation2D"

    parsed = []
    skipped = []

    for file in files:
        timestep = timestep_from_vtu(file)
        if timestep is None:
            skipped.append(file)
            continue
        parsed.append((timestep, file.resolve()))

    parsed.sort(key=lambda item: item[0])

    return source, parsed, skipped


def main():
    args = parse_args()

    simdir = Path(args.root).expanduser().resolve()

    if not simdir.is_dir():
        raise SystemExit(f"ERROR: Simulation2D directory does not exist: {simdir}")

    final_output = simdir / args.output_name
    tmp_output = simdir / f".{args.output_name}.tmp.vtkhdf"
    backup_output = simdir / f"{args.output_name}.bak"

    source, parsed_files, skipped = collect_vtu_files(simdir)

    if skipped:
        print("[warning] Some files could not be parsed as Output_t<number>.vtu:")
        for file in skipped[:10]:
            print(f"          {file}")

    if not parsed_files:
        raise SystemExit(
            "ERROR: No valid Output_t*.vtu files found in either:\n"
            f"       {simdir / 'merged'}\n"
            f"       {simdir}"
        )

    timesteps = [t for t, _ in parsed_files]
    files = [str(path) for _, path in parsed_files]

    print("[info] Rebuilding VTKHDF with ParaView")
    print(f"[info] Simulation2D directory : {simdir}")
    print(f"[info] Source location        : {source}")
    print(f"[info] Number of VTU files    : {len(files)}")
    print(f"[info] First timestep         : {timesteps[0]}")
    print(f"[info] Last timestep          : {timesteps[-1]}")
    print(f"[info] First file             : {files[0]}")
    print(f"[info] Last file              : {files[-1]}")
    print(f"[info] Temporary output       : {tmp_output}")
    print(f"[info] Final output           : {final_output}")

    if timesteps[0] != 0 and final_output.exists():
        raise SystemExit(
            "ERROR: Existing Simulation.vtkhdf exists, but the available VTU series "
            f"starts at timestep {timesteps[0]}, not 0.\n"
            "This would overwrite old timesteps if --replace is used.\n"
            "Do not replace unless you intentionally want a partial VTKHDF."
        )

    if tmp_output.exists():
        print(f"[info] Removing old temporary file: {tmp_output}")
        tmp_output.unlink()

    reader = XMLUnstructuredGridReader(FileName=files)

    try:
        reader.TimestepValues = timesteps
        print("[info] Assigned explicit timestep values from filenames.")
    except Exception as exc:
        print(f"[warning] Could not assign TimestepValues: {exc}")
        print("[warning] Continuing with ParaView-inferred timesteps.")

    print("[info] Writing temporary VTKHDF...")
    SaveData(str(tmp_output), proxy=reader)
    Delete(reader)

    if not tmp_output.is_file():
        raise SystemExit(f"ERROR: Temporary VTKHDF was not created: {tmp_output}")

    if tmp_output.stat().st_size == 0:
        raise SystemExit(f"ERROR: Temporary VTKHDF is empty: {tmp_output}")

    print(f"[success] Temporary VTKHDF created: {tmp_output}")
    print(f"[success] Temporary size: {tmp_output.stat().st_size / 1e9:.3f} GB")

    if not args.replace:
        print("[info] --replace was not used.")
        print("[info] Existing Simulation.vtkhdf was not modified.")
        return 0

    print("[info] --replace enabled.")

    if final_output.exists():
        if backup_output.exists():
            print(f"[info] Removing previous backup: {backup_output}")
            backup_output.unlink()

        print(f"[info] Backing up old file to: {backup_output}")
        final_output.rename(backup_output)

    tmp_output.rename(final_output)

    if not final_output.is_file() or final_output.stat().st_size == 0:
        raise SystemExit(f"ERROR: Final VTKHDF is missing or empty: {final_output}")

    print(f"[success] Final VTKHDF written: {final_output}")

    if backup_output.exists():
        print(f"[info] Backup kept here: {backup_output}")

    return 0


if __name__ == "__main__":
    sys.exit(main())