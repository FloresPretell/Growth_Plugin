#!/usr/bin/env python
from __future__ import annotations

import argparse
import os
import sys
import time
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path
from statistics import mean

import numpy as np
import paraview.simple as pv
from paraview import servermanager

REPO_ROOT = Path(__file__).resolve().parents[1]
POSTANALYSIS_DIR = REPO_ROOT / "Postanalysis"
if str(POSTANALYSIS_DIR) not in sys.path:
    sys.path.insert(0, str(POSTANALYSIS_DIR))

import util_postanalisis_export_png as export_util


@dataclass
class TimestepProfile:
    index: int
    time_value: float
    fetch_seconds: float
    lsf_interp_seconds: float
    field_interp_seconds: float
    inhibitor_interp_seconds: float
    png_write_seconds: float
    total_seconds: float
    point_count: int
    cell_count: int


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Profile the local_safe PNG export path.")
    parser.add_argument("--simdir", required=True, help="Simulation merged directory containing Output_t*.vtu files")
    parser.add_argument("--fields", nargs="+", required=True, help="Fields to export, e.g. u_t")
    parser.add_argument("--max-files", type=int, default=2, help="Maximum number of timestep files to profile")
    parser.add_argument("--outdir", required=True, help="Temporary output directory for profiler PNGs")
    parser.add_argument(
        "--image-size",
        nargs=2,
        type=int,
        metavar=("WIDTH", "HEIGHT"),
        default=(1024, 1024),
        help="Image size used by the local_safe exporter",
    )
    return parser.parse_args()


def file_size_mb(path: Path) -> float:
    return path.stat().st_size / (1024.0 * 1024.0)


def format_seconds(value: float) -> str:
    return f"{value:.3f}s"


def fetch_arrays_with_metadata(reader, time_value: float, array_names: list[str]):
    start = time.perf_counter()
    pv.UpdatePipeline(time=time_value, proxy=reader)
    data = servermanager.Fetch(reader)
    fetch_seconds = time.perf_counter() - start
    if data is None or data.GetPoints() is None:
        raise RuntimeError(f"Failed to fetch ParaView data at time {time_value}")

    points = export_util.vtk_to_numpy(data.GetPoints().GetData())
    point_data = data.GetPointData()
    arrays = {}
    available_point_arrays = []
    for idx in range(point_data.GetNumberOfArrays()):
        array = point_data.GetArray(idx)
        if array is not None and array.GetName():
            available_point_arrays.append(array.GetName())

    cell_data = data.GetCellData()
    available_cell_arrays = []
    for idx in range(cell_data.GetNumberOfArrays()):
        array = cell_data.GetArray(idx)
        if array is not None and array.GetName():
            available_cell_arrays.append(array.GetName())

    for name in array_names:
        array = point_data.GetArray(name)
        if array is None:
            raise RuntimeError(f"Missing required array '{name}' at time {time_value}")
        arrays[name] = export_util.vtk_to_numpy(array)

    return {
        "fetch_seconds": fetch_seconds,
        "x": points[:, 0],
        "y": points[:, 1],
        "arrays": arrays,
        "point_count": int(data.GetNumberOfPoints()),
        "cell_count": int(data.GetNumberOfCells()),
        "available_point_arrays": available_point_arrays,
        "available_cell_arrays": available_cell_arrays,
    }


def classify_bottlenecks(
    total_dir_files: int,
    profiled_files: int,
    fields: list[str],
    used_pvtu: bool,
    avg_fetch: float,
    avg_interpolate: float,
    avg_write: float,
    point_count: int,
    avg_file_size_mb: float,
    overwrote_existing_outputs: bool,
) -> list[str]:
    classes: list[str] = []

    if total_dir_files > profiled_files:
        classes.append(
            f"B. It is processing too many timesteps in the default local run: {total_dir_files} available vs {profiled_files} profiled."
        )
    if len(fields) < 5:
        classes.append(
            f"C. Default local export still processes more fields than this profile: {len(fields)} profiled vs 5 scientific fields plus inhibitor in the full run."
        )
    if used_pvtu:
        classes.append("D. The reader is using partitioned .pvtu inputs.")
    if avg_fetch > max(avg_interpolate, avg_write):
        classes.append("E. ParaView UpdatePipeline/servermanager.Fetch is the dominant measured cost per timestep.")
    if avg_interpolate >= avg_fetch and avg_interpolate >= avg_write:
        classes.append("F. Python-side scattered interpolation is the dominant measured cost per timestep.")
    if avg_file_size_mb > 5 or point_count > 100000:
        classes.append(
            f"A. The dataset is non-trivial even for local profiling: mean file size {avg_file_size_mb:.2f} MB, point count about {point_count}."
        )
    if avg_write > 0.5:
        classes.append("H. PNG writing itself is a measurable part of the cost.")
    if overwrote_existing_outputs:
        classes.append("G. Existing profiler PNGs are overwritten; current local_safe export does not skip already-generated outputs.")
    return classes or ["J. Another concrete cause was not isolated beyond the measured fetch/interpolation/write split."]


def markdown_table(rows: list[list[str]]) -> str:
    if not rows:
        return ""
    header = rows[0]
    divider = ["---"] * len(header)
    lines = [
        "| " + " | ".join(header) + " |",
        "| " + " | ".join(divider) + " |",
    ]
    for row in rows[1:]:
        lines.append("| " + " | ".join(row) + " |")
    return "\n".join(lines)


def main() -> int:
    args = parse_args()
    simdir = Path(args.simdir).resolve()
    outdir = Path(args.outdir).resolve()
    report_dir = REPO_ROOT / "reports_IA"
    report_dir.mkdir(parents=True, exist_ok=True)
    outdir.mkdir(parents=True, exist_ok=True)
    os.environ.setdefault("MPLCONFIGDIR", str(outdir / ".matplotlib"))
    Path(os.environ["MPLCONFIGDIR"]).mkdir(parents=True, exist_ok=True)

    if not simdir.is_dir():
        raise SystemExit(f"ERROR: simdir does not exist: {simdir}")

    all_vtu_files = sorted(simdir.glob("Output_t*.vtu"))
    all_pvtu_files = sorted(simdir.glob("Output_t*.pvtu"))
    selected_files = [Path(p) for p in export_util._sorted_vtu_files(str(simdir))]
    if not selected_files:
        raise SystemExit(f"ERROR: current local_safe file selection found no VTU files in {simdir}")

    profiled_files = selected_files[: args.max_files]
    if not profiled_files:
        raise SystemExit("ERROR: no files selected after applying --max-files")

    print(f"simdir            : {simdir}", flush=True)
    print(f"profiled files    : {len(profiled_files)} / {len(selected_files)} selected by current local_safe code", flush=True)
    print(f"directory VTU     : {len(all_vtu_files)}", flush=True)
    print(f"directory PVTU    : {len(all_pvtu_files)}", flush=True)
    print(f"fields            : {', '.join(args.fields)}", flush=True)
    print(f"image_size        : {tuple(args.image_size)}", flush=True)
    print(f"outdir            : {outdir}", flush=True)

    pv.ResetSession()
    print("creating reader   : pv.XMLUnstructuredGridReader", flush=True)
    construct_start = time.perf_counter()
    reader = pv.XMLUnstructuredGridReader(FileName=[str(path) for path in profiled_files])
    reader.PointArrayStatus = args.fields + ["lsf", "inhibitor"]
    reader.TimeArray = "None"
    construct_seconds = time.perf_counter() - construct_start

    print("loading bounds    : first metadata load", flush=True)
    bounds_start = time.perf_counter()
    bounds = export_util._initial_bounds(reader)
    bounds_seconds = time.perf_counter() - bounds_start

    print("building grid     : numpy meshgrid", flush=True)
    mesh_start = time.perf_counter()
    x_lin = np.linspace(bounds[0], bounds[1], args.image_size[0])
    y_lin = np.linspace(bounds[2], bounds[3], args.image_size[1])
    grid_x, grid_y = np.meshgrid(x_lin, y_lin)
    mesh_seconds = time.perf_counter() - mesh_start

    time_values = list(reader.TimestepValues)
    selected_time_values = time_values[: len(profiled_files)]

    timestep_profiles: list[TimestepProfile] = []
    available_point_arrays: list[str] = []
    available_cell_arrays: list[str] = []
    overwrote_existing_outputs = False

    for field in args.fields + ["inhibitor"]:
        (outdir / field).mkdir(parents=True, exist_ok=True)

    for idx, time_value in enumerate(selected_time_values):
        print(f"timestep {idx}     : fetch start", flush=True)
        total_start = time.perf_counter()
        fetched = fetch_arrays_with_metadata(reader, time_value, args.fields + ["lsf", "inhibitor"])
        print(f"timestep {idx}     : fetch done in {fetched['fetch_seconds']:.3f}s", flush=True)
        if not available_point_arrays:
            available_point_arrays = fetched["available_point_arrays"]
            available_cell_arrays = fetched["available_cell_arrays"]

        print(f"timestep {idx}     : interpolate lsf", flush=True)
        lsf_start = time.perf_counter()
        lsf_grid = export_util._interpolate_to_grid(fetched["x"], fetched["y"], fetched["arrays"]["lsf"], grid_x, grid_y)
        lsf_interp_seconds = time.perf_counter() - lsf_start
        print(f"timestep {idx}     : lsf done in {lsf_interp_seconds:.3f}s", flush=True)

        field_interp_seconds = 0.0
        png_write_seconds = 0.0

        for field in args.fields:
            print(f"timestep {idx}     : interpolate field {field}", flush=True)
            interp_start = time.perf_counter()
            scalar_grid = export_util._interpolate_to_grid(
                fetched["x"], fetched["y"], fetched["arrays"][field], grid_x, grid_y
            )
            field_duration = time.perf_counter() - interp_start
            field_interp_seconds += field_duration
            print(f"timestep {idx}     : field {field} interp done in {field_duration:.3f}s", flush=True)

            output_path = outdir / field / f"{idx:03d}.png"
            overwrote_existing_outputs = overwrote_existing_outputs or output_path.exists()
            print(f"timestep {idx}     : write field {field} -> {output_path.name}", flush=True)
            write_start = time.perf_counter()
            export_util._save_local_safe_png(
                str(output_path),
                scalar_grid,
                lsf_grid < 0,
                bounds,
                0,
                export_util.DEFAULT_MAX_RANGES[field],
                image_size=tuple(args.image_size),
            )
            write_duration = time.perf_counter() - write_start
            png_write_seconds += write_duration
            print(f"timestep {idx}     : field {field} write done in {write_duration:.3f}s", flush=True)

        print(f"timestep {idx}     : interpolate inhibitor", flush=True)
        inhibitor_interp_start = time.perf_counter()
        inhibitor_grid = export_util._interpolate_to_grid(
            fetched["x"], fetched["y"], fetched["arrays"]["inhibitor"], grid_x, grid_y
        )
        inhibitor_interp_seconds = time.perf_counter() - inhibitor_interp_start
        print(f"timestep {idx}     : inhibitor interp done in {inhibitor_interp_seconds:.3f}s", flush=True)

        inhibitor_output = outdir / "inhibitor" / f"{idx:03d}.png"
        overwrote_existing_outputs = overwrote_existing_outputs or inhibitor_output.exists()
        print(f"timestep {idx}     : write inhibitor -> {inhibitor_output.name}", flush=True)
        inhibitor_write_start = time.perf_counter()
        export_util._save_local_safe_png(
            str(inhibitor_output),
            inhibitor_grid,
            lsf_grid > 0,
            bounds,
            0,
            export_util.DEFAULT_MAX_RANGES["inhibitor"],
            image_size=tuple(args.image_size),
        )
        inhibitor_write_duration = time.perf_counter() - inhibitor_write_start
        png_write_seconds += inhibitor_write_duration
        print(f"timestep {idx}     : inhibitor write done in {inhibitor_write_duration:.3f}s", flush=True)

        timestep_total_seconds = time.perf_counter() - total_start
        print(f"timestep {idx}     : total done in {timestep_total_seconds:.3f}s", flush=True)
        timestep_profiles.append(
            TimestepProfile(
                index=idx,
                time_value=float(time_value),
                fetch_seconds=fetched["fetch_seconds"],
                lsf_interp_seconds=lsf_interp_seconds,
                field_interp_seconds=field_interp_seconds,
                inhibitor_interp_seconds=inhibitor_interp_seconds,
                png_write_seconds=png_write_seconds,
                total_seconds=timestep_total_seconds,
                point_count=fetched["point_count"],
                cell_count=fetched["cell_count"],
            )
        )

    avg_fetch = mean(item.fetch_seconds for item in timestep_profiles)
    avg_lsf = mean(item.lsf_interp_seconds for item in timestep_profiles)
    avg_field_interp = mean(item.field_interp_seconds for item in timestep_profiles)
    avg_inhib_interp = mean(item.inhibitor_interp_seconds for item in timestep_profiles)
    avg_write = mean(item.png_write_seconds for item in timestep_profiles)
    avg_total = mean(item.total_seconds for item in timestep_profiles)
    avg_interpolate = avg_lsf + avg_field_interp + avg_inhib_interp
    avg_file_size = mean(file_size_mb(path) for path in profiled_files)
    point_count = timestep_profiles[0].point_count if timestep_profiles else 0
    cell_count = timestep_profiles[0].cell_count if timestep_profiles else 0
    used_pvtu = any(path.suffix == ".pvtu" for path in profiled_files)

    classifications = classify_bottlenecks(
        total_dir_files=len(selected_files),
        profiled_files=len(profiled_files),
        fields=args.fields,
        used_pvtu=used_pvtu,
        avg_fetch=avg_fetch,
        avg_interpolate=avg_interpolate,
        avg_write=avg_write,
        point_count=point_count,
        avg_file_size_mb=avg_file_size,
        overwrote_existing_outputs=overwrote_existing_outputs,
    )

    print("")
    print("Summary", flush=True)
    print(f"reader construct   : {format_seconds(construct_seconds)}", flush=True)
    print(f"initial bounds     : {format_seconds(bounds_seconds)}", flush=True)
    print(f"grid build         : {format_seconds(mesh_seconds)}", flush=True)
    print(f"avg fetch          : {format_seconds(avg_fetch)}", flush=True)
    print(f"avg lsf interp     : {format_seconds(avg_lsf)}", flush=True)
    print(f"avg field interp   : {format_seconds(avg_field_interp)}", flush=True)
    print(f"avg inhib interp   : {format_seconds(avg_inhib_interp)}", flush=True)
    print(f"avg png write      : {format_seconds(avg_write)}", flush=True)
    print(f"avg total/timestep : {format_seconds(avg_total)}", flush=True)
    print(f"points/cells       : {point_count} / {cell_count}", flush=True)
    print(f"available arrays   : {', '.join(available_point_arrays)}", flush=True)

    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    report_path = report_dir / f"local_safe_performance_report_{timestamp}.md"

    file_rows = [["file", "size_mb", "selected_for_profile"]]
    profiled_set = {path.resolve() for path in profiled_files}
    for path in selected_files:
        file_rows.append(
            [
                str(path.name),
                f"{file_size_mb(path):.3f}",
                "yes" if path.resolve() in profiled_set else "no",
            ]
        )

    profile_rows = [
        [
            "idx",
            "time_value",
            "fetch_s",
            "lsf_interp_s",
            "field_interp_s",
            "inhibitor_interp_s",
            "png_write_s",
            "total_s",
            "points",
            "cells",
        ]
    ]
    for item in timestep_profiles:
        profile_rows.append(
            [
                str(item.index),
                f"{item.time_value:.6g}",
                f"{item.fetch_seconds:.3f}",
                f"{item.lsf_interp_seconds:.3f}",
                f"{item.field_interp_seconds:.3f}",
                f"{item.inhibitor_interp_seconds:.3f}",
                f"{item.png_write_seconds:.3f}",
                f"{item.total_seconds:.3f}",
                str(item.point_count),
                str(item.cell_count),
            ]
        )

    report_lines = [
        "# local_safe Performance Report",
        "",
        f"- Timestamp: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}",
        f"- Command tested: python tools/profile_local_safe_export.py --simdir {simdir} --fields {' '.join(args.fields)} --max-files {args.max_files} --outdir {outdir}",
        "- Config used: standalone profiler against the current `local_safe` implementation in `Postanalysis/util_postanalisis_export_png.py`",
        f"- Case tested: `{simdir}`",
        f"- Fields exported in profile: `{', '.join(args.fields)}` plus `inhibitor`",
        f"- Timesteps selected for profile: {len(selected_time_values)}",
        f"- Timesteps available to current local_safe selection: {len(selected_files)}",
        f"- Reader used: `pv.XMLUnstructuredGridReader`",
        f"- Data source type selected by current code: `{profiled_files[0].suffix}` files from the merged directory",
        f"- Uses ParaView Fetch: yes (`servermanager.Fetch` per timestep)",
        f"- Uses Python scattered interpolation: yes (`scipy.interpolate.griddata` for `lsf`, each field, and `inhibitor`)",
        f"- Existing outputs regenerated if re-run: {'yes' if overwrote_existing_outputs else 'no during this run, but code overwrites when files exist'}",
        "",
        "## File Inventory",
        "",
        markdown_table(file_rows),
        "",
        "## Available Arrays",
        "",
        f"- Point arrays: {', '.join(available_point_arrays)}",
        f"- Cell arrays: {', '.join(available_cell_arrays) if available_cell_arrays else '<none>'}",
        f"- Points / cells in profiled timestep: {point_count} / {cell_count}",
        "",
        "## Timings",
        "",
        f"- Reader construction time: {construct_seconds:.3f}s",
        f"- Initial metadata/bounds load time: {bounds_seconds:.3f}s",
        f"- Grid build time: {mesh_seconds:.3f}s",
        f"- Average fetch time per timestep: {avg_fetch:.3f}s",
        f"- Average `lsf` interpolation time per timestep: {avg_lsf:.3f}s",
        f"- Average field interpolation time per timestep: {avg_field_interp:.3f}s",
        f"- Average inhibitor interpolation time per timestep: {avg_inhib_interp:.3f}s",
        f"- Average PNG write time per timestep: {avg_write:.3f}s",
        f"- Average total time per timestep: {avg_total:.3f}s",
        "",
        markdown_table(profile_rows),
        "",
        "## Bottleneck Classification",
        "",
    ]
    report_lines.extend(f"- {item}" for item in classifications)
    report_lines.extend(
        [
            "",
            "## Interpretation",
            "",
            "- The current `local_safe` implementation does one `UpdatePipeline` + `servermanager.Fetch` per timestep.",
            f"- It then does {1 + len(args.fields) + 1} scattered-grid interpolations per profiled timestep: one for `lsf`, {len(args.fields)} for the requested field set, and one for `inhibitor`.",
            "- In the default local runner, the field list is larger than this profile, so full local runtime scales roughly with both timestep count and field count.",
            "- The current code does not skip existing PNG outputs; reruns regenerate files.",
            "",
            "## Files Changed",
            "",
            "- `tools/profile_local_safe_export.py`",
            "",
            "## Shaheen Impact",
            "",
            "- None. This profiler is external and does not change the `render` export path or `run_shaheen.sh` behavior.",
        ]
    )

    report_path.write_text("\n".join(report_lines), encoding="utf-8")
    print(f"report            : {report_path}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
