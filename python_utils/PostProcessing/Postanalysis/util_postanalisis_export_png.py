import paraview.simple as pv
import glob
import os
from pathlib import Path


DEFAULT_MAX_RANGES = {
    "u_t": 2,
    "u_u": 1.5,
    "u_b": 1.5,
    "u_p": 0.002,
    "u_ca_cyt": 0.9,
    "inhibitor": 0.01,
}


def Create_view_white(image_size=(2048, 2048)):
    """Create ParaView render view with black background"""
    try:
        renderView = pv.GetActiveViewOrCreate('RenderView')
        renderView.ViewSize = list(image_size)
        renderView.Background = [0, 0, 0]  # Black background
        renderView.BackgroundColorMode = 'Single Color'
        renderView.OrientationAxesVisibility = 0  # No axes
        renderView.CameraParallelProjection = 1
        return renderView
    except Exception as e:
        print(f"❌ Error in Create_view_white: {e}")
        return None


def calculate_perfect_parallel_scale(width, height, image_size, margin_factor=0):
    """Calculate perfect parallel scale to fit geometry in image"""
    geom_ratio = width / height
    image_ratio = image_size[0] / image_size[1]

    if geom_ratio > image_ratio:
        # Geometry wider than image → width limits
        scale = (width / 2) * (1 + margin_factor)
    else:
        # Geometry taller than image → height limits
        scale = (height / 2) * (1 + margin_factor)

    return scale


def Camara_position(reader, renderView, image_size=(2048, 2048)):
    """Configure camera to center and fit geometry"""
    try:
        t0 = reader.TimestepValues[0]
        pv.UpdatePipeline(time=t0, proxy=reader)

        # Get bounds
        data_info = reader.GetDataInformation()
        bounds = data_info.GetBounds()

        width = bounds[1] - bounds[0]
        height = bounds[3] - bounds[2]
        center_x = (bounds[0] + bounds[1]) / 2
        center_y = (bounds[2] + bounds[3]) / 2
        center_z = (bounds[4] + bounds[5]) / 2

        # Configure camera
        camera = renderView.GetActiveCamera()
        camera.SetPosition([center_x, center_y, center_z + max(width, height) * 1.5])
        camera.SetFocalPoint([center_x, center_y, center_z])
        camera.SetViewUp([0, 1, 0])

        # Perfect parallel scale
        optimal_scale = calculate_perfect_parallel_scale(width, height, image_size, margin_factor=0.0)
        camera.SetParallelScale(optimal_scale)

        return camera

    except Exception as e:
        print(f"❌ Error in Camara_position: {e}")
        return None


def export_component_pngs(variable, reader, renderView, display_resample, base_output_dir, min_range=0.0, max_range=2.0):
    """Export PNG frames for a single variable"""
    # Set coloring
    pv.ColorBy(display_resample, ('POINTS', variable))

    # Colormap setup
    lut = pv.GetColorTransferFunction(variable, display_resample)

    # Try multiple colormap names (different ParaView versions have different names)
    colormap_names = [
        'Inferno (matplotlib)',  # Try this first
        'Inferno',               # Fallback
        'jet',                   # Generic fallback
        'viridis',               # Another option
    ]
    
    colormap_applied = False
    for cmap_name in colormap_names:
        try:
            lut.ApplyPreset(cmap_name, True)
            print(f"  [info] Applied colormap: {cmap_name}")
            colormap_applied = True
            break
        except RuntimeError:
            continue
    
    if not colormap_applied:
        print(f"  [warning] Could not apply any colormap preset; using default")
        # Just use the default without ApplyPreset
    
    
    lut.InvertTransferFunction()
    lut.RescaleTransferFunction(min_range, max_range)
    display_resample.LookupTable = lut

    # Output folder
    out_directory_component = os.path.join(base_output_dir, variable)
    os.makedirs(out_directory_component, exist_ok=True)

    # Save screenshots for each timestep
    for idx, t in enumerate(reader.TimestepValues):
        print(f"  Saving {variable} timestep {idx} / {len(reader.TimestepValues)}")
        renderView.ViewTime = t
        pv.Render()
        pv.SaveScreenshot(
            os.path.join(out_directory_component, f"{idx:03d}.png"),
            renderView,
            ImageResolution=renderView.ViewSize,
            TransparentBackground=1,
        )


def _sorted_vtu_files(vtu_file):
    """Collect VTU files from directory or pattern"""
    if os.path.isdir(vtu_file):
        return sorted(glob.glob(os.path.join(vtu_file, "Output_t*.vtu")))
    return sorted(glob.glob(vtu_file))


def _detect_and_open_file(file_path):
    """
    Open either VTKHDF or VTU files.
    
    Args:
        file_path: Path to Simulation.vtkhdf OR directory/pattern for VTU files
    
    Returns:
        ParaView reader object
    """
    file_path = os.path.abspath(file_path)

    # === Check for VTKHDF ===
    if file_path.endswith(".vtkhdf") and os.path.isfile(file_path):
        print(f"[info] Detected VTKHDF file: {file_path}")
        try:
            reader = pv.OpenDataFile(file_path)
            print(f"[success] Opened VTKHDF: {file_path}")
            return reader
        except Exception as e:
            print(f"[error] Failed to open VTKHDF: {e}")
            raise SystemExit(f"ERROR: Cannot open VTKHDF: {e}") from e

    # === Check for VTU files (backward compatibility) ===
    vtu_files = _sorted_vtu_files(file_path)
    if vtu_files:
        print(f"[info] Detected VTU file series: {len(vtu_files)} files")
        reader = pv.XMLUnstructuredGridReader(FileName=vtu_files)
        print(f"[success] Opened VTU series")
        return reader

    # === No valid input found ===
    raise SystemExit(
        f"ERROR: No valid input found at {file_path}\n"
        f"       Expected: Simulation.vtkhdf OR Output_t*.vtu files"
    )


def export_png_from_vtu(vtu_file, outdir, fields_to_export=["u_u", "u_p", "u_t", "u_b", "u_ca_cyt"], 
                        min_range=0, max_ranges=DEFAULT_MAX_RANGES):
    """
    Export PNGs from VTKHDF or VTU files.
    """
    pv.ResetSession()
    renderView = Create_view_white()

    # === Open file (VTKHDF or VTU) ===
    reader = _detect_and_open_file(vtu_file)

    # Configure reader (handle both VTKHDF and XMLUnstructuredGridReader)
    reader.PointArrayStatus = fields_to_export + ["lsf", "inhibitor"]
    
    # TimeArray only exists for XMLUnstructuredGridReader, NOT VTKHDFReader
    reader_type = type(reader).__name__
    if reader_type == "XMLUnstructuredGridReader":
        reader.TimeArray = "None"
    else:
        print(f"[info] Reader type: {reader_type} (no TimeArray property)")

    # Show and configure camera
    tmp_display = pv.Show(reader, renderView)
    pv.Render()
    camara = Camara_position(reader, renderView)
    pv.Hide(reader, renderView)
    pv.Render()

    # Create output directory
    os.makedirs(outdir, exist_ok=True)

    # === Rest of code unchanged ===
    contour = pv.Contour(Input=reader)
    contour.ContourBy = ['POINTS', 'lsf']
    contour.Isosurfaces = [0.0]
    display_lsf = pv.Show(contour, renderView)
    display_lsf.ColorArrayName = ['POINTS', '']
    display_lsf.Representation = 'Surface'
    display_lsf.DiffuseColor = [0, 0, 0]
    display_lsf.LineWidth = 2
    pv.Hide(display_lsf, renderView)

    thresh = pv.Threshold(Input=reader)
    thresh.Scalars = ['POINTS', 'lsf']
    thresh.LowerThreshold = 0
    thresh.ThresholdMethod = 'Below Lower Threshold'

    display_th = pv.Show(thresh, renderView)
    display_th.Representation = 'Surface'
    pv.Hide(display_th, renderView)

    resample = pv.ResampleWithDataset(SourceDataArrays=reader, DestinationMesh=thresh)
    resample.CellLocator = 'Static Cell Locator'
    display_resample = pv.Show(resample, renderView)

    # === Export intracellular fields ===
    for comp in fields_to_export:
        print(f"Exporting {comp}...")
        export_component_pngs(
            variable=comp,
            reader=reader,
            renderView=renderView,
            display_resample=display_resample,
            base_output_dir=outdir,
            min_range=min_range,
            max_range=max_ranges[comp],
        )

    pv.Hide(resample, renderView)

    # === Setup threshold for extracellular region (lsf > 0) ===
    thresh_inhib = pv.Threshold(Input=reader)
    thresh_inhib.Scalars = ['POINTS', 'lsf']
    thresh_inhib.UpperThreshold = 0
    thresh_inhib.ThresholdMethod = 'Above Upper Threshold'

    display_th_inhib = pv.Show(thresh_inhib, renderView)
    display_th_inhib.Representation = 'Surface'
    pv.Hide(thresh_inhib, renderView)

    resample_inhib = pv.ResampleWithDataset(SourceDataArrays=reader, DestinationMesh=thresh_inhib)
    resample_inhib.CellLocator = 'Static Cell Locator'
    display_resample_inhib = pv.Show(resample_inhib, renderView)
    pv.Hide(resample_inhib, renderView)

    # === Export inhibitor field ===
    print(f"Exporting inhibitor...")
    export_component_pngs(
        variable="inhibitor",
        reader=reader,
        renderView=renderView,
        display_resample=display_resample_inhib,
        base_output_dir=outdir,
        min_range=min_range,
        max_range=max_ranges["inhibitor"],
    )


def export_one(path, outdir, fields_to_export=["u_u", "u_p", "u_t", "u_b", "u_ca_cyt"]):
    """Wrapper: export from path to outdir"""
    print(f"🚀 Exporting PNGs from: {path}")
    export_png_from_vtu(path, outdir, fields_to_export=fields_to_export, min_range=0, max_ranges=DEFAULT_MAX_RANGES)
    print(f"✅ Done: {path}")


def export_one_local_safe(path, outdir, fields_to_export=["u_u", "u_p", "u_t", "u_b", "u_ca_cyt"], max_timesteps=None):
    """Placeholder: local_safe not available without numpy"""
    raise SystemExit(
        "ERROR: local_safe mode requires numpy and scipy.\n"
        "       Use --export-mode render instead (works without numpy)."
    )