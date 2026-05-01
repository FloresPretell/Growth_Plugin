def Create_view_white(image_size=(2048, 2048)):
    try:
        renderView = pv.GetActiveViewOrCreate('RenderView')
        renderView.ViewSize = list(image_size)
        renderView.Background = [0,0,0]  # Fondo negro
        renderView.BackgroundColorMode = 'Single Color'  
        renderView.OrientationAxesVisibility = 0  # Sin ejes
        renderView.CameraParallelProjection = 1
        ##pv.LoadPalette(paletteName='WhiteBackground')
        return renderView

    except Exception as e:
        print(f" Error en BackgroundWhite: {e}")
    
def calculate_perfect_parallel_scale(width, height, image_size, margin_factor=0):
    """Calcula escala paralela perfecta para llenar la imagen"""
    geom_ratio = width / height
    image_ratio = image_size[0] / image_size[1]
    
    if geom_ratio > image_ratio:
        # Geometría más ancha que imagen → width limita
        scale = (width / 2) * (1 + margin_factor)
        limiting_dimension = "WIDTH"
    else:
        # Geometría más alta que imagen → height limita  
        scale = (height / 2) * (1 + margin_factor)
        limiting_dimension = "HEIGHT"
    
    return scale  

def Camara_position(reader, renderView,image_size=(2048, 2048)):
    try:

        t0 = reader.TimestepValues[0]
        pv.UpdatePipeline(time=t0, proxy=reader) # fix the camara to the initial one before any resample data
                
        # Obtener bounds sin mostrar
        data_info = reader.GetDataInformation()
        bounds = data_info.GetBounds()

    
        width = bounds[1] - bounds[0]
        height = bounds[3] - bounds[2]
        center_x = (bounds[0] + bounds[1]) / 2
        center_y = (bounds[2] + bounds[3]) / 2
        center_z = (bounds[4] + bounds[5]) / 2
        
        # Configurar cámara perfecta
        camera = renderView.GetActiveCamera()
        camera.SetPosition([center_x, center_y, center_z + max(width, height) * 1.5])
        camera.SetFocalPoint([center_x, center_y, center_z])
        camera.SetViewUp([0, 1, 0])
        
        

        # ¡CLAVE! Escala paralela perfecta
        optimal_scale = calculate_perfect_parallel_scale(width, height, image_size, margin_factor=0.0)
        camera.SetParallelScale(optimal_scale)
        
        return camera

    except Exception as e:
        print(f"❌ Error en Camara_position: {e}")
        return None


def export_component_pngs(variable,reader,renderView,display_resample,base_output_dir,min_range=0.0,max_range=2.0):
    # Set coloring array
    pv.ColorBy(display_resample, ('POINTS', variable))

    # LUT
    lut = pv.GetColorTransferFunction(variable, display_resample)
    lut.ApplyPreset('Inferno (matplotlib)', True)
    lut.InvertTransferFunction()
    lut.RescaleTransferFunction(min_range, max_range)

    display_resample.LookupTable = lut

    # Output folder for this component
    out_directory_component = os.path.join(base_output_dir, variable)
    os.makedirs(out_directory_component, exist_ok=True)

    # Save screenshots
    for idx, t in enumerate(reader.TimestepValues):
        print(f"Saving {variable} timestep {idx} / {len(reader.TimestepValues)}")
        renderView.ViewTime = t
        pv.Render()
        pv.SaveScreenshot(
            os.path.join(out_directory_component, f"{idx:03d}.png"),
            renderView,
            ImageResolution=renderView.ViewSize,
            TransparentBackground=1)
        
############
import paraview.simple as pv
import glob
import os


def export_one(path, fields_to_export = ["u_u", "u_p", "u_t", "u_b", "u_ca_cyt"] ):
    print(f"🚀 Exporting PNGs from: {path}")
    export_png_from_vtu(path, fields_to_export = fields_to_export , min_range=0, max_ranges={'u_t': 2, 'u_u': 1.5, 'u_b': 1.5, 'u_p': 0.002, 'u_ca_cyt': 0.9, 'inhibitor': 0.01})
    print(f"✅ Done: {path}")
    


def export_png_from_vtu(vtu_file, fields_to_export = ["u_u", "u_p", "u_t", "u_b", "u_ca_cyt"], min_range=0, max_ranges= {'u_t': 2, 'u_u': 1.5, 'u_b': 1.5, 'u_p': 0.002, 'u_ca_cyt': 0.9, 'inhibitor': 0.01}):
    pv.ResetSession()
    renderView = Create_view_white()

   #look for the pvtu files
    if os.path.isdir(vtu_file):
        files = sorted(glob.glob(os.path.join(vtu_file, "Output_t*.vtu"))) ## change to vtu
    else:
        files = sorted(glob.glob(vtu_file))
        
    # Leer archivo VTU
    reader = pv.XMLUnstructuredGridReader(FileName=files) 
    
    reader.PointArrayStatus = fields_to_export +  ["lsf", "inhibitor"]
    #reader.UpdatePipeline()
    reader.TimeArray = "None"

    #display = pv.Show(reader, renderView)
    #pv.ColorBy(display, ('POINTS', variable))
    #display.Representation = 'Surface'
    
    tmp_display = pv.Show(reader, renderView)
    pv.Render()

    camara = Camara_position(reader,renderView) ## center the image
    pv.Hide(reader, renderView)
    pv.Render()

    #save pngs:
    out_directory = os.path.join(vtu_file, "png")
    os.makedirs(out_directory, exist_ok=True)

    contour = pv.Contour(Input=reader)
    contour.ContourBy = ['POINTS', 'lsf']
    contour.Isosurfaces = [0.0]
    display_lsf = pv.Show(contour, renderView)
    #pv.ColorBy(display_lsf, None)
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

    
    for comp in fields_to_export:
        export_component_pngs(variable=comp,reader=reader,renderView=renderView,display_resample=display_resample,base_output_dir=out_directory,min_range=0.0,max_range=max_ranges[comp])

    pv.Hide(resample, renderView)
    
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

    export_component_pngs(variable="inhibitor",reader=reader,renderView=renderView,display_resample=display_resample_inhib,base_output_dir=out_directory,min_range=0.0,max_range=max_ranges["inhibitor"])