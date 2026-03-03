

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
    
############
import paraview.simple as pv
import glob
import os


def export_one(path,variable_input='concentration'):
    print(f"🚀 Exporting PNGs from: {path}")
    export_png_from_vtu(path, variable=variable_input, min_range=0.3, max_range=0.9)
    print(f"✅ Done: {path}")
    
        
def export_png_from_vtu(vtu_file, variable = 'concentration',min_range=0, max_range=1.0):
    pv.ResetSession()
    renderView = Create_view_white()

   #look for the pvtu files
    if os.path.isdir(vtu_file):
        files = sorted(glob.glob(os.path.join(vtu_file, "Results_*.pvtu")))
    else:
        files = sorted(glob.glob(vtu_file))
        
    # Leer archivo VTU
    reader = pv.XMLPartitionedUnstructuredGridReader(FileName=files)
    reader.PointArrayStatus = [variable]
    #reader.UpdatePipeline()
    reader.TimeArray = "None"

    display = pv.Show(reader, renderView)
    pv.ColorBy(display, ('POINTS', variable))
    display.Representation = 'Surface'

    lut = pv.GetColorTransferFunction(variable, display)
    lut.ApplyPreset('Inferno (matplotlib)', True)
    lut.InvertTransferFunction()
    lut.RescaleTransferFunction(min_range, max_range)
    
    camara = Camara_position(reader,renderView) ## center the image
    pv.Render()

    #save pngs:
    out_directory = os.path.join(vtu_file, "png")
    os.makedirs(out_directory, exist_ok=True)

    """     #Aditional for contours of level set function in 0 -------------------------
    #look for the pvtu files
    if os.path.isdir(vtu_file):
        files = sorted(glob.glob(os.path.join(vtu_file, "LSF_t*.pvtu")))
    else:
        files = sorted(glob.glob(vtu_file))
        
    # Leer archivo VTU
    lsf = pv.XMLPartitionedUnstructuredGridReader(FileName=files)
    lsf.PointArrayStatus = ["ca_cyt"]
    lsf.TimeArray = "None"
    
    contour = pv.Contour(Input=lsf)
    contour.ContourBy = ['POINTS', 'ca_cyt']
    contour.Isosurfaces = [0.0]
    display_lsf = pv.Show(contour, renderView)
    pv.ColorBy(display_lsf, None)
    display_lsf.Representation = 'Surface'
    display_lsf.DiffuseColor = [0, 0, 0]
    display_lsf.LineWidth = 2
    pv.Hide(display_lsf, renderView) """

    # Guardar PNGs
    for idx, t in enumerate(reader.TimestepValues):
        print (f" Saving timestep {idx} / {len(reader.TimestepValues)}")
        renderView.ViewTime = t
        pv.Render()
        pv.SaveScreenshot(os.path.join(out_directory, f"{idx:03d}.png"), 
                        renderView,
                        #quality=100,
                        ImageResolution=renderView.ViewSize, 
                        TransparentBackground=1)
