from paraview.simple import *

# disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()

# set the number of frames
num_frames = 167 # adjust based on number of VTK files

# setup view
view1 = GetActiveViewOrCreate('RenderView')
view1.ViewSize = [2000, 1500] 
view1.InteractionMode = '2D'
view1.CameraPosition = [0.2, 1.0, 0.2]  # isometric view, ish
view1.CameraFocalPoint = [0.0, 0.0, 0.0]
view1.CameraViewUp = [0.0, 1.0, 0.0]  # y-axis is up

# iterate over each VTK file
for i in range(num_frames):
    # create a new 'Legacy VTK Reader'
    vtk_file = f'/Users/cchen/Desktop/simpleopt/src/test0011_08/{i:03}.vtk' #path to VTK files
    legacyVTKReader1 = LegacyVTKReader(FileNames=[vtk_file])

    # Show data in view
    Show(legacyVTKReader1, view1)

    # Render
    Render()
    
    # reset camera to show all data
    ResetCamera()

    # save frame as PNG
    WriteImage(f'/Users/cchen/Desktop/simpleopt/src/test0011_08/images_all_2/frame_{i:04}.png')

    # Once done with a file, you can remove the data from memory using the Delete function
    Delete(legacyVTKReader1)
    del legacyVTKReader1
