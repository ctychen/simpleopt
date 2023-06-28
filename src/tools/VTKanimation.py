from paraview.simple import *

# disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()

# set the number of frames
num_frames = 100 # adjust based on your number of VTK files

# setup view
view1 = GetActiveViewOrCreate('RenderView')
view1.ViewSize = [800, 600] # adjust based on your preferred frame size
view1.InteractionMode = '2D'
view1.CameraPosition = [1.0, -1.0, 1.0]  # isometric view
view1.CameraFocalPoint = [0.0, 0.0, 0.0]
view1.CameraViewUp = [0.0, 1.0, 0.0]  # y-axis is up

# iterate over each VTK file
for i in range(num_frames):
    # create a new 'Legacy VTK Reader'
    vtk_file = f'/Users/cchen/Desktop/simpleopt/src/test0011_08/{i:03}.vtk' # adjust path to your VTK files
    legacyVTKReader1 = LegacyVTKReader(FileNames=[vtk_file])

    # Render
    Render()
    # reset camera to show all data
    ResetCamera()
    # save frame as PNG
    WriteImage(f'/Users/cchen/Desktop/simpleopt/src/test0011_08/images_all/frame_{i:04}.png')

# # iterate over each VTK file
# for i in range(num_frames):
#     # create a new 'Legacy VTK Reader'
#     vtk_file = f'/Users/cchen/Desktop/simpleopt/src/test0011_08/{i:03}.vtk' # adjust path to your VTK files
#     legacyVTKReader1 = LegacyVTKReader(FileNames=[vtk_file])

#     # Render
#     Render()
#     # save frame as PNG
#     WriteImage(f'/Users/cchen/Desktop/simpleopt/src/test0011_08/images_all/frame_{i:04}.png')
