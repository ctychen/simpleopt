import vtk
from vtk.util import numpy_support

for i in range(100): # replace 100 with the number of VTK files you have
    # Set file name
    filename = f"/Users/cchen/Desktop/simpleopt/src/test0011_08/{i:03}.vtk"  # adjust path to your VTK files

    # Read the source file.
    reader = vtk.vtkUnstructuredGridReader()
    reader.SetFileName(filename)
    reader.Update()  # Needed because of GetScalarRange
    output = reader.GetOutput()

    # Create the mapper that corresponds the objects of the vtk.vtk file
    # into graphics elements
    mapper = vtk.vtkDataSetMapper()
    mapper.SetInputData(output)

    # Create the Actor
    actor = vtk.vtkActor()
    actor.SetMapper(mapper)

    # Create the Renderer
    renderer = vtk.vtkRenderer()
    renderer.AddActor(actor)
    renderer.SetBackground(1, 1, 1)  # Set background to white

    # Create the RendererWindow
    renderer_window = vtk.vtkRenderWindow()
    renderer_window.AddRenderer(renderer)

    # Create the RendererWindowInteractor and display the vtk_file
    interactor = vtk.vtkRenderWindowInteractor()
    interactor.SetRenderWindow(renderer_window)

    # Setup camera
    camera = renderer.GetActiveCamera()
    camera.SetPosition(0, -1, 0)
    camera.SetFocalPoint(0, 0, 0)
    camera.SetViewUp(0, 0, 1)

    # Render
    interactor.Initialize()
    renderer_window.Render()

    # Save image
    window_to_image_filter = vtk.vtkWindowToImageFilter()
    window_to_image_filter.SetInput(renderer_window)
    window_to_image_filter.Update()

    writer = vtk.vtkPNGWriter()
    writer.SetFileName(f"/Users/cchen/Desktop/simpleopt/src/test0011_08/images_all/test{i:04}.png") # adjust path to save your frames
    writer.SetInputConnection(window_to_image_filter.GetOutputPort())
    writer.Write()

# Close window
interactor.TerminateApp()
