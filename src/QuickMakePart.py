import sys
import numpy as np
import os

try:
    runMode = os.environ["runMode"]
except:
    runMode = 'local'
    os.environ["runMode"] = runMode

if (runMode == 'docker'):
    FreeCADPath = '/usr/lib/freecad-daily/lib'
else:
    FreeCADPath = '/usr/lib/freecad-python3/lib'

sys.path.append(FreeCADPath)
print(sys.path)

import FreeCAD
import Part

# Create a new document
doc = FreeCAD.newDocument()

# Define the box dimensions
length = 10  # Length of the box in mm
width = 10   # Width of the box in mm
height = 10  # Height of the box in mm

# Create a box shape
box_shape = Part.makeBox(length, width, height)

# Create a solid object from the shape
box_solid = Part.Solid(box_shape)

# Add the solid to the document
box_solid = doc.addObject("Part::Feature", "Box").Shape 

# Export the document to an STP file
output_path = "box.stp"
doc.exportStep(output_path)

# Print the path to the exported STP file
print(f"Box exported to: {output_path}")

