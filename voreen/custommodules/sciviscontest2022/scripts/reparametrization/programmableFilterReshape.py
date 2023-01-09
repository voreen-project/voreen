import numpy as np
import os

# Get Input
input = self.GetInput()
numPoints = input.GetNumberOfPoints()

# Create Output Dataset
output = self.GetOutput()
output.ShallowCopy(input)

minZ = {}

path = "/home/m_ever14/Desktop/minZ.npy"
if(not os.path.exists(path)):
    for i in range(0,numPoints):
        coord = input.GetPoint(i)
        x,y,z = coord[:3]
        if(int(x) in minZ):
            minZ[int(x)] = np.min([z, minZ[int(x)]])
        else:
            minZ[int(x)] = z

    np.save(path, minZ)
else:
    minZ = np.load(path, allow_pickle=True).item()
 
newPoints = vtk.vtkPoints()

for i in range(0,numPoints):
    coord = input.GetPoint(i)
    x,y,z = coord[:3]
    z = z-minZ[int(x)]
    newPoints.InsertPoint(i,x,y,z)
output.SetPoints(newPoints)
