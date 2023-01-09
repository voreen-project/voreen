import os
import numpy as np
#### import the simple module from the paraview
from paraview.simple import *
#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()

ENSEMBLE_PATH = "/media/m_ever14/Seagate Expansion Drive/SciVis22/originalEnsemble"
DATA_PATH = "/media/m_ever14/Seagate Expansion Drive/SciVis22/savedData"
OUTPUT_PATH = "/media/m_ever14/Seagate Expansion Drive/SciVis22/flattenVTI"

def zMinSingleRun(path, ensembleName):
    # create a new 'XML Structured Grid Reader'
    output1000vts = XMLStructuredGridReader(FileName=[path+'/output.1000.vts'])

    # get active view
    renderView1 = GetActiveViewOrCreate('RenderView')
    # uncomment following to set a specific view size
    # renderView1.ViewSize = [1424, 937]

    # get layout
    #layout1 = GetLayout()

    # show data in view
    output1000vtsDisplay = Show(output1000vts, renderView1, 'StructuredGridRepresentation')

    # trace defaults for the display properties.
    output1000vtsDisplay.Representation = 'Outline'

    # reset view to fit data
    renderView1.ResetCamera()

    # update the view to ensure updated data information
    renderView1.Update()

    # create a new 'Programmable Filter'
    programmableFilter1 = ProgrammableFilter(Input=output1000vts)

    # Properties modified on programmableFilter1
    programmableFilter1.OutputDataSetType = 'vtkStructuredGrid'
    programmableFilter1.Script = """
    import numpy as np
    import os
    
    # Get Input
    input = self.GetInput()
    numPoints = input.GetNumberOfPoints()
    
    # Create Output Dataset
    output = self.GetOutput()
    output.ShallowCopy(input)
    
    minZ = {}
    
    path = """ + '"' + os.path.join(DATA_PATH, ensembleName+"_minZ.npy")+'"'+"""
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
        minZ = np.load(path, allow_pickle=True).item()"""
    programmableFilter1.RequestInformationScript = ''
    programmableFilter1.RequestUpdateExtentScript = ''
    programmableFilter1.PythonPath = ''

    # show data in view
    programmableFilter1Display = Show(programmableFilter1, renderView1, 'StructuredGridRepresentation')

    # trace defaults for the display properties.
    programmableFilter1Display.Representation = 'Outline'

    # hide data in view
    Hide(output1000vts, renderView1)

    # update the view to ensure updated data information
    renderView1.Update()

def flattenAndResample(path, ensembleName, filename, height):
    # create a new 'XML Structured Grid Reader'
    output1000vts = XMLStructuredGridReader(FileName=[os.path.join(path, filename)])

    # get active view
    renderView1 = GetActiveViewOrCreate('RenderView')
    # uncomment following to set a specific view size
    # renderView1.ViewSize = [1424, 937]

    # get layout
    #layout1 = GetLayout()

    # show data in view
    output1000vtsDisplay = Show(output1000vts, renderView1, 'StructuredGridRepresentation')

    # trace defaults for the display properties.
    output1000vtsDisplay.Representation = 'Outline'

    # reset view to fit data
    renderView1.ResetCamera()

    # update the view to ensure updated data information
    renderView1.Update()

    # create a new 'Programmable Filter'
    programmableFilter1 = ProgrammableFilter(Input=output1000vts)

    # Properties modified on programmableFilter1
    programmableFilter1.OutputDataSetType = 'vtkStructuredGrid'
    programmableFilter1.Script = """
    import numpy as np
    import os

    # Get Input
    input = self.GetInput()
    numPoints = input.GetNumberOfPoints()

    # Create Output Dataset
    output = self.GetOutput()
    output.ShallowCopy(input)

    minZ = {}

    path = """ + '"' + os.path.join(DATA_PATH, ensembleName+"_minZ.npy")+'"'+"""
    minZ = np.load(path, allow_pickle=True).item()

    newPoints = vtk.vtkPoints()

    for i in range(0,numPoints):
        coord = input.GetPoint(i)
        x,y,z = coord[:3]
        z = z-minZ[int(x)]
        newPoints.InsertPoint(i,x,y,z)
    output.SetPoints(newPoints)"""
    programmableFilter1.RequestInformationScript = ''
    programmableFilter1.RequestUpdateExtentScript = ''
    programmableFilter1.PythonPath = ''

    # show data in view
    programmableFilter1Display = Show(programmableFilter1, renderView1, 'StructuredGridRepresentation')

    # trace defaults for the display properties.
    programmableFilter1Display.Representation = 'Outline'

    # hide data in view
    Hide(output1000vts, renderView1)

    # update the view to ensure updated data information
    renderView1.Update()

    # create a new 'Resample To Image'
    resampleToImage1 = ResampleToImage(Input=programmableFilter1)

    # set active source
    SetActiveSource(programmableFilter1)

    # set active source
    SetActiveSource(resampleToImage1)

    # Properties modified on resampleToImage1
    resampleToImage1.UseInputBounds = 0
    resampleToImage1.SamplingDimensions = [600, 500, 50]
    resampleToImage1.SamplingBounds = [-498.0, 700.0, -500.0, 498.0, 0.0, height]

    # show data in view
    resampleToImage1Display = Show(resampleToImage1, renderView1, 'UniformGridRepresentation')

    # trace defaults for the display properties.
    resampleToImage1Display.Representation = 'Outline'

    # hide data in view
    Hide(programmableFilter1, renderView1)

    # update the view to ensure updated data information
    renderView1.Update()

    # save data
    #'/home/m_ever14/Desktop/resampledTest.vti'
    outpath = os.path.join(os.path.join(OUTPUT_PATH, ensembleName),filename).replace("vts", "vti")
    print(outpath)
    SaveData(outpath, proxy=resampleToImage1,
             PointDataArrays=['O2', 'convht_1', 'frhosiesrad_1', 'rhof_1', 'rhowatervapor', 'theta', 'u', 'v',
                              'vtkGhostType', 'vtkValidPointMask', 'w'],
             CellDataArrays=['vtkGhostType'])
    # destroy output1000vts
    Delete(resampleToImage1)
    del resampleToImage1
    Delete(programmableFilter1)
    del programmableFilter1
    Delete(output1000vts)
    del output1000vts

# Find minimal z values
#for member in sorted(os.listdir(ENSEMBLE_PATH)):
#    memberPath = os.path.join(ENSEMBLE_PATH, member)
#    zMinSingleRun(memberPath, member)
#    print("Member " + member + " done")

# Determine joint height
height = 890.0
for file in sorted(os.listdir(DATA_PATH)):
    zMins = np.load(os.path.join(DATA_PATH, file), allow_pickle = True).item()
    zMinMax = np.max(list(zMins.values()))
    height = min(height, 890.0-zMinMax*2)

# Flatten data and resample
#for member in sorted(os.listdir(ENSEMBLE_PATH)):
#    memberPath = os.path.join(ENSEMBLE_PATH, member)
#    if(not os.path.exists(os.path.join(OUTPUT_PATH, member))):
#        os.mkdir(os.path.join(OUTPUT_PATH, member))
#    for timestep in sorted(os.listdir(memberPath)):
#        flattenAndResample(memberPath, member, timestep, height)
#        print("Timestep " + timestep + " for " + member + " done")

# Only previously corrupted timestep
flattenAndResample(os.path.join(ENSEMBLE_PATH, "mountain_backcurve40"), "mountain_backcurve40", "output.2000.vts", height)
        
        
        
        
        
        
