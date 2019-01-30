import voreen
import voreenqt
import os
import time

startTime = time.time()
inputPath = "/home/spider/Desktop/SciVis/input/"
outputPath = "/home/spider/Desktop/SciVis/output/"

files = []
for file in os.listdir(inputPath):
    if file.endswith(".vti"):
        files.append(file)

fileCounter = 1

for vtiFile in files:
    subStartTime = time.time()

    fullInputPath = inputPath + vtiFile
    fullOutputPath = outputPath + vtiFile.replace(".vti", ".hdf5")

    print ("Converting %s to %s (%s of %s)" % fullInputPath, fullOutputPath, fileCounter, len(files))
    fileCounter += 1

    if os.path.isfile(fullOutputPath):
        print ("File %s already exists. Continuing" % (fullOutputPath))
        voreenqt.processEvents()
        continue

    voreenqt.processEvents()

    voreen.setPropertyValue("VolumeListSave", "fileNameHDF5", fullOutputPath)
    voreen.loadVolumes(fullInputPath, True, True, "VolumeListSource")
    
    print ("Converting file took %s seconds" % (time.time() - subStartTime))
    voreenqt.processEvents()

print ("The whole dataset was converted in %s seconds" % (time.time() - startTime))
