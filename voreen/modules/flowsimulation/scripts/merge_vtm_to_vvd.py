import voreen
import voreenqt
import os
import time

startTime = time.time()
inputPath = "/home/spider/Desktop/OpenLB/PALMAII/OpenLB-intel/simulations/aorta3d/vtkData/data/"
outputPath = "/home/spider/Desktop/OpenLB/PALMAII/OpenLB-intel/simulations/aorta3d/converted/"

files = []
for file in os.listdir(inputPath):
    if file.endswith(".vtm"):
        files.append(file)

fileCounter = 0

names = ["physVelocity", "physPressure"]

for vtmFile in files:
    subStartTime = time.time()
    fileCounter += 1

    fullInputPath = inputPath + vtmFile
    print ("Merging %s (%s of %s)" % (fullInputPath, fileCounter, len(files)))

    for name in names:

        fullOutputPath = outputPath + name + "_" + vtmFile.replace(".vtm", ".vvd")

        if os.path.isfile(fullOutputPath):
            print ("File %s already exists. Continuing" % (fullOutputPath))
            voreenqt.processEvents()
            continue

        voreenqt.processEvents()

        voreen.setPropertyValue("VolumeSave", "outputFilename", fullOutputPath)
        voreen.loadVolumes(fullInputPath + "?name=" + name, True, True, "VolumeListSource")

print ("The whole dataset was converted in %s seconds" % (time.time() - startTime))
