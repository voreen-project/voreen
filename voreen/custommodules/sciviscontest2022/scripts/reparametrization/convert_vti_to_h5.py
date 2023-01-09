import voreen
import voreenqt
import os
import time

startTime = time.time()
inputPath = "/media/m_ever14/Seagate Expansion Drive/SciVis22/flattenVTI/"
outputPath = "/media/m_ever14/Seagate Expansion Drive/SciVis22/HDF5/"

fields = ["O2", "convht_1", "frhosiesrad_1", "rhof_1", "rhowatervapor", "theta", "u", "v", "w"]

def convertRun(inPath, outPath):
    #inPath = inputPath#inputPath+inPath+"/"
    #outPath = outputPath#outputPath+outPath+"/"
    files = []
    for file in os.listdir(inPath):
        if file.endswith(".vti"):
            files.append(file)

    fileCounter = 1

    for vtiFile in files:
        subStartTime = time.time()

        fullInputPath = os.path.join(inPath, vtiFile)
        fullOutputPath = os.path.join(outPath, vtiFile.replace(".vti", ".hdf5"))

        print ("Converting %s to %s (%s of %s)" % (fullInputPath, fullOutputPath, fileCounter, len(files)))
        fileCounter += 1

        if os.path.isfile(fullOutputPath):
            print ("File %s already exists. Continuing" % (fullOutputPath))
            voreenqt.processEvents()
            continue

        voreenqt.processEvents()

        firstField = True
        for field in fields:
            if(firstField):
                voreen.loadVolumes(fullInputPath+"?name="+field, True, True, "VolumeListSource")
                firstField = False
            voreen.loadVolumes(fullInputPath+"?name="+field, True, False, "VolumeListSource")
        voreen.setPropertyValue("VolumeListSave", "fileNameHDF5", fullOutputPath)
    
        print ("Converting file took %s seconds" % (time.time() - subStartTime))
        voreenqt.processEvents()

    print ("The whole run was converted after %s seconds" % (time.time() - startTime)) 

runs = []
for run in os.listdir(inputPath):
    if(not os.path.exists(os.path.join(outputPath, run))):
        #print("mkdir \""+os.path.join(outputPath, run)+"\"")
        os.system("mkdir \""+os.path.join(outputPath, run)+"\"")
    convertRun(os.path.join(inputPath, run), os.path.join(outputPath, run))
