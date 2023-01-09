import voreen
import os
import time

startTime = time.time()
inputPath = "/DATA/groupdrive/scivis_contest/2022/vts/"
#outputPath = "/DATA/groupdrive/scivis_contest/2022/H5/"
outputPath = "/DATA/owncloud/shares/SciVis Contest 22/ensemble"

def convert_directory(directory):

    files = []
    for file in os.listdir(os.path.join(inputPath, directory)):
        if file.endswith(".vts"):
            files.append(file)

    files.sort()

    fileCounter = 1

    for vtiFile in files:
        subStartTime = time.time()

        fullInputPath = os.path.join(inputPath, directory, vtiFile)
        fullOutputPath = os.path.join(outputPath, directory)
        os.makedirs(fullOutputPath, exist_ok = True)
        fullOutputPath = os.path.join(fullOutputPath, vtiFile.replace(".vts", ".h5"))

        print ("Converting %s to %s (%s of %s)" % (fullInputPath, fullOutputPath, fileCounter, len(files)))
        fileCounter += 1

        if os.path.isfile(fullOutputPath):
            print ("File %s already exists. Continuing" % (fullOutputPath))
            continue

        voreen.loadVolumes(fullInputPath, True, True)
        voreen.render()
        voreen.setPropertyValue("VolumeListSave", "fileNameHDF5", fullOutputPath)
        voreen.setPropertyValue("VolumeListSave", "save", 1)
        voreen.render()
        print ("Converting file took %s seconds" % (time.time() - subStartTime))

    print ("The whole run was converted in %s seconds" % (time.time() - startTime))

directories = []
for directory in os.listdir(inputPath):
    if os.path.isdir(os.path.join(inputPath, directory)):
        convert_directory(directory)
