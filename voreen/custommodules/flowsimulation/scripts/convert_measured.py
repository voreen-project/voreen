# This script assumes a workspace containing
#   * PBReader
#   * VolumeListSource named magnitude_save connected to PBReader magnitude output
#   * VolumeListSource named velocity_save connected to PBReader velocity output

import voreen
import voreenqt
import os

# Name of the file indicating a valid directory.
parameterFile = "FOV_parameters.txt"

# Fields the PBReader outputs.
fields = ["magnitude", "velocity"]

# Input path, folder containing series of measurements.
inputPath = "/home/spider/ownCloud/data/4D-PC-MRI/RAW/"

# Output path, containing a VVD file for each field specified above.
outputPath = "/home/spider/ownCloud/data/4D-PC-MRI/converted/"

# Determine if all converted runs should be stored in a flattened hierarchy.
flatten = True

for root, subdirs, files in os.walk(inputPath):
    if parameterFile in files:
        # Load dataset.
        voreen.setPropertyValue("PBReader", "folderProp", root)
        voreen.setPropertyValue("PBReader", "loadButtonProp", 1)
        # Store it as VVD file.
        for field in fields:
            processorName = field + "_save"
            voreen.setPropertyValue(processorName, "baseName", field + "_")
            fieldOutputPath = os.path.relpath(root, inputPath)
            if flatten:
                fieldOutputPath = "-".join(os.path.normpath(fieldOutputPath).split(os.sep))
            fieldOutputPath = os.path.join(outputPath, fieldOutputPath)
            os.makedirs(fieldOutputPath,  exist_ok=True)
            voreen.setPropertyValue(processorName, "folderNameVVD", fieldOutputPath)
            voreen.setPropertyValue(processorName, "save", 1)

        # Print progress.
        print("Finished %s" % (root))
        voreenqt.processEvents()
