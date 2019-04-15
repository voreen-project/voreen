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
inputPath = "/home/spider/ownCloud/data/RAW/Aneurysma/Aneurysma I/3D-Druck"

# Output path.
outputPath = "/mnt/data_drive/data/FlowEnsembles"

for ensemble in os.listdir(inputPath):
    for run in os.listdir(os.path.join(inputPath, ensemble)):
        runPath = os.path.join(inputPath, ensemble, run)
        if parameterFile in os.listdir(runPath):
            # Load dataset.
            voreen.setPropertyValue("PBReader", "folderProp", runPath)
            voreen.setPropertyValue("PBReader", "loadButtonProp", 1)
            # Store it as VVD file.
            for field in fields:
                processorName = field + "_save"
                voreen.setPropertyValue(processorName, "baseName", field + "_")
                fieldOutputPath = os.path.join(outputPath, field, ensemble, run)
                os.makedirs(fieldOutputPath,  exist_ok=True)
                voreen.setPropertyValue(processorName, "folderNameVVD", fieldOutputPath)
                voreen.setPropertyValue(processorName, "save", 1)

            # Print progress.
            print("Finished %s/%s" % (ensemble, run))
            voreenqt.processEvents()
