import voreen
import voreenqt
import os

inputs = ["velocity"]
measures = ["magnitude", "vorticity"]
inputPath = "/mnt/data_drive/data/FlowEnsembles"

for property in os.listdir(inputPath):
    if property in inputs:
        for ensemble in os.listdir(os.path.join(inputPath, property)):
            for run in os.listdir(os.path.join(inputPath, property, ensemble)):
                for it in os.listdir(os.path.join(inputPath, property, ensemble, run)):
                    if it.endswith(".vvd"):
                        file = os.path.join(inputPath, property, ensemble, run, it)
                        voreen.loadVolume(file, property)
                        for measure in measures:
                            outputPath = os.path.join(inputPath, measure, ensemble, run)
                            os.makedirs(outputPath, exist_ok=True)
                            outputPath = os.path.join(outputPath, it.replace(property, "t"))
                            print(outputPath)
                            voreen.setPropertyValue(measure, "outputFilename", outputPath)
                            voreen.render()

                print("Finished %s/%s" % (ensemble, run))
                voreenqt.processEvents()
