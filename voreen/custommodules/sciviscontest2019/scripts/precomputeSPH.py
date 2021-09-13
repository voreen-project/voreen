# Voreen Python script
import voreen
import voreenqt

numberOfFiles = 265

resolution = 16
voreen.setPropertyValue("CosmologyVolumeConverter", "spreadMode", "SPH")
voreen.setPropertyValue("CosmologyVolumeConverter", "particleProperty", "uu")
voreen.setPropertyValue("CosmologyVolumeConverter", "volumeDimensions", [resolution, resolution, resolution])
pathToFolder = '/home/m_stei44/Desktop/SciVis2019/data/uu16metaTS1/E1/'

for i in range(0, numberOfFiles):
   timestep = (1.0 * i)/numberOfFiles
   voreen.setPropertyValue("CosmologyVolumeConverter", "timeStep", timestep)
   index = str(i).zfill(len(str(numberOfFiles)))
   filename = "SPH_Internal_Energy_t_{}.h5".format(index)
   voreen.setPropertyValue("VolumeSave", "outputFilename", pathToFolder + filename)
   voreen.render();
