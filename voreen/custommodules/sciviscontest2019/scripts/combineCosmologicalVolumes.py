# Voreen Python script
import voreen
import voreenqt

numberOfFiles = 625

pathToFolder = '/home/m_stei44/Desktop/SciVis2019/data/64entropy/E1/'
for i in range(0, numberOfFiles):
   timestep = i
   index = str(i).zfill(len(str(numberOfFiles)))
   filename = "SPH_Internal_Energy_t_{}.h5".format(index)
   voreen.setPropertyValue("VolumeSelector", "volumeID", i)
   voreen.setPropertyValue("VolumeSelector 2", "volumeID", i)
   voreen.setPropertyValue("VolumeSave", "outputFilename", pathToFolder + filename)
   voreen.render();
