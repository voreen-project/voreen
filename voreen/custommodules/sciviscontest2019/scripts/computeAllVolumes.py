# Voreen Python script
import voreen
import voreenqt
import os

numberOfFiles = 265
resolution = 64

voreen.setPropertyValue("CosmologyVolumeConverter", "volumeDimensions", [resolution, resolution, resolution])
pathToBaseFolder = '/home/m_stei44/Desktop/SciVis2019/data/'

nV_properties = ["mass", "phi", "velocity", "velocityx", "velocityy", "velocityz"]
#SPH_properties = ["uu", "temperature", "sphDensity"]
#density_types = ["baryon", "dark matter", "wind", "star", "gas", "agn", "all"]


#####################
#  nearestVoxel
#####################
voreen.setPropertyValue("CosmologyVolumeConverter", "spreadMode", "nearestVoxel")
nv_path = pathToBaseFolder + "nearestVoxel/"
os.mkdir(nv_path)

for prop in nV_properties:
   voreen.setPropertyValue("CosmologyVolumeConverter", "particleProperty", prop) 
   nv_prop_path = nv_path + prop + "/"
   os.mkdir(nv_prop_path)
   for i in range(0, numberOfFiles):
      print(i)
      timestep = i
      voreen.setPropertyValue("CosmologyVolumeConverter", "timeStep", timestep)
      index = str(i).zfill(len(str(numberOfFiles)))
      filename = "nV_{}_t_{}.h5".format(prop, index)
      voreen.setPropertyValue("VolumeSave", "outputFilename", nv_prop_path + filename)
      voreen.render();

#####################
#  Densities
#####################
#voreen.setPropertyValue("CosmologyVolumeConverter", "spreadMode", "AMOUNT")
#densities_path = pathToBaseFolder + "Densities/"
#os.mkdir(densities_path)

#for type in density_types:
#   voreen.setPropertyValue("CosmologyVolumeConverter", "particleType", type) 
#   densities_type_path = densities_path + type + "/"
#   os.mkdir(densities_type_path)
#   for i in range(0, numberOfFiles):
#      timestep = i
#      voreen.setPropertyValue("CosmologyVolumeConverter", "timeStep", timestep)
#      index = str(i).zfill(len(str(numberOfFiles)))
#      filename = "density_{}_t_{}.h5".format(type, index)
#      voreen.setPropertyValue("VolumeSave", "outputFilename", densities_type_path + filename)
#      voreen.render();

#####################
#  SPH
#####################
#voreen.setPropertyValue("CosmologyVolumeConverter", "spreadMode", "SPH")
#SPH_path = pathToBaseFolder + "SPH/"
#os.mkdir(SPH_path)

#for prop in SPH_properties:
#   voreen.setPropertyValue("CosmologyVolumeConverter", "particleProperty", prop) 
#   SPH_prop_path = SPH_path + prop + "/"
#   os.mkdir(SPH_prop_path)
#
#   for i in range(0, numberOfFiles):
#      timestep = i
#      voreen.setPropertyValue("CosmologyVolumeConverter", "timeStep", timestep)
#      index = str(i).zfill(len(str(numberOfFiles)))
#      filename = "SPH_{}_t_{}.h5".format(prop, index)
#      voreen.setPropertyValue("VolumeSave", "outputFilename", SPH_prop_path + filename)
#      voreen.render();



