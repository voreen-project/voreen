# Voreen Python script
import voreen

# Fetch port data.
volume = voreen.getPortData("DynamicPythonProcessor", "VolumePort Inport #1")
print ("Format: %s " % (volume.format))
print ("Dimensions: (%i, %i, %i)" % (volume.dimX, volume.dimY, volume.dimZ))
print ("Spacing: (%f, %f, %f)" % (volume.spacingX, volume.spacingY, volume.spacingZ))
print ("Offset: (%f, %f, %f)" % (volume.offsetX, volume.offsetY, volume.offsetZ))

# Perform thresholding
lowerThreshold = 0.0
upperThreshold = 1.0

for i in range(0, len(volume.data)):
    if (volume.data[i] < lowerThreshold):
        volume.data[i] = lowerThreshold
    if (volume.data[i] > upperThreshold):
        volume.data[i] = upperThreshold

# Set output format.
volume.format = "float"

# Set data to port.
voreen.setPortData("DynamicPythonProcessor", "VolumePort Outport #2", volume);