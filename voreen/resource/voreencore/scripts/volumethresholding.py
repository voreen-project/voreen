# Voreen Python script
import voreen

# Fetch port data.
volume = voreen.getPortData("DynamicPythonProcessor", "VolumePort-in (1)")
print ("Format: %s " % (volume.format))
print ("Dimensions: (%i, %i, %i)" % (volume.dimX, volume.dimY, volume.dimZ))
print ("Spacing: (%f, %f, %f)" % (volume.spacingX, volume.spacingY, volume.spacingZ))
print ("Offset: (%f, %f, %f)" % (volume.offsetX, volume.offsetY, volume.offsetZ))
print ("RealWorldMapping (Offset, Scale): (%f, %f)" % (volume.rwmOffset, volume.rwmScale))

# Perform thresholding
lowerThreshold = 0.0
upperThreshold = 0.7

for i in range(0, len(volume.data)):
    if (volume.data[i] < lowerThreshold):
        volume.data[i] = lowerThreshold
    if (volume.data[i] > upperThreshold):
        volume.data[i] = upperThreshold

# Set output format.
volume.format = "float"

# Set data to port.
voreen.setPortData("DynamicPythonProcessor", "VolumePort-out (2)", volume);