# Sample script for loading a transfer function from disc.
import voreen
import voreenqt

# usage: 
# voreen.loadTransferFunction(processor name, 
#                             property id, 
#                             transfunc path)
voreen.loadTransferFunction(
    "SingleVolumeRaycaster", 
    "transferFunction", 
    voreen.getBasePath() + "/resource/transferfuncs/nucleon.tfi")
