# Sample script for loading volume data.
import voreen

# usage: voreen.loadVolume(filepath, [name of VolumeSource processor])
voreen.loadVolume(voreen.getBasePath() + "/resource/voreenve/volumes/nucleon.vvd", "VolumeSource")
