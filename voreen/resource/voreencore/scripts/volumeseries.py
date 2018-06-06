# Consecutively loads the ten time steps of the golf ball data set
# (http://www.voreen.org/files/volumedynamics.zip), 
# and takes a snapshot of each frame.

import voreen
import voreenqt

golfBallPath = voreen.getBasePath() + "/data/volumes/golfball"
snapshotPath = voreen.getBasePath() + "/data/screenshots"

for i in range(0, 10):
    # usage: voreen.loadVolume(filepath, [VolumeSource processor]) 
    voreen.loadVolume(golfBallPath + "/golfball0_%d-512x256x256.dat" % (i))
    
    # alternative: use VolumeSeriesSource as data supplier
    #voreen.setPropertyValue("VolumeSeriesSource", "step", i) 

    voreen.repaint()
    voreenqt.processEvents() # to see each frame on screen
    voreen.snapshot(snapshotPath + "/golfball0_%d.png" % (i))

voreenqt.messageBox("Wrote ten frames to %s/golfball0_*.png" % snapshotPath)
