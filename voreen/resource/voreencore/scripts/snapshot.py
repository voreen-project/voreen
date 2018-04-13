# Sample script for taking canvas snapshots. 
# Writes a snapshot of each canvas to disc.
import voreen
import time

filename = voreen.getBasePath() + "/data/screenshots/snapshot-%d.png"

start = time.time()
count = voreen.canvasCount()
print "canvas count is %d:" % count
for i in range(0, count):
    voreen.snapshotCanvas(i, filename % (i+1))
    print " * saved snapshot to file %s" % (filename % (i+1))

end = time.time()
runtime = (end - start)
print 'total runtime: %fs' % runtime
