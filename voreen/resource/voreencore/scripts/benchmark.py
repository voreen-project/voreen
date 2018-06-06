# Voreen Python script for measuring frame rates. 
import voreen
import voreenqt
import time
import math

# --- configuration ---
frames    = 100
canvasDim = (1024, 768)

processor  = "SingleVolumeRaycaster"
cameraProp = "camera"
canvas     = "Canvas"
initCam    = ((0.0, 0.0, 2.75), (0.0, 0.0, 0.0), (0.0, 1.0, 0.0))
rotAxis    = (0.0, 1.0, 0.0)
# --- configuration ---

def benchmark():
    # store current canvas size and camera settings
    prevDim = voreen.getPropertyValue(canvas, "canvasSize")
    prevCam = voreen.getPropertyValue(processor, cameraProp)

    # resize canvas and initialize camera
    voreen.setPropertyValue(canvas, "canvasSize", canvasDim)
    voreen.setPropertyValue(processor, cameraProp, initCam)
    #voreen.repaint() # not necessary any longer

    # make sure all Qt events have been processed before starting
    voreenqt.processEvents()  

    # start loop for animation
    counter = 0
    angleIncr = 2*math.pi / frames  # do a full rotation
    start = time.time()
    while (counter < frames):
        counter = counter + 1

        # update camera for current frame
        voreen.rotateCamera(processor, cameraProp, angleIncr, rotAxis)

        # render network state
        #voreen.repaint() # not necessary any longer
        voreenqt.processEvents() # this will trigger repaint.
        
    end = time.time()

    # restore canvas and camera
    voreen.setPropertyValue(canvas, "canvasSize", prevDim)
    voreen.setPropertyValue(processor, cameraProp, prevCam)
    #voreen.repaint() # not necessary any longer

    # evaluate results
    runtime = (end - start)
    if runtime > 0:
        fps = frames / runtime
    else:
        fps = -1

    message =  'frames: %i\n' % (frames)
    message += 'total runtime: %f\n' % (runtime)
    message += 'fps: %f' %(fps)
    print message
    voreenqt.messageBox("Finished!\n" + message)


# Start
if voreenqt.questionBox(
        "Welcome to benchmark.py!\n" + 
        "Will run for %i frames with %ix%i pixels.\n" % 
        (frames, canvasDim[0], canvasDim[1]) +
        "Start?"):
    benchmark()
else:
    print "aborting on user request"
