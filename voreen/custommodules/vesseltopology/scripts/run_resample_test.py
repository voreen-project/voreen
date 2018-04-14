import voreen
import voreenqt
import time
import math
import os
import subprocess

# --- configuration ---
iterations = 36
output_folder_base = "~/nosnapshot/resample_results"

def update_workspace_and_gui():
    voreen.invalidateProcessors()
    voreenqt.processEvents()

def run_test(parameters, output_folder):
    voreen.setPropertyValue("SegmentationValidation", "autoExport", False)

    asyncprocessors = ["Resample", "ResampleAndTransform", "ResampleAndTransformInverse"]
    for p in asyncprocessors:
        voreen.setPropertyValue(p, "synchronousComputation", True)

    update_workspace_and_gui()


    update_workspace_and_gui()

    start = time.time()
    for i in range(iterations):
        rotation = i*2*math.pi/iterations
        elapsed = time.time() - start
        approx_total = 0
        if i > 0:
            approx_total = elapsed*iterations/i
        approx_remaining = approx_total - elapsed
        print("[{}/{}, {:5.1f}s / {:5.1f}s -- {:5.1f}s] Starting rotation={:7.5f}".format(i, iterations, elapsed, approx_total, approx_remaining, rotation))

        output_file_path = output_folder + "/rot{:7.5f}.csv".format(rotation)
        if os.path.isfile(output_file_path):
            os.remove(output_file_path)

        voreen.setPropertyValue("RotationParameters", "rotationAngle", rotation)
        voreen.setPropertyValue("VolumeInputSwitch", "inputselect", 1)
        voreen.setPropertyValue("SegmentationValidation", "exportfile", output_file_path)
        voreen.setPropertyValue("SegmentationValidation", "autoExport", True)

        voreen.repaint()
        voreenqt.processEvents()
    end = time.time()

    voreenqt.processEvents()
    voreen.setPropertyValue("SegmentationValidation", "exportfile", "")
    voreen.setPropertyValue("SegmentationValidation", "autoExport", False)
    voreen.setPropertyValue("VolumeInputSwitch", "inputselect", 2)

    for p in asyncprocessors:
        voreen.setPropertyValue(p, "synchronousComputation", False)

    print("Finished conf {} in {:5.1f}s\n".format(parameters, end-start))

def run_tests():
    voreen.repaint()

    # make sure all Qt events have been processed before starting
    voreenqt.processEvents()


    dataset_url = voreen.getPropertyValue("VolumeSource", "volumeURL")
    parameter_index = dataset_url.find("?")

    if parameter_index != -1:
        url_without_parameters = dataset_url[:parameter_index]
    else:
        url_without_parameters = dataset_url

    input_file_identifier = os.path.basename(url_without_parameters)

    commit = subprocess.check_output(["git", "describe", "--always"]).rstrip()
    output_name = input_file_identifier + ":" + commit
    output_folder = os.path.expanduser(os.path.join(output_folder_base, output_name))

    if not os.path.isdir(output_folder):
        print("Creating output folder: {}".format(output_folder))
        os.makedirs(output_folder)

    voreen.setPropertyValue("RotationParameters", "transformationType", "rotation")

    configurations = [
            (1, 2, 3, 4),
            ]

    for configuration in configurations:
        run_test(configuration, output_folder)

run_tests()
