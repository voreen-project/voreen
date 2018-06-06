import voreen
import voreenqt
import time
import math
import os
import subprocess

# --- configuration ---
iterations = 10
total_shift_amount = 0.01625
#total_shift_amount = 1.0
shift_direction = (1.0,1.0,1.0)
output_folder_base = "~/nosnapshot/GERoMe_shift_results"

def update_workspace_and_gui():
    voreen.invalidateProcessors()
    voreenqt.processEvents()

def run_test(parameters, output_folder):
    num_refinement_iterations, min_voxel_length, min_elongation, min_bulge = parameters #unpack parameters

    output_file_name = "it{}_len{:0>2}_elo{:.2f}_bul{:.2f}.csv".format(num_refinement_iterations, min_voxel_length, min_elongation, min_bulge)
    output_file_path = os.path.join(output_folder, output_file_name)

    voreen.setPropertyValue("VesselGraphComparison", "enabled", False)

    update_workspace_and_gui()

    voreen.setPropertyValue("VesselGraphComparison", "statExportFile", output_file_path)
    if os.path.isfile(output_file_path):
        os.remove(output_file_path)

    update_workspace_and_gui()
    voreen.setPropertyValue("VesselGraphComparison", "enabled", True)

    voreen.setPropertyValue("VesselGraphCreator", "numRefinementIterations", num_refinement_iterations)
    voreen.setPropertyValue("VesselGraphCreator", "minVoxelLength", min_voxel_length)
    voreen.setPropertyValue("VesselGraphCreator", "minElongation", min_elongation)
    voreen.setPropertyValue("VesselGraphCreator", "minBulgeSize", min_bulge)

    start = time.time()
    for i in range(iterations):
        translation_amount = i*total_shift_amount/iterations
        translation_vector = (shift_direction[0]*translation_amount, shift_direction[1]*translation_amount, shift_direction[2]*translation_amount)
        elapsed = time.time() - start
        approx_total = 0
        if i > 0:
            approx_total = elapsed*iterations/i
        approx_remaining = approx_total - elapsed
        print("[{}/{}, {:5.1f}s / {:5.1f}s -- {:5.1f}s] Starting translation={:7.5f} for conf {}".format(i, iterations, elapsed, approx_total, approx_remaining, translation_amount, parameters))
        voreen.setPropertyValue("RotationParameters", "translationAmount", translation_vector)
        voreen.setPropertyValue("VesselGraphComparison", "datasetIdentifier", "{}".format(translation_amount))
        voreen.setPropertyValue("VolumeInputSwitch", "inputselect", 1)

        voreen.repaint()
        voreenqt.processEvents()
    end = time.time()

    voreenqt.processEvents()
    voreen.setPropertyValue("VesselGraphComparison", "statExportFile", "")
    voreen.setPropertyValue("VolumeInputSwitch", "inputselect", 2)

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

    # Yeah, this does not make sense semantically, but we don't want to create another workspace...
    voreen.setPropertyValue("RotationParameters", "transformationType", "translation")

    configurations = [
            #(0,0,0,0),
            #(1,0,0,1.0),
            #(2,0,0,1.0),
            #(3,0,0,1.0),
            #(2,0,0,1.1),
            #(1,0,0,1.5),
            #(2,0,0,1.5),
            #(3,0,0,1.5),
            #(1,0,0,2.0),
            (2,0,0,2.0),
            #(3,0,0,2.0),
            #(1,15,0,0),
            #(2,15,0,0),
            #(1,0,1,0),
            #(2,0,1,0),
            ]

    for configuration in configurations:
        run_test(configuration, output_folder)

# Start
#if voreenqt.questionBox(
#        "Welcome to benchmark.py!\n" +
#        "Start?"):
#else:
#    print "aborting on user request"
run_tests()
