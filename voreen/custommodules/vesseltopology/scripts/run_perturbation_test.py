import voreen
import voreenqt
import time
import math
import os
import subprocess
import time

# --- configuration ---
iterations = 100
output_folder_base = "~/nosnapshot/perturbation_results_bgm"

def update_workspace_and_gui():
    voreen.invalidateProcessors()
    voreenqt.processEvents()


def run_test(parameters, output_folder):
    (perturbation_method, ) = parameters #unpack parameters

    output_file_name = "perturbation_{}.csv".format(perturbation_method)
    output_file_path = os.path.join(output_folder, output_file_name)

    voreen.setPropertyValue("VesselGraphCreator", "synchronousComputation", False)
    voreen.setPropertyValue("VesselGraphCreator", "continuousUpdate_", False)
    voreen.setPropertyValue("VesselGraphComparison", "enabled", False)

    update_workspace_and_gui()

    voreen.setPropertyValue("VesselGraphComparison", "statExportFile", output_file_path)
    if os.path.isfile(output_file_path):
        os.remove(output_file_path)

    update_workspace_and_gui()

    voreen.setPropertyValue("VesselGraphPerturbation", "perturbationMethod", perturbation_method)

    voreen.setPropertyValue("VesselGraphCreator", "synchronousComputation", True)
    voreen.setPropertyValue("VesselGraphCreator", "continuousUpdate_", True)
    voreen.setPropertyValue("VesselGraphComparison", "enabled", True)

    start = time.time()
    for i in range(iterations):
        rotation = i*2*math.pi/iterations
        elapsed = time.time() - start
        approx_total = 0
        if i > 0:
            approx_total = elapsed*iterations/i
        approx_remaining = approx_total - elapsed
        perturbation_amount = float(i)/(iterations-1);
        print("[{}/{}, {:5.1f}s / {:5.1f}s -- {:5.1f}s] Starting perturbation={:7.5f} for conf {}".format(i, iterations, elapsed, approx_total, approx_remaining, perturbation_amount, parameters))
        voreen.setPropertyValue("VesselGraphPerturbation", "perturbationAmount", perturbation_amount)

        voreen.repaint()
        voreenqt.processEvents()
    end = time.time()

    voreenqt.processEvents()
    voreen.setPropertyValue("VesselGraphComparison", "statExportFile", "")

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

    configurations = [
            ("add_edges",),
            ("split_nodes",),
            ("subdivide_edges",),
            ("split_edges",),
            ("move_nodes",),
            ("change_properties",),
            ]

    for configuration in configurations:
        run_test(configuration, output_folder)

run_tests()
