import os 
import re
import shutil
import voreen
import voreenqt

# read original ensemble path
simulation_ensemble_path = voreen.getPropertyValue("EnsembleDataSource", "ensemblepath")

# name of the new subolder
name = "use_case_c_ensemble"

def copy_files(simulation_ensemble_path, target_ensemble_path):
    # regex-pattern
    pattern = re.compile(r"results_iT\d+iC\d+\.vti$")
    # loop for all folders
    for main_folder in os.listdir(simulation_ensemble_path):
        main_path = os.path.join(simulation_ensemble_path, main_folder)
        vtk_data_path = os.path.join(main_path, "vtkData", "data")

        if os.path.isdir(vtk_data_path):
            target_folder = os.path.join(target_ensemble_path, main_folder)
            os.makedirs(target_folder, exist_ok=True)

            for filename in os.listdir(vtk_data_path):
                if pattern.match(filename):
                    src_file = os.path.join(vtk_data_path, filename)
                    dest_file = os.path.join(target_folder, filename)
                    shutil.copy2(src_file, dest_file)


# prevent recursive nesting
if name in os.path.basename(simulation_ensemble_path).lower():
    voreenqt.messageBox("It looks as if the ensemble has already been created, since the target path already exists. Skipping.")
else:
    # define shifted path
    target_ensemble_path = os.path.join(simulation_ensemble_path, name)
    os.makedirs(target_ensemble_path, exist_ok=True)

    # do the copy
    copy_files(simulation_ensemble_path, target_ensemble_path)
    
    # update ensemble path in voreen
    voreen.setPropertyValue("EnsembleDataSource", "ensemblepath", target_ensemble_path)

    #print(f"Copied files from {simulation_ensemble_path} to {target_ensemble_path}")
    voreenqt.messageBox("File transfer is done. You can load the ensemble now")

