/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2018 University of Muenster, Germany,                        *
 * Department of Computer Science.                                                 *
 * For a list of authors please refer to the file "CREDITS.txt".                   *
 *                                                                                 *
 * This file is part of the Voreen software package. Voreen is free software:      *
 * you can redistribute it and/or modify it under the terms of the GNU General     *
 * Public License version 2 as published by the Free Software Foundation.          *
 *                                                                                 *
 * Voreen is distributed in the hope that it will be useful, but WITHOUT ANY       *
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR   *
 * A PARTICULAR PURPOSE. See the GNU General Public License for more details.      *
 *                                                                                 *
 * You should have received a copy of the GNU General Public License in the file   *
 * "LICENSE.txt" along with this file. If not, see <http://www.gnu.org/licenses/>. *
 *                                                                                 *
 * For non-commercial academic use see the license exception specified in the file *
 * "LICENSE-academic.txt". To get information about commercial licensing please    *
 * contact the authors.                                                            *
 *                                                                                 *
 ***********************************************************************************/

#include "flowsimulationcluster.h"

#include "voreen/core/datastructures/geometry/glmeshgeometry.h"
#include "voreen/core/datastructures/volume/volumeatomic.h"
#include "voreen/core/ports/conditions/portconditionvolumelist.h"

#include "modules/core/io/rawvolumereader.h"

namespace voreen {

const std::string FlowSimulationCluster::loggerCat_("voreen.flowreen.FlowSimulationCluster");

FlowSimulationCluster::FlowSimulationCluster()
    : Processor()
    // ports
    , geometryDataPort_(Port::INPORT, "geometryDataPort", "Geometry Input", false)
    , measuredDataPort_(Port::INPORT, "measuredDataPort", "Measured Data Input", false)
    , parameterPort_(Port::INPORT, "parameterPort", "Parameterization", false)
    , useLocalInstance_("useLocalInstance", "Use local Instance", false)
    , localInstancePath_("localInstancePath", "Local Instance Path", "Path", "", "EXE (*.exe)", FileDialogProperty::OPEN_FILE, Processor::VALID, Property::LOD_DEFAULT, VoreenFileWatchListener::ALWAYS_OFF)
    , institution_("institution", "Institution")
    , username_("username", "Username", "s_leis06")
    , emailAddress_("emailAddress", "E-Mail Address", "s_leis06@uni-muenster.de")
    , clusterAddress_("clusterAddress", "Cluster Address", "palma2c.uni-muenster.de")
    , programPath_("programPath", "Program Path", "~/OpenLB")
    , dataPath_("dataPath", "Data Path", "/scratch/tmp")
    , toolchain_("toolchain", "Toolchain")
    , simulationType_("simulationType", "Simulation Type")
    , configNodes_("configNodes", "Nodes", 1, 1, 2)
    , configTasks_("configTasks", "Tasks", 18, 1, 72)
    , configTasksPerNode_("configTasksPerNode", "Tasks per Node", 18, 1, 72)
    , configCPUsPerTask_("configCPUsPerTask", "CPUs per Task", 4, 1, 72)
    , configMemory_("configMemory", "Memory (GB/Node)", 16, 1, 32)
    , configPartition_("configPartition", "Partition")
    , configTimeDays_("configTimeDays", "Max. Time Days", 0, 0, 6)
    , configTimeQuarters_("configTimeQuarters", "Max. Quarters (1/4h)", 1, 1, 4*24)
    , simulationResults_("simulationResults", "Simulation Results", "Simulation Results", VoreenApplication::app()->getTemporaryPath("simulations"), "", FileDialogProperty::DIRECTORY, Processor::VALID, Property::LOD_DEFAULT, VoreenFileWatchListener::ALWAYS_OFF)
    , uploadDataPath_("uploadDataPath", "Upload Data Path", "Upload Data Path", VoreenApplication::app()->getTemporaryPath(), "", FileDialogProperty::DIRECTORY, Processor::VALID, Property::LOD_DEFAULT)
    , compileOnUpload_("compileOnUpload", "Compile on Upload", false)
    , deleteOnDownload_("deleteOnDownload", "Delete original Data", false)
    , triggerEnqueueSimulations_("triggerEnqueueSimulations", "Enqueue Simulations", Processor::VALID)
    , triggerFetchResults_("triggerFetchResults", "Fetch Results", Processor::VALID)
    , progress_("progress", "Progress")
{
    addPort(geometryDataPort_);
    addPort(measuredDataPort_); // Currently ignored
    measuredDataPort_.addCondition(new PortConditionVolumeListEnsemble());
    measuredDataPort_.addCondition(new PortConditionVolumeListAdapter(new PortConditionVolumeType3xFloat()));
    addPort(parameterPort_);

    addProperty(useLocalInstance_);
    //ON_CHANGE_LAMBDA(useLocalInstance_, [this] {
    //    setPropertyGroupVisible("cluster-general", !useLocalInstance_.get());
    //});
    useLocalInstance_.setGroupID("local-instance");
    addProperty(localInstancePath_);
    localInstancePath_.setGroupID("local-instance");
    setPropertyGroupGuiName("local-instance", "Local Instance");
//#ifndef WIN32
    // This is essentially a workaround for OpenLB not being compilable using MSVC.
    // So we compile it, e. g. using cygwin, and treat it as a local compute node.
    // TODO: see below!
    setPropertyGroupVisible("local-instance", false);
//#endif

    addProperty(institution_);
    ON_CHANGE(institution_, FlowSimulationCluster, institutionChanged);
    institution_.addOption("jena", "Jena");
    institution_.addOption("wwu", "WWU");
    institution_.setGroupID("cluster-general");
    addProperty(username_);
    username_.setGroupID("cluster-general");
    addProperty(emailAddress_);
    emailAddress_.setGroupID("cluster-general");
    addProperty(clusterAddress_);
    clusterAddress_.setGroupID("cluster-general");
    addProperty(programPath_);
    programPath_.setGroupID("cluster-general");
    addProperty(dataPath_);
    dataPath_.setGroupID("cluster-general");
    addProperty(toolchain_);
    toolchain_.setGroupID("cluster-general");
    toolchain_.addOption("foss", "foss");
    toolchain_.addOption("intel", "intel");
    addProperty(simulationType_);
    simulationType_.setGroupID("cluster-general");
    simulationType_.addOption("default", "default");
    //simulationType_.addOption("steered", "Steered"); // TODO: implement!
    simulationType_.addOption("aorta3d", "aorta3d");
    setPropertyGroupGuiName("cluster-general", "General Cluster Config");

    addProperty(configNodes_);
    configNodes_.setGroupID("cluster-resources");
    addProperty(configTasks_);
    configTasks_.setGroupID("cluster-resources");
    addProperty(configTasksPerNode_);
    configTasksPerNode_.setGroupID("cluster-resources");
    addProperty(configCPUsPerTask_);
    configCPUsPerTask_.setGroupID("cluster-resources");
    addProperty(configMemory_);
    configMemory_.setGroupID("cluster-resources");
    addProperty(configPartition_);
    configPartition_.addOption("normal", "normal");
    configPartition_.addOption("express", "express");
    configPartition_.setGroupID("cluster-resources");
    addProperty(configTimeDays_);
    configTimeDays_.setGroupID("cluster-resources");
    addProperty(configTimeQuarters_);
    configTimeQuarters_.setGroupID("cluster-resources");
    setPropertyGroupGuiName("cluster-resources", "Cluster Resource Config");

    addProperty(simulationResults_);
    simulationResults_.setGroupID("results");
    addProperty(uploadDataPath_);
    uploadDataPath_.setGroupID("results");
    addProperty(compileOnUpload_);
    compileOnUpload_.setGroupID("results");
    addProperty(deleteOnDownload_);
    deleteOnDownload_.setGroupID("results");
    addProperty(triggerEnqueueSimulations_);
    triggerEnqueueSimulations_.setGroupID("results");
    ON_CHANGE(triggerEnqueueSimulations_, FlowSimulationCluster, enqueueSimulations);
    addProperty(triggerFetchResults_);
    triggerFetchResults_.setGroupID("results");
    ON_CHANGE(triggerFetchResults_, FlowSimulationCluster, fetchResults);
    addProperty(progress_);
    progress_.setGroupID("results");
    setPropertyGroupGuiName("results", "Results");

    // Set visibility according to institution.
    institutionChanged();
}

FlowSimulationCluster::~FlowSimulationCluster() {
}

bool FlowSimulationCluster::isReady() const {
    if(!isInitialized()) {
        setNotReadyErrorMessage("Not initialized.");
        return false;
    }
    if(!geometryDataPort_.isReady()) {
        setNotReadyErrorMessage("Geometry Port not ready.");
        return false;
    }

    // Note: measuredDataPort ist optional!
    if(measuredDataPort_.hasData() && !measuredDataPort_.isReady()) {
        setNotReadyErrorMessage("Measured Data Port not ready.");
        return false;
    }

    if(!parameterPort_.isReady()) {
        setNotReadyErrorMessage("Parameter Port not ready.");
        return false;
    }

    return true;
}

void FlowSimulationCluster::process() {
    // Everything to process is done using callbacks.
}

void FlowSimulationCluster::institutionChanged() {
    if(institution_.get() == "jena") {
        setPropertyGroupVisible("cluster-resources", false);
    }
    else if(institution_.get() == "wwu") {
        setPropertyGroupVisible("cluster-resources", true);
    }
    else {
        tgtAssert(false, "unhandled institution");
    }
}

void FlowSimulationCluster::enqueueSimulations() {

    const GlMeshGeometryBase* geometryData = dynamic_cast<const GlMeshGeometryBase*>(geometryDataPort_.getData());
    if (!geometryData) {
        VoreenApplication::app()->showMessageBox("Error", "No simulation geometry. Did you perform the segmentation?", true);
        LERROR("Invalid simulation geometry");
        return;
    }

    const FlowParametrizationList* flowParametrization = parameterPort_.getData();
    if (!flowParametrization || flowParametrization->empty()) {
        VoreenApplication::app()->showMessageBox("Error", "No parametrization. Did you add one?", true);
        LERROR("No parametrization");
        return;
    }

    VoreenApplication::app()->showMessageBox("Information", "Start enqueuing runs, this may take a while...");
    LINFO("Configuring and enqueuing Simulations '" << simulationType_.get() << "'");

    // (Re-)compile code on cluster, if desired.
#ifdef WIN32
    if (!useLocalInstance_.get() && compileOnUpload_.get()) {
#else
    if (compileOnUpload_.get()) {
#endif

        // Execute compile script on the cluster.
        std::string command = "ssh " + username_.get() + "@" + clusterAddress_.get() + " " + generateCompileScript();
        int ret = executeCommand(command);

        // Check if compilation was successful.
        if (ret != EXIT_SUCCESS) {
            VoreenApplication::app()->showMessageBox("Error", "Could not compile program", true);
            LERROR("Could not compile program");
            return;
        }
    }

    std::string simulationPathSource = uploadDataPath_.get() + "/" + flowParametrization->getName() + "/";
    tgt::FileSystem::createDirectoryRecursive(simulationPathSource);
    std::string simulationPathDest = username_.get() + "@" + clusterAddress_.get() + ":" + programPath_.get() + "/" +
                                     toolchain_.get() + "/simulations/" + simulationType_.get() + "/";

    // Copy simulation geometry.
    // TODO: support changing/multiple geometries.
    tgt::FileSystem::createDirectory(simulationPathSource + "geometry/");
    std::string geometryFilename = simulationPathSource + "geometry/" + "geometry.stl";
    try {
        std::ofstream file(geometryFilename);
        geometryData->exportAsStl(file);
        file.close();
    }
    catch (std::exception& e) {
        VoreenApplication::app()->showMessageBox("Error", "Could not write geometry file", true);
        LERROR("Could not write geometry file: " << e.what());
        return;
    }

    // Copy velocity data.
    if(const VolumeList* volumeList = measuredDataPort_.getData()) {
        tgt::FileSystem::createDirectory(simulationPathSource + "velocity/");
        int nrLength = static_cast<int>(std::to_string(volumeList->size() - 1).size());
        for(size_t i=0; i<volumeList->size(); i++) {

            // Enumerate volumes.
            std::ostringstream suffix;
            suffix << std::setw(nrLength) << std::setfill('0') << i;
            std::string volumeName = "velocity" + suffix.str() + ".raw";

            std::string velocityFilename = simulationPathSource + "velocity/ " + volumeName;
            std::fstream velocityFile(velocityFilename.c_str(), std::ios::out | std::ios::binary);

            const VolumeRAM* volume = volumeList->at(i)->getRepresentation<VolumeRAM>();
            velocityFile.write(reinterpret_cast<const char*>(volume->getData()), volume->getNumBytes());
            if (!velocityFile.good()) {
                LERROR("Could not write velocity file");
                continue;
            }

            velocityFile.close();
        }
    }

    // Create configurations.
    for(size_t i=0; i<flowParametrization->size(); i++) {

        std::string parameterPathSource = simulationPathSource + flowParametrization->at(i).getName() + "/";
        tgt::FileSystem::createDirectoryRecursive(parameterPathSource);

        // Create parameter configuration file and commit to cluster.
        std::string parameterFilename = parameterPathSource + "config.xml";
        std::ofstream parameterFile(parameterFilename);
        if (!parameterFile.good()) {
            LERROR("Could not write parameter file");
            continue;
        }
        parameterFile << flowParametrization->toXMLString(i);
        parameterFile.close();

        std::string submissionScriptFilename = parameterPathSource + "submit.cmd";
        std::ofstream submissionScriptFile(submissionScriptFilename);
        if (!submissionScriptFile.good()) {
            LERROR("Could not write submission script file");
            continue;
        }
        submissionScriptFile << generateSubmissionScript(flowParametrization->at(i).getName());
        submissionScriptFile.close();
    }

#ifdef WIN32
    // Allow to use a local simulation program.
    if (useLocalInstance_.get()) {

        // TODO: make platform independent.
        // TODO: figure out and set working directory according to upload path.

        // Move directory to the very same place, the local instance is located.
        simulationPathSource = tgt::FileSystem::cleanupPath(simulationPathSource, true);
        simulationPathDest = tgt::FileSystem::cleanupPath(tgt::FileSystem::dirName(localInstancePath_.get()) + "/" + flowParametrization->getName(), true);
        std::string command = "move " + simulationPathSource + " " + simulationPathDest;
        if (executeCommand(command) != EXIT_SUCCESS) {
            VoreenApplication::app()->showMessageBox("Error", "Could not move ensemble", true);
            LERROR("Could not move ensemble");
            return;
        }

        // Run simulations one after the other.
        std::vector<std::string> failed;
        for (size_t i = 0; i < flowParametrization->size(); i++) {

            // Start job.
            std::string command = localInstancePath_.get() + " " + flowParametrization->getName() + " " + flowParametrization->at(i).getName() + " " + simulationResults_.get() + "/"; // Add a trailing '/' !
            int ret = executeCommand(command);
            if (ret != EXIT_SUCCESS) {
                failed.push_back(flowParametrization->at(i).getName());
                continue;
            }
        }

        if (!failed.empty()) {
            VoreenApplication::app()->showMessageBox("Error", "Some runs could not be enqueued. See log for details.", true);
            LERROR("Could not enqueue runs: \n* " << strJoin(failed, "\n* "));
        }
        else {
            VoreenApplication::app()->showMessageBox("Information", "All runs successfully enqueued!");
        }

        progress_.setProgress(1.0f);
        return;
    }
#endif

    // Copy data to cluster.
    std::string command = "scp -r " + simulationPathSource + " " + simulationPathDest;
    int ret = executeCommand(command);
    if (ret != EXIT_SUCCESS) {
        VoreenApplication::app()->showMessageBox("Error", "Data could not be copied!", true);
        LERROR("Data could not be copied");
        return;
    }

    // Don't delete data, will be handled by Voreen Temp. data!
    //tgt::FileSystem::deleteDirectoryRecursive(simulationPathSource);

    // Enqueue simulations.
    std::vector<std::string> failed;
    for(size_t i=0; i<flowParametrization->size(); i++) {

        progress_.setProgress(i * 1.0f / flowParametrization->size());

        std::string localSimulationPath = programPath_.get() + "/" + toolchain_.get() + "/simulations/" + simulationType_.get() + "/" +
                                          flowParametrization->getName() + "/" + flowParametrization->at(i).getName();

        // Enqueue job.
        command = "ssh " + username_.get() + "@" + clusterAddress_.get() + " " + generateEnqueueScript(localSimulationPath);
        ret = executeCommand(command);
        if (ret != EXIT_SUCCESS) {
            failed.push_back(flowParametrization->at(i).getName());
            continue;
        }
    }

    if (!failed.empty()) {
        VoreenApplication::app()->showMessageBox("Error", "Some runs could not be enqueued. See log for details.", true);
        LERROR("Could not enqueue runs: \n* " << strJoin(failed, "\n* "));
    }
    else {
        VoreenApplication::app()->showMessageBox("Information", "All runs successfully enqueued!");
    }

    progress_.setProgress(1.0f);
}

void FlowSimulationCluster::fetchResults() {
    std::string username = username_.get();
    std::string clusterAddress = clusterAddress_.get();

    VoreenApplication::app()->showMessageBox("Information", "Start fetching data, this may take a while...");

    // Create output directory.
    std::string directory = simulationResults_.get();
    if (!tgt::FileSystem::createDirectoryRecursive(directory)) {
        LERROR("Could not create output directory: " << directory);
        return;
    }

    std::string address = username + "@" + clusterAddress;
    std::string simulationPath = dataPath_.get() + "/" + username + "/simulations/" + simulationType_.get();
    std::string source = address + ":" + simulationPath;
    const FlowParametrizationList* parametrizationList = parameterPort_.getData();
    if (parametrizationList) {
        source += "/" + parametrizationList->getName();

        std::string dest = directory + "/" + parametrizationList->getName();
        if (!tgt::FileSystem::dirExists(dest) && !tgt::FileSystem::createDirectory(dest)) {
            LERROR("Could not create ensemble directory: " << dest);
            return;
        }

        std::vector<std::string> failed;
        bool deletionFailed = false;
        for(const FlowParameters& parameters : parametrizationList->getFlowParametrizations()) {
            std::string command = "scp -r ";
            command += source + "/" + parameters.getName() + " ";
            command += dest;

            int ret = executeCommand(command);
            if (ret != EXIT_SUCCESS) {
                failed.push_back(parametrizationList->getName() + "/" + parameters.getName());
            }

            if(deleteOnDownload_.get()) {
                command = "ssh " + address + " \"rm -rf " + simulationPath + "/" + parametrizationList->getName() + "/" + parameters.getName() + "\"";
                ret = executeCommand(command);
                if (ret != EXIT_SUCCESS) {
                    deletionFailed = true;
                }
            }
        }

        if(deletionFailed && failed.empty()) {
            VoreenApplication::app()->showMessageBox("Error", "Some data could not be deleted. You may clean up manually.", true);
        }

        if(!failed.empty()) {
            std::string message = "Some data could not be fetched. See log for details.";
            if(deletionFailed) {
                message += "\n\nAlso, some data could not be deleted. You may clean up manually.";
            }
            VoreenApplication::app()->showMessageBox("Error", message, true);
            LERROR("Could not fetch results of: \n* " << strJoin(failed, "\n* "));
            return;
        }
    }
    else {
        std::string command = "scp -r " + source + " " + simulationResults_.get();
        int ret = executeCommand(command);
        if (ret != EXIT_SUCCESS) {
            VoreenApplication::app()->showMessageBox("Error", "Data could not be fetched!", true);
            LERROR("Could not fetch results");
            return;
        }

        if(deleteOnDownload_.get()) {
            command = "ssh " + address + " \" rm -rf " + simulationPath + "\"";
            ret = executeCommand(command);
            if (ret != EXIT_SUCCESS) {
                VoreenApplication::app()->showMessageBox("Warning", "Data could not be deleted!", true);
                LERROR("Could not delete data");
                return;
            }
        }
    }

    VoreenApplication::app()->showMessageBox("Information", "Fetching data was successful!");
}

int FlowSimulationCluster::executeCommand(const std::string& command) const {
#ifndef WIN32
    std::string result;
    std::array<char, 128> buffer;

    auto pipe = popen(command.c_str(), "r");
    if (!pipe) {
        LERROR("Failed to open pipe");
        return EXIT_FAILURE;
    }

    while (!feof(pipe)) {
        if (fgets(buffer.data(), buffer.size(), pipe) != nullptr) {
            result += buffer.data();
        }
    }

    if(!result.empty()) {
        LINFO(result);
    }

    return pclose(pipe);
#else
    return system(command.c_str());
#endif
}

std::string FlowSimulationCluster::generateCompileScript() const {
    std::stringstream script;
    script << "\"";

    if(institution_.get() == "jena") {
        script << "cd " << programPath_.get() << "/" << toolchain_.get() << "/simulations/" << simulationType_.get();
        script << " && make -j4"; // Assumes that using 4 threads is okay, which should be the case on any system.
    }
    else if(institution_.get() == "wwu") {
        script << "module load " << toolchain_.get();
        script << " && cd " << programPath_.get() << "/" << toolchain_.get() << "/simulations/" << simulationType_.get();
        script << " && make -j4"; // Assumes that using 4 threads is okay, which should be the case on any system.
    }

    script << "\"";
    return script.str();
}

std::string FlowSimulationCluster::generateEnqueueScript(const std::string& parametrizationPath) const {
    std::stringstream script;
    script << "\"";

    if(institution_.get() == "jena") {
        script << "cd " << parametrizationPath;
#ifdef WIN32
        // Using windows, the script's DOS line breaks need to be converted to UNIX line breaks.
        script << " && sed -i 's/\\r$//' submit.cmd";
#endif
        script << " && qsub submit.cmd";
    }
    else if(institution_.get() == "wwu") {
        script << "module load " << toolchain_.get();
        script << " && cd " << parametrizationPath;
#ifdef WIN32
        // Using windows, the script's DOS line breaks need to be converted to UNIX line breaks.
        script << " && sed -i 's/\\r$//' submit.cmd";
#endif
        script << " && sbatch submit.cmd";
    }

    script << "\"";
    return script.str();
}

std::string FlowSimulationCluster::generateSubmissionScript(const std::string& parametrizationName) const {
    tgtAssert(parameterPort_.hasData(), "no data");
    std::stringstream script;

    if(institution_.get() == "jena") {

        script << "#!/bin/bash" << std::endl;
        script << "# The job should be placed into the queue 'all.q'." << std::endl;
        script << "#$ -q short.q" << std::endl;
        script << "# E-Mail notification to...." << std::endl;
        script << "#$ -M " << emailAddress_.get() << std::endl;
        script << "# Report on finished and abort." << std::endl;
        script << "#$ -m bea" << std::endl;
        script << "#$ -N " << parameterPort_.getData()->getName() << "_" << parametrizationName << std::endl;
        script << "cd " << programPath_.get() << "/" << toolchain_.get() << "/simulations/" << simulationType_.get()
               << "/" << parameterPort_.getData()->getName() << "/" << parametrizationName << std::endl;
        script << "# This is the file to be executed." << std::endl;
        script << programPath_.get() << "/" << toolchain_.get() << "/simulations/" << simulationType_.get()
               << "/" << simulationType_.get();
        // First argument: ensemble name
        script << " " << parameterPort_.getData()->getName();
        // Second argument: run name
        script << " " << parametrizationName;
        // Third argument: output directory
        script << " "
               << dataPath_.get() + "/" + username_.get() + "/simulations/" // Cluster code needs a trailing '/' !
               << " > output.log" // forward output to output.log.
               << std::endl;

    }
    else if(institution_.get() == "wwu") {

        int quarters = configTimeQuarters_.get();
        int minutes = quarters * 15;
        int hours = minutes / 60;
        minutes = minutes % 60;

        script << "#!/bin/bash" << std::endl;
        script << std::endl;
        script << "# set the number of nodes" << std::endl;
        script << "#SBATCH --nodes=" << configNodes_.get() << std::endl;
        script << std::endl;
        script << "# MPI/OMP Hybrid config" << std::endl;
        script << "#SBATCH --ntasks=" << configTasks_.get() << std::endl;
        script << "#SBATCH --ntasks-per-node=" << configTasksPerNode_.get() << std::endl;
        script << "#SBATCH --cpus-per-task=" << configCPUsPerTask_.get() << std::endl;
        script << std::endl;
        // Don't use exclusive access, it might block other jobs if number of tasks is low.
        //script << "# set the number of CPU cores per node" << std::endl;
        //script << "#SBATCH --exclusive" << std::endl;
        //script << std::endl;
        script << "# How much memory is needed (per node)" << std::endl;
        script << "#SBATCH --mem=" << configMemory_.get() << "G" << std::endl;
        script << std::endl;
        script << "# set a partition" << std::endl;
        script << "#SBATCH --partition " << configPartition_.get() << std::endl;
        script << std::endl;
        script << "# set max wallclock time" << std::endl;
        script << "#SBATCH --time=" << configTimeDays_.get() << "-" << std::setw(2) << std::setfill('0') << hours << ":"
               << std::setw(2) << std::setfill('0') << minutes << ":00" << std::endl;
        script << std::endl;
        script << "# set name of job" << std::endl;
        script << "#SBATCH --job-name=" << parameterPort_.getData()->getName() << "-" << parametrizationName
               << std::endl;
        script << std::endl;
        script << "# mail alert at start, end and abortion of execution" << std::endl;
        script << "#SBATCH --mail-type=ALL" << std::endl;
        script << std::endl;
        script << "# set an output file" << std::endl;
        script << "#SBATCH --output output.log" << std::endl;
        script << std::endl;
        script << "# send mail to this address" << std::endl;
        script << "#SBATCH --mail-user=" << emailAddress_.get() << std::endl;
        script << std::endl;
        script << "# run the application" << std::endl;
        if (configCPUsPerTask_.get() > 1) {
            script << "OMP_NUM_THREADS=" << configCPUsPerTask_.get() << " ";
        }
        if (configTasks_.get() > 1 || configTasksPerNode_.get() > 1) {
            script << "mpirun ";
        }

        // Add executable.
        script << programPath_.get() << "/" << toolchain_.get() << "/simulations/" << simulationType_.get()
               << "/" << simulationType_.get();
        // First argument: ensemble name
        script << " " << parameterPort_.getData()->getName();
        // Second argument: run name
        script << " " << parametrizationName;
        // Third argument: output directory
        script << " " << dataPath_.get() + "/" + username_.get() + "/simulations/"; // Cluster code needs a trailing '/'

        script << std::endl;

    }


    return script.str();
}

}   // namespace
