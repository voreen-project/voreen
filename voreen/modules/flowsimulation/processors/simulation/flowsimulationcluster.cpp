/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2024 University of Muenster, Germany,                        *
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
#include "voreen/core/datastructures/volume/volumedecorator.h"
#include "voreen/core/ports/conditions/portconditionvolumelist.h"

#include "modules/core/io/rawvolumereader.h"
#include "modules/core/io/vvdvolumewriter.h"
#ifdef VRN_MODULE_VTK
#include "modules/vtk/io/vtivolumewriter.h"
#endif

#include <boost/process.hpp>


namespace {

    /**
     * Executes a system command.
     * On Linux, it also logs the output to the debug console.
     */
    int executeCommand(const std::string& command) {
#ifndef WIN32
        std::string result;
        std::array<char, 128> buffer;

        auto pipe = popen(command.c_str(), "r");
        if (!pipe) {
            LERRORC("System", "Failed to open pipe");
            return EXIT_FAILURE;
        }

        while (!feof(pipe)) {
            if (fgets(buffer.data(), buffer.size(), pipe) != nullptr) {
                result += buffer.data();
            }
        }

        if (!result.empty()) {
            LINFOC("System", result);
        }

        return pclose(pipe);
#else
        return system(command.c_str());
#endif
    }

    /**
     * Copies an entire directory from src to dst.
     * If abortOnError is set, the copy process is being stopped as soon as an error occurs.
     */
    bool copyDirectory(const boost::filesystem::path& src, const boost::filesystem::path& dst, bool abortOnError = false) {
        if (boost::filesystem::exists(dst)) {
            return false;
        }

        if (boost::filesystem::is_directory(src)) {
            bool success = true;
            boost::filesystem::create_directories(dst);
            for (auto& item : boost::filesystem::directory_iterator(src)) {
                success &= copyDirectory(item.path(), dst / item.path().filename(), abortOnError);
                if (!success && abortOnError) {
                    return false;
                }
            }
            return success;
        }
        else if (boost::filesystem::is_regular_file(src)) {
            boost::system::error_code ec;
            boost::filesystem::copy(src, dst, ec);
            if (ec.value() == 0) {
                return true;
            }
        }

        return false;
    }

}


namespace voreen {


/**
 * Simple background thread that allows to execute a system command asyncronously.
 */
class ExecutorProcess {
public:

    ExecutorProcess(const std::string& cd, const std::string& command, const std::string& name, bool detach)
        : cd_(cd)
        , command_(command)
        , name_(name)
    {}

    ~ExecutorProcess() {
        process_.terminate();
    }

    bool isFinished() {
        return !process_.running();
    }

    void run() {
        // We need to change the current working directory accordingly,
        // but also need to restore the old one afterward.
        auto path = boost::filesystem::current_path();
        boost::filesystem::current_path(cd_);
        process_ = boost::process::child(command_);
        if (detach_) {
            process_.detach();
        }
        boost::filesystem::current_path(path);
    }

    const std::string& getName() const {
        return name_;
    }

    bool successful() const {
        return process_.exit_code() == EXIT_SUCCESS;
    }

private:

    std::string cd_;
    std::string command_;
    std::string name_;
    bool detach_;

    boost::process::child process_;
};


const std::string FlowSimulationCluster::loggerCat_("voreen.flowsimulation.FlowSimulationCluster");

FlowSimulationCluster::FlowSimulationCluster()
    : Processor()
    // ports
    , geometryDataPort_(Port::INPORT, "geometryDataPort", "Geometry Input", false)
    , geometryVolumeDataPort_(Port::INPORT, "geometryVolumeDataPort", "Segmentation Input", false)
    , measuredDataPort_(Port::INPORT, "measuredDataPort", "Measured Data Input", false)
    , configPort_(Port::INPORT, "parameterPort", "Simulation Config", false)
    , useLocalInstance_("useLocalInstance", "Use local Instance", false)
    , localInstancePath_("localInstancePath", "Local Instance Path", "Path", "", "", FileDialogProperty::OPEN_FILE, Processor::INVALID_RESULT, Property::LOD_DEFAULT, VoreenFileWatchListener::ALWAYS_OFF)
    , detachProcesses_("detachProcesses", "Detach Processes", true)
    , overwriteExistingConfig_("overwriteExistingConfig", "Overwrite Existing Config", true)
    , stopProcesses_("stopProcesses", "Stop Runs")
    , workloadManager_("institution", "Institution")
    , username_("username", "Username", "s_leis06")
    , emailAddress_("emailAddress", "E-Mail Address", "s_leis06@uni-muenster.de")
    , clusterAddress_("clusterAddress", "Cluster Address", "palma.uni-muenster.de")
    , programPath_("programPath", "Program Path", "~")
    , dataPath_("dataPath", "Data Path", "/scratch/tmp")
    , simulationType_("simulationType", "Simulation Type")
    , configToolchain_("toolchain", "Toolchain")
    , configNodes_("configNodes", "Nodes", 1, 1, 2)
    , configNumGPUs_("configNumGPUs", "Num GPUs", 0, 0, 8)
    , configTasksPerNode_("configTasksPerNode", "Tasks per Node", 9, 1, 36)
    , configCPUsPerTask_("configCPUsPerTask", "CPUs per Task", 4, 1, 36)
    , configMemory_("configMemory", "Memory (GB/Node)", 16, 1, 92)
    , configPartition_("configPartition", "Partition")
    , configTimeDays_("configTimeDays", "Max. Time Days", 0, 0, 6)
    , configTimeQuarters_("configTimeQuarters", "Max. Quarters (1/4h)", 1, 1, 4*24)
    , refreshClusterCode_("refreshClusterCode", "Refresh Cluster Code")
    , simulationResults_("simulationResults", "Simulation Results", "Simulation Results", VoreenApplication::app()->getTemporaryPath("simulations"), "", FileDialogProperty::DIRECTORY, Processor::VALID, Property::LOD_DEFAULT, VoreenFileWatchListener::ALWAYS_OFF)
    , uploadDataPath_("uploadDataPath", "Upload Data Path", "Upload Data Path", VoreenApplication::app()->getTemporaryPath(), "", FileDialogProperty::DIRECTORY, Processor::VALID, Property::LOD_DEFAULT)
    , deleteOnDownload_("deleteOnDownload", "Delete original Data", false)
    , triggerEnqueueSimulations_("triggerEnqueueSimulations", "Enqueue Simulations")
    , triggerFetchResults_("triggerFetchResults", "Fetch Results")
    , progress_("progress", "Progress")
    , numEnqueuedThreads_(0)
    , numFinishedThreads_(0)
{
    addPort(geometryDataPort_);
    addPort(geometryVolumeDataPort_);
    geometryVolumeDataPort_.addCondition(new PortConditionVolumeListEnsemble());
    geometryVolumeDataPort_.addCondition(new PortConditionVolumeListAdapter(new PortConditionVolumeChannelCount(1)));
    addPort(measuredDataPort_);
    measuredDataPort_.addCondition(new PortConditionVolumeListEnsemble());
    measuredDataPort_.addCondition(new PortConditionVolumeListAdapter(new PortConditionVolumeChannelCount(3)));
    addPort(configPort_);

    addProperty(useLocalInstance_);
    ON_CHANGE_LAMBDA(useLocalInstance_, [this] {
        setPropertyGroupVisible("cluster-general", !useLocalInstance_.get());
        setPropertyGroupVisible("cluster-resources", !useLocalInstance_.get());
        if (!useLocalInstance_.get()) {
            workloadManagerChanged();
        }
        deleteOnDownload_.setVisibleFlag(!useLocalInstance_.get());
        triggerFetchResults_.setVisibleFlag(!useLocalInstance_.get());
    });
    useLocalInstance_.setGroupID("local-instance");
    addProperty(localInstancePath_);
    localInstancePath_.setGroupID("local-instance");
    addProperty(detachProcesses_);
    detachProcesses_.setGroupID("local-instance");
    addProperty(overwriteExistingConfig_);
    overwriteExistingConfig_.setGroupID("local-instance");
    addProperty(stopProcesses_);
    ON_CHANGE(stopProcesses_, FlowSimulationCluster, threadsStopped);
    stopProcesses_.setGroupID("local-instance");
    setPropertyGroupGuiName("local-instance", "Local Instance");

    addProperty(workloadManager_);
    ON_CHANGE(workloadManager_, FlowSimulationCluster, workloadManagerChanged);
    workloadManager_.addOption("slurm", "slurm");
    workloadManager_.addOption("qsub", "qsub");
    workloadManager_.setGroupID("cluster-general");
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
    //addProperty(simulationType_);
    //simulationType_.setGroupID("cluster-general");
    //simulationType_.addOption("default", "default");
    setPropertyGroupGuiName("cluster-general", "General Cluster Config");

    addProperty(configToolchain_);
    configToolchain_.setGroupID("cluster-resources");
    configToolchain_.addOption("foss", "foss");
    configToolchain_.addOption("intel", "intel");
    addProperty(configNodes_);
    configNodes_.setGroupID("cluster-resources");
    addProperty(configNumGPUs_);
    configNumGPUs_.setGroupID("cluster-resources");
    addProperty(configTasksPerNode_);
    configTasksPerNode_.setGroupID("cluster-resources");
    addProperty(configCPUsPerTask_);
    configCPUsPerTask_.setGroupID("cluster-resources");
    addProperty(configMemory_);
    configMemory_.setGroupID("cluster-resources");
    addProperty(configPartition_);
    configPartition_.addOption("normal", "normal");
    configPartition_.addOption("express", "express");
    configPartition_.addOption("gpu2080", "gpu2080");
    configPartition_.setGroupID("cluster-resources");
    addProperty(configTimeDays_);
    configTimeDays_.setGroupID("cluster-resources");
    addProperty(configTimeQuarters_);
    configTimeQuarters_.setGroupID("cluster-resources");
    addProperty(refreshClusterCode_);
    ON_CHANGE(refreshClusterCode_, FlowSimulationCluster, refreshClusterCode);
    refreshClusterCode_.setGroupID("cluster-resources");
    setPropertyGroupGuiName("cluster-resources", "Cluster Resource Config");

    addProperty(simulationResults_);
    simulationResults_.setGroupID("results");
    addProperty(uploadDataPath_);
    uploadDataPath_.setGroupID("results");
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
    workloadManagerChanged();
}

FlowSimulationCluster::~FlowSimulationCluster() {
    threadsStopped();
}

bool FlowSimulationCluster::isReady() const {
    if(!isInitialized()) {
        setNotReadyErrorMessage("Not initialized.");
        return false;
    }

    if(!configPort_.isReady()) {
        setNotReadyErrorMessage("Parameter Port not ready.");
        return false;
    }

    return true;
}

void FlowSimulationCluster::process() {
    // Most of the logic is handled by callbacks.
    // We just use the process method to start threads, if necessary.

    for (auto iter = runningThreads_.begin(); iter != runningThreads_.end();) {
        auto run = iter->get();
        if (run->isFinished()) {
            LINFO("Run " << run->getName() << " finished " << (run->successful() ? "successfully" : "unsuccessfully"));
            iter = runningThreads_.erase(iter);
            progress_.setProgress(static_cast<float>(++numFinishedThreads_) / static_cast<float>(numEnqueuedThreads_));
        }
        else {
            iter++;
        }
    }

    size_t maxNumThreads = boost::thread::hardware_concurrency();
    while (!waitingThreads_.empty() && runningThreads_.size() < maxNumThreads) {
        auto thread = std::move(waitingThreads_.front());
        waitingThreads_.pop_front();

        thread->run();
        LINFO("Run " << thread->getName() << " started");
        runningThreads_.push_back(std::move(thread));
    }

    if (!runningThreads_.empty()) {
        invalidate(INVALID_RESULT);
    }
}

void FlowSimulationCluster::threadsStopped() {
    runningThreads_.clear();
    waitingThreads_.clear();
    numEnqueuedThreads_ = 0;
    numFinishedThreads_ = 0;
}

void FlowSimulationCluster::workloadManagerChanged() {
    if(workloadManager_.get() == "qsub") {
        setPropertyGroupVisible("cluster-resources", false);
    }
    else if(workloadManager_.get() == "slurm") {
        setPropertyGroupVisible("cluster-resources", true);
    }
    else {
        tgtAssert(false, "unhandled workload manager");
    }
}

void FlowSimulationCluster::refreshClusterCode() {
    LINFO("Uploading program, this may take a while...");

    // Copy code to cluster.
    //std::string modulePath = getModulePath(); // Does not work!
    std::string modulePath = VoreenApplication::app()->getModulePath("flowsimulation");
    std::string simulationPathSource = modulePath + "/ext/openlb";
    std::string simulationPathDest = username_.get() + "@" + clusterAddress_.get() + ":" + programPath_.get();
    std::string command = "scp -r " + simulationPathSource + " " + simulationPathDest;
    int ret = executeCommand(command);
    if (ret != EXIT_SUCCESS) {
        VoreenApplication::app()->showMessageBox("Error", "Code could not be copied", true);
        LERROR("Code could not be copied");
        return;
    }

    LINFO("Compiling program, this may take a while...");

    // Execute compile script on the cluster.
    command = "ssh " + username_.get() + "@" + clusterAddress_.get() + " " + generateCompileScript();
    ret = executeCommand(command);

    // Check if compilation was successful.
    if (ret != EXIT_SUCCESS) {
        VoreenApplication::app()->showMessageBox("Error", "Error compiling program", true);
        LERROR("Error compiling program");
    }
    else {
        VoreenApplication::app()->showMessageBox("Finished", "Compilation finished", false);
        LINFO("Compilation finished");
    }
}

void FlowSimulationCluster::stepCopyGeometryData(FlowSimulationConfig& config, const std::string& simulationPathSource) {
    if(const auto* geometryData = dynamic_cast<const GlMeshGeometryBase*>(geometryDataPort_.getData())) {
        std::unique_ptr<GlMeshGeometryBase> geometry(dynamic_cast<GlMeshGeometryBase*>(geometryData->clone().release()));
        tgt::FileSystem::createDirectory(simulationPathSource + "geometry/");
        std::string geometryFilename = simulationPathSource + "geometry/" + "geometry.stl";

        try {
            std::ofstream file(geometryFilename);
            geometry->setTransformationMatrix(config.getTransformationMatrix() * geometry->getTransformationMatrix());
            geometry->exportAsStl(file);
            file.close();
        }
        catch (std::exception& e) {
            std::string errorMessage = "Could not write geometry file: ";
            throw VoreenException(errorMessage + e.what());
        }
        auto geometryFiles = std::map<float, std::string>();
        geometryFiles[0.0f] = "geometry/geometry.stl";
        config.setGeometryFiles(geometryFiles, true);
    }
    else if (const auto* segmentationData = geometryVolumeDataPort_.getData()) {
        auto geometryFiles = checkAndConvertVolumeList(segmentationData, config.getTransformationMatrix(), simulationPathSource, "geometry");
        config.setGeometryFiles(geometryFiles, false);
    }

    if (config.getGeometryFiles().empty()) {
        throw VoreenException("No GlMeshGeometry data.");
    }
}

void FlowSimulationCluster::stepCopyMeasurementData(const VolumeList* volumeList, FlowSimulationConfig& config, const std::string& simulationPathSource) {
    auto measurementData = checkAndConvertVolumeList(volumeList, config.getTransformationMatrix(), simulationPathSource, "measured");
    config.setMeasuredDataFiles(measurementData);
}

void FlowSimulationCluster::stepCreateSimulationConfigs(FlowSimulationConfig& config, const std::string& simulationPathSource) {
    for(size_t i=0; i<config.size(); i++) {

        std::string parameterPathSource = simulationPathSource + config.at(i).name_ + "/";
        tgt::FileSystem::createDirectoryRecursive(parameterPathSource);

        // Create parameter configuration file and commit to cluster.
        std::string parameterFilename = parameterPathSource + "config.xml";
        std::ofstream parameterFile(parameterFilename);
        if (!parameterFile.good()) {
            LERROR("Could not write parameter file");
            continue;
        }
        parameterFile << config.toXMLString(i);
        parameterFile.close();

        std::string submissionScriptFilename = parameterPathSource + "submit.cmd";
        std::ofstream submissionScriptFile(submissionScriptFilename);
        if (!submissionScriptFile.good()) {
            LERROR("Could not write submission script file");
            continue;
        }
        submissionScriptFile << generateSubmissionScript(config.at(i).name_);
        submissionScriptFile.close();
    }
}

void FlowSimulationCluster::runLocal(FlowSimulationConfig& config, std::string simulationPathSource, std::string simulationPathDest) {
    // Move directory to the very same place, the local instance is located.
    simulationPathSource = tgt::FileSystem::cleanupPath(simulationPathSource, true);
    simulationPathDest = tgt::FileSystem::cleanupPath(tgt::FileSystem::dirName(localInstancePath_.get()) + "/" + config.getName(), true);

    if(tgt::FileSystem::dirExists(simulationPathDest) && overwriteExistingConfig_.get()) {
        tgt::FileSystem::deleteDirectoryRecursive(simulationPathDest);
    }

    if(!copyDirectory(simulationPathSource, simulationPathDest, true)) {
        VoreenApplication::app()->showMessageBox("Error", "Could not copy configurations, they might already exist", true);
        return;
    }

    // Enqueue jobs.
    for (size_t i = 0; i < config.size(); i++) {

        std::string ensemble = config.getName();
        std::string run = config.at(i).name_;

        std::string workingDirectory = tgt::FileSystem::cleanupPath(tgt::FileSystem::dirName(localInstancePath_.get()) + "/" + ensemble + "/" + run, true);
        std::string runCommand = localInstancePath_.get() + " " + ensemble + " " + run + " " + simulationResults_.get() + "/"; // Add a trailing '/' !;
        std::string name = ensemble + "-" + run;
        waitingThreads_.emplace_back(
                std::unique_ptr<ExecutorProcess>(
                        new ExecutorProcess(workingDirectory, runCommand, name, detachProcesses_.get())
                )
        );
        numEnqueuedThreads_++;
        LINFO("Enqueued run " << name);
    }
}
void FlowSimulationCluster::runCluster(FlowSimulationConfig& config, std::string simulationPathSource, std::string simulationPathDest) {
    // Copy data to cluster.
    std::string command = "scp -r " + simulationPathSource + " " + simulationPathDest;
    int ret = executeCommand(command);
    if (ret != EXIT_SUCCESS) {
        throw VoreenException("Data could not be copied");
    }

    // Don't delete data, will be handled by Voreen Temp. data!
    //tgt::FileSystem::deleteDirectoryRecursive(simulationPathSource);

    // Enqueue simulations.
    std::vector<std::string> failed;
    for(size_t i=0; i<config.size(); i++) {

        progress_.setProgress(i * 1.0f / config.size());

        std::string localSimulationPath = programPath_.get() + "/openlb/voreen/" + config.getName() + "/" + config.at(i).name_;

        // Enqueue job.
        command = "ssh " + username_.get() + "@" + clusterAddress_.get() + " " + generateEnqueueScript(localSimulationPath);
        ret = executeCommand(command);
        if (ret != EXIT_SUCCESS) {
            failed.push_back(config.at(i).name_);
        }
    }

    progress_.setProgress(1.0f);

    if (!failed.empty()) {
        VoreenApplication::app()->showMessageBox("Error", "Some runs could not be enqueued. See log for details.", true);
        LERROR("Could not enqueue runs: \n* " << strJoin(failed, "\n* "));
    }
    else {
        VoreenApplication::app()->showMessageBox("Information", "All runs successfully enqueued!");
    }
}

void FlowSimulationCluster::enqueueSimulations() {

    const FlowSimulationConfig* originalConfig = configPort_.getData();
    if (!originalConfig || originalConfig->empty()) {
        VoreenApplication::app()->showMessageBox("Error", "No parametrization. Did you add one?", true);
        LERROR("No parametrization");
        return;
    }

    // We need to create a copy here, since we modify the paths of the individual runs.
    FlowSimulationConfig config(*originalConfig);

    if(config.getFlowFeatures() == FF_NONE) {
        VoreenApplication::app()->showMessageBox("Error", "No flow feature selected. Did you add one?", true);
        LERROR("No flow feature selected");
        return;
    }

    // Make sure we have a binary selected.
    if(useLocalInstance_.get() && localInstancePath_.get().empty()) {
        VoreenApplication::app()->showMessageBox("Error",
                                                 "No local instance binary selected. It is normally located here: "
                                                 "'voreen/modules/flowsimulation/ext/openlb/voreen/simulation_cluster'", true);
        return;
    }

    LINFO("Configuring and enqueuing Simulations...");

    std::string simulationPathSource = uploadDataPath_.get() + "/" + config.getName() + "/";
    tgt::FileSystem::deleteDirectoryRecursive(simulationPathSource);
    tgt::FileSystem::createDirectoryRecursive(simulationPathSource);
    std::string simulationPathDest = username_.get() + "@" + clusterAddress_.get() + ":" + programPath_.get() + "/openlb/voreen/";

    try {

        // Copy simulation geometry.
        stepCopyGeometryData(config, simulationPathSource);

        // Copy measurement data.
        stepCopyMeasurementData(measuredDataPort_.getData(), config, simulationPathSource);

        // Create configurations.
        stepCreateSimulationConfigs(config, simulationPathSource);

        // Allow to use a local simulation program.
        if (useLocalInstance_.get()) {
            runLocal(config, simulationPathSource, simulationPathDest);
        } else {
            runCluster(config, simulationPathSource, simulationPathDest);
        }
    }
    catch(VoreenException& e) {
        LERROR(e.what());
        VoreenApplication::app()->showMessageBox("Error", e.what(), true);
    }
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
    const FlowSimulationConfig* config = configPort_.getData();
    if (config) {
        source += "/" + config->getName();

        std::string dest = directory + "/" + config->getName();
        if (!tgt::FileSystem::dirExists(dest) && !tgt::FileSystem::createDirectory(dest)) {
            LERROR("Could not create ensemble directory: " << dest);
            return;
        }

        std::vector<std::string> failed;
        bool deletionFailed = false;
        for(const Parameters& parameters : config->getFlowParameterSets()) {
            std::string command = "scp -r ";
            command += source + "/" + parameters.name_ + " ";
            command += dest;

            int ret = executeCommand(command);
            if (ret != EXIT_SUCCESS) {
                failed.push_back(config->getName() + "/" + parameters.name_);
            }

            if(deleteOnDownload_.get()) {
                command = "ssh " + address + " \"rm -rf " + simulationPath + "/" + config->getName() + "/" + parameters.name_ + "\"";
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

std::string FlowSimulationCluster::generateCompileScript() const {
    std::stringstream script;
    script << "\"";

    std::string options = "PARALLEL_MODE=HYBRID ";
    if(configCPUsPerTask_.get() > 0) {
        options += "OMPFLAGS=-fopenmp ";
    }

    std::string platforms = "CPU_SISD";// CPU_SIMD"; // Might be possible in a later OpenLB release.
    if(configPartition_.get() == "gpu2080" && configNumGPUs_.get() > 0) {
        options += " GPU_CUDA";
    }
    options += "PLATFORMS='" + platforms + "' ";

    if(workloadManager_.get() == "qsub") {
        script << "cd " << programPath_.get() << "/openlb/voreen";
        script << " && make clean && " << options << "make";
    }
    else if(workloadManager_.get() == "slurm") {
        script << "module load palma/2021a ";

        if(configToolchain_.get() == "intel") {
            script << "intel/2021a";
            options += "CXX=mpicxx CC=icc CXXFLAGS='-std=c++17 -O3 -Wall -xHost -ipo' ";
        }
        else if(configToolchain_.get() == "foss"){
            script << "foss/2021a";
            options += "CXX=mpic++ CC=gcc CXXFLAGS='-std=c++17 -Wall -march=native -mtune=native' ";
        }

        script << " && cd " << programPath_.get() << "/openlb/voreen";
        script << " && make clean && " << options << "make";
    }

    script << "\"";
    return script.str();
}

std::string FlowSimulationCluster::generateEnqueueScript(const std::string& parametrizationPath) const {
    std::stringstream script;
    script << "\"";

    if(workloadManager_.get() == "qsub") {
        script << "cd " << parametrizationPath;
#ifdef WIN32
        // Using windows, the script's DOS line breaks need to be converted to UNIX line breaks.
        script << " && sed -i 's/\\r$//' submit.cmd";
#endif
        script << " && qsub submit.cmd";
    }
    else if(workloadManager_.get() == "slurm") {
        script << "cd " << parametrizationPath;
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
    tgtAssert(configPort_.hasData(), "no data");
    std::stringstream script;

    if(workloadManager_.get() == "qsub") {

        script << "#!/bin/bash" << std::endl;
        script << "# The job should be placed into the queue 'all.q'." << std::endl;
        script << "#$ -q short.q" << std::endl;
        script << "# E-Mail notification to...." << std::endl;
        script << "#$ -M " << emailAddress_.get() << std::endl;
        script << "# Report on finished and abort." << std::endl;
        script << "#$ -m bea" << std::endl;
        script << "#$ -N " << configPort_.getData()->getName() << "_" << parametrizationName << std::endl;
        script << "cd " << programPath_.get() << "/" << configToolchain_.get() << "/simulations/" << simulationType_.get()
               << "/" << configPort_.getData()->getName() << "/" << parametrizationName << std::endl;
        script << "# This is the file to be executed." << std::endl;
        script << programPath_.get() << "/" << configToolchain_.get() << "/simulations/" << simulationType_.get()
               << "/" << simulationType_.get();
        // First argument: ensemble name
        script << " " << configPort_.getData()->getName();
        // Second argument: run name
        script << " " << parametrizationName;
        // Third argument: output directory
        script << " "
               << dataPath_.get() + "/" + username_.get() + "/simulations/" // Cluster code needs a trailing '/' !
               << " > output.log" // forward output to output.log.
               << std::endl;

    }
    else if(workloadManager_.get() == "slurm") {

        int quarters = configTimeQuarters_.get();
        int minutes = quarters * 15;
        int hours = minutes / 60;
        minutes = minutes % 60;

        script << "#!/bin/bash" << std::endl;
        script << std::endl;
        script << "# set the number of nodes" << std::endl;
        script << "#SBATCH --nodes=" << configNodes_.get() << std::endl;
        script << std::endl;
        script << "# set the number of gpus" << std::endl;
        script << "#SBATCH --gres=gpu:" << configNumGPUs_.get() << std::endl;
        script << std::endl;
        script << "# MPI/OMP Hybrid config" << std::endl;
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
        script << "#SBATCH --job-name=" << configPort_.getData()->getName() << "-" << parametrizationName
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
        script << "module load palma/2021a foss/2021a" << std::endl; // TODO: make adjustable.
        script << std::endl;
        script << "# run the application" << std::endl;
        if (configCPUsPerTask_.get() > 1) {
            script << "OMP_NUM_THREADS=" << configCPUsPerTask_.get() << " ";
        }
        if (configTasksPerNode_.get() > 1) {
            script << "mpirun ";
        }
        if (configNumGPUs_.get() > 1) {
            script << "bash -c 'export CUDA_VISIBLE_DEVICES=${OMPI_COMM_WORLD_LOCAL_RANK}; ";
        }

        // Add executable.
        script << programPath_.get() << "/openlb/voreen/simulation_cluster";
        // First argument: ensemble name
        script << " " << configPort_.getData()->getName();
        // Second argument: run name
        script << " " << parametrizationName;
        // Third argument: output directory
        script << " " << dataPath_.get() + "/" + username_.get() + "/simulations/"; // Cluster code needs a trailing '/'

        script << std::endl;

    }


    return script.str();
}

#ifdef VRN_MODULE_VTK
std::map<float, std::string> FlowSimulationCluster::checkAndConvertVolumeList(const VolumeList* volumes, tgt::mat4 transformation, const std::string& simulationPathSource, const std::string& subdirectory) const {
    tgt::FileSystem::createDirectory(simulationPathSource + subdirectory + "/");
    int nrLength = static_cast<int>(std::to_string(volumes->size() - 1).size());

    std::map<float, std::string> result;
    for (size_t i = 0; i < volumes->size(); i++) {

        VolumeBase* volume = volumes->at(i);

        // Enumerate volumes.
        std::ostringstream suffix;
        suffix << std::setw(nrLength) << std::setfill('0') << i;
        std::string volumeName = subdirectory + suffix.str() + ".vti";
        std::string path = simulationPathSource + subdirectory + "/" + volumeName;

        if(transformation != tgt::mat4::identity) {
            std::unique_ptr<VolumeDecoratorReplace> transformedVolume(
                    new VolumeDecoratorReplaceTransformation(volume,
                                                             transformation * volume->getPhysicalToWorldMatrix())
            );
            VTIVolumeWriter().writeVectorField(path, transformedVolume.get());
        }
        else {
            VTIVolumeWriter().writeVectorField(path, volume);
        }

        std::string url = path + "?fieldName=" + volume->getModality().getName();
        result[volume->getTimestep()] = url;
    }

    return result;
}
#else
std::map<float, std::string> FlowSimulationCluster::checkAndConvertVolumeList(const VolumeList*, tgt::mat4, const std::string&, const std::string&) const {
    throw VoreenException("Need to convert to vti, but VTI module not enabled");
}
#endif

}   // namespace
