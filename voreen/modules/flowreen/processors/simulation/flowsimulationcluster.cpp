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
#include "voreen/core/datastructures/volume/volumefactory.h"
#include "voreen/core/datastructures/volume/volumeminmaxmagnitude.h"
#include "voreen/core/ports/conditions/portconditionvolumetype.h"

#include "modules/core/io/rawvolumereader.h"

#include "modules/flowreen/utils/geometryconverter.h"

#ifdef VRN_MODULE_ENSEMBLEANALYSIS
#include "voreen/core/io/volumeserializerpopulator.h"
#include "voreen/core/io/volumeserializer.h"

#include "custommodules/ensembleanalysis/utils/colorpool.h"
#endif

namespace voreen {

const std::string FlowSimulationCluster::loggerCat_("voreen.flowreen.FlowSimulationCluster");

FlowSimulationCluster::FlowSimulationCluster()
    : Processor()
    // ports
    , geometryDataPort_(Port::INPORT, "geometryDataPort", "Geometry Input", false)
    , measuredDataPort_(Port::INPORT, "measuredDataPort", "Measured Data Input", false)
    , parameterPort_(Port::INPORT, "parameterPort", "Parameterization", false)
#ifdef VRN_MODULE_ENSEMBLEANALYSIS
    , ensemblePort_(Port::OUTPORT, "ensembleport", "Ensemble Output", false)
#endif
    , username_("username", "Username", "s_leis06")
    , clusterAddress_("clusterAddress", "Cluster Address", "palma2c.uni-muenster.de")
    , simulationPath_("simulationPath", "Simulation Path", "~/OpenLB-intel/simulations")
    , simulationType_("simulationType", "Simulation Type")
    , configNodes_("configNodes", "Nodes", 1, 1, 2)
    , configTasks_("configTasks", "Tasks", 18, 1, 72)
    , configTasksPerNode_("configTasksPerNode", "Tasks per Node", 18, 1, 72)
    , configCPUsPerTask_("configCPUsPerTask", "CPUs per Task", 4, 1, 72)
    , configMemory_("configMemory", "Memory (GB/Node)", 16, 1, 32)
    , configPartition_("configPartition", "Partition")
    , simulationResults_("simulationResults", "Simulation Results", "Simulation Results", VoreenApplication::app()->getTemporaryPath("simulations"), "", FileDialogProperty::DIRECTORY, Processor::VALID)
    , triggerEnqueueSimulations_("triggerEnqueueSimulations", "Enqueue Simulations", Processor::VALID)
    , triggerFetchResults_("triggerFetchResults", "Fetch Results", Processor::VALID)
    , progress_("progress", "Progress")
{
    addPort(geometryDataPort_);
    addPort(measuredDataPort_); // Currently ignored
    //measuredDataPort_.addCondition(new PortConditionVolumeType3xFloat());
    addPort(parameterPort_);
#ifdef VRN_MODULE_ENSEMBLEANALYSIS
    addPort(ensemblePort_);
#endif

    addProperty(username_);
    username_.setGroupID("cluster-general");
    addProperty(clusterAddress_);
    clusterAddress_.setGroupID("cluster-general");
    addProperty(simulationPath_);
    simulationPath_.setGroupID("cluster-general");
    addProperty(simulationType_);
    simulationType_.setGroupID("cluster-general");
    simulationType_.addOption("default", "Default");
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
    setPropertyGroupGuiName("cluster-resources", "Cluster Resource Config");

    addProperty(simulationResults_);
    addProperty(triggerEnqueueSimulations_);
    ON_CHANGE(triggerEnqueueSimulations_, FlowSimulationCluster, enqueueSimulations);
    addProperty(triggerFetchResults_);
    ON_CHANGE(triggerFetchResults_, FlowSimulationCluster, fetchResults);
    addProperty(progress_);
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

    // Note: measuredDataPort is currently ignored!

    if(!parameterPort_.isReady()) {
        setNotReadyErrorMessage("Parameter Port not ready.");
        return false;
    }

    return true;
}

void FlowSimulationCluster::process() {
}

void FlowSimulationCluster::enqueueSimulations() {

    const Geometry* geometryData = geometryDataPort_.getData();
    if (!geometryData) {
        LERROR("No simulation geometry");
        return;
    }

    const FlowParametrizationList* flowParametrization = parameterPort_.getData();
    if (!flowParametrization || flowParametrization->empty()) {
        LERROR("No parameterization");
        return;
    }

    VoreenApplication::app()->showMessageBox("Information", "Start enqueuing data, this may take a while...");
    LINFO("Configuring and enqueuing Simulations '" << simulationType_.get() << "'");

    // This processor assumes:
    //      * on the cluster, there is a folder in home called 'OpenLB-intel/simulations'
    //      * all simulations have already been compiled
    //      * the result will be saved to /scratch/tmp/<user>/simulations/<simulation_name>/<run_name>

    std::string simulationPathSource = VoreenApplication::app()->getTemporaryPath(flowParametrization->getName()) + "/";
    tgt::FileSystem::createDirectoryRecursive(simulationPathSource);
    std::string simulationPathDest = username_.get() + "@" + clusterAddress_.get() + ":" + simulationPath_.get() + "/" +
                                     simulationType_.get() + "/";

    // Copy simulation geometry.
    // TODO: support changing/multiple geometries.
    tgt::FileSystem::createDirectory(simulationPathSource + "geometry/");
    std::string geometryFilename = simulationPathSource + "geometry/" + "geometry.stl";
    if (!exportGeometryToSTL(geometryData, geometryFilename)) {
        LERROR("Could not write geometry file");
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
            velocityFile.write(static_cast<const char*>(volume->getData()), volume->getNumBytes());
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

    // Copy data to cluster.
    std::string command = "scp -r " + simulationPathSource + " " + simulationPathDest;
    int ret = executeCommand(command);
    if (ret != EXIT_SUCCESS) {
        VoreenApplication::app()->showMessageBox("Error", "Data could not be copied!", true);
        LERROR("Could not fetch results");
        return;
    }

    // Delete local data.
    tgt::FileSystem::deleteDirectoryRecursive(simulationPathSource);

    // Enqueue simulations.
    for(size_t i=0; i<flowParametrization->size(); i++) {

        std::string enqueueScriptFilename = VoreenApplication::app()->getTemporaryPath("submit.sh");
        std::ofstream enqueueScriptFile(enqueueScriptFilename);
        if (!enqueueScriptFile.good()) {
            LERROR("Could not write enqueue script file");
            continue;
        }
        std::string localSimulationPath = simulationPath_.get() + "/" + simulationType_.get() + "/" +
                                          flowParametrization->getName() + "/" + flowParametrization->at(i).getName();
        enqueueScriptFile << generateEnqueueScript(localSimulationPath);
        enqueueScriptFile.close();

        // Enqueue job.
        command = "ssh " + username_.get() + "@" + clusterAddress_.get() + " 'bash -s' < " + enqueueScriptFilename;
        ret = executeCommand(command);
        if (ret != EXIT_SUCCESS) {
            LERROR("Could not enqueue job");
            continue;
        }
    }

    progress_.setProgress(1.0f);

#ifdef VRN_MODULE_ENSEMBLEANALYSIS
    ensemblePort_.clear();
#endif
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

    std::string source = username + "@" + clusterAddress + ":/scratch/tmp/" + username + "/simulations";
    const FlowParametrizationList* parametrizationList = parameterPort_.getData();
    if (parametrizationList) {
        std::string paramPath = "/" + parametrizationList->getName() + "/";
        source += paramPath;
        std::vector<std::string> failed;
        for(const FlowParameters& parameters : parametrizationList->getFlowParametrizations()) {
            std::string command = "scp -r ";
            command += source + parameters.getName() + " ";
            command += simulationResults_.get() + paramPath + parameters.getName();

            int ret = executeCommand(command);
            if (ret != EXIT_SUCCESS) {
                failed.push_back(parametrizationList->getName() + "/" + parameters.getName());
            }
        }

        if(!failed.empty()) {
            VoreenApplication::app()->showMessageBox("Error", "Some data could not be fetched. See log for details.", true);
            LERROR("Could not fetch results of: \n* " << strJoin(failed, "\n* "));
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
        VoreenApplication::app()->showMessageBox("Fetching Done", "Fetching data was successful!");
    }


#if 0//#ifdef VRN_MODULE_ENSEMBLEANALYSIS
    std::unique_ptr<RawVolumeReader> reader(new RawVolumeReader());
    RawVolumeReader::ReadHints hints;
    reader->setReadHints(hints);

    EnsembleDataset* dataset = new EnsembleDataset();

    std::vector<std::string> runs = tgt::FileSystem::listSubDirectories(simulationResults_.get(), true);
    for(const std::string& run : runs) {
        std::string runPath = simulationResults_.get() + "/" + run;
        std::vector<std::string> fileNames = tgt::FileSystem::readDirectory(runPath, true, false);

        std::vector<EnsembleDataset::TimeStep> timeSteps;
        for(const std::string& fileName : fileNames) {
            std::string url = runPath + "/" + fileName;

            EnsembleDataset::TimeStep timeStep;
            timeStep.path_ = fileName;
            timeStep.time_ = 0.0f; // TODO;
            timeStep.duration_ = 0.0f; // TODO:
            bool timeSet = false;

            const std::vector<VolumeURL>& subURLs = reader->listVolumes(url);
            for(const VolumeURL& subURL : subURLs) {
                VolumeBase* volumeHandle = reader->read(subURL);
                if(!volumeHandle)
                    break;

                if (!volumeHandle->hasDerivedData<VolumeMinMax>())
                    LWARNING("Volume does not contain min max information - needs to be calculated");

                // TODO: run simulation on cluster and find out meta data.
                const FloatMetaData* time = dynamic_cast<const FloatMetaData*>(volumeHandle->getMetaData(SIMULATED_TIME_NAME));
                if(!time) {
                    delete volumeHandle;
                    LERROR("Meta data '" << SIMULATED_TIME_NAME << "' not present for " << subURL.getPath());
                    break;
                }

                const MetaDataBase* scalar = volumeHandle->getMetaData(SCALAR_FIELD_NAME);
                if(!scalar) {
                    delete volumeHandle;
                    LERROR("Meta data '" << SCALAR_FIELD_NAME << "' not present for " << subURL.getPath());
                    break;
                }

                // Add additional information gained reading the file structure.
                Volume* volume = dynamic_cast<Volume*>(volumeHandle);
                tgtAssert(volume, "volumeHandle must be volume");
                volume->getMetaDataContainer().addMetaData("run_name", new StringMetaData(run));

                timeStep.channels_[scalar->toString()] = volumeHandle;
            }

            // Calculate duration the current timeStep is valid.
            // Note that the last time step has a duration of 0.
            if (!timeSteps.empty())
                timeSteps.back().duration_ = timeStep.time_ - timeSteps.back().time_;

            timeSteps.push_back(timeStep);
        }

        // Update dataset.
        tgt::vec3 color = ColorPool::getDistinctColor(dataset->getRuns().size());
        dataset->addRun(EnsembleDataset::Run{ run, color, timeSteps });
    }

    ensemblePort_.setData(dataset, true);
#endif
}

int FlowSimulationCluster::executeCommand(const std::string& command) const {

    std::string result;
    std::array<char, 128> buffer;

    auto pipe = popen(command.c_str(), "r");
    if (!pipe) {
        LERROR("Failed to open pipe");
        return EXIT_FAILURE;
    }

    while (!feof(pipe)) {
        if (fgets(buffer.data(), 128, pipe) != nullptr) {
            result += buffer.data();
        }
    }

    if(!result.empty()) {
        LINFO(result);
    }

    return pclose(pipe);
}

std::string FlowSimulationCluster::generateEnqueueScript(const std::string& parametrizationPath) const {
    std::stringstream script;

    script << "#!/bin/bash" << std::endl;
    script << "module add intel/2018a" << std::endl;
    script << "cd " << parametrizationPath << std::endl;
    script << "sbatch submit.cmd" << std::endl;

    return script.str();
}

std::string FlowSimulationCluster::generateSubmissionScript(const std::string& parametrizationName) const {
    tgtAssert(parameterPort_.hasData(), "no data");
    std::stringstream script;

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
    script << "# set the number of CPU cores per node" << std::endl;
    script << "#SBATCH --exclusive" << std::endl;
    script << std::endl;
    script << "# How much memory is needed (per node)" << std::endl;
    script << "#SBATCH --mem=" << configMemory_.get() << "G" << std::endl;
    script << std::endl;
    script << "# set a partition" << std::endl;
    script << "#SBATCH --partition " << configPartition_.get() << std::endl;
    script << std::endl;
    script << "# set max wallclock time" << std::endl;
    script << "#SBATCH --time=00:20:00" << std::endl;
    script << std::endl;
    script << "# set name of job" << std::endl;
    script << "#SBATCH --job-name=" + parameterPort_.getData()->getName() + "-" + parametrizationName << std::endl;
    script << std::endl;
    script << "# mail alert at start, end and abortion of execution" << std::endl;
    script << "#SBATCH --mail-type=ALL" << std::endl;
    script << std::endl;
    script << "# set an output file" << std::endl;
    script << "#SBATCH --output output.dat" << std::endl;
    script << std::endl;
    script << "# send mail to this address" << std::endl;
    script << "#SBATCH --mail-user=" + username_.get() + "@uni-muenster.de" << std::endl;
    script << std::endl;
    script << "# run the application" << std::endl;
    script << "OMP_NUM_THREADS=" << configCPUsPerTask_.get() << " "
           << "mpirun ~/" + simulationPath_.get() + "/" + simulationType_.get() + "/" + simulationType_.get() << " "
           << parametrizationName << std::endl;

    return script.str();
}

}   // namespace
