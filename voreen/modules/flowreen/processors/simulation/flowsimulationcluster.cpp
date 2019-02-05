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

#include "modules/hdf5/io/hdf5volumereader.h"
#include "modules/hdf5/io/hdf5volumewriter.h"

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
    , simulationResults_("simulationResults", "Simulation Results", "Simulation Results", VoreenApplication::app()->getTemporaryPath("simulations"), "", FileDialogProperty::DIRECTORY, Processor::VALID)
    , triggerFetchResults_("triggerFetchResults", "Fetch Results", Processor::VALID)
{
    addPort(geometryDataPort_);
    //addPort(measuredDataPort_); // Currently ignored
    addPort(parameterPort_);
#ifdef VRN_MODULE_ENSEMBLEANALYSIS
    addPort(ensemblePort_);
#endif

    addProperty(username_);
    addProperty(clusterAddress_);
    addProperty(simulationResults_);
    addProperty(triggerFetchResults_);
    ON_CHANGE(triggerFetchResults_, FlowSimulationCluster, fetchResults);
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
    const Geometry *geometryData = geometryDataPort_.getData();
    if (!geometryData) {
        LERROR("No simulation geometry");
        return;
    }

    const FlowParametrizationList* flowParametrization = parameterPort_.getData();
    if (!flowParametrization || flowParametrization->empty()) {
        LERROR("No parameterization");
        return;
    }

    std::string username = username_.get();
    std::string clusterAddress = clusterAddress_.get();

    for(const FlowParameters& flowParameters : flowParametrization->getFlowParametrizations()) {

        std::string directory = simulationResults_.get() + "/" + flowParametrization->getName() + "/" + flowParameters.getName();
        if(!tgt::FileSystem::createDirectoryRecursive(directory)) {
            LERROR("Could not create output directory: " << directory);
            continue;
        }

        // 1. Create parameter configuration file and commit to cluster.
        std::string parameterFilename = VoreenApplication::app()->getUniqueTmpFilePath(".txt");
        std::ofstream parameterFile(parameterFilename);
        if (!parameterFile.good()) {
            LERROR("Could not write parameter file");
            continue;
        }
        parameterFile << flowParameters.toCSVString();
        parameterFile.close();

        std::string destination = username + "@" + clusterAddress + ":/scratch/tmp/" + username +
                                  "/simulation_test/parameter.txt";
        std::string command = "scp -r " + parameterFilename + " " + destination;
        int ret = executeCommand(command);
        if (ret != EXIT_SUCCESS) {
            LERROR("Could not copy parameter file");
            continue;
        }

        // 2. Copy simulation geometry.
        std::string geometryFilename = VoreenApplication::app()->getUniqueTmpFilePath(".stl");
        if (!exportGeometryToSTL(geometryData, geometryFilename)) {
            LERROR("Could not write geometry file");
            continue;
        }

        destination = username + "@" + clusterAddress + ":/scratch/tmp/" + username + "/simulation_test/geometry.stl";
        command = "scp -r " + geometryFilename + " " + destination;
        ret = executeCommand(command);
        if (ret != EXIT_SUCCESS) {
            LERROR("Could not copy geometry file");
            continue;
        }

        // 3. execute script
        command = "ssh " + username + "@" + clusterAddress + " 'bash -s' < local_script.sh";
        ret = executeCommand(command);
        if (ret != EXIT_SUCCESS) {
            LERROR("Could not execute script");
            continue;
        }
    }

#ifdef VRN_MODULE_ENSEMBLEANALYSIS
    ensemblePort_.clear();
#endif
}

void FlowSimulationCluster::fetchResults() {
    std::string username = username_.get();
    std::string clusterAddress = clusterAddress_.get();

    // 1. Create parameter configuration file and commit to cluster.
    VoreenApplication::app()->showMessageBox("Information", "Start fetching data, this may take a while...");
    std::string source = username + "@" + clusterAddress + ":/scratch/tmp/" + username +
                              "/simulations";
    std::string command = "scp -r " + source + " " + simulationResults_.get();
    int ret = executeCommand(command);
    if (ret != EXIT_SUCCESS) {
        VoreenApplication::app()->showMessageBox("Error", "Data could not be fetched!", true);
        LERROR("Could not fetch results");
        return;
    }
    VoreenApplication::app()->showMessageBox("Fetching Done", "Fetching data was successful!");


#if 0// #ifdef VRN_MODULE_ENSEMBLEANALYSIS
    VolumeSerializerPopulator populator;

    EnsembleDataset* dataset = new EnsembleDataset();

    std::vector<std::string> runs = tgt::FileSystem::listSubDirectories(simulationResults_.get(), true);
    for(const std::string& run : runs) {
        std::string runPath = simulationResults_.get() + "/" + run;
        std::vector<std::string> fileNames = tgt::FileSystem::readDirectory(runPath, true, false);

        std::vector<EnsembleDataset::TimeStep> timeSteps;
        for(const std::string& fileName : fileNames) {
            std::string url = runPath + "/" + fileName;
            std::vector<VolumeReader*> readers = populator.getVolumeSerializer()->getReaders(url);
            if(readers.empty()) {
                LERROR("No valid volume reader found for " << url);
                break;
            }

            VolumeReader* reader = readers.front();
            tgtAssert(reader, "Reader was null");

            EnsembleDataset::TimeStep timeStep;
            timeStep.path_ = fileName;
            timeStep.time_ = 0.0f;
            timeStep.duration_ = 0.0f;
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

                if (!timeSet) {
                    timeStep.time_ = time->getValue();
                    timeSet = true;
                }
                else if (timeStep.time_ != time->getValue())
                    LWARNING("Meta data '" << SIMULATED_TIME_NAME << "' not equal channel-wise in t=" << timeSteps.size() << " of run" << run);

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

}   // namespace
