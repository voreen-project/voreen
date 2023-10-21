/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2021 University of Muenster, Germany,                        *
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

#include "flowsimulationresult.h"

#include "modules/vtk/io/vtmvolumereader.h"
#include "custommodules/flowsimulation/processors/volume/volumemerger.h"

#include <vtkXMLDataParser.h>
#include <vtkXMLUtilities.h>



namespace voreen {



const std::string FlowSimulationResult::loggerCat_("voreen.flowsimulation.FlowSimulationResult");

FlowSimulationResult::FlowSimulationResult()
    : AsyncComputeProcessor()
    , outport_(Port::OUTPORT, "outport", "Outport")
    , inputFile_("inputFile", "Input File", "Load PVD file", "", "*.pvd", FileDialogProperty::OPEN_FILE, VALID, Property::LOD_DEFAULT, VoreenFileWatchListener::OPTIONAL_ON)
    , fields_("fields", "Fields")
    , timeStep_("timestep", "Timestep", 0, 0, 1000, INVALID_RESULT, IntProperty::DYNAMIC)
    , selectMostRecentTimeStep_("selectMostRecentTimeStep", "Select Most Recent Timestep", true)
{
    addPort(outport_);
    addProperty(inputFile_);
    ON_CHANGE(inputFile_, FlowSimulationResult, onFileChange);
    addProperty(fields_);
    addProperty(timeStep_);
    timeStep_.setTracking(false);
    addProperty(selectMostRecentTimeStep_);
}

FlowSimulationResult::~FlowSimulationResult() {
}

bool FlowSimulationResult::isReady() const {
    return true;
}

std::vector<std::string> FlowSimulationResult::loadPvdFile(std::string file) const {
    vtkXMLDataElement* root = vtkXMLUtilities::ReadElementFromFile(file.c_str());
    if (root == nullptr) {
        return {};
    }

    vtkXMLDataElement* collection = root->FindNestedElementWithName("Collection");
    if (collection == nullptr) {
        return {};
    }

    std::vector<std::string> files;
    for (int i = 0; i < collection->GetNumberOfNestedElements(); i++) {
        vtkXMLDataElement* dataSet = collection->GetNestedElement(i);
        if (dataSet) {
            const char* timestepStr = dataSet->GetAttribute("timestep");
            if (timestepStr) {
                int timestep = std::stoi(timestepStr);
                std::cout << "Timestep: " << timestep << std::endl; // TODO: do we need this information?
            }
            const char* fileStr = dataSet->GetAttribute("file");
            if (fileStr) {
                std::string path = tgt::FileSystem::dirName(file) + "/" + fileStr;
                files.emplace_back(path);
            }
        }
    }

    return files;
}

void FlowSimulationResult::onFileChange() {

    auto pvdFile = inputFile_.get();

    auto files = loadPvdFile(pvdFile);

    // In case we have less time steps than before, this indicates that a different simulation has been loaded.
    if (files.size() < timeStepPaths_.size()) {
        timeStep_.set(0);
    }

    timeStepPaths_ = std::move(files);
    timeStep_.setMaxValue(timeStepPaths_.size() - 1);

    if (selectMostRecentTimeStep_.get()) {
        timeStep_.set(timeStep_.getMaxValue());
    }

    // We query the first time step and first block as reference, to extract the fields.
    VolumeURL url(timeStepPaths_[0]);
    url.addSearchParameter("block", "0");
    auto volumes = VTMVolumeReader().listVolumes(url.getURL());

    std::deque<Option<std::string>> options;
    for (const auto& volume : volumes) {
        auto name = volume.getMetaDataContainer().getMetaData("name")->toString();
        options.push_back(Option<std::string>(name, name, name));
    }
    fields_.setOptions(options);
    fields_.invalidate();
}

FlowSimulationResult::ComputeInput FlowSimulationResult::prepareComputeInput() {
    auto timeStep = timeStep_.get();
    auto file = timeStepPaths_[timeStep];
    auto field = fields_.get();

    return {file, field};
}

FlowSimulationResult::ComputeOutput FlowSimulationResult::compute(ComputeInput input, ProgressReporter& progressReporter) const {
    VolumeURL url(input.path);
    url.addSearchParameter("name", input.field);
    std::unique_ptr<VolumeList> volumes(VTMVolumeReader().read(url.getURL()));

    VolumeMerger merger;
    VolumeListPort* port = new VolumeListPort(Port::OUTPORT, "volumelist");
    port->connect(merger.getPort("volumelist.input"));
    port->setData(volumes.get(), false);
    //dynamic_cast<VolumeListPort*>(merger.getPort("volumelist.input"))->setData(volumes.get(), false);
    auto computeInput = merger.prepareComputeInput();
    auto output = merger.compute(std::move(computeInput), progressReporter);

    auto ref = std::move(output.outputVolume);
    // auto ref = std::unique_ptr<VolumeBase>(volumes->first());
    return { std::move(ref) };
}

void FlowSimulationResult::processComputeOutput(ComputeOutput output) {
    outport_.setData(output.volume.release());
}


}   // namespace
