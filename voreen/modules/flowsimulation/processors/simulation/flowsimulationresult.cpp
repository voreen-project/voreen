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

#include "flowsimulationresult.h"

#include "modules/flowsimulation/processors/volume/volumemerger.h"
#include "modules/vtk/io/vtivolumereader.h"
#include "modules/vtk/io/vtmvolumereader.h"

#include <vtkXMLDataParser.h>
#include <vtkXMLUtilities.h>
#include <vtkSmartPointer.h>
#include <vtkMultiBlockDataSet.h>
#include <vtkXMLMultiBlockDataReader.h>
#include <vtkImageData.h>
#include <vtkImageResample.h>
#include <vtkCompositeDataIterator.h>
#include <vtkImageAppend.h>


namespace voreen {



const std::string FlowSimulationResult::loggerCat_("voreen.flowsimulation.FlowSimulationResult");

FlowSimulationResult::FlowSimulationResult()
    : AsyncComputeProcessor()
    , outport_(Port::OUTPORT, "outport", "Outport")
    , inputFile_("inputFile", "Input File", "Load PVD file", "", "*.pvd", FileDialogProperty::OPEN_FILE, VALID, Property::LOD_DEFAULT, VoreenFileWatchListener::OPTIONAL_ON)
    , fields_("fields", "Fields")
    , timeStep_("timestep", "Timestep", 0, 0, 1, INVALID_RESULT, IntProperty::DYNAMIC)
    , selectMostRecentTimeStep_("selectMostRecentTimeStep", "Select Most Recent Timestep", true)
{
    addPort(outport_);
    addProperty(inputFile_);
    ON_CHANGE(inputFile_, FlowSimulationResult, onFileChange);
    addProperty(fields_);
    ON_CHANGE_LAMBDA(fields_, [this] { selectedField_ = fields_.get(); });
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

    auto dirName = tgt::FileSystem::dirName(file);

    std::vector<std::string> files;
    for (int i = 0; i < collection->GetNumberOfNestedElements(); i++) {
        vtkXMLDataElement* dataSet = collection->GetNestedElement(i);
        if (dataSet) {
            /*
            const char* timestepStr = dataSet->GetAttribute("timestep");
            if (timestepStr) {
                int timestep = std::stoi(timestepStr);
                std::cout << "Timestep: " << timestep << std::endl; // TODO: do we need this information?
            }
            */
            const char* fileStr = dataSet->GetAttribute("file");
            if (fileStr) {
                std::string path = dirName + "/" + fileStr;
                files.emplace_back(path);
            }
        }
    }

    return files;
}

void FlowSimulationResult::onFileChange() {

    // Don't do anything unless deserialization has finished.
    if(firstProcessAfterDeserialization()) {
        return;
    }

    fields_.setOptions(std::deque<Option<std::string>>());

    if(!tgt::FileSystem::fileExists(inputFile_.get())) {
        inputFile_.set("");
        fields_.invalidate();
        return;
    }

    auto pvdFile = inputFile_.get();

    auto files = loadPvdFile(pvdFile);
    if (files.empty()) {
        LWARNING("File contains no time steps");
    }

    // In case we have less time steps than before, this indicates that a different simulation has been loaded.
    if (files.size() < timeStepPaths_.size()) {
        timeStep_.set(0);
    }

    timeStepPaths_ = std::move(files);

    if (!timeStepPaths_.empty()) {
        timeStep_.setMinValue(0);
        timeStep_.setMaxValue(timeStepPaths_.size() - 1);
    }

    if (selectMostRecentTimeStep_.get()) {
        // HACK: Since the .pvd file might be updated before all volumes are written to disk
        //  this might crash. As a workaround, we select the second to last time step.
        int timeStep = std::max(timeStep_.getMaxValue() - 1, 0);
        timeStep_.set(timeStep_.getMaxValue());
    }

    auto queryFieldsFromFile = [&] (const std::string& file) {
        std::vector<std::string> fields;

        VolumeURL url(file);
        url.addSearchParameter("block", "0");
        auto volumes = VTMVolumeReader().listVolumes(url.getURL());

        for (const auto& volume : volumes) {
            auto name = volume.getMetaDataContainer().getMetaData("name")->toString();
            fields.push_back(name);
        }

        return fields;
    };

    // Begin editing options.
    fields_.blockCallbacks(true);

    // We probe all vtm files in that directory, such as cuboid and material file.
    // These are not mentioned in the pvd file.
    auto baseDir = tgt::FileSystem::dirName(pvdFile);
    auto additionalFiles = tgt::FileSystem::listFiles(baseDir);
    for (auto& file : additionalFiles) {
        if (tgt::FileSystem::fileExtension(file, true) == "vtm") {
            auto fullPath = baseDir + "/" + file;
            auto staticFields = queryFieldsFromFile(fullPath);
            for (const auto& field : staticFields) {
                fields_.addOption(field, field, fullPath);
            }
        }
    }

    // We also probe the first time step contained in the pvd file, as it contains all the fields.
    if (!timeStepPaths_.empty()) {
        auto timeVaryingFields = queryFieldsFromFile(timeStepPaths_[0]);
        for (const auto& field: timeVaryingFields) {
            fields_.addOption(field, field, "");
        }
    }

    // End editing options.
    fields_.blockCallbacks(false);

    // Select previously selected field.
    if (fields_.hasKey(selectedField_)) {
        fields_.select(selectedField_);
    }
}

FlowSimulationResult::ComputeInput FlowSimulationResult::prepareComputeInput() {
    if(inputFile_.get().empty()) {
        throw InvalidInputException("No input file specified.", InvalidInputException::S_WARNING);
    }

    std::string path = fields_.getValue();
    if(path.empty()) {

        if(timeStepPaths_.empty()) {
            throw InvalidInputException("No time steps found.", InvalidInputException::S_WARNING);
        }

        auto timeStep = timeStep_.get();
        if(timeStep < 0 || timeStep >= static_cast<int>(timeStepPaths_.size())) {
            throw InvalidInputException("Invalid time step index.", InvalidInputException::S_WARNING);
        }
        path = timeStepPaths_[timeStep];
    }

    VolumeURL url(path);
    url.addSearchParameter("name", fields_.getDescription());

    return { url };
}

#if 0
FlowSimulationResult::ComputeOutput FlowSimulationResult::compute(ComputeInput input, ProgressReporter& progressReporter) const {
    vtkSmartPointer<vtkXMLMultiBlockDataReader> reader = vtkSmartPointer<vtkXMLMultiBlockDataReader>::New();
    reader->SetFileName(input.url.getPath().c_str());
    reader->Update();

    vtkSmartPointer<vtkImageAppend> imageAppend = vtkSmartPointer<vtkImageAppend>::New();

    vtkMultiBlockDataSet* blockData = vtkMultiBlockDataSet::SafeDownCast(reader->GetOutput());
    vtkCompositeDataIterator* iter = blockData->NewIterator();
    for (iter->InitTraversal(); !iter->IsDoneWithTraversal(); iter->GoToNextItem()) {
        vtkImageData* block = vtkImageData::SafeDownCast(iter->GetCurrentDataObject());
        if(!block) // We only can handle image data.
            continue;

        imageAppend->AddInputData(block);
    }
    iter->Delete();
    imageAppend->Update();

    // Step 3: Optional - Resample the Data
    vtkSmartPointer<vtkImageResample> resampler = vtkSmartPointer<vtkImageResample>::New();
    resampler->SetInputData(imageAppend->GetOutput());
    resampler->SetInterpolationModeToLinear();  // Set your preferred interpolation mode
    //resampler->SetOutputSpacing(newSpacing);    // Set your desired spacing
    //resampler->SetOutputOrigin(newOrigin);      // Set your desired origin
    int extend[] = {100, 100, 100};
    resampler->SetOutputExtent(extend);      // Set your desired extent
    resampler->Update();
    vtkSmartPointer<vtkImageData> resampledImageData = resampler->GetOutput();

    auto volume = std::unique_ptr<Volume>(createVolumeFromVtkImageData(input.path, resampledImageData));

    return { std::move(volume) };
}

#else

FlowSimulationResult::ComputeOutput FlowSimulationResult::compute(ComputeInput input, ProgressReporter& progressReporter) const {

    std::unique_ptr<VolumeList> volumes;
    try {
        volumes = std::unique_ptr<VolumeList>(VTMVolumeReader().read(input.url.getURL()));
    }
    catch (std::exception& e) {
        LERRORC("FlowSimulationResult", e.what());
        return {nullptr};
    }

    VolumeMerger merger;
    merger.setPadding(1);
    merger.setIntersectionResolutionStrategy(VolumeMerger::IRS_LAST);
    VolumeListPort port(Port::OUTPORT, "volumelist", "", false, Processor::VALID);
    port.setProcessor(const_cast<FlowSimulationResult*>(this));
    port.connect(merger.getPort("volumelist.input"));
    port.setData(volumes.get(), false);

    std::unique_ptr<Volume> result;
    try {
        auto computeInput = merger.prepareComputeInput();
        auto output = merger.compute(std::move(computeInput), progressReporter);

        result = std::move(output.outputVolume);
        result->setSpacing(result->getSpacing());
        result->setOffset(result->getOffset());
    }
    catch (InvalidInputException& e) {
        LERRORC("VolumeMerger", e.what());
        return {nullptr};
    }

    // Clean up.
    while (!volumes->empty()) {
        delete volumes->first();
    }

    return { std::move(result) };
}
#endif

void FlowSimulationResult::processComputeOutput(ComputeOutput output) {
    outport_.setData(output.volume.release());
}

void FlowSimulationResult::serialize(Serializer& s) const {
    AsyncComputeProcessor::serialize(s);
    s.serialize("selectedField", selectedField_);
}

void FlowSimulationResult::deserialize(Deserializer& d) {
    d.optionalDeserialize("selectedField", selectedField_, {});
    AsyncComputeProcessor::deserialize(d);
}

}   // namespace
