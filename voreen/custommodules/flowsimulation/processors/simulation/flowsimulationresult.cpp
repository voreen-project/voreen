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
#include "modules/vtk/io/vtivolumereader.h"

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

    std::vector<std::string> filesToQuery;
    auto baseDir = tgt::FileSystem::dirName(pvdFile);
    auto additionalFiles = tgt::FileSystem::listFiles(baseDir);
    for (auto& file : additionalFiles) {
        if (tgt::FileSystem::fileExtension(file, true) == "vtm") {
            filesToQuery.emplace_back(baseDir + "/" + file);
        }
    }
    filesToQuery.emplace_back(timeStepPaths_[0]);

    std::deque<Option<std::string>> options;
    for (auto& file : filesToQuery) {
        VolumeURL url(file);
        url.addSearchParameter("block", "0");
        auto volumes = VTMVolumeReader().listVolumes(url.getURL());

        for (const auto& volume : volumes) {
            auto name = volume.getMetaDataContainer().getMetaData("name")->toString();
            options.emplace_back(Option<std::string>(name, name, name));
            // TODO: there might be an error here
            if (fieldToPath_.size() < additionalFiles.size()) {
                fieldToPath_[name] = baseDir + "/" + file;
            }
        }
    }

    fields_.setOptions(options);
    fields_.invalidate();
}

std::string FlowSimulationResult::getFileForConfig() const {
    auto path = fieldToPath_.find(fields_.get());
    if(path != fieldToPath_.end()) {
        return path->second;
    }
    else {
        auto timeStep = timeStep_.get();
        auto file = timeStepPaths_[timeStep];
        return file;
    }
}

FlowSimulationResult::ComputeInput FlowSimulationResult::prepareComputeInput() {
    auto file = getFileForConfig();
    auto field = fields_.get();
    return {file, field};
}

#if 0
FlowSimulationResult::ComputeOutput FlowSimulationResult::compute(ComputeInput input, ProgressReporter& progressReporter) const {
    vtkSmartPointer<vtkXMLMultiBlockDataReader> reader = vtkSmartPointer<vtkXMLMultiBlockDataReader>::New();
    reader->SetFileName(input.path.c_str());
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
    VolumeURL url(input.path);
    url.addSearchParameter("name", input.field);
    std::unique_ptr<VolumeList> volumes(VTMVolumeReader().read(url.getURL()));

    VolumeMerger merger;
    merger.setPadding(1);
    merger.setAllowIntersections(true);
    VolumeListPort port(Port::OUTPORT, "volumelist", "", false, Processor::VALID);
    port.setProcessor(const_cast<FlowSimulationResult*>(this));
    port.connect(merger.getPort("volumelist.input"));
    port.setData(volumes.get(), false);

    auto computeInput = merger.prepareComputeInput();
    auto output = merger.compute(std::move(computeInput), progressReporter);

    // Clean up.
    while (!volumes->empty()) {
        delete volumes->first();
    }

    auto volume = std::move(output.outputVolume);
    volume->setSpacing(volume->getSpacing() / VOREEN_LENGTH_TO_SI);
    volume->setOffset(volume->getOffset() / VOREEN_LENGTH_TO_SI);

    return { std::move(volume) };
}
#endif

void FlowSimulationResult::processComputeOutput(ComputeOutput output) {
    outport_.setData(output.volume.release());
}


}   // namespace
