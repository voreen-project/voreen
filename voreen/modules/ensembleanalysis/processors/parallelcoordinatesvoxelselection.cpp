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

#include "parallelcoordinatesvoxelselection.h"

#include "modules/ensembleanalysis/utils/ensemblehash.h"
#include "modules/ensembleanalysis/utils/utils.h"

namespace voreen {

const std::string ParallelCoordinatesVoxelSelection::loggerCat_("voreen.EnsembleAnalysis.ParallelCoordinatesVoxelSelection");

ParallelCoordinatesVoxelSelection::ParallelCoordinatesVoxelSelection()
    : AsyncComputeProcessor<ComputeInput, ComputeOutput>()
    , ensembleport_(Port::INPORT, "port_ensemble", "Ensemble Input" )
    , volumeport_(Port::OUTPORT, "port_volume", "Binary Volume Output" )
    , propertySections_("property_sections", "Sections")
    , outputDimensions_("outputDimensions", "Output Dimensions", tgt::ivec3(200), tgt::ivec3(2), tgt::ivec3(1000))
{
    // --- Initialize Ports --- //
    addPort(ensembleport_ );
    addPort(volumeport_ );

    // --- Initialize Properties --- //
    addProperty(propertySections_ );
    propertySections_.setVisibleFlag(false );
    addProperty(outputDimensions_);
}

ParallelCoordinatesVoxelSelectionInput ParallelCoordinatesVoxelSelection::prepareComputeInput() {

    // Clear output first!
    volumeport_.clear();

    PortDataPointer<EnsembleDataset> ensemble = ensembleport_.getThreadSafeData();
    if (!ensemble) {
        throw InvalidInputException("No input", InvalidInputException::S_WARNING);
    }

    auto sectionData = propertySections_.get();
    if(EnsembleHash(*ensemble).getHash() != sectionData.ensembleHash) {
        LERROR("Ensemble hash does not match");
        throw InvalidInputException("Ensemble hash does not match", InvalidInputException::S_IGNORE);
    }

    const EnsembleMember* member = nullptr;
    for( const auto& m : ensemble->getMembers() ) {
        if(m.getName() == sectionData.member ) {
            member = &m;
            break;
        }
    }

    if(!member) {
        if(sectionData.member == "Aggregated") {
            throw InvalidInputException("Selection on Aggregated plots can't be selected, switch to non-aggregated mode", InvalidInputException::S_WARNING);
        }
        else {
            throw InvalidInputException("Could not find member " + sectionData.member, InvalidInputException::S_ERROR);
        }
    }

    // Collect Volumes.
    auto volumes = std::vector<std::pair<const VolumeBase*, int>>();
    for( const auto& field : sectionData.fields ) {
        float time = mapRange(sectionData.time, 0.0f, 1.0f, ensemble->getCommonTimeInterval().x, ensemble->getCommonTimeInterval().y);
        size_t timeStep = member->getTimeStep(time);
        const VolumeBase* volume = member->getTimeSteps()[timeStep].getVolume(field.first);
        if(!volume) {
            std::string url = member->getTimeSteps()[timeStep].getURL(field.first).getURL();
            throw InvalidInputException("Could not load volume " + url, InvalidInputException::S_ERROR);
        }
        volumes.emplace_back(std::pair<const VolumeBase*, int>(volume, field.second));
    }

    tgtAssert(volumes.size() == sectionData.sections.size(), "size mismatch");

    // Bounds have been added later but old behavior was using common bounds.
    tgt::Bounds bounds = sectionData.bounds;
    if(!bounds.isDefined()) {
        bounds = ensemble->getCommonBounds();
    }

    // Create output volume.
    tgt::ivec3 newDims = outputDimensions_.get();
    std::unique_ptr<VolumeRAM_UInt8> outputVolume(new VolumeRAM_UInt8(newDims));
    outputVolume->fill(0);

    // Clear old data.
    volumeport_.clear();

    return ParallelCoordinatesVoxelSelectionInput{
            std::move(ensemble),
            bounds,
            std::move(sectionData),
            std::move(volumes),
            std::move(outputVolume)
    };
}

ParallelCoordinatesVoxelSelectionOutput ParallelCoordinatesVoxelSelection::compute(ParallelCoordinatesVoxelSelectionInput input, ProgressReporter& progress) const {

    auto ensemble = std::move(input.ensemble);
    auto ensembleBounds = input.bounds;
    auto sectionData = std::move(input.sectionData);
    auto volumes = std::move(input.inputVolumes);

    auto output = std::move(input.outputVolume);
    const tgt::svec3 newDims = output->getDimensions();

    // First we consider all voxels.
    // The idea is that we thin out the voxel indices by testing for each of them,
    // if it is within the selected value ranges of each selected (field) volume.
    std::deque<tgt::svec3> voxels;
    for(size_t z=0; z<newDims.z; z++) {
        for(size_t y=0; y<newDims.y; y++) {
            for(size_t x=0; x<newDims.x; x++) {
                voxels.emplace_back(tgt::vec3(x, y, z));
            }
        }
    }

    // In favor of memory caching, we iterate volumes first.
    for(size_t i=0; i<volumes.size(); i++) {

        const VolumeBase* volumeHandle = volumes[i].first;
        int channel = volumes[i].second;

        tgt::Bounds volumeBounds = volumeHandle->getBoundingBox().getBoundingBox();
        RealWorldMapping rwm = volumeHandle->getRealWorldMapping();
        tgt::mat4 worldToVoxel = volumeHandle->getWorldToVoxelMatrix();
        VolumeRAMRepresentationLock volume(volumeHandle);

        const auto& sections = sectionData.sections[i];

        // Progress shall be displayed in detail.
        SubtaskProgressReporter voxelProgress(progress, tgt::vec2(i, i+1) / tgt::vec2(volumes.size()));
        size_t numProcessedVoxels = 0;
        size_t numVoxels = voxels.size();

        for(auto iter = voxels.begin(); iter != voxels.end();) {
            voxelProgress.setProgress(static_cast<float>(numProcessedVoxels++) / static_cast<float>(numVoxels));

            // Map sample position to world space.
            tgt::vec3 pos = mapRange(tgt::vec3(*iter), tgt::vec3::zero, tgt::vec3(newDims), ensembleBounds.getLLF(), ensembleBounds.getURB());

            // Convert voxel to world pos.
            if(!volumeBounds.containsPoint(pos)) {
                iter++; // The voxel is not contained in this volume but might in another.
                continue;
            }

            const auto value = rwm.normalizedToRealWorld(volume->getVoxelNormalized( worldToVoxel * pos, channel ));

            bool selected = sections.empty();
            for( const auto& section : sections ) {
                if( value >= section.x && value <= section.y ) {
                    selected = true;
                    break;
                }
            }

            if( !selected ) {
                iter = voxels.erase(iter);
            }
            else {
                iter++;
            }
        }

        if(voxels.empty()) {
            break;
        }
    }

    // Enable remaining voxels.
    for(const tgt::svec3& p : voxels) {
        output->voxel(p) = std::numeric_limits<uint8_t>::max();
    }

    progress.setProgress(1.0f);

    tgt::vec3 spacing = ensembleBounds.diagonal() / tgt::vec3(newDims);
    std::unique_ptr<Volume> volume(new Volume(output.release(), spacing, ensembleBounds.getLLF()));
    volume->setTimestep(sectionData.time);

    progress.setProgress(1.0f);

    return ParallelCoordinatesVoxelSelectionOutput{
            std::move(volume)
    };
}

void ParallelCoordinatesVoxelSelection::processComputeOutput(ParallelCoordinatesVoxelSelectionOutput output) {
    volumeport_.setData(output.volume.release(), true);
}

}