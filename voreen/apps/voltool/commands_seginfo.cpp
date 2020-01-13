/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2020 University of Muenster, Germany,                        *
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

#include "commands_seginfo.h"
#include "voreen/core/io/volumeserializer.h"
#include "voreen/core/io/volumeserializerpopulator.h"

using tgt::ivec3;
using tgt::vec3;

namespace voreen {

CommandSegInfo::CommandSegInfo() /*:
    Command("--seginfo", "", "", "<IN OUT>", 2) */
{
//    loggerCat_ += "." + name_;
}

bool CommandSegInfo::checkParameters(const std::vector<std::string>& parameters) {
    return (parameters.size() == 2);
}

bool CommandSegInfo::execute(const std::vector<std::string>& parameters) {

    VolumeSerializerPopulator volLoadPop;
    const VolumeSerializer* serializer = volLoadPop.getVolumeSerializer();

    //load volume dataset
    const VolumeRAM_UInt8* sourceDataset_;
    VolumeCollection* volumeCollection = serializer->read(parameters[0]);
    sourceDataset_ = dynamic_cast<const VolumeRAM_UInt8*>(volumeCollection->first()->getRepresentation<VolumeRAM>());
    tgt::ivec3 dimensions = sourceDataset_->getDimensions();

    std::map<uint8_t, std::vector<tgt::ivec3> > segments = std::map<uint8_t, std::vector<tgt::ivec3> >();

    for (int voxel_z=0; voxel_z<dimensions.z; voxel_z++) {
        for (int voxel_y=0; voxel_y<dimensions.y; voxel_y++) {
            for (int voxel_x=0; voxel_x<dimensions.x; voxel_x++) {

                tgt::ivec3 pos = tgt::ivec3(voxel_x, voxel_y, voxel_z);
                uint8_t v = sourceDataset_->voxel(pos);
                std::map<uint8_t, std::vector<tgt::ivec3> >::iterator it = segments.find(v);
                if(it == segments.end()) {
                    std::vector<tgt::ivec3> positions = std::vector<tgt::ivec3>();
                    positions.push_back(pos);
                    segments.insert(std::make_pair(v, positions));
                } else {
                    it->second.push_back(pos);
                }
            }
        }
    }

    std::cout << "No of Segments: " << segments.size() << std::endl;

    for(std::map<uint8_t, std::vector<tgt::ivec3> >::iterator it = segments.begin(); it != segments.end(); it++) {
        tgt::dvec3 avg = tgt::dvec3(0.0);
        std::vector<tgt::ivec3> positions = it->second;
        for(size_t i = 0; i < positions.size(); i++) {
            avg += tgt::dvec3(positions[i]);
        }
        std::cout << "Segment " << (int)it->first << ": " << std::endl;
        std::cout << "\tSize " << positions.size() << std::endl;
        std::cout << "\tPercent " << (100.0 * (double)positions.size()) / (double)(dimensions.x*dimensions.y*dimensions.z) << std::endl;
        std::cout << "\tAvg Pos " << avg / (double)positions.size() << std::endl;
        //std::cout << "\tnormalized " << 2.0*(avg / (double)positions.size()) / tgt::max(tgt::dvec3(dimensions)) - tgt::dvec3(1.0) << std::endl;
        std::cout << "\tnormalized " << 2.0*(avg / (double)positions.size()) / tgt::dvec3(dimensions) - tgt::dvec3(1.0) << std::endl;
    }

    return true;
}

/////////////////////////////////////////////////////////////////////////////7

CommandCutToSegs::CommandCutToSegs() /*:
    Command("--segcut", "", "", "<INORIG INSEG OUT>", 3) */
{
//    loggerCat_ += "." + name_;
}

bool CommandCutToSegs::checkParameters(const std::vector<std::string>& parameters) {
    return (parameters.size() == 3);
}

bool CommandCutToSegs::execute(const std::vector<std::string>& parameters) {

    VolumeSerializerPopulator volLoadPop;
    const VolumeSerializer* serializer = volLoadPop.getVolumeSerializer();

    //load original dataset
    const VolumeRAM* origDataset_;
    VolumeCollection* volumeCollection = serializer->read(parameters[0]);
    origDataset_ = volumeCollection->first()->getRepresentation<VolumeRAM>();
    VolumeRAM* targetDataset = origDataset_->clone();

    //load segmentation dataset
    const VolumeRAM_UInt8* segDataset_;
    volumeCollection = serializer->read(parameters[1]);
    segDataset_ = dynamic_cast<const VolumeRAM_UInt8*>(volumeCollection->first()->getRepresentation<VolumeRAM>());
    tgt::ivec3 dimensions = segDataset_->getDimensions();

    for (int voxel_z=0; voxel_z<dimensions.z; voxel_z++) {
        for (int voxel_y=0; voxel_y<dimensions.y; voxel_y++) {
            for (int voxel_x=0; voxel_x<dimensions.x; voxel_x++) {

                tgt::ivec3 pos = tgt::ivec3(voxel_x, voxel_y, voxel_z);
                uint8_t v = segDataset_->voxel(pos);
                if(v == 0)
                    targetDataset->setVoxelNormalized(0.f, pos);
            }
        }
    }

    Volume h(targetDataset, volumeCollection->first());
    serializer->write(parameters.back(), &h);
    //delete targetDataset;
    delete origDataset_;

    return true;
}

/////////////////////////////////////////////////////////////////////////////7

CommandCubify::CommandCubify() /*:
    Command("--cubify", "", "", "<IN OUT>", 2) */
{
//   loggerCat_ += "." + name_;
}

bool CommandCubify::checkParameters(const std::vector<std::string>& parameters) {
    return (parameters.size() == 2);
}

bool CommandCubify::execute(const std::vector<std::string>& parameters) {

    VolumeSerializerPopulator volLoadPop;
    const VolumeSerializer* serializer = volLoadPop.getVolumeSerializer();

    //load original dataset
    const VolumeRAM* origDataset;
    VolumeCollection* volumeCollection = serializer->read(parameters[0]);
    origDataset = volumeCollection->first()->getRepresentation<VolumeRAM>();

    tgt::ivec3 oldims = origDataset->getDimensions();
    tgt::ivec3 newdims = tgt::ivec3(tgt::max(oldims));
    tgt::ivec3 llf = (newdims - oldims) / 2;
    tgt::ivec3 urb = (newdims + oldims) / 2;

    VolumeRAM* targetDataset;
    if(origDataset->getBitsAllocated() == 8)
        targetDataset = new VolumeRAM_UInt8(newdims);
    else if(origDataset->getBitsAllocated() == 16)
        targetDataset = new VolumeRAM_UInt16(newdims);
    else
        return false;

    for (int voxel_z=0; voxel_z < newdims.z; voxel_z++) {
        for (int voxel_y=0; voxel_y < newdims.y; voxel_y++) {
            for (int voxel_x=0; voxel_x < newdims.x; voxel_x++) {

                tgt::ivec3 pos = tgt::ivec3(voxel_x, voxel_y, voxel_z);

                if(tgt::hor(tgt::lessThan(pos, llf)) || tgt::hor(tgt::greaterThanEqual(pos, urb)))
                    targetDataset->setVoxelNormalized(0.f, pos);
                else {
                    tgt::ivec3 oldPos = pos - llf;
                    targetDataset->setVoxelNormalized(origDataset->getVoxelNormalized(oldPos), pos);
                }
            }
        }
    }

    Volume h(targetDataset, volumeCollection->first()->getSpacing(), vec3(0.0f));
    oldVolumePosition(&h);
    serializer->write(parameters.back(), &h);
    delete origDataset;
    //delete targetDataset;

    return true;
}

} //namespace

