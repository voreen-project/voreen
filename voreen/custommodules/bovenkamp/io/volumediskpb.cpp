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

#include "volumediskpb.h"

#include "voreen/core/datastructures/volume/volumeatomic.h"
#include "voreen/core/datastructures/volume/volumefactory.h"
#include "voreen/core/utils/hashing.h"

namespace voreen {

const std::string VolumeDiskPB::loggerCat_("voreen.bovenkamp.VolumeDiskPB");

VolumeDiskPB::VolumeDiskPB(const std::string& magnitudeFilename,
                           const tgt::bvec3& invertPosition,
                           const tgt::svec3& dimensions,
                           int timeStep)
    : VolumeDisk(VolumeGeneratorFloat().getFormat(), dimensions)
    , invertPosition_(invertPosition)
    , invertVelocity_(false)
    , timeStep_(timeStep)
{
    filenames_.push_back(magnitudeFilename);
}

VolumeDiskPB::VolumeDiskPB(const std::string& velocityXFilename,
                           const std::string& velocityYFilename,
                           const std::string& velocityZFilename,
                           const tgt::bvec3& invertPosition,
                           const tgt::bvec3& invertVelocity,
                           const tgt::svec3& dimensions,
                           int timeStep)
    : VolumeDisk(VolumeGenerator3xFloat().getFormat(), dimensions)
    , invertPosition_(invertPosition)
    , invertVelocity_(invertVelocity)
    , timeStep_(timeStep)
{
    filenames_.push_back(velocityXFilename);
    filenames_.push_back(velocityYFilename);
    filenames_.push_back(velocityZFilename);
}

VolumeDiskPB::~VolumeDiskPB() {
}

std::string VolumeDiskPB::getHash() const {
    std::string configStr;

    for (const std::string& filename : filenames_)
        configStr += filename + "#";
    configStr += std::to_string(timeStep_) + "#";

    configStr += VoreenHash::getHash(&invertPosition_[0], invertPosition_.size) + "#";
    configStr += VoreenHash::getHash(&invertVelocity_[0], invertVelocity_.size) + "#";

    return VoreenHash::getHash(configStr);
}

VolumeRAM* VolumeDiskPB::loadVolume() const {
    return loadBrick(tgt::svec3::zero, dimensions_);
}

VolumeRAM* VolumeDiskPB::loadSlices(const size_t firstZSlice, const size_t lastZSlice) const {
    if (firstZSlice > lastZSlice)
        throw VoreenException("last slice must be behind first slice");

    return loadBrick(tgt::svec3(0, 0, firstZSlice), tgt::svec3(dimensions_.x, dimensions_.y, lastZSlice-firstZSlice+1));
}

VolumeRAM* VolumeDiskPB::loadBrick(const tgt::svec3& offset, const tgt::svec3& dimensions) const {
    // check parameters
    if (tgt::hmul(dimensions) == 0)
        throw VoreenException("requested brick dimensions are zero");
    if (!tgt::hand(tgt::lessThanEqual(offset + dimensions, dimensions_)))
        throw VoreenException("requested brick (at least partially) outside volume dimensions");

    std::unique_ptr<VolumeRAM> result;
    if (getNumChannels() == 1)
        result.reset(VolumeGeneratorFloat().create(dimensions));
    else if (getNumChannels() == 3)
        result.reset(VolumeGenerator3xFloat().create(dimensions));

    tgtAssert(result, "Unhandled channel count");

    for (size_t i = 0; i < filenames_.size(); i++)
        readFile(filenames_[i], result.get(), i, offset, dimensions);
    
    return result.release();
}

size_t toLinear(const tgt::svec3& dim, size_t x, size_t y, size_t z) {
    return dim.x*dim.y*z + dim.x*y + x;
}

void VolumeDiskPB::readFile(const std::string& filename, VolumeRAM* volume, size_t channel, const tgt::svec3& brickOffset, const tgt::svec3& brickDimensions) const {

    // TODO: fix arbitrary slice/brick reading!
    if(brickOffset != tgt::svec3::zero || brickDimensions != dimensions_) {
        throw VoreenException("Reading anything other than the whole volume at once from disk might yield wrong results!");
        //LERROR("Reading anything other than the whole volume at once from disk might yield wrong results!");
    }

    // Open file.
    std::ifstream ifs(filename.c_str());
    if (ifs.fail())
        throw VoreenException("File: '" + filename + "' could not be opened");

    // Seek required offset. This operation is quite expensive.
    const size_t timeStepOffset = brickDimensions.x*brickDimensions.z*timeStep_;
    const size_t xzOffset = brickOffset.x * brickOffset.z;

    for (size_t i = 0; i < timeStepOffset + xzOffset; i++)
        ifs.ignore(std::numeric_limits<std::streamsize>::max(), ifs.widen('\n'));

    //TODO: checks
    const size_t numChannels = getNumChannels();
    float* voxels = reinterpret_cast<float*>(volume->getData());
    std::string tmp;

    for (size_t x = 0; x < brickDimensions.x; x++) {

        // Seek required z offset.
        for (size_t i = 0; i < brickOffset.z; i++)
            ifs.ignore(std::numeric_limits<std::streamsize>::max(), ifs.widen('\n'));

        for (size_t z = 0; z < brickDimensions.z; z++) {

            getline(ifs, tmp);
            std::stringstream line(tmp);

            // Seek required y offset.
            for (size_t i = 0; i < brickOffset.y; i++)
                line.ignore(std::numeric_limits<std::streamsize>::max(), line.widen('\t'));

            for (size_t y = 0; y < brickDimensions.y; y++) {

                // retrieve
                getline(line, tmp, '\t');
                float value = static_cast<float>(atof(tmp.c_str()));

                // Invert value if desired.
                if (numChannels == 3 &&
                  ((channel == 0 && invertVelocity_.x) ||
                  (channel == 1 && invertVelocity_.y) ||
                  (channel == 2 && invertVelocity_.z)))
                    value = -value;

                // Calculate index.
                size_t index = toLinear(brickDimensions,
                    (invertPosition_.x ? brickDimensions.x - 1 - x : x),
                    (invertPosition_.y ? brickDimensions.y - 1 - y : y),
                    (invertPosition_.z ? brickDimensions.z - 1 - z : z))
                    * numChannels + channel;

                // Set actual voxel value.
                voxels[index] = value;
            }
        }
    }

    ifs.close();
}

} // namespace voreen
