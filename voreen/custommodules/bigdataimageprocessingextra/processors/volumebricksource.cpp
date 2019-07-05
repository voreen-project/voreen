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

#include "volumebricksource.h"

#include "voreen/core/datastructures/volume/volume.h"
#include "voreen/core/datastructures/volume/volumeram.h"
#include "voreen/core/datastructures/volume/volumedisk.h"
#include "voreen/core/datastructures/octree/volumeoctree.h"
#include "voreen/core/voreenapplication.h"
#include "voreen/core/datastructures/callback/lambdacallback.h"

#include "modules/hdf5/io/hdf5volumereader.h"
#include "modules/hdf5/io/volumediskhdf5.h"

#include "tgt/filesystem.h"

#include <memory>

namespace voreen {

const std::string VolumeBrickSource::loggerCat_("voreen.bigdataimageprocessingextra.volumebricksource");

VolumeBrickSource::VolumeBrickSource()
    : VolumeProcessor()
    , outport_(Port::OUTPORT, "connectedcomponentanalysis.outport", "Volume Output")
    , volumeFilePath_("volumeFilePath", "Volume file input", "Path", "", "HDF5 (*.h5)", FileDialogProperty::OPEN_FILE, Processor::INVALID_RESULT, Property::LOD_DEFAULT)
    , volumeInfo_("volumeInfo","Volume info")
    , volumeDimensions_("volumeDimensions", "Volume dimensions", tgt::ivec3(0), tgt::ivec3(0), tgt::ivec3(std::numeric_limits<int>::max()))
    , brickDimensions_("brickDimensions", "Brick dimensions", tgt::ivec3(1), tgt::ivec3(1), tgt::ivec3(std::numeric_limits<int>::max()))
    , currentBrickOffset_("currentBrickOffset", "Current brick offset", tgt::ivec3(0), tgt::ivec3(0), tgt::ivec3(std::numeric_limits<int>::max()))
    , numBricks_("numBricks", "Number of bricks", tgt::ivec3(1), tgt::ivec3(1), tgt::ivec3(std::numeric_limits<int>::max()))
    , brickToLoad_("brickToLoad", "Brick to load", 0, 0, std::numeric_limits<int>::max())
    , loadButton_("loadButton_", "Load brick", Processor::INVALID_RESULT, Property::LOD_DEFAULT)
    , brickInfo_("brickInfo","Brick info")
    , currentVolume_(nullptr)
{
    addPort(outport_);

    addProperty(volumeFilePath_);
        ON_CHANGE(volumeFilePath_, VolumeBrickSource, tryOpenVolume);

    addProperty(volumeInfo_);

    addProperty(volumeDimensions_);
        volumeDimensions_.setReadOnlyFlag(true);
        ON_CHANGE(volumeDimensions_, VolumeBrickSource, adaptToChangedNumBricks);

    addProperty(brickDimensions_);
        brickDimensions_.setReadOnlyFlag(true);
        ON_CHANGE(brickDimensions_, VolumeBrickSource, updateCurrentBrickOffset);

    addProperty(currentBrickOffset_);
        currentBrickOffset_.setReadOnlyFlag(true);

    addProperty(numBricks_);
        ON_CHANGE(numBricks_, VolumeBrickSource, adaptToChangedNumBricks);

    addProperty(brickToLoad_);
        ON_CHANGE_LAMBDA(brickToLoad_, [this] () {
                setOutput(nullptr);
                });
        ON_CHANGE(brickToLoad_, VolumeBrickSource, updateCurrentBrickOffset);

    addProperty(loadButton_);
        ON_CHANGE(loadButton_, VolumeBrickSource, loadBrick);

    addProperty(brickInfo_);
}

VolumeBrickSource::~VolumeBrickSource() {
}
VoreenSerializableObject* VolumeBrickSource::create() const {
    return new VolumeBrickSource();
}

void VolumeBrickSource::initialize() {
}
void VolumeBrickSource::deinitialize() {
}
tgt::svec3 VolumeBrickSource::getCurrentBrickDimensions() const {
    const tgt::ivec3 offset(currentBrickOffset_.get());
    const tgt::ivec3 maxBrickSize(brickDimensions_.get());
    const tgt::ivec3 volSize(volumeDimensions_.get());

    tgt::ivec3 overrun = offset + maxBrickSize - volSize;
    tgt::ivec3 result = maxBrickSize - tgt::max(tgt::ivec3::zero, overrun);
    tgtAssert(tgt::hand(tgt::greaterThan(result, tgt::ivec3::zero)), "Invalid brick dimensions");
    return tgt::svec3(result);
}

void VolumeBrickSource::process() {
}
void VolumeBrickSource::loadBrick() {
    if(!currentVolume_) {
        LWARNING("No volume file opened");
        return;
    }
    tgt::svec3 currentBrickSize = getCurrentBrickDimensions();
    tgt::svec3 currentBrickOffset = currentBrickOffset_.get();

    VolumeRAM* brick = currentVolumeDisk_->loadBrick(currentBrickOffset, currentBrickSize);
    const tgt::vec3 volSpacing = currentVolume_->getSpacing();
    const tgt::vec3 volOffset = currentVolume_->getOffset();

    // Now create the VolumeDiskHDF5 by handing over fileVolume...
    Volume* vol = new Volume(brick, volSpacing, volOffset + tgt::vec3(currentBrickOffset)*volSpacing);

    // Set real world mapping...
    if(currentVolume_->hasMetaData(VolumeBase::META_DATA_NAME_REAL_WORLD_MAPPING)) {
        // ... from file if available
        vol->setRealWorldMapping(currentVolume_->getRealWorldMapping());
    }  else {
        // ... otherwise a static denomalizing mapping.
        vol->setRealWorldMapping(RealWorldMapping::createDenormalizingMapping(currentVolume_->getBaseType()));
    }

    // Add origin information to the volume.
    vol->setOrigin(volumeFilePath_.get());

    setOutput(vol);
}

void VolumeBrickSource::tryOpenVolume() {
    volumeInfo_.reset();
    setOutput(nullptr);

    std::unique_ptr<HDF5FileVolume> hdf5vol = nullptr;

    // See which volumes are present in the selected file.
    std::vector<VolumeURL> volumes = HDF5VolumeReader().listVolumes(volumeFilePath_.get());
    try {
        if(volumes.empty()) {
            throw tgt::IOException("No volumes in file");
        }
        // If there are any: Pick the first
        VolumeURL url = volumes[0];

        // If the implementation did not change, there should be a 'path' parameter
        tgtAssert(url.getSearchParameter("path") != "", "Could not find path in VolumeURL");

        hdf5vol = HDF5FileVolume::openVolume(url.getPath(), url.getSearchParameter("path"), true /*readonly*/);
    } catch(tgt::IOException& e) {
        VoreenApplication::app()->showMessageBox("Could not open input volume",
                std::string("Could not open input volume:\n\n") + e.what(),
                true);
        currentVolumeDisk_ = nullptr;
        currentVolume_ = nullptr;
        return;
    }
    // Update properties
    numBricks_.setMaxValue(hdf5vol->getDimensions());
    if(Processor::isInitialized()) {
        numBricks_.set(hdf5vol->getDimensions()/hdf5vol->getChunkSize());
    }
    brickToLoad_.setMaxValue(tgt::hmul(numBricks_.get())-1);

    std::unique_ptr<tgt::vec3> offset(hdf5vol->tryReadOffset());
    std::unique_ptr<tgt::vec3> spacing(hdf5vol->tryReadSpacing());
    tgtAssert(hdf5vol, "Volume load failed silently");
    currentVolumeDisk_ = new VolumeDiskHDF5(std::move(hdf5vol));
    currentVolume_ = std::unique_ptr<VolumeBase>(
        new Volume(
            currentVolumeDisk_,
            spacing.get() ? *spacing : tgt::vec3::one,
            offset.get() ? *offset : tgt::vec3::zero)
        );
    currentVolume_->setOrigin(volumeFilePath_.get());

    volumeDimensions_.set(currentVolume_->getDimensions());
    volumeInfo_.setVolume(currentVolume_.get());
}

void VolumeBrickSource::updateBrickDimensions() {
    brickDimensions_.set(tgt::iceil(tgt::vec3(volumeDimensions_.get())/tgt::vec3(numBricks_.get())));
}

void VolumeBrickSource::updateCurrentBrickOffset() {
    int brickId = brickToLoad_.get();
    const tgt::ivec3 numBricks = numBricks_.get();
    tgt::ivec3 brickOffset;
    brickOffset.x = brickId%numBricks.x;
    brickId /= numBricks.x;
    brickOffset.y = brickId%numBricks.y;
    brickId /= numBricks.y;
    brickOffset.z = brickId%numBricks.z;
    brickId /= numBricks.z;

    currentBrickOffset_.set(brickDimensions_.get()*brickOffset);
}
void VolumeBrickSource::setOutput(const VolumeBase* volume) {
    outport_.setData(volume);
    brickInfo_.setVolume(volume);
}
bool VolumeBrickSource::isValidNumBricks(const tgt::ivec3& v) const {
    return isValidNumBricksInDim<0>(v)
        && isValidNumBricksInDim<1>(v)
        && isValidNumBricksInDim<2>(v);
}

template<int dim>
bool VolumeBrickSource::isValidNumBricksInDim(const tgt::ivec3& v) const {
    if(v[dim]<=0) {
        return false;
    }
    if(v[dim]>volumeDimensions_.get()[dim]) {
        return false;
    }
    int nBrickDim = tgt::iceil(static_cast<float>(volumeDimensions_.get()[dim])/v[dim]);
    return nBrickDim*v[dim]-volumeDimensions_.get()[dim] < nBrickDim;
}

template<int dim>
tgt::ivec3 VolumeBrickSource::makeValidNumBricks(const tgt::ivec3& numBricks) const {
    for(int i = 0; i < volumeDimensions_.get()[dim]; ++i) {
        tgt::ivec3 candidate;
        tgt::ivec3 delta(tgt::ivec3::zero);
        delta[dim] = i;
        candidate = numBricks + delta;
        if(isValidNumBricksInDim<dim>(candidate)) {
            return candidate;
        }
        candidate = numBricks - delta;
        if(isValidNumBricksInDim<dim>(candidate)) {
            return candidate;
        }
    }
    tgtAssert(false, "no suitable numBricks");
    return numBricks;
}
void VolumeBrickSource::adaptToChangedNumBricks() {
    tgt::ivec3 numBricks = numBricks_.get();
    if(volumeDimensions_.get() != tgt::ivec3::zero && !isValidNumBricks(numBricks_.get())) {
        numBricks = makeValidNumBricks<0>(numBricks);
        numBricks = makeValidNumBricks<1>(numBricks);
        numBricks = makeValidNumBricks<2>(numBricks);
        numBricks_.set(numBricks);
    }    
    updateBrickDimensions();
    brickToLoad_.setMaxValue(tgt::hmul(numBricks)-1);
    updateCurrentBrickOffset();
}

} // namespace voreen
