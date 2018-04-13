/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2018 University of Muenster, Germany.                        *
 * Visualization and Computer Graphics Group <http://viscg.uni-muenster.de>        *
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

#include "voreen/core/datastructures/volume/volumebase.h"

#include "voreen/core/datastructures/volume/volumehash.h"
#include "voreen/core/utils/hashing.h"
#include "voreen/core/datastructures/volume/volume.h"
#include "voreen/core/datastructures/volume/volumedisk.h"
#include "voreen/core/datastructures/octree/volumeoctree.h"

namespace voreen {

const std::string VolumeBase::loggerCat_("voreen.VolumeBase");
const std::string VolumeBase::META_DATA_NAME_OFFSET("Offset");
const std::string VolumeBase::META_DATA_NAME_SPACING("Spacing");
const std::string VolumeBase::META_DATA_NAME_TRANSFORMATION("Transformation");
const std::string VolumeBase::META_DATA_NAME_TIMESTEP("Timestep");
const std::string VolumeBase::META_DATA_NAME_MODALITY("Modality");
const std::string VolumeBase::META_DATA_NAME_REAL_WORLD_MAPPING("RealWorldMapping");

//---- constructor / destructor and clone ----

VolumeBase::VolumeBase()
    // set the format strings and OpenGL information to default values that indicate that something went wrong when initializing the subclass (or we have a decorator without a decorated volume)
    : format_("none")
    , baseType_("none")
    , glInternalFormat_(0)
    , glFormat_(static_cast<GLenum>(0))
    , glType_(static_cast<GLenum>(0))
    , dataInvalidator_(*this)
{
    Observable<VolumeObserver>::addObserver(&dataInvalidator_);
    // we only want to register Volume objects, i.e., volumes with actual data
    //if (VolumeMemoryManager::isInited())
    //    VolumeMemoryManager::getRef().registerVolume(this);
}

VolumeBase::~VolumeBase() {
    // Volume and VolumeDecorator already notify the observer!
    //notifyDelete();
    stopRunningThreads();
    clearFinishedThreads();
    clearDerivedData();
    Observable<VolumeObserver>::removeObserver(&dataInvalidator_);
}

Volume* VolumeBase::clone() const {
    VolumeRAM* v = getRepresentation<VolumeRAM>()->clone();
    return new Volume(v, this);
}

//---- hash functionality ----

std::string VolumeBase::getHash() const {
    return getRawDataHash() + "-" + getMetaDataHash();
}

std::string VolumeBase::getRawDataHash() const {
    return getDerivedData<VolumeHash>()->getHash();
}

std::string VolumeBase::getMetaDataHash() const {
    MetaDataContainer metaData;

    std::vector<std::string> keys = getMetaDataKeys();
    for(size_t i=0; i<keys.size(); i++) {
        const MetaDataBase* md = getMetaData(keys[i]);
        if(md) {
            MetaDataBase* cl = md->clone();
            if(cl)
                metaData.addMetaData(keys[i], cl);
            else {
                LERROR("Failed to clone metadata!");
            }
        }
    }

    XmlSerializer s;

    s.serialize("metaData", metaData);
    //p->serializeValue(s);

    std::stringstream stream;
    s.write(stream);

    return VoreenHash::getHash(stream.str());
}

//---- derived data ----

void VolumeBase::addDerivedData(VolumeDerivedData* data) {
    derivedDataMutex_.lock();
    for (std::set<VolumeDerivedData*>::const_iterator it=derivedData_.begin(); it!=derivedData_.end(); ++it) {
        if (typeid(**it) == typeid(*data)) {
            derivedData_.erase(it);
            break;
        }
    }
    derivedData_.insert(data);
    derivedDataMutex_.unlock();
}
//---- data and representations

VolumeRAM* VolumeBase::getSlice(size_t sliceNumber) const {
    // Although we do not use the VolumeMemoryManager directy in this function,
    // it may be locked in functions called here. As we need to lock representationMutex_,
    // and since vmmMutex always has to be locked before the representationMutex (this
    // prevents deadlocks), we lock it now (if available).
    boost::recursive_mutex* vmmMutex = 0;
    if (VolumeMemoryManager::isInited())
        vmmMutex = VolumeMemoryManager::getRef().getMutex();
    VolumeLockGuard vmmGuard(vmmMutex);

    boost::lock_guard<boost::recursive_mutex> lock(representationMutex_);
    tgtAssert(sliceNumber < getDimensions().z, "Invalid slice number");
    if(hasRepresentation<VolumeRAM>()) {
        tgt::svec3 offset(0,0,sliceNumber);
        tgt::svec3 dimensions = getDimensions();
        dimensions.z = 1;
        return getRepresentation<VolumeRAM>()->getSubVolume(dimensions, offset);
    } else if(hasRepresentation<VolumeOctree>()) {
        return getRepresentation<VolumeOctree>()->createSlice(voreen::SliceAlignment::XY_PLANE, sliceNumber);
    } else if(hasRepresentation<VolumeDisk>()) {
        return getRepresentation<VolumeDisk>()->loadSlices(sliceNumber,sliceNumber);
    } else {
        tgtAssert(false, "Found no VolumeRepresentation to read slice data from.");
        return nullptr;
    }
}

//---- other functions ----

tgt::vec3 VolumeBase::getCubeSize() const {
    // this is ok, as it takes into account 0.5 voxel before and after the first and last sample
    tgt::vec3 cubeSize = tgt::vec3(getDimensions()) * getSpacing();
    return cubeSize;
}

tgt::vec3 VolumeBase::getLLF() const {
    // our LLF is 0.5 voxel before the first sample to consistently account for 3d texture sampling and bounding boxes
    return getOffset() - tgt::vec3(0.5) * getSpacing(); // offset points to the first voxel
}

tgt::vec3 VolumeBase::getURB() const {
    // LLF is 0.5 voxel before and cube size contains an additional voxel, so that the URB is now correct
    return getLLF()+getCubeSize();
}

tgt::mat4 VolumeBase::getVoxelToWorldMatrix() const {
    return getPhysicalToWorldMatrix() * getVoxelToPhysicalMatrix();
}

tgt::mat4 VolumeBase::getWorldToVoxelMatrix() const {
    tgt::mat4 result = tgt::mat4::identity;
    getVoxelToWorldMatrix().invert(result);
    return result;
}

tgt::mat4 VolumeBase::getWorldToTextureMatrix() const {
    tgt::mat4 result = tgt::mat4::identity;
    getTextureToWorldMatrix().invert(result);
    return result;
}

tgt::mat4 VolumeBase::getTextureToWorldMatrix() const {
    return getVoxelToWorldMatrix() * getTextureToVoxelMatrix();
}

tgt::mat4 VolumeBase::getVoxelToPhysicalMatrix() const {
    // 1. Multiply by spacing 2. Apply offset
    return tgt::mat4::createTranslation(getOffset()) * tgt::mat4::createScale(getSpacing());
}

tgt::mat4 VolumeBase::getPhysicalToVoxelMatrix() const {
    return tgt::mat4::createScale(1.0f/getSpacing()) * tgt::mat4::createTranslation(-getOffset());
}

tgt::mat4 VolumeBase::getPhysicalToWorldMatrix() const {
    return getMetaDataValue<Mat4MetaData>(META_DATA_NAME_TRANSFORMATION, tgt::mat4::identity);
}

tgt::mat4 VolumeBase::getWorldToPhysicalMatrix() const {
    tgt::mat4 result = tgt::mat4::identity;
    getPhysicalToWorldMatrix().invert(result);
    return result;
}

tgt::mat4 VolumeBase::getTextureToPhysicalMatrix() const {
    return getVoxelToPhysicalMatrix() * getTextureToVoxelMatrix();
}

tgt::mat4 VolumeBase::getPhysicalToTextureMatrix() const {
    return getVoxelToTextureMatrix() * getPhysicalToVoxelMatrix();
}

tgt::mat4 VolumeBase::getTextureToVoxelMatrix() const {
    // translation accounts for OpenGL texture sampling (see getVoxelToTextureMatrix() function)
    return tgt::mat4::createTranslation(tgt::vec3(-0.5f)) * tgt::mat4::createScale(getDimensions());
}

tgt::mat4 VolumeBase::getVoxelToTextureMatrix() const {
    // OpenGL texture sampling puts the texel in the center -> the texture contains 0.5 texel more at the border
    return tgt::mat4::createScale(1.0f/tgt::vec3(getDimensions())) * tgt::mat4::createTranslation(tgt::vec3(0.5f));
}

std::vector<tgt::vec3> VolumeBase::getCubeVertices() const {
    std::vector<tgt::vec3> cubeVertices;
    tgt::vec3 llf = getLLF();
    tgt::vec3 urb = getURB();

    cubeVertices.push_back(tgt::vec3(llf.x, llf.y, urb.z));// llb 0
    cubeVertices.push_back(tgt::vec3(urb.x, llf.y, urb.z));// lrb 1
    cubeVertices.push_back(tgt::vec3(urb.x, urb.y, urb.z));// urb 2
    cubeVertices.push_back(tgt::vec3(llf.x, urb.y, urb.z));// ulb 3

    cubeVertices.push_back(tgt::vec3(llf.x, llf.y, llf.z));// llf 4
    cubeVertices.push_back(tgt::vec3(urb.x, llf.y, llf.z));// lrf 5
    cubeVertices.push_back(tgt::vec3(urb.x, urb.y, llf.z));// urf 6
    cubeVertices.push_back(tgt::vec3(llf.x, urb.y, llf.z));// ulf 7
    return cubeVertices;
}

MeshGeometry VolumeBase::getBoundingBox(bool applyTransformation) const {
    // create cube with correct positions and texture coordinates
    MeshGeometry boundingBox = MeshGeometry::createCube(getLLF(), getURB(), tgt::vec3(0.f), tgt::vec3(1.f));
    // transform to world coordinates
    if (applyTransformation)
        boundingBox.transform(getPhysicalToWorldMatrix());
    return boundingBox;
}

void VolumeBase::stopRunningThreads() {
    // copy set of threads because they are removed from derivedDataThreads when they finish
    derivedDataThreadMutex_.lock();
    std::set<VolumeDerivedDataThreadBase*> copy = derivedDataThreads_;
    derivedDataThreadMutex_.unlock();

    for (std::set<VolumeDerivedDataThreadBase*>::iterator it=copy.begin(); it!=copy.end(); ++it) {
        VolumeDerivedDataThreadBase* tmp = *it;
        if(tmp->isRunning()) {
            tmp->interrupt();
            tmp->join();
        }
    }
}

void VolumeBase::clearFinishedThreads() {
    derivedDataThreadMutex_.lock();
    for (std::set<VolumeDerivedDataThreadBase*>::iterator it=derivedDataThreadsFinished_.begin(); it!=derivedDataThreadsFinished_.end(); ++it) {
        if((*it)->isRunning()) {
            (*it)->interrupt();
            (*it)->join();
        }
        delete *it;
    }
    derivedDataThreadsFinished_.clear();
    derivedDataThreadMutex_.unlock();
}

void VolumeBase::clearDerivedData() {
    notifyPendingDataInvalidation();
    derivedDataMutex_.lock();
    for (std::set<VolumeDerivedData*>::iterator it=derivedData_.begin(); it!=derivedData_.end(); ++it) {
        delete *it;
    }
    derivedData_.clear();
    derivedDataMutex_.unlock();
}

const VolumeURL& VolumeBase::getOrigin() const {
    return origin_;
}

VolumeURL& VolumeBase::getOrigin() {
    return origin_;
}

void VolumeBase::setOrigin(const VolumeURL& origin) {
    origin_ = origin;
}

float VolumeBase::getTimestep() const {
    return getMetaDataValue<FloatMetaData>(META_DATA_NAME_TIMESTEP, 0.0f);
}

std::string VolumeBase::getFormat() const {
    return format_;
}

std::string VolumeBase::getBaseType() const {
    return baseType_;
}

GLint VolumeBase::getOpenGLInternalFormat() const {
    return glInternalFormat_;
}

GLenum VolumeBase::getOpenGLFormat() const {
    return glFormat_;
}

GLenum VolumeBase::getOpenGLType() const {
    return glType_;
}

size_t VolumeBase::getNumChannels() const {
    boost::lock_guard<boost::recursive_mutex> lock(representationMutex_);
    if (getNumRepresentations() > 0) {
        return getRepresentation(0)->getNumChannels();
    }
    else {
        tgtAssert(false, "Volume has no representation!");
        return 0;
    }
}

tgt::svec3 VolumeBase::getDimensions() const {
    boost::lock_guard<boost::recursive_mutex> lock(representationMutex_);
    if (getNumRepresentations() > 0) {
        return getRepresentation(0)->getDimensions();
    }
    else {
        tgtAssert(false, "Volume has no representation!");
        return tgt::svec3::zero;
    }
}

size_t VolumeBase::getNumVoxels() const {
    boost::lock_guard<boost::recursive_mutex> lock(representationMutex_);
    if (getNumRepresentations() > 0) {
        return getRepresentation(0)->getNumVoxels();
    }
    else {
        tgtAssert(false, "Volume has no representation!");
        return 0;
    }
}

size_t VolumeBase::getBytesPerVoxel() const {
    boost::lock_guard<boost::recursive_mutex> lock(representationMutex_);
    if (getNumRepresentations() > 0) {
        return getRepresentation(0)->getBytesPerVoxel();
    }
    else {
        tgtAssert(false, "Volume has no representation!");
        return 0;
    }
}

tgt::vec3 VolumeBase::getSpacing() const {
    return getMetaDataValue<Vec3MetaData>(META_DATA_NAME_SPACING, tgt::vec3(1.0f));
}

RealWorldMapping VolumeBase::getRealWorldMapping() const {
    return getMetaDataValue<RealWorldMappingMetaData>(META_DATA_NAME_REAL_WORLD_MAPPING, RealWorldMapping());
}

tgt::vec3 VolumeBase::getOffset() const {
    return getMetaDataValue<Vec3MetaData>(META_DATA_NAME_OFFSET, tgt::vec3(0.0f));
}

Modality VolumeBase::getModality() const {
    return Modality(getMetaDataValue<StringMetaData, std::string>(META_DATA_NAME_MODALITY, "unknown"));
}

void VolumeBase::notifyDelete() {
    std::vector<VolumeObserver*> observers = Observable<VolumeObserver>::getObservers();
    for (size_t i=0; i<observers.size(); ++i)
        observers[i]->volumeDelete(this);
}

void VolumeBase::notifyChanged() {
    std::vector<VolumeObserver*> observers = Observable<VolumeObserver>::getObservers();
    for (size_t i=0; i<observers.size(); ++i)
        observers[i]->volumeChange(this);
}

/// ----------------------------------------------------------------------------
/// VolumeBaseDataInvalidator --------------------------------------------------
/// ----------------------------------------------------------------------------
VolumeBaseDataInvalidator::VolumeBaseDataInvalidator(VolumeBase& base)
    : base_(base)
{
}

void VolumeBaseDataInvalidator::volumeDelete(const VolumeBase* source) {
    tgtAssert(source == &base_, "Invalid volume observation");
    base_.notifyPendingDataInvalidation();
}
void VolumeBaseDataInvalidator::volumeChange(const VolumeBase* source) {
    tgtAssert(source == &base_, "Invalid volume observation");
    base_.notifyPendingDataInvalidation();
}
void VolumeBaseDataInvalidator::volumeRepresentationDelete(const VolumeBase* source, const VolumeRepresentation* rep) {
    tgtAssert(source == &base_, "Invalid volume observation");
    base_.notifyPendingDataInvalidation();
}
void VolumeBaseDataInvalidator::volumeDataDelete(const VolumeBase* source) {
    tgtAssert(source == &base_, "Invalid volume observation");
    base_.notifyPendingDataInvalidation();
}

} // namespace
