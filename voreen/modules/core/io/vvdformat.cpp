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

#include "vvdformat.h"

#include "voreen/core/datastructures/volume/volumedisk.h"
#include "voreen/core/datastructures/octree/volumeoctree.h"

#include "voreen/core/datastructures/volume/volumehash.h"
#include "voreen/core/datastructures/volume/volumeminmax.h"
#include "voreen/core/datastructures/volume/volumeminmaxmagnitude.h"
#include "voreen/core/datastructures/volume/histogram.h"
#include "voreen/core/datastructures/volume/volumepreview.h"
#include "voreen/core/datastructures/volume/volumefactory.h"
#include "voreen/core/datastructures/meta/templatemetadata.h"

#include "voreen/core/utils/voreenfilepathhelper.h"

#include "tgt/filesystem.h"

namespace voreen {

const std::string VvdObject::loggerCat_("voreen.datastructure.VvdObject");

//------------------------VvdDataObjectBase-------------------------------------
VvdDataObjectBase::VvdDataObjectBase(const VolumeBase* volume)
    : dimensions_(volume->getDimensions())
    , format_(volume->getFormat())
{}

//------------------------VvdRawDataObject--------------------------------------
VvdRawDataObject::VvdRawDataObject(const VolumeBase* volume, std::string absRawPath)
    : VvdDataObjectBase(volume)
    , absRawPath_(absRawPath)
{}

void VvdRawDataObject::serialize(Serializer& s) const {
    s.serialize("Paths", VoreenFilePathHelper(absRawPath_));
    s.serialize("format", format_);
    int x = dimensions_.x;
    int y = dimensions_.y;
    int z = dimensions_.z;
    s.serialize("x", x);
    s.serialize("y", y);
    s.serialize("z", z);
}

void VvdRawDataObject::deserialize(Deserializer& s) {
    try {
        VoreenFilePathHelper tmp;
        s.deserialize("Paths",tmp);
        absRawPath_ = tmp.getPath();
    } catch (...) {
        // old deserialization
        s.removeLastError();
        s.deserialize("filename", absRawPath_);
        absRawPath_ = tgt::FileSystem::cleanupPath(tgt::FileSystem::dirName(s.getDocumentPath()) + "/" + absRawPath_,true);
    }
    s.deserialize("format", format_);
    int x;
    int y;
    int z;
    s.deserialize("x", x);
    s.deserialize("y", y);
    s.deserialize("z", z);
    dimensions_ = tgt::ivec3(x,y,z);
}

//------------------------VvdRawDataNoExternalFileObject--------------------------------------
VvdRawDataNoExternalFileObject::VvdRawDataNoExternalFileObject(const VolumeBase* volume)
    : VvdDataObjectBase(volume)
    , rawData_(&(reinterpret_cast<const char*>(volume->getRepresentation<VolumeRAM>()->getData())[0]),
               &(reinterpret_cast<const char*>(volume->getRepresentation<VolumeRAM>()->getData())[volume->getBytesPerVoxel()*volume->getNumVoxels()]))
{
}

void VvdRawDataNoExternalFileObject::serialize(Serializer& s) const {
    s.serialize("format", format_);
    int x = dimensions_.x;
    int y = dimensions_.y;
    int z = dimensions_.z;
    s.serialize("x", x);
    s.serialize("y", y);
    s.serialize("z", z);
    s.serializeBinaryBlob("data",rawData_);
}

void VvdRawDataNoExternalFileObject::deserialize(Deserializer& s) {
    s.deserialize("format", format_);
    int x;
    int y;
    int z;
    s.deserialize("x", x);
    s.deserialize("y", y);
    s.deserialize("z", z);
    dimensions_ = tgt::ivec3(x,y,z);
    s.deserializeBinaryBlob("data",rawData_);
}


//-----------------------------------------------------------------------------------------------------

VvdObject::~VvdObject() {
    delete representationData_;
    representationData_ = 0;
    // no need to delete: derived data pointers are either copied from existing volume or assigned to new volume
    //while(!derivedData_.empty()) {
        //delete *derivedData_.begin();
        //derivedData_.erase(derivedData_.begin());
    //}
}

VvdObject::VvdObject(const VolumeBase* vb, std::string absRAWPath, bool useExternalFile)
    : representationData_(0)
{
    tgtAssert(vb,"no volume");
    //create representation data
    bool error = true;
        //try raw
    if(!useExternalFile) {
        representationData_ = new VvdRawDataNoExternalFileObject(vb);
        error = false;
    }

    if(error && !absRAWPath.empty() && tgt::FileSystem::isAbsolutePath(absRAWPath) && tgt::FileSystem::fileExtension(absRAWPath) == "raw") {
        representationData_ = new VvdRawDataObject(vb, absRAWPath);
        error = false;
    }
        //try octree
    if(error && vb->hasRepresentation<VolumeOctree>()) {
        const VolumeOctree* octree = vb->getRepresentation<VolumeOctree>();
        // Why even create that octree?!

        error = false;
    }
        //error
    if(error)
        throw std::invalid_argument("absRawPath does not exist and VolumeBase has no VolumeOctree representation");

    //copy meta and derived data
    std::vector<std::string> keys = vb->getMetaDataKeys();
    for(size_t i=0; i<keys.size(); i++) {
        const MetaDataBase* md = vb->getMetaData(keys[i]);
        if(md) {

            if(keys[i] == VolumeBase::META_DATA_NAME_TRANSFORMATION) {
                const TemplateMetaData<tgt::mat4>* matrixMD = dynamic_cast<const TemplateMetaData<tgt::mat4>*>(md);
                if(matrixMD) {
                    // write transformation matrix unless it is the identity matrix (default)
                    tgt::mat4 transformation = matrixMD->getValue();
                    if (transformation == tgt::mat4::createIdentity())
                        continue;
                }
            }

            if(keys[i] == VolumeBase::META_DATA_NAME_MODALITY) {
                const TemplateMetaData<std::string>* stringMD = dynamic_cast<const TemplateMetaData<std::string>*>(md);
                if(stringMD) {
                    // write modality unless it is unknown (default)

                    std::string modality = stringMD->getValue();
                    if (modality == "unknown")
                        continue;
                }
            }

            if(keys[i] == VolumeBase::META_DATA_NAME_TIMESTEP) {
                const TemplateMetaData<float>* floatMD = dynamic_cast<const TemplateMetaData<float>*>(md);
                if(floatMD) {
                    // write timestep unless it is 0.0f (default)
                    float timestep = floatMD->getValue();
                    if (timestep == 0.0f)
                        continue;
                }
            }

            MetaDataBase* cl = md->clone();
            if(cl)
                metaData_.addMetaData(keys[i], cl);
            else {
                std::cout << "Failed to clone!";
            }
        }
    }

    // Be careful to avoid adding nullptr elements in case any of the derived data computations fail
    // (e.g., if the required volume representation is not available).
    //derivedData_.insert(vb->getDerivedData<VolumeHash>());
    if(auto derived = vb->hasDerivedData<VolumeMinMax>()) {
        derivedData_.insert(derived);
    }
    if(auto derived = vb->hasDerivedData<VolumeMinMaxMagnitude>()) {
        derivedData_.insert(derived);
    }
    if(auto derived = vb->hasDerivedData<VolumePreview>()) {
        derivedData_.insert(derived);
    }
    if(auto derived = vb->hasDerivedData<VolumeHistogramIntensity>()) {
        derivedData_.insert(derived);
    }
    if(auto derived = vb->hasDerivedData<VolumeHistogramIntensityGradient>()) {
        derivedData_.insert(derived);
    }
    //TODO: removed hard-coded classes
}

VvdObject::VvdObject(const VvdObject& other) {
    metaData_ = other.metaData_;
    derivedData_ = other.derivedData_;

    if(VvdRawDataNoExternalFileObject* raw = dynamic_cast<VvdRawDataNoExternalFileObject*>(other.representationData_)) {
        representationData_ = new VvdRawDataNoExternalFileObject(*raw);
    } else
        if(VvdRawDataObject* raw = dynamic_cast<VvdRawDataObject*>(other.representationData_)) {
            representationData_ = new VvdRawDataObject(*raw);
        } else
            if(VvdOctreeDataObject* octree = dynamic_cast<VvdOctreeDataObject*>(other.representationData_)) {
                representationData_ = new VvdOctreeDataObject(*octree);
            } else {
                LERROR("Unknown VvdDataObject class!");
                tgtAssert(false,"Unknown VvdDataObject class!");
            }
}

Volume* VvdObject::createVolume() {
    VolumeRepresentation* volumeRep = 0;
    if(VvdRawDataNoExternalFileObject* raw = dynamic_cast<VvdRawDataNoExternalFileObject*>(representationData_)) {
        VolumeFactory vf;
        volumeRep = vf.create(raw->getFormat(), raw->getDimensions());
        std::copy(raw->rawData_.begin(), raw->rawData_.end(),reinterpret_cast<char*>(static_cast<VolumeRAM*>(volumeRep)->getData()));
    } else
        if(VvdRawDataObject* raw = dynamic_cast<VvdRawDataObject*>(representationData_)) {
            volumeRep = new VolumeDiskRaw(raw->getRawPath(), raw->getFormat(), raw->getDimensions());
        } else
            if(VvdOctreeDataObject* octree = dynamic_cast<VvdOctreeDataObject*>(representationData_)) {

            } else {
                LERROR("Unknown VvdDataObject class!");
                tgtAssert(false,"Unknown VvdDataObject class!");
            }

    Volume* volume = new Volume(volumeRep, &metaData_, derivedData_);
    return volume;
}


void VvdObject::serialize(Serializer& s) const {

    if(VvdRawDataObject* raw = dynamic_cast<VvdRawDataObject*>(representationData_))
        s.serialize("RawData", raw);
    else
        if(VvdRawDataNoExternalFileObject* raw = dynamic_cast<VvdRawDataNoExternalFileObject*>(representationData_))
            s.serialize("RawDataNoExternalFile", raw);
        else
            if(VvdOctreeDataObject* octree = dynamic_cast<VvdOctreeDataObject*>(representationData_))
                s.serialize("OctreeData", octree);
            else {
                LERROR("Unknown VvdDataObject class!");
                tgtAssert(false,"Unknown VvdDataObject class!");
            }

    //serialize meta and derived data
    metaData_.serialize(s);
    tgtAssert(derivedData_.count(nullptr) == 0, "derivedData_ containers nullptr");
    s.serialize("DerivedData", derivedData_, "DerivedItem");
}

void VvdObject::deserialize(Deserializer& s) {
    try {
        VvdRawDataObject* tmp = 0;
        s.deserialize("RawData", tmp);
        representationData_ = tmp;
    } catch (...) {
        s.removeLastError();
    }

    try {
        VvdRawDataNoExternalFileObject* tmp = 0;
        s.deserialize("RawDataNoExternalFile", tmp);
        representationData_ = tmp;
    } catch (...) {
        s.removeLastError();
    }

    try {
        VvdOctreeDataObject* tmp = 0;
        s.deserialize("OctreeData", tmp);
        representationData_ = tmp;
    } catch (...) {
        s.removeLastError();
    }

    metaData_.deserialize(s);

    //renaming if old vvd files are used
    metaData_.renameMetaData("transformation", VolumeBase::META_DATA_NAME_TRANSFORMATION);
    metaData_.renameMetaData("timestep", VolumeBase::META_DATA_NAME_TIMESTEP);
    metaData_.renameMetaData("modality", VolumeBase::META_DATA_NAME_MODALITY);
    metaData_.renameMetaData("spacing", VolumeBase::META_DATA_NAME_SPACING);
    metaData_.renameMetaData("offset", VolumeBase::META_DATA_NAME_OFFSET);


    try {
        s.deserialize("DerivedData", derivedData_, "DerivedItem");
    }
    catch (SerializationNoSuchDataException& /*e*/) {
        s.removeLastError();
    }
    // Handle invalid vvd files gracefully:
    derivedData_.erase(nullptr);
}

} // namespace voreen

