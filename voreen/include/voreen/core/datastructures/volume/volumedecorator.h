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

#ifndef VRN_VOLUMEHANDLEDECORATOR_H
#define VRN_VOLUMEHANDLEDECORATOR_H

#include "voreen/core/datastructures/volume/volumebase.h"

#include "voreen/core/datastructures/meta/templatemetadata.h"
#include "voreen/core/datastructures/meta/realworldmappingmetadata.h"

#include "tgt/matrix.h"
#include "tgt/vector.h"

namespace voreen {

/**
 * Basic dectorator for volumes. Implemented according to decorator pattern: it is derived from VolumeBase class and contains the decorated VolumeBase as a member (using a pointer).
 */
class VRN_CORE_API VolumeDecoratorIdentity : public VolumeBase, public VolumeObserver {
public:
    VolumeDecoratorIdentity(const VolumeBase* vhb);

    virtual ~VolumeDecoratorIdentity();

    /// returns the decorated VolumeBase object
    virtual const VolumeBase* getDecorated() const;


    //---- basic volume information ----

    virtual std::string getFormat() const;
    virtual std::string getBaseType() const;
    virtual GLint getOpenGLInternalFormat() const;
    virtual GLenum getOpenGLFormat() const;
    virtual GLenum getOpenGLType() const;


    //---- representations (overloaded to work on decorated volume) ----

    // FIXME: cannot overload these methods as they are virtual :( however, they only call virtual methods internally and should work without problems
    /*template <class T>
    bool hasRepresentation() const;

    template <class T>
    const T* getRepresentation() const;

    template <class T>
    bool canConvertToRepresentation() const;

    template <class T>
    void removeRepresentation();*/

    virtual size_t getNumRepresentations() const;

    virtual const VolumeRepresentation* useConverter(const RepresentationConverterBase* converter) const;

    virtual void addRepresentation(VolumeRepresentation* rep);


    //---- meta data ----

    virtual std::vector<std::string> getMetaDataKeys() const;

    virtual const MetaDataBase* getMetaData(const std::string& key) const;

    virtual bool hasMetaData(const std::string& key) const;

    virtual Modality getModality() const;


    //VolumeObserver implementation:
    virtual void volumeDelete(const VolumeBase* source);

    virtual void volumeChange(const VolumeBase* /*source*/);

    virtual void volumeRepresentationDelete(const VolumeBase* volume, const VolumeRepresentation* rep);

    virtual void volumeDataDelete(const VolumeBase* source);

protected:

    virtual const VolumeRepresentation* getRepresentation(size_t i) const;

    virtual void removeRepresentation(size_t i);

    const VolumeBase* base_;
};

//-------------------------------------------------------------------------------------------------

/// Decorates a Volume, replacing a metadata item.
class VRN_CORE_API VolumeDecoratorReplace : public VolumeDecoratorIdentity {
public:
    /**
     * Decorates a volumehandle, replacing a metadata item.
     *
     * @param vhb The Volume to decorate
     * @param key Key of the MetaData item to replace.
     * @param value New Value. The decorator takes ownership.
     */
    VolumeDecoratorReplace(const VolumeBase* vhb, const std::string& key, MetaDataBase* value, bool keepDerivedData);
    virtual ~VolumeDecoratorReplace() {
        stopRunningThreads();
        delete value_;
    }

    virtual std::vector<std::string> getMetaDataKeys() const;
    virtual const MetaDataBase* getMetaData(const std::string& key) const;
    virtual bool hasMetaData(const std::string& key) const;

    MetaDataBase* getValue() const;
    void setValue(MetaDataBase* value);

protected:
    std::string key_;
    MetaDataBase* value_;
};

//-------------------------------------------------------------------------------------------------

/// Decorates a Volume, replacing its physical-to-world transformation matrix.
class VRN_CORE_API VolumeDecoratorReplaceTransformation : public VolumeDecoratorReplace {
public:
    VolumeDecoratorReplaceTransformation(const VolumeBase* vhb, tgt::mat4 matrix) :
      VolumeDecoratorReplace(vhb, META_DATA_NAME_TRANSFORMATION, new Mat4MetaData(matrix), true) {}
};

//-------------------------------------------------------------------------------------------------

/// Decorates a Volume, replacing its voxel spacing.
class VRN_CORE_API VolumeDecoratorReplaceSpacing : public VolumeDecoratorReplace {
public:
    VolumeDecoratorReplaceSpacing(const VolumeBase* vhb, tgt::vec3 spacing) :
      VolumeDecoratorReplace(vhb, META_DATA_NAME_SPACING, new Vec3MetaData(spacing), true) {}
};

//-------------------------------------------------------------------------------------------------

/// Decorates a Volume, replacing its offset in physical coordinates.
class VRN_CORE_API VolumeDecoratorReplaceOffset : public VolumeDecoratorReplace {
public:
    VolumeDecoratorReplaceOffset(const VolumeBase* vhb, tgt::vec3 offset) :
      VolumeDecoratorReplace(vhb, META_DATA_NAME_OFFSET, new Vec3MetaData(offset), true) {}
};

//-------------------------------------------------------------------------------------------------

/// Decorates a Volume, replacing its real-world mapping.
class VRN_CORE_API VolumeDecoratorReplaceRealWorldMapping : public VolumeDecoratorReplace {
public:
    VolumeDecoratorReplaceRealWorldMapping(const VolumeBase* vhb, RealWorldMapping rwm) :
      VolumeDecoratorReplace(vhb, META_DATA_NAME_REAL_WORLD_MAPPING, new RealWorldMappingMetaData(rwm), false) {}
};

//-------------------------------------------------------------------------------------------------

/// Decorates a Volume, replacing its timestep.
class VRN_CORE_API VolumeDecoratorReplaceTimestep : public VolumeDecoratorReplace {
public:
    VolumeDecoratorReplaceTimestep(const VolumeBase* vhb, float timestep) :
      VolumeDecoratorReplace(vhb, META_DATA_NAME_TIMESTEP, new FloatMetaData(timestep), true) {}
};


/////////////////////////////////////////////////////////////////////////////////////////////////////
// Template implementations
/////////////////////////////////////////////////////////////////////////////////////////////////////

/*
template <class T>
bool VolumeDecoratorIdentity::hasRepresentation() const {
    if (base_)
        return base_->hasRepresentation<T>();
    else
        return false;
}

template <class T>
const T* VolumeDecoratorIdentity::getRepresentation() const {
    // VolumeDecorator is handled as an individual VolumeBase by the memory manager, so it has to be notified separately
    if (VolumeMemoryManager::isInited())
        VolumeMemoryManager::getRef().notifyUse(const_cast<VolumeDecoratorIdentity*>(this));

    if (base_)
        return base_->getRepresentation<T>();
    else
        return 0;
}

template <class T>
void VolumeDecoratorIdentity::removeRepresentation() {
    stopRunningThreads();
    if (base_)
        const_cast<VolumeBase*>(base_)->removeRepresentation<T>();
    else
        LERROR("Removing representation from VolumeDecorator with no decorated volume not possible!");
}

template <class T>
bool VolumeDecoratorIdentity::canConvertToRepresentation() const {
    if (base_)
        return base_->canConvertToRepresentation<T>();
    else
        return false;
}*/

} // namespace

#endif
