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

#ifndef VRN_TRANSFUNCPROPERTYBASE_H
#define VRN_TRANSFUNCPROPERTYBASE_H

#include "voreen/core/properties/property.h"
#include "voreen/core/datastructures/volume/volume.h"

namespace voreen {

    class TransFuncBase;

/**
 * Base class of all transfer function properties.
 * It handles the volume interaction and the domain fitting strategy.
 *
 * @NOTE pointers are NOT deleted and deletion has to be done in sub class
 */
class VRN_CORE_API TransFuncPropertyBase : public Property, public VolumeObserver {

    friend class TransFuncPropertyWidgetBase;
    friend class TransFuncPropertyEditorBase;

public:
    /// Automatic domain fitting behavior of the transfunc on volume changes.
    enum DomainAutoFittingStrategy {
        FIT_DOMAIN_NEVER    = 0,    ///< Never automatically fit domain.
        FIT_DOMAIN_INITIAL  = 1,    ///< Fit domain to volume, if the standard domain [0.0:1.0] is currently assigned.
        FIT_DOMAIN_ALWAYS   = 2     ///< Always fit domain to the assigned volume.
    };

    /**
     * Constructor
     *
     * @param ident identifier that is used in serialization
     * @param guiText text that is shown in the gui
     * @param invalidationLevel The owner is invalidated with this InvalidationLevel upon change.
     * @param lod level of detail in the gui representation
     */
    TransFuncPropertyBase(const std::string& ident, const std::string& guiText, int invalidationLevel = Processor::INVALID_RESULT,
                          Property::LevelOfDetail lod = Property::LOD_DEFAULT);
    TransFuncPropertyBase();
    virtual ~TransFuncPropertyBase();

    //-------------------------------------------------
    //  Property Functions
    //-------------------------------------------------
    /** @see Property::serialize */
    virtual void serialize(Serializer& s) const;
    /** @see Property::deserialize */
    virtual void deserialize(Deserializer& s);
    /** @see Property::deserialize */
    virtual void deinitialize();
protected:
    /**
     * Sets the stored transfer function to the given one.
     * @note should only be called from subclasses
     */
    virtual void set(TransFuncBase* tf);
    /** Reset to default function */
    virtual void reset();
    //-------------------------------------------------
    //  Domain Fitting Handling
    //-------------------------------------------------
public:
    /** Determines the domain fitting behavior of the TransFunc on volume changes. */
    void setDomainFittingStrategy(DomainAutoFittingStrategy strategy);
    /** Returns the transfunc's domain fitting behavior. */
    DomainAutoFittingStrategy getDomainFittingStrategy() const;

    //-------------------------------------------------
    //  Volume Handling
    //-------------------------------------------------
public:
    /**
     * Assigns the given volumehandle to this property. It is tested whether the given volumehandle
     * is different to the already stored one. If this is the case, the texture of the transfer
     * function is resized according to the bitdepth of the volume
     *
     * @param handle volumehandle that should be assigned to this property
     * @param channel the channel used by the property
     */
    void setVolume(const VolumeBase* volume, size_t channel = 0);

    /**
     * Returns the volume that is assigned to this property.
     *
     * @return volume that is associated with this property
     */
    const VolumeBase* getVolume() const;

    /**
     * Returns the channel of the volume that is assigned to this property.
     *
     * @return channel of the volume that is associated with this property
     */
    const size_t getVolumeChannel() const;

    /** @see VolumeObserver */
    virtual void volumeDelete(const VolumeBase* source);
    /** @see VolumeObserver */
    virtual void volumeChange(const VolumeBase* source);
    /** @see VolumeObserver */
    virtual void derivedDataThreadFinished(const VolumeBase* source, const VolumeDerivedData* derivedData);

    /**
     * Returns the flag if the histogram of the volume should be computed automatically using a background thread (default is true).
     */
    virtual bool getComputeHistogram() const;

    /**
     * Sets the flag if the histogram of the volume should be computed automatically using a background thread (default is true).
     */
    virtual void setComputeHistogram(bool autoHist = true);

    //-------------------------------------------------
    //  Functions to override
    //-------------------------------------------------
public:
    /** Returns the transfer function. */
    virtual TransFuncBase* get() const = 0;
    /** Sets the tf domain bounds from the current volume. */
    virtual void applyDomainFromData() = 0;
protected:
    /** appls a texture resize according to the requeste size. */
    virtual void applyTextureResizeFromData(int requestedSize) = 0;
    /** Checks, if the domain should be updated from volume */
    virtual bool isDomainFittingEnabled() const = 0;
    /** Checks, if TF meta data is present. */
    virtual bool isTFMetaPresent() const = 0;
    /** Applys TFMeta from data. */
    virtual void applyTFMetaFromData() = 0;

    //-------------------------------------------------
    //  Member
    //-------------------------------------------------
protected:
    const VolumeBase* volume_; ///< volumehandle that is associated with the transfer function property
    size_t channel_;           ///< channel of the volumehandle that is associated with the transfer function property
    bool computeHistogram_;    ///< determines if the histogram of the volumehandle should be requested by the transfer function property

    DomainAutoFittingStrategy domainFittingStrategy_; ///< current used fitting strategy

    TransFuncBase* baseFunction_; ///< pointer to the transfer function of this property
};

} // namespace voreen

#endif // VRN_TRANSFUNCPROPERTYBASE_H
