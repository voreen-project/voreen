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

#ifndef VRN_VOLUMEPORT_H
#define VRN_VOLUMEPORT_H

#include "voreen/core/ports/genericport.h"
#include "voreen/core/datastructures/volume/volume.h"

#include "voreen/core/properties/floatproperty.h"
#include "voreen/core/properties/optionproperty.h"

namespace voreen {

#ifdef DLL_TEMPLATE_INST
template class VRN_CORE_API GenericPort<VolumeBase>;
#endif

class VRN_CORE_API VolumePort : public GenericPort<VolumeBase>, public VolumeObserver {
public:
    // By default, outports have an invalidationLevel of VALID, while inports will have an
    // invalidation level of INVALID_RESULT
    VolumePort(PortDirection direction, const std::string& id, const std::string& guiName = "",
               bool allowMultipleConnections = false);

    VolumePort(PortDirection direction, const std::string& id, const std::string& guiName,
               bool allowMultipleConnections,
               Processor::InvalidationLevel invalidationLevel);



    virtual Port* create(PortDirection direction, const std::string& id, const std::string& guiName = "") const {return new VolumePort(direction,id,guiName);}
    virtual std::string getClassName() const {return "VolumePort";}
    virtual std::string getContentDescription() const;
    virtual std::string getContentDescriptionHTML() const;

    /**
     * Returns true, if the port contains a Volume object
     * and the Volume has a valid volume.
     */
    virtual bool isReady() const;

    /**
     * Assigns the passed volume to the port.
     * Also removes/adds VolumeObservers to data
     *
     */
    virtual void setData(const VolumeBase* handle, bool takeOwnership = true);

    /**
     * Implementation of VolumeObserver interface.
     */
    virtual void volumeDelete(const VolumeBase* source);

    /**
     * Implementation of VolumeObserver interface.
     *
     * VolumePort needs to overwrite this method, since it should not call volumeDelete (as in the oiginal implementation).
     */
    virtual void volumeDataDelete(const VolumeBase* source);

    /**
     * Implementation of VolumeObserver interface.
     */
    virtual void volumeChange(const VolumeBase* source);

    /**
     * Shows or hides the port's texture access properties in the user interface.
     *
     * These properties (texFilterMode_, texClampMode_, texBorderIntensity_) determine how
     * a processor accesses the volume data it receives via an inport. By default,
     * the texture access properties are hidden.
     *
     * @note Texture access properties can only be shown for inports.
     */
    void showTextureAccessProperties(bool show);

    /**
     * Returns the port's texture filter mode property, which may be used
     * to determine how a processor accesses its volume input data.
     * The property is only available for inports.
     *
     * @see showTextureAccessProperties
     */
    IntOptionProperty& getTextureFilterModeProperty();

    /**
     * Returns the port's texture clamp mode property, which may be used
     * to determine how a processor accesses its volume input data.
     * The property is only available for inports.
     *
     * @see showTextureAccessProperties
     */
    GLEnumOptionProperty& getTextureClampModeProperty();

    /**
     * Returns the port's texture border intensity property, which may be used
     * to determine how a processor accesses its volume input data.
     * The property is only available for inports.
     *
     * @see showTextureAccessProperties
     */
    FloatProperty& getTextureBorderIntensityProperty();

    /// This port type supports caching.
    virtual bool supportsCaching() const;

    /**
     * Returns an MD5 hash of the stored volume,
     * or and empty string, if no volume is assigned.
     */
    virtual std::string getHash() const;

    /**
     * Saves the assigned volume to the given path.
     * If a filename without extension is passed,
     * ".vvd" is appended to it.
     *
     * @throws VoreenException If saving failed or
     *      no volume is assigned.
     */
    virtual void saveData(const std::string& path) const;

    /**
     * Loads a volume from the given path and assigns it
     * to the port. If a filename without extension is passed,
     * ".vvd" is appended to it.
     *
     * @throws VoreenException if loading failed.
     */
    virtual void loadData(const std::string& path);

    virtual tgt::col3 getColorHint() const;

protected:
    IntOptionProperty texFilterMode_;         ///< texture filtering mode to use for volume access
    GLEnumOptionProperty texClampMode_;       ///< texture clamp mode to use for the volume
    FloatProperty texBorderIntensity_;        ///< clamp border intensity
private:
    void construct(std::string id, std::string guiName); // internally used by constructors
};

} // namespace

#endif // VRN_VOLUMEPORT_H
