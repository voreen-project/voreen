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

#ifndef VRN_VOLUMEURLPROPERTY_H
#define VRN_VOLUMEURLPROPERTY_H

#include "voreen/core/properties/stringproperty.h"
#include "voreen/core/properties/condition.h"
#include "voreen/core/utils/voreenfilewatcher.h"

namespace voreen {

class VolumeBase;
class VolumeInfoProperty;
class ProgressBar;

/**
 * Property encapsulating the source URL of a volume.
 *
 * Additionally, the property is able to load the volume from the assigned URL and
 * and performs the memory management.
 */
class VRN_CORE_API VolumeURLProperty : public StringProperty, public VoreenFileWatchListener {
public:
    /**
     * Constructor.
     *
     * @param ident identifier that is used in serialization
     * @param guiText text that is shown in the gui
     * @param value the initial url to assign
     */
    VolumeURLProperty(const std::string& id, const std::string& guiText, const std::string& url = "",
       int invalidationLevel = Processor::INVALID_PARAMETERS, Property::LevelOfDetail lod = Property::LOD_DEFAULT,
       VoreenFileWatchListener::WatchMode watchMode = VoreenFileWatchListener::OPTIONAL_OFF);
    VolumeURLProperty();

    virtual Property* create() const;

    virtual std::string getClassName() const       { return "VolumeURLProperty"; }
    virtual std::string getTypeDescription() const { return "VolumeURL"; }

    /**
     * Assigns the passed URL to the property, and clears the currently referenced volume.
     * @param url the volume URL to assign, may be empty
     */
    void set(const std::string& url);

    /// @see set
    void setURL(const std::string& url);

    /// Returns the stored source URL, may be empty.
    std::string getURL() const;

    /**
     * Assigns the passed volume to the property
     * and updates the source URL accordingly.
     *
     * @param handle the volume to assign
     * @param owner determines whether the property should take ownership of the passed handle
     *  (i.e., whether the property should delete the volume, when a new one is assigned or loaded)
     */
    void setVolume(VolumeBase* handle, bool owner = false);

    /**
     * Returns the referenced volume, may be null.
     */
    VolumeBase* getVolume() const;

    /**
     * Loads the volume from the currently assigned source URL.
     *
     * @see setURL
     *
     * @note The property takes ownership of the loaded
     *       volume and deletes it on its own
     *       destruction or when a new handle is assigned.
     *
     * @throws FileException, bad_alloc if the volume could not be loaded
     */
    void loadVolume();

    /**
     * Clears the assigned source URL and deletes the
     * corresponding volume, if set and owned by the property.
     */
    void clear();

    /**
     * Adds a VolumeInfoProperty to this. The VolumeInfoProperty will be
     * synchronized with the laoded volume.
     * @see VolumeInfoProperty
     * @param pointer pointer to the property.
     */
    void addInfoProperty(VolumeInfoProperty* pointer);

    /// @see Property::serialize
    virtual void serialize(Serializer& s) const;

    /// @see Property::deserialize
    virtual void deserialize(Deserializer& s);

    /**
     * Inherited from VoreenFileWatchListener.
     * This implementation calls loadVolume()
     */
    virtual void fileActionCallback() override;

protected:
    /// Deletes the assigned volume, if it is owned by the property.
    virtual void deinitialize();

private:
    /// Returns the property's progress bar and generates it on first access.
    ProgressBar* getProgressBar();

    VolumeBase* volume_;             ///< the volume belonging to the assigned URL
    bool volumeOwner_;               ///< determines, if the property owns the volume_
    VolumeInfoProperty* infoProp_;   ///< pointer to an info property, can be NULL
    bool recursiveSet_;              ///< determines if set was called recursively

    ProgressBar* progressBar_;

    static const std::string loggerCat_;
};

} // namespace voreen

#endif // VRN_VOLUMEURLPROPERTY_H
