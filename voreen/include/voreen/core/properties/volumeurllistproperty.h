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

#ifndef VRN_VOLUMEURLLISTPROPERTY_H
#define VRN_VOLUMEURLLISTPROPERTY_H

#include "voreen/core/properties/templateproperty.h"
#include "voreen/core/properties/condition.h"
#include "voreen/core/utils/voreenfilewatcher.h"

#include <vector>
#include <map>

namespace voreen {

class VolumeList;
class VolumeBase;
class ProgressBar;

#ifdef DLL_TEMPLATE_INST
template class VRN_CORE_API TemplateProperty<std::vector<std::string> >;
#endif

/**
 * Property that stores a list of volume source URLs.
 *
 * Additionally, the property provides functionality to load volumes from the assigned URLs
 * and to maintain their selection state.
 */
class VRN_CORE_API VolumeURLListProperty : public TemplateProperty<std::vector<std::string> >, public VoreenFileWatchListener {
public:
    /**
     * Constructor.
     *
     * @param ident identifier that is used in serialization
     * @param guiText text that is shown in the gui
     * @param value the initial source URLs to assign
     */
    VolumeURLListProperty(const std::string& id, const std::string& guiText, const std::vector<std::string>& value,
        int invalidationLevel = Processor::INVALID_PARAMETERS, Property::LevelOfDetail lod = Property::LOD_DEFAULT,
        VoreenFileWatchListener::WatchMode watchMode = VoreenFileWatchListener::OPTIONAL_OFF);
    VolumeURLListProperty();

    virtual Property* create() const;

    virtual std::string getClassName() const       { return "VolumeURLListProperty"; }
    virtual std::string getTypeDescription() const { return "URL List"; }

    /**
     * Assigns the passed URL list to the property, and clears the currently loaded volumes.
     * All passed URLs are marked as not selected. No volumes are loaded.
     *
     * @see loadVolumes
     */
    void set(const std::vector<std::string>& URLs);

    /**
     * Assigns the passed URL list to the property, and clears the currently loaded volumes.
     * No volumes are loaded.
     *
     * @param URLs the URLs to assign to the property
     * @param selected if true, the assigned URLs are marked as selected
     *
     * @see loadVolumes
     */
    void setURLs(const std::vector<std::string>& URLs, bool selected = false);

    /**
     * Adds the passed URL, if not already contained by the property.
     *
     * @param url the URL to add to the property
     * @param selected if true, the added URL is marked as selected
     */
    void addURL(const std::string& url, bool selected = false);

    /**
     * Removes the passed URL and deletes the corresponding volume,
     * if the handle is owned by the property.
     */
    void removeURL(const std::string& url);

    /**
     * Returns the stored list of source URLs.
     */
    const std::vector<std::string>& getURLs() const;

    /**
     * Returns whether the property contains the passed URL.
     */
    bool containsURL(const std::string& url) const;

    /**
     * Adds the passed handle's URL to the property and stores the handle.
     * If the handle's URL is already registered and has a volume assigned,
     * the volume is replaced by the passed one.
     *
     * @param handle the volume to add to the property
     * @param owner if true, the property takes ownership of the volume
     * @param selected if true, the handle's URL is marked as selected
     */
    void addVolume(VolumeBase* handle, bool owner = false, bool selected = false);

    /**
     * Removes the passed volume and deletes it, if it is owned by the property.
     * The volume's URL is not removed from the property.
     *
     * @see removeURL
     */
    void removeVolume(VolumeBase* handle);

    /**
     * Returns a collection containing the volumes that have
     * been loaded from the assigned source URLs.
     *
     * @param selectedOnly if true, only the selected volumes are returned
     * @return a new VolumeList, deleting it is up to the caller
     */
    VolumeList* getVolumes(bool selectedOnly = false) const;

    /**
     * Returns the volume corresponding to the passed URL,
     * or null if the volume has not been loaded, yet.
     */
    VolumeBase* getVolume(const std::string& url) const;

    /// Sets the selection state of the passed URL.
    void setSelected(const std::string& url, bool selected);

    /// Sets the selection state of all assigned URLs.
    void setAllSelected(bool selected);

    /// Returns whether the passed URL is selected.
    bool isSelected(const std::string& url) const;

    /**
     * Loads a single volume from passed source URL.
     * If the volume has already been loaded, it is reloaded.
     *
     * @note The passed URL must be contained by the property.
     *
     * @note The property takes ownership of the loaded volume
     *       and deletes it when the URL is removed from the property
     *       or the property is destructed.
     *
     * @param invalidate set to false to prevent loading the volume from invalidating the property, which is necessary during deserialization
     *
     * @throws FileException, bad_alloc if the volume could not be loaded
     */
    void loadVolume(const std::string& url, bool invalidateProperty = true);

    /**
     * Loads volumes from the currently assigned source URLs. Already
     * loaded volumes are ignored.
     *
     * @param selectedOnly if true, only the selected volumes are loaded
     * @param removeOnFailure if true, URLs of volumes that could not be loaded are removed
     * @param forceReload if true, volumes are reloaded even if they have been loaded successfully before
     *
     * @note The property takes ownership of the loaded
     *       volumes and deletes them on its own
     *       destruction or when the corresponding URL is removed.
     */
    void loadVolumes(bool selectedOnly = false, bool removeOnFailure = false, bool forceReload = false);

    /**
     * Clears the URL list and deletes all volumes that are owned
     * by the property.
     */
    void clear();

    /// @see Property::serialize
    virtual void serialize(Serializer& s) const;

    /// @see Property::deserialize
    virtual void deserialize(Deserializer& s);

    /**
     * Inherited from VoreenFileWatchListener.
     * This implementation calls loadVolumes()
     */
    virtual void fileActionCallback() override;

    bool getPreviewsVisible();
    void setPreviewsVisible(bool previewsVisible);

protected:
    /// Clears the property, thereby deleting all volume that are owned by it.
    virtual void deinitialize();

private:
    /// Returns whether the volume corresponding to the passed URL is owned by the property.
    bool isOwner(const std::string& url) const;

    /// Returns the property's progress bar and generates it on first access.
    ProgressBar* getProgressBar();

    std::map<std::string, VolumeBase*> handleMap_;  ///< maps from URL to volume (transient)
    std::map<std::string, bool> selectionMap_;      ///< maps from URL to selection state (persisted)
    std::map<std::string, bool> ownerMap_;          ///< maps from URL to owner state

    bool previewsVisible_;

    ProgressBar* progressBar_;

    static const std::string loggerCat_;
};

} // namespace voreen

#endif // VRN_VOLUMEURLLISTPROPERTY_H
