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

#ifndef VRN_WORKSPACE_H
#define VRN_WORKSPACE_H

#include "voreen/core/network/networkserializer.h"
#include "voreen/core/io/serialization/serialization.h"
#include "voreen/core/network/processornetworkobserver.h"

#include <string>
#include <vector>
#include <map>

namespace tgt {
    class GLCanvas;
}

namespace voreen {

class ProcessorNetwork;
class Animation;
class Workspace;

/**
 * Stores information about the visibility and grouping of a workspace's
 * properties in application mode.
 *
 * The configuration consists of an ordered list of property groups,
 * to which properties are assigned. Each property can only be member
 * of a single group at a time. If a property is not assigned to a group,
 * it is not shown in application mode (default).
 * The property group membership info includes a non-negative priority value
 * that determines the in-group ordering of property widgets in the user interface:
 * a property with a lower priority value is placed above one with a larger value.
 *
 * Also stores the names of all selected/unselected menu entitites.
 *
 * @see PropertyListWidget
 *
 */
class VRN_CORE_API ApplicationModeConfiguration : public Serializable {

public:
    ApplicationModeConfiguration(Workspace* workspace);

    /// Appends a new property group with the passed group name. The name must be unique, but may be empty.
    void addPropertyGroup(const std::string& groupName);

    /// Removes the property group with the specified group name.
    void removePropertyGroup(const std::string& groupName);

    /// Removes the property group at the specified index.
    void removePropertyGroup(size_t index);

    /// Renames the property group at the specified index to the specified name.
    void renamePropertyGroup(size_t index, const std::string& newGroupName);

    /// Swaps the position of the property groups at the specified positions.
    void swapPropertyGroups(size_t posA, size_t posB);

    /// Returns whether a group with the specified name exists.
    bool hasPropertyGroup(const std::string& groupName) const;

    /// Returns the names of all registered property groups.
    std::vector<std::string> getPropertyGroups() const;


    /**
     * Assigns the passed property to the specified group with the specified priority.
     * The property group must exist. A negative priority value causes the property to be removed from the group.
     * If the property is currently assigned to another group, it is implicitly removed from that group.
     */
    void setPropertyGroupMembership(const Property* property, const std::string& groupName, int priority);

    /**
     * Returns the group membership information for the passed property.
     *
     * @return Pair consisting of the group the property belongs to and its in-group priority.
     *         If the property does not belong to any group, <"", -1> is returned.
     */
    std::pair<std::string, int> getPropertyGroupMembership(const Property* property) const;

    /// Returns the properties belonging to the specified group, ordered by their priority.
    std::vector<Property*> getGroupProperties(const std::string& groupName) const;

    /**
     * Returns the priority value of the property in the specified group.
     * If the property does not belong to the group, -1 is returned.
     */
    int getPropertyGroupPriority(const Property* property, const std::string& groupName) const;

    /// Returns whether the passed property is member of any property group, i.e., whether it is shown in application mode.
    bool isPropertyVisible(const Property* property) const;

    /// Removes all properties of the passed processor from the property groups.
    void removeProcessorProperties(const Processor* processor);

    /// Clears all property groups.
    void clearGroups();

    ///
    void setMenuEntityVisibility(std::string& name, bool visible);

    /**
     * checks, if the entity is visible.
     * @ return first value, if entity is has been set. second value is stored visibility. by default false.
     */
    std::pair<bool, bool> isMenuEntityVisible(const std::string& name) const;

    /// returns the main canvas in app mode or ""
    std::string getMainCanvas() const;
    /// sets the main canvas in app mode
    void setMainCanvas(std::string mainCanvas);

    virtual void serialize(Serializer& s) const;
    virtual void deserialize(Deserializer& s);

private:
    std::vector<std::string> propertyGroups_;   ///< names of the property groups

    /// maps from property to <groupName, in-group priority>
    std::map<const Property*, std::pair<std::string, int> > propertyGroupMembership_;

    //map of all un/visible entities
    std::map<std::string, bool> menuEntityVisiblity_;

    //main canvas
    std::string mainCanvas_;

    Workspace* workspace_;
};

//-----------------------------------------------------------------------------

/**
 * Interface for workspace observers. Objects of this type can be registered at a Workspace.
 *
 * @see ProcessorNetworkObserver
 */
class VRN_CORE_API WorkspaceObserver : public Observer {
public:

    /**
     * Is triggered when the application mode configuration of the workspace has changed.
     */
    virtual void applicationModeConfigurationChanged() {};
    /**
     * Is triggered when ever a part of the workspace has been modified
     */
    virtual void workspaceHasBeenModified() {};
};

//-----------------------------------------------------------------------------

/**
 * Marker interface to be implemented by classes that operate on a Workspace object.
 */
class VRN_CORE_API UsesWorkspace {

public:
    virtual ~UsesWorkspace() {};
    virtual void setWorkspace(Workspace* workspace) = 0;

};

//-----------------------------------------------------------------------------

class VRN_CORE_API Workspace : public Serializable, public ProcessorNetworkObserver, public Observable<WorkspaceObserver> {

    friend class ApplicationModeConfiguration;

public:

    /// Current Network version. Used to check for incompatible versions.
    static const int WORKSPACE_VERSION;

    /**
     * Constructor.
     */
    Workspace();

    /**
     * Calls clear().
     */
    ~Workspace();

    /**
     * Deletes the current network and the resource containers,
     * and clears the serialization error collector.
     *
     * The workspace's filename is not cleared.
     */
    void clear();
private:
    /** private function called in clear() and delete*/
    void clearResources();
public:
    /**
     * Updates the workspace from the specified file.
     *
     * Non-fatal errors occurring during workspace load are saved
     * and can be requested using @c getErrors().
     *
     * @param filename the workspace file to load
     * @param workDir Absolute working directory of the workspace, used for relative-to-absolute path conversions.
     *      Is passed as document path to the XMLDeserializer. If an empty string is passed, the location
     *      of the workspace file is used as working directory.
     *
     * @throw SerializationException on serialization errors
     */
    void load(const std::string& filename, const std::string& workDir = "");

    /**
     * Saves the workspace to the specified file.
     *
     * @param filename the file the workspace will be written to
     * @param if true, an existing file will be overwritten, otherwise an exception is thrown when the file already exists
     * @param workDir Absolute working directory of the workspace, used for absolute-to-relative path conversions.
     *      Is passed as document path to the XMLSerializer. If an empty string is passed, the output location
     *      of the workspace file is used as working directory.
     *
     * @throw SerializationException on serialization errors
     */
    void save(const std::string& filename, bool overwrite = true, const std::string& workDir = "");

    /**
     * Returns the errors that have occurred during serialization.
     */
    std::vector<std::string> getErrors() const;

    /// Flag prompting the application to not overwrite the workspace file, not used by serialize()
    bool readOnly() const;

    ProcessorNetwork* getProcessorNetwork() const;
    void setProcessorNetwork(ProcessorNetwork* network);

    void setFilename(const std::string& filename);
    std::string getFilename() const;

    Animation* getAnimation() const;
    void setAnimation(Animation* anim);

    const std::string& getDescription() const;
    bool hasDescription() const;
    void setDescription(const std::string& description);

    const ApplicationModeConfiguration& getApplicationModeConfig() const;
    ApplicationModeConfiguration& getApplicationModeConfig();

    bool isModified() const;
    void setModified(bool isModified = true);

    /// @see Serializable::serialize
    virtual void serialize(Serializer& s) const;

    /// @see Serializable::deserialize
    virtual void deserialize(Deserializer& s);

    //-------------------------------------------------------------------------------
    //  Observer related functions
    //-------------------------------------------------------------------------------
    /// notify all registered workspace observer
    void notifyApplicationModeConfigChanged() const;
    void notifyWorkspaceHasBeenModified() const;
    /// network observer functions
    virtual void networkChanged() { setModified(true); }
private:

    ProcessorNetwork* network_;
    Animation* animation_;
    std::string filename_;
    int version_;
    bool readOnly_;
    std::string description_;
    bool modified_;

    ApplicationModeConfiguration applicationModeConfig_;

    std::vector<std::string> errorList_;

    static const std::string loggerCat_;
};

} // namespace

#endif //VRN_WORKSPACE_H
