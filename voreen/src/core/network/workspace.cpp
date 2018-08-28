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

#include "voreen/core/network/workspace.h"

#include "voreen/core/voreenapplication.h"
#include "voreen/core/network/processornetwork.h"
#include "voreen/core/animation/animation.h"
#include "voreen/core/properties/link/linkevaluatorhelper.h"
#include "voreen/core/animation/animatedprocessor.h"
#include "voreen/core/utils/stringutils.h"

#include "tgt/filesystem.h"
#include "tgt/glcanvas.h"

#ifdef WIN32
#include <windows.h>
#else
#include <unistd.h>
#endif

#include <cstdio>
#include <sstream>

namespace voreen {

//-------------------------------------------------------------------------------------------------

const std::string Workspace::loggerCat_("voreen.Workspace");
const int Workspace::WORKSPACE_VERSION = 2;


ApplicationModeConfiguration::ApplicationModeConfiguration(Workspace* workspace)
    : workspace_(workspace)
    , mainCanvas_("")
{
    tgtAssert(workspace_, "no workspace passed");
}

void ApplicationModeConfiguration::addPropertyGroup(const std::string& groupName) {
    tgtAssert(workspace_, "no workspace set");

    if (hasPropertyGroup(groupName)) {
        LWARNINGC("voreen.Workspace.ApplicationModeConfiguration", "Property group '" << groupName << "' already present");
        return;
    }

    propertyGroups_.push_back(groupName);

    workspace_->notifyApplicationModeConfigChanged();
}

void ApplicationModeConfiguration::removePropertyGroup(const std::string& groupName) {
    tgtAssert(workspace_, "no workspace set");

    std::vector<std::string>::iterator it = std::find(propertyGroups_.begin(), propertyGroups_.end(), groupName);
    if (it != propertyGroups_.end()) {
        propertyGroups_.erase(it);
    }
    tgtAssert(std::find(propertyGroups_.begin(), propertyGroups_.end(), groupName) == propertyGroups_.end(), "group not removed");

    // remove membership entries for this group
    for (std::map<const Property*, std::pair<std::string, int> >::iterator it=propertyGroupMembership_.begin(); it != propertyGroupMembership_.end(); ) {
        if (it->second.first == groupName)
            propertyGroupMembership_.erase(it++); //< advances the iterator before it becomes invalid by the deletion
        else
            it++;
    }

    workspace_->notifyApplicationModeConfigChanged();
}

void ApplicationModeConfiguration::removePropertyGroup(size_t i) {
    tgtAssert(workspace_, "no workspace set");

    if (i >= propertyGroups_.size())
        return;

    removePropertyGroup(std::string(propertyGroups_.at(i)));

    workspace_->notifyApplicationModeConfigChanged();
}

void ApplicationModeConfiguration::renamePropertyGroup(size_t index, const std::string& newGroupName) {
    tgtAssert(workspace_, "no workspace set");

    if (index >= propertyGroups_.size())
        return;

    // replace group name
    std::string prevName = propertyGroups_[index];
    propertyGroups_[index] = newGroupName;

    // update membership entries for this group
    for (std::map<const Property*, std::pair<std::string, int> >::iterator it=propertyGroupMembership_.begin(); it != propertyGroupMembership_.end(); it++) {
        if (it->second.first == prevName)
            it->second.first = newGroupName;
    }

    workspace_->notifyApplicationModeConfigChanged();
}

void ApplicationModeConfiguration::swapPropertyGroups(size_t indexA, size_t indexB) {
    tgtAssert(workspace_, "no workspace set");

    if (indexA == indexB || indexA >= propertyGroups_.size() || indexB >= propertyGroups_.size())
        return;

    std::swap(propertyGroups_.at(indexA), propertyGroups_.at(indexB));

    workspace_->notifyApplicationModeConfigChanged();
}

bool ApplicationModeConfiguration::hasPropertyGroup(const std::string& groupName) const {
    return std::find(propertyGroups_.begin(), propertyGroups_.end(), groupName) != propertyGroups_.end();
}

std::vector<std::string> ApplicationModeConfiguration::getPropertyGroups() const {
    return propertyGroups_;
}

void ApplicationModeConfiguration::setPropertyGroupMembership(const Property* property, const std::string& groupName, int priority) {
    tgtAssert(workspace_, "no workspace set");
    tgtAssert(property, "null pointer passed");

    if (!hasPropertyGroup(groupName)) {
        tgtAssert(false, "property group does not exist");
        return;
    }

    if (priority < 0) { // remove property from group
        std::map<const Property*, std::pair<std::string, int> >::iterator it = propertyGroupMembership_.find(property);
        if (it != propertyGroupMembership_.end() && it->second.first == groupName)
            propertyGroupMembership_.erase(it);
    }
    else { // insert/overwrite entry
        propertyGroupMembership_[property] = std::make_pair(groupName, priority);
    }

    workspace_->notifyApplicationModeConfigChanged();
}

std::pair<std::string, int> ApplicationModeConfiguration::getPropertyGroupMembership(const Property* property) const {
    std::map<const Property*, std::pair<std::string, int> >::const_iterator it = propertyGroupMembership_.find(property);
    if (it == propertyGroupMembership_.end()) { // property is not member of a group
        return std::pair<std::string, int>("", -1);
    }
    else {
        tgtAssert(hasPropertyGroup(it->second.first), "invalid property group membership entry: group ID does not exist");
        return it->second;
    }
}

int ApplicationModeConfiguration::getPropertyGroupPriority(const Property* property, const std::string& groupName) const {
    std::pair<std::string, int> membership = getPropertyGroupMembership(property);
    if (membership.first == groupName) //< property belongs to specified group
        return membership.second;
    else
        return -1;
}

std::vector<Property*> ApplicationModeConfiguration::getGroupProperties(const std::string& groupName) const {
    std::vector<Property*> result;
    if (!hasPropertyGroup(groupName))
        return result;

    // collect group properties in temp vector, ordered by priority
    std::vector<std::pair<Property*, int> > groupProperties;
    std::map<const Property*, std::pair<std::string, int> >::const_iterator it;
    for (it = propertyGroupMembership_.begin(); it != propertyGroupMembership_.end(); it++) {
        Property* property = const_cast<Property*>(it->first);
        std::string tGroupName = it->second.first;
        int groupPriority = it->second.second;

        if (tGroupName == groupName && groupPriority >= 0) { //< property does belong to group => insert into temp vector (sorted)
            size_t index = 0;
            while (index < groupProperties.size() && groupProperties.at(index).second < groupPriority)
                index++;
            groupProperties.insert(groupProperties.begin() + index, std::make_pair(property, groupPriority));
        }
    }

    // copy group properties from temp vector to result vector
    for (size_t i=0; i<groupProperties.size(); i++)
        result.push_back(groupProperties.at(i).first);

    return result;
}

bool ApplicationModeConfiguration::isPropertyVisible(const Property* property) const {
    tgtAssert(property, "null pointer passed");
    std::pair<std::string, int> groupMembership = getPropertyGroupMembership(property);
    return (groupMembership.second >= 0);
}

void ApplicationModeConfiguration::removeProcessorProperties(const Processor* processor) {
    tgtAssert(workspace_, "no workspace set");
    tgtAssert(processor, "null pointer passed");
    const std::vector<Property*>& procProperties = processor->getProperties();
    for (size_t i=0; i<procProperties.size(); i++) {
        propertyGroupMembership_.erase(procProperties.at(i));
    }

    workspace_->notifyApplicationModeConfigChanged();
}

void ApplicationModeConfiguration::clearGroups() {
    tgtAssert(workspace_, "no workspace set");

    propertyGroups_.clear();
    propertyGroupMembership_.clear();

    workspace_->notifyApplicationModeConfigChanged();
}

void ApplicationModeConfiguration::setMenuEntityVisibility(std::string& name, bool visible) {
    std::pair<std::map<std::string, bool>::iterator,bool> tmp = menuEntityVisiblity_.insert(std::pair<std::string,bool>(name,visible));
    //if second = false, name already existed => update visible state
    if(!tmp.second)
        tmp.first->second = visible;
    workspace_->notifyApplicationModeConfigChanged();
}

std::pair<bool, bool> ApplicationModeConfiguration::isMenuEntityVisible(const std::string& name) const{
    std::map<std::string, bool>::const_iterator it = menuEntityVisiblity_.find(name);
    if(it != menuEntityVisiblity_.end())
        return std::pair<bool, bool>(true, it->second);
    else
        return std::pair<bool, bool>(false,false);
}

std::string ApplicationModeConfiguration::getMainCanvas() const {
    return mainCanvas_;
}

void ApplicationModeConfiguration::setMainCanvas(std::string mainCanvas) {
    mainCanvas_ = mainCanvas;
    workspace_->notifyApplicationModeConfigChanged();
}

void ApplicationModeConfiguration::serialize(Serializer& s) const  {
    s.serialize("PropertyGroups", propertyGroups_);

    // create serialization map, which uses the property ID instead of the property itself as key
    std::map<std::string, std::pair<std::string, int> > tempMembershipMap;
    if (!propertyGroupMembership_.empty()) {

        // collect currently present properties from network, in order to filter out invalid property pointers (i.e., of deleted processors)
        std::set<const Property*> properties;
        std::map<std::string, Property*> idToProperty;
        if (workspace_ && workspace_->getProcessorNetwork()) {
            const std::vector<Processor*> processors = workspace_->getProcessorNetwork()->getProcessors();
            for (size_t i=0; i<processors.size(); i++) {
                Processor* proc = processors.at(i);
                for (size_t j=0; j<proc->getProperties().size(); j++)
                    properties.insert(proc->getProperties().at(j));
            }
        }

        // actually create serialization map
        for (std::map<const Property*, std::pair<std::string, int> >::const_iterator it=propertyGroupMembership_.begin(); it!=propertyGroupMembership_.end(); it++) {
            if (properties.find(it->first) != properties.end()) //< property is in network
                tempMembershipMap[it->first->getFullyQualifiedID()] = it->second;
            else
                LWARNINGC("voreen.ApplicationModeConfiguration", "property not in network");
        }
    }
    s.serialize("GroupMembership", tempMembershipMap);
    s.serialize("MenuEntityVisibility",menuEntityVisiblity_);
    s.serialize("MainCanvas", mainCanvas_);
}

void ApplicationModeConfiguration::deserialize(Deserializer& s) {
    clearGroups();

    std::map<std::string, bool> defaultMap;
    s.optionalDeserialize("MenuEntityVisibility", menuEntityVisiblity_, defaultMap);
    std::string defaultStr("");
    s.optionalDeserialize("MainCanvas", mainCanvas_, defaultStr);

    s.deserialize("PropertyGroups", propertyGroups_);

    // deserialize temporary group membership map and transform it to final map (i.e. replace propertyID -> property object)
    std::map<std::string, std::pair<std::string, int> > tempMembershipMap;
    s.deserialize("GroupMembership", tempMembershipMap);
    if (!tempMembershipMap.empty()) {

        // create auxiliary map from property ids to Property objects
        std::map<std::string, Property*> idToProperty;
        if (workspace_ && workspace_->getProcessorNetwork()) {
            const std::vector<Processor*> processors = workspace_->getProcessorNetwork()->getProcessors();
            for (size_t i=0; i<processors.size(); i++) {
                Processor* proc = processors.at(i);
                for (size_t j=0; j<proc->getProperties().size(); j++) {
                    Property* prop = proc->getProperties().at(j);
                    idToProperty[prop->getFullyQualifiedID()] = prop;
                }
            }
        }

        // create final map
        for (std::map<std::string, std::pair<std::string, int> >::iterator it=tempMembershipMap.begin(); it!=tempMembershipMap.end(); it++) {
            std::string propertyID = it->first;
            if (idToProperty.count(propertyID)) {
                const std::pair<std::string, int>& membershipEntry = it->second;
                if (std::find(propertyGroups_.begin(), propertyGroups_.end(), membershipEntry.first) != propertyGroups_.end())
                    propertyGroupMembership_[idToProperty[propertyID]] = membershipEntry;
                else
                    LWARNINGC("voreen.ApplicationModeConfiguration", "Property group '" << membershipEntry.first << "' does not exist");
            }
            else {
                LWARNINGC("voreen.ApplicationModeConfiguration", "Property with id '" << propertyID << "' not found in network");
            }
        }
    }

}

//-----------------------------------------------------------------------------

Workspace::Workspace()
    : version_(WORKSPACE_VERSION)
    , network_(new ProcessorNetwork())
    , animation_(0)
    , filename_("")
    , readOnly_(false)
    , description_("")
    , applicationModeConfig_(this)
    , modified_(false)
{
    network_->addObserver(static_cast<ProcessorNetworkObserver*>(this));
}

Workspace::~Workspace() {
    clearResources();
}

std::vector<std::string> Workspace::getErrors() const {
    return errorList_;
}

void Workspace::load(const std::string& filename, const std::string& workDir) {
    // open file for reading
    std::fstream fileStream(filename.c_str(), std::ios_base::in);
    if (fileStream.fail()) {
        //LERROR("Failed to open file '" << tgt::FileSystem::absolutePath(filename) << "' for reading.");
        throw SerializationException("Failed to open workspace file '" + tgt::FileSystem::absolutePath(filename) + "' for reading.");
    }

    std::string documentPath;
    if (!workDir.empty())
        documentPath = workDir;
    else
        documentPath = filename;

    // read data stream into deserializer
    XmlDeserializer d(documentPath);
    d.setUseAttributes(true);
    NetworkSerializer ser;
    try {
        d.read(fileStream, &ser);
    }
    catch (SerializationException& e) {
        throw SerializationException("Failed to read serialization data stream from workspace file '"
                                     + filename + "': " + e.what());
    }
    catch (...) {
        throw SerializationException("Failed to read serialization data stream from workspace file '"
                                     + filename + "' (unknown exception).");
    }

    // deserialize workspace from data stream
    try {
        d.deserialize("Workspace", *this);
        errorList_ = d.getErrors();
        setFilename(filename);
    }
    catch (std::exception& e) {
        throw SerializationException("Deserialization from workspace file '" + filename + "' failed: " + e.what());
    }
    catch (...) {
        throw SerializationException("Deserialization from workspace file '" + filename + "' failed (unknown exception).");
    }
}

void Workspace::save(const std::string& filename, bool overwrite, const std::string& workDir) {
    // check if file is already present
    if (!overwrite && tgt::FileSystem::fileExists(filename))
        throw SerializationException("File '" + filename + "' already exists.");

    std::string documentPath;
    if (!workDir.empty())
        documentPath = workDir;
    else
        documentPath = filename;

    // serialize workspace
    XmlSerializer s(documentPath);
    s.setUseAttributes(true);
    s.serialize("Workspace", *this);
    errorList_ = s.getErrors();

    // write serialization data to temporary string stream
    std::ostringstream textStream;

    try {
        s.write(textStream);
        if (textStream.fail())
            throw SerializationException("Failed to write serialization data to string stream.");
    }
    catch (std::exception& e) {
        throw SerializationException("Failed to write serialization data to string stream: " + std::string(e.what()));
    }
    catch (...) {
        throw SerializationException("Failed to write serialization data to string stream (unknown exception).");
    }

    // Now we have a valid StringStream containing the serialization data.
    // => Open output file and write it to the file.
    // For added data security we write to a temporary file and afterwards move it into place
    // (which should be an atomic operation).
    const std::string tmpfilename = filename + ".tmp";
    std::fstream fileStream(tmpfilename.c_str(), std::ios_base::out);
    if (fileStream.fail())
        throw SerializationException("Failed to open file '" + tmpfilename + "' for writing.");

    try {
        fileStream << textStream.str();
    }
    catch (std::exception& e) {
        throw SerializationException("Failed to write serialization data stream to file '"
                                     + tmpfilename + "': " + std::string(e.what()));
    }
    catch (...) {
        throw SerializationException("Failed to write serialization data stream to file '"
                                     + tmpfilename + "' (unknown exception).");
    }
    fileStream.close();

    // Finally move the temporary file into place. It is important that this happens in-place,
    // without deleting the old file first.
    bool success;
#ifdef WIN32
    // rename() does not replace existing files on Windows, so we have to use this
    success = (MoveFileEx(tmpfilename.c_str(), filename.c_str(),
                          MOVEFILE_REPLACE_EXISTING | MOVEFILE_COPY_ALLOWED) != 0);

#else
    // atomic replace
    success = (rename(tmpfilename.c_str(), filename.c_str()) == 0);
#endif
    if (!success) {
#ifdef WIN32
        _unlink(tmpfilename.c_str()); // ignore failure here
#else
        unlink(tmpfilename.c_str()); // ignore failure here
#endif
        throw SerializationException("Failed to rename temporary file '" + tmpfilename + "' to '"
                                     + filename + "'");
    }

    // saving successful
    setFilename(filename);

    //a saved workspace in no longer modified
    setModified(false);
}

void Workspace::clear() {
    clearResources();
    setModified(true);
}

void Workspace::clearResources() {
    filename_ = "";

    // animation
    if (animation_) {
        delete animation_;
        animation_ = 0;
    }

    // application mode settings
    applicationModeConfig_.clearGroups();

    // network
    if (network_) {
        delete network_;
        network_ = 0;
    }

    //clear description
    description_ = "";

    errorList_.clear();
    readOnly_ = false;
}

bool Workspace::readOnly() const {
    return readOnly_;
}

ProcessorNetwork* Workspace::getProcessorNetwork() const {
    return network_;
}

void Workspace::setFilename(const std::string& filename) {
    filename_ = filename;

    // replace backslashes
    std::string::size_type pos = filename_.find("\\");
    while (pos != std::string::npos) {
        filename_[pos] = '/';
        pos = filename_.find("\\");
    }

}

std::string Workspace::getFilename() const {
    return filename_;
}

bool Workspace::isModified() const {
    return modified_;
}

void Workspace::setModified(bool isModified) {
    bool oldValue = modified_;
    if (modified_ != isModified)
        modified_ = isModified;
    if(isModified || oldValue) //notify, if it has been modified, or has changed to no modification
        notifyWorkspaceHasBeenModified();
}

void Workspace::setProcessorNetwork(ProcessorNetwork* network) {
    network_ = network;
    if (network_) {
        network_->setWorkspace(this);
        if(!network_->isObservedBy(static_cast<ProcessorNetworkObserver*>(this)))
            network_->addObserver(static_cast<ProcessorNetworkObserver*>(this));
    }
}

void Workspace::serialize(Serializer& s) const {
    // Always serialize newest version.
    s.serialize("version", WORKSPACE_VERSION);

    //HACK: new workspaces are always NOT read only
    s.serialize("readonly", false);

    // Serialize network...
    s.serialize("ProcessorNetwork", network_);

    // Serialize animation...
    s.serialize("Animation", animation_);

    // serialize application mode description
    s.serialize("ApplicationModeConfig", applicationModeConfig_);

    s.serialize("GlobalDescription", description_);
}

void Workspace::deserialize(Deserializer& s) {
    // clear existing network and containers
    clearResources();

    // Deserialize version number for incompatibility checks.
    s.deserialize("version", version_);

    try {
        s.deserialize("readonly", readOnly_);
    }
    catch (SerializationNoSuchDataException&) {
        s.removeLastError();
        readOnly_ = false;
    }

    // Deserialize network...
    s.deserialize("ProcessorNetwork", network_);

    // Deserialize animation if present...
    try {
        s.deserialize("Animation", animation_);
    } catch (SerializationNoSuchDataException&) {
        s.removeLastError();
    }

    if (animation_) {
        animation_->setNetwork(network_);

        // register this as observer on all propertytimelines for undo/redo
        const std::vector<AnimatedProcessor*> animproc = this->animation_->getAnimatedProcessors();
        std::vector<AnimatedProcessor*>::const_iterator it;
        for (it = animproc.begin();it != animproc.end();it++)
        {
            const std::vector<PropertyTimeline*> timelines = (*it)->getPropertyTimelines();
            std::vector<PropertyTimeline*>::const_iterator it2;
            for (it2 = timelines.begin();it2 != timelines.end(); ++it2) {
                (*it2)->registerUndoObserver(this->animation_);
            }
        }

        // register this as observer in the processornetwork to register added and removed processors
        ProcessorNetwork* net = const_cast<ProcessorNetwork*>(network_);
        net->addObserver(this->animation_);
    }

    // deserialize application mode config
    try {
        s.deserialize("ApplicationModeConfig", applicationModeConfig_);
    }
    catch (SerializationNoSuchDataException&) {
        s.removeLastError();
        applicationModeConfig_.clearGroups();
    }

    try {
        s.deserialize("GlobalDescription", description_);
    }
    catch (SerializationNoSuchDataException&) {
        s.removeLastError();
        description_ = "";
    }

    if (network_) {
        network_->setWorkspace(this);
        if(!network_->isObservedBy(static_cast<ProcessorNetworkObserver*>(this)))
            network_->addObserver(static_cast<ProcessorNetworkObserver*>(this));
    }
}

Animation* Workspace::getAnimation() const {
    return animation_;
}

void Workspace::setAnimation(Animation* anim) {
    if(animation_ != anim) {
        animation_ = anim;
        /*setModified(true); //FIXME: new animation should trigger modified
        notifyWorkspaceHasBeenModified();*/
    }
}

const std::string& Workspace::getDescription() const {
    return description_;
}

bool Workspace::hasDescription() const {
    return !description_.empty();
}

void Workspace::setDescription(const std::string& description) {
    if(description.compare(description_) != 0) {
        description_ = description;
        setModified(true);
    }
}

const ApplicationModeConfiguration& Workspace::getApplicationModeConfig() const {
    return applicationModeConfig_;
}
ApplicationModeConfiguration& Workspace::getApplicationModeConfig() {
    return applicationModeConfig_;
}

void Workspace::notifyApplicationModeConfigChanged() const {
    std::vector<WorkspaceObserver*> observers = getObservers();
    for (size_t i=0; i<observers.size(); i++)
        observers.at(i)->applicationModeConfigurationChanged();
}

void Workspace::notifyWorkspaceHasBeenModified() const {
    std::vector<WorkspaceObserver*> observers = getObservers();
    for (size_t i=0; i<observers.size(); i++)
        observers.at(i)->workspaceHasBeenModified();
}

} // namespace
