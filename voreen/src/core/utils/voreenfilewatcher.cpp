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

#include "voreen/core/utils/voreenfilewatcher.h"

#include "voreen/core/voreenapplication.h"
#include "voreen/core/utils/commandqueue.h"

#include "tgt/assert.h"
#include "tgt/filesystem.h"

namespace voreen {

//=========================================================================================
// VoreenFileWatcher

VoreenFileWatcher::VoreenFileWatcher()
#ifdef VRN_USE_GENERIC_FILE_WATCHER
    : fileWatcher_(true) // Using the generic file watcher prevents racing conditions under windows and allows remote file watching.
#endif
{
    fileWatcher_.followSymlinks(true);
    //fileWatcher_.allowOutOfScopeLinks(true); // Could lead to crashes due to endless recursion if symlinks build recursive trees.
    fileWatcher_.watch();
}

VoreenFileWatcher::~VoreenFileWatcher() {
    tgtAssert(metaListeners_.empty(), "Some file watches still registered");
}

void VoreenFileWatcher::requestAddWatch(VoreenFileWatchListener* listener, const std::string& directory) {
    tgtAssert(listener, "listener must not be null");

    // Create Meta Listener if not already existing.
    if (!metaListeners_[directory].get())
        metaListeners_[directory].reset(new MetaFileWatchListener(&fileWatcher_, directory));

    // Register listener.
    metaListeners_[directory]->listeners_.insert(listener);
}

void VoreenFileWatcher::requestRemoveWatch(VoreenFileWatchListener* listener, const std::string& directory) {
    tgtAssert(listener, "listener must not be null");
    tgtAssert(metaListeners_.find(directory) != metaListeners_.end(), "no watch registered for this directory");
    tgtAssert(metaListeners_[directory]->listeners_.find(listener) != metaListeners_[directory]->listeners_.end(), "listener not registered for the specified directory");

    // Deregister listener.
    metaListeners_[directory]->listeners_.erase(listener);
    
    // Remove file watch if no listener is left.
    if (metaListeners_[directory]->listeners_.empty())
        metaListeners_.erase(directory);
}

void VoreenFileWatcher::requestRemoveAllWatches(VoreenFileWatchListener* listener) {
    tgtAssert(listener, "listener must not be null");

    // Deregister listener.
    for (auto iter = metaListeners_.begin(); iter != metaListeners_.end();) {
        std::set<VoreenFileWatchListener*>& listeners = iter->second->listeners_;
        auto it = listeners.find(listener);
        if (it != listeners.end()) {
            listeners.erase(it);

            // Remove file watch if no listener is left.
            if (listeners.empty()) {
                iter = metaListeners_.erase(iter);
                continue;
            }
        }

        iter++;
    }
}

//=========================================================================================
// VoreenFileWatcher::FileWatchListener

VoreenFileWatcher::MetaFileWatchListener::MetaFileWatchListener(efsw::FileWatcher* parent, const std::string& directory)
    : parent_(parent)
{
    // Request file watch ID.
    watchID_ = parent_->addWatch(directory, this, false);

    // Check for errors.
    switch (watchID_) {
    case efsw::Error::FileNotFound:
        LWARNINGC("VoreenFileWatcher", "A FileNotFound error occurd watching: " + directory);
        break;
    case efsw::Error::FileNotReadable:
        LWARNINGC("VoreenFileWatcher", "A FileNotReadable error occurd watching: " + directory);
        break;
    case efsw::Error::FileOutOfScope:
        LWARNINGC("VoreenFileWatcher", "A FileOutOfScope error occurd watching: " + directory);
        break;
    case efsw::Error::FileRemote:
        LWARNINGC("VoreenFileWatcher", "A FileRemote error occurd watching: " + directory);
        break;
    case efsw::Error::FileRepeated:
        LWARNINGC("VoreenFileWatcher", "A FileRepeated error occurd watching: " + directory);
        break;
    case efsw::Error::Unspecified:
        LWARNINGC("VoreenFileWatcher", "An unspecified error occurd watching: " + directory);
        break;
    default:
        // Everything went fine.
        break;
    }
}

VoreenFileWatcher::MetaFileWatchListener::~MetaFileWatchListener() {
    if (watchID_ >= 0)
        parent_->removeWatch(watchID_);
}

void VoreenFileWatcher::MetaFileWatchListener::handleFileAction(efsw::WatchID watchid, const std::string& dir, const std::string& filename, efsw::Action action, std::string oldFilename) {

#ifdef VRN_DEBUG
    switch (action)
    {
    case efsw::Actions::Add:
        LDEBUGC("MetaFileWatchListener", "DIR (" << dir << ") FILE (" << filename << ") has event Added");
        break;
    case efsw::Actions::Delete:
        LDEBUGC("MetaFileWatchListener", "DIR (" << dir << ") FILE (" << filename << ") has event Delete");
        break;
    case efsw::Actions::Modified:
        LDEBUGC("MetaFileWatchListener", "DIR (" << dir << ") FILE (" << filename << ") has event Modified");
        break;
    case efsw::Actions::Moved:
        LDEBUGC("MetaFileWatchListener", "DIR (" << dir << ") FILE (" << filename << ") has event Moved from (" << oldFilename << ")");
        break;
    default:
        tgtAssert(false, "Unhandled File Action!");
        LERRORC("MetaFileWatchListener", "Unhandled File Action!");
        return;
    }
#endif

    // Enqueue file action.
    if (VoreenApplication::app()) {
        for (VoreenFileWatchListener* listener : listeners_) {
            std::string path = tgt::FileSystem::cleanupPath(dir + "/" + filename);
            if (listener->isWatching(path))
            {
                tgtAssert(listener->isFileWatchEnabled(), "File Watch Callback called unintentionally");

                VoreenApplication::app()->getCommandQueue()->enqueue(listener, MemberFunctionCallback<VoreenFileWatchListener>(listener, &VoreenFileWatchListener::fileActionCallback));
            }
        }
    }
}


//=========================================================================================
// VoreenFileWatchListener

VoreenFileWatchListener::VoreenFileWatchListener(WatchMode watchMode)
    : watchMode_(watchMode)
{
    fileWatchEnabled_ = (watchMode_ == WatchMode::OPTIONAL_ON || watchMode_ == ALWAYS_ON);
}

VoreenFileWatchListener::~VoreenFileWatchListener() {
    removeAllWatches();

    // Remove all Enqueued tasks.
    if (VoreenApplication::app())
        VoreenApplication::app()->getCommandQueue()->removeAll(this);
}

bool VoreenFileWatchListener::addWatch(const std::string& path) {

    // Ignore empty paths.
    if (path.empty())
        return true;

    // Retrieve parent dir.
    std::string parentDir = tgt::FileSystem::parentDir(path);
    if (parentDir.empty() && VoreenApplication::app())
        parentDir = VoreenApplication::app()->getBasePath();

    if (!tgt::FileSystem::dirExists(parentDir))
        return false;

    // Get actual file or folder name (also works for directories).
    std::string name = tgt::FileSystem::fileName(path);

    // Convert to absolute path.
    parentDir = tgt::FileSystem::cleanupPath(tgt::FileSystem::absolutePath(parentDir));
    std::string absolutePath = tgt::FileSystem::cleanupPath(parentDir + "/" + name);
    if (absolutePath != path) {
        // Overwrite old (identical) path, if set already.
        absolutePaths_[path] = absolutePath;
    }

    // Insert path (may occur multiple times).
    watchedPaths_.insert(std::make_pair(parentDir, name));

    // If parent directory is already being watched, skip request.
    if (watchedPaths_.count(parentDir) == 1 && isFileWatchEnabled() && VoreenFileWatcher::isInited())
        VoreenFileWatcher::getRef().requestAddWatch(this, parentDir);

    return true;
}

bool VoreenFileWatchListener::removeWatch(const std::string& path) {

    // Ignore empty paths.
    if (path.empty())
        return true;

    // Retrieve parent dir.
    std::string parentDir = getParentDir(path);

    // Get actual file or folder name (also works for directories).
    std::string name = tgt::FileSystem::fileName(path);

    // Check, if path is being watched.
    auto range = watchedPaths_.equal_range(parentDir);
    bool found = false, single = false;
    for (auto iter = range.first; iter != range.second; iter++) {
        if (iter->second == name) {
            if (!found) {
                found = true;
                single = true;
            }
            else {
                single = false;
                break;
            }
        }
    }

    if (!found)
        return false;

    if (single)
        absolutePaths_.erase(path);

    // If parent directory is still being watched due to another contained file, skip request.
    if (watchedPaths_.count(parentDir) == 0 && isFileWatchEnabled() && VoreenFileWatcher::isInited())
        VoreenFileWatcher::getRef().requestRemoveWatch(this, parentDir);

    return true;
}

void VoreenFileWatchListener::removeAllWatches() {
    if (VoreenFileWatcher::isInited())
        VoreenFileWatcher::getRef().requestRemoveAllWatches(this);

    absolutePaths_.clear();
    watchedPaths_.clear();
}

bool VoreenFileWatchListener::isWatching(const std::string& path) const {

    // Fast check for parent directory.
    std::string parentDir = getParentDir(path);
    if (watchedPaths_.find(parentDir) == watchedPaths_.end())
        return false;

    // Get actual file or folder name (also works for directories).
    std::string name = tgt::FileSystem::fileName(path);

    auto range = watchedPaths_.equal_range(parentDir);
    for (auto iter = range.first; iter != range.second; iter++) {
        if (iter->second == name) {
            return true;
        }
    }

    return false;
}

std::string VoreenFileWatchListener::getParentDir(std::string path) const {
    auto iter = absolutePaths_.find(path);
    if (iter != absolutePaths_.end())
        path = iter->second;
    std::string parentDir = tgt::FileSystem::parentDir(path);
    if (parentDir.empty() && VoreenApplication::app())
        parentDir = VoreenApplication::app()->getBasePath();
    return parentDir;
}

std::vector<std::string> VoreenFileWatchListener::getWatches() const {
    std::set<std::string> result;
    for (const auto& path : watchedPaths_) {
        result.insert(tgt::FileSystem::cleanupPath(path.first + "/" + path.second));
    }
    return std::vector<std::string>(result.begin(), result.end());
}

VoreenFileWatchListener::WatchMode VoreenFileWatchListener::getWatchMode() const {
    return watchMode_;
}

bool VoreenFileWatchListener::isFileWatchEditable() const {
    return watchMode_ == VoreenFileWatchListener::OPTIONAL_OFF || watchMode_ == VoreenFileWatchListener::OPTIONAL_ON;
}

void VoreenFileWatchListener::setFileWatchEnabled(bool enabled) {
    if (!isFileWatchEditable() || enabled == fileWatchEnabled_)
        return;

    // Toggle now since addWatch/removeWatch rely on this state.
    fileWatchEnabled_ = enabled;

    // When toggling file watching either remove or (re-)add listeners.
    if (VoreenFileWatcher::isInited()) {
        if (enabled) {
            for (const auto& path : watchedPaths_) {
                std::string parentDir = path.first;
                VoreenFileWatcher::getRef().requestAddWatch(this, parentDir);
            }
        }
        else
            VoreenFileWatcher::getRef().requestRemoveAllWatches(this);
    }
}

bool VoreenFileWatchListener::isFileWatchEnabled() const {
    return fileWatchEnabled_;
}

void VoreenFileWatchListener::serialize(Serializer& s) const {
    s.serialize("watchMode", watchMode_);
    s.serialize("fileWatchEnabled", fileWatchEnabled_);
}

void VoreenFileWatchListener::deserialize(Deserializer& s) {
    int watchMode = watchMode_;
    s.optionalDeserialize("watchMode", watchMode, watchMode);
    watchMode_ = static_cast<WatchMode>(watchMode);
    s.optionalDeserialize("fileWatchEnabled", fileWatchEnabled_, fileWatchEnabled_);
}

} // namespace
