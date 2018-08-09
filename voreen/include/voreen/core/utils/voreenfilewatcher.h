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

#ifndef VRN_VOREENFILEWATCHER_H
#define VRN_VOREENFILEWATCHER_H

#include "voreen/core/utils/exception.h"
#include "voreen/core/properties/property.h"
#include "voreen/core/voreencoreapi.h"

#include "tgt/singleton.h"

#include "efsw/include/efsw/efsw.hpp"

#include <string>
#include <set>
#include <map>

namespace voreen {

class VoreenFileWatchListener;

/**
 * This class watches the file system. It manages listeners and informs them
 * if any change was made.
 */
class VRN_CORE_API VoreenFileWatcher : public tgt::Singleton<VoreenFileWatcher> {
    friend class VoreenFileWatchListener;
public:

    VoreenFileWatcher();
    ~VoreenFileWatcher();

    /**
     * Registers a directory watch for the specified listener.
     * @param listener listener to register a directory for
     * @param directory directory to be watched by listener
     */
    void requestAddWatch(VoreenFileWatchListener* listener, const std::string& directory);

    /**
     * Degeristers a directory watch for the specified listener.
     * @param listener listener to deregister a directory for
     * @param directory directory to be no longer watched by listener
     */
    void requestRemoveWatch(VoreenFileWatchListener* listener, const std::string& directory);

    /**
     * Deregisters all file watches registered for the specified listener.
     * @param listener listener to remove all watched for
     */
    void requestRemoveAllWatches(VoreenFileWatchListener* listener);

private:

    // Delete copy constructor and assignment operator, which enables use of unique_ptr for members.
    VoreenFileWatcher(const VoreenFileWatcher&);
    VoreenFileWatcher& operator=(const VoreenFileWatcher&);

    /**
     * The MetaFileWatchListener acts as the connection between the filewatch (being implemented as directory watch)
     * and the actual VoreenFileWatchListeners. Since multiple listeners can watch a file or directory inside the same
     * directory but only one listener at a time can be registered per directory, this class forwards the callback
     * to the actual receiver.
     */
    class MetaFileWatchListener : public efsw::FileWatchListener {
    public:
        MetaFileWatchListener(efsw::FileWatcher* parent, const std::string& directory);
        ~MetaFileWatchListener();

        /**
         * This function is called by the async file watch thread.
         * For each registered listener the implemented fileActionCallback
         * will be enqueued into the CommandQueue.
         *
         * Note, that all old commands added by the registered listeners will be removed beforehand.
         */
        virtual void handleFileAction(efsw::WatchID watchid, const std::string& dir, const std::string& filename, efsw::Action action, std::string oldFilename = "");

        std::set<VoreenFileWatchListener*> listeners_; ///< all registered VoreenFileWatchListener for this watch

    private:

        efsw::WatchID watchID_; ///< unique watch id
        efsw::FileWatcher* parent_;
    };

    /// Maps directories to their listeners. 
    std::map<std::string, std::unique_ptr<MetaFileWatchListener> > metaListeners_; 
    efsw::FileWatcher fileWatcher_;
};

/**
 * This class should be used within properties which do rely on files such as:
 * FileDialogProperty, VolumeURLProperty, etc.
 * The owning property gets invalidated as soon as the underlying file was either:
 * moved, renamed, deleted, readded
 */
class VRN_CORE_API VoreenFileWatchListener /*: public Serializable*/ {
    friend class VoreenFileWatcher::MetaFileWatchListener;
public:

    /**
    * Modes a property can be initialized by to determine
    * if the file watcher should handle file changes
    * for any particular property.
    *
    * This is used by:
    *   - FileDialogProperty
    *   - VolumeURLProperty
    *   - VolumeURLListProperty
    */
    enum WatchMode {
        OPTIONAL_ON,  // selectable, on by default
        OPTIONAL_OFF, // selectable, off by default
        ALWAYS_ON,    // not selectable, on
        ALWAYS_OFF    // not selectable, off
    };

    /**
     * Constructor.
     *
     * @param watchMode 
     */
    VoreenFileWatchListener(WatchMode watchMode = OPTIONAL_OFF);

    /**
     * Destructor.
     */
    ~VoreenFileWatchListener();

    /**
     * Adds a path to be watched for changes.
     * Empty paths will succeed but be ignored.
     * Path may be watched already.
     *
     * Note that the file or directory itself does not need to exist
     * when function is called, but the containing directory does.
     * This is necessary since the path needs to be resolved
     * to it's absolute variant.
     *
     * Returns true if successful, false otherwise.
     */
    bool addWatch(const std::string& path);

    /**
     * Removes a path.
     * Empty paths will succeed but be ignored.
     * Path must be watched already.
     *
     * Returns true if successful, false otherwise.
     */
    bool removeWatch(const std::string& path);

    /**
     * Removes all watches.
     */
    void removeAllWatches();

    /**
     * Determines, if the specified path is currently being watched.
     */
    bool isWatching(const std::string& path) const;

    /**
     * Returns the currently watched paths.
     *
     * Note that even if a path is being watched multiple times
     * it will only be contained once in the returned vector.
     */
    std::vector<std::string> getWatches() const;

    /**
     * Returns the watch mode the listener was initialized with.
     */
    WatchMode getWatchMode() const;

    /**
     * Determines if file watch can be disabled/enabled.
     * This is true, if mode set to OPTION_ON or OPTIONAL_OFF.
     */
    bool isFileWatchEditable() const;

    /**
     * Sets the file watch policy.
     * If watch mode is ALWAYS_ON or ALWAYS_OFF, this has no effect.
     */
    void setFileWatchEnabled(bool enabled);
    bool isFileWatchEnabled() const;

    /// @see Serializable::serialize
    virtual void serialize(Serializer& s) const;

    /// @see Serializable::deserialize
    virtual void deserialize(Deserializer& s);

    /**
     * This function will be enqueued into the CommandQueue by a MetaFileWatchListener.
     * Reimplement in subclass to specify behavior.
     *
     * Note that any old command originating from this callback being
     * enqueued but not yet executed will be removed beforehand.
     */
    virtual void fileActionCallback() = 0;

private:

    std::string getParentDir(std::string path) const;

    std::map<std::string, std::string> absolutePaths_;      ///< Maps incoming path to it's absolute parent path
    std::multimap<std::string, std::string> watchedPaths_;  ///< Maps parent directories to contained watches

    bool fileWatchEnabled_;
    WatchMode watchMode_; ///< Watch mode
};

} // namespace

#endif // VRN_VOREENFILEWATCHER_H
