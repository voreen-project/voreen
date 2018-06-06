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

#ifndef VRN_VOLUMELIST_H
#define VRN_VOLUMELIST_H

#include "voreen/core/utils/observer.h"
#include "voreen/core/datastructures/volume/volume.h"

#include <vector>
#include <string>

namespace voreen {

class VolumeList;
class Modality;

/**
 * Interface for volume list observers.
 */
class VRN_CORE_API VolumeListObserver : public Observer {
public:

    /**
     * This method is called by the observed list after
     * a volume has been added.
     *
     * @param source the calling volume list
     * @param handle the volume that has been added
     */
    virtual void volumeAdded(const VolumeList* /*source*/, const VolumeBase* /*handle*/) {};

    /**
     * This method is called by the observed list after
     * a volume has been removed.
     *
     * @param source the calling volume list
     * @param handle the volume that has been removed
     */
    virtual void volumeRemoved(const VolumeList* /*source*/, const VolumeBase* /*handle*/) {};

    /**
    * This method is called by the observed list after
    * a volume has been changed. This usually happens when
    * a volume reloads the wrapped volume.
    *
    * @param source the calling list
    * @param handle the volume that has been changed
    */
    virtual void volumeChanged(const VolumeList* /*source*/, const VolumeBase* /*handle*/) {};

};

#ifdef DLL_TEMPLATE_INST
template class VRN_CORE_API Observable<VolumeListObserver>;
#endif

/**
 * Ordered list of volumes, which can be observed.
 *
 * The volume list does not perform memory management, i.e., added volumes
 * are not deleted on removal or destruction of the list.
 *
 * @see VolumeContainer
 */
class VRN_CORE_API VolumeList : public Observable<VolumeListObserver>,
                                public VolumeObserver,
                                public DataInvalidationObservable,
                                public DataInvalidationObserver {

public:

    VolumeList();

    virtual ~VolumeList();

    /**
     * Adds the passed volume to the list,
     * if it is not already contained.
     *
     * The VolumeList does not take ownership of the
     * added Volume and does therefore not delete it
     * on its own destruction.
     */
    virtual void add(VolumeBase* volumeHandle);

    /**
     * Adds volume contained by the passed VolumeList
     * to the list, if they are not already contained.
     *
     * The VolumeList does not take ownership of the
     * added volume and does therefore not delete them
     * on its own destruction.
     */
    virtual void add(const VolumeList* volumeList);

    /**
     * Removes the passed volume from the VolumeList
     * without deleting it.
     */
    virtual void remove(const VolumeBase* volume);

    /**
     * Removes all volume contained by the passed list
     * without deleting them.
     */
    virtual void remove(const VolumeList* volumeList);

    /**
     * Returns whether the passed volume is contained
     * by the list.
     */
    virtual bool contains(const VolumeBase* volume) const;

    /**
     * Returns the VolumeElement at a specified index position.
     *
     * @param The index of the volume to return. Must be valid, i.e. i < size().
     */
    virtual VolumeBase* at(size_t i) const;

    /**
     * Returns the first volume of the list, or null
     * if the list is empty.
     */
    virtual VolumeBase* first() const;

    /**
     * Clears the list without deleting the volumes.
     */
    virtual void clear();

    /**
     * Returns a list containing all volumes of the
     * specified modality.
     */
    virtual VolumeList* selectModality(const Modality& name) const;

    /**
     * Returns a list containing all volumes with the
     * specified origin.
     */
    virtual VolumeList* selectOrigin(const VolumeURL& origin) const;

    /**
     * Returns a new list containing the volumes between the passed start and end indices (inclusive).
     */
    virtual VolumeList* subList(size_t start, size_t end) const;

    /**
     * Returns a new list containing the volumes at the passed indices.
     */
    virtual VolumeList* subList(const std::vector<size_t>& indices) const;

    /**
     * Returns the number of volumes contained by the list.
     */
    virtual size_t size() const;

    /**
     * Returns whether the list is currently empty.
     */
    virtual bool empty() const;

    /**
     * Returns the hash of the list based on the hashes of the contained volumes.
     */
    virtual std::string getHash() const;

    /// @see VolumeObserver::volumeChange
    virtual void volumeChange(const VolumeBase* volume);

    /// @see VolumeObserver::volumeDelete
    virtual void volumeDelete(const VolumeBase* volume);

    /// @see DataInvalidationObserver::dataAboutToInvalidate
    virtual void dataAboutToInvalidate(const DataInvalidationObservable* data);

protected:
    /**
     * Returns an iterator pointing to the position of the passed Volume
     * within the volumeHandles_ vector. Returns volumeHandles_.end(), if
     * the handle is not contained by the list.
     */
    std::vector<VolumeBase*>::iterator find(const VolumeBase* volume);
    std::vector<VolumeBase*>::const_iterator find(const VolumeBase* volume) const;

    /// Notifies all VolumeListObservers that a volume has been added.
    void notifyVolumeAdded(const VolumeBase* volume);
    /// Notifies all VolumeListObservers that a volumes has been removed.
    void notifyVolumeRemoved(const VolumeBase* volume);
    /// Notifies all VolumeListObservers that a volume has been changed.
    void notifyVolumeChanged(const VolumeBase* volume);

    /// Vector storing the volumes contained by the list.
    std::vector<VolumeBase*> volumes_;

    /// category for logging.
    static const std::string loggerCat_;
};

}   // namespace

#endif
