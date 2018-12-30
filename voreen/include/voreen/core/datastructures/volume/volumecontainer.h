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

#ifndef VRN_VOLUMECONTAINER_H
#define VRN_VOLUMECONTAINER_H

#include "voreen/core/datastructures/volume/volumelist.h"
#include "voreen/core/datastructures/volume/volume.h"

#include <vector>
#include <string>

namespace voreen {

/**
 * Extension of the VolumeList that performs memory management,
 * i.e., the volume container takes ownership of passed volumes and
 * deletes them on removal or its own destruction.
 *
 * Additionally, this class provides convenience functions for
 * volume loading.
 *
 * @see VolumeList
 */
class VRN_CORE_API VolumeContainer : public VolumeList {
public:

    VolumeContainer();

    /**
     * Deletes all contained VolumeHandles.
     *
     * @see releaseAll
     */
    virtual ~VolumeContainer();

    /**
     * Deletes the passed Volume and removes it
     * from the container. If the passed handle is not contained
     * by the container, this function has no effect.
     *
     * @see releaseVolume
     */
    virtual void remove(const VolumeBase* handle);

    /**
     * Deletes all VolumeHandles contained by the passed collection
     * and removes them from the container. Handles that are not
     * contained by the container are not touched.
     *
     * @see releaseVolume
     */
    virtual void remove(const VolumeList* volumeList);

    /**
     * Removes all VolumeHandles objects from the container and deletes them.
     *
     * @see releaseAll
     */
    virtual void clear();

    /**
     * Removes the passed volume from the container
     * without deleting it.
     */
    virtual void release(const VolumeBase* handle);

    /**
     * Removes all contained volumes from the container
     * without deleting them.
     */
    virtual void releaseAll();

protected:

    /// category for logging.
    static const std::string loggerCat_;
};

}   // namespace

#endif
