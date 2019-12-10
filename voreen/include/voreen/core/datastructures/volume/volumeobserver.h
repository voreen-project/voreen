/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2019 University of Muenster, Germany,                        *
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

#ifndef VRN_VOLUMEOBSERVER_H
#define VRN_VOLUMEOBSERVER_H


#include "voreen/core/utils/observer.h"

namespace voreen {

class VolumeBase;

/**
 * Interface for volume observers.
 */
class VRN_CORE_API VolumeObserver : public Observer {
public:
    /**
     * This method is called by the observed Volume's destructor.
     *
     * @param source the calling Volume
     */
    virtual void volumeDelete(const VolumeBase* source) = 0;

    /**
     * This method is called by the observed Volume
     * after its member Volume object has changed.
     *
     * When this function is called, the new Volume object has
     * already been assigned. The former Volume object is still
     * valid at this point, but it is deleted immediately after
     * this function has been called.
     *
     * @param source the calling Volume
     */
    virtual void volumeChange(const VolumeBase* source) = 0;

    /**
     * This method is called by the observed Volume after a VolumeDerivedDataThread has finished.
     * Then the observer itself needs to check if the derived data has changed.
     */
    virtual void derivedDataThreadFinished(const VolumeBase* source) {}

    /**
     * This method is called when removing a representation (e.g., VolumeRAM) from the volume.
     *
     * Is not called by the destructor, where volumeDelete is called instead before deleting the representations.
     */
    virtual void volumeRepresentationDelete(const VolumeBase* source, const VolumeRepresentation* rep) { }

    /**
     * Special case of volumeDelete, which is only called by VolumeDecorator as a workaround for missing virtual template methods.
     * Is called when the decorated volume, i.e. the actual data is deleted, but the VolumeDecorator still exists.
     *
     * The default implementation will call volumeDelete, to make sure that data access is stopped (overwrite for other special cases).
     */
    virtual void volumeDataDelete(const VolumeBase* source) { volumeDelete(source); }
};

} // namespace

#endif
