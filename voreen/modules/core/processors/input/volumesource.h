/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2020 University of Muenster, Germany,                        *
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

#ifndef VRN_VOLUMESOURCEPROCESSOR_H
#define VRN_VOLUMESOURCEPROCESSOR_H

#include "voreen/core/processors/processor.h"
#include "voreen/core/ports/volumeport.h"
#include "voreen/core/datastructures/volume/volume.h"
#include "voreen/core/properties/volumeurlproperty.h"
#include "voreen/core/properties/volumeinfoproperty.h"

namespace voreen {

class Volume;

/**
 * Volume data set supplier in the network.
 */
class VRN_CORE_API VolumeSource : public Processor {

public:
    VolumeSource();
    virtual Processor* create() const;

    virtual std::string getClassName() const    { return "VolumeSource";    }
    virtual std::string getCategory() const     { return "Input";           }
    virtual CodeState getCodeState() const      { return CODE_STATE_STABLE; }

    virtual void invalidate(int inv = INVALID_RESULT);

    /**
     * Loads the volume specified by filename. The loading is
     * delegated to the processor's VolumeURLProperty.
     *
     * If the data set could not be successfully loaded, an exception is thrown.
     *
     * @param filename the volume to load
     */
    void loadVolume(const std::string& filename);

    /**
     * Clears the loaded volume.
     */
    void clearVolume();

    /**
     * Assigns a volume to this processor.
     * The processor does \e not take ownership of the assigned volume.
     *
     * @param volume The volume to assign, is written to the processor's outport
     * @param owner if true, the VolumeSource takes ownership of the passed volume
     */
    void setVolume(VolumeBase* volume, bool owner = false);

    /**
     * Returns a reference to the volume of the loaded volume. May be null.
     */
    VolumeBase* getVolume() const;

protected:
    virtual void setDescriptions() {
        setDescription("Provides a single volume.");
    }

    virtual void process();
    virtual void initialize();
    virtual void deinitialize();

    VolumeURLProperty volumeURL_;
    VolumeInfoProperty volumeInfo_;

    /// The volume port the loaded data set is written to.
    VolumePort outport_;

    bool forceVolumeLoad_; ///< if true, the volume will be loaded in process

    static const std::string loggerCat_;
};

} // namespace

#endif
