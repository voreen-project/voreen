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

#ifndef VRN_VOLUMELISTSOURCE_H
#define VRN_VOLUMELISTSOURCE_H

#include "voreen/core/processors/processor.h"
#include "voreen/core/ports/genericport.h"
#include "voreen/core/properties/volumeurllistproperty.h"

namespace voreen {

class Volume;
class VolumeList;

/**
 * Loads multiple volumes and provides them
 * as VolumeList through its outport.
 */
class VRN_CORE_API VolumeListSource : public Processor {

public:
    VolumeListSource();
    virtual ~VolumeListSource();
    virtual Processor* create() const;

    virtual std::string getClassName() const  { return "VolumeListSource"; }
    virtual std::string getCategory() const   { return "Input";            }
    virtual CodeState getCodeState() const    { return CODE_STATE_STABLE;  }

    virtual void invalidate(int inv = INVALID_RESULT);

    /**
     * Assigns a volume list to this processor.
     *
     * @param owner if true, the processor takes ownership of the passed volumes
     * @param selected if true, the passed volumes will be selected
     */
    void setVolumeList(VolumeList* collection, bool owner = false, bool selected = true);

    /**
     * Returns the currently assigned volume list.
     */
    VolumeList* getVolumeList() const;

    /**
     * Loads volumes from the passed URL and adds them to the output list.
     * If the data set could not be successfully loaded, an exception is thrown.
     *
     * @param selected if true, the loaded volumes will be selected
     */
    void loadVolumes(const std::string& url, bool selected = true);

protected:
    virtual void setDescriptions() {
        setDescription("Loads multiple volumes and provides them as VolumeList.");
    }

    virtual void process();
    virtual void initialize();

    /// The volume port the loaded data set is written to.
    VolumeListPort outport_;

    /// Property storing the loaded volume collection.
    VolumeURLListProperty volumeURLList_;

    static const std::string loggerCat_;
};

} // namespace

#endif
