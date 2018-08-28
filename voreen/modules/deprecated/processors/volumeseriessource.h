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

#ifndef VRN_VOLUMESERIESSOURCE_H
#define VRN_VOLUMESERIESSOURCE_H

#include "voreen/core/processors/processor.h"
#include "voreen/core/ports/volumeport.h"
#include "voreen/core/datastructures/volume/volume.h"
#include "voreen/core/properties/intproperty.h"
#include "voreen/core/properties/filedialogproperty.h"

namespace voreen {

/**
 * Supplies a single volume out of a series of (time-varying) volume files. The series is
 * defined in a .sdat file, which is like a .dat file, just with multiple "ObjectFileName:"
 * entries, each corresponding to a single volume file. The data is loaded on demand based on
 * the "time step" property.
 */
class VolumeSeriesSource : public Processor {
public:
    VolumeSeriesSource();
    ~VolumeSeriesSource();
    virtual Processor* create() const;

    virtual std::string getCategory() const  { return "Data Source";        }
    virtual std::string getClassName() const { return "VolumeSeriesSource"; }
    virtual CodeState getCodeState() const   { return CODE_STATE_OBSOLETE;   }

protected:
    virtual void setDescriptions() {
        setDescription("");
    }

    virtual void initialize();
    virtual void process();

    void openSeries();
    void loadStep();

    Volume* volume_;

    FileDialogProperty filename_;
    IntProperty step_;

    VolumePort outport_;

    std::vector<std::string> files_;
    bool needUpload_;
    float spreadMin_, spreadMax_;

    static const std::string loggerCat_;
};

} // namespace

#endif // VRN_VOLUMESERIESSOURCE
