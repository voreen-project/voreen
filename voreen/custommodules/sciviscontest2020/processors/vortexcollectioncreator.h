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

#ifndef VRN_VORTEXCOLLECTIONCREATOR_H
#define VRN_VORTEXCOLLECTIONCREATOR_H

#include "voreen/core/processors/processor.h"
#include "voreen/core/properties/numeric/intervalproperty.h"
#include "voreen/core/properties/string/stringlistproperty.h"

#include "modules/ensembleanalysis/ports/ensembledatasetport.h"
#include "custommodules/sciviscontest2020/ports/vortexport.h"

#include <chrono>

namespace voreen {

class VortexCollectionCreator : public Processor {
public:
    VortexCollectionCreator();

    Processor* create() const override          { return new VortexCollectionCreator(); }
    std::string getClassName() const override   { return "VortexCollectionCreator";     }
    std::string getCategory() const override    { return "Vortex Processing";           }
    bool isReady() const override {
        return _inportEnsemble.isReady();
    }

private:
    void process() override {}
    void updateButton();

    EnsembleDatasetPort _inportEnsemble;
    VortexCollectionPort _outportVortexCollection;

    StringListProperty _propertySelectedMembers;
    IntIntervalProperty _propertyTimestepInterval;
    IntProperty _propertyCorelineLength;
    ButtonProperty _propertyUpdateButton;

    FileDialogProperty _propertyFileDialog;
    ButtonProperty _propertySaveButton;
};

}

#endif // VRN_VORTEXCOLLECTIONCREATOR_H
