/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2021 University of Muenster, Germany,                        *
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

#ifndef VRN_VORTEXLISTSELECTOR_H
#define VRN_VORTEXLISTSELECTOR_H

#include "voreen/core/processors/processor.h"
#include "voreen/core/ports/geometryport.h"
#include "voreen/core/properties/optionproperty.h"
#include "voreen/core/properties/numeric/intervalproperty.h"
#include "voreen/core/properties/string/stringlistproperty.h"

#include "custommodules/sciviscontest2020/ports/vortexport.h"

namespace voreen {
class VortexListSelector : public Processor {
private:
    enum RotationOptions {
        OPTION_CW, //Clockwise
        OPTION_CCW, //Counterclockwise
        OPTION_B, //Both
    };

public:
    VortexListSelector();

    Processor* create() const override
    {
        return new VortexListSelector();
    }
    std::string getClassName() const override
    {
        return "VortexListSelector";
    }
    std::string getCategory() const override
    {
        return "Vortex Processing";
    }

    bool isReady() const override
    {
        return _inportVortexCollection.isReady();
    }

    static void Process( const VortexCollection& vortices, const std::vector<int>& runs, int firstTimestep, int lastTimestep, int minLength,RotationOptions rot, std::vector<Vortex>& outVortexList );

private:
    void process() override;
    void updatePropertyCorelineLength();

    VortexCollectionPort _inportVortexCollection;
    VortexListPort _outportVortexList;
    GeometryPort _outportGeometry;

    StringListProperty _propertyMembers;
    IntIntervalProperty _propertyTimesteps;
    IntProperty _propertyCorelineLength;
    OptionProperty<RotationOptions> _Rotation;
};

}

#endif // VRN_VORTEXLISTSELECTOR_H
