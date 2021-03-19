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

#ifndef VRN_ENSEMBLEINFORMATION_H
#define VRN_ENSEMBLEINFORMATION_H

#include "voreen/core/processors/processor.h"
#include "voreen/core/properties/intproperty.h"
#include "voreen/core/properties/floatproperty.h"

#include "../ports/ensembledatasetport.h"

namespace voreen {

/**
 * Displays meta data of the incoming ensemble as read-only properties.
 */
class VRN_CORE_API EnsembleInformation : public Processor {
public:
    EnsembleInformation();
    Processor* create() const;

    std::string getClassName() const    { return "EnsembleInformation"; }
    std::string getCategory() const     { return "Utility"; };
    CodeState getCodeState() const      { return CODE_STATE_TESTING; }
    bool isUtility() const              { return true; }

protected:
    virtual void setDescriptions() {
        setDescription("Displays meta data of the incoming ensemble as read-only properties");
    }

    virtual void process();
    virtual void adjustPropertiesToInput();

private:

    EnsembleDatasetPort ensemblePort_;

    StringOptionProperty fieldNames_;
    FloatProperty minValue_;
    FloatProperty maxValue_;
    FloatProperty minMagnitude_;
    FloatProperty maxMagnitude_;
    IntProperty numChannels_;

};

} // namespace

#endif // VRN_ENSEMBLEINFORMATION_H
