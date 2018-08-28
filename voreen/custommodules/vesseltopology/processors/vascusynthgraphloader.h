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

#ifndef VRN_VASCUSYNTHGRAPHLOADER_H
#define VRN_VASCUSYNTHGRAPHLOADER_H

#include "voreen/core/processors/processor.h"

#include "../ports/vesselgraphport.h"
#include "voreen/core/properties/filedialogproperty.h"

namespace voreen {

class VascuSynthGraphLoader : public Processor {
public:
    VascuSynthGraphLoader();
    virtual ~VascuSynthGraphLoader();
    virtual std::string getCategory() const { return "Geometry"; }
    virtual std::string getClassName() const { return "VascuSynthGraphLoader"; }
    virtual CodeState getCodeState() const { return Processor::CODE_STATE_EXPERIMENTAL; }
    virtual Processor* create() const { return new VascuSynthGraphLoader(); }

protected:
    virtual void setDescriptions() {
        setDescription("This processor can be used to load ground truth graph data provided by VascuSynth in an xml file.");
    }

    virtual void process();

    VesselGraphPort outport_;

    // properties
    FileDialogProperty graphFilePath_;


    static const std::string loggerCat_;
};

} // namespace voreen
#endif // VRN_VASCUSYNTHGRAPHLOADER_H
