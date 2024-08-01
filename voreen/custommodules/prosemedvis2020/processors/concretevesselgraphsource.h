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

#ifndef VRN_CONCRETEVESSELGRAPHSOURCE_H
#define VRN_CONCRETEVESSELGRAPHSOURCE_H

#include "voreen/core/processors/processor.h"

#include "voreen/core/properties/filedialogproperty.h"
#include "voreen/core/properties/buttonproperty.h"

#include "../ports/concretevesselgraphport.h"

namespace voreen {

class ConcreteVesselGraphSource : public Processor {
public:
    ConcreteVesselGraphSource();
    virtual ~ConcreteVesselGraphSource();
    virtual std::string getCategory() const { return "Input"; }
    virtual std::string getClassName() const { return "ConcreteVesselGraphSource"; }
    virtual CodeState getCodeState() const { return Processor::CODE_STATE_TESTING; }
    virtual Processor* create() const { return new ConcreteVesselGraphSource(); }

protected:
    virtual void setDescriptions() {
        setDescription("This processor can be used to load concrete vessel graph files that have been previously saved using <b>ConcreteVesselGraphCreator</b>. "
            "ConcreteVesselgraphs are serialized in a custom (but simple) json format."
        );
        graphFilePath_.setDescription("Path to the *.json file to be loaded.");
    }

    virtual void process();

    ConcreteVesselGraphPort outport_;

    // properties
    FileDialogProperty graphFilePath_;
    ButtonProperty reload_;

};

} // namespace voreen
#endif 
