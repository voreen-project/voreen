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

#ifndef VRN_BINARYGEOMETRYSAVE_H
#define VRN_BINARYGEOMETRYSAVE_H

#include "voreen/core/processors/processor.h"

#include "voreen/core/ports/geometryport.h"
#include "voreen/core/properties/boolproperty.h"
#include "voreen/core/properties/buttonproperty.h"
#include "voreen/core/properties/filedialogproperty.h"

namespace voreen {

class VRN_CORE_API BinaryGeometrySave : public Processor {
public:
	BinaryGeometrySave();
    Processor* create() const override;

    std::string getClassName() const override	{ return "BinaryGeometrySave";    }
    std::string getCategory() const override 	{ return "Output";          }
    CodeState getCodeState() const override  	{ return CODE_STATE_EXPERIMENTAL; }
    bool isEndProcessor() const override     	{ return true;              }

    bool isReady() const override { return true; }
    void invalidate(int inv = 1) override;

protected:
    void setDescriptions() override {
        setDescription("Writes the input geometry to a Voreen Binary Geometry file (.vbge).");
    }

    void process() override;

    void saveFile();

    GeometryPort inport_;

    FileDialogProperty fileProp_;
    ButtonProperty saveButton_;
    BoolProperty continousSave_;

    static const std::string loggerCat_;
};

} // namespace voreen

#endif // VRN_BINARYGEOMETRYSAVE_H
