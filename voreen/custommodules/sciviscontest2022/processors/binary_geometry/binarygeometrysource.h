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

#ifndef VRN_BINARYGEOMETRYSOURCE_H
#define VRN_BINARYGEOMETRYSOURCE_H

#include "voreen/core/processors/processor.h"
#include "voreen/core/ports/geometryport.h"
#include "voreen/core/properties/filedialogproperty.h"
#include "voreen/core/properties/optionproperty.h"
#include "voreen/core/properties/intproperty.h"
#include "voreen/core/properties/buttonproperty.h"

#include "voreen/core/properties/colorproperty.h"

namespace voreen {

/**
 * Reads a Voreen Binary Geometry from a file
 * and provides it as geometry through its outport.
 */
class VRN_CORE_API BinaryGeometrySource : public Processor {
public:
	BinaryGeometrySource();
    Processor* create() const override;

    std::string getClassName() const override { return "BinaryGeometrySource";   }
    std::string getCategory() const override  { return "Input";            }
    CodeState getCodeState() const override   { return CODE_STATE_EXPERIMENTAL;  }
    bool usesExpensiveComputation() const override { return true; }

protected:
    void setDescriptions() override {
        setDescription("Loads a serialized Voreen Binary Geometry (.vbge) from a file.");
    }

    void process() override;
    void initialize() override;

private:
    /**
     * Delegates to readVoreenGeometry() and assigns the returned geometry to the outport.
     */
    void readGeometry();

    /**
     * Deserializes a geometry from a Voreen Geometry file (*.vge).
     *
     * @return the read geometry, if deserialization succeeded
     *
     * @throw VoreenException if deserialization failed
     */
    Geometry* readVoreenGeometry(const std::string& filename);

    /**
     * Removed the geometry from the output and clears the file property.
     */
    void clearGeometry();

    /**
     * Triggers reloading the geometry file.
     */
    void forceReload();

    /// Adjusts the visibility of the skipItemCount_ property.
    void updatePropertyVisibility();

    FileDialogProperty geometryFile_;   ///< filename of the file containing the point information
    ButtonProperty loadGeometry_;
    ButtonProperty clearGeometry_;

    GeometryPort outport_;

    bool forceReload_;

    static const std::string loggerCat_;
};

}   //namespace

#endif // VRN_BINARYGEOMETRYSOURCE_H
