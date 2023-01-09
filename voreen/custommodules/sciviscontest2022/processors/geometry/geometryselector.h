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

#ifndef VRN_GEOMETRYSELECTOR_H
#define VRN_GEOMETRYSELECTOR_H

#include "voreen/core/datastructures/geometry/geometrysequence.h"
#include "voreen/core/ports/geometryport.h"
#include "voreen/core/processors/processor.h"
#include "voreen/core/properties/intproperty.h"
#include "voreen/core/properties/vectorproperty.h"

namespace voreen {

/**
 * Selects a single geometry out of an input geometry sequence, if the input is just a geometry it is simply returned.
 */
class VRN_CORE_API GeometrySelector : public Processor {

public:
    GeometrySelector();
    virtual Processor* create() const;

    virtual std::string getClassName() const { return "GeometrySelector"; }
    virtual std::string getCategory() const { return "Input"; }
    virtual CodeState getCodeState() const { return CODE_STATE_EXPERIMENTAL; }
    virtual bool isUtility() const { return true; }

protected:
    virtual void setDescriptions() {
        setDescription("Selects a single geometry from the input geometry sequence.");
    }

    virtual void process();
    virtual void initialize();
    virtual void deinitialize();

    virtual void adjustPropertiesToInput();

    GeometryPort inport_;   ///< Inport for the geometry sequence.
    GeometryPort outport_;  ///< The geometry port the selected geometry.

    IntProperty geometryID_;   ///< id of the selected geometry

    static const std::string loggerCat_;
};

} // namespace

#endif
