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

#ifndef VRN_GEOMETRYMERGE_H
#define VRN_GEOMETRYMERGE_H

#include "voreen/core/processors/processor.h"

#include "voreen/core/ports/geometryport.h"

#include "voreen/core/properties/boolproperty.h"

namespace voreen {

/**
 * Closes all holes of the input geometry.
 */
class VRN_CORE_API GeometryMerge : public Processor {
public:
    GeometryMerge();
    virtual ~GeometryMerge();
    virtual Processor* create() const;

    virtual std::string getClassName() const { return "GeometryMerge";          }
    virtual std::string getCategory() const  { return "Geometry";               }
    virtual CodeState getCodeState() const   { return CODE_STATE_EXPERIMENTAL;  }

protected:
    virtual void setDescriptions() {
        setDescription("Merges a geometry sequence into one geometry.");
    }

    virtual void process();

    GeometryPort inport_;        ///< Inport for a geometry sequence to merge.
    GeometryPort outport_;       ///< Outport for a geometry sequence that was merged.

    BoolProperty enabled_;       ///< Determines whether the merging should be performed.

    /// category used in logging
    static const std::string loggerCat_;
};

} //namespace

#endif // VRN_GEOMETRYCLOSE_H
