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

#ifndef VRN_GEOMETRYSEQUENCECREATOR_H
#define VRN_GEOMETRYSEQUENCECREATOR_H

#include "voreen/core/processors/processor.h"

#include "voreen/core/ports/geometryport.h"

#include "voreen/core/properties/boolproperty.h"
#include "voreen/core/properties/boundingboxproperty.h"

namespace voreen {

/**
 * Can be used to convert N geometries into a GeometrySequence.
 * @see GeoemtrySequence
 */
class VRN_CORE_API GeometrySequenceCreator : public Processor {
public:
    GeometrySequenceCreator();
    virtual ~GeometrySequenceCreator();
    virtual Processor* create() const;

    virtual std::string getClassName() const { return "GeometrySequenceCreator"; }
    virtual std::string getCategory() const  { return "Geometry";         }
    virtual CodeState getCodeState() const   { return CODE_STATE_STABLE;  }

protected:
    virtual void setDescriptions() {
        setDescription("Can be used to convert N geometries into a GeometrySequence. <br>" \
                       "If a GeometrySequence is added, it will be split into its components. <br>" \
                       "The clipping takes place in world-coordinates.");
    }

    virtual void process();
    virtual void adjustPropertiesToInput();
    void adjustPropertyVisibility();
    void invalidateOutput();

    GeometryPort inport_;        ///< Inport for the N geometries
    GeometryPort outport_;       ///< Outport for the created geometrysequence

    bool invalidOutput_;     ///< is set in adjustPropertiesToInput()

    BoolProperty enableClippingProp_;               ///< enable clipping
    FloatBoundingBoxProperty clippingBBoxProp_;     ///< the clipping itself

    /// category used in logging
    static const std::string loggerCat_;
};

} //namespace

#endif // VRN_GEOMETRYSEQUENCECREATOR_H
