/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2018 University of Muenster, Germany.                        *
 * Visualization and Computer Graphics Group <http://viscg.uni-muenster.de>        *
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

#ifndef VRN_MESHFRUSTUMCLIPPING_H
#define VRN_MESHFRUSTUMCLIPPING_H

#include "voreen/core/processors/processor.h"
#include "voreen/core/ports/geometryport.h"

#include "voreen/core/datastructures/geometry/meshlistgeometry.h"
#include "voreen/core/properties/cameraproperty.h"
#include "voreen/core/properties/boolproperty.h"
#include "voreen/core/properties/floatproperty.h"

namespace voreen {

class MeshFrustumClipping : public Processor {
public:
    MeshFrustumClipping();
    virtual ~MeshFrustumClipping();
    virtual Processor* create() const;

    virtual std::string getClassName() const { return "MeshFrustumClipping"; }
    virtual std::string getCategory() const  { return "Geometry"; }

protected:
    virtual void setDescriptions() {
        setDescription("Makes it possible to perform clipping of geometry against the camera frustum planes.");
    }

    virtual void process();

    GeometryPort inport_;       ///< Inport for a list of mesh geometries to clip.
    GeometryPort outport_;      ///< Outport for a list of mesh geometries that were clipped.

    BoolProperty clipLeft_;       ///< Clipping against left frustum plane
    FloatProperty offsetLeft_;  ///< Left plane offset

    BoolProperty clipRight_;       ///< Clipping against right frustum plane
    FloatProperty offsetRight_; ///< Right plane offset

    BoolProperty clipTop_;       ///< Clipping against top frustum plane
    FloatProperty offsetTop_;   ///< Top plane offset

    BoolProperty clipBottom_;   ///< Clipping against bottom frustum plane
    FloatProperty offsetBottom_;///< Bottom plane offset

    BoolProperty clipNear_;       ///< Clipping against near frustum plane
    FloatProperty offsetNear_;  ///< Near plane offset

    BoolProperty clipFar_;       ///< Clipping against far frustum plane
    FloatProperty offsetFar_;   ///< Far plane offset

    CameraProperty camera_;

    MeshListGeometry geometry_;  ///< Clipped input geometry.

    /// category used in logging
    static const std::string loggerCat_;
};

} //namespace

#endif // VRN_MESHFRUSTUMCLIPPING_H
