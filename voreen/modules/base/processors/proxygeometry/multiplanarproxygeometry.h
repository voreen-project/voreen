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

#ifndef VRN_MULTIPLANARPROXYGEOMETRY_H
#define VRN_MULTIPLANARPROXYGEOMETRY_H

#include "voreen/core/processors/processor.h"
#include "voreen/core/ports/volumeport.h"
#include "voreen/core/ports/genericport.h"
#include "voreen/core/ports/geometryport.h"
#include "voreen/core/properties/floatproperty.h"
#include "voreen/core/datastructures/geometry/meshgeometry.h"

namespace voreen {

/**
 * Provides a volume proxy geometry consisting of three orthogonal, axis-aligned slices.
 * Each slice is represented by a front and a back face. The slice positions are specified in
 * normalized coordinates.
 */
class VRN_CORE_API MultiPlanarProxyGeometry : public Processor {
public:
    MultiPlanarProxyGeometry();
    virtual Processor* create() const;

    virtual std::string getClassName() const { return "MultiplanarProxyGeometry"; }
    virtual std::string getCategory() const  { return "Volume Proxy Geometry";    }
    virtual CodeState getCodeState() const   { return CODE_STATE_STABLE;          }

protected:
    virtual void setDescriptions() {
        setDescription("Provides a proxy geometry consisting of three orthogonal, axis-aligned slices. Each slice is represented by a front and a back face. The slice positions are specified in (normalized) texture coordinates.\
<p><span style=\"font-weight: bold\">Note:</span> Since generated geometry is non-convex, it cannot be used in a raycasting pipeline.</p>");
    }

    void process();

    VolumePort inport_;      ///< Inport for the dataset the proxy is generated for.
    GeometryPort outport_;   ///< Outport for the generated mesh proxy geometry.

    FloatProperty slicePosX_;
    FloatProperty slicePosY_;
    FloatProperty slicePosZ_;
};

} // namespace voreen

#endif // VRN_MULTIPLANARPROXYGEOMETRY_H
