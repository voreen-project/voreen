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

#ifndef VRN_STREAMLINETOGEOMETRY_H
#define VRN_STREAMLINETOGEOMETRY_H

#include "voreen/core/processors/processor.h"

#include "../../ports/streamlinelistport.h"
#include "voreen/core/ports/geometryport.h"

namespace voreen {

/**
 * This processor converts a StreamlineList into a bounding box which can be rendered via the BoundingBoxRenderer.
 * The output geometry consists of only two points (LLF and URB) of the StreamlineList bounds.
 *
 * @see StreamlineCreator
 * @see Streamline
 */
class VRN_CORE_API StreamlineToGeometry : public Processor {
public:
    StreamlineToGeometry();
    virtual ~StreamlineToGeometry();
    virtual Processor* create() const         {return new StreamlineToGeometry(); }

    virtual std::string getClassName() const  { return "StreamlineToGeometry";  }
    virtual std::string getCategory() const   { return "Converter";          }
    virtual CodeState getCodeState() const    { return CODE_STATE_EXPERIMENTAL; }

protected:
    virtual void setDescriptions() {
        setDescription("This processor converts a StreamlineList into a geometry that can be rendered by the GeometryRenderer. "
                       "Alternatively, the output can be used as foreground seeds for the random walker. "
                       "In order to do so, choose PointSegmentListGeometry and Voxel Coordinates.");
    }

    virtual void process();

    //--------------
    //  Member
    //--------------
private:

    enum GeometryType {
        GEOMETRY_GLMESH,
        GEOMETRY_POINTSEGMENTLIST
    };

    enum TargetCoordinateSystem {
        VOXEL_COORDINATES,
        WORLD_COORDINATES,
    };

    StreamlineListPort inport_;
    GeometryPort outport_;
    OptionProperty<GeometryType> geometryType_;
    OptionProperty<TargetCoordinateSystem> targetCoordinateSystem_;

    static const std::string loggerCat_;
};

}   //namespace

#endif
