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

#ifndef VRN_POINTSEGMENTLISTRENDERER_H
#define VRN_POINTSEGMENTLISTRENDERER_H

#include "voreen/core/processors/geometryrendererbase.h"
#include "voreen/core/properties/propertyvector.h"
#include "voreen/core/datastructures/geometry/glmeshgeometry.h"

#include "voreen/core/properties/optionproperty.h"
#include "voreen/core/properties/boolproperty.h"
#include "voreen/core/properties/intproperty.h"
#include "voreen/core/properties/floatproperty.h"
#include "voreen/core/properties/vectorproperty.h"
#include "voreen/core/properties/colorproperty.h"

namespace voreen {

class Geometry;

template<typename T>
class PointSegmentListGeometry;

/**
 * Renders a PointSegmentListGeometry<tgt::vec3> which is passed through a geometry inport.
 * The point list can be rendered as points, spheres or illuminated spheres.
 */
class VRN_CORE_API PointSegmentListRenderer : public GeometryRendererBase {
public:
    PointSegmentListRenderer();
    virtual ~PointSegmentListRenderer();
    virtual std::string getCategory() const { return "Geometry"; }
    virtual std::string getClassName() const { return "PointSegmentListRenderer"; }
    virtual CodeState getCodeState() const { return Processor::CODE_STATE_STABLE; }
    virtual Processor* create() const { return new PointSegmentListRenderer(); }

    virtual void serialize(Serializer& s) const override;
    virtual void deserialize(Deserializer& s) override;

protected:
    virtual void setDescriptions() {
        setDescription("Renders a list of segments, each consisting of points (PointSegmentListGeometryVec3), as point primitives, spheres or illuminated spheres. Each segment can be assigned a different color.\
<p>see GeometrySource</p>");
    }

    /// Actual rendering is implemented in render()-method
    virtual void process();

    virtual void initialize();
    virtual void deinitialize();

    /**
     * Actually renders the geometry.
     */
    virtual void render();

    /**
     * Returns the enclosing Bounding Box in world coordinates.
     */
    virtual tgt::Bounds getBoundingBox() const;

    /**
     * Property callback for invalidating the triangle mesh geometry on property changes.
     */
    void invalidateGeometry();

    /**
     * Updates the list of segments with enabled status and color
     */
    void updateSegmentList();

    /**
     * Updates the propertys for the configuration of segments
     */
    void updateSegmentListUI();

    /**
     * Updates the propertys for the configuration of segments color and enabled status
     */
    void updateSegmentUI();

    /// mesh for speeding up rendering
    std::unique_ptr<GlMeshGeometryBase> mesh_;

    tgt::Shader* sphereShader_; //< sphere rendering shader with illumination
    tgt::Shader* circleShader_; //< rendering points as circles (old GL_POINT_SMOOTH simulation)

    // properties
    StringOptionProperty coordinateSystem_;
    StringOptionProperty renderingPrimitiveProp_;
    BoolProperty applyUniformColor_;
    ColorProperty color_;
    BoolProperty depthTest_;
    FloatProperty pointSize_;
    BoolProperty pointSmooth_;
    FloatProperty sphereDiameter_;
    IntProperty sphereSlicesStacks_;

    OptionProperty<int> segmentList_;
    ColorProperty segmentColor_;
    BoolProperty segmentEnabled_;

    GeometryPort geometryInport_;

    struct Segment :public Serializable {
        tgt::vec4 color_;
        bool enabled_;

        virtual void serialize(Serializer& s) const override;
        virtual void deserialize(Deserializer& s) override;
    };

    std::vector<Segment> segments_;
    const PointSegmentListGeometry<tgt::vec3>* segmentListGeom_;

    static const std::string loggerCat_;
};

} // namespace voreen

#endif //VRN_POINTSEGMENTLISTRENDERER
