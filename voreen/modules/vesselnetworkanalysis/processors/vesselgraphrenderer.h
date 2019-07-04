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

#ifndef VRN_VESSELGRAPHRENDERER_H
#define VRN_VESSELGRAPHRENDERER_H

#include "voreen/core/processors/geometryrendererbase.h"

#include "../ports/vesselgraphport.h"
#include "voreen/core/datastructures/geometry/glmeshgeometry.h"
#include "voreen/core/properties/colorproperty.h"
#include "voreen/core/properties/shaderproperty.h"
#include "voreen/core/properties/floatproperty.h"
#include "voreen/core/properties/intproperty.h"
#include "voreen/core/properties/boolproperty.h"
#include "voreen/core/properties/optionproperty.h"
#include "voreen/core/properties/lightsourceproperty.h"
#include "voreen/core/properties/transfunc/1d/1dkeys/transfunc1dkeysproperty.h"
#include "voreen/core/datastructures/volume/volumeatomic.h"

#include <functional>

namespace voreen {

// Helper class that collects data to display it in a transfer function property
template<typename T>
struct ArbitraryHistogramTF {
    ArbitraryHistogramTF(const std::string& id, const std::string& guiName, std::function<float(const T& val)> getVal)
        : transfunc_(id, guiName)
        , valueVol_(nullptr)
        , getValue(getVal)
    {
    }

    ArbitraryHistogramTF(const ArbitraryHistogramTF&) = delete;
    ArbitraryHistogramTF(ArbitraryHistogramTF&& other)
        : transfunc_(other.transfunc_)
        , valueVol_(other.valueVol_.release())
        , getValue(other.getValue)
    {
    }


    // Set new data using the supplied iterators and the getVal function.
    template<class Iterator>
    void setData(Iterator first, Iterator last) {
        VolumeRAM* volram = new VolumeRAM_Float(tgt::svec3(std::distance(first, last),1,1));
        //float* data = static_cast<float*>(volram->getData());
        for(size_t i = 0; first != last; ++first, ++i) {
            //data[i] = static_cast<float>(nodes[i].getConnectivity());
            float val = getValue(*first);
            tgtAssert(!std::isnan(val), "nan edge value");
            volram->setVoxelNormalized(val, i, 0, 0);
        }
        transfunc_.setVolume(nullptr);
        valueVol_ = std::unique_ptr<Volume>(new Volume(volram, tgt::vec3::one, tgt::vec3::zero));
        transfunc_.setVolume(valueVol_.get());
    }


    TransFunc1DKeysProperty& getProperty() {
        return transfunc_;
    }

    // Conversion function from arbitrary values to float (which can be displayed in the TF)
    const std::function<float(const T& val)> getValue;

private:
    // Holds the values to create the histogram from
    std::unique_ptr<Volume> valueVol_;
    TransFunc1DKeysProperty transfunc_;
};

// Renders a VesselGraph as a geometry using spheres for nodes and cylinders for edges between them.
// Features of nodes and edges can be displayed using color maps.
class VesselGraphRenderer : public GeometryRendererBase {
public:
    VesselGraphRenderer();
    virtual ~VesselGraphRenderer();
    virtual std::string getCategory() const { return "Geometry"; }
    virtual std::string getClassName() const { return "VesselGraphRenderer"; }
    virtual CodeState getCodeState() const { return Processor::CODE_STATE_EXPERIMENTAL; }
    virtual Processor* create() const { return new VesselGraphRenderer(); }

    virtual tgt::Bounds getBoundingBox() const;

    virtual void render();

protected:
    virtual void setDescriptions() {
        setDescription("Renders Vessel-Graphs.");
    }

    void adjustLightingPropertyVisibility();
    void setLightingUniforms(tgt::Shader& shader) const;

    /// Actual rendering is implemented in render()-method
    virtual void process();

    virtual void initialize();
    virtual void deinitialize();

    enum EdgeProperty {
        LENGTH,
        DISTANCE,
        CURVENESS,
        STRAIGHTNESS,
        VOLUME,
        AVG_CROSS_SECTION,
        MIN_RADIUS_AVG,
        MIN_RADIUS_STD_DEV,
        AVG_RADIUS_AVG,
        AVG_RADIUS_STD_DEV,
        MAX_RADIUS_AVG,
        MAX_RADIUS_STD_DEV,
        ROUNDNESS_AVG,
        ROUNDNESS_STD_DEV,
        ELONGATION,
        RELATIVE_BULGE_SIZE
    };

    void adaptToNewInput();
    void showActiveEdgeTF();

    /// triangle mesh for rendering a node
    GlMeshGeometryUInt16Normal nodeTriangleMesh_;

    /// triangle mesh for rendering a node at the border of the sample
    GlMeshGeometryUInt16Normal borderNodeTriangleMesh_;

    /// triangle mesh for rendering an edge
    GlMeshGeometryUInt16Normal edgeTriangleMesh_;

    // properties
    BoolProperty enabled_;
    ShaderProperty shader_;
    ShaderProperty activeEdgeShader_;
    BoolProperty renderNodeRadii_;
    FloatProperty nodeRadiusMultiplier_;
    FloatProperty edgeCrossSectionMultiplier_;
    ArbitraryHistogramTF<VesselGraphNode> nodeDegreeTF_;
    OptionProperty<EdgeProperty> activeEdgeProperty_;
    std::map<EdgeProperty, ArbitraryHistogramTF<VesselGraphEdge>> edgeTFs_;
    IntProperty activeEdgeID_;
    ColorProperty activeEdgeColor_;

    // Lighting related properties
    BoolProperty enableLighting_;
    LightSourceProperty lightPosition_; ///< The position of the light source in world coordinates
    ColorProperty lightAmbient_;        ///< The light source's ambient color
    ColorProperty lightDiffuse_;        ///< The light source's diffuse color
    ColorProperty lightSpecular_;       ///< The light source's specular color
    FloatProperty materialShininess_;   ///< The material's specular exponent

    VesselGraphPort graphInport_;

    static const std::string loggerCat_;

};

} // namespace voreen
#endif // VRN_VESSELGRAPHRENDERER_H
