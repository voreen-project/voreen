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

#include "vesselgraphrenderer.h"

#include "voreen/core/datastructures/callback/lambdacallback.h"
#include "tgt/quaternion.h"

#include <boost/iterator/filter_iterator.hpp>

namespace voreen {



const std::string VesselGraphRenderer::loggerCat_("voreen.vesselgraphrenderer");


#define ADD_EDGE_PROPERTY(Type, name, displayname) \
    activeEdgeProperty_.addOption("edge" #name, displayname, EdgeProperty::Type); \
    edgeTFs_.emplace(EdgeProperty::Type, ArbitraryHistogramTF<VesselGraphEdge>("edge" #name "TF", displayname " TF", [] (const VesselGraphEdge& edge) { \
                return static_cast<float>(edge.get ## name ()); \
            }));\

VesselGraphRenderer::VesselGraphRenderer()
    : GeometryRendererBase()
    , enabled_("enabled", "Enabled", true)
    , shader_("shader", "Shader", "vesselgraphrenderer.frag", "vesselgraphrenderer.vert")
    , activeEdgeShader_("activeEdgeShader", "Active Edge Shader", "vesselgraphrenderer.frag", "vesselgraphrenderer_activeedge.vert")
    , nodeRadiusMultiplier_("nodeRadiusMultiplier_", "Node Radius Multiplier", 0.5, 0, 10)
    , edgeCrossSectionMultiplier_("edgeCrossSectionMultiplier_", "Edge Cross Section Multiplier", 0.3, 0, 10)
    , nodeDegreeTF_("nodeDegreeTf", "Node Degree", [] (const VesselGraphNode& node) {
                return static_cast<float>(node.getDegree());
            })
    , activeEdgeProperty_("activeEdgeProperty", "Active edge property")
    , edgeTFs_()
    , renderNodeRadii_("renderNodeRadii", "Render Node Radii", true)
    , activeEdgeID_("activeEdgeID", "Active Edge ID", -1, -1, std::numeric_limits<int>::max()/2) // This is stupid. Somehow Qt divides by 0 if we set INT_MAX and causes a SIGFPE...
    , activeEdgeColor_("activeEdgeColor", "Active Edge Color", tgt::Color(0.25f, 0.75f, 0.f, 1.f))
    , enableLighting_("enableLighting", "Enable Lighting", false)
    , lightPosition_("lightPosition", "Light Source Position", tgt::vec4(2.3f, 1.5f, 1.5f, 1.f), tgt::vec4(-10000), tgt::vec4(10000))
    , lightAmbient_("lightAmbient", "Ambient Light", tgt::Color(0.4f, 0.4f, 0.4f, 1.f))
    , lightDiffuse_("lightDiffuse", "Diffuse Light", tgt::Color(0.8f, 0.8f, 0.8f, 1.f))
    , lightSpecular_("lightSpecular", "Specular Light", tgt::Color(0.6f, 0.6f, 0.6f, 1.f))
    , materialShininess_("materialShininess", "Shininess", 60.f, 0.1f, 128.f)
    , graphInport_(Port::INPORT, "graph.input", "Graph Input")
{
    ADD_EDGE_PROPERTY(LENGTH, Length, "Length")
    ADD_EDGE_PROPERTY(DISTANCE, Distance, "Distance")
    ADD_EDGE_PROPERTY(CURVENESS, Curveness, "Curveness")
    ADD_EDGE_PROPERTY(STRAIGHTNESS, Straightness, "Straightness")
    ADD_EDGE_PROPERTY(VOLUME, Volume, "Volume")
    ADD_EDGE_PROPERTY(AVG_CROSS_SECTION, AvgCrossSection, "Average Cross Section")
    ADD_EDGE_PROPERTY(MIN_RADIUS_AVG, MinRadiusAvg, "MinRadius Mean")
    ADD_EDGE_PROPERTY(MIN_RADIUS_STD_DEV, MinRadiusStdDeviation, "MinRadius Standard Deviation")
    ADD_EDGE_PROPERTY(AVG_RADIUS_AVG, AvgRadiusAvg, "AvgRadius Mean")
    ADD_EDGE_PROPERTY(AVG_RADIUS_STD_DEV, AvgRadiusStdDeviation, "AvgRadius Standard Deviation")
    ADD_EDGE_PROPERTY(MAX_RADIUS_AVG, MaxRadiusAvg, "MaxRadius Mean")
    ADD_EDGE_PROPERTY(MAX_RADIUS_STD_DEV, MaxRadiusStdDeviation, "MaxRadius Standard Deviation")
    ADD_EDGE_PROPERTY(ROUNDNESS_AVG, RoundnessAvg, "Roundness Mean")
    ADD_EDGE_PROPERTY(ROUNDNESS_STD_DEV, RoundnessStdDeviation, "Roundness Standard Deviation")
    ADD_EDGE_PROPERTY(ELONGATION, Elongation, "Elongation")
    ADD_EDGE_PROPERTY(RELATIVE_BULGE_SIZE, RelativeBulgeSize, "Relative Bulge Size")

    addProperty(enabled_);
    addProperty(renderNodeRadii_);
    addProperty(nodeRadiusMultiplier_);
    addProperty(edgeCrossSectionMultiplier_);
    addProperty(nodeDegreeTF_.getProperty());
    addProperty(activeEdgeProperty_);
        ON_CHANGE(activeEdgeProperty_, VesselGraphRenderer, showActiveEdgeTF);
    for(auto& pair : edgeTFs_) {
        addProperty(pair.second.getProperty());
    }
    addProperty(activeEdgeID_);
    addProperty(activeEdgeColor_);
    addProperty(shader_);
    addProperty(activeEdgeShader_);

        addProperty(enableLighting_);
            enableLighting_.setGroupID("lighting");
            ON_CHANGE(enableLighting_, VesselGraphRenderer, adjustLightingPropertyVisibility);
        addProperty(lightPosition_);
            lightPosition_.setGroupID("lighting");
        addProperty(lightAmbient_);
            lightAmbient_.setGroupID("lighting");
        addProperty(lightDiffuse_);
            lightDiffuse_.setGroupID("lighting");
        addProperty(lightSpecular_);
            lightSpecular_.setGroupID("lighting");
        addProperty(materialShininess_);
            materialShininess_.setGroupID("lighting");
    setPropertyGroupGuiName("lighting", "Lighting Parameters");
    adjustLightingPropertyVisibility();

    addPort(graphInport_);

    showActiveEdgeTF();
}

VesselGraphRenderer::~VesselGraphRenderer() {
}

void VesselGraphRenderer::adjustLightingPropertyVisibility() {
    bool visible = enableLighting_.get();
    lightPosition_.setVisibleFlag(visible);
    lightAmbient_.setVisibleFlag(visible);
    lightDiffuse_.setVisibleFlag(visible);
    lightSpecular_.setVisibleFlag(visible);
    materialShininess_.setVisibleFlag(visible);
}

void VesselGraphRenderer::process() {
    tgtAssert(graphInport_.isReady(), "inport not ready");
    if(graphInport_.hasChanged()) {
        adaptToNewInput();
    }
}

void VesselGraphRenderer::initialize() {
    GeometryRendererBase::initialize();
    nodeTriangleMesh_.setSphereGeometry(1.0f /*radius*/, tgt::vec3::zero, tgt::vec4::one, 8);
    borderNodeTriangleMesh_.setCuboidGeometry(1,1,1);
    edgeTriangleMesh_.setCylinderGeometry(tgt::vec4::one, 1, 1, 1, 8, 8, false, false);
    shader_.rebuild();
    activeEdgeShader_.rebuild();
}

void VesselGraphRenderer::deinitialize() {
    GeometryRendererBase::deinitialize();
}

void VesselGraphRenderer::setLightingUniforms(tgt::Shader& shader) const {
    // enable or disable lighting
    shader.setUniform("lightingEnabled_", enableLighting_.get());

    // set lighting and material parameters
    std::string prefix = "lightSource_.";
    shader.setUniform(prefix+"position_", lightPosition_.get());
    shader.setUniform(prefix+"ambientColor_", lightAmbient_.get().xyz());
    shader.setUniform(prefix+"diffuseColor_", lightDiffuse_.get().xyz());
    shader.setUniform(prefix+"specularColor_", lightSpecular_.get().xyz());
    shader.setUniform("material_.shininess_", materialShininess_.get());

    LGL_ERROR;
}

void VesselGraphRenderer::render() {
    tgtAssert(graphInport_.isReady(), "inport not ready");

    if(!enabled_.get()) {
        return;
    }

    const VesselGraph* graph = graphInport_.getData();
    if(!graph) {
        LWARNING("No input graph");
        return;
    }

    MatStack.matrixMode(tgt::MatrixStack::MODELVIEW);

    tgt::Shader* shader = shader_.getShader();
    if(!shader) {
        LWARNING("No shader");
        return;
    }

    // activate shader
    shader->activate();

    tgt::TextureUnit transferUnit;
    transferUnit.activate();
    if(!nodeDegreeTF_.getProperty().get()) {
        LWARNING("No TF");
        return;
    }
    nodeDegreeTF_.getProperty().get()->getTexture()->bind();
    nodeDegreeTF_.getProperty().get()->setUniform(shader, "transferFunc_", "transferFuncTex_", transferUnit.getUnitNumber());

    float boundsDiagonal = tgt::length(graph->getBounds().diagonal());
    if(boundsDiagonal == 0) {
        boundsDiagonal = 1;
    }
    float nodeRadius = boundsDiagonal * 0.01 * nodeRadiusMultiplier_.get();
    float edgeCrossSection = boundsDiagonal * 0.01 * edgeCrossSectionMultiplier_.get();

    setLightingUniforms(*shader);

    for(const VesselGraphNode& node : graph->getNodes()) {
        MatStack.pushMatrix();
        MatStack.translate(node.pos_);
        float r;
        if(renderNodeRadii_.get()) {
            r = std::max(node.getRadius(), nodeRadius);
        } else {
            r = nodeRadius;
        }
        MatStack.scale(tgt::vec3(r));
        // set matrix stack uniforms
        shader->setUniform("modelViewMatrixStack_", MatStack.getModelViewMatrix());
        shader->setUniform("modelViewProjectionMatrixStack_", MatStack.getProjectionMatrix() * MatStack.getModelViewMatrix());
        tgt::mat4 viewInverse;
        MatStack.getModelViewMatrix().invert(viewInverse);
        shader->setUniform("modelViewMatrixInverseStack_", viewInverse);
        shader->setUniform("property_", nodeDegreeTF_.getValue(node));

        if(node.isAtSampleBorder_) {
            borderNodeTriangleMesh_.render(GL_TRIANGLE_STRIP);
        } else {
            nodeTriangleMesh_.render(GL_TRIANGLE_STRIP);
        }

        MatStack.popMatrix();
    }

    auto it = edgeTFs_.find(activeEdgeProperty_.getValue());
    tgtAssert(it != edgeTFs_.end(), "Missing tf in map");
    ArbitraryHistogramTF<VesselGraphEdge>& currentEdgeTF = it->second;
    if(!currentEdgeTF.getProperty().get()) {
        LWARNING("No TF");
        return;
    }
    currentEdgeTF.getProperty().get()->getTexture()->bind();
    currentEdgeTF.getProperty().get()->setUniform(shader, "transferFunc_", "transferFuncTex_", transferUnit.getUnitNumber());

    int currentEdgeID = -1;
    const int activeEdgeID = activeEdgeID_.get();
    for(const VesselGraphEdge& edge : graph->getEdges()) {
        ++currentEdgeID;
        if(currentEdgeID == activeEdgeID) {
            continue;
        }
        const tgt::vec3 axis = edge.getNode2().pos_ - edge.getNode1().pos_;

        MatStack.pushMatrix();

        MatStack.translate(edge.getNode2().pos_);
        MatStack.multMatrix(generateMatrixFromQuat(tgt::generateQuaternionFromTo( tgt::normalize(axis), tgt::vec3(0,0,1))));
        MatStack.scale(tgt::vec3(edgeCrossSection, edgeCrossSection, -tgt::length(axis)));

        // set matrix stack uniforms
        shader->setUniform("modelViewMatrixStack_", MatStack.getModelViewMatrix());
        shader->setUniform("modelViewProjectionMatrixStack_", MatStack.getProjectionMatrix() * MatStack.getModelViewMatrix());
        tgt::mat4 viewInverse;
        MatStack.getModelViewMatrix().invert(viewInverse);
        shader->setUniform("modelViewMatrixInverseStack_", viewInverse);

        shader->setUniform("property_", currentEdgeTF.getValue(edge));

        edgeTriangleMesh_.render(GL_TRIANGLES);

        MatStack.popMatrix();
    }
    shader->deactivate();

    if(activeEdgeID != -1) {
        const VesselGraphEdge& edge = graph->getEdges()[activeEdgeID];

        tgt::Shader* shader = activeEdgeShader_.getShader();
        if(!shader) {
            LWARNING("No active edge shader");
            return;
        }
        glDisable(GL_DEPTH_TEST);
        shader->activate();

        const tgt::vec3 axis = edge.getNode2().pos_ - edge.getNode1().pos_;

        MatStack.pushMatrix();

        MatStack.translate(edge.getNode2().pos_);
        MatStack.multMatrix(generateMatrixFromQuat(tgt::generateQuaternionFromTo( tgt::normalize(axis), tgt::vec3(0,0,1))));
        MatStack.scale(tgt::vec3(edgeCrossSection, edgeCrossSection, -tgt::length(axis)));

        // set matrix stack uniforms
        shader->setUniform("modelViewMatrixStack_", MatStack.getModelViewMatrix());
        shader->setUniform("modelViewProjectionMatrixStack_", MatStack.getProjectionMatrix() * MatStack.getModelViewMatrix());
        tgt::mat4 viewInverse;
        MatStack.getModelViewMatrix().invert(viewInverse);
        shader->setUniform("modelViewMatrixInverseStack_", viewInverse);

        shader->setUniform("color_", activeEdgeColor_.get());

        setLightingUniforms(*shader);

        edgeTriangleMesh_.render(GL_TRIANGLES);

        MatStack.popMatrix();

        shader->deactivate();
        glEnable(GL_DEPTH_TEST);
    }

}

tgt::Bounds VesselGraphRenderer::getBoundingBox() const {
    const VesselGraph* graph = graphInport_.getData();
    if(!graph) {
        return tgt::Bounds();
    }
    return graph->getBounds();
}

void VesselGraphRenderer::adaptToNewInput() {
    const VesselGraph* graph = graphInport_.getData();
    if(!graph) {
        return;
    }

    activeEdgeID_.setMaxValue(static_cast<int>(graph->getEdges().size()) - 1);


    auto nodes = graph->getNodes();
    auto edges = graph->getEdges();

    nodeDegreeTF_.setData(nodes.begin(), nodes.end());

    // Only use valid edges for TFs (otherwise e.g. a mean maxRadius of -1 (symbolic)
    // will mess up the TF range and thus its precision)
    auto is_valid_edge = [] (const VesselGraphEdge& e) {
        return e.getMinRadiusAvg() >= 0;
    };

    for(auto& pair : edgeTFs_) {
        pair.second.setData(
                boost::make_filter_iterator(is_valid_edge, edges.begin(), edges.end()),
                boost::make_filter_iterator(is_valid_edge, edges.end(), edges.end())
                );
    }
}

void VesselGraphRenderer::showActiveEdgeTF() {
    for(auto& pair : edgeTFs_) {
        pair.second.getProperty().setVisibleFlag(pair.first == activeEdgeProperty_.getValue());
    }
}
} // namespace voreen
