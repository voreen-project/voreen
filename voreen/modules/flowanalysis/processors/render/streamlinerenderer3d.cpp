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

#include "streamlinerenderer3d.h"

#include "voreen/core/datastructures/volume/volumeatomic.h"
#include "voreen/core/datastructures/volume/volumeminmax.h"
#include "voreen/core/datastructures/callback/lambdacallback.h"

#include "../../utils/flowutils.h"

namespace voreen {

// TODO: Should use Shaders!
static const uint32_t GEOMETRY_TESSELATION = 8;

std::string StreamlineRenderer3D::loggerCat_("flowanalysis.StreamlineRenderer3D");

StreamlineRenderer3D::StreamlineRenderer3D()
    : RenderProcessor()
    //ports
    , streamlineInport_(Port::INPORT, "streamlineInport", "Streamline Input")
    , imgOutport_(Port::OUTPORT, "image.streamlines", "Streamline Image Output", true, Processor::INVALID_RESULT, RenderPort::RENDERSIZE_RECEIVER)
    //properties
        //style
    , streamlineStyle_("streamlineStyle", "Streamline Style:")
    , lineWidth_("lineWidth", "Line Width:", 1.0f, 1.0f, 10.0f)
        //color
    , color_("colorProp", "Color:")
    , transferFunction_("tfProp", "Color Map:")
    , colorRotationMatrix_("colorRotationMatrix", "To be linked with FlowDirectionOverlay", tgt::mat4::identity,
                           tgt::mat4(-1.1f), tgt::mat4(1.1f), Processor::INVALID_RESULT, NumericProperty<tgt::mat4>::STATIC, Property::LOD_DEBUG)
    , rotateAroundX_("rotationaroundx", "x-Axis Rotation (degrees)")
    , rotateAroundY_("rotationaroundy", "y-Axis Rotation (degrees)")
    , rotateAroundZ_("rotationaroundz", "z-Axis Rotation (degrees)")
    , enableMaximumIntensityProjection_("maximumIntensityProjection", "Enable Maximum Intensity Projection (MIP)", false)
    , timeWindowStart_("timeWindowStart", "Time Window Start", 0.0f, 0.0f, 1000.0f, Processor::INVALID_RESULT, FloatProperty::DYNAMIC)
    , timeWindowSize_("timeWindowSize", "Time Window Size", 0.0f, 0.0f, 1000.0f, Processor::INVALID_RESULT, FloatProperty::DYNAMIC)
        //lighting
    , enableLighting_("enableLighting", "Enable Lighting", false)
    , lightPosition_("lightPosition", "Light Source Position", tgt::vec4(2.3f, 1.5f, 1.5f, 1.f), tgt::vec4(-10000), tgt::vec4(10000))
    , lightAmbient_("lightAmbient", "Ambient Light", tgt::Color(0.4f, 0.4f, 0.4f, 1.f))
    , lightDiffuse_("lightDiffuse", "Diffuse Light", tgt::Color(0.8f, 0.8f, 0.8f, 1.f))
    , lightSpecular_("lightSpecular", "Specular Light", tgt::Color(0.6f, 0.6f, 0.6f, 1.f))
    , materialShininess_("materialShininess", "Shininess", 60.f, 0.1f, 128.f)
        //must haves
    , streamlineShader_("streamlineShaderProp", "Shader:", "streamlinerenderer3d.frag", "streamlinerenderer3d.vert", ""/*"streamlinerenderer3d.geom"*/, Processor::INVALID_PROGRAM, Property::LOD_DEBUG)
    , requiresRecompileShader_(true)
    , camera_("camera", "Camera", tgt::Camera(tgt::vec3(0.f, 0.f, 3.5f), tgt::vec3(0.f, 0.f, 0.f), tgt::vec3(0.f, 1.f, 0.f)))
    , cameraHandler_(nullptr)
    , requiresRebuild_(false)
{
    //ports
    addPort(streamlineInport_);
        ON_CHANGE(streamlineInport_, StreamlineRenderer3D, onStreamlineDataChange);
    addPort(imgOutport_);
    //properties
        // style
    addProperty(streamlineStyle_);
        streamlineStyle_.addOption("lines", "Lines", STYLE_LINES);
        streamlineStyle_.addOption("tubes", "Tubes", STYLE_TUBES);
        streamlineStyle_.onChange(MemberFunctionCallback<StreamlineRenderer3D>(this, &StreamlineRenderer3D::onStyleChange));
    addProperty(lineWidth_);
        ON_CHANGE_LAMBDA(lineWidth_, [this] {
            requiresRebuild_ |= streamlineStyle_.getValue() != STYLE_LINES;
        });
        // color
    addProperty(color_);
        color_.addOption("velocity" , "Velocity" , COLOR_VELOCITY);
        color_.addOption("direction" , "Direction" , COLOR_DIRECTION);
        color_.onChange(MemberFunctionCallback<StreamlineRenderer3D>(this, &StreamlineRenderer3D::onColorChange));
    addProperty(transferFunction_);
    addProperty(rotateAroundX_);
        rotateAroundX_.addOption("0", "0", 0.f);
        rotateAroundX_.addOption("90", "90", 90.f);
        rotateAroundX_.addOption("180", "180", 180.f);
        rotateAroundX_.addOption("270", "270", 270.f);
        rotateAroundX_.setVisibleFlag(false);
        rotateAroundX_.onChange(MemberFunctionCallback<StreamlineRenderer3D>(this, &StreamlineRenderer3D::computeDirectionColorRotationMatrix));
    addProperty(rotateAroundY_);
        rotateAroundY_.addOption("0", "0", 0.f);
        rotateAroundY_.addOption("90", "90", 90.f);
        rotateAroundY_.addOption("180", "180", 180.f);
        rotateAroundY_.addOption("270", "270", 270.f);
        rotateAroundY_.setVisibleFlag(false);
        rotateAroundY_.onChange(MemberFunctionCallback<StreamlineRenderer3D>(this, &StreamlineRenderer3D::computeDirectionColorRotationMatrix));
    addProperty(rotateAroundZ_);
        rotateAroundZ_.addOption("0", "0", 0.f);
        rotateAroundZ_.addOption("90", "90", 90.f);
        rotateAroundZ_.addOption("180", "180", 180.f);
        rotateAroundZ_.addOption("270", "270", 270.f);
        rotateAroundZ_.setVisibleFlag(false);
        rotateAroundZ_.onChange(MemberFunctionCallback<StreamlineRenderer3D>(this, &StreamlineRenderer3D::computeDirectionColorRotationMatrix));
    addProperty(colorRotationMatrix_);
        colorRotationMatrix_.setVisibleFlag(false);
        colorRotationMatrix_.setReadOnlyFlag(true);
    addProperty(enableMaximumIntensityProjection_);
    addProperty(timeWindowStart_);
        timeWindowStart_.setNumDecimals(4);
        timeWindowStart_.setReadOnlyFlag(true);
    addProperty(timeWindowSize_);
        timeWindowSize_.setNumDecimals(4);
        timeWindowSize_.setReadOnlyFlag(true);
        //lighting
    addProperty(enableLighting_);
    ON_CHANGE_LAMBDA(enableLighting_, [this] {
        bool lighting = enableLighting_.get();
        setPropertyGroupVisible("lighting", lighting);
    });
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
        //must have
    addProperty(streamlineShader_);
    addProperty(camera_);
        cameraHandler_ = new CameraInteractionHandler("cameraHandler", "Camera Handler", &camera_);
    addInteractionHandler(cameraHandler_);
}

StreamlineRenderer3D::~StreamlineRenderer3D() {
    delete cameraHandler_;
}

void StreamlineRenderer3D::initialize() {
    RenderProcessor::initialize();
    requiresRecompileShader_ = true;
}

void StreamlineRenderer3D::deinitialize() {
    meshes_.clear();
    RenderProcessor::deinitialize();
}

void StreamlineRenderer3D::beforeProcess() {
    RenderProcessor::beforeProcess();

    if (requiresRecompileShader_)
        compile();

    if (requiresRebuild_)
        rebuild();
}

void StreamlineRenderer3D::process() {

    //activate output
    imgOutport_.activateTarget();
    imgOutport_.clearTarget();

    if(tgt::Shader* shader = streamlineShader_.getShader()) {
        shader->activate();

        const tgt::Camera& camera = camera_.get();

        // set transformation uniforms
        setGlobalShaderParameters(shader, &camera, imgOutport_.getSize());
        shader->setUniform("velocityTransformMatrix_",streamlineInport_.getData()->getVelocityTransformMatrix());
        shader->setUniform("timeWindowStart_", timeWindowStart_.get());
        shader->setUniform("timeWindowSize_", timeWindowSize_.get());

        shader->setUniform("enableLighting_", enableLighting_.get());
        if(enableLighting_.get()) {
            shader->setUniform("cameraPosition_", camera.getPosition());
            shader->setUniform("lightSource_.position_", lightPosition_.get().xyz());
            shader->setUniform("lightSource_.ambientColor_", lightAmbient_.get().xyz());
            shader->setUniform("lightSource_.diffuseColor_", lightDiffuse_.get().xyz());
            shader->setUniform("lightSource_.specularColor_", lightSpecular_.get().xyz());
            //prog->setUniform("lightSource_.attenuation_", tgt::vec3(1.f, 0.f, 0.f));
            shader->setUniform("shininess_", materialShininess_.get());
        }

        switch(color_.getValue()) {
        case COLOR_VELOCITY:
        {
            tgt::TextureUnit transFuncUnit;
            transFuncUnit.activate();
            transferFunction_.get()->getTexture()->bind();
            transferFunction_.get()->setUniform(shader, "transFuncParam_", "transFuncTex_", transFuncUnit.getUnitNumber());

            glEnable(GL_BLEND);

            // In case of MIP rendering, we change the blend equation accordingly.
            if(enableMaximumIntensityProjection_.get()) {
                glBlendFuncSeparate(GL_ONE, GL_ONE, GL_ONE, GL_ONE);
                glBlendEquation(GL_MAX);
                glDisable(GL_DEPTH_TEST);
            }
            else {
                glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
                if(transferFunction_.get()->getAlphaMode() == TransFuncBase::TF_USE_ALPHA) {
                    glDisable(GL_DEPTH_TEST);
                }
            }
            break;
        }
        case COLOR_DIRECTION:
            shader->setUniform("colorRotationMatrix_", colorRotationMatrix_.get(), true);
            break;
        default:
            tgtAssert(false,"Unknown Color Coding");
            LERROR("Unknown Color Coding");
        }

        glLineWidth(lineWidth_.get());

        for(size_t i = 0; i < meshes_.size(); i++) {
            meshes_[i]->render();
        }

        // Restore state.
        glLineWidth(1.0f);
        glBlendEquation(GL_FUNC_ADD);
        glBlendFunc(GL_ONE, GL_ZERO);
        glDisable(GL_BLEND);
        glEnable(GL_DEPTH_TEST);
        shader->deactivate();
        LGL_ERROR;
    }

    imgOutport_.deactivateTarget();
}

//--------------------------------------------------------------------------------------
//      Callbacks
//--------------------------------------------------------------------------------------
void StreamlineRenderer3D::onStreamlineDataChange() {

    const StreamlineListBase* streamlines = streamlineInport_.getData();
    if (!streamlines)
        return;

    // Fit camera.
    camera_.adaptInteractionToScene(streamlines->getOriginalWorldBounds(), 0.0f, true);

    // Update transfer function.
    float* data = new float[2];
    data[0] = streamlines->getMinMagnitude();
    data[1] = streamlines->getMaxMagnitude();
    VolumeRAM_Float* representation = new VolumeRAM_Float(data, tgt::svec3(2,1,1));
    tfVolume_.reset(new Volume(representation, tgt::vec3::one, tgt::vec3::zero));
    tfVolume_->addDerivedData(new VolumeMinMax(data[0], data[1], data[0], data[1])); // To save time by not triggering a background thread.
    transferFunction_.setVolume(tfVolume_.get());

    // Update time window.
    bool timeWindowEmpty = streamlines->getTemporalRange().x == streamlines->getTemporalRange().y;
    timeWindowStart_.setReadOnlyFlag(timeWindowEmpty);
    timeWindowSize_.setReadOnlyFlag(timeWindowEmpty);
    std::deque<Option<StreamlineStyle>> styleOptions;
    styleOptions.emplace_back(Option<StreamlineStyle>("lines", "Lines", STYLE_LINES));
    if(timeWindowEmpty) {
        styleOptions.emplace_back(Option<StreamlineStyle>("tubes", "Tubes", STYLE_TUBES));
        // styleOptions.emplace_back(Option<StreamlineStyle>("arrows", "Arrows", STYLE_ARROWS)); // TODO: should use shaders!
    }
    streamlineStyle_.setOptions(styleOptions);
    streamlineStyle_.invalidate();

    bool timeWindowChanged = streamlines->getTemporalRange().x != timeWindowStart_.getMinValue()
            || (streamlines->getTemporalRange().y - streamlines->getTemporalRange().x) != timeWindowSize_.getMaxValue();
    if(timeWindowChanged) {
        timeWindowStart_.setMinValue(streamlines->getTemporalRange().x);
        timeWindowStart_.setMaxValue(streamlines->getTemporalRange().y);
        timeWindowStart_.set(timeWindowStart_.getMinValue()); // default: Start at front.

        timeWindowSize_.setMinValue(0);
        timeWindowSize_.setMaxValue(streamlines->getTemporalRange().y - streamlines->getTemporalRange().x);
        timeWindowSize_.set(timeWindowSize_.getMaxValue()); // default: End at back.
    }

    // Force rebuild
    requiresRebuild_ = true;
}

void StreamlineRenderer3D::onStyleChange() {
    requiresRebuild_ = true;
    lineWidth_.setTracking(streamlineStyle_.getValue() == STYLE_LINES);
    /*
    //update shaders
    switch(streamlineStyleProp_.getValue()) {
    case STYLE_LINES:
    case STYLE_TUBES:
    case STYLE_ARROWS:
        requiresRecompileShader_ = true;
        break;

    default:
        LERROR("Unsupported style");
    }
     */
}

void StreamlineRenderer3D::onColorChange() {
    //update visibility
    switch(color_.getValue()) {
    case COLOR_VELOCITY:
        transferFunction_.setVisibleFlag(true);
        rotateAroundX_.setVisibleFlag(false);
        rotateAroundY_.setVisibleFlag(false);
        rotateAroundZ_.setVisibleFlag(false);
        colorRotationMatrix_.setVisibleFlag(false);
        enableMaximumIntensityProjection_.setVisibleFlag(true);
        break;
    case COLOR_DIRECTION:
        transferFunction_.setVisibleFlag(false);
        rotateAroundX_.setVisibleFlag(true);
        rotateAroundY_.setVisibleFlag(true);
        rotateAroundZ_.setVisibleFlag(true);
        colorRotationMatrix_.setVisibleFlag(true);
        enableMaximumIntensityProjection_.setVisibleFlag(false);
        break;
    default:
        LERROR("Unsupported color coding");
    }
    // shader header must be updated
    requiresRecompileShader_ = true;
}

void StreamlineRenderer3D::computeDirectionColorRotationMatrix() {
    // start with identity
    tgt::mat4 rotation = tgt::mat4::createIdentity();

    // apply the three rotations if necessary
    if (rotateAroundX_.getValue() != 0.f)
        rotation = tgt::mat4::createRotationXDegree(rotateAroundX_.getValue()) * rotation;

    if (rotateAroundY_.getValue() != 0.f)
        rotation = tgt::mat4::createRotationYDegree(rotateAroundY_.getValue()) * rotation;

    if (rotateAroundZ_.getValue() != 0.f)
        rotation = tgt::mat4::createRotationZDegree(rotateAroundZ_.getValue()) * rotation;

    // set the transformation
    colorRotationMatrix_.set(rotation);
}


//--------------------------------------------------------------------------------------
//      Helpers
//--------------------------------------------------------------------------------------
std::string StreamlineRenderer3D::generateHeader(const tgt::GpuCapabilities::GlVersion*) {
    //generate basic header
    std::string header = RenderProcessor::generateHeader(&tgt::GpuCapabilities::GlVersion::SHADER_VERSION_330);
    //add define for fragment shader
    switch(color_.getValue()) {
    case COLOR_VELOCITY:
        header += "#define COLOR_VELOCITY\n";
        header += transferFunction_.get()->getShaderDefines();
        break;
    case COLOR_DIRECTION:
        header += "#define COLOR_DIRECTION\n";
        break;
    default:
        tgtAssert(false,"Should not get here");
        LERROR("Unknown Color Coding");
    }
    return header;
}

void StreamlineRenderer3D::compile() {
    streamlineShader_.setHeader(generateHeader());
    streamlineShader_.rebuild();
    requiresRecompileShader_ = false;
}

void StreamlineRenderer3D::rebuild() {

    // clear old mesh data
    meshes_.clear();

    // create new meshes according to selected option
    switch (streamlineStyle_.getValue()) {
    case STYLE_LINES:
        meshes_.push_back(std::unique_ptr<GlMeshGeometryBase>(createLineGeometry(streamlineInport_.getData()->getStreamlines())));
        break;
    case STYLE_TUBES:
        for(const Streamline& streamline : streamlineInport_.getData()->getStreamlines()) {
            meshes_.push_back(std::unique_ptr<GlMeshGeometryBase>(createTubeGeometry(streamline)));
        }
        break;
    case STYLE_ARROWS:
        for(const Streamline& streamline : streamlineInport_.getData()->getStreamlines()) {
            meshes_.push_back(std::unique_ptr<GlMeshGeometryBase>(createArrowGeometry(streamline)));
        }
        break;
    default:
        LERROR("Unsupported option");
    }

    requiresRebuild_ = false;
}

GlMeshGeometryBase* StreamlineRenderer3D::createLineGeometry(const std::vector<Streamline>& streamlines) {

    GlMeshGeometryUInt32Color* mesh = new GlMeshGeometryUInt32Color();
    mesh->setPrimitiveType(GL_LINE_STRIP);
    mesh->enablePrimitiveRestart();

    for (const Streamline& streamline : streamlines) {
        for (size_t i = 0; i < streamline.getNumElements(); i++) {
            const Streamline::StreamlineElement& element = streamline.getElementAt(i);
            mesh->addIndex(static_cast<uint32_t>(mesh->getVertices().size()));
            mesh->addVertex(VertexColor(element.position_, tgt::vec4(element.velocity_, element.time_)));
            // TODO: we abuse color.a as time and radius is currently unused.
        }
        mesh->addIndex(mesh->getPrimitiveRestartIndex());
    }

    return mesh;
}

GlMeshGeometryBase* StreamlineRenderer3D::createTubeGeometry(const Streamline& streamline) const {

    GlMeshGeometryUInt32Color* mesh = new GlMeshGeometryUInt32Color();
    mesh->setPrimitiveType(GL_TRIANGLES);

    const uint32_t tesselation = GEOMETRY_TESSELATION;
    const float angleStep = (tgt::PIf * 2.0f) / tesselation;

    for (size_t k = 0; k < streamline.getNumElements(); k++) {

        const Streamline::StreamlineElement& element = streamline.getElementAt(k);
        tgt::mat4 transformation = createTransformationMatrix(element.position_, element.velocity_);

        uint32_t offset = static_cast<uint32_t>(k) * tesselation;

        float radius = std::max(element.radius_, lineWidth_.get());
        tgt::vec3 color = element.velocity_;
        tgt::vec4 v(0.0f, 0.0f, 0.0f, 1.0f);
        for (uint32_t i = 0; i < tesselation; i++) {

            // Calculate Vertex and add it to the mesh.
            float angle = angleStep * i;
            v.x = radius * cosf(angle);
            v.y = radius * sinf(angle);
            mesh->addVertex(VertexColor((transformation * v).xyz(), tgt::vec4(color, element.time_)));

            // Calculate indices and add them to the mesh.
            if (k < streamline.getNumElements() - 1) {

                unsigned int j = (i + 1) % tesselation;

                mesh->addIndex(offset + i);
                mesh->addIndex(offset + j);
                mesh->addIndex(offset + i + tesselation);
                mesh->addIndex(offset + j);
                mesh->addIndex(offset + j + tesselation);
                mesh->addIndex(offset + i + tesselation);
            }
        }
    }

    return mesh;
}

GlMeshGeometryBase* StreamlineRenderer3D::createArrowGeometry(const Streamline& streamline) const {

    GlMeshGeometryUInt32Color* mesh = new GlMeshGeometryUInt32Color();
    mesh->setPrimitiveType(GL_TRIANGLES);

    // Define some constants for easier access.
    const uint32_t tesselation = GEOMETRY_TESSELATION;
    const float angleStep = (tgt::PIf * 2.0f) / tesselation;

    float radius = std::max(streamline.getFirstElement().radius_, lineWidth_.get()); // FIXME: radius is only a rough approximation.
    float length = streamline.getPhysicalLength();

    // Calculate the number of arrows the bundle will be split into.
    const float desiredLengthOfArrow = std::min(radius * 4.0f, length);
    const size_t samples = std::max<size_t>(2, length / desiredLengthOfArrow);

    // Resample the bundle.
    Streamline centroid = streamline.resample(samples);

    // Specify the dimensions of each arrow-submesh.
    const float arrowRatio = 0.9f;
    const float lengthOfArrow = (length / static_cast<float>(samples - 1)) * arrowRatio;
    const float lengthOfCylinder = lengthOfArrow * 0.7f;

    for (size_t k = 0; k < centroid.getNumElements() - 1; k++) {

        const float radiusOfCylinder = centroid.getElementAt(k).radius_ * 0.7f;
        const float radiusOfCone     = centroid.getElementAt(k).radius_ * 1.1f;

        const tgt::vec3& p0 = centroid.getElementAt(k).position_;
        const tgt::vec3& p1 = centroid.getElementAt(k+1).position_;

        tgt::mat4 transformation = createTransformationMatrix(p0, p1 - p0);

        uint32_t offset = static_cast<uint32_t>(mesh->getVertices().size());

        tgt::vec3 color = centroid.getElementAt(k).velocity_;
        tgt::vec4 v(0.0f, 0.0f, 0.0f, 1.0f);

        // Cylinder
        mesh->addVertex(VertexColor(centroid.getElementAt(k).position_, tgt::vec4(color, 1.0f)));
        for (uint32_t i = 0; i < tesselation; i++) {

            // Calculate Vertex and add it to the mesh.
            float angle = angleStep * i;
            v.x = radiusOfCylinder * cosf(angle);
            v.y = radiusOfCylinder * sinf(angle);
            v.z = 0.0f;
            mesh->addVertex(VertexColor((transformation * v).xyz(), tgt::vec4(color, 1.0f)));

            v.z = lengthOfCylinder;
            mesh->addVertex(VertexColor((transformation * v).xyz(), tgt::vec4(color, 1.0f)));

            // Calculate indices and add them to the mesh.

            // Bottom of the Cylinder
            mesh->addIndex(offset + (i * 2 + 3) % (2 * tesselation));
            mesh->addIndex(offset + i * 2 + 1);
            mesh->addIndex(offset);

            // Hull of the Cylinder
            if(i < tesselation - 1) {
                mesh->addIndex(offset + i * 2 + 1);
                mesh->addIndex(offset + i * 2 + 3);
                mesh->addIndex(offset + i * 2 + 2);
                mesh->addIndex(offset + i * 2 + 3);
                mesh->addIndex(offset + i * 2 + 4);
                mesh->addIndex(offset + i * 2 + 2);
            } else { // Handle last quad separately
                mesh->addIndex(offset + i * 2 + 1);
                mesh->addIndex(offset + 1);
                mesh->addIndex(offset + i * 2 + 2);
                mesh->addIndex(offset + 1);
                mesh->addIndex(offset + 2);
                mesh->addIndex(offset + i * 2 + 2);
            }
        }

        // Cone
        offset = static_cast<uint32_t>(mesh->getVertices().size());
        v = tgt::vec4(0.0f, 0.0f, lengthOfCylinder, 1.0f);
        mesh->addVertex(VertexColor((transformation * v).xyz(), tgt::vec4(color, 1.0f)));
        for (uint32_t i = 0; i < tesselation; i++) {

            // Calculate Vertex and add it to the mesh.
            float angle = angleStep * i;
            v.x = radiusOfCone * cosf(angle);
            v.y = radiusOfCone * sinf(angle);
            mesh->addVertex(VertexColor((transformation * v).xyz(), tgt::vec4(color, 1.0f)));

            // Calculate indices and add them to the mesh.

            // Bottom of the Cone
            mesh->addIndex(offset + (i + 1) % tesselation + 1);
            mesh->addIndex(offset + i + 1);
            mesh->addIndex(offset);

            // Top of the Cone
            mesh->addIndex(offset + i + 1);
            mesh->addIndex(offset + (i + 1) % tesselation + 1);
            mesh->addIndex(offset + tesselation + 1);
        }

        v = tgt::vec4(0.0f, 0.0f, lengthOfArrow, 1.0f);
        mesh->addVertex(VertexColor((transformation * v).xyz(), tgt::vec4(color, 1.0f)));
    }

    return mesh;
}

}   // namespace
