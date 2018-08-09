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

#include "streamlinerenderer3d.h"

#include "voreen/core/datastructures/volume/volumeatomic.h"
#include "voreen/core/datastructures/volume/volumeminmax.h"
#include "voreen/core/datastructures/callback/lambdacallback.h"

namespace voreen {

    // TODO: Could be adjustable via Property ( debug ).
    static const uint32_t GEOMETRY_TESSELATION = 16;
    static const float CENTROID_LINE_WIDTH = 3.0f;

StreamlineRenderer3D::StreamlineRenderer3D()
    : RenderProcessor()
    //ports
    , streamlineInport_(Port::INPORT, "streamlineInport", "Streamline Input")
    , imgOutport_(Port::OUTPORT, "image.streamlines", "Streamline Image Output", true, Processor::INVALID_RESULT, RenderPort::RENDERSIZE_RECEIVER)
    //properties
        //rendering option
    , renderingOptionProp_("renderingOption", "Render:")
        //style
    , streamlineStyleProp_("streamlineStyle", "Streamline Style:")
    , streamlineBundleStyleProp_("streamlineBundleStyle", "Streamline Bundle Style:")
        //color
    , colorProp_("colorProp","Color:")
    , tfProp_("tfProp","Color Map:")
    , tfVolume_(0)
    , colorRotationMatrix_("colorRotationMatrix", "To be linked with FlowDirectionOverlay", tgt::mat4::identity,
                           tgt::mat4(-1.1f), tgt::mat4(1.1f), Processor::INVALID_RESULT, NumericProperty<tgt::mat4>::STATIC, Property::LOD_DEBUG)
    , rotateAroundX_("rotationaroundx", "x-Axis Rotation (degrees)")
    , rotateAroundY_("rotationaroundy", "y-Axis Rotation (degrees)")
    , rotateAroundZ_("rotationaroundz", "z-Axis Rotation (degrees)")
        //must haves
    , streamlineShaderProp_("streamlineShaderProp", "Shader:","streamlinerenderer3d.frag","streamlinerenderer3d.vert",""/*"streamlinerenderer3d.geom"*/,Processor::INVALID_PROGRAM,Property::LOD_DEBUG)
    ,       requiresRecompileShader_(true)
    , cameraProp_("camera", "Camera: ", tgt::Camera(),true,true,500.f,Processor::INVALID_RESULT,Property::LOD_DEBUG)
    ,       cameraHandler_(0)
    , bundleMeshStartIndex_(0)
    , requiresRebuild_(false)
{
    //ports
    addPort(streamlineInport_);
        streamlineInport_.onNewData(MemberFunctionCallback<StreamlineRenderer3D>(this, &StreamlineRenderer3D::onStreamlineDataChange));
    addPort(imgOutport_);
    //properties
        // rendering option
    addProperty(renderingOptionProp_);
        renderingOptionProp_.addOption("streamlines", "Streamlines (No Bundles available)", OPTION_STREAMLINES);
        renderingOptionProp_.addOption("streamlinebundles noise off", "Streamline Bundles - Noise off", OPTION_STREAMLINEBUNDLES_NOISE_OFF);
        renderingOptionProp_.addOption("streamlinebundles noise on", "Streamline Bundles - Noise on", OPTION_STREAMLINEBUNDLES_NOISE_ON);
        renderingOptionProp_.setReadOnlyFlag(true);
        renderingOptionProp_.onChange(LambdaFunctionCallback([this] { requiresRebuild_ = true; }));
    // style
    addProperty(streamlineStyleProp_);
        streamlineStyleProp_.addOption("lines", "Lines", STYLE_LINES);
        streamlineStyleProp_.addOption("arrows", "Arrows", STYLE_ARROWS);
        streamlineStyleProp_.onChange(MemberFunctionCallback<StreamlineRenderer3D>(this, &StreamlineRenderer3D::onStyleChange));
    addProperty(streamlineBundleStyleProp_);
        streamlineBundleStyleProp_.addOption("lines", "Lines", BUNDLE_STYLE_LINES);
        streamlineBundleStyleProp_.addOption("tubes", "Tubes", BUNDLE_STYLE_TUBES);
        streamlineBundleStyleProp_.addOption("arrows", "Arrows", BUNDLE_STYLE_ARROWS);
        streamlineBundleStyleProp_.onChange(LambdaFunctionCallback([this] { requiresRebuild_ = true; })); // force mesh rebuild
        // color
    addProperty(colorProp_);
        colorProp_.addOption("velocity" , "Velocity" , COLOR_VELOCITY);
        colorProp_.addOption("direction" , "Direction" , COLOR_DIRECTION);
        colorProp_.onChange(MemberFunctionCallback<StreamlineRenderer3D>(this, &StreamlineRenderer3D::onColorChange));
    addProperty(tfProp_);
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

        //must have
    addProperty(streamlineShaderProp_);
    addProperty(cameraProp_);
        cameraHandler_ = new CameraInteractionHandler("cameraHandler", "Camera Handler", &cameraProp_);
    addInteractionHandler(cameraHandler_);
}

StreamlineRenderer3D::~StreamlineRenderer3D() {
    delete cameraHandler_;
    delete tfVolume_;
}

void StreamlineRenderer3D::initialize() {
    RenderProcessor::initialize();

    //init shader
    if(streamlineStyleProp_.getValue() == STYLE_LINES)
        requiresRecompileShader_ = true;
    else
        LERROR("Unsupported style");
}

void StreamlineRenderer3D::deinitialize() {

    // clear mesh data
    for(GlMeshGeometryUInt32Color* mesh : meshes_)
        delete mesh;
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

    if(tgt::Shader* shader = streamlineShaderProp_.getShader()) {
        shader->activate();

        // set transformation uniforms
        setGlobalShaderParameters(shader,&cameraProp_.get(),imgOutport_.getSize());
        shader->setUniform("voxelToWorldMatrix_", streamlineInport_.getData()->getVoxelToWorldMatrix());
        shader->setUniform("velocityTransformMatrix_",streamlineInport_.getData()->getVelocityTransformMatrix());

        //set color uniforms
        switch(colorProp_.getValue()) {
        case COLOR_VELOCITY:
        {
            tgt::TextureUnit transFuncUnit;
            transFuncUnit.activate();
            tfProp_.get()->getTexture()->bind();
            tfProp_.get()->setUniform(shader, "transFuncParam_", "transFuncTex_", transFuncUnit.getUnitNumber());
            break;
        }
        case COLOR_DIRECTION:
            shader->setUniform("colorRotationMatrix_", colorRotationMatrix_.get(), true);
            break;
        default:
            tgtAssert(false,"Should not get here!");
            LERROR("Unknown Color Coding");
        }

        // Render every Mesh being created before.
        for(size_t i = 0; i < meshes_.size(); i++) {
            glLineWidth((i >= bundleMeshStartIndex_) ? CENTROID_LINE_WIDTH : 1.0f);
            meshes_[i]->render();
        }

        // Restore state.
        glLineWidth(1.0f);
        shader->deactivate();
        LGL_ERROR;
    }

    imgOutport_.deactivateTarget();
}

//--------------------------------------------------------------------------------------
//      Callbacks
//--------------------------------------------------------------------------------------
void StreamlineRenderer3D::onStreamlineDataChange() {

    //update camera
    cameraProp_.adaptInteractionToScene(streamlineInport_.getData()->getOriginalWorldBounds());

    //update tf
    delete tfVolume_;
    float* data = new float[2];
    data[0] = streamlineInport_.getData()->getMinMagnitude();
    data[1] = streamlineInport_.getData()->getMaxMagnitude();
    VolumeRAM_Float* rep = new VolumeRAM_Float(data, tgt::svec3(2,1,1));
    tfVolume_ = new Volume(rep, tgt::vec3::one, tgt::vec3::zero);
    tfVolume_->addDerivedData(new VolumeMinMax(data[0], data[1], data[0], data[1])); //to save time by not triggering an background thread
    tfProp_.setVolume(tfVolume_, 0);

    // Update Rendering Option selection.
    if (streamlineInport_.getData()->getStreamlineBundles().empty()) {
        renderingOptionProp_.setReadOnlyFlag(true);
        renderingOptionProp_.setOptionDescription("streamlines", "Streamlines (No Bundles available)");
        renderingOptionProp_.selectByValue(OPTION_STREAMLINES);
    }
    else {
        renderingOptionProp_.setReadOnlyFlag(false);
        renderingOptionProp_.setOptionDescription("streamlines", "Streamlines");
    }

    // Force rebuild
    requiresRebuild_ = true;
}

void StreamlineRenderer3D::onStyleChange() {
    //update shaders
    switch(streamlineStyleProp_.getValue()) {
    case STYLE_LINES:
        requiresRecompileShader_ = true;
        break;
    default:
        LERROR("Unsupported style");
    }
}

void StreamlineRenderer3D::onColorChange() {
    //update visibility
    switch(colorProp_.getValue()) {
    case COLOR_VELOCITY:
        tfProp_.setVisibleFlag(true);
        rotateAroundX_.setVisibleFlag(false);
        rotateAroundY_.setVisibleFlag(false);
        rotateAroundZ_.setVisibleFlag(false);
        colorRotationMatrix_.setVisibleFlag(false);
        break;
    case COLOR_DIRECTION:
        tfProp_.setVisibleFlag(false);
        rotateAroundX_.setVisibleFlag(true);
        rotateAroundY_.setVisibleFlag(true);
        rotateAroundZ_.setVisibleFlag(true);
        colorRotationMatrix_.setVisibleFlag(true);
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
    tgt::GpuCapabilities::GlVersion* version = new tgt::GpuCapabilities::GlVersion(3,3,0);
    std::string header = RenderProcessor::generateHeader(version);
    delete version;
    //add define for fragment shader
    switch(colorProp_.getValue()) {
    case COLOR_VELOCITY:
        header += "#define COLOR_VELOCITY\n";
        header += tfProp_.get()->getShaderDefines();
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
    streamlineShaderProp_.setHeader(generateHeader());
    streamlineShaderProp_.rebuild();
    requiresRecompileShader_ = false;
}

void StreamlineRenderer3D::rebuild() {

    // clear old mesh data
    for (GlMeshGeometryUInt32Color* mesh : meshes_)
        delete mesh;
    meshes_.clear();

    // create new meshes according to selected option
    switch (renderingOptionProp_.getValue()) {
    case OPTION_STREAMLINES:
        buildStreamlineData(streamlineInport_.getData()->getStreamlines());
        bundleMeshStartIndex_ = 1;
        break;
    case OPTION_STREAMLINEBUNDLES_NOISE_OFF:
        buildStreamlineBundleData(streamlineInport_.getData()->getStreamlineBundles());
        bundleMeshStartIndex_ = 0;
        break;
    case OPTION_STREAMLINEBUNDLES_NOISE_ON:
    {
        std::vector<Streamline> noise;
        noise.reserve(streamlineInport_.getData()->getStreamlineNoise().size());
        for (size_t i : streamlineInport_.getData()->getStreamlineNoise())
            noise.push_back(streamlineInport_.getData()->getStreamlines()[i]);

        buildStreamlineData(noise);
        bundleMeshStartIndex_ = meshes_.size();
        buildStreamlineBundleData(streamlineInport_.getData()->getStreamlineBundles());

        break;
    }
    default:
        LERROR("Unsupported option");
    }

    requiresRebuild_ = false;
}

void StreamlineRenderer3D::buildStreamlineData(const std::vector<Streamline>& streamlines) {

    if(streamlines.empty())
        return;

    GlMeshGeometryUInt32Color* mesh = new GlMeshGeometryUInt32Color();
    mesh->setPrimitiveType(GL_LINE_STRIP);
    mesh->enablePrimitiveRestart();

    for (const Streamline& streamline : streamlines) {
        for (size_t l = 0; l < streamline.getNumElements(); l++) {
            mesh->addIndex(static_cast<uint32_t>(mesh->getVertices().size()));
            mesh->addVertex(VertexColor(streamline.getElementAt(l).position_, tgt::vec4(streamline.getElementAt(l).velocity_, 1.f)));
        }
        mesh->addIndex(mesh->getPrimitiveRestartIndex());
    }

    meshes_.push_back(mesh);
}

void StreamlineRenderer3D::buildStreamlineBundleData(const std::vector<StreamlineBundle>& bundles) {

    if(bundles.empty())
        return;

    switch(streamlineBundleStyleProp_.getValue()) {
    case BUNDLE_STYLE_LINES:
    {
        std::vector<Streamline> streamlines;
        streamlines.reserve(bundles.size());
        for (const StreamlineBundle& bundle : bundles) {
            streamlines.push_back(bundle.getCentroid());
        }

        buildStreamlineData(streamlines);

        break;
    }
    case BUNDLE_STYLE_TUBES:

        // Create one mesh for each bundle
        for (const StreamlineBundle& bundle : bundles) {
            meshes_.push_back(createStreamTube(bundle));
        }

        break;

    case BUNDLE_STYLE_ARROWS:

        // Create one mesh for each bundle
        for (const StreamlineBundle& bundle : bundles) {
            meshes_.push_back(createArrowPath(bundle));
        }

        break;

    default:
        LERROR("Unsupported style");
    }
}

GlMeshGeometryUInt32Color* StreamlineRenderer3D::createStreamTube(const StreamlineBundle& bundle) const {

    GlMeshGeometryUInt32Color* mesh = new GlMeshGeometryUInt32Color();
    mesh->setPrimitiveType(GL_TRIANGLES);

    const uint32_t tesselation = GEOMETRY_TESSELATION;
    const float angleStep = (tgt::PIf * 2.0f) / tesselation;
    const float radius = bundle.getRadius();

    // Return an empty mesh in case of a zero radius.
    if(radius == 0.0f)
        return mesh;

    Streamline centroid = bundle.getCentroid();

    for (size_t k = 0; k < centroid.getNumElements(); k++) {

        const Streamline::StreamlineElement& element = centroid.getElementAt(k);
        tgt::mat4 transformation = createTransformationMatrix(element.position_, element.velocity_);

        uint32_t offset = static_cast<uint32_t>(k) * tesselation;

        tgt::vec3 color = element.velocity_;
        tgt::vec4 v(0.0f, 0.0f, 0.0f, 1.0f);
        for (unsigned int i = 0; i < tesselation; i++) {

            // Calculate Vertex and add it to the mesh.
            float angle = angleStep * i;
            v.x = radius * cosf(angle);
            v.y = radius * sinf(angle);
            mesh->addVertex(VertexColor((transformation * v).xyz(), tgt::vec4(color, 1.0f)));

            // Calculate indices and add them to the mesh.
            if (k < centroid.getNumElements() - 1) {

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

GlMeshGeometryUInt32Color* StreamlineRenderer3D::createArrowPath(const StreamlineBundle& bundle) const {

    GlMeshGeometryUInt32Color* mesh = new GlMeshGeometryUInt32Color();
    mesh->setPrimitiveType(GL_TRIANGLES);

    // Define some constants for easier access.
    const uint32_t tesselation = GEOMETRY_TESSELATION;
    const float angleStep = (tgt::PIf * 2.0f) / tesselation;
    const float radius = bundle.getRadius();

    // Return an empty mesh in case of a zero radius.
    if(radius == 0.0f)
        return mesh;

    // Get the streamline-representation of the bundle in order to use the resample function.
    Streamline centroid = bundle.getCentroid();

    // Calculate the number of arrows the bundle will be split into.
    const float desiredLengthOfArrow = std::min(radius * 4.0f, centroid.getLength());
    const size_t samples = std::max<size_t>(2, centroid.getLength() / desiredLengthOfArrow);

    // Resample the bundle.
    centroid = centroid.resample(samples);

    // Specify the dimensions of each arrow-submesh.
    const float arrowRatio = 0.9f;
    const float lengthOfArrow = (centroid.getLength() / static_cast<float>(samples - 1)) * arrowRatio;
    const float lengthOfCylinder = lengthOfArrow * 0.7f;

    const float radiusOfCylinder = radius * 0.7f;
    const float radiusOfCone = radius * 1.1f;

    for (size_t k = 0; k < centroid.getNumElements() - 1; k++) {

        const tgt::vec3& p0 = centroid.getElementAt(k).position_;
        const tgt::vec3& p1 = centroid.getElementAt(k+1).position_;

        tgt::mat4 transformation = createTransformationMatrix(p0, p1 - p0);

        uint32_t offset = static_cast<uint32_t>(mesh->getVertices().size());

        tgt::vec3 color = centroid.getElementAt(k).velocity_;
        tgt::vec4 v(0.0f, 0.0f, 0.0f, 1.0f);

        // Cylinder
        mesh->addVertex(VertexColor(centroid.getElementAt(k).position_, tgt::vec4(color, 1.0f)));
        for (unsigned int i = 0; i < tesselation; i++) {

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
        for (unsigned int i = 0; i < tesselation; i++) {

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

tgt::mat4 StreamlineRenderer3D::createTransformationMatrix(const tgt::vec3& position, const tgt::vec3& velocity) const {

    tgt::vec3 tangent(tgt::normalize(velocity));

    tgt::vec3 temp(0.0f, 0.0f, 1.0f);
    if(1.0f - std::abs(tgt::dot(temp, tangent)) <= std::numeric_limits<float>::epsilon())
        temp = tgt::vec3(0.0f, 1.0f, 0.0f);

    tgt::vec3 binormal(tgt::normalize(tgt::cross(temp, tangent)));
    tgt::vec3 normal(tgt::normalize(tgt::cross(tangent, binormal)));

    return tgt::mat4(normal.x, binormal.x, tangent.x, position.x,
                     normal.y, binormal.y, tangent.y, position.y,
                     normal.z, binormal.z, tangent.z, position.z,
                     0.0f, 0.0f, 0.0f, 1.0f);
}

}   // namespace
