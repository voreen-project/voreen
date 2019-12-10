/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2019 University of Muenster, Germany,                        *
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

#include "geometryslicerenderer.h"

namespace voreen {

GeometrySliceRenderer::GeometrySliceRenderer()
    : RenderProcessor()
    , geometryInport_(Port::INPORT, "inputgeometry", "Input Geometry")
    , volumeInport_(Port::INPORT, "inputvolume", "Input Volume (for slice rendering)")
    , renderOutport_(Port::OUTPORT, "renderoutport", "Output Rendering (geometry slice view)", false, Processor::INVALID_RESULT, RenderPort::RENDERSIZE_RECEIVER)
    , sliceAlignment_("sliceAlignmentProp", "Slice Alignment (link with SliceViewer)")
    , sliceIndex_("sliceIndex", "Slice Number (link with SliceViewer)", 0, 0, 10000)
    , pickingMatrix_("pickingMatrix", "Picking Matrix (link with SliceViewer)", tgt::mat4::createIdentity(), tgt::mat4(-1e6f), tgt::mat4(1e6f))
    , enable_("enable", "Enable", true)
    , solidColor_("solidcolor", "Color")
    , alwaysFreeMemory_("alwaysfreememory", "Always free GPU memory", false, Processor::INVALID_RESULT, Property::LOD_ADVANCED)
    , clipShader_(0)
    , renderShader_(0)
    , numInputTriangles_(0)
    , writtenEdges_(0)
    , vertexArrayID_(0)
    , tfo_(0)
    , edgeVbo_(0)
    , edgeQuery_(0)
{
    ON_CHANGE(geometryInport_, GeometrySliceRenderer, adjustPropertyVisibility);
    addPort(geometryInport_);
    addPort(volumeInport_);
    addPort(renderOutport_);

    addProperty(enable_);

    sliceAlignment_.addOption("xy-plane", "XY-Plane (axial)", XY_PLANE);
    sliceAlignment_.addOption("xz-plane", "XZ-Plane (coronal)", XZ_PLANE);
    sliceAlignment_.addOption("yz-plane", "YZ-Plane (sagittal)", YZ_PLANE);

    addProperty(sliceAlignment_);
    addProperty(sliceIndex_);
    addProperty(pickingMatrix_);

    addProperty(solidColor_);

    addProperty(alwaysFreeMemory_);
}

Processor* GeometrySliceRenderer::create() const {
    return new GeometrySliceRenderer();
}

void GeometrySliceRenderer::initialize() {
    RenderProcessor::initialize();

    // load clip shader and set it up for the default vertex layout (color)
    currentVertexLayout_ = VertexBase::COLOR;
    std::vector<std::string> filenames(2);
    filenames[0] = "geometryslicerenderer_clipping.vert";
    filenames[1] = "geometryslicerenderer_clipping.geom";

    std::vector<tgt::ShaderObject::ShaderType> types(2);
    types[0] = tgt::ShaderObject::VERTEX_SHADER;
    types[1] = tgt::ShaderObject::GEOMETRY_SHADER;

    std::string shaderDefines = "#version 330 core \n #extension GL_ARB_gpu_shader5 : enable \n #define USE_COLOR \n";
    try {
        clipShader_ = ShdrMgr.loadSeparate(filenames, types, shaderDefines, false, false);
    } catch (tgt::Exception& e) {
        LERROR("Failed to load shader: " << e.what());
        clipShader_ = 0;
        return;
    }

    // prepare program for xfb
    static const char* xfbVaryings[] = { "edgeData.position1", "edgeData.color1", "edgeData.position2", "edgeData.color2" };
    glTransformFeedbackVaryings(clipShader_->getID(), 4, xfbVaryings, GL_INTERLEAVED_ATTRIBS);
    glLinkProgram(clipShader_->getID());

    GLint result = GL_FALSE;
    glGetProgramiv(clipShader_->getID(), GL_LINK_STATUS, &result);
    if (!result) {
        LERROR(clipShader_->getLinkerLog());
    }

    // load render shader
    try {
        renderShader_ = ShdrMgr.load("geometryslicerenderer_render", "", false);
    } catch (tgt::Exception& e) {
        LERROR("Failed to load shader: " << e.what());
        renderShader_ = 0;
        return;
    }

    glGenVertexArrays(1, &vertexArrayID_);
    glGenTransformFeedbacks(1, &tfo_);
    glGenBuffers(1, &edgeVbo_);
    glGenQueries(1, &edgeQuery_);
}

void GeometrySliceRenderer::deinitialize() {

    if(vertexArrayID_) {
        glDeleteVertexArrays(1, &vertexArrayID_);
    }
    if(tfo_) {
        glDeleteTransformFeedbacks(1, &tfo_);
    }
    if(edgeVbo_) {
        glDeleteBuffers(1, &edgeVbo_);
    }
    if(edgeQuery_) {
        glDeleteQueries(1, &edgeQuery_);
    }

    // dispose clip shader and render shader
    ShdrMgr.dispose(clipShader_);
    ShdrMgr.dispose(renderShader_);

    RenderProcessor::deinitialize();
}

void GeometrySliceRenderer::process() {

    tgtAssert(geometryInport_.getData(), "no input data");
    tgtAssert(volumeInport_.getData(), "no volume data");

    if (!enable_.get()) {
        renderOutport_.activateTarget();
        renderOutport_.clearTarget();
        renderOutport_.deactivateTarget();
        return;
    }

    if (!clipShader_ || !renderShader_) {
        LERROR("Invalid shader(s)!");
        return;
    }

    // check the input geometry, as only triangles are supported
    const GlMeshGeometryBase* glmgb = dynamic_cast<const GlMeshGeometryBase*>(geometryInport_.getData());
    bool primitiveType = glmgb && ( (glmgb->getPrimitiveType() == GL_TRIANGLES) ||  (glmgb->getPrimitiveType() == GL_TRIANGLE_STRIP) || (glmgb->getPrimitiveType() == GL_TRIANGLE_FAN) );
    if (!primitiveType) {
        LERROR("Only GlMeshGeometry with triangle types is currently supported");
        if (alwaysFreeMemory_.get())
            freeBufferMemory();
        renderOutport_.activateTarget();
        renderOutport_.clearTarget();
        renderOutport_.deactivateTarget();
        return;
    }
    else if (glmgb->isEmpty()) {
        // empty geometry -> do nothing but clear the outport
        if (alwaysFreeMemory_.get())
            freeBufferMemory();
        renderOutport_.activateTarget();
        renderOutport_.clearTarget();
        renderOutport_.deactivateTarget();
        return;
    }

    //check if the vertex layout has changed and update shader defines
    if (glmgb->getVertexLayout() != currentVertexLayout_) {
        currentVertexLayout_ = glmgb->getVertexLayout();
        updateShader(glmgb);
    }

    //get the number of triangles in the input geometry
    size_t meshTriangles = 0;
    if (glmgb->getPrimitiveType() == GL_TRIANGLES && !glmgb->usesIndexedDrawing())
        meshTriangles = glmgb->getNumVertices() / 3; // exactly (vertices / 3) triangles in geometry
    else
        meshTriangles = glmgb->usesIndexedDrawing() ? glmgb->getNumIndices() : glmgb->getNumVertices();  // FIXME: very bad approximation for number of triangles...
    if (numInputTriangles_ != meshTriangles) {
        numInputTriangles_ = meshTriangles;
        setupTransformFeedbackSizes();
    }

    // activate program
    clipShader_->activate();

    //disable rasterization and activate transform feedback
    glEnable(GL_RASTERIZER_DISCARD);
    glBindTransformFeedback(GL_TRANSFORM_FEEDBACK, tfo_);
    glBeginTransformFeedback(GL_POINTS);

    glBeginQueryIndexed(GL_TRANSFORM_FEEDBACK_PRIMITIVES_WRITTEN, 0, edgeQuery_);

    //get clipping plane and transform it to model space
    tgt::vec3 normal(0.f);
    switch (sliceAlignment_.getValue()) {
        case XY_PLANE:
            normal = tgt::vec3(0.f, 0.f, 1.f);
            break;
        case XZ_PLANE:
            normal = tgt::vec3(0.f, 1.f, 0.f);
            break;
        case YZ_PLANE:
            normal = tgt::vec3(1.f, 0.f, 0.f);
            break;
        default:
            break;
    }

    // define plane in voxel space and transform it first to world and then to model space
    tgt::plane plane = tgt::plane(normal, static_cast<float>(sliceIndex_.get()));
    plane = plane.transform(volumeInport_.getData()->getVoxelToWorldMatrix());
    plane = plane.transform(glmgb->getInvertedTransformationMatrix());

    clipShader_->setUniform("plane", plane.toVec4());

    if (!glmgb->supportsColors())
        clipShader_->setUniform("color", solidColor_.get());

    LGL_ERROR;
    glmgb->render();
    LGL_ERROR;

    glEndTransformFeedback();

    // each edge is writteb as a single vertex -> we just have to request this number
    glEndQueryIndexed(GL_TRANSFORM_FEEDBACK_PRIMITIVES_WRITTEN, 0);
    glGetQueryObjectuiv(edgeQuery_, GL_QUERY_RESULT, &writtenEdges_);

    //clean up OpenGL status
    glBindTransformFeedback(GL_TRANSFORM_FEEDBACK, 0);
    glDisable(GL_RASTERIZER_DISCARD);
    clipShader_->deactivate();

    LDEBUG("Recorded edges in clip edge buffer: " << writtenEdges_);

    renderOutport_.activateTarget();
    renderOutport_.clearTarget();

    LGL_ERROR;

    if (writtenEdges_ > 0) {
        renderShader_->activate();

        glBindVertexArray(vertexArrayID_);
        glBindBuffer(GL_ARRAY_BUFFER, edgeVbo_);

        // setup attributes
        glEnableVertexAttribArray(0);
        glVertexAttribPointer(0, 4, GL_FLOAT, GL_FALSE, 8 * sizeof(GLfloat), 0);
        glEnableVertexAttribArray(1);
        glVertexAttribPointer(1, 4, GL_FLOAT, GL_FALSE, 8 * sizeof(GLfloat), (void*) (4 * sizeof(GLfloat)));

        // compute the transformation matrix to project the world space geometry to the slice view
        float canvasWidth = static_cast<float>(renderOutport_.getSize().x);
        float canvasHeight = static_cast<float>(renderOutport_.getSize().y);

        // geometry has been computed in model space -> transform
        // model space -> world space
        tgt::mat4 transform = glmgb->getTransformationMatrix();

        // world space -> voxel space
        tgt::mat4 model = tgt::mat4::identity;
        volumeInport_.getData()->getVoxelToWorldMatrix().invert(model);
        transform = model * transform;

        // voxel space to viewport
        tgt::mat4 view = tgt::mat4::identity;
        pickingMatrix_.get().invert(view);
        transform = view * transform;

        // translate to correct z-coordinate
        tgt::mat4 translate = tgt::mat4::createTranslation(tgt::vec3(0.f, 0.f, static_cast<float>(-sliceIndex_.get())));
        transform = translate * transform;

        // project to canvas
        tgt::mat4 proj = tgt::mat4::createOrtho(0.f, canvasWidth, 0.f, canvasHeight, -1.f, 1.f);
        transform = proj * transform;

        renderShader_->setUniform("trafoMatrix_", transform);

        // The edges have been recorded into the transform feedback buffer -> they can be rendered directly
        glDrawArrays(GL_LINES, 0, writtenEdges_ * 2);    // multiply by 2 since each edge consists of two vertices

        // disable attributes
        glDisableVertexAttribArray(0);
        glDisableVertexAttribArray(1);

        glBindBuffer(GL_ARRAY_BUFFER, 0);
        glBindVertexArray(0);

        LGL_ERROR;

        renderShader_->deactivate();
    }

    LGL_ERROR;

    renderOutport_.deactivateTarget();

    LGL_ERROR;

    if (alwaysFreeMemory_.get())
        freeBufferMemory();
}

void GeometrySliceRenderer::updateShader(const GlMeshGeometryBase* mesh) {
    // set the defines
    clipShader_->setHeaders("#version 330 core \n #extension GL_ARB_gpu_shader5 : enable \n" + mesh->getShaderDefines());

    static const char* xfbVaryings[] = { "edgeData.position1", "edgeData.color1", "edgeData.position2", "edgeData.color2" };
    glTransformFeedbackVaryings(clipShader_->getID(), 4, xfbVaryings, GL_INTERLEAVED_ATTRIBS);

    clipShader_->rebuild();

    // check the link status after setting the transform feedback varyings
    GLint result = GL_FALSE;
    glGetProgramiv(clipShader_->getID(), GL_LINK_STATUS, &result);
    if (!result) {
        LERROR(clipShader_->getLinkerLog());
    }
}

void GeometrySliceRenderer::setupTransformFeedbackSizes() {

    // set up the buffer objects
    glBindTransformFeedback(GL_TRANSFORM_FEEDBACK, tfo_);
    glBindBuffer(GL_TRANSFORM_FEEDBACK_BUFFER, edgeVbo_);
    // this is the absolute maximum that will be recorded (3 edges per triangle can occur if the triangle lies within the clipping plane)
    glBufferData(GL_TRANSFORM_FEEDBACK_BUFFER, 3 * numInputTriangles_ * 16 * sizeof(GLfloat), NULL, GL_STREAM_DRAW);
    glBindBufferBase(GL_TRANSFORM_FEEDBACK_BUFFER, 0, edgeVbo_);
}

void GeometrySliceRenderer::freeBufferMemory() {
    //TODO: is there a better way to do this? Set the buffer size to 0 using glBufferData?
    glDeleteBuffers(1, &edgeVbo_);
    glGenBuffers(1, &edgeVbo_);
    numInputTriangles_ = 0;
}

void GeometrySliceRenderer::adjustPropertyVisibility() {
    if (!geometryInport_.getData())
        return;

    const GlMeshGeometryBase* glmgb = dynamic_cast<const GlMeshGeometryBase*>(geometryInport_.getData());

    if (glmgb && !glmgb->supportsColors())
        solidColor_.setVisibleFlag(true);
    else
        solidColor_.setVisibleFlag(false);
}

/*void GeometrySliceRenderer::retrieveClipEdges(GLuint edgeVbo, size_t edgeVertices, std::vector<tgt::vec4>& edges) {
    // retrieve the data from the GPU
    glBindBuffer(GL_TRANSFORM_FEEDBACK_BUFFER, edgeVbo);
    float* cpuData = new float[edgeVertices * 8];
    glGetBufferSubData(GL_TRANSFORM_FEEDBACK_BUFFER, 0, sizeof(GLfloat) * 8 * edgeVertices, &cpuData[0]);
    glBindBuffer(GL_TRANSFORM_FEEDBACK_BUFFER, 0);

    // convert the floats and store the edges
    for (size_t i = 0; i < edgeVertices; ++i) {

        tgt::vec4 pos = tgt::vec4(cpuData[i * 8], cpuData[i * 8 + 1], cpuData[i * 8 + 2], cpuData[i * 8 + 3]);
        tgt::vec4 color = tgt::vec4(cpuData[i * 8 + 4], cpuData[i * 8 + 5], cpuData[i * 8 + 6], cpuData[i * 8 + 7]);

        edges.push_back(pos);
    }

    delete[] cpuData;
}*/

} // namespace
