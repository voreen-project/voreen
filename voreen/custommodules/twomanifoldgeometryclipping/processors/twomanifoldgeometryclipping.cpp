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

#include "twomanifoldgeometryclipping.h"

namespace voreen {

using tgt::ShaderManager;    
using tgt::ShaderObject;
using tgt::Shader;

using tgt::mat4;
using tgt::mat3;
using tgt::transpose;
using tgt::vec3;
using tgt::vec2;
using tgt::vec4;
using tgt::ivec2;
using tgt::dot;
using tgt::cross;

const std::string TwoManifoldGeometryClipping::loggerCat_("voreen.TwoManifoldGeometryClipping.TwoManifoldGeometryClipping");

TwoManifoldGeometryClipping::TwoManifoldGeometryClipping()
    : Processor()
    , inport_(Port::INPORT, "geometry.input", "Geometry Input")
    , outport_(Port::OUTPORT, "geometry.output", "Geometry Output", false)
    , program_(0)  
    , enableClipping_("enable.clipping", "Enable Clipping", false)
    , setFixedClipColor_("setfixedcolor", "Set Fixed Color to Clipped Areas", false)
    , clipColor_("clipcolor", "Clip Area Color", tgt::vec4(1.f))
    , epsilon_("epsilon", "Epsilon Value (* 10e-4)", 0.f, 0.f, 10.f)   
    , onlyTriangulateTopLevel_("onlytriangulatetoplevel", "Only Triangulate Outer CCW Loops", false) 
    , alwaysFreeMemory_("alwaysfreememory", "Always free GPU memory", false) 
    , planeNormal_("plane.normal", "Clipping plane normal", tgt::vec3(0.f, 0.f, 1.f), tgt::vec3(-1.f), tgt::vec3(1.f))
    , planePosition_("plane.position", "Clipping plane position", 0.f, -10000.f, 10000.f)
    , clipMeshVbo_(0)
    , edgeVbo_(0)
    , tfo_(0)
    , clipQuery_(0)
    , edgeQuery_(0) 
    , writtenMeshVertices_(0)
    , writtenEdgeVertices_(0) 
    , outputGeometry_(0)  
    , currentVertexLayout_(VertexBase::COLOR_NORMAL) 
    , sizePerVertex_(44)  
    , numInputTriangles_(0)
{
    addPort(inport_);
    addPort(outport_);

    addProperty(enableClipping_);
    addProperty(planeNormal_);
    addProperty(planePosition_);
    addProperty(setFixedClipColor_);
    addProperty(clipColor_);
    addProperty(epsilon_);
    addProperty(onlyTriangulateTopLevel_);
    addProperty(alwaysFreeMemory_);

    //enableClipping_.onChange(CallMemberAction<TwoManifoldGeometryClipping>(this, &TwoManifoldGeometryClipping::updatePropertyVisibilities));
}

TwoManifoldGeometryClipping::~TwoManifoldGeometryClipping() {}

Processor* TwoManifoldGeometryClipping::create() const {
    return new TwoManifoldGeometryClipping();
}

void TwoManifoldGeometryClipping::initialize() {
    Processor::initialize();

    //load shader program and set it up for the default vertex layout (color and normals)
    std::vector<std::string> filenames(2);
    filenames[0] = "clipping.vert";
    filenames[1] = "clipping.geom";
    
    std::vector<ShaderObject::ShaderType> types(2);
    types[0] = ShaderObject::VERTEX_SHADER;
    types[1] = ShaderObject::GEOMETRY_SHADER;

    std::string shaderDefines = "#version 430 core\n #define USE_COLOR \n #define USE_NORMAL \n";
    program_ = ShdrMgr.loadSeparate(filenames, types, shaderDefines, false, false);

    // prepare program for xfb
    static const char* xfbVaryings[] = { "triangleData.position", "triangleData.color", "triangleData.normal", "gl_NextBuffer", "edgeData.position", "edgeData.color", "edgeData.normal", "edgeData.pos2d" };
    glTransformFeedbackVaryings(program_->getID(), 8, xfbVaryings, GL_INTERLEAVED_ATTRIBS);
    glLinkProgram(program_->getID());

    GLint result = GL_FALSE;
    glGetProgramiv(program_->getID(), GL_LINK_STATUS, &result);
    if (!result) {
        LERROR(program_->getLinkerLog());
    }

    glGenTransformFeedbacks(1, &tfo_);
    glGenBuffers(1, &clipMeshVbo_);
    glGenBuffers(1, &edgeVbo_);

    glGenQueries(1, &clipQuery_);
    glGenQueries(1, &edgeQuery_);
}

void TwoManifoldGeometryClipping::deinitialize() {

    glDeleteTransformFeedbacks(1, &tfo_);
    glDeleteBuffers(1, &clipMeshVbo_);
    glDeleteBuffers(1, &edgeVbo_);

    glDeleteQueries(1, &clipQuery_);
    glDeleteQueries(1, &edgeQuery_);

    ShdrMgr.dispose(program_);
    
    Processor::deinitialize();
}

void TwoManifoldGeometryClipping::updateShader(const GlMeshGeometryBase* mesh) {
    // set the defines
    program_->setHeaders("#version 430 core\n" + mesh->getShaderDefines());
    //program_->rebuild();

    // set the transform feedback varyings
    switch (currentVertexLayout_) {
        case VertexBase::SIMPLE : 
                    { 
                        static const char* xfbVaryings[] = { "triangleData.position", "gl_NextBuffer", "edgeData.position", "edgeData.pos2d" };
                        glTransformFeedbackVaryings(program_->getID(), 4, xfbVaryings, GL_INTERLEAVED_ATTRIBS);
                        //glLinkProgram(program_->getID());
                        break;
                    }
        case VertexBase::COLOR : 
                    { 
                        static const char* xfbVaryings[] = { "triangleData.position", "triangleData.color", "gl_NextBuffer", "edgeData.position", "edgeData.color", "edgeData.pos2d" };
                        glTransformFeedbackVaryings(program_->getID(), 6, xfbVaryings, GL_INTERLEAVED_ATTRIBS);
                        //glLinkProgram(program_->getID());
                        break;
                    }
        case VertexBase::NORMAL : 
                    { 
                        static const char* xfbVaryings[] = { "triangleData.position", "triangleData.normal", "gl_NextBuffer", "edgeData.position", "edgeData.normal", "edgeData.pos2d" };
                        glTransformFeedbackVaryings(program_->getID(), 6, xfbVaryings, GL_INTERLEAVED_ATTRIBS);
                        //glLinkProgram(program_->getID());
                        break;
                    }
        case VertexBase::COLOR_NORMAL : 
                    { 
                        static const char* xfbVaryings[] = { "triangleData.position", "triangleData.color", "triangleData.normal", "gl_NextBuffer", "edgeData.position", "edgeData.color", "edgeData.normal", "edgeData.pos2d" };
                        glTransformFeedbackVaryings(program_->getID(), 8, xfbVaryings, GL_INTERLEAVED_ATTRIBS);
                        //glLinkProgram(program_->getID());
                        break;
                    }
        case VertexBase::TEXCOORD : 
                    { 
                        static const char* xfbVaryings[] = { "triangleData.position", "triangleData.texCoord", "triangleData.texIndex", "gl_NextBuffer", 
                                                             "edgeData.position", "edgeData.texCoord", "edgeData.texIndex", "edgeData.pos2d" };
                        glTransformFeedbackVaryings(program_->getID(), 8, xfbVaryings, GL_INTERLEAVED_ATTRIBS);
                        //glLinkProgram(program_->getID());
                        break;
                    }
        case VertexBase::COLOR_TEXCOORD : 
                    { 
                        static const char* xfbVaryings[] = { "triangleData.position", "triangleData.color", "triangleData.texCoord", "triangleData.texIndex", "gl_NextBuffer", 
                                                             "edgeData.position", "edgeData.color", "edgeData.texCoord", "edgeData.texIndex", "edgeData.pos2d" };
                        glTransformFeedbackVaryings(program_->getID(), 10, xfbVaryings, GL_INTERLEAVED_ATTRIBS);
                        //glLinkProgram(program_->getID());
                        break;
                    }
        case VertexBase::NORMAL_TEXCOORD : 
                    { 
                        static const char* xfbVaryings[] = { "triangleData.position", "triangleData.normal", "triangleData.texCoord", "triangleData.texIndex", 
                                                             "gl_NextBuffer", "edgeData.position", "edgeData.normal", "edgeData.texCoord", "edgeData.texIndex", "edgeData.pos2d" };
                        glTransformFeedbackVaryings(program_->getID(), 10, xfbVaryings, GL_INTERLEAVED_ATTRIBS);
                        //glLinkProgram(program_->getID());
                        break;
                    }
        case VertexBase::COLOR_NORMAL_TEXCOORD: 
                    { 
                        static const char* xfbVaryings[] = { "triangleData.position", "triangleData.color", "triangleData.normal", "triangleData.texCoord", "triangleData.texIndex", 
                                                             "gl_NextBuffer", "edgeData.position", "edgeData.color", "edgeData.normal", "edgeData.texCoord", "edgeData.texIndex", "edgeData.pos2d" };
                        glTransformFeedbackVaryings(program_->getID(), 12, xfbVaryings, GL_INTERLEAVED_ATTRIBS);
                        //glLinkProgram(program_->getID());
                        break;
                    }
        default: break;
    }

    program_->rebuild();

    // check the link status after setting the transform feedback varyings
    GLint result = GL_FALSE;
    glGetProgramiv(program_->getID(), GL_LINK_STATUS, &result);
    if (!result) {
        LERROR(program_->getLinkerLog());
    }
}

void TwoManifoldGeometryClipping::updateSizePerVertex(const GlMeshGeometryBase* mesh) {
    size_t numFloats = 4; //position
    size_t numInts = 0;
    if (mesh->supportsColors())
        numFloats += 4;
    if (mesh->supportsNormals())
        numFloats += 3;
    if (mesh->supportsTextureData()) {
        numFloats += 2;
        numInts += 2;
    }

    sizePerVertex_ = numFloats * sizeof(GLfloat) + numInts * sizeof(GLint);
}

void TwoManifoldGeometryClipping::setupTransformFeedbackSizes() {

    // set up the buffer objects
    glBindTransformFeedback(GL_TRANSFORM_FEEDBACK, tfo_);
    glBindBuffer(GL_TRANSFORM_FEEDBACK_BUFFER, clipMeshVbo_);
    glBufferData(GL_TRANSFORM_FEEDBACK_BUFFER, 6 * numInputTriangles_ * sizePerVertex_, NULL, GL_STREAM_READ);
    glBindBufferBase(GL_TRANSFORM_FEEDBACK_BUFFER, 0, clipMeshVbo_);

    glBindBuffer(GL_TRANSFORM_FEEDBACK_BUFFER, edgeVbo_);
    // edges have additional 2D positions
    glBufferData(GL_TRANSFORM_FEEDBACK_BUFFER, 2 * numInputTriangles_ * (sizePerVertex_ + 2 * sizeof(GLfloat)), NULL, GL_STREAM_READ);
    glBindBufferBase(GL_TRANSFORM_FEEDBACK_BUFFER, 1, edgeVbo_);
}

void TwoManifoldGeometryClipping::freeBufferMemory() {
    //TODO: is there a better way to do this? Set the buffer size to 0 using glBufferData
    glDeleteBuffers(1, &clipMeshVbo_);
    glDeleteBuffers(1, &edgeVbo_);
    glGenBuffers(1, &clipMeshVbo_);
    glGenBuffers(1, &edgeVbo_);
    numInputTriangles_ = 0;
}

template<>
GlMeshGeometry<uint32_t, VertexBase>* TwoManifoldGeometryClipping::createEmptyMesh<VertexBase>(const GlMeshGeometryBase* inputMesh) {
    GlMeshGeometryUInt32Simple* mesh = new GlMeshGeometryUInt32Simple();
    mesh->setTransformationMatrix(inputMesh->getTransformationMatrix());
    return mesh;
}
   
template<>
GlMeshGeometry<uint32_t, VertexColor>* TwoManifoldGeometryClipping::createEmptyMesh<VertexColor>(const GlMeshGeometryBase* inputMesh) {
    GlMeshGeometryUInt32Color* mesh = new GlMeshGeometryUInt32Color();
    mesh->setTransformationMatrix(inputMesh->getTransformationMatrix());
    return mesh;
}

template<>
GlMeshGeometry<uint32_t, VertexNormal>* TwoManifoldGeometryClipping::createEmptyMesh<VertexNormal>(const GlMeshGeometryBase* inputMesh) {
    GlMeshGeometryUInt32Normal* mesh = new GlMeshGeometryUInt32Normal();
    mesh->setTransformationMatrix(inputMesh->getTransformationMatrix());
    return mesh;
}
    
template<>
GlMeshGeometry<uint32_t, VertexTexCoord>* TwoManifoldGeometryClipping::createEmptyMesh<VertexTexCoord>(const GlMeshGeometryBase* inputMesh) {
    GlMeshGeometryUInt32TexCoord* mesh = new GlMeshGeometryUInt32TexCoord();
    mesh->setTransformationMatrix(inputMesh->getTransformationMatrix());
    mesh->setTextureData(inputMesh->getTextureData());
    return mesh;
}

template<>
GlMeshGeometry<uint32_t, VertexColorNormal>* TwoManifoldGeometryClipping::createEmptyMesh<VertexColorNormal>(const GlMeshGeometryBase* inputMesh) {
    GlMeshGeometryUInt32ColorNormal* mesh = new GlMeshGeometryUInt32ColorNormal();
    mesh->setTransformationMatrix(inputMesh->getTransformationMatrix());
    return mesh;
}

template<>
GlMeshGeometry<uint32_t, VertexColorTexCoord>* TwoManifoldGeometryClipping::createEmptyMesh<VertexColorTexCoord>(const GlMeshGeometryBase* inputMesh) {
    GlMeshGeometryUInt32ColorTexCoord* mesh = new GlMeshGeometryUInt32ColorTexCoord();
    mesh->setTransformationMatrix(inputMesh->getTransformationMatrix());
    mesh->setTextureData(inputMesh->getTextureData());
    return mesh;
}
    
template<>
GlMeshGeometry<uint32_t, VertexNormalTexCoord>* TwoManifoldGeometryClipping::createEmptyMesh<VertexNormalTexCoord>(const GlMeshGeometryBase* inputMesh) {
    GlMeshGeometryUInt32NormalTexCoord* mesh = new GlMeshGeometryUInt32NormalTexCoord();
    mesh->setTransformationMatrix(inputMesh->getTransformationMatrix());
    mesh->setTextureData(inputMesh->getTextureData());
    return mesh;
}
    
template<>
GlMeshGeometry<uint32_t, VertexColorNormalTexCoord>* TwoManifoldGeometryClipping::createEmptyMesh<VertexColorNormalTexCoord>(const GlMeshGeometryBase* inputMesh) {
    GlMeshGeometryUInt32ColorNormalTexCoord* mesh = new GlMeshGeometryUInt32ColorNormalTexCoord();
    mesh->setTransformationMatrix(inputMesh->getTransformationMatrix());
    mesh->setTextureData(inputMesh->getTextureData());
    return mesh;
}

template<>
void TwoManifoldGeometryClipping::setFixedColors<VertexColor>(std::list<Loop<VertexColor>* >& polygons) {
    for (std::list<Loop<VertexColor>* >::iterator i = polygons.begin(); i != polygons.end(); ++i) {
        for (std::list<ClipVertex<VertexColor> >::iterator j = (*i)->vertices_.begin(); j != (*i)->vertices_.end(); ++j) 
            j->vertex_.setColor(clipColor_.get());
    }
}


template<>
void TwoManifoldGeometryClipping::setFixedColors<VertexColorNormal>(std::list<Loop<VertexColorNormal>* >& polygons) {
    for (std::list<Loop<VertexColorNormal>* >::iterator i = polygons.begin(); i != polygons.end(); ++i) {
        for (std::list<ClipVertex<VertexColorNormal> >::iterator j = (*i)->vertices_.begin(); j != (*i)->vertices_.end(); ++j) 
            j->vertex_.setColor(clipColor_.get());
    }
}


template<>
void TwoManifoldGeometryClipping::setFixedColors<VertexColorTexCoord>(std::list<Loop<VertexColorTexCoord>* >& polygons) {
    for (std::list<Loop<VertexColorTexCoord>* >::iterator i = polygons.begin(); i != polygons.end(); ++i) {
        for (std::list<ClipVertex<VertexColorTexCoord> >::iterator j = (*i)->vertices_.begin(); j != (*i)->vertices_.end(); ++j) 
            j->vertex_.setColor(clipColor_.get());
    }
}

template<>
void TwoManifoldGeometryClipping::setFixedColors<VertexColorNormalTexCoord>(std::list<Loop<VertexColorNormalTexCoord>* >& polygons) {
    for (std::list<Loop<VertexColorNormalTexCoord>* >::iterator i = polygons.begin(); i != polygons.end(); ++i) {
        for (std::list<ClipVertex<VertexColorNormalTexCoord> >::iterator j = (*i)->vertices_.begin(); j != (*i)->vertices_.end(); ++j) 
            j->vertex_.setColor(clipColor_.get());
    }
}

template<>
void TwoManifoldGeometryClipping::buildTriangleMeshFromBufferObject<VertexBase>(GLuint vbo, size_t meshVertices, GlMeshGeometry<uint32_t, VertexBase>* mesh) {
    // retrieve the data from the GPU
    glBindBuffer(GL_TRANSFORM_FEEDBACK_BUFFER, vbo);
    float* cpuData = new float[meshVertices * 4];
    glGetBufferSubData(GL_TRANSFORM_FEEDBACK_BUFFER, 0, sizeof(GLfloat) * 4 * meshVertices, &cpuData[0]);
    glBindBuffer(GL_TRANSFORM_FEEDBACK_BUFFER, 0);

    // convert the floats to triangles and add them to the geometry
    VertexBase v[3];
    for (size_t i = 0; i < meshVertices; ++i) {  
        vec4 pos = vec4(cpuData[i * 4], cpuData[i * 4 + 1], cpuData[i * 4 + 2], cpuData[i * 4 + 3]);
        v[i % 3] = VertexBase(pos.xyz());
        //check if three consecutive vertices have been read in 
        if ((i % 3) == 2) {
            //new triangle complete
            mesh->addVertex(v[0]);
            mesh->addVertex(v[1]);
            mesh->addVertex(v[2]);
        }
    }

    delete[] cpuData;
}

template<>
void TwoManifoldGeometryClipping::buildTriangleMeshFromBufferObject<VertexColor>(GLuint vbo, size_t meshVertices, GlMeshGeometry<uint32_t, VertexColor>* mesh) {
    // retrieve the data from the GPU
    glBindBuffer(GL_TRANSFORM_FEEDBACK_BUFFER, vbo);
    float* cpuData = new float[meshVertices * 8];
    glGetBufferSubData(GL_TRANSFORM_FEEDBACK_BUFFER, 0, sizeof(GLfloat) * 8 * meshVertices, &cpuData[0]);
    glBindBuffer(GL_TRANSFORM_FEEDBACK_BUFFER, 0);

    // convert the floats to triangles and add them to the geometry
    VertexColor v[3];
    for (size_t i = 0; i < meshVertices; ++i) {
        
        vec4 pos = vec4(cpuData[i * 8], cpuData[i * 8 + 1], cpuData[i * 8 + 2], cpuData[i * 8 + 3]);
        vec4 color = vec4(cpuData[i * 8 + 4], cpuData[i * 8 + 5], cpuData[i * 8 + 6], cpuData[i * 8 + 7]);
        v[i % 3] = VertexColor(pos.xyz(), color);

        //check if three consecutive vertices have been read in 
        if ((i % 3) == 2) {
            //new triangle complete
            mesh->addVertex(v[0]);
            mesh->addVertex(v[1]);
            mesh->addVertex(v[2]);
        }
    }

    delete[] cpuData;
}

template<>
void TwoManifoldGeometryClipping::buildTriangleMeshFromBufferObject<VertexNormal>(GLuint vbo, size_t meshVertices, GlMeshGeometry<uint32_t, VertexNormal>* mesh) {
    // retrieve the data from the GPU
    glBindBuffer(GL_TRANSFORM_FEEDBACK_BUFFER, vbo);
    float* cpuData = new float[meshVertices * 7];
    glGetBufferSubData(GL_TRANSFORM_FEEDBACK_BUFFER, 0, sizeof(GLfloat) * 7 * meshVertices, &cpuData[0]);
    glBindBuffer(GL_TRANSFORM_FEEDBACK_BUFFER, 0);

    // convert the floats to triangles and add them to the geometry
    VertexNormal v[3];
    for (size_t i = 0; i < meshVertices; ++i) {
        
        vec4 pos = vec4(cpuData[i * 7], cpuData[i * 7 + 1], cpuData[i * 7 + 2], cpuData[i * 7 + 3]);
        vec3 normal = vec3(cpuData[i * 7 + 4], cpuData[i * 7 + 5], cpuData[i * 7 + 6]);
        v[i % 3] = VertexNormal(pos.xyz(), normal);

        //check if three consecutive vertices have been read in 
        if ((i % 3) == 2) {
            //new triangle complete
            mesh->addVertex(v[0]);
            mesh->addVertex(v[1]);
            mesh->addVertex(v[2]);
        }
    }

    delete[] cpuData;
}


template<>
void TwoManifoldGeometryClipping::buildTriangleMeshFromBufferObject<VertexColorNormal>(GLuint vbo, size_t meshVertices, GlMeshGeometry<uint32_t, VertexColorNormal>* mesh) {
    // retrieve the data from the GPU
    glBindBuffer(GL_TRANSFORM_FEEDBACK_BUFFER, vbo);
    float* cpuData = new float[meshVertices * 11];
    glGetBufferSubData(GL_TRANSFORM_FEEDBACK_BUFFER, 0, sizeof(GLfloat) * 11 * meshVertices, &cpuData[0]);
    glBindBuffer(GL_TRANSFORM_FEEDBACK_BUFFER, 0);

    // convert the floats to triangles and add them to the geometry
    VertexColorNormal v[3];
    for (size_t i = 0; i < meshVertices; ++i) {
        
        vec4 pos = vec4(cpuData[i * 11], cpuData[i * 11 + 1], cpuData[i * 11 + 2], cpuData[i * 11 + 3]);
        vec4 color = vec4(cpuData[i * 11 + 4], cpuData[i * 11 + 5], cpuData[i * 11 + 6], cpuData[i * 11 + 7]);
        vec3 normal = vec3(cpuData[i * 11 + 8], cpuData[i * 11 + 9], cpuData[i * 11 + 10]);
        v[i % 3] = VertexColorNormal(pos.xyz(), color, normal);

        //check if three consecutive vertices have been read in 
        if ((i % 3) == 2) {
            //new triangle complete
            mesh->addVertex(v[0]);
            mesh->addVertex(v[1]);
            mesh->addVertex(v[2]);
        }
    }

    delete[] cpuData;
}

template<>
void TwoManifoldGeometryClipping::buildTriangleMeshFromBufferObject<VertexTexCoord>(GLuint vbo, size_t meshVertices, GlMeshGeometry<uint32_t, VertexTexCoord>* mesh) {
    // retrieve the data from the GPU -> void* because there are also integers present
    glBindBuffer(GL_TRANSFORM_FEEDBACK_BUFFER, vbo);
    void* cpuData = operator new(meshVertices * sizePerVertex_);
    glGetBufferSubData(GL_TRANSFORM_FEEDBACK_BUFFER, 0, meshVertices * sizePerVertex_, cpuData);
    glBindBuffer(GL_TRANSFORM_FEEDBACK_BUFFER, 0);

    // convert the data to triangles and add them to the geometry
    VertexTexCoord v[3];
    char* current;
    size_t i = 0;
    for (current = reinterpret_cast<char*>(cpuData); current < (reinterpret_cast<char*>(cpuData) + meshVertices * sizePerVertex_); current += sizePerVertex_) {
        float* floatCurrent = reinterpret_cast<float*>(current);  
        vec4 pos = vec4(*floatCurrent, *(floatCurrent + 1), *(floatCurrent + 2), *(floatCurrent + 3));
        vec2 texCoord = vec2(*(floatCurrent + 4), *(floatCurrent + 5));
        int* intCurrent = reinterpret_cast<int*>(floatCurrent + 6);
        ivec2 texIndex = ivec2(*intCurrent, *(intCurrent + 1));
        v[i % 3] = VertexTexCoord(pos.xyz(), texCoord, texIndex);
        //check if three consecutive vertices have been read in 
        if ((i % 3) == 2) {
            //new triangle complete
            mesh->addVertex(v[0]);
            mesh->addVertex(v[1]);
            mesh->addVertex(v[2]);
        }
        i++;
    }

    operator delete(cpuData);
}

template<>
void TwoManifoldGeometryClipping::buildTriangleMeshFromBufferObject<VertexColorTexCoord>(GLuint vbo, size_t meshVertices, GlMeshGeometry<uint32_t, VertexColorTexCoord>* mesh) {
    // retrieve the data from the GPU -> void* because there are also integers present
    glBindBuffer(GL_TRANSFORM_FEEDBACK_BUFFER, vbo);
    void* cpuData = operator new(meshVertices * sizePerVertex_);
    glGetBufferSubData(GL_TRANSFORM_FEEDBACK_BUFFER, 0, meshVertices * sizePerVertex_, cpuData);
    glBindBuffer(GL_TRANSFORM_FEEDBACK_BUFFER, 0);

    // convert the data to triangles and add them to the geometry
    VertexColorTexCoord v[3];
    char* current;
    size_t i = 0;
    for (current = reinterpret_cast<char*>(cpuData); current < (reinterpret_cast<char*>(cpuData) + meshVertices * sizePerVertex_); current += sizePerVertex_) {
        float* floatCurrent = reinterpret_cast<float*>(current);  
        vec4 pos = vec4(*floatCurrent, *(floatCurrent + 1), *(floatCurrent + 2), *(floatCurrent + 3));
        vec4 color = vec4(*(floatCurrent + 4), *(floatCurrent + 5), *(floatCurrent + 6), *(floatCurrent + 7));
        vec2 texCoord = vec2(*(floatCurrent + 8), *(floatCurrent + 9));

    /*

    //TODO: free vertex buffer memory?
    
    
    //addCapsToMesh(edgeVbo_, writtenEdgeVertices_, outputMesh);
    
    */
        int* intCurrent = reinterpret_cast<int*>(floatCurrent + 10);
        ivec2 texIndex = ivec2(*intCurrent, *(intCurrent + 1));
        v[i % 3] = VertexColorTexCoord(pos.xyz(), color, texCoord, texIndex);
        //check if three consecutive vertices have been read in 
        if ((i % 3) == 2) {
            //new triangle complete
            mesh->addVertex(v[0]);
            mesh->addVertex(v[1]);
            mesh->addVertex(v[2]);
        }
        i++;
    }

    operator delete(cpuData);
}

template<>
void TwoManifoldGeometryClipping::buildTriangleMeshFromBufferObject<VertexNormalTexCoord>(GLuint vbo, size_t meshVertices, GlMeshGeometry<uint32_t, VertexNormalTexCoord>* mesh) {
    // retrieve the data from the GPU -> void* because there are also integers present
    glBindBuffer(GL_TRANSFORM_FEEDBACK_BUFFER, vbo);
    void* cpuData = operator new(meshVertices * sizePerVertex_);
    glGetBufferSubData(GL_TRANSFORM_FEEDBACK_BUFFER, 0, meshVertices * sizePerVertex_, cpuData);
    glBindBuffer(GL_TRANSFORM_FEEDBACK_BUFFER, 0);

    // convert the data to triangles and add them to the geometry
    VertexNormalTexCoord v[3];
    char* current;
    size_t i = 0;
    for (current = reinterpret_cast<char*>(cpuData); current < (reinterpret_cast<char*>(cpuData) + meshVertices * sizePerVertex_); current += sizePerVertex_) {
        float* floatCurrent = reinterpret_cast<float*>(current);  
        vec4 pos = vec4(*floatCurrent, *(floatCurrent + 1), *(floatCurrent + 2), *(floatCurrent + 3));
        vec3 normal = vec3(*(floatCurrent + 4), *(floatCurrent + 5), *(floatCurrent + 6));
        vec2 texCoord = vec2(*(floatCurrent + 7), *(floatCurrent + 8));
        int* intCurrent = reinterpret_cast<int*>(floatCurrent + 9);
        ivec2 texIndex = ivec2(*intCurrent, *(intCurrent + 1));
        v[i % 3] = VertexNormalTexCoord(pos.xyz(), normal, texCoord, texIndex);
        //check if three consecutive vertices have been read in 
        if ((i % 3) == 2) {
            //new triangle complete
            mesh->addVertex(v[0]);
            mesh->addVertex(v[1]);
            mesh->addVertex(v[2]);
        }
        i++;
    }

    operator delete(cpuData);
}

template<>
void TwoManifoldGeometryClipping::buildTriangleMeshFromBufferObject<VertexColorNormalTexCoord>(GLuint vbo, size_t meshVertices, GlMeshGeometry<uint32_t, VertexColorNormalTexCoord>* mesh) {
    // retrieve the data from the GPU -> void* because there are also integers present
    glBindBuffer(GL_TRANSFORM_FEEDBACK_BUFFER, vbo);
    void* cpuData = operator new(meshVertices * sizePerVertex_);
    glGetBufferSubData(GL_TRANSFORM_FEEDBACK_BUFFER, 0, meshVertices * sizePerVertex_, cpuData);
    glBindBuffer(GL_TRANSFORM_FEEDBACK_BUFFER, 0);

    // convert the data to triangles and add them to the geometry
    VertexColorNormalTexCoord v[3];
    char* current;
    size_t i = 0;
    for (current = reinterpret_cast<char*>(cpuData); current < (reinterpret_cast<char*>(cpuData) + meshVertices * sizePerVertex_); current += sizePerVertex_) {
        float* floatCurrent = reinterpret_cast<float*>(current);  
        vec4 pos = vec4(*floatCurrent, *(floatCurrent + 1), *(floatCurrent + 2), *(floatCurrent + 3));
        vec4 color = vec4(*(floatCurrent + 4), *(floatCurrent + 5), *(floatCurrent + 6), *(floatCurrent + 7));
        vec3 normal = vec3(*(floatCurrent + 8), *(floatCurrent + 9), *(floatCurrent + 10));
        vec2 texCoord = vec2(*(floatCurrent + 11), *(floatCurrent + 12));
        int* intCurrent = reinterpret_cast<int*>(floatCurrent + 13);
        ivec2 texIndex = ivec2(*intCurrent, *(intCurrent + 1));
        v[i % 3] = VertexColorNormalTexCoord(pos.xyz(), color, normal, texCoord, texIndex);
        //check if three consecutive vertices have been read in 
        if ((i % 3) == 2) {
            //new triangle complete
            mesh->addVertex(v[0]);
            mesh->addVertex(v[1]);
            mesh->addVertex(v[2]);
        }
        i++;
    }

    operator delete(cpuData);
}

template<>
void TwoManifoldGeometryClipping::retrieveClipEdges<VertexBase>(GLuint edgeVbo, size_t edgeVertices, std::vector<ClipEdge<VertexBase> >& edges) {
    // retrieve the data from the GPU
    glBindBuffer(GL_TRANSFORM_FEEDBACK_BUFFER, edgeVbo);
    float* cpuData = new float[edgeVertices * 6];
    glGetBufferSubData(GL_TRANSFORM_FEEDBACK_BUFFER, 0, sizeof(GLfloat) * 6 * edgeVertices, &cpuData[0]);
    glBindBuffer(GL_TRANSFORM_FEEDBACK_BUFFER, 0);

    // convert the floats and store the edges
    ClipVertex<VertexBase> v[2];
    for (size_t i = 0; i < edgeVertices; ++i) {
        
        vec4 pos = vec4(cpuData[i * 6], cpuData[i * 6 + 1], cpuData[i * 6 + 2], cpuData[i * 6 + 3]);
        VertexBase vcn(pos.xyz());
        vec2 pos2d = vec2(cpuData[i * 6 + 4], cpuData[i * 6 + 5]);
        ClipVertex<VertexBase> cv; 
        cv.vertex_ = vcn;
        cv.pos2d_ = pos2d;
        cv.earTipStatus_ = false;
        v[i % 2] = cv;

        //check if two consecutive vertices have been read in to get rid of degenerated edges
        if ((i % 2) == 1 && !(v[0] == v[1])) {
            //new edge complete
            ClipEdge<VertexBase> e;
            e.start_ = v[0];
            e.end_ = v[1];
            edges.push_back(e);
        }
    }

    delete[] cpuData;
}

template<>
void TwoManifoldGeometryClipping::retrieveClipEdges<VertexColor>(GLuint edgeVbo, size_t edgeVertices, std::vector<ClipEdge<VertexColor> >& edges) {
    // retrieve the data from the GPU
    glBindBuffer(GL_TRANSFORM_FEEDBACK_BUFFER, edgeVbo);
    float* cpuData = new float[edgeVertices * 10];
    glGetBufferSubData(GL_TRANSFORM_FEEDBACK_BUFFER, 0, sizeof(GLfloat) * 10 * edgeVertices, &cpuData[0]);
    glBindBuffer(GL_TRANSFORM_FEEDBACK_BUFFER, 0);

    // convert the floats and store the edges
    ClipVertex<VertexColor> v[2];
    for (size_t i = 0; i < edgeVertices; ++i) {
        
        vec4 pos = vec4(cpuData[i * 10], cpuData[i * 10 + 1], cpuData[i * 10 + 2], cpuData[i * 10 + 3]);
        vec4 color = vec4(cpuData[i * 10 + 4], cpuData[i * 10 + 5], cpuData[i * 10 + 6], cpuData[i * 10 + 7]);
        VertexColor vcn(pos.xyz(), color);
        vec2 pos2d = vec2(cpuData[i * 10 + 8], cpuData[i * 10 + 9]);
        ClipVertex<VertexColor> cv; 
        cv.vertex_ = vcn;
        cv.pos2d_ = pos2d;
        cv.earTipStatus_ = false;
        v[i % 2] = cv;

        //check if two consecutive vertices have been read in to get rid of degenerated edges
        if ((i % 2) == 1 && !(v[0] == v[1])) {
            //new edge complete
            ClipEdge<VertexColor> e;
            e.start_ = v[0];
            e.end_ = v[1];
            edges.push_back(e);
        }
    }

    delete[] cpuData;
}

template<>
void TwoManifoldGeometryClipping::retrieveClipEdges<VertexNormal>(GLuint edgeVbo, size_t edgeVertices, std::vector<ClipEdge<VertexNormal> >& edges) {
    // retrieve the data from the GPU
    glBindBuffer(GL_TRANSFORM_FEEDBACK_BUFFER, edgeVbo);
    float* cpuData = new float[edgeVertices * 9];
    glGetBufferSubData(GL_TRANSFORM_FEEDBACK_BUFFER, 0, sizeof(GLfloat) * 9 * edgeVertices, &cpuData[0]);
    glBindBuffer(GL_TRANSFORM_FEEDBACK_BUFFER, 0);

    // convert the floats and store the edges
    ClipVertex<VertexNormal> v[2];
    for (size_t i = 0; i < edgeVertices; ++i) {
        
        vec4 pos = vec4(cpuData[i * 9], cpuData[i * 9 + 1], cpuData[i * 9 + 2], cpuData[i * 9 + 3]);
        //vec3 normal = vec3(cpuData[i * 9 + 4], cpuData[i * 9 + 5], cpuData[i * 9 + 6]);
        VertexNormal vcn(pos.xyz(), capNormal_);
        vec2 pos2d = vec2(cpuData[i * 9 + 7], cpuData[i * 9 + 8]);
        ClipVertex<VertexNormal> cv; 
        cv.vertex_ = vcn;
        cv.pos2d_ = pos2d;
        cv.earTipStatus_ = false;
        v[i % 2] = cv;

        //check if two consecutive vertices have been read in to get rid of degenerated edges
        if ((i % 2) == 1 && !(v[0] == v[1])) {
            //new edge complete
            ClipEdge<VertexNormal> e;
            e.start_ = v[0];
            e.end_ = v[1];
            edges.push_back(e);
        }
    }

    delete[] cpuData;
}

template<>
void TwoManifoldGeometryClipping::retrieveClipEdges<VertexColorNormal>(GLuint edgeVbo, size_t edgeVertices, std::vector<ClipEdge<VertexColorNormal> >& edges) {
    // retrieve the data from the GPU
    glBindBuffer(GL_TRANSFORM_FEEDBACK_BUFFER, edgeVbo);
    float* cpuData = new float[edgeVertices * 13];
    glGetBufferSubData(GL_TRANSFORM_FEEDBACK_BUFFER, 0, sizeof(GLfloat) * 13 * edgeVertices, &cpuData[0]);
    glBindBuffer(GL_TRANSFORM_FEEDBACK_BUFFER, 0);

    // convert the floats and store the edges
    ClipVertex<VertexColorNormal> v[2];
    for (size_t i = 0; i < edgeVertices; ++i) {
        
        vec4 pos = vec4(cpuData[i * 13], cpuData[i * 13 + 1], cpuData[i * 13 + 2], cpuData[i * 13 + 3]);
        vec4 color = vec4(cpuData[i * 13 + 4], cpuData[i * 13 + 5], cpuData[i * 13 + 6], cpuData[i * 13 + 7]);
        //vec3 normal = vec3(cpuData[i * 13 + 8], cpuData[i * 13 + 9], cpuData[i * 13 + 10]);
        //this sets the inverted plane normal
        VertexColorNormal vcn(pos.xyz(), color, capNormal_);
        vec2 pos2d = vec2(cpuData[i * 13 + 11], cpuData[i * 13 + 12]);
        ClipVertex<VertexColorNormal> cv; 
        cv.vertex_ = vcn;
        cv.pos2d_ = pos2d;
        cv.earTipStatus_ = false;
        v[i % 2] = cv;

        //check if two consecutive vertices have been read in to get rid of degenerated edges
        if ((i % 2) == 1 && !(v[0] == v[1])) {
            //new edge complete
            ClipEdge<VertexColorNormal> e;
            e.start_ = v[0];
            e.end_ = v[1];
            edges.push_back(e);
        }
    }

    delete[] cpuData;
}

template<>
void TwoManifoldGeometryClipping::retrieveClipEdges<VertexTexCoord>(GLuint edgeVbo, size_t edgeVertices, std::vector<ClipEdge<VertexTexCoord> >& edges) {
    // retrieve the data from the GPU -> void* because there are also integers present
    glBindBuffer(GL_TRANSFORM_FEEDBACK_BUFFER, edgeVbo);
    size_t actualVertexSize = sizePerVertex_ + 2 * sizeof(GLfloat);
    void* cpuData = operator new(edgeVertices * actualVertexSize);
    glGetBufferSubData(GL_TRANSFORM_FEEDBACK_BUFFER, 0, edgeVertices * actualVertexSize, cpuData);
    glBindBuffer(GL_TRANSFORM_FEEDBACK_BUFFER, 0);

    // convert the data and store the edges
    ClipVertex<VertexTexCoord> v[2];
    char* current;
    size_t i = 0;
    for (current = reinterpret_cast<char*>(cpuData); current < (reinterpret_cast<char*>(cpuData) + edgeVertices * actualVertexSize); current += actualVertexSize) {
        float* floatCurrent = reinterpret_cast<float*>(current);  
        vec4 pos = vec4(*floatCurrent, *(floatCurrent + 1), *(floatCurrent + 2), *(floatCurrent + 3));
        vec2 texCoord = vec2(*(floatCurrent + 4), *(floatCurrent + 5));
        int* intCurrent = reinterpret_cast<int*>(floatCurrent + 6);
        ivec2 texIndex = ivec2(*intCurrent, *(intCurrent + 1));
        VertexTexCoord vcn(pos.xyz(), texCoord, texIndex);
        floatCurrent = reinterpret_cast<float*>(intCurrent + 2);
        vec2 pos2d = vec2(*floatCurrent, *(floatCurrent + 1));
        ClipVertex<VertexTexCoord> cv;
        cv.vertex_ = vcn;
        cv.pos2d_ = pos2d;
        cv.earTipStatus_ = false;
        v[i % 2] = cv;
        //check if two consecutive vertices have been read in and get rid of degenerated edges
        if ((i % 2) == 1 && !(v[0] == v[1])) {
            //new edge complete
            ClipEdge<VertexTexCoord> e;
            e.start_ = v[0];
            e.end_ = v[1];
            edges.push_back(e);
        }
        i++;
    }

    operator delete(cpuData);
}

template<>
void TwoManifoldGeometryClipping::retrieveClipEdges<VertexColorTexCoord>(GLuint edgeVbo, size_t edgeVertices, std::vector<ClipEdge<VertexColorTexCoord> >& edges) {
    // retrieve the data from the GPU -> void* because there are also integers present
    glBindBuffer(GL_TRANSFORM_FEEDBACK_BUFFER, edgeVbo);
    size_t actualVertexSize = sizePerVertex_ + 2 * sizeof(GLfloat);
    void* cpuData = operator new(edgeVertices * actualVertexSize);
    glGetBufferSubData(GL_TRANSFORM_FEEDBACK_BUFFER, 0, edgeVertices * actualVertexSize, cpuData);
    glBindBuffer(GL_TRANSFORM_FEEDBACK_BUFFER, 0);

    // convert the data and store the edges
    ClipVertex<VertexColorTexCoord> v[2];
    char* current;
    size_t i = 0;
    for (current = reinterpret_cast<char*>(cpuData); current < (reinterpret_cast<char*>(cpuData) + edgeVertices * actualVertexSize); current += actualVertexSize) {
        float* floatCurrent = reinterpret_cast<float*>(current);  
        vec4 pos = vec4(*floatCurrent, *(floatCurrent + 1), *(floatCurrent + 2), *(floatCurrent + 3));
        vec4 color = vec4(*(floatCurrent + 4), *(floatCurrent + 5), *(floatCurrent + 6), *(floatCurrent + 7));
        vec2 texCoord = vec2(*(floatCurrent + 8), *(floatCurrent + 9));
        int* intCurrent = reinterpret_cast<int*>(floatCurrent + 10);
        ivec2 texIndex = ivec2(*intCurrent, *(intCurrent + 1));
        VertexColorTexCoord vcn(pos.xyz(), color, texCoord, texIndex);
        floatCurrent = reinterpret_cast<float*>(intCurrent + 2);
        vec2 pos2d = vec2(*floatCurrent, *(floatCurrent + 1));
        ClipVertex<VertexColorTexCoord> cv;
        cv.vertex_ = vcn;
        cv.pos2d_ = pos2d;
        cv.earTipStatus_ = false;
        v[i % 2] = cv;
        //check if two consecutive vertices have been read in and get rid of degenerated edges
        if ((i % 2) == 1 && !(v[0] == v[1])) {
            //new edge complete
            ClipEdge<VertexColorTexCoord> e;
            e.start_ = v[0];
            e.end_ = v[1];
            edges.push_back(e);
        }
        i++;
    }

    operator delete(cpuData);
}

template<>
void TwoManifoldGeometryClipping::retrieveClipEdges<VertexNormalTexCoord>(GLuint edgeVbo, size_t edgeVertices, std::vector<ClipEdge<VertexNormalTexCoord> >& edges) {
    // retrieve the data from the GPU -> void* because there are also integers present
    glBindBuffer(GL_TRANSFORM_FEEDBACK_BUFFER, edgeVbo);
    size_t actualVertexSize = sizePerVertex_ + 2 * sizeof(GLfloat);
    void* cpuData = operator new(edgeVertices * actualVertexSize);
    glGetBufferSubData(GL_TRANSFORM_FEEDBACK_BUFFER, 0, edgeVertices * actualVertexSize, cpuData);
    glBindBuffer(GL_TRANSFORM_FEEDBACK_BUFFER, 0);

    // convert the data and store the edges
    ClipVertex<VertexNormalTexCoord> v[2];
    char* current;
    size_t i = 0;
    for (current = reinterpret_cast<char*>(cpuData); current < (reinterpret_cast<char*>(cpuData) + edgeVertices * actualVertexSize); current += actualVertexSize) {
        float* floatCurrent = reinterpret_cast<float*>(current);  
        vec4 pos = vec4(*floatCurrent, *(floatCurrent + 1), *(floatCurrent + 2), *(floatCurrent + 3));
        // normal not needed
        //vec3 normal = vec3(*(floatCurrent + 4), *(floatCurrent + 5), *(floatCurrent + 6));
        vec2 texCoord = vec2(*(floatCurrent + 7), *(floatCurrent + 8));
        int* intCurrent = reinterpret_cast<int*>(floatCurrent + 9);
        ivec2 texIndex = ivec2(*intCurrent, *(intCurrent + 1));
        VertexNormalTexCoord vcn(pos.xyz(), capNormal_, texCoord, texIndex);
        floatCurrent = reinterpret_cast<float*>(intCurrent + 2);
        vec2 pos2d = vec2(*floatCurrent, *(floatCurrent + 1));
        ClipVertex<VertexNormalTexCoord> cv;
        cv.vertex_ = vcn;
        cv.pos2d_ = pos2d;
        cv.earTipStatus_ = false;
        v[i % 2] = cv;
        //check if two consecutive vertices have been read in and get rid of degenerated edges
        if ((i % 2) == 1 && !(v[0] == v[1])) {
            //new edge complete
            ClipEdge<VertexNormalTexCoord> e;
            e.start_ = v[0];
            e.end_ = v[1];
            edges.push_back(e);
        }
        i++;
    }

    operator delete(cpuData);
}


template<>
void TwoManifoldGeometryClipping::retrieveClipEdges<VertexColorNormalTexCoord>(GLuint edgeVbo, size_t edgeVertices, std::vector<ClipEdge<VertexColorNormalTexCoord> >& edges) {
    // retrieve the data from the GPU -> void* because there are also integers present
    glBindBuffer(GL_TRANSFORM_FEEDBACK_BUFFER, edgeVbo);
    size_t actualVertexSize = sizePerVertex_ + 2 * sizeof(GLfloat);
    void* cpuData = operator new(edgeVertices * actualVertexSize);
    glGetBufferSubData(GL_TRANSFORM_FEEDBACK_BUFFER, 0, edgeVertices * actualVertexSize, cpuData);
    glBindBuffer(GL_TRANSFORM_FEEDBACK_BUFFER, 0);

    // convert the data and store the edges
    ClipVertex<VertexColorNormalTexCoord> v[2];
    char* current;
    size_t i = 0;
    for (current = reinterpret_cast<char*>(cpuData); current < (reinterpret_cast<char*>(cpuData) + edgeVertices * actualVertexSize); current += actualVertexSize) {
        float* floatCurrent = reinterpret_cast<float*>(current);  
        vec4 pos = vec4(*floatCurrent, *(floatCurrent + 1), *(floatCurrent + 2), *(floatCurrent + 3));
        vec4 color = vec4(*(floatCurrent + 4), *(floatCurrent + 5), *(floatCurrent + 6), *(floatCurrent + 7));
        // normal not needed
        //vec3 normal = vec3(*(floatCurrent + 8), *(floatCurrent + 9), *(floatCurrent + 10));
        vec2 texCoord = vec2(*(floatCurrent + 11), *(floatCurrent + 12));
        int* intCurrent = reinterpret_cast<int*>(floatCurrent + 13);
        ivec2 texIndex = ivec2(*intCurrent, *(intCurrent + 1));
        VertexColorNormalTexCoord vcn(pos.xyz(), color, capNormal_, texCoord, texIndex);
        floatCurrent = reinterpret_cast<float*>(intCurrent + 2);
        vec2 pos2d = vec2(*floatCurrent, *(floatCurrent + 1));
        ClipVertex<VertexColorNormalTexCoord> cv;
        cv.vertex_ = vcn;
        cv.pos2d_ = pos2d;
        cv.earTipStatus_ = false;
        v[i % 2] = cv;
        //check if two consecutive vertices have been read in and get rid of degenerated edges
        if ((i % 2) == 1 && !(v[0] == v[1])) {
            //new edge complete
            ClipEdge<VertexColorNormalTexCoord> e;
            e.start_ = v[0];
            e.end_ = v[1];
            edges.push_back(e);
        }
        i++;
    }

    operator delete(cpuData);
}

void TwoManifoldGeometryClipping::process() {
    const Geometry* inputGeometry = inport_.getData();
    tgtAssert(inputGeometry, "no input geometry");

    outport_.setData(0);
    outputGeometry_ = 0;

    const GlMeshGeometryBase* tmgb = dynamic_cast<const GlMeshGeometryBase*>(inport_.getData());

    //check if clipping is enabled, input geometry and shader program are valid
    if (!enableClipping_.get()) {
        outport_.setData(inputGeometry, false);
        if (alwaysFreeMemory_.get())
            freeBufferMemory();
        return;
    }
    else if (!tmgb) {
        LWARNING("Input geometry is not of type GlMeshGeometry, bypassing clipping.");
        outport_.setData(inputGeometry, false);
        if (alwaysFreeMemory_.get())
            freeBufferMemory();
        return;
    }
    else if (!(tmgb->getPrimitiveType() == GL_TRIANGLES || tmgb->getPrimitiveType() == GL_TRIANGLE_STRIP || tmgb->getPrimitiveType() == GL_TRIANGLE_FAN)) {
        LWARNING("Input geometry has to be of type GL_TRIANGLES, GL_TRIANGLE_STRIP, or GL_TRIANGLE_FAN, bypassing clipping.");
        outport_.setData(inputGeometry, false);
        if (alwaysFreeMemory_.get())
            freeBufferMemory();
        return;
    }
    else if (tmgb->isEmpty()) {
        // empty geometry -> if storage is allocated by the buffers, free it
        if (alwaysFreeMemory_.get())
            freeBufferMemory();
        outport_.setData(inputGeometry, false);
        return;
    }

    //check if the vertex layout has changed and update shader defines
    if (tmgb->getVertexLayout() != currentVertexLayout_) {
        currentVertexLayout_ = tmgb->getVertexLayout();
        updateShader(tmgb);
        updateSizePerVertex(tmgb);
        //set numInputTriangles_ to 0 to indicate that buffer storage has to be re-allocated
        numInputTriangles_ = 0;
    }

    // estimate the number of triangles in the geometry (TODO: compute precisely)
    size_t numVertices = 0;
    if (tmgb->usesIndexedDrawing())
        numVertices = tmgb->getNumIndices();
    else
        numVertices = tmgb->getNumVertices();
    size_t numTriangles = 0;
    if (tmgb->getPrimitiveType() == GL_TRIANGLES) {
        numTriangles = numVertices / 3;
    }
    else if (tmgb->getPrimitiveType() == GL_TRIANGLE_STRIP || tmgb->getPrimitiveType() == GL_TRIANGLE_FAN) {
        numTriangles = numVertices - 2;
    }

    //get the number of triangles in the input geometry
    if (numInputTriangles_ != numTriangles) { 
        numInputTriangles_ = numTriangles;
        setupTransformFeedbackSizes();
    }

    // activate program
    program_->activate();
    program_->setIgnoreUniformLocationError(true);

    //disable rasterization and activate transform feedback
    glEnable(GL_RASTERIZER_DISCARD);
    glBindTransformFeedback(GL_TRANSFORM_FEEDBACK, tfo_);
    glBeginTransformFeedback(GL_POINTS);

    glBeginQueryIndexed(GL_TRANSFORM_FEEDBACK_PRIMITIVES_WRITTEN, 0, clipQuery_);
    glBeginQueryIndexed(GL_TRANSFORM_FEEDBACK_PRIMITIVES_WRITTEN, 1, edgeQuery_);
   
    //get clipping plane equation and compute 2D transformation 
    //vec4 plane = vec4(normalize(planeNormal_.get()), planePosition_.get());i
    
    //get clipping plane and transform it to model space
    tgt::plane plane = tgt::plane(normalize(planeNormal_.get()), planePosition_.get());
    plane = plane.transform(tmgb->getInvertedTransformationMatrix());

    //find a transformation to rotate the edges from world coordinates into the xy-plane
    /*mat4 rotationMatrix;
    float angle = acos(dot(-plane.n, vec3(0.f,0.f,1.f)));
    if (std::abs(angle) < 0.0005f) {
        rotationMatrix = mat4(1.f);
    }
    else {
        vec3 rotationAxis;
        if (std::abs(std::abs(angle) - tgt::PIf) < 0.0005f)
            rotationAxis = vec3(1.f,0.f,0.f);
        else
           rotationAxis = cross(-plane.n, vec3(0.f,0.f,1.f));
        
        rotationMatrix = mat4::createRotation(angle, rotationAxis);  
    }*/

    // find the best plane to rotate to
    vec3 destinationNormal;
    ivec2 planeIndices;
    size_t elementIndex = tgt::maxElem(tgt::abs(plane.n));
    switch(elementIndex) {
        case 0: {
                    destinationNormal = vec3(-1.f, 0.f,0.f);
                    planeIndices = ivec2(2,1);
                    break;
                }
        case 1: {
                    destinationNormal = vec3(0.f, 1.f, 0.f);
                    planeIndices = ivec2(2,0);
                    break;
                }
        default: {
                     destinationNormal = vec3(0.f, 0.f, 1.f);
                     planeIndices = ivec2(0,1);
                 }
    }

    mat4 rotationMatrix;
    float angle = acos(dot(-plane.n, destinationNormal));
    if (std::abs(angle) < 0.0005f) {
        rotationMatrix = mat4(1.f);
    }
    else {
        vec3 rotationAxis;
        if (std::abs(std::abs(angle) - tgt::PIf) < 0.0005f) {
            if (elementIndex == 0)
                rotationAxis = vec3(0.f,0.f,1.f);
            else 
                rotationAxis = vec3(1.f,0.f,0.f);
        }
        else
           rotationAxis = cross(-plane.n, destinationNormal);
        
        rotationMatrix = mat4::createRotation(angle, rotationAxis);  
    }
    
    
    //set uniforms
    //program_->setUniform("modelMatrix", inport_.getData()->getTransformationMatrix());
    program_->setUniform("plane", plane.toVec4());
    program_->setUniform("xyTransform", rotationMatrix);
    program_->setUniform("xyIndices", planeIndices);

    tgt::Stopwatch bmTimer;
    bmTimer.start();

    LGL_ERROR;
    tmgb->render();
    LGL_ERROR;

    glEndTransformFeedback();

    glEndQueryIndexed(GL_TRANSFORM_FEEDBACK_PRIMITIVES_WRITTEN, 0);
    glEndQueryIndexed(GL_TRANSFORM_FEEDBACK_PRIMITIVES_WRITTEN, 1);

    glGetQueryObjectuiv(clipQuery_, GL_QUERY_RESULT, &writtenMeshVertices_);        
    glGetQueryObjectuiv(edgeQuery_, GL_QUERY_RESULT, &writtenEdgeVertices_);

    bmTimer.stop();
    LDEBUG("Clipping time: " << bmTimer.getRuntime() << " ms");
    bmTimer.reset();

    //clean up OpenGL status
    glBindTransformFeedback(GL_TRANSFORM_FEEDBACK, 0);
    glDisable(GL_RASTERIZER_DISCARD);
    program_->deactivate();

    LDEBUG("Recorded vertices in clipped input mesh: " << writtenMeshVertices_);
    LDEBUG("Recorded vertices in clip edge buffer: " << writtenEdgeVertices_);

    //if nothing has been clipped: just copy the input
    /*if (writtenMeshVertices_ > 0 && writtenEdgeVertices_ == 0) {
        outport_.setData(inputGeometry, false);      
        if (alwaysFreeMemory_.get())
            freeBufferMemory();
        return;
    }*/

    //compute the model space plane normal for the vertex normals if necessary
    if (tmgb->supportsNormals()) {
        mat3 tmpMat = inputGeometry->getInvertedTransformationMatrix().getRotationalPartMat3();
        mat3 normalTransform;
        tmpMat.invert(normalTransform);
        transpose(normalTransform);
        capNormal_ = normalTransform * (-planeNormal_.get()); 
    }

    //build output geometry depending on vertex layout
    switch (currentVertexLayout_) {
        case VertexBase::SIMPLE : 
                    { 
                        outputGeometry_ = buildOutputGeometry<VertexBase>(tmgb); 
                        break;
                    }
        case VertexBase::COLOR : 
                    { 
                        outputGeometry_ = buildOutputGeometry<VertexColor>(tmgb);
                        break; 
                    }
        case VertexBase::NORMAL : 
                    { 
                        outputGeometry_ = buildOutputGeometry<VertexNormal>(tmgb); 
                        break;
                    }
        case VertexBase::COLOR_NORMAL : 
                    { 
                        outputGeometry_ = buildOutputGeometry<VertexColorNormal>(tmgb); 
                        break;
                    }
        case VertexBase::TEXCOORD : 
                    { 
                        outputGeometry_ = buildOutputGeometry<VertexTexCoord>(tmgb); 
                        break;
                    }
        case VertexBase::COLOR_TEXCOORD : 
                    { 
                        outputGeometry_ = buildOutputGeometry<VertexColorTexCoord>(tmgb); 
                        break;
                    }
        case VertexBase::NORMAL_TEXCOORD : 
                    { 
                        outputGeometry_ = buildOutputGeometry<VertexNormalTexCoord>(tmgb); 
                        break;
                    }
        case VertexBase::COLOR_NORMAL_TEXCOORD: 
                    { 
                        outputGeometry_ = buildOutputGeometry<VertexColorNormalTexCoord>(tmgb); 
                        break;
                    }
        default: break;
    }

    outport_.setData(outputGeometry_);
    
    if (alwaysFreeMemory_.get())
        freeBufferMemory();
}

}   // namespace
