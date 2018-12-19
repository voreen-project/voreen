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

#include "geometrysource.h"

#include "voreen/core/voreenapplication.h"
#include "voreen/core/datastructures/geometry/meshlistgeometry.h"
#include "voreen/core/datastructures/geometry/pointlistgeometry.h"
#include "voreen/core/datastructures/geometry/pointsegmentlistgeometry.h"
#include "voreen/core/datastructures/geometry/trianglemeshgeometryindexed.h"
#include "voreen/core/datastructures/geometry/glmeshgeometry.h"
#include "voreen/core/datastructures/meta/serializablevectormetadata.h"

#include "voreen/core/datastructures/callback/memberfunctioncallback.h"
#include "tgt/filesystem.h"
#include "tgt/exception.h"
#include "ext/tinyobj/tiny_obj_loader.h"

#include <algorithm>
#include <vector>
#include <fstream>

using tgt::vec3;
using tgt::ivec3;
using tgt::ivec2;
using std::vector;
using std::string;

namespace voreen {

const std::string GeometrySource::loggerCat_("voreen.core.GeometrySource");

GeometrySource::GeometrySource()
    : Processor()
    , geometryFile_("geometryFile", "Geometry File", "Open Geometry File", VoreenApplication::app()->getUserDataPath(), "Geometry (*.vge *.ply *.obj *.stl);;Voreen Geometry (*.vge);;PLY mesh (*.ply);;OBJ mesh (*.obj);;STL mesh (*.stl);;ASCII point lists (*.txt)")
    , geometryType_("geometryType", "Geometry Type")
    , skipItemCount_("skipItems", "Items to skip after each point", 0, 0, 100)
    , loadGeometry_("loadGeometry", "Load Geometry")
    , clearGeometry_("clearGeometry", "Clear Geometry")
    //, useIndexedGeometry_("useIndexed", "Prefer indexed geometry for external formats", false)
    , calculateNormals_("calcNormals", "Calculate missing normals for external formats", false)
    , addColorToOBJ_("addobjcolor", "Add Color to OBJ", false)
    , objColor_("objcolor", "OBJ color", tgt::vec4::one)
    , outport_(Port::OUTPORT, "geometry.pointlist", "PointList Output", false)
    , forceReload_(false)
{
    addPort(outport_);

    geometryType_.addOption("geometry", "Voreen Geometry (.vge)");
    geometryType_.addOption("pointlist", "Pointlist");
    geometryType_.addOption("segmentlist", "Segmented Pointlist");

    geometryFile_.onChange(MemberFunctionCallback<GeometrySource>(this, &GeometrySource::forceReload));
    geometryType_.onChange(MemberFunctionCallback<GeometrySource>(this, &GeometrySource::updatePropertyVisibility));
    loadGeometry_.onChange(MemberFunctionCallback<GeometrySource>(this, &GeometrySource::forceReload));
    clearGeometry_.onChange(MemberFunctionCallback<GeometrySource>(this, &GeometrySource::clearGeometry));
    //useIndexedGeometry_.onChange(MemberFunctionCallback<GeometrySource>(this, &GeometrySource::forceReload));

    addProperty(geometryFile_);
    addProperty(geometryType_);
    addProperty(skipItemCount_);
    addProperty(loadGeometry_);
    addProperty(clearGeometry_);

    //addProperty(useIndexedGeometry_);
    addProperty(calculateNormals_);

    addProperty(addColorToOBJ_);
    addProperty(objColor_);
}

Processor* GeometrySource::create() const {
    return new GeometrySource();
}

void GeometrySource::process() {
    if (geometryFile_.get() != "" && forceReload_) {
        try {
            readGeometry();
        }
        catch (tgt::FileNotFoundException& f) {
            LERROR(f.what());
        }
        forceReload_ = false;
        updatePropertyVisibility();
    }
}

void GeometrySource::initialize() {
    Processor::initialize();
    forceReload_ = true;
}

void GeometrySource::readGeometry() {
    std::string filename = geometryFile_.get();
    setProgress(0.f);

    if (geometryFile_.get() == "")
        return;

    if (!tgt::FileSystem::fileExists(geometryFile_.get()))
        throw tgt::FileNotFoundException("File does not exists", geometryFile_.get());

    if (geometryType_.isSelected("geometry")) {
        LINFO("Reading geometry file: " << geometryFile_.get());
        if(endsWith(filename, ".ply")) {
            LINFO("Reading PLY file.");
            try {
                Geometry* geometry = readPLYGeometry(filename);
                tgtAssert(geometry, "null pointer returned (exception expected)");
                outport_.setData(geometry);
                setProgress(1.f);
            }
            catch (VoreenException& e) {
                LERROR(e.what());
                setProgress(0.f);
            }
            catch (tgt::CorruptedFileException& e) {
                LERROR(e.what());
                setProgress(0.f);
            }
        }
        else if(endsWith(filename, ".obj")) {
            LINFO("Reading Wavefront obj file.");
            try {
                Geometry *geometry = 0;
                if (addColorToOBJ_.get())
                    geometry = readOBJGeometryWithColor(filename, objColor_.get());
                else
                    geometry = readOBJGeometry(filename);
                tgtAssert(geometry, "null pointer returned (exception expected)");
                outport_.setData(geometry);
                setProgress(1.f);
            }
            catch (VoreenException &e) {
                LERROR(e.what());
                setProgress(0.f);
            }
            catch (tgt::CorruptedFileException &e) {
                LERROR(e.what());
                setProgress(0.f);
            }
        }
        else if(endsWith(filename, ".stl")) {
                LINFO("Reading STL file.");
                try {
                    Geometry* geometry = readSTLGeometry(filename);
                    tgtAssert(geometry, "null pointer returned (exception expected)");
                    outport_.setData(geometry);
                    setProgress(1.f);
                }
                catch (VoreenException& e) {
                    LERROR(e.what());
                    setProgress(0.f);
                }
                catch (tgt::CorruptedFileException& e) {
                    LERROR(e.what());
                    setProgress(0.f);
                }
        } else {
            try {
                Geometry* geometry = readVoreenGeometry(geometryFile_.get());
                tgtAssert(geometry, "null pointer returned (exception expected)");
                if(TriangleMeshGeometryBase* tmgb = dynamic_cast<TriangleMeshGeometryBase*>(geometry))
                    tmgb->loadTextureData();
                else if (GlMeshGeometryBase* glmgb = dynamic_cast<GlMeshGeometryBase*>(geometry))
                    glmgb->loadTextureData();
                outport_.setData(geometry);
                setProgress(1.f);
            }
            catch (VoreenException& e) {
                LERROR(e.what());
                setProgress(0.f);
            }
        }
    }
    else if (geometryType_.isSelected("pointlist") || geometryType_.isSelected("segmentlist")) {
        PointListType listType;
        if (geometryType_.isSelected("pointlist")) {
            LINFO("Reading point list from file: " << geometryFile_.get());
            listType = PointList;
        }
        else {
            LINFO("Reading segmented point list from file: " << geometryFile_.get());
            listType = SegmentedPointList;
        }

        try {
            Geometry* geometry = readPointList(geometryFile_.get(), listType, skipItemCount_.get());
            tgtAssert(geometry, "null pointer returned (exception expected)");
            outport_.setData(geometry);
            setProgress(1.f);
        }
        catch (VoreenException& e) {
            LERROR(e.what());
            setProgress(0.f);
        }
    }
    else {
        LWARNING("Unknown geometry type: " << geometryType_.get());
    }

    updatePropertyVisibility();
}

Geometry* GeometrySource::readPLYGeometry(const std::string& filename) {
    // read PLY file (very incomplete reader!)
    tgt::File* file = FileSys.open(filename);
    std::string line = file->getLine();
    if(line != "ply")
        throw tgt::CorruptedFileException("Expected 'ply' at first", filename);

    enum {
        ASCII,
        BINARY_BIG_ENDIAN,
        BINARY_LITTLE_ENDIAN
    };

    // read header:
    line = trim(file->getLine());
    int numVertices = -1;
    int numFaces = -1;
    int format = -1;

    while(line != "end_header") {
        vector<std::string> expl = strSplit(line, ' ');
        if(expl.size() > 0) {
            string command = expl[0];

            if(command == "format") {
                if(expl.size() != 3)
                    throw tgt::CorruptedFileException("Could not parse format", filename);
                else {
                    if(expl[2] != "1.0")
                        throw tgt::CorruptedFileException("Unknown format: " + line, filename);

                    if(expl[1] == "ascii")
                        format = ASCII;
                    else if(expl[1] == "binary_big_endian")
                        format = BINARY_BIG_ENDIAN;
                    else if(expl[1] == "binary_little_endian")
                        format = BINARY_LITTLE_ENDIAN;
                    else
                        throw tgt::CorruptedFileException("Unknown format: " + line, filename);
                }
            }
            else if(command == "comment") {
                LINFO(line);
            }
            else if(command == "element") {
                if(expl.size() != 3)
                    throw tgt::CorruptedFileException("Could not parse element" + line, filename);
                else {
                    string type = expl[1];
                    if(type == "vertex") {
                        numVertices = stoi(expl[2]);
                        LINFO(numVertices << " vertices");
                    }
                    else if(type == "face") {
                        numFaces = stoi(expl[2]);
                        LINFO(numFaces << " faces");
                    }
                }
            }
            else if(command == "property") {
            }
        }

        line = trim(file->getLine());
    }

    // read data:
    vector<VertexNormal> vertices;
    if(format == ASCII) {
        for(int i=0; i<numVertices; i++) {
            line = trim(file->getLine());
            vector<std::string> expl = strSplit(line, ' ');

            if(expl.size() >= 3) {
                vec3 pos;
                pos.x = stof(expl[0]);
                pos.y = stof(expl[1]);
                pos.z = stof(expl[2]);
                vertices.push_back(VertexNormal(pos, vec3(1.0f)));
                //vertices.push_back(VertexNormal(pos, vec3(stof(expl[4]))));
            }
            setProgress(0.5f * static_cast<float>(i) / numVertices);
        }
    } else if(format == BINARY_BIG_ENDIAN || format == BINARY_LITTLE_ENDIAN) {
        std::vector<tgt::vec3> positions;
        positions.resize(numVertices);
        if(!file->read(&(*positions.begin()), numVertices * sizeof(tgt::vec3))) {
            throw tgt::CorruptedFileException("Could not parse binary vertex data", filename);
        }
        vertices.resize(numVertices);
        int counter = 0;

#ifdef VRN_MODULE_OPENMP
    #pragma omp parallel for
#endif
        for(int i = 0; i < numVertices; i++) {
            tgt::vec3 pos = positions[i];
            if(format == BINARY_BIG_ENDIAN) {
                unsigned char *memp = reinterpret_cast<unsigned char*>(&pos);
                std::reverse(memp, memp + sizeof(float));
                std::reverse(memp + sizeof(float), memp + 2 * sizeof(float));
                std::reverse(memp + 2 * sizeof(float), memp + 3 * sizeof(float));
            }
            vertices[i] = (VertexNormal(pos, vec3(1.0f)));
            counter++;
            setProgress(0.5f * static_cast<float>(counter) / numVertices);
        }
    }

    std::vector<tgt::ivec3> faces;
    faces.resize(numFaces);
    std::vector<tgt::vec3> vertexNormals;
    if(calculateNormals_.get())
        vertexNormals.resize(numVertices, tgt::vec3(0.f));

    if(format == BINARY_BIG_ENDIAN || format == BINARY_LITTLE_ENDIAN) {
        int counter = 0;
        const int bytesPerFace = 1 + 3 * sizeof(int);
        size_t numBytes = numFaces * bytesPerFace;
        std::vector<unsigned char> faceBuffer(numBytes);
        file->read(&(*faceBuffer.begin()), numBytes);

//not: do not use OpenMP here as the code might throw an exception that might then be unhandled!
//#ifdef VRN_MODULE_OPENMP
//    #pragma omp parallel for
//#endif
        for(int i=0; i<numFaces; i++) {
            unsigned char* curOffset = &(*faceBuffer.begin()) + i * bytesPerFace;
            unsigned char faceVertexNum = *curOffset;
            if(faceVertexNum != 3)
                throw tgt::CorruptedFileException("Faces with number of vertices != 3 currently not supported", filename);

            int* elms = reinterpret_cast<int*>(curOffset + 1);
            tgt::ivec3 faceIndices(elms[0], elms[1], elms[2]);
            if(format == BINARY_BIG_ENDIAN) {
                unsigned char *memp = reinterpret_cast<unsigned char*>(&faceIndices);
                std::reverse(memp, memp + sizeof(int));
                std::reverse(memp + sizeof(int), memp + 2 * sizeof(int));
                std::reverse(memp + 2 * sizeof(int), memp + 3 * sizeof(int));
            }
            faces[i] = faceIndices;

            if(calculateNormals_.get()) {
                tgt::vec3 norm = normalize(cross(vertices[faceIndices.y].pos_ - vertices[faceIndices.x].pos_, vertices[faceIndices.z].pos_ - vertices[faceIndices.x].pos_));
                vertexNormals[faceIndices.x] += norm;
                vertexNormals[faceIndices.y] += norm;
                vertexNormals[faceIndices.z] += norm;
            }
            counter++;
            setProgress(0.5f + 0.5f * static_cast<float>(counter) / numFaces);
        }
    } else {
        for(int i=0; i<numFaces; i++) {
            tgt::ivec3 faceIndices;
            if(format == ASCII) {
                line = trim(file->getLine());
                vector<std::string> expl = strSplit(line, ' ');

                if(expl.size() == 4 && expl[0] == "3") {
                    faceIndices.x = stoi(expl[1]);
                    faceIndices.y = stoi(expl[2]);
                    faceIndices.z = stoi(expl[3]);
                }
            }

            faces[i] = faceIndices;

            if(calculateNormals_.get()) {
                tgt::vec3 norm = normalize(cross(vertices[faceIndices.y].pos_ - vertices[faceIndices.x].pos_, vertices[faceIndices.z].pos_ - vertices[faceIndices.x].pos_));
                vertexNormals[faceIndices.x] += norm;
                vertexNormals[faceIndices.y] += norm;
                vertexNormals[faceIndices.z] += norm;
            }
            setProgress(0.5f + 0.5f * static_cast<float>(i) / numFaces);
        }
    }

    if(calculateNormals_.get()) {
#ifdef VRN_MODULE_OPENMP
    #pragma omp parallel for
#endif
        for(int i = 0; i < numVertices; i++)
            vertices[i].normal_ = normalize(vertexNormals.at(i));
    }

    GlMeshGeometryUInt32Normal* mesh = new GlMeshGeometryUInt32Normal();
    mesh->setVertices(vertices);
    for(size_t i = 0; i < static_cast<size_t>(numFaces); ++i)
        for (size_t v = 0; v < 3; ++v)
            mesh->addIndex(faces.at(i)[v]);

    LINFO("Read " << vertices.size() << " vertices and " << numFaces << " triangles " << mesh->getBoundingBox());

    file->close();
    delete file;
    return mesh;
}

Geometry* GeometrySource::readOBJGeometry(const std::string& filename) {
    std::vector<tinyobj::shape_t> shapes;
    std::string basePath = tgt::FileSystem::dirName(filename) + "/";

    // read obj file
    std::string err = tinyobj::LoadObj(shapes, filename.c_str(), basePath.c_str());
    if(!err.empty())
        throw tgt::CorruptedFileException(std::string("Failed to read obj: ") + err, filename);

    std::vector<VertexNormalTexCoord> vertices;
    std::vector<tgt::ivec3> faces;

    int baseIndex = 0;
    std::vector<std::string> textureNames;

    int textureCounter = 0;
    std::map<std::string, int> texCounterMap;

    for (size_t i = 0; i < shapes.size(); i++) {

        unsigned short requiredTextures = 0;
        std::map<unsigned short, unsigned char> textureFunctions;

        std::string texNameA = trim(shapes[i].material.ambient_texname, " \n\r\t");
        if(!texNameA.empty()){
            bool validTexture = false;
            if(texCounterMap.count(texNameA) == 0) {
                if(textureCounter < 16) {
                    validTexture = true;
                    texCounterMap[texNameA] = textureCounter++;
                    textureFunctions[texCounterMap[texNameA]] = 0;
                    textureNames.push_back(tgt::FileSystem::cleanupPath(basePath + texNameA));
                }
            } else
                validTexture = true;

            if(validTexture) {
                requiredTextures |= 1 << texCounterMap[texNameA];
                textureFunctions[texCounterMap[texNameA]] |= 1;
            }
        }

        std::string texNameD = trim(shapes[i].material.diffuse_texname, " \n\r\t");
        if(!texNameD.empty()){
            bool validTexture = false;
            if(texCounterMap.count(texNameD) == 0) {
                if(textureCounter < 16) {
                    validTexture = true;
                    texCounterMap[texNameD] = textureCounter++;
                    textureFunctions[texCounterMap[texNameD]] = 0;
                    textureNames.push_back(tgt::FileSystem::cleanupPath(basePath + texNameD));
                }
            } else
                validTexture = true;

            if(validTexture) {
                requiredTextures |= 1 << texCounterMap[texNameD];
                textureFunctions[texCounterMap[texNameD]] |= 2;
            }
        }

        std::string texNameS = trim(shapes[i].material.specular_texname, " \n\r\t");
        if(!texNameS.empty()){
            bool validTexture = false;
            if(texCounterMap.count(texNameS) == 0) {
                if(textureCounter < 16) {
                    validTexture = true;
                    texCounterMap[texNameS] = textureCounter++;
                    textureFunctions[texCounterMap[texNameS]] = 0;
                    textureNames.push_back(tgt::FileSystem::cleanupPath(basePath + texNameS));
                }
            } else
                validTexture = true;

            if(validTexture) {
                requiredTextures |= 1 << texCounterMap[texNameS];
                textureFunctions[texCounterMap[texNameS]] |= 4;
            }
        }

        bool bumpMapFound = false;
        std::string texNameB = trim(shapes[i].material.bump_texname, " \n\r\t");
        if(!texNameB.empty()){
            bool validTexture = false;
            if(texCounterMap.count(texNameB) == 0) {
                if(textureCounter < 16) {
                    validTexture = true;
                    texCounterMap[texNameB] = textureCounter++;
                    textureFunctions[texCounterMap[texNameB]] = 0;
                    textureNames.push_back(tgt::FileSystem::cleanupPath(basePath + texNameB));
                }
            } else
                validTexture = true;

            if(validTexture) {
                requiredTextures |= 1 << texCounterMap[texNameB];
                textureFunctions[texCounterMap[texNameB]] |= 8;
                bumpMapFound = true;
            }
        }

        std::string texNameN = trim(shapes[i].material.normal_texname, " \n\r\t");
        if(!bumpMapFound && !texNameN.empty()){
            bool validTexture = false;
            if(texCounterMap.count(texNameN) == 0) {
                if(textureCounter < 16) {
                    validTexture = true;
                    texCounterMap[texNameN] = textureCounter++;
                    textureFunctions[texCounterMap[texNameN]] = 0;
                    textureNames.push_back(tgt::FileSystem::cleanupPath(basePath + texNameN));
                }
            } else
                validTexture = true;

            if(validTexture) {
                requiredTextures |= 1 << texCounterMap[texNameN];
                //textureFunctions |= 8;
                textureFunctions[texCounterMap[texNameN]] |= 8;
            }
        }

        bool useTextures = false;
        if(requiredTextures > 0 && !shapes[i].mesh.texcoords.empty())
            useTextures = true;

        unsigned short actualTextureFunctions = 0;
        int useCounter = 0;
        for(std::map<unsigned short, unsigned char>::const_iterator it = textureFunctions.begin(); it != textureFunctions.end(); it++) {
            actualTextureFunctions |=  it->second << (useCounter * 4);
            useCounter++;
        }

        bool hasNormals   = !shapes[i].mesh.normals.empty();
        for (size_t v = 0; v < shapes[i].mesh.positions.size() / 3; v++) {
            vertices.push_back(VertexNormalTexCoord(tgt::vec3::fromPointer(&shapes[i].mesh.positions[3*v]),
                                              hasNormals  ? tgt::vec3::fromPointer(&(shapes[i].mesh.normals[3*v])) : tgt::vec3(0.0),
                                              useTextures ? tgt::vec2::fromPointer(&(shapes[i].mesh.texcoords[2*v])) : tgt::vec2(0.f),
                                              useTextures ? tgt::ivec2(requiredTextures, actualTextureFunctions) : tgt::ivec2(0)
                                              )
                              );
        }

        for (size_t v = 0; v < shapes[i].mesh.indices.size() / 3; v++) {
            tgt::ivec3 curFace = tgt::ivec3(tgt::Vector3<unsigned int>::fromPointer(&shapes[i].mesh.indices[3*v])) + baseIndex;
            faces.push_back(curFace);
        }

        if(!hasNormals && calculateNormals_.get()) {
            std::vector<tgt::vec3> vertexNormals = std::vector<tgt::vec3>(vertices.size(), tgt::vec3(0.f));
            for(size_t i=0; i < faces.size(); i++) {
                tgt::ivec3 faceIndices = faces.at(i);
                tgt::vec3 norm = normalize(cross(vertices[faceIndices.y].pos_ - vertices[faceIndices.x].pos_, vertices[faceIndices.z].pos_ - vertices[faceIndices.x].pos_));
                vertexNormals[faceIndices.x] += norm;
                vertexNormals[faceIndices.y] += norm;
                vertexNormals[faceIndices.z] += norm;
            }
            for(size_t i = 0; i < vertexNormals.size(); i++)
                vertices.at(i).normal_ = normalize(vertexNormals.at(i));
        }

        baseIndex = static_cast<int>(vertices.size());
    }

    GlMeshGeometryUInt32NormalTexCoord* mesh = new GlMeshGeometryUInt32NormalTexCoord();
    mesh->setVertices(vertices);
    for(size_t i = 0; i < faces.size(); ++i)
        for (size_t v = 0; v < 3; ++v)
            mesh->addIndex(faces.at(i)[v]);

    LINFO("Model requires the following textures:");
    //for(std::map<std::string, int>::const_iterator it = texCounterMap.begin(); it != texCounterMap.end(); it++) {
    for(size_t i = 0; i < textureNames.size(); i++)
        LINFO("    " + textureNames.at(i));

    if(!textureNames.empty())
        mesh->loadTextureData(textureNames);

    LINFO("Successfully read " << shapes.size() << " shapes with a total of " << vertices.size() << " vertices and " << faces.size() << " triangles.");

    return mesh;
}

Geometry* GeometrySource::readOBJGeometryWithColor(const std::string& filename, tgt::vec4 color) {
    std::vector<tinyobj::shape_t> shapes;
    std::string basePath = tgt::FileSystem::dirName(filename) + "/";

    // read obj file
    std::string err = tinyobj::LoadObj(shapes, filename.c_str(), basePath.c_str());
    if(!err.empty())
        throw tgt::CorruptedFileException(std::string("Failed to read obj: ") + err, filename);

    std::vector<VertexColorNormalTexCoord> vertices;
    std::vector<tgt::ivec3> faces;

    int baseIndex = 0;
    std::vector<std::string> textureNames;

    int textureCounter = 0;
    std::map<std::string, int> texCounterMap;

    for (size_t i = 0; i < shapes.size(); i++) {

        unsigned short requiredTextures = 0;
        std::map<unsigned short, unsigned char> textureFunctions;

        std::string texNameA = trim(shapes[i].material.ambient_texname, " \n\r\t");
        if(!texNameA.empty()){
            bool validTexture = false;
            if(texCounterMap.count(texNameA) == 0) {
                if(textureCounter < 16) {
                    validTexture = true;
                    texCounterMap[texNameA] = textureCounter++;
                    textureFunctions[texCounterMap[texNameA]] = 0;
                    textureNames.push_back(tgt::FileSystem::cleanupPath(basePath + texNameA));
                }
            } else
                validTexture = true;

            if(validTexture) {
                requiredTextures |= 1 << texCounterMap[texNameA];
                textureFunctions[texCounterMap[texNameA]] |= 1;
            }
        }

        std::string texNameD = trim(shapes[i].material.diffuse_texname, " \n\r\t");
        if(!texNameD.empty()){
            bool validTexture = false;
            if(texCounterMap.count(texNameD) == 0) {
                if(textureCounter < 16) {
                    validTexture = true;
                    texCounterMap[texNameD] = textureCounter++;
                    textureFunctions[texCounterMap[texNameD]] = 0;
                    textureNames.push_back(tgt::FileSystem::cleanupPath(basePath + texNameD));
                }
            } else
                validTexture = true;

            if(validTexture) {
                requiredTextures |= 1 << texCounterMap[texNameD];
                textureFunctions[texCounterMap[texNameD]] |= 2;
            }
        }

        std::string texNameS = trim(shapes[i].material.specular_texname, " \n\r\t");
        if(!texNameS.empty()){
            bool validTexture = false;
            if(texCounterMap.count(texNameS) == 0) {
                if(textureCounter < 16) {
                    validTexture = true;
                    texCounterMap[texNameS] = textureCounter++;
                    textureFunctions[texCounterMap[texNameS]] = 0;
                    textureNames.push_back(tgt::FileSystem::cleanupPath(basePath + texNameS));
                }
            } else
                validTexture = true;

            if(validTexture) {
                requiredTextures |= 1 << texCounterMap[texNameS];
                textureFunctions[texCounterMap[texNameS]] |= 4;
            }
        }

        bool bumpMapFound = false;
        std::string texNameB = trim(shapes[i].material.bump_texname, " \n\r\t");
        if(!texNameB.empty()){
            bool validTexture = false;
            if(texCounterMap.count(texNameB) == 0) {
                if(textureCounter < 16) {
                    validTexture = true;
                    texCounterMap[texNameB] = textureCounter++;
                    textureFunctions[texCounterMap[texNameB]] = 0;
                    textureNames.push_back(tgt::FileSystem::cleanupPath(basePath + texNameB));
                }
            } else
                validTexture = true;

            if(validTexture) {
                requiredTextures |= 1 << texCounterMap[texNameB];
                textureFunctions[texCounterMap[texNameB]] |= 8;
                bumpMapFound = true;
            }
        }

        std::string texNameN = trim(shapes[i].material.normal_texname, " \n\r\t");
        if(!bumpMapFound && !texNameN.empty()){
            bool validTexture = false;
            if(texCounterMap.count(texNameN) == 0) {
                if(textureCounter < 16) {
                    validTexture = true;
                    texCounterMap[texNameN] = textureCounter++;
                    textureFunctions[texCounterMap[texNameN]] = 0;
                    textureNames.push_back(tgt::FileSystem::cleanupPath(basePath + texNameN));
                }
            } else
                validTexture = true;

            if(validTexture) {
                requiredTextures |= 1 << texCounterMap[texNameN];
                //textureFunctions |= 8;
                textureFunctions[texCounterMap[texNameN]] |= 8;
            }
        }

        bool useTextures = false;
        if(requiredTextures > 0 && !shapes[i].mesh.texcoords.empty())
            useTextures = true;

        unsigned short actualTextureFunctions = 0;
        int useCounter = 0;
        for(std::map<unsigned short, unsigned char>::const_iterator it = textureFunctions.begin(); it != textureFunctions.end(); it++) {
            actualTextureFunctions |=  it->second << (useCounter * 4);
            useCounter++;
        }

        bool hasNormals   = !shapes[i].mesh.normals.empty();
        for (size_t v = 0; v < shapes[i].mesh.positions.size() / 3; v++) {
            vertices.push_back(VertexColorNormalTexCoord(tgt::vec3::fromPointer(&shapes[i].mesh.positions[3*v]),
                                              objColor_.get(),
                                              hasNormals  ? tgt::vec3::fromPointer(&(shapes[i].mesh.normals[3*v])) : tgt::vec3(0.0),
                                              useTextures ? tgt::vec2::fromPointer(&(shapes[i].mesh.texcoords[2*v])) : tgt::vec2(0.f),
                                              useTextures ? tgt::ivec2(requiredTextures, actualTextureFunctions) : tgt::ivec2(0)
                                              )
                              );
        }

        for (size_t v = 0; v < shapes[i].mesh.indices.size() / 3; v++) {
            tgt::ivec3 curFace = tgt::ivec3(tgt::Vector3<unsigned int>::fromPointer(&shapes[i].mesh.indices[3*v])) + baseIndex;
            faces.push_back(curFace);
        }

        if(!hasNormals && calculateNormals_.get()) {
            std::vector<tgt::vec3> vertexNormals = std::vector<tgt::vec3>(vertices.size(), tgt::vec3(0.f));
            for(size_t i=0; i < faces.size(); i++) {
                tgt::ivec3 faceIndices = faces.at(i);
                tgt::vec3 norm = normalize(cross(vertices[faceIndices.y].pos_ - vertices[faceIndices.x].pos_, vertices[faceIndices.z].pos_ - vertices[faceIndices.x].pos_));
                vertexNormals[faceIndices.x] += norm;
                vertexNormals[faceIndices.y] += norm;
                vertexNormals[faceIndices.z] += norm;
            }
            for(size_t i = 0; i < vertexNormals.size(); i++)
                vertices.at(i).normal_ = normalize(vertexNormals.at(i));
        }

        baseIndex = static_cast<int>(vertices.size());
    }

    GlMeshGeometryUInt32ColorNormalTexCoord* mesh = new GlMeshGeometryUInt32ColorNormalTexCoord();
    mesh->setVertices(vertices);
    for(size_t i = 0; i < faces.size(); ++i)
        for (size_t v = 0; v < 3; ++v)
            mesh->addIndex(faces.at(i)[v]);

    LINFO("Model requires the following textures:");
    //for(std::map<std::string, int>::const_iterator it = texCounterMap.begin(); it != texCounterMap.end(); it++) {
    for(size_t i = 0; i < textureNames.size(); i++)
        LINFO("    " + textureNames.at(i));

    if(!textureNames.empty())
        mesh->loadTextureData(textureNames);

    LINFO("Successfully read " << shapes.size() << " shapes with a total of " << vertices.size() << " vertices and " << faces.size() << " triangles.");

    return mesh;
}

Geometry* GeometrySource::readSTLGeometry(const std::string& filename) {

    std::ifstream f(filename.c_str(), std::ios::in);
    if (!f.good()) {
        throw std::runtime_error("STL File not valid.");
    }

    std::unique_ptr<GlMeshGeometryUInt32Normal> mesh(new GlMeshGeometryUInt32Normal());

    char buf[6];
    buf[5] = 0;
    f.read(buf, 5);
    const std::string asciiHeader = "solid";
    if (std::string(buf) == asciiHeader) {
        f.seekg(0, std::ios::beg);
        if (f.good()) {
            std::string s0, s1;
            while (!f.eof()) {
                f >> s0;
                if (s0 == "facet") {
                    tgt::vec3 normal;
                    GlMeshGeometryUInt16Normal::VertexType v1, v2, v3;
                    f >> s1 >> normal.x >> normal.y >> normal.z;
                    f >> s0 >> s1;
                    f >> s0 >> v1.pos_.x >> v1.pos_.y >> v1.pos_.z;
                    f >> s0 >> v2.pos_.x >> v2.pos_.y >> v2.pos_.z;
                    f >> s0 >> v3.pos_.x >> v3.pos_.y >> v3.pos_.z;
                    f >> s0;
                    f >> s0;

                    // TODO: calculate smooth normals.
                    v1.normal_ = v2.normal_ = v3.normal_ = normal;

                    mesh->addVertex(v1);
                    mesh->addVertex(v2);
                    mesh->addVertex(v3);

                } else if (s0 == "endsolid") {
                    break;
                }
            }
        }
    } else {
        f.close();
        f.open(filename.c_str(), std::ios::in | std::ios::binary);
        char comment[80];
        f.read(comment, 80);

        if (!f.good()) {
            throw tgt::CorruptedFileException("STL File not valid.");
        }

        comment[79] = 0;
        int32_t nFacets;
        f.read(reinterpret_cast<char *>(&nFacets), sizeof(int32_t));

        if (!f.good()) {
            throw tgt::CorruptedFileException("STL File not valid.");
        }

        float v[12];
        unsigned short uint16;
        for (int32_t i = 0; i < nFacets; ++i) {
            for (unsigned int j = 0; j < 12; ++j) {
                f.read(reinterpret_cast<char *>(&v[j]), sizeof(float));
            }
            f.read(reinterpret_cast<char *>(&uint16), sizeof(unsigned short));

            tgt::vec3 normal;
            GlMeshGeometryUInt16Normal::VertexType v1, v2, v3;

            normal.x = v[0];
            normal.y = v[1];
            normal.z = v[2];
            v1.pos_.x = v[3];
            v1.pos_.y = v[4];
            v1.pos_.z = v[5];
            v2.pos_.x = v[6];
            v2.pos_.y = v[7];
            v2.pos_.z = v[8];
            v3.pos_.x = v[9];
            v3.pos_.y = v[10];
            v3.pos_.z = v[11];

            // TODO: calculate smooth normals.
            v1.normal_ = v2.normal_ = v3.normal_ = normal;

            mesh->addVertex(v1);
            mesh->addVertex(v2);
            mesh->addVertex(v3);
        }
    }
    f.close();

    return mesh.release();
}

Geometry* GeometrySource::readVoreenGeometry(const std::string& filename) {
    // read Voreen geometry serialization (.vge)
    std::ifstream stream;
    stream.open(filename.c_str(), std::ios_base::in);
    if (stream.fail())
        throw VoreenException("Failed to open file " + geometryFile_.get() + " for reading");

    XmlDeserializer deserializer;
    try {
        deserializer.read(stream);
        Geometry* geometry = 0;
        deserializer.deserialize("Geometry", geometry);
        return geometry;
    }
    catch (SerializationException &e) {
        throw VoreenException("Failed to deserialize Voreen Geometry from " + geometryFile_.get() + ": " + e.what());
    }
}

Geometry* GeometrySource::readPointList(const std::string& filename, PointListType listType, int skipItems) {
    tgtAssert(skipItems >= 0, "skipItems must be non-negative");

    std::ifstream inFile;
    inFile.open(filename.c_str());
    if (inFile.fail())
        throw VoreenException("Failed to open file for reading: " + filename);

    PointListGeometryVec3* pointListGeometry = 0;

    PointSegmentListGeometry<tgt::vec3>* pointSegmentListGeometry = 0;
    std::vector<std::vector<tgt::vec3> > segments;

    if (listType == PointList)
        pointListGeometry = new PointListGeometryVec3();
    else if (listType == SegmentedPointList)
        pointSegmentListGeometry = new PointSegmentListGeometryVec3();
    else {
        tgtAssert(false, "unknown point list type");
        throw VoreenException("Unknown point list type: " + genericToString(listType));
    }

    int lastSegID = -1;
    while (true) {

        tgt::vec3 point;
        float dummy;
        int segID;

        if (!(inFile >> point.x))
            break;

        if (!(inFile >> point.y))
            break;

        if (!(inFile >> point.z))
            break;

        if (listType == SegmentedPointList) {
            if (!(inFile >> segID))
                break;
        }

        // skip items according to property
        bool error = false;
        for (int i=0; i<skipItems && !error; ++i) {
            error = !(inFile >> dummy);
        }

        if(error)
            break;

        if (listType == SegmentedPointList) {
            // if new segment started, append vector for it
            if (segID > lastSegID) {
                segments.push_back(std::vector<tgt::vec3>());
                lastSegID = segID;
            }

            // append next skelPoint to last segment vector
            segments.back().push_back(point);
        }
        else if (listType == PointList) {
            tgtAssert(pointListGeometry, "no pointListGeometry");
            pointListGeometry->addPoint(point);
        }
    }

    if (listType == SegmentedPointList) {
        size_t numPoints = 0;
        for (size_t i=0; i<segments.size(); ++i)
            numPoints += segments[i].size();
        LINFO("Read " << segments.size() << " segments consisting of " << numPoints << " points.");
        tgtAssert(pointSegmentListGeometry, "no pointSegmentListGeometry");
        pointSegmentListGeometry->setData(segments);
    }
    else if (listType == PointList) {
        tgtAssert(pointListGeometry, "no pointListGeometry");
        LINFO("Read " << pointListGeometry->getNumPoints() << " points.");
    }

    if (pointListGeometry)
        return pointListGeometry;
    else if (pointSegmentListGeometry)
        return pointSegmentListGeometry;
    else {
        tgtAssert(false, "no geometry created");
        return 0;
    }
}

void GeometrySource::clearGeometry() {
    outport_.setData(0);
    geometryFile_.set("");
    updatePropertyVisibility();
}

void GeometrySource::forceReload() {
    forceReload_ = true;
    invalidate();
}

void GeometrySource::updatePropertyVisibility() {
    loadGeometry_.setReadOnlyFlag(geometryFile_.get() == "");
    clearGeometry_.setReadOnlyFlag(!outport_.getData());
    skipItemCount_.setVisibleFlag(!geometryType_.isSelected("geometry"));
}

} // namespace voreen
