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

#include "voreen/core/datastructures/geometry/glmeshgeometry.h"

namespace voreen {

/////////////// BASE CLASS IMPLEMENTATION /////////////////////////////////
GlMeshGeometryBase::GlMeshGeometryBase()
    : Geometry()
    , hasTransparentComponents_(false)
    , textureData_(0)
    , textureOwner_(false)
    , primitiveType_(GL_TRIANGLES)
    , usesIndexedDrawing_(false)
{}

GlMeshGeometryBase::~GlMeshGeometryBase() {
    clearTextureData();
}

bool GlMeshGeometryBase::supportsNormals() const {
    return supportsNormals(getVertexLayout());
}

bool GlMeshGeometryBase::supportsColors() const {
    return supportsColors(getVertexLayout());
}

bool GlMeshGeometryBase::supportsTextureData() const {
    return supportsTextureData(getVertexLayout());
}

std::string GlMeshGeometryBase::getShaderDefines(bool useTransparency) const {
    return getShaderDefines(getVertexLayout(), useTransparency);
}

bool GlMeshGeometryBase::supportsNormals(VertexBase::VertexLayout layout) {
    return layout & VertexBase::NORMAL;
}

bool GlMeshGeometryBase::supportsColors(VertexBase::VertexLayout layout) {
    return layout & VertexBase::COLOR;
}

bool GlMeshGeometryBase::supportsTextureData(VertexBase::VertexLayout layout) {
    return layout & VertexBase::TEXCOORD;
}

std::string GlMeshGeometryBase::getShaderDefines(VertexBase::VertexLayout layout, bool useTransparency) {
    std::string header = "";
    if(supportsColors(layout))
        header += "#define USE_COLOR\n";
    if(supportsNormals(layout))
        header += "#define USE_NORMAL\n";
    if(supportsTextureData(layout))
        header += "#define USE_TEXCOORD\n";
    if(useTransparency)
        header += "#define USE_TRANSPARENCY\n";
    return header;
}

GLenum GlMeshGeometryBase::getPrimitiveType() const {
    return primitiveType_;
}

void GlMeshGeometryBase::setPrimitiveType(GLenum primitiveType) {

    // use assertions to check type
    tgtAssert(  (primitiveType == GL_POINTS) || (primitiveType == GL_LINE_STRIP) || (primitiveType == GL_LINE_LOOP)
             || (primitiveType == GL_LINES) || (primitiveType == GL_LINE_STRIP_ADJACENCY) || (primitiveType == GL_LINES_ADJACENCY)
             || (primitiveType == GL_TRIANGLE_STRIP) || (primitiveType == GL_TRIANGLE_FAN) || (primitiveType == GL_TRIANGLES)
             || (primitiveType == GL_TRIANGLE_STRIP_ADJACENCY) || (primitiveType == GL_TRIANGLES_ADJACENCY ), "invalid primitive type" );

    primitiveType_ = primitiveType;
}

bool GlMeshGeometryBase::usesIndexedDrawing() const {
    return usesIndexedDrawing_;
}

void GlMeshGeometryBase::enableIndexedDrawing() {
    if (getNumIndices() == 0) {
        LERROR("Cannot enable indexed drawing: no indices present in mesh.");
        usesIndexedDrawing_ = false;
    }
    else
        usesIndexedDrawing_ = true;
}

void GlMeshGeometryBase::disableIndexedDrawing() {
    usesIndexedDrawing_ = false;
}

void GlMeshGeometryBase::serialize(Serializer& s) const {
    Geometry::serialize(s);
    s.serialize("textureFilenames", textureFilenames_, "filename");
    s.serialize("usesIndexedDrawing", usesIndexedDrawing_);

    // TODO: serialize this as string?
    s.serialize("primitiveType", static_cast<int>(primitiveType_));
}

void GlMeshGeometryBase::deserialize(Deserializer& d) {
    Geometry::deserialize(d);
    textureFilenames_.clear();
    try{
        std::vector<std::string> filenames;
        d.deserialize("textureFilenames", filenames, "filename");
        textureFilenames_ = filenames;
    } catch(SerializationException&) {
        d.removeLastError();
    }
    d.deserialize("usesIndexedDrawing", usesIndexedDrawing_);

    // TODO: use strings for ptimitive types?
    int primType;
    d.deserialize("primitiveType", primType);
    primitiveType_ = static_cast<GLenum>(primType);
}

//----------------------------------- texture stuff

void GlMeshGeometryBase::loadTextureData() {
    if(textureOwner_)
        delete textureData_;
    textureData_ = 0;
    if(textureFilenames_.empty())
        return;
    textureData_ = GLSL::createArrayTexture(textureFilenames_);
    if(!textureData_) {
        LERRORC("voreen.GlMeshGeometryBase", "Failed to load textures! Texturing will be disabled for this mesh.");
        return;
    }
}

void GlMeshGeometryBase::clearTextureData() {
    if(textureOwner_)
        delete textureData_;
    textureFilenames_.clear();
    textureData_ = 0;
}

void GlMeshGeometryBase::setTextureData(tgt::Texture* tex, bool owner) {
    clearTextureData();
    textureData_ = tex;
    textureOwner_ = owner;
}

tgt::Texture* GlMeshGeometryBase::getTextureData() const {
    if(!textureFilenames_.empty() && !textureData_) {
        LERRORC("voreen.GlMeshGeometryBase", "Texture filenames found but no texture(s). Please call loadTextureData() before using mesh texture(s), returning 0 at the moment.");
        return 0;
    }
    return textureData_;
}

bool GlMeshGeometryBase::hasTextureData() const {
    return textureData_ != 0;
}

void GlMeshGeometryBase::loadTextureData(const std::vector<std::string>& filenames) {
    clearTextureData();
    if(filenames.empty())
        return;
    textureData_ = GLSL::createArrayTexture(filenames);
    if(!textureData_) {
        LERRORC("voreen.GlMeshGeometryBase", "Failed to load textures! Texturing will be disabled for this mesh.");
        return;
    }
    textureFilenames_ = filenames;
    textureOwner_ = true;
}

void GlMeshGeometryBase::loadTextureData(const std::string& filename) {
    std::vector<std::string> tmp;
    tmp.push_back(filename);
    loadTextureData(tmp);
}

void GlMeshGeometryBase::setTextureDataAndFilenames(tgt::Texture* tex, const std::vector<std::string>& filenames, bool owner) {
    setTextureData(tex, owner);
    textureFilenames_ = filenames;
}

//----------------------------------- transparency

bool GlMeshGeometryBase::isTransparent() const {
    return hasTransparentComponents_;
}
void GlMeshGeometryBase::setHasTransparentComponents(bool hasTransparentComponents) {
    hasTransparentComponents_ = hasTransparentComponents;
}


} // namespace
