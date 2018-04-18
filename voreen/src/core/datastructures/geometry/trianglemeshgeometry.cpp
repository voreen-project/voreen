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

#include "tgt/texturemanager.h"
#include "voreen/core/datastructures/geometry/trianglemeshgeometry.h"
#include "voreen/core/utils/glsl.h"

namespace voreen {

using tgt::vec3;
using tgt::vec4;

TriangleMeshGeometryBase::TriangleMeshGeometryBase()
    : Geometry(), textureData_(0), textureOwner_(false), hasTransparentComponents_(false)
{}

TriangleMeshGeometryBase::~TriangleMeshGeometryBase() {
    clearTextureData();
}

bool TriangleMeshGeometryBase::supportsNormals() const {
    return getVertexLayout() & NORMAL;
}

bool TriangleMeshGeometryBase::supportsColors() const {
    return getVertexLayout() & COLOR;
}

bool TriangleMeshGeometryBase::supportsTextureData() const {
    return getVertexLayout() & TEXCOORD;
}

void TriangleMeshGeometryBase::loadTextureData() {
    if(textureOwner_)
        delete textureData_;
    textureData_ = 0;
    if(textureFilenames_.empty())
        return;
    textureData_ = GLSL::createArrayTexture(textureFilenames_);
    if(!textureData_) {
        LERRORC("voreen.TriangleMeshGeometryBase", "Failed to load textures! Texturing will be disabled for this mesh.");
        return;
    }
}

void TriangleMeshGeometryBase::clearTextureData() {
    if(textureOwner_)
        delete textureData_;
    textureFilenames_.clear();
    textureData_ = 0;
}

void TriangleMeshGeometryBase::setTextureData(tgt::Texture* tex, bool owner) {
    clearTextureData();
    textureData_ = tex;
    textureOwner_ = owner;
}

tgt::Texture* TriangleMeshGeometryBase::getTextureData() const {
    if(!textureFilenames_.empty() && !textureData_) {
        LERRORC("voreen.TriangleMeshGeometryBase", "Texture filenames found but no texture(s). Please call loadTextureData() before using mesh texture(s), returning 0 at the moment.");
        return 0;
    }
    return textureData_;
}

bool TriangleMeshGeometryBase::hasTextureData() const {
    return textureData_ != 0;
}

void TriangleMeshGeometryBase::loadTextureData(const std::vector<std::string>& filenames) {
    clearTextureData();
    if(filenames.empty())
        return;
    textureData_ = GLSL::createArrayTexture(filenames);
    if(!textureData_) {
        LERRORC("voreen.TriangleMeshGeometryBase", "Failed to load textures! Texturing will be disabled for this mesh.");
        return;
    }
    textureFilenames_ = filenames;
    textureOwner_ = true;
}

void TriangleMeshGeometryBase::loadTextureData(const std::string& filename) {
    std::vector<std::string> tmp;
    tmp.push_back(filename);
    loadTextureData(tmp);
}

void TriangleMeshGeometryBase::setTextureDataAndFilenames(tgt::Texture* tex, const std::vector<std::string>& filenames, bool owner) {
    setTextureData(tex, owner);
    textureFilenames_ = filenames;
}

void TriangleMeshGeometryBase::serialize(Serializer& s) const {
    Geometry::serialize(s);
    s.serialize("textureFilenames", textureFilenames_, "filename");
}

void TriangleMeshGeometryBase::deserialize(Deserializer& s) {
    Geometry::deserialize(s);
    textureFilenames_.clear();
    try{
        std::vector<std::string> filenames;
        s.deserialize("textureFilenames", filenames, "filename");
        textureFilenames_ = filenames;
    } catch(SerializationException&) {
        s.removeLastError();
    }
}

std::string TriangleMeshGeometryBase::getShaderDefines() const {
    std::string header = "";
    if(getVertexLayout() & COLOR)
        header += "#define USE_COLOR\n";
    if(getVertexLayout() & NORMAL)
        header += "#define USE_NORMAL\n";
    if(getVertexLayout() & TEXCOORD)
        header += "#define USE_TEXCOORD\n";
    return header;
}
bool TriangleMeshGeometryBase::isTransparent() const {
    return hasTransparentComponents_;
}
void TriangleMeshGeometryBase::setHasTransparentComponents(bool hasTransparentComponents) {
    hasTransparentComponents_ = hasTransparentComponents;
}

// --------------------------------------------------------------------

TriangleMeshGeometryNormal* TriangleMeshGeometryNormal::createCube(VertexType llfVertex, VertexType urbVertex) {
    TriangleMeshGeometryNormal* ret = new TriangleMeshGeometryNormal();
    ret->addCube(llfVertex, urbVertex);
    return ret;
}

void TriangleMeshGeometryNormal::addCube(VertexType llfVertex, VertexType urbVertex) {
    // expecting llfVertex.pos_ < urbVertex.pos_
    if (llfVertex.pos_.x > urbVertex.pos_.x) {
        std::swap(llfVertex.pos_.x, urbVertex.pos_.x);
        std::swap(llfVertex.normal_.x, urbVertex.normal_.x);
    }
    if (llfVertex.pos_.y > urbVertex.pos_.y) {
        std::swap(llfVertex.pos_.y, urbVertex.pos_.y);
        std::swap(llfVertex.normal_.y, urbVertex.normal_.y);
    }
    if (llfVertex.pos_.z > urbVertex.pos_.z) {
        std::swap(llfVertex.pos_.z, urbVertex.pos_.z);
        std::swap(llfVertex.normal_.z, urbVertex.normal_.z);
    }

    VertexNormal llf(tgt::vec3(llfVertex.pos_.x, llfVertex.pos_.y, llfVertex.pos_.z), tgt::vec3(llfVertex.normal_.x, llfVertex.normal_.y, llfVertex.normal_.z));
    VertexNormal lrf(tgt::vec3(urbVertex.pos_.x, llfVertex.pos_.y, llfVertex.pos_.z), tgt::vec3(urbVertex.normal_.x, llfVertex.normal_.y, llfVertex.normal_.z));
    VertexNormal lrb(tgt::vec3(urbVertex.pos_.x, llfVertex.pos_.y, urbVertex.pos_.z), tgt::vec3(urbVertex.normal_.x, llfVertex.normal_.y, urbVertex.normal_.z));
    VertexNormal llb(tgt::vec3(llfVertex.pos_.x, llfVertex.pos_.y, urbVertex.pos_.z), tgt::vec3(llfVertex.normal_.x, llfVertex.normal_.y, urbVertex.normal_.z));

    VertexNormal ulb(tgt::vec3(llfVertex.pos_.x, urbVertex.pos_.y, urbVertex.pos_.z), tgt::vec3(llfVertex.normal_.x, urbVertex.normal_.y, urbVertex.normal_.z));
    VertexNormal ulf(tgt::vec3(llfVertex.pos_.x, urbVertex.pos_.y, llfVertex.pos_.z), tgt::vec3(llfVertex.normal_.x, urbVertex.normal_.y, llfVertex.normal_.z));
    VertexNormal urf(tgt::vec3(urbVertex.pos_.x, urbVertex.pos_.y, llfVertex.pos_.z), tgt::vec3(urbVertex.normal_.x, urbVertex.normal_.y, llfVertex.normal_.z));
    VertexNormal urb(tgt::vec3(urbVertex.pos_.x, urbVertex.pos_.y, urbVertex.pos_.z), tgt::vec3(urbVertex.normal_.x, urbVertex.normal_.y, urbVertex.normal_.z));

    addTriangle(TriangleType(urb, urf, ulf));
    addTriangle(TriangleType(urb, ulf, ulb));

    addTriangle(TriangleType(llf, ulf, urf));
    addTriangle(TriangleType(llf, urf, lrf));

    addTriangle(TriangleType(llf, llb, ulb));
    addTriangle(TriangleType(llf, ulb, ulf));

    addTriangle(TriangleType(urb, ulb, llb));
    addTriangle(TriangleType(urb, llb, lrb));

    addTriangle(TriangleType(urb, lrb, lrf));
    addTriangle(TriangleType(urb, lrf, urf));

    addTriangle(TriangleType(llf, lrf, lrb));
    addTriangle(TriangleType(llf, lrb, llb));

    invalidate();
}

//-------------------------------------------------------------------------------------------------

TriangleMeshGeometryColorNormal* TriangleMeshGeometryColorNormal::createCube(tgt::vec3 coordLlf, tgt::vec3 coordUrb, tgt::vec3 colorLlf, tgt::vec3 colorUrb, float alpha, tgt::vec3 texLlf, tgt::vec3 texUrb) {
    // expecting coordLlf < coordUrb
    if (coordLlf.x > coordUrb.x) {
        std::swap(coordLlf.x, coordUrb.x);
        std::swap(texLlf.x, texUrb.x);
        std::swap(colorLlf.x, colorUrb.x);
    }
    if (coordLlf.y > coordUrb.y) {
        std::swap(coordLlf.y, coordUrb.y);
        std::swap(texLlf.y, texUrb.y);
        std::swap(colorLlf.y, colorUrb.y);
    }
    if (coordLlf.z > coordUrb.z) {
        std::swap(coordLlf.z, coordUrb.z);
        std::swap(texLlf.z, texUrb.z);
        std::swap(colorLlf.z, colorUrb.z);
    }

    TriangleMeshGeometryColorNormal* ret = new TriangleMeshGeometryColorNormal();

    VertexColorNormal llf(tgt::vec3(coordLlf.x, coordLlf.y, coordLlf.z), tgt::vec4(colorLlf.x, colorLlf.y, colorLlf.z, alpha), tgt::vec3(texLlf.x, texLlf.y, texLlf.z));
    VertexColorNormal lrf(tgt::vec3(coordUrb.x, coordLlf.y, coordLlf.z), tgt::vec4(colorUrb.x, colorLlf.y, colorLlf.z, alpha), tgt::vec3(texUrb.x, texLlf.y, texLlf.z));
    VertexColorNormal lrb(tgt::vec3(coordUrb.x, coordLlf.y, coordUrb.z), tgt::vec4(colorUrb.x, colorLlf.y, colorUrb.z, alpha), tgt::vec3(texUrb.x, texLlf.y, texUrb.z));
    VertexColorNormal llb(tgt::vec3(coordLlf.x, coordLlf.y, coordUrb.z), tgt::vec4(colorLlf.x, colorLlf.y, colorUrb.z, alpha), tgt::vec3(texLlf.x, texLlf.y, texUrb.z));

    VertexColorNormal ulb(tgt::vec3(coordLlf.x, coordUrb.y, coordUrb.z), tgt::vec4(colorLlf.x, colorUrb.y, colorUrb.z, alpha), tgt::vec3(texLlf.x, texUrb.y, texUrb.z));
    VertexColorNormal ulf(tgt::vec3(coordLlf.x, coordUrb.y, coordLlf.z), tgt::vec4(colorLlf.x, colorUrb.y, colorLlf.z, alpha), tgt::vec3(texLlf.x, texUrb.y, texLlf.z));
    VertexColorNormal urf(tgt::vec3(coordUrb.x, coordUrb.y, coordLlf.z), tgt::vec4(colorUrb.x, colorUrb.y, colorLlf.z, alpha), tgt::vec3(texUrb.x, texUrb.y, texLlf.z));
    VertexColorNormal urb(tgt::vec3(coordUrb.x, coordUrb.y, coordUrb.z), tgt::vec4(colorUrb.x, colorUrb.y, colorUrb.z, alpha), tgt::vec3(texUrb.x, texUrb.y, texUrb.z));

    ret->addTriangle(TriangleType(urb, urf, ulf));
    ret->addTriangle(TriangleType(urb, ulf, ulb));

    ret->addTriangle(TriangleType(llf, ulf, urf));
    ret->addTriangle(TriangleType(llf, urf, lrf));

    ret->addTriangle(TriangleType(llf, llb, ulb));
    ret->addTriangle(TriangleType(llf, ulb, ulf));

    ret->addTriangle(TriangleType(urb, ulb, llb));
    ret->addTriangle(TriangleType(urb, llb, lrb));

    ret->addTriangle(TriangleType(urb, lrb, lrf));
    ret->addTriangle(TriangleType(urb, lrf, urf));

    ret->addTriangle(TriangleType(llf, lrf, lrb));
    ret->addTriangle(TriangleType(llf, lrb, llb));

    return ret;
}

TriangleMeshGeometryColorNormal* TriangleMeshGeometryColorNormal::createCube(tgt::vec3 coordLlf, tgt::vec3 coordUrb, tgt::vec3 colorLlf, tgt::vec3 colorUrb, float alpha) {
    // expecting coordLlf < coordUrb
    if (coordLlf.x > coordUrb.x) {
        std::swap(coordLlf.x, coordUrb.x);
        std::swap(colorLlf.x, colorUrb.x);
    }
    if (coordLlf.y > coordUrb.y) {
        std::swap(coordLlf.y, coordUrb.y);
        std::swap(colorLlf.y, colorUrb.y);
    }
    if (coordLlf.z > coordUrb.z) {
        std::swap(coordLlf.z, coordUrb.z);
        std::swap(colorLlf.z, colorUrb.z);
    }

    TriangleMeshGeometryColorNormal* ret = new TriangleMeshGeometryColorNormal();

    VertexColorNormal llf(tgt::vec3(coordLlf.x, coordLlf.y, coordLlf.z), tgt::vec4(colorLlf.x, colorLlf.y, colorLlf.z, alpha), tgt::vec3(0.0f));
    VertexColorNormal lrf(tgt::vec3(coordUrb.x, coordLlf.y, coordLlf.z), tgt::vec4(colorUrb.x, colorLlf.y, colorLlf.z, alpha), tgt::vec3(0.0f));
    VertexColorNormal lrb(tgt::vec3(coordUrb.x, coordLlf.y, coordUrb.z), tgt::vec4(colorUrb.x, colorLlf.y, colorUrb.z, alpha), tgt::vec3(0.0f));
    VertexColorNormal llb(tgt::vec3(coordLlf.x, coordLlf.y, coordUrb.z), tgt::vec4(colorLlf.x, colorLlf.y, colorUrb.z, alpha), tgt::vec3(0.0f));

    VertexColorNormal ulb(tgt::vec3(coordLlf.x, coordUrb.y, coordUrb.z), tgt::vec4(colorLlf.x, colorUrb.y, colorUrb.z, alpha), tgt::vec3(0.0f));
    VertexColorNormal ulf(tgt::vec3(coordLlf.x, coordUrb.y, coordLlf.z), tgt::vec4(colorLlf.x, colorUrb.y, colorLlf.z, alpha), tgt::vec3(0.0f));
    VertexColorNormal urf(tgt::vec3(coordUrb.x, coordUrb.y, coordLlf.z), tgt::vec4(colorUrb.x, colorUrb.y, colorLlf.z, alpha), tgt::vec3(0.0f));
    VertexColorNormal urb(tgt::vec3(coordUrb.x, coordUrb.y, coordUrb.z), tgt::vec4(colorUrb.x, colorUrb.y, colorUrb.z, alpha), tgt::vec3(0.0f));

    llf.normal_ = lrf.normal_ = lrb.normal_ = llb.normal_ = ulb.normal_ = ulf.normal_ = urf.normal_ = urb.normal_ = tgt::vec3(0.0f, 1.0f, 0.0f);
    ret->addTriangle(TriangleType(urb, urf, ulf));
    ret->addTriangle(TriangleType(urb, ulf, ulb));

    llf.normal_ = lrf.normal_ = lrb.normal_ = llb.normal_ = ulb.normal_ = ulf.normal_ = urf.normal_ = urb.normal_ = tgt::vec3(0.0f, 0.0f, -1.0f);
    ret->addTriangle(TriangleType(llf, ulf, urf));
    ret->addTriangle(TriangleType(llf, urf, lrf));

    llf.normal_ = lrf.normal_ = lrb.normal_ = llb.normal_ = ulb.normal_ = ulf.normal_ = urf.normal_ = urb.normal_ = tgt::vec3(-1.0f, 0.0f, 0.0f);
    ret->addTriangle(TriangleType(llf, llb, ulb));
    ret->addTriangle(TriangleType(llf, ulb, ulf));

    llf.normal_ = lrf.normal_ = lrb.normal_ = llb.normal_ = ulb.normal_ = ulf.normal_ = urf.normal_ = urb.normal_ = tgt::vec3(0.0f, 0.0f, 1.0f);
    ret->addTriangle(TriangleType(urb, ulb, llb));
    ret->addTriangle(TriangleType(urb, llb, lrb));

    llf.normal_ = lrf.normal_ = lrb.normal_ = llb.normal_ = ulb.normal_ = ulf.normal_ = urf.normal_ = urb.normal_ = tgt::vec3(1.0f, 0.0f, 0.0f);
    ret->addTriangle(TriangleType(urb, lrb, lrf));
    ret->addTriangle(TriangleType(urb, lrf, urf));

    llf.normal_ = lrf.normal_ = lrb.normal_ = llb.normal_ = ulb.normal_ = ulf.normal_ = urf.normal_ = urb.normal_ = tgt::vec3(0.0f, -1.0f, 0.0f);
    ret->addTriangle(TriangleType(llf, lrf, lrb));
    ret->addTriangle(TriangleType(llf, lrb, llb));

    return ret;
}

void TriangleMeshGeometryColorNormal::addCube(VertexNormal llfVertex, VertexNormal urbVertex) {
    // expecting llfVertex.pos_ < urbVertex.pos_
    if (llfVertex.pos_.x > urbVertex.pos_.x) {
        std::swap(llfVertex.pos_.x, urbVertex.pos_.x);
        std::swap(llfVertex.normal_.x, urbVertex.normal_.x);
    }
    if (llfVertex.pos_.y > urbVertex.pos_.y) {
        std::swap(llfVertex.pos_.y, urbVertex.pos_.y);
        std::swap(llfVertex.normal_.y, urbVertex.normal_.y);
    }
    if (llfVertex.pos_.z > urbVertex.pos_.z) {
        std::swap(llfVertex.pos_.z, urbVertex.pos_.z);
        std::swap(llfVertex.normal_.z, urbVertex.normal_.z);
    }

    float alpha = 1.0f;

    VertexColorNormal llf(tgt::vec3(llfVertex.pos_.x, llfVertex.pos_.y, llfVertex.pos_.z), tgt::vec4(llfVertex.normal_.x, llfVertex.normal_.y, llfVertex.normal_.z, alpha), tgt::vec3(0.0f));
    VertexColorNormal lrf(tgt::vec3(urbVertex.pos_.x, llfVertex.pos_.y, llfVertex.pos_.z), tgt::vec4(urbVertex.normal_.x, llfVertex.normal_.y, llfVertex.normal_.z, alpha), tgt::vec3(0.0f));
    VertexColorNormal lrb(tgt::vec3(urbVertex.pos_.x, llfVertex.pos_.y, urbVertex.pos_.z), tgt::vec4(urbVertex.normal_.x, llfVertex.normal_.y, urbVertex.normal_.z, alpha), tgt::vec3(0.0f));
    VertexColorNormal llb(tgt::vec3(llfVertex.pos_.x, llfVertex.pos_.y, urbVertex.pos_.z), tgt::vec4(llfVertex.normal_.x, llfVertex.normal_.y, urbVertex.normal_.z, alpha), tgt::vec3(0.0f));

    VertexColorNormal ulb(tgt::vec3(llfVertex.pos_.x, urbVertex.pos_.y, urbVertex.pos_.z), tgt::vec4(llfVertex.normal_.x, urbVertex.normal_.y, urbVertex.normal_.z, alpha), tgt::vec3(0.0f));
    VertexColorNormal ulf(tgt::vec3(llfVertex.pos_.x, urbVertex.pos_.y, llfVertex.pos_.z), tgt::vec4(llfVertex.normal_.x, urbVertex.normal_.y, llfVertex.normal_.z, alpha), tgt::vec3(0.0f));
    VertexColorNormal urf(tgt::vec3(urbVertex.pos_.x, urbVertex.pos_.y, llfVertex.pos_.z), tgt::vec4(urbVertex.normal_.x, urbVertex.normal_.y, llfVertex.normal_.z, alpha), tgt::vec3(0.0f));
    VertexColorNormal urb(tgt::vec3(urbVertex.pos_.x, urbVertex.pos_.y, urbVertex.pos_.z), tgt::vec4(urbVertex.normal_.x, urbVertex.normal_.y, urbVertex.normal_.z, alpha), tgt::vec3(0.0f));

    llf.normal_ = lrf.normal_ = lrb.normal_ = llb.normal_ = ulb.normal_ = ulf.normal_ = urf.normal_ = urb.normal_ = tgt::vec3(0.0f, 1.0f, 0.0f);
    addTriangle(TriangleType(urb, urf, ulf));
    addTriangle(TriangleType(urb, ulf, ulb));

    llf.normal_ = lrf.normal_ = lrb.normal_ = llb.normal_ = ulb.normal_ = ulf.normal_ = urf.normal_ = urb.normal_ = tgt::vec3(0.0f, 0.0f, -1.0f);
    addTriangle(TriangleType(llf, ulf, urf));
    addTriangle(TriangleType(llf, urf, lrf));

    llf.normal_ = lrf.normal_ = lrb.normal_ = llb.normal_ = ulb.normal_ = ulf.normal_ = urf.normal_ = urb.normal_ = tgt::vec3(-1.0f, 0.0f, 0.0f);
    addTriangle(TriangleType(llf, llb, ulb));
    addTriangle(TriangleType(llf, ulb, ulf));

    llf.normal_ = lrf.normal_ = lrb.normal_ = llb.normal_ = ulb.normal_ = ulf.normal_ = urf.normal_ = urb.normal_ = tgt::vec3(0.0f, 0.0f, 1.0f);
    addTriangle(TriangleType(urb, ulb, llb));
    addTriangle(TriangleType(urb, llb, lrb));

    llf.normal_ = lrf.normal_ = lrb.normal_ = llb.normal_ = ulb.normal_ = ulf.normal_ = urf.normal_ = urb.normal_ = tgt::vec3(1.0f, 0.0f, 0.0f);
    addTriangle(TriangleType(urb, lrb, lrf));
    addTriangle(TriangleType(urb, lrf, urf));

    llf.normal_ = lrf.normal_ = lrb.normal_ = llb.normal_ = ulb.normal_ = ulf.normal_ = urf.normal_ = urb.normal_ = tgt::vec3(0.0f, -1.0f, 0.0f);
    addTriangle(TriangleType(llf, lrf, lrb));
    addTriangle(TriangleType(llf, lrb, llb));

    invalidate();
}

void TriangleMeshGeometryColorNormal::addQuad(const VertexNormal& v1, const VertexNormal& v2, const VertexNormal& v3, const VertexNormal& v4) {
    vec3 a = v2.pos_ - v1.pos_;
    vec3 b = v3.pos_ - v1.pos_;
    vec3 n = normalize(cross(a, b));
    float alpha = 1.0f;

    addTriangle(TriangleType(VertexColorNormal(v1.pos_, vec4(v1.normal_, alpha), n), VertexColorNormal(v2.pos_, vec4(v2.normal_, alpha), n), VertexColorNormal(v3.pos_, vec4(v3.normal_, alpha), n)));
    addTriangle(TriangleType(VertexColorNormal(v1.pos_, vec4(v1.normal_, alpha), n), VertexColorNormal(v3.pos_, vec4(v3.normal_, alpha), n), VertexColorNormal(v4.pos_, vec4(v4.normal_, alpha), n)));

    invalidate();
}

//-------------------------------------------------------------------------------------------------

TriangleMeshGeometrySimple* TriangleMeshGeometrySimple::createCube(tgt::vec3 coordLlf, tgt::vec3 coordUrb) {
    // expecting coordLlf < coordUrb
    if (coordLlf.x > coordUrb.x) {
        std::swap(coordLlf.x, coordUrb.x);
    }
    if (coordLlf.y > coordUrb.y) {
        std::swap(coordLlf.y, coordUrb.y);
    }
    if (coordLlf.z > coordUrb.z) {
        std::swap(coordLlf.z, coordUrb.z);
    }

    TriangleMeshGeometrySimple* ret = new TriangleMeshGeometrySimple();

    VertexBase llf(tgt::vec3(coordLlf.x, coordLlf.y, coordLlf.z));
    VertexBase lrf(tgt::vec3(coordUrb.x, coordLlf.y, coordLlf.z));
    VertexBase lrb(tgt::vec3(coordUrb.x, coordLlf.y, coordUrb.z));
    VertexBase llb(tgt::vec3(coordLlf.x, coordLlf.y, coordUrb.z));

    VertexBase ulb(tgt::vec3(coordLlf.x, coordUrb.y, coordUrb.z));
    VertexBase ulf(tgt::vec3(coordLlf.x, coordUrb.y, coordLlf.z));
    VertexBase urf(tgt::vec3(coordUrb.x, coordUrb.y, coordLlf.z));
    VertexBase urb(tgt::vec3(coordUrb.x, coordUrb.y, coordUrb.z));

    ret->addTriangle(TriangleType(urb, urf, ulf));
    ret->addTriangle(TriangleType(urb, ulf, ulb));

    ret->addTriangle(TriangleType(llf, ulf, urf));
    ret->addTriangle(TriangleType(llf, urf, lrf));

    ret->addTriangle(TriangleType(llf, llb, ulb));
    ret->addTriangle(TriangleType(llf, ulb, ulf));

    ret->addTriangle(TriangleType(urb, ulb, llb));
    ret->addTriangle(TriangleType(urb, llb, lrb));

    ret->addTriangle(TriangleType(urb, lrb, lrf));
    ret->addTriangle(TriangleType(urb, lrf, urf));

    ret->addTriangle(TriangleType(llf, lrf, lrb));
    ret->addTriangle(TriangleType(llf, lrb, llb));

    return ret;
}

//-----------------------------------------------------------------------------

TriangleMeshGeometryColor* TriangleMeshGeometryColor::createCube(tgt::vec3 coordLlf, tgt::vec3 coordUrb, tgt::vec3 color, float alpha) {
    // expecting coordLlf < coordUrb
    if (coordLlf.x > coordUrb.x) {
        std::swap(coordLlf.x, coordUrb.x);
    }
    if (coordLlf.y > coordUrb.y) {
        std::swap(coordLlf.y, coordUrb.y);
    }
    if (coordLlf.z > coordUrb.z) {
        std::swap(coordLlf.z, coordUrb.z);
    }

    TriangleMeshGeometryColor* ret = new TriangleMeshGeometryColor();

    VertexColor llf(tgt::vec3(coordLlf.x, coordLlf.y, coordLlf.z), tgt::vec4(color, alpha));
    VertexColor lrf(tgt::vec3(coordUrb.x, coordLlf.y, coordLlf.z), tgt::vec4(color, alpha));
    VertexColor lrb(tgt::vec3(coordUrb.x, coordLlf.y, coordUrb.z), tgt::vec4(color, alpha));
    VertexColor llb(tgt::vec3(coordLlf.x, coordLlf.y, coordUrb.z), tgt::vec4(color, alpha));

    VertexColor ulb(tgt::vec3(coordLlf.x, coordUrb.y, coordUrb.z), tgt::vec4(color, alpha));
    VertexColor ulf(tgt::vec3(coordLlf.x, coordUrb.y, coordLlf.z), tgt::vec4(color, alpha));
    VertexColor urf(tgt::vec3(coordUrb.x, coordUrb.y, coordLlf.z), tgt::vec4(color, alpha));
    VertexColor urb(tgt::vec3(coordUrb.x, coordUrb.y, coordUrb.z), tgt::vec4(color, alpha));

    ret->addTriangle(TriangleType(urb, urf, ulf));
    ret->addTriangle(TriangleType(urb, ulf, ulb));

    ret->addTriangle(TriangleType(llf, ulf, urf));
    ret->addTriangle(TriangleType(llf, urf, lrf));

    ret->addTriangle(TriangleType(llf, llb, ulb));
    ret->addTriangle(TriangleType(llf, ulb, ulf));

    ret->addTriangle(TriangleType(urb, ulb, llb));
    ret->addTriangle(TriangleType(urb, llb, lrb));

    ret->addTriangle(TriangleType(urb, lrb, lrf));
    ret->addTriangle(TriangleType(urb, lrf, urf));

    ret->addTriangle(TriangleType(llf, lrf, lrb));
    ret->addTriangle(TriangleType(llf, lrb, llb));

    return ret;
}

} // namespace
