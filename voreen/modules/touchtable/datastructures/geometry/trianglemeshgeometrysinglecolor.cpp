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

#include "trianglemeshgeometrysinglecolor.h"

namespace voreen{

TriangleMeshGeometrySingleColor::TriangleMeshGeometrySingleColor()
    : TriangleMeshGeometry<VertexNormal>(){

}

void TriangleMeshGeometrySingleColor::render() const{
    glColor4f(getColor().r, getColor().g, getColor().b, getColor().a);
    TriangleMeshGeometry<VertexNormal>::render();
}

TriangleMeshGeometryBase::VertexLayout TriangleMeshGeometrySingleColor::getVertexLayout() const{
    return TriangleMeshGeometryBase::NORMAL;
}

void TriangleMeshGeometrySingleColor::setColor(const tgt::vec4& color){
    color_ = color;
}

tgt::vec4 TriangleMeshGeometrySingleColor::getColor()const{
    return color_;
}

TriangleMeshGeometryBodyParts3D::TriangleMeshGeometryBodyParts3D()
    : TriangleMeshGeometrySingleColor()
    , highlight_(false){

}

void TriangleMeshGeometryBodyParts3D::setHighlight(bool b)const{
    highlight_ = b;
}

bool TriangleMeshGeometryBodyParts3D::getHighlight() const{
    return highlight_;
}

void TriangleMeshGeometryBodyParts3D::setTag(const std::string& str)const{
    tag_ = str;
}

std::string TriangleMeshGeometryBodyParts3D::getTag() const{
    return tag_;
}

void TriangleMeshGeometryBodyParts3D::renderUncolored() const{
    TriangleMeshGeometry<VertexNormal>::render();
}

tgt::vec4 TriangleMeshGeometryBodyParts3D::getColor() const{
    if(highlight_)
        return 1.2f*color_;
    else
        return color_;
}

//-------------------------------------------------------------------------

TriangleMeshGeometryCollectionBodyParts3D::TriangleMeshGeometryCollectionBodyParts3D()
    : TriangleMeshGeometryBase()
    , boundingBoxValid_(false)
    , numTriangles_(0)
    , numTrianglesValid_(true){

}

TriangleMeshGeometryCollectionBodyParts3D::~TriangleMeshGeometryCollectionBodyParts3D(){
    clear();
}

TriangleMeshGeometryBase::VertexLayout TriangleMeshGeometryCollectionBodyParts3D::getVertexLayout() const{
    return TriangleMeshGeometryBase::NORMAL;
}

size_t TriangleMeshGeometryCollectionBodyParts3D::getNumMeshes() const{
    return meshes_.size();
}

size_t TriangleMeshGeometryCollectionBodyParts3D::getNumTriangles() const{
    if(numTrianglesValid_)
        return numTriangles_;

    numTriangles_ = 0;
    for(size_t i = 0; i < meshes_.size() ; ++i){
        numTriangles_ += meshes_[i]->getNumTriangles();
    }
    numTrianglesValid_ = true;
    return numTriangles_;
}

bool TriangleMeshGeometryCollectionBodyParts3D::isEmpty() const{
    return (meshes_.size() == 0);
}

void TriangleMeshGeometryCollectionBodyParts3D::clear(){
    for(size_t i = 0; i < meshes_.size() ; ++i){
        if(meshes_[i])
            delete meshes_[i];
    }
    meshes_.clear();

    invalidate();
}

void TriangleMeshGeometryCollectionBodyParts3D::addMesh(TriangleMeshGeometryBodyParts3D* m){
    meshes_.push_back(m);
    invalidate();
}

TriangleMeshGeometryBodyParts3D* TriangleMeshGeometryCollectionBodyParts3D::removeMesh(size_t i){
    tgtAssert( i < meshes_.size(), "Invalid Mesh Index");
    TriangleMeshGeometryBodyParts3D* m = meshes_[i];
    meshes_.erase(meshes_.begin() + i);
    return m;
}

void TriangleMeshGeometryCollectionBodyParts3D::deleteMesh(size_t i){
    delete removeMesh(i);
}

const TriangleMeshGeometryBodyParts3D* TriangleMeshGeometryCollectionBodyParts3D::getMesh(size_t i) const{
    tgtAssert( i < meshes_.size(), "Invalid Mesh Index");
    return meshes_[i];
}

const TriangleMeshGeometryBodyParts3D* TriangleMeshGeometryCollectionBodyParts3D::operator[](size_t i) const{
    tgtAssert( i < meshes_.size(), "Invalid Mesh Index");
    return meshes_[i];
}

tgt::Bounds TriangleMeshGeometryCollectionBodyParts3D::getBoundingBox(bool transformed) const{
    if(!boundingBoxValid_)
        updateBoundingBox();

    if(transformed)
        return boundingBox_.transform(getTransformationMatrix());
    else
        return boundingBox_;
}

void TriangleMeshGeometryCollectionBodyParts3D::invalidate(){
    boundingBoxValid_ = false;
    numTrianglesValid_ = false;
}

void TriangleMeshGeometryCollectionBodyParts3D::render()const{
    for(size_t i = 0 ; i < meshes_.size() ; ++i){
        meshes_[i]->render();
    }
}

void TriangleMeshGeometryCollectionBodyParts3D::updateBoundingBox() const{
    boundingBox_ = tgt::Bounds();
    for(size_t i = 0; i < meshes_.size() ; ++i){
        boundingBox_.addVolume(meshes_[i]->getBoundingBox());
    }
    boundingBoxValid_ = true;
}

}
