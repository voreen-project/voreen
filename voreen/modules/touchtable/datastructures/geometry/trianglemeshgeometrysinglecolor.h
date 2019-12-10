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

#ifndef VRN_TRIANGLECLIPMESH_H
#define VRN_TRIANGLECLIPMESH_H

#include <vector>
#include "tgt/tgt_gl.h"
#include "tgt/glmath.h"
#include "voreen/core/datastructures/geometry/trianglemeshgeometry.h"
#include "voreen/core/datastructures/geometry/geometry.h"
#include "voreen/core/datastructures/geometry/vertex.h"

namespace voreen{

class VRN_CORE_API TriangleMeshGeometrySingleColor : public TriangleMeshGeometry<VertexNormal>{
public:
    TriangleMeshGeometrySingleColor();
    virtual Geometry* create() const { return new TriangleMeshGeometrySingleColor();}
    virtual std::string getClassName() const { return "TriangleMeshGeometrySingleColor";}
    VertexLayout getVertexLayout() const;
    virtual void render() const;

    void setColor(const tgt::vec4& color);
    virtual tgt::vec4 getColor()const;
protected:
    tgt::vec4 color_;
};

class VRN_CORE_API TriangleMeshGeometryBodyParts3D : public TriangleMeshGeometrySingleColor{
public:
    TriangleMeshGeometryBodyParts3D();
    virtual Geometry* create() const {return new TriangleMeshGeometryBodyParts3D();}
    virtual std::string getClassName() const { return "TriangleMeshGeometryBodyParts3D"; }

    void setHighlight( bool b ) const;
    bool getHighlight() const;

    void setTag(const std::string& t) const;
    std::string getTag() const;

    void renderUncolored() const;

    tgt::vec4 getColor() const;
private:
    mutable bool highlight_;
    mutable std::string tag_;
};

class VRN_CORE_API TriangleMeshGeometryCollectionBodyParts3D : public TriangleMeshGeometryBase{
public:
    TriangleMeshGeometryCollectionBodyParts3D();

    ~TriangleMeshGeometryCollectionBodyParts3D();

    virtual Geometry* create() const {return new TriangleMeshGeometryCollectionBodyParts3D();}
    virtual std::string getClassName() const { return "TriangleMeshGeometryCollectionBodyParts3D";}
    VertexLayout getVertexLayout() const;

    size_t getNumMeshes() const;

    size_t getNumTriangles() const;

    virtual bool isEmpty() const;

    virtual void clear();

    virtual void calculateNormals() {}; //TODO

    void addMesh(TriangleMeshGeometryBodyParts3D* m);

    TriangleMeshGeometryBodyParts3D* removeMesh(size_t i);

    void deleteMesh(size_t i);

    const TriangleMeshGeometryBodyParts3D* getMesh(size_t i) const;

    const TriangleMeshGeometryBodyParts3D* operator[](size_t i) const;

    virtual tgt::Bounds getBoundingBox(bool transformed = true) const;

    virtual void invalidate();

    virtual void render() const;

private:
    void updateBoundingBox() const;

    mutable tgt::Bounds boundingBox_;
    mutable bool boundingBoxValid_;

    mutable size_t numTriangles_;
    mutable bool numTrianglesValid_;

    std::vector<TriangleMeshGeometryBodyParts3D*> meshes_;
};

//--------------- More Generic Mesh Collections ----------------------------//

template <class V>
class TriangleMeshGeometryCollection : public TriangleMeshGeometryBase{
public:
    typedef TriangleMeshGeometry<V> MeshType;
    typedef V VertexType;

    TriangleMeshGeometryCollection();

    virtual ~TriangleMeshGeometryCollection();

    virtual size_t getNumMeshes() const;

    virtual size_t getNumTriangles() const;

    virtual bool isEmpty() const;

    virtual void clear();

    virtual void calculateNormals() {} //TODO

    void addMesh(MeshType* m);

    MeshType* removeMesh(size_t i);

    void deleteMesh(size_t i);

    const MeshType* getMesh(size_t i) const;

    const MeshType* operator[](size_t i) const;

    virtual tgt::Bounds getBoundingBox(bool transformed = true) const;

    void invalidate();

    virtual void render() const;

protected:
    virtual void updateBoundingBox() const;

    mutable tgt::Bounds boundingBox_;
    mutable bool boundingBoxValid_;

    mutable size_t numTriangles_;
    mutable bool numTrianglesValid_;

    std::vector<MeshType*> meshes_;
};

template <class M>
TriangleMeshGeometryCollection<M>::TriangleMeshGeometryCollection()
    : boundingBox_()
    , boundingBoxValid_(false)
    , numTriangles_(0)
    , numTrianglesValid_(true){

}

template <class M>
TriangleMeshGeometryCollection<M>::~TriangleMeshGeometryCollection(){
    clear();
}

template <class M>
size_t TriangleMeshGeometryCollection<M>::getNumMeshes() const{
    return meshes_.size();
}

template <class M>
size_t TriangleMeshGeometryCollection<M>::getNumTriangles() const{
    if(numTrianglesValid_)
        return numTriangles_;

    numTriangles_ = 0; //shouldnt overflow (practically)
    for(size_t i = 0; i < meshes_.size() ; ++i){
        numTriangles_ += meshes_[i]->getNumTriangles();
    }
    numTrianglesValid_ = true;
    return numTriangles_;
}

template <class M>
bool TriangleMeshGeometryCollection<M>::isEmpty() const {
    return (meshes_.size() == 0);
}

template <class M>
void TriangleMeshGeometryCollection<M>::clear(){
    for(size_t i = 0; i < meshes_.size() ; ++i){
        delete meshes_[i];
    }
    meshes_.clear();

    invalidate();
}

template <class M>
void TriangleMeshGeometryCollection<M>::addMesh(MeshType* m){
    meshes_.push_back(m);
    invalidate();
}

template <class M>
TriangleMeshGeometry<M>* TriangleMeshGeometryCollection<M>::removeMesh(size_t i){
    tgtAssert( i < meshes_.size(), "Invalid Mesh Index");
    MeshType* m = meshes_[i];
    meshes_.erase(meshes_.begin()+i);
    return m;
}

template <class M>
void TriangleMeshGeometryCollection<M>::deleteMesh(size_t i){
    delete removeMesh(i);
}

template <class M>
const typename TriangleMeshGeometryCollection<M>::MeshType* TriangleMeshGeometryCollection<M>::getMesh(size_t i) const{
    tgtAssert( i < meshes_.size(), "Invalid Mesh Index");
    return meshes_[i];
}

template <class M>
const typename TriangleMeshGeometryCollection<M>::MeshType* TriangleMeshGeometryCollection<M>::operator[](size_t i) const{
    tgtAssert( i < meshes_.size(), "Invalid Mesh Index");
    return meshes_[i];
}

template <class M>
tgt::Bounds TriangleMeshGeometryCollection<M>::getBoundingBox(bool transformed) const{
    if(!boundingBoxValid_)
        updateBoundingBox();

    if(transformed)
        return boundingBox_.transform(getTransformationMatrix());
    else
        return boundingBox_;
}

template <class M>
void TriangleMeshGeometryCollection<M>::invalidate() {
    boundingBoxValid_ = false;
    numTrianglesValid_ = false;
}

template <class M>
void TriangleMeshGeometryCollection<M>::render() const{
    for(size_t i = 0; i < meshes_.size() ; ++i){
        meshes_[i]->render();
    }
}

template <class M>
void TriangleMeshGeometryCollection<M>::updateBoundingBox() const{
    boundingBox_ = tgt::Bounds();
    for(size_t i = 0; i < meshes_.size() ; ++i){
        boundingBox_.addVolume(meshes_[i]->getBoundingBox());
    }
    boundingBoxValid_ = true;
}

//--------------------------------------------------------------------------------------------------------------------

/// A Collection of triangle meshes with vertex type VertexBase
class VRN_CORE_API TriangleMeshGeometryCollectionSimple : public TriangleMeshGeometryCollection<VertexBase>{
    virtual Geometry* create() const { return new TriangleMeshGeometryCollectionSimple();}
    virtual std::string getClassName() const { return "TriangleMeshGeometryCollectionSimple";}
    virtual TriangleMeshGeometryBase::VertexLayout getVertexLayout() const { return SIMPLE;}
};

//--------------------------------------------------------------------------------------------------------------------

/// A Collection of triangle meshes with vertex type VertexNormal
class VRN_CORE_API TriangleMeshGeometryCollectionVec3 : public TriangleMeshGeometryCollection<VertexNormal>{
    virtual Geometry* create() const { return new TriangleMeshGeometryCollectionVec3();}
    virtual std::string getClassName() const { return "TriangleMeshGeometryCollectionVec3"; }
    virtual TriangleMeshGeometryBase::VertexLayout getVertexLayout() const { return NORMAL; }
};

//--------------------------------------------------------------------------------------------------------------------

/// A Collection of triangle meshes with vertex type VertexColorNormal
class VRN_CORE_API TriangleMeshGeometryCollectionColorNormal : public TriangleMeshGeometryCollection<VertexColorNormal>{
    virtual Geometry* create() const { return new TriangleMeshGeometryCollectionColorNormal();}
    virtual std::string getClassName() const { return "TriangleMeshGeometryCollectionColorNormal"; }
    virtual TriangleMeshGeometryBase::VertexLayout getVertexLayout() const { return COLOR_NORMAL; }
};


}

#endif //VRN_TRIANGLECLIPMESH_H
