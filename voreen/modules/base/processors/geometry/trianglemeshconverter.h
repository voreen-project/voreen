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

#ifndef VRN_TRIANGLEMESHCONVERTER_H
#define VRN_TRIANGLEMESHCONVERTER_H

#include "voreen/core/processors/processor.h"

#include "voreen/core/datastructures/geometry/trianglemeshgeometry.h"
#include "voreen/core/datastructures/geometry/trianglemeshgeometryindexed.h"

#include "voreen/core/ports/geometryport.h"

#include "voreen/core/properties/vectorproperty.h"
#include "voreen/core/properties/boolproperty.h"
#include "voreen/core/properties/optionproperty.h"
#include "voreen/core/properties/colorproperty.h"
#include "voreen/core/properties/stringproperty.h"

namespace voreen {

/**
 * Can be used to transform between different types of triangle mesh.
 */
class VRN_CORE_API TriangleMeshConverter : public Processor {
public:
    TriangleMeshConverter();
    virtual ~TriangleMeshConverter();
    virtual Processor* create() const;

    virtual std::string getClassName() const { return "TriangleMeshConverter"; }
    virtual std::string getCategory() const  { return "Geometry";         }
    virtual CodeState getCodeState() const   { return CODE_STATE_STABLE;  }

protected:
    virtual void setDescriptions() {
        setDescription("Can be used to transform between different types of triangle mesh");
    }

    virtual void process();

    virtual void adjustPropertiesToInput();

    void updatePropVisibility();

    GeometryPort inport_;        ///< Inport for a triangle mesh geometry to convert.
    GeometryPort outport_;       ///< Outport for a triangle mesh geometry that was converted.

    BoolProperty enabled_;       ///< Determines whether the conversion is performed.

    StringProperty inputMesh_;          ///< Used as output to show the input geometry vertex type.
    OptionProperty<TriangleMeshGeometryBase::VertexLayout> targetMesh_;   ///< Into what kind of mesh should the source mesh be converted?

    ColorProperty defaultColor_;        ///< Used to user-define the color of converted vertices

    /// category used in logging
    static const std::string loggerCat_;

private:
    template<class I>
    Geometry* switchCaseTriangleMesh(const TriangleMeshGeometry<I>* input);
    template<class I, class O>
    void convertTriangleMeshGeometry(const TriangleMeshGeometry<I>* inputGeometry, TriangleMeshGeometry<O>* outputGeometry);

    template<class I_Vertex>
    Geometry* switchCaseTriangleMeshUInt16Indexed(const TriangleMeshGeometryIndexed<uint16_t, I_Vertex>* input);
    template<class I_Vertex>
    Geometry* switchCaseTriangleMeshUInt32Indexed(const TriangleMeshGeometryIndexed<uint32_t, I_Vertex>* input);
    template<class Index, class I_Vertex, class O_Vertex>
    void convertTriangleMeshGeometryIndexed(const TriangleMeshGeometryIndexed<Index, I_Vertex>* inputGeometry, TriangleMeshGeometryIndexed<Index, O_Vertex>* outputGeometry);
};

template<class I>
Geometry* TriangleMeshConverter::switchCaseTriangleMesh(const TriangleMeshGeometry<I>* input) {
    switch(targetMesh_.getValue()) {
        case TriangleMeshGeometryBase::SIMPLE: {
            TriangleMeshGeometrySimple* out = new TriangleMeshGeometrySimple();
            convertTriangleMeshGeometry<I,VertexBase>(input,out);
            return out; }
            break;
        case TriangleMeshGeometryBase::NORMAL: {
            TriangleMeshGeometryNormal* out = new TriangleMeshGeometryNormal();
            convertTriangleMeshGeometry<I,VertexNormal>(input,out);
            return out; }
            break;
        case TriangleMeshGeometryBase::COLOR: {
            TriangleMeshGeometryColor* out = new TriangleMeshGeometryColor();
            convertTriangleMeshGeometry<I,VertexColor>(input,out);
            return out; }
            break;
        case TriangleMeshGeometryBase::TEXCOORD: {
            TriangleMeshGeometryTexCoord* out = new TriangleMeshGeometryTexCoord();
            convertTriangleMeshGeometry<I,VertexTexCoord>(input,out);
            return out; }
            break;
        case TriangleMeshGeometryBase::COLOR_NORMAL: {
            TriangleMeshGeometryColorNormal* out = new TriangleMeshGeometryColorNormal();
            convertTriangleMeshGeometry<I,VertexColorNormal>(input,out);
            return out; }
            break;
        case TriangleMeshGeometryBase::NORMAL_TEXCOORD: {
            TriangleMeshGeometryNormalTexCoord* out = new TriangleMeshGeometryNormalTexCoord();
            convertTriangleMeshGeometry<I,VertexNormalTexCoord>(input,out);
            return out; }
            break;
        case TriangleMeshGeometryBase::COLOR_TEXCOORD: {
            TriangleMeshGeometryColorTexCoord* out = new TriangleMeshGeometryColorTexCoord();
            convertTriangleMeshGeometry<I,VertexColorTexCoord>(input,out);
            return out; }
            break;
        case TriangleMeshGeometryBase::COLOR_NORMAL_TEXCOORD: {
            TriangleMeshGeometryColorNormalTexCoord* out = new TriangleMeshGeometryColorNormalTexCoord();
            convertTriangleMeshGeometry<I,VertexColorNormalTexCoord>(input,out);
            return out; }
        default:
            LWARNING("Unknown source triangle mesh layout; could not convert.");
            break;
    }
    //should not get here
    return 0;
}

template<class I, class O>
void TriangleMeshConverter::convertTriangleMeshGeometry(const TriangleMeshGeometry<I>* inputGeometry, TriangleMeshGeometry<O>* outputGeometry) {
    tgtAssert(inputGeometry,"no input geometry");
    tgtAssert(outputGeometry,"no output geometry");
    for(size_t i = 0; i < inputGeometry->getNumTriangles(); i++) {
        Triangle<I> inTriangle = inputGeometry->getTriangle(i);
        Triangle<O> outTriangle(O(inTriangle.v_[0].pos_,(inputGeometry->supportsColors() ? inTriangle.v_[0].getColor() : defaultColor_.get()),
                                                inTriangle.v_[0].getNormal(),inTriangle.v_[0].getTexCoord(),inTriangle.v_[0].getTexIndex()),
                                O(inTriangle.v_[1].pos_,(inputGeometry->supportsColors() ? inTriangle.v_[1].getColor() : defaultColor_.get()),
                                                inTriangle.v_[1].getNormal(),inTriangle.v_[1].getTexCoord(),inTriangle.v_[1].getTexIndex()),
                                O(inTriangle.v_[2].pos_,(inputGeometry->supportsColors() ? inTriangle.v_[2].getColor() : defaultColor_.get()),
                                                inTriangle.v_[2].getNormal(),inTriangle.v_[2].getTexCoord(),inTriangle.v_[2].getTexIndex()));
        outputGeometry->addTriangle(outTriangle);
    }
    //calc normals if not set from input
    if(!inputGeometry->supportsNormals())
        outputGeometry->calculateNormals();
}

template<class I_Vertex>
Geometry* TriangleMeshConverter::switchCaseTriangleMeshUInt16Indexed(const TriangleMeshGeometryIndexed<uint16_t, I_Vertex>* input) {
    switch(targetMesh_.getValue()) {
        case TriangleMeshGeometryBase::SIMPLE: {
                TriangleMeshGeometryUInt16IndexedSimple* out = new TriangleMeshGeometryUInt16IndexedSimple();
                convertTriangleMeshGeometryIndexed<uint16_t,I_Vertex,VertexBase>(input,out);
                return out;
            } break;
        case TriangleMeshGeometryBase::NORMAL: {
                TriangleMeshGeometryUInt16IndexedNormal* out = new TriangleMeshGeometryUInt16IndexedNormal();
                convertTriangleMeshGeometryIndexed<uint16_t,I_Vertex,VertexNormal>(input,out);
                return out;
          } break;
        case TriangleMeshGeometryBase::COLOR: {
                TriangleMeshGeometryUInt16IndexedColor* out = new TriangleMeshGeometryUInt16IndexedColor();
                convertTriangleMeshGeometryIndexed<uint16_t,I_Vertex,VertexColor>(input,out);
                return out;
            } break;
        case TriangleMeshGeometryBase::TEXCOORD: {
                TriangleMeshGeometryUInt16IndexedTexCoord* out = new TriangleMeshGeometryUInt16IndexedTexCoord();
                convertTriangleMeshGeometryIndexed<uint16_t,I_Vertex,VertexTexCoord>(input,out);
                return out;
            } break;
        case TriangleMeshGeometryBase::COLOR_NORMAL: {
                TriangleMeshGeometryUInt16IndexedColorNormal* out = new TriangleMeshGeometryUInt16IndexedColorNormal();
                convertTriangleMeshGeometryIndexed<uint16_t,I_Vertex,VertexColorNormal>(input,out);
                return out;
            } break;
        case TriangleMeshGeometryBase::NORMAL_TEXCOORD: {
                TriangleMeshGeometryUInt16IndexedNormalTexCoord* out = new TriangleMeshGeometryUInt16IndexedNormalTexCoord();
                convertTriangleMeshGeometryIndexed<uint16_t,I_Vertex,VertexNormalTexCoord>(input,out);
                return out;
            } break;
        case TriangleMeshGeometryBase::COLOR_TEXCOORD: {
                TriangleMeshGeometryUInt16IndexedColorTexCoord* out = new TriangleMeshGeometryUInt16IndexedColorTexCoord();
                convertTriangleMeshGeometryIndexed<uint16_t,I_Vertex,VertexColorTexCoord>(input,out);
                return out;
            } break;
        case TriangleMeshGeometryBase::COLOR_NORMAL_TEXCOORD: {
                TriangleMeshGeometryUInt16IndexedColorNormalTexCoord* out = new TriangleMeshGeometryUInt16IndexedColorNormalTexCoord();
                convertTriangleMeshGeometryIndexed<uint16_t,I_Vertex,VertexColorNormalTexCoord>(input,out);
                return out;
            } break;
        default:
            LWARNING("Unknown source triangle mesh layout; could not convert.");
            break;
    }
    //should not get here
    return 0;
}

template<class I_Vertex>
Geometry* TriangleMeshConverter::switchCaseTriangleMeshUInt32Indexed(const TriangleMeshGeometryIndexed<uint32_t, I_Vertex>* input) {
    switch(targetMesh_.getValue()) {
        case TriangleMeshGeometryBase::SIMPLE: {
                TriangleMeshGeometryUInt32IndexedSimple* out = new TriangleMeshGeometryUInt32IndexedSimple();
                convertTriangleMeshGeometryIndexed<uint32_t,I_Vertex,VertexBase>(input,out);
                return out;
            } break;
        case TriangleMeshGeometryBase::NORMAL: {
                TriangleMeshGeometryUInt32IndexedNormal* out = new TriangleMeshGeometryUInt32IndexedNormal();
                convertTriangleMeshGeometryIndexed<uint32_t,I_Vertex,VertexNormal>(input,out);
                return out;
          } break;
        case TriangleMeshGeometryBase::COLOR: {
                TriangleMeshGeometryUInt32IndexedColor* out = new TriangleMeshGeometryUInt32IndexedColor();
                convertTriangleMeshGeometryIndexed<uint32_t,I_Vertex,VertexColor>(input,out);
                return out;
            } break;
        case TriangleMeshGeometryBase::TEXCOORD: {
                TriangleMeshGeometryUInt32IndexedTexCoord* out = new TriangleMeshGeometryUInt32IndexedTexCoord();
                convertTriangleMeshGeometryIndexed<uint32_t,I_Vertex,VertexTexCoord>(input,out);
                return out;
            } break;
        case TriangleMeshGeometryBase::COLOR_NORMAL: {
                TriangleMeshGeometryUInt32IndexedColorNormal* out = new TriangleMeshGeometryUInt32IndexedColorNormal();
                convertTriangleMeshGeometryIndexed<uint32_t,I_Vertex,VertexColorNormal>(input,out);
                return out;
            } break;
        case TriangleMeshGeometryBase::NORMAL_TEXCOORD: {
                TriangleMeshGeometryUInt32IndexedNormalTexCoord* out = new TriangleMeshGeometryUInt32IndexedNormalTexCoord();
                convertTriangleMeshGeometryIndexed<uint32_t,I_Vertex,VertexNormalTexCoord>(input,out);
                return out;
            } break;
        case TriangleMeshGeometryBase::COLOR_TEXCOORD: {
                TriangleMeshGeometryUInt32IndexedColorTexCoord* out = new TriangleMeshGeometryUInt32IndexedColorTexCoord();
                convertTriangleMeshGeometryIndexed<uint32_t,I_Vertex,VertexColorTexCoord>(input,out);
                return out;
            } break;
        case TriangleMeshGeometryBase::COLOR_NORMAL_TEXCOORD: {
                TriangleMeshGeometryUInt32IndexedColorNormalTexCoord* out = new TriangleMeshGeometryUInt32IndexedColorNormalTexCoord();
                convertTriangleMeshGeometryIndexed<uint32_t,I_Vertex,VertexColorNormalTexCoord>(input,out);
                return out;
            } break;
        default:
            LWARNING("Unknown source triangle mesh layout; could not convert.");
            break;
    }
    //should not get here
    return 0;
}

template<class Index, class I_Vertex, class O_Vertex>
void TriangleMeshConverter::convertTriangleMeshGeometryIndexed(const TriangleMeshGeometryIndexed<Index, I_Vertex>* inputGeometry, TriangleMeshGeometryIndexed<Index, O_Vertex>* outputGeometry) {
    tgtAssert(inputGeometry,"no input geometry");
    tgtAssert(outputGeometry,"no output geometry");
    //set vertices
    std::vector<O_Vertex> vec;
    std::vector<I_Vertex> vertices = inputGeometry->getVertices();
    vec.resize(vertices.size());
    for(size_t i = 0; i < vertices.size(); i++) {
        vec[i] = O_Vertex(vertices[i].pos_,(inputGeometry->supportsColors() ? vertices[i].getColor() : defaultColor_.get()),
                          vertices[i].getNormal(),vertices[i].getTexCoord(),vertices[i].getTexIndex());
    }
    outputGeometry->setVertices(vec);
    //copy triangles
    for(size_t i = 0; i < inputGeometry->getNumTriangles(); i++) {
        outputGeometry->addTriangle(inputGeometry->getTriangle(i));
    }
    //calc normals if not set from input
    if(!inputGeometry->supportsNormals())
        outputGeometry->calculateNormals();
}


} //namespace

#endif // VRN_TRIANGLEMESHCONVERTER_H
