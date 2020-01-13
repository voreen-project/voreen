/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2020 University of Muenster, Germany,                        *
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

#include "trianglemeshconverter.h"
#include "voreen/core/datastructures/geometry/trianglemeshgeometry.h"
#include "voreen/core/datastructures/callback/memberfunctioncallback.h"

namespace voreen {

const std::string TriangleMeshConverter::loggerCat_("voreen.base.TriangleMeshConverter");

TriangleMeshConverter::TriangleMeshConverter()
    : Processor()
    , inport_(Port::INPORT, "geometry.geometry", "Geometry Input")
    , outport_(Port::OUTPORT, "geometry.clippedgeometry", "Converted Geometry Output")
    , enabled_("enabled", "Enabled", false)
    , inputMesh_("inputMesh", "Input Type:","no input")
    , targetMesh_("targetmesh", "Convert into...", Processor::INVALID_RESULT)
    , defaultColor_("defaultColor", "Geometry Color")
{
    addPort(inport_);
    addPort(outport_);

    addProperty(enabled_);

    addProperty(inputMesh_);
        inputMesh_.setEditable(false);
    addProperty(targetMesh_);
        targetMesh_.addOption("simple", "Simple", TriangleMeshGeometryBase::SIMPLE);
        targetMesh_.addOption("vertexvec3", "VertexNormal", TriangleMeshGeometryBase::NORMAL);
        targetMesh_.addOption("vertexcolor", "VertexColor", TriangleMeshGeometryBase::COLOR);
        targetMesh_.addOption("vertextexcoord", "VertexTexCoord", TriangleMeshGeometryBase::TEXCOORD);
        targetMesh_.addOption("vertexvec4vec3", "VertexColorNormal", TriangleMeshGeometryBase::COLOR_NORMAL);
        targetMesh_.addOption("vertexnormaltexcoord", "VertexNormalTexCoord", TriangleMeshGeometryBase::NORMAL_TEXCOORD);
        targetMesh_.addOption("vertexcolortexcoord", "VertexColorTexCoord", TriangleMeshGeometryBase::COLOR_TEXCOORD);
        targetMesh_.addOption("vertexcolornormaltexcoord", "VertexColorNormalTexCoord", TriangleMeshGeometryBase::COLOR_NORMAL_TEXCOORD);
        targetMesh_.selectByKey("vertexvec4vec3");
    addProperty(defaultColor_);

    ON_PROPERTY_CHANGE(targetMesh_,TriangleMeshConverter,updatePropVisibility);
}

TriangleMeshConverter::~TriangleMeshConverter()
{}

Processor* TriangleMeshConverter::create() const {
    return new TriangleMeshConverter();
}

void TriangleMeshConverter::updatePropVisibility() {
    if(targetMesh_.getValue() == TriangleMeshGeometryBase::COLOR) {
        defaultColor_.setVisibleFlag(true);
    } else {
        defaultColor_.setVisibleFlag(false);
    }
}

void TriangleMeshConverter::adjustPropertiesToInput() {
    const Geometry* inputGeometry = inport_.getData();
    if (!inputGeometry) {
        inputMesh_.set("no input");
    } else {
        const TriangleMeshGeometryBase* tm = dynamic_cast<const TriangleMeshGeometryBase*>(inputGeometry);
        if (!tm) {
            inputMesh_.set("no triangle mesh");
        } else {
            switch(tm->getVertexLayout()) {
            case TriangleMeshGeometryBase::SIMPLE:
                inputMesh_.set("VertexBase");
                break;
            case TriangleMeshGeometryBase::NORMAL:
                inputMesh_.set("VertexNormal");
                break;
            case TriangleMeshGeometryBase::COLOR:
                inputMesh_.set("VertexColor");
                break;
            case TriangleMeshGeometryBase::TEXCOORD:
                inputMesh_.set("VertexTexCoord");
                break;
            case TriangleMeshGeometryBase::COLOR_NORMAL:
                inputMesh_.set("VertexColorNormal");
                break;
            case TriangleMeshGeometryBase::NORMAL_TEXCOORD:
                inputMesh_.set("VertexNormalTexCoord");
                break;
            case TriangleMeshGeometryBase::COLOR_TEXCOORD:
                inputMesh_.set("VertexColorTexCoord");
                break;
            case TriangleMeshGeometryBase::COLOR_NORMAL_TEXCOORD:
                inputMesh_.set("VertexColorNormalTexCoord");
                break;
            default:
                inputMesh_.set("unsupported type");
                tgtAssert(false,"unknown type");
            break;
            }
        }
    }
}


void TriangleMeshConverter::process() {
    const Geometry* inputGeometry = inport_.getData();
    tgtAssert(inport_.getData(), "no geometry");
    //if disabled, copy geometry into outport
    if (!enabled_.get()) {
        outport_.setData(inputGeometry, false);
        return;
    }
    //check, if we have a triangle mesh
    const TriangleMeshGeometryBase* tm = dynamic_cast<const TriangleMeshGeometryBase*>(inputGeometry);
    if (!tm) {
        LWARNING("This processor can only be used with triangle meshes; skipping execution.");
        outport_.setData(inputGeometry, false);
        return;
    }

    Geometry* outputGeometry = 0;
    //case indexed mesh
    if(const TriangleMeshGeometryIndexedBase* tmi = dynamic_cast<const TriangleMeshGeometryIndexedBase*>(tm)) {
        if(tmi->getIndexType() == TriangleMeshGeometryIndexedBase::INDEX_UINT16) {
            switch(tmi->getVertexLayout()) {
                case TriangleMeshGeometryBase::SIMPLE:
                    outputGeometry = switchCaseTriangleMeshUInt16Indexed<VertexBase>(static_cast<const TriangleMeshGeometryUInt16IndexedSimple*>(inputGeometry));
                    break;
                case TriangleMeshGeometryBase::NORMAL:
                    outputGeometry = switchCaseTriangleMeshUInt16Indexed<VertexNormal>(static_cast<const TriangleMeshGeometryUInt16IndexedNormal*>(inputGeometry));
                    break;
                case TriangleMeshGeometryBase::COLOR:
                    outputGeometry = switchCaseTriangleMeshUInt16Indexed<VertexColor>(static_cast<const TriangleMeshGeometryUInt16IndexedColor*>(inputGeometry));
                    break;
                case TriangleMeshGeometryBase::TEXCOORD:
                    outputGeometry = switchCaseTriangleMeshUInt16Indexed<VertexTexCoord>(static_cast<const TriangleMeshGeometryUInt16IndexedTexCoord*>(inputGeometry));
                    break;
                case TriangleMeshGeometryBase::COLOR_NORMAL:
                    outputGeometry = switchCaseTriangleMeshUInt16Indexed<VertexColorNormal>(static_cast<const TriangleMeshGeometryUInt16IndexedColorNormal*>(inputGeometry));
                    break;
                case TriangleMeshGeometryBase::NORMAL_TEXCOORD:
                    outputGeometry = switchCaseTriangleMeshUInt16Indexed<VertexNormalTexCoord>(static_cast<const TriangleMeshGeometryUInt16IndexedNormalTexCoord*>(inputGeometry));
                    break;
                case TriangleMeshGeometryBase::COLOR_TEXCOORD:
                    outputGeometry = switchCaseTriangleMeshUInt16Indexed<VertexColorTexCoord>(static_cast<const TriangleMeshGeometryUInt16IndexedColorTexCoord*>(inputGeometry));
                    break;
                case TriangleMeshGeometryBase::COLOR_NORMAL_TEXCOORD:
                    outputGeometry = switchCaseTriangleMeshUInt16Indexed<VertexColorNormalTexCoord>(static_cast<const TriangleMeshGeometryUInt16IndexedColorNormalTexCoord*>(inputGeometry));
                    break;
                default:
                    LWARNING("Unknown source triangle mesh layout; could not convert.");
                    break;
            }
        } else { //unit32
            switch(tmi->getVertexLayout()) {
                case TriangleMeshGeometryBase::SIMPLE:
                    outputGeometry = switchCaseTriangleMeshUInt32Indexed<VertexBase>(static_cast<const TriangleMeshGeometryUInt32IndexedSimple*>(inputGeometry));
                    break;
                case TriangleMeshGeometryBase::NORMAL:
                    outputGeometry = switchCaseTriangleMeshUInt32Indexed<VertexNormal>(static_cast<const TriangleMeshGeometryUInt32IndexedNormal*>(inputGeometry));
                    break;
                case TriangleMeshGeometryBase::COLOR:
                    outputGeometry = switchCaseTriangleMeshUInt32Indexed<VertexColor>(static_cast<const TriangleMeshGeometryUInt32IndexedColor*>(inputGeometry));
                    break;
                case TriangleMeshGeometryBase::TEXCOORD:
                    outputGeometry = switchCaseTriangleMeshUInt32Indexed<VertexTexCoord>(static_cast<const TriangleMeshGeometryUInt32IndexedTexCoord*>(inputGeometry));
                    break;
                case TriangleMeshGeometryBase::COLOR_NORMAL:
                    outputGeometry = switchCaseTriangleMeshUInt32Indexed<VertexColorNormal>(static_cast<const TriangleMeshGeometryUInt32IndexedColorNormal*>(inputGeometry));
                    break;
                case TriangleMeshGeometryBase::NORMAL_TEXCOORD:
                    outputGeometry = switchCaseTriangleMeshUInt32Indexed<VertexNormalTexCoord>(static_cast<const TriangleMeshGeometryUInt32IndexedNormalTexCoord*>(inputGeometry));
                    break;
                case TriangleMeshGeometryBase::COLOR_TEXCOORD:
                    outputGeometry = switchCaseTriangleMeshUInt32Indexed<VertexColorTexCoord>(static_cast<const TriangleMeshGeometryUInt32IndexedColorTexCoord*>(inputGeometry));
                    break;
                case TriangleMeshGeometryBase::COLOR_NORMAL_TEXCOORD:
                    outputGeometry = switchCaseTriangleMeshUInt32Indexed<VertexColorNormalTexCoord>(static_cast<const TriangleMeshGeometryUInt32IndexedColorNormalTexCoord*>(inputGeometry));
                    break;
                default:
                    LWARNING("Unknown source triangle mesh layout; could not convert.");
                    break;
            }
        }
    } else { //case normal triangle mesh
        switch(tm->getVertexLayout()) {
            case TriangleMeshGeometryBase::SIMPLE:
                outputGeometry = switchCaseTriangleMesh<VertexBase>(static_cast<const TriangleMeshGeometrySimple*>(inputGeometry));
                break;
            case TriangleMeshGeometryBase::NORMAL:
                outputGeometry = switchCaseTriangleMesh<VertexNormal>(static_cast<const TriangleMeshGeometryNormal*>(inputGeometry));
                break;
            case TriangleMeshGeometryBase::COLOR:
                outputGeometry = switchCaseTriangleMesh<VertexColor>(static_cast<const TriangleMeshGeometryColor*>(inputGeometry));
                break;
            case TriangleMeshGeometryBase::TEXCOORD:
                outputGeometry = switchCaseTriangleMesh<VertexTexCoord>(static_cast<const TriangleMeshGeometryTexCoord*>(inputGeometry));
                break;
            case TriangleMeshGeometryBase::COLOR_NORMAL:
                outputGeometry = switchCaseTriangleMesh<VertexColorNormal>(static_cast<const TriangleMeshGeometryColorNormal*>(inputGeometry));
                break;
            case TriangleMeshGeometryBase::NORMAL_TEXCOORD:
                outputGeometry = switchCaseTriangleMesh<VertexNormalTexCoord>(static_cast<const TriangleMeshGeometryNormalTexCoord*>(inputGeometry));
                break;
            case TriangleMeshGeometryBase::COLOR_TEXCOORD:
                outputGeometry = switchCaseTriangleMesh<VertexColorTexCoord>(static_cast<const TriangleMeshGeometryColorTexCoord*>(inputGeometry));
                break;
            case TriangleMeshGeometryBase::COLOR_NORMAL_TEXCOORD:
                outputGeometry = switchCaseTriangleMesh<VertexColorNormalTexCoord>(static_cast<const TriangleMeshGeometryColorNormalTexCoord*>(inputGeometry));
                break;
            default:
                LWARNING("Unknown source triangle mesh layout; could not convert.");
                break;
        }
    }

    //set trafo and output
    if(outputGeometry)
        outputGeometry->setTransformationMatrix(tm->getTransformationMatrix());

    outport_.setData(outputGeometry);
}

}  //namespace
