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

#ifndef VRN_GEOMETRYSLICERENDERER_H
#define VRN_GEOMETRYSLICERENDERER_H

#include "voreen/core/processors/renderprocessor.h"

#include "voreen/core/ports/geometryport.h"
#include "voreen/core/ports/volumeport.h"
#include "voreen/core/ports/renderport.h"

#include "voreen/core/properties/boolproperty.h"
#include "voreen/core/properties/optionproperty.h"
#include "voreen/core/properties/colorproperty.h"
#include "voreen/core/properties/intproperty.h"
#include "voreen/core/properties/matrixproperty.h"

#include "voreen/core/datastructures/geometry/glmeshgeometry.h"

#include "tgt/shadermanager.h"

#include "voreen/core/datastructures/volume/slice/slicehelper.h"

namespace voreen {

/**
 * This processor renders the intersection between the plane defined by a SliceViewer and a GlMeshGeometry input.
 *
 * It requires OpenGL 4.3 since the intersection is computed by a geometry shader and recorded using transform feedback.
 */
class VRN_CORE_API GeometrySliceRenderer : public RenderProcessor {

public:

    GeometrySliceRenderer();
    //~GeometrySliceRenderer();

    virtual Processor* create() const;

    virtual std::string getClassName() const { return "GeometrySliceRenderer"; }
    virtual std::string getCategory() const { return "Geometry"; }
    virtual Processor::CodeState getCodeState() const { return CODE_STATE_TESTING; }

protected:
    virtual void setDescriptions() {
        setDescription("Renders the intersection of an input geometry consisting of triangles and a plane which is defined by a slice view. Requires OpenGL 4.3 features since the intersection is computed on the GPU.");
    }

    virtual void initialize();
    virtual void deinitialize();

    virtual void process();

    virtual void adjustPropertyVisibility();

    void updateShader(const GlMeshGeometryBase* mesh);

    void setupTransformFeedbackSizes();

    void freeBufferMemory();

    /// just for debugging purposes
    //void retrieveClipEdges(GLuint edgeVbo, size_t edgeVertices, std::vector<tgt::vec4>& edges);

    GeometryPort geometryInport_;       ///< input geometry
    VolumePort volumeInport_;   ///< volume is necessary when using the processor in combination with a slice viewer to obtain the voxelToWorld-matrix

    RenderPort renderOutport_;      ///< rendering output for overlaying with slice viewer

    BoolProperty enable_;

    ColorProperty solidColor_;  ///< color for geometries without vertex color

    /// Property containing the available alignments: xy (axial), xz (coronal), yz (sagittal)
    OptionProperty<SliceAlignment> sliceAlignment_;
    IntProperty sliceIndex_;                ///< Property containing the currently selected slice
    FloatMat4Property pickingMatrix_;       ///< Picking matrix from SliceViewer

    BoolProperty alwaysFreeMemory_;     ///< if set, the OpenGL buffer objects are destroyed at the end of the process()-method to free the GPU memory

    // OpenGL resources
    // shader
    tgt::Shader* clipShader_;       ///< shader for computing the intersection
    tgt::Shader* renderShader_;     ///< shader for rendering the geometry
    // objects
    GLuint vertexArrayID_;
    GLuint edgeVbo_;                ///< vertex buffer object for recording clip edges using transform feedback
    GLuint tfo_;                    ///< transform feedback object
    // meta information
    VertexBase::VertexLayout currentVertexLayout_;    ///< vertex layout of the current input geometry (needed for setting shader defines)
    //size_t sizePerVertex_;                                    ///< size per transform feedback vertex in bytes
    size_t numInputTriangles_;                                ///< number of triangles in the input geometry
    // query
    GLuint edgeQuery_;              ///< transform feedback query for retrieving the number of written vertices for the clip edges
    GLuint writtenEdges_;           ///< number of egdes written during transform feedback
};

}   // namespace

#endif
