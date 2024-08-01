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

#ifndef VRN_TWOMANIFOLDGEOMETRYCLIPPING_H
#define VRN_TWOMANIFOLDGEOMETRYCLIPPING_H

#include "voreen/core/processors/processor.h"

#include "voreen/core/ports/geometryport.h"
#include "voreen/core/properties/boolproperty.h"
#include "voreen/core/properties/floatproperty.h"
#include "voreen/core/properties/shaderproperty.h"
#include "voreen/core/properties/vectorproperty.h"
#include "voreen/core/properties/colorproperty.h"
#include "voreen/core/datastructures/geometry/glmeshgeometry.h"

#include "tgt/shadermanager.h"

#include "../util/twomanifoldclippingutil.h"
#include "../util/loopconstructor.h"
#include "../util/polygoncreator.h"
#include "../util/polygontriangulator.h"

#include "tgt/stopwatch.h"

namespace voreen {

class VRN_CORE_API TwoManifoldGeometryClipping : public Processor {
public:
    TwoManifoldGeometryClipping();
    virtual ~TwoManifoldGeometryClipping();
    virtual Processor* create() const;

    virtual std::string getCategory() const   { return "Geometry";              }
    virtual std::string getClassName() const  { return "TwoManifoldGeometryClipping";   }
    virtual CodeState getCodeState() const    { return CODE_STATE_EXPERIMENTAL; }

protected:
    virtual void setDescriptions() {
        setDescription("Clips the input geometry by a plane on the GPU and afterwards adds cap polygons to the description of the mesh to close the geometry. \
                For a input triangle mesh that is a closed two-manifold, the output of the processor is again a closed two-manifold mesh. \
                The processor requires OpenGL version 4.3 for the clipping process.");
        onlyTriangulateTopLevel_.setDescription("For input meshes that are not closed two-manifolds, this heuristic might yield aceptable results when closing the mesh.");
        setFixedClipColor_.setDescription("If the input mesh contains color information, the color of the cap polygons can be set to a fixed value.");
        clipColor_.setDescription("Color of the cap polygons.");
    }

    virtual void process();

    virtual void initialize();

    virtual void deinitialize();

    /// updates defines and transform feedback varyings in the shader program according to the geometry
    virtual void updateShader(const GlMeshGeometryBase* mesh);

    /// update the size of each vertex in bytes for transform feedback
    virtual void updateSizePerVertex(const GlMeshGeometryBase* mesh);

    /// allocates memory for the xfb buffer objects and sets the transform feedback object state
    virtual void setupTransformFeedbackSizes();

    /// deletes the buffer objects and generates them again to free the GPU memory
    virtual void freeBufferMemory();

    template<typename T>
    GlMeshGeometry<uint32_t, T>* buildOutputGeometry(const GlMeshGeometryBase* inputMesh) {

        float epsilon = epsilon_.get() * 0.0001f; 
        //build a new output mesh and copy the trafo matrix and texture data
        GlMeshGeometry<uint32_t, T>* mesh = createEmptyMesh<T>(inputMesh);

        //complete mesh has been clipped: output empty geometry
        if (writtenMeshVertices_ == 0 && writtenEdgeVertices_ == 0)
            return mesh;
        
        // retrieve the clipped geometry and put it into the output mesh
        // TODO: this is inefficient and should be changed as soon as GPU geometry classes are available
        if (clipMeshVbo_ && writtenMeshVertices_ && mesh)
            buildTriangleMeshFromBufferObject<T>(clipMeshVbo_, writtenMeshVertices_, mesh);
        else 
            return mesh;
        
        // clipped mesh is ready, now add cap mesh
        // 1. download the edges
        std::vector<ClipEdge<T> > edges;
        if (edgeVbo_ && writtenEdgeVertices_)
            retrieveClipEdges<T>(edgeVbo_, writtenEdgeVertices_, edges);
        else 
            return mesh;

        LDEBUG("Retrieved " << edges.size() << " clip edges");

        tgt::Stopwatch bmTimer;
        bmTimer.start();

        // 2. construct loops
        SimpleLoopConstructor<T> lc(/*epsilon*/0.f);
        std::list<Loop<T> > loops;
        lc.constructLoops(edges,loops);

        bmTimer.stop();
        LDEBUG("Loop construction time: " << bmTimer.getRuntime() << " ms");
        bmTimer.reset();

        LDEBUG("Found " << loops.size() << " Loops");

        bmTimer.start();
        // 3. create polygons
        SimplePolygonCreator<T> polygonCreator(onlyTriangulateTopLevel_.get());
        std::list<Loop<T>* > polygons;
        polygonCreator.createPolygons(loops, polygons);

        bmTimer.stop();
        LDEBUG("Polygon creation time: " << bmTimer.getRuntime() << " ms");
        bmTimer.reset();
        
        
        //TODO: set fixed color, fixed or normalized texture coordinates, or fixed texture indices if necessary
        if (setFixedClipColor_.get())
            setFixedColors<T>(polygons);

        bmTimer.start();
        
        // 4. triangulate polygons and add to mesh
        EarClippingTriangulator<T> triangulator(epsilon);
        
        for(typename std::list<Loop<T>* >::iterator i = polygons.begin(); i != polygons.end(); ++i) 
            if ((*i)->ccw_)
                triangulator.triangulate(*(*i), mesh);
        
        bmTimer.stop();
        LDEBUG("Triangulation time: " << bmTimer.getRuntime() << " ms");
        bmTimer.reset();

        return mesh;
    }

    /// do nothing for meshes without colors -> overwrite for color vertex classes
    template<typename T>
    void setFixedColors(std::list<Loop<T>* >& polygons) {
        return;
    }

    /// overwrite for specialized subclasses
    template<typename T>
    GlMeshGeometry<uint32_t, T>* createEmptyMesh(const GlMeshGeometryBase* inputMesh) {
        return 0;
    }

    /// downloads the data of the vertex buffer, interprets it as triangles, and adds those triangles to the mesh (overwrite for specialized subclasses)
    template<typename T>
    void buildTriangleMeshFromBufferObject(GLuint vbo, size_t meshVertices, GlMeshGeometry<uint32_t, T>* mesh) {
        return;
    }

    /**
     * Retrieves the clip edges from a vertex buffer object and appends them to the vector (overwrite for specialized subclasses).
     * Caution: vector should be empty, might otherwise lead to problems in later steps!
     */
    template<typename T>
    void retrieveClipEdges(GLuint edgeVbo, size_t edgeVertices, std::vector<ClipEdge<T> >& edges) {
        edges.clear();
    }

private:
    GeometryPort inport_;
    GeometryPort outport_;

    BoolProperty enableClipping_;       ///< disables or enables the clipping
    FloatVec3Property planeNormal_;     ///< normal of the clipping plane (in world space)
    FloatProperty planePosition_;       ///< position of the clipping plane (in world space)

    FloatProperty epsilon_;             ///< configure the epsilon for loop construction and polygon triangulation

    BoolProperty setFixedClipColor_;    ///< if set to true, the clipped areas are set to a fixed color   
    ColorProperty clipColor_;           ///< specifies the color for the clipped areas

    BoolProperty onlyTriangulateTopLevel_;  ///< "robust" mode, only triangulate top level of loop hierarchy ?

    BoolProperty alwaysFreeMemory_;     ///< if set, the OpenGL buffer objects are destroyed at the end of the process()-method to free the GPU memory



    tgt::Shader* program_;          ///< clipping shader program
    GLuint clipMeshVbo_;            ///< vertex buffer object for recording clipped mesh using transform feedback
    GLuint edgeVbo_;                ///< vertex buffer object for recording clip edges using transform feedback
    GLuint tfo_;                    ///< transform feedback object

    VertexBase::VertexLayout currentVertexLayout_;    ///< vertex layout of the current input geometry (needed for setting shader defines)
    size_t sizePerVertex_; ///< size per transform feedback vertex in bytes (for edge vertices, 2* sizeof(GLfloat) has to be added for the 2d position)
    size_t numInputTriangles_; ///< number of triangles in the input geometry

    GLuint clipQuery_;              ///< transform feedback query for retrieving the number of written vertices in the clip mesh
    GLuint edgeQuery_;              ///< transform feedback query for retrieving the number of written vertices for the clip edges

    GLuint writtenMeshVertices_;    ///< number of vertices written during transform feedback for the clipped mesh
    GLuint writtenEdgeVertices_;    ///< number of vertices written during transform feedback for the edges

    tgt::vec3 capNormal_;           ///< inverted normal of the clipping plane in model space

    GlMeshGeometryBase* outputGeometry_;    ///< the output geometry

    static const std::string loggerCat_; ///< category used in logging
};

}   //namespace

#endif // VRN_GPUGEOMETRYCLIPPING_H
