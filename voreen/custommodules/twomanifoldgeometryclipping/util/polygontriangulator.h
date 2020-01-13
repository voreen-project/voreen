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

#ifndef VRN_GPUGEOMETRYCLIPPING_POLYGONTRIANGULATOR
#define VRN_GPUGEOMETRYCLIPPING_POLYGONTRIANGULATOR

#include "voreen/core/datastructures/geometry/trianglemeshgeometry.h"
#include "voreen/core/datastructures/geometry/vertex.h"

#include "twomanifoldclippingutil.h"

#include "ext/tgt/logmanager.h"

namespace voreen {

/**
 * Abstract base class for triangulator classes.
 */  
template<typename T>
class PolygonTriangulatorBase {

public:

    virtual void triangulate(Loop<T>& polygon, GlMeshGeometry<uint32_t, T>* mesh) = 0;    
   
};  

/**
 * Simple ear-clipping implementation.
 * Runs in O(n^2) and always selects ear with minimal inner angle to improve triangulation quality.
 * Prior to the ear-clipping process an O(n) convexity test is performed and an O(n) triangulation is used for convex polygons.
 */
template<typename T>
class EarClippingTriangulator : public PolygonTriangulatorBase<T> 
{

public:

    EarClippingTriangulator(float epsilon = 0.f)
        : PolygonTriangulatorBase<T>()
        , epsilon_(epsilon)
    {  }      
    
    virtual void triangulate(Loop<T>& polygon, GlMeshGeometry<uint32_t, T>* mesh) {

        LDEBUGC("voreen.gpugeometryclipping.earclippingtriangulator", "Triangulating polygon with " << polygon.vertices_.size() << " vertices");
         
        //perform O(n) convexity test (and remember the largest interior angle)
        typename std::list<ClipVertex<T> >::iterator firstVertex = polygon.vertices_.end(); --firstVertex; --firstVertex;
        typename std::list<ClipVertex<T> >::iterator secondVertex = polygon.vertices_.end(); --secondVertex;
        typename std::list<ClipVertex<T> >::iterator thirdVertex = polygon.vertices_.begin();
        bool convex = true;
        float largestAngleCos = 2.f;
        typename std::list<ClipVertex<T> >::iterator largestAngleVertex;
        while (thirdVertex != polygon.vertices_.end()) {
            //check three subsequent vertices if they make a left turn
            if (!leftTurn(firstVertex->pos2d_, secondVertex->pos2d_, thirdVertex->pos2d_)) {
                convex = false;
                break;
            }
            else {
                //compute interior angle -> min cos = max angle
                float cos = tgt::dot(tgt::normalize(firstVertex->pos2d_ - secondVertex->pos2d_), tgt::normalize(thirdVertex->pos2d_ - secondVertex->pos2d_));
                if (cos < largestAngleCos) {
                    largestAngleCos = cos;
                    largestAngleVertex = secondVertex;
                }
            }

            // set to the next triple
            firstVertex++;
            if (firstVertex == polygon.vertices_.end())
                firstVertex = polygon.vertices_.begin();

            secondVertex++;
            if (secondVertex == polygon.vertices_.end())
                secondVertex = polygon.vertices_.begin();

            thirdVertex++;
        }

        //if the polygon is convex: O(n) triangulation starting from the largest angle vertex
        if (convex) {
            LDEBUGC("voreen.gpugeometryclipping.earclippingtriangulator", "Found convex polygon... using O(n) triangulation");

            //change the list so that the largest angle vertex is at the front
            if (largestAngleVertex != polygon.vertices_.begin()) {
                typename std::list<ClipVertex<T> >::iterator lastVertex = largestAngleVertex;
                lastVertex--;
                polygon.vertices_.splice(polygon.vertices_.end(), polygon.vertices_, polygon.vertices_.begin(), lastVertex);
            }

            //now insert diagonals starting at the largestAngleVertex until the polygon is triangulated
            size_t numTriangles = polygon.vertices_.size() - 2;
            while (numTriangles > 0) {
                // vertex with largest angle is always at the beginning 
                typename std::list<ClipVertex<T> >::iterator firstVertex = polygon.vertices_.begin();
                typename std::list<ClipVertex<T> >::iterator secondVertex = ++polygon.vertices_.begin();
                typename std::list<ClipVertex<T> >::iterator thirdVertex = ++(++polygon.vertices_.begin());
                //add triangle and remove middle vertex from polygon
                //mesh->addTriangle(Triangle<T>(firstVertex->vertex_, secondVertex->vertex_, thirdVertex->vertex_));
                mesh->addVertex(firstVertex->vertex_);
                mesh->addVertex(secondVertex->vertex_);
                mesh->addVertex(thirdVertex->vertex_);
                polygon.vertices_.erase(secondVertex);
                numTriangles--;
            }

            return;        
        }

        // polygon is not convex -> apply ear clipping algorithm
        // initialize ear tip status and find the first ear
        float maxCos = -2.f; 
        typename std::list<ClipVertex<T> >::iterator bestEar = polygon.vertices_.end();;
        typename std::list<ClipVertex<T> >::iterator currentVertex; 
        for (currentVertex = polygon.vertices_.begin(); currentVertex != polygon.vertices_.end(); ++currentVertex) {
            initEarTipStatus(polygon, currentVertex);
            if (currentVertex->earTipStatus_ && currentVertex->minCos_ > maxCos) {
                maxCos = currentVertex->minCos_;
                bestEar = currentVertex;
            }
        }

        if (maxCos == -2.f || bestEar == polygon.vertices_.end()) {
            LERRORC("voreen.gpugeometryclipping.earclippingtriangulator", "Triangulation failed, found no valid ear in the polygon!");
            return;
        }

        while (polygon.vertices_.size() > 3) {
            // add triangle and remove best ear
            typename std::list<ClipVertex<T> >::iterator currentVertex = bestEar;
            // set predecessor and successor vertices
            typename std::list<ClipVertex<T> >::iterator predVertex;
            if (currentVertex == polygon.vertices_.begin())
                predVertex = polygon.vertices_.end();
            else
                predVertex = currentVertex;
            --predVertex;

            typename std::list<ClipVertex<T> >::iterator succVertex = currentVertex;
            ++succVertex;
            if (succVertex == polygon.vertices_.end())
                succVertex = polygon.vertices_.begin();

            //mesh->addTriangle(Triangle<T>(predVertex->vertex_, currentVertex->vertex_, succVertex->vertex_));
            mesh->addVertex(predVertex->vertex_);
            mesh->addVertex(currentVertex->vertex_);
            mesh->addVertex(succVertex->vertex_);
            succVertex = polygon.vertices_.erase(currentVertex);
            if (succVertex == polygon.vertices_.end()) {
                succVertex = polygon.vertices_.begin();
                predVertex = polygon.vertices_.end();
                --predVertex;
            }
            else if (succVertex == polygon.vertices_.begin()){
                predVertex = polygon.vertices_.end();
                --predVertex;
            }
            else {
                predVertex = succVertex;
                --predVertex;
            }

            // check ear tip status of predecessor and successor again
            initEarTipStatus(polygon, predVertex);
            initEarTipStatus(polygon, succVertex);

            //find new best ear
            bestEar = polygon.vertices_.end();
            float maxCos = -2.f;
            for (currentVertex = polygon.vertices_.begin(); currentVertex != polygon.vertices_.end(); ++currentVertex) {
                if (currentVertex->earTipStatus_ && currentVertex->minCos_ > maxCos) {
                    maxCos = currentVertex->minCos_;
                    bestEar = currentVertex;
                }
            }

            if (maxCos == -2.f || bestEar == polygon.vertices_.end()) {
                LERRORC("voreen.gpugeometryclipping.earclippingtriangulator", "Triangulation failed, found no ear anymore during triangulation!");
                return;
            }

        }

        //write out last triangle
        currentVertex = polygon.vertices_.begin();
        T v0 = currentVertex->vertex_;
        currentVertex++;
        T v1 = currentVertex->vertex_;
        currentVertex++;
        T v2 = currentVertex->vertex_;

        //mesh->addTriangle(Triangle<T>(v0,v1,v2));
        mesh->addVertex(v0);
        mesh->addVertex(v1);
        mesh->addVertex(v2);

        polygon.vertices_.clear(); 
    }

protected:

    /// set the ear tip status (including the cosine of the max interior angles if the vertex is a convex vertex)
    virtual void initEarTipStatus(Loop<T>& polygon, typename std::list<ClipVertex<T> >::iterator& currentVertex) {
        // set predecessor and successor vertices
        typename std::list<ClipVertex<T> >::iterator predVertex;
        if (currentVertex == polygon.vertices_.begin())
            predVertex = polygon.vertices_.end();
        else
            predVertex = currentVertex;
        --predVertex;

        typename std::list<ClipVertex<T> >::iterator succVertex = currentVertex;
        ++succVertex;
        if (succVertex == polygon.vertices_.end())
            succVertex = polygon.vertices_.begin();

        currentVertex->earTipStatus_ = false;

        //check if currentVertex is convex (include collinear case for degenerated cases)
        if (epsilon_ == 0.f) {
            if (!leftTurnCollinear(predVertex->pos2d_, currentVertex->pos2d_, succVertex->pos2d_)) 
                return;
        }
        else {  // use epsilon
            if (!leftTurnEpsilon(predVertex->pos2d_, currentVertex->pos2d_, succVertex->pos2d_, -epsilon_))
                return;
        }

        // convex vertex -> check triangle against all points in polygon
        for (typename std::list<ClipVertex<T> >::iterator checkPoint = polygon.vertices_.begin(); checkPoint != polygon.vertices_.end(); ++checkPoint) {
            // no need to account for degenerated edges, as those should have been removed during loop construction
            if (!(*checkPoint == *predVertex) && !(*checkPoint == *currentVertex) && !(*checkPoint == *succVertex)) {
                tgt::vec2 checkPos = checkPoint->pos2d_;
                if (epsilon_ == 0.f) {
                    if (leftTurn(predVertex->pos2d_, currentVertex->pos2d_, checkPos) && leftTurn(currentVertex->pos2d_, succVertex->pos2d_, checkPos) && leftTurn(succVertex->pos2d_, predVertex->pos2d_, checkPos)) {
                        //found point inside -> no ear tip
                        return;
                    }
                }   
                else { //use epsilon
                    if (leftTurnEpsilon(predVertex->pos2d_, currentVertex->pos2d_, checkPos, epsilon_) 
                            && leftTurnEpsilon(currentVertex->pos2d_, succVertex->pos2d_, checkPos, epsilon_) 
                            && leftTurnEpsilon(succVertex->pos2d_, predVertex->pos2d_, checkPos, epsilon_)) {
                        //found point inside -> no ear tip
                        return;
                    }
                }
            }
        }
        
        // set ear tip status and compute the cos of all interior angles to find the min (if necessary)
        currentVertex->earTipStatus_ = true;
        float cos1 = tgt::dot(tgt::normalize(succVertex->pos2d_ - predVertex->pos2d_), tgt::normalize(currentVertex->pos2d_ - predVertex->pos2d_)); 
        float cos2 = tgt::dot(tgt::normalize(predVertex->pos2d_ - currentVertex->pos2d_), tgt::normalize(succVertex->pos2d_ - currentVertex->pos2d_));
        float cos3 = tgt::dot(tgt::normalize(predVertex->pos2d_ - succVertex->pos2d_), tgt::normalize(currentVertex->pos2d_ - succVertex->pos2d_));
        currentVertex->minCos_ = std::min(cos1, std::min(cos2, cos3)); 

    }

    float epsilon_;
};

} //namespace

#endif
