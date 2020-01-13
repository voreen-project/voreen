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

#include "twomanifoldclippingutil.h"

namespace voreen {

float orientationTest(const tgt::vec2& pa, const tgt::vec2& pb, const tgt::vec2& pc) {
    float acx = pa.x - pc.x;
    float bcx = pb.x - pc.x;
    float acy = pa.y - pc.y;
    float bcy = pb.y - pc.y;

    return acx * bcy - acy * bcx;
}

bool leftTurnEpsilon(const tgt::vec2& pa, const tgt::vec2& pb, const tgt::vec2& pc, float epsilon) {
    return (orientationTest(pa, pb, pc) >= epsilon);
}

bool leftTurn(const tgt::vec2& pa, const tgt::vec2& pb, const tgt::vec2& pc) {
    return (orientationTest(pa, pb, pc) > 0);
}

bool leftTurnCollinear(const tgt::vec2& pa, const tgt::vec2& pb, const tgt::vec2& pc) {
    return (orientationTest(pa, pb, pc) >= 0);
}

bool linesProperIntersect(const tgt::vec2& p0, const tgt::vec2& p1, const tgt::vec2& p2, const tgt::vec2& p3) {
    //compute bounding boxes
    tgt::vec2 r0 = tgt::vec2(std::min(p0.x, p1.x), std::min(p0.y, p1.y));
    tgt::vec2 r1 = tgt::vec2(std::max(p0.x, p1.x), std::max(p0.y, p1.y));
    tgt::vec2 r2 = tgt::vec2(std::min(p2.x, p3.x), std::min(p2.y, p3.y));
    tgt::vec2 r3 = tgt::vec2(std::max(p2.x, p3.x), std::max(p2.y, p3.y));

    return ((r1.x >= r2.x) && (r3.x >= r0.x) && (r1.y >= r2.y) && (r3.y >= r0.y)
            && (orientationTest(p0, p2, p3) * orientationTest(p1, p2, p3) < 0.f)
            && (orientationTest(p2, p0, p1) * orientationTest(p3, p0, p1) < 0.f));
}

bool linesIntersect(const tgt::vec2& p0, const tgt::vec2& p1, const tgt::vec2& p2, const tgt::vec2& p3) {
    //compute bounding boxes
    tgt::vec2 r0 = tgt::vec2(std::min(p0.x, p1.x), std::min(p0.y, p1.y));
    tgt::vec2 r1 = tgt::vec2(std::max(p0.x, p1.x), std::max(p0.y, p1.y));
    tgt::vec2 r2 = tgt::vec2(std::min(p2.x, p3.x), std::min(p2.y, p3.y));
    tgt::vec2 r3 = tgt::vec2(std::max(p2.x, p3.x), std::max(p2.y, p3.y));

    return ((r1.x >= r2.x) && (r3.x >= r0.x) && (r1.y >= r2.y) && (r3.y >= r0.y)
            && (orientationTest(p0, p2, p3) * orientationTest(p1, p2, p3) <= 0.f)
            && (orientationTest(p2, p0, p1) * orientationTest(p3, p0, p1) <= 0.f));
}

/*
    template<typename T>

//computes the winding order of the (transformed) vertices (cw or ccw) and sets the coorresponding flag
void setLoopOrientation(Loop& loop) {
    float area = 0;
    //begin with edge to the first vertex and then iterate through all edges
    tgt::vec2 begin = loop.vertices_.back().pos2d_;
    tgt::vec2 end;
    for (std::list<ClipVertex>::iterator i = loop.vertices_.begin(); i != loop.vertices_.end(); ++i) {
        end = i->pos2d_;
        area += begin.x * end.y - begin.y * end.x;
        begin = end;
    }

    loop.ccw_ = (area > 0);
}


//triangulates the loop by naive(!) ear-clipping - n^3 method for n vertices.
//vertices are removed from the loop every time a triangle is added, triangles is assumed to have enough memory already!!!
//uses the positions of the vertices that have been transformed to 2d
bool triangulateLoop(Loop& loop, std::vector<ClipVertex>& triangles) {
    size_t trianglesIndex = 0;

    //outer loop: repeat until fully triangulated
    while (loop.vertices_.size() > 2) {
        //second loop: loop over triangles
        ClipVertex first = loop.vertices_.back();
        ClipVertex second = loop.vertices_.front();
        std::list<ClipVertex>::iterator currentVertex = loop.vertices_.begin();
        currentVertex++;
        bool found = false;
        size_t checkedVertices = 0;     // for robustness: prevent endless loop
        while (!found) {
            ClipVertex third = *currentVertex;

            //transform triangle positions into xy-plane
            tgt::vec2 firstPos = first.pos2d_;
            tgt::vec2 secondPos = second.pos2d_;
            tgt::vec2 thirdPos = third.pos2d_;
            
            //check if the second vertex is convex or a reflex vertex
            if (leftTurn(firstPos, secondPos, thirdPos)) {
                bool inside = false;
                //convex vertex -> now check if any of the other points lies within the triangle
                for (std::list<ClipVertex>::iterator checkPoint = loop.vertices_.begin(); checkPoint != loop.vertices_.end(); ++checkPoint) {
                    // no need to account for degenerated edges, as those should have been removed during loop construction
                    if (!(*checkPoint == first) && !(*checkPoint == second) && !(*checkPoint == third)) {
                        tgt::vec2 checkPos = checkPoint->pos2d_;
                        if (leftTurn(firstPos, secondPos, checkPos) && leftTurn(secondPos, thirdPos, checkPos) && leftTurn(thirdPos, firstPos, checkPos)) {
                            //found point inside -> cannot clip triangle
                            inside = true;
                            break;
                        }
                    }
                }

                if (!inside) {
                    found = true;
                    //convex and no point inside -> clip triangle :-)
                    //triangles.at(trianglesIndex++) = first;
                    //triangles.at(trianglesIndex++) = second;
                    //triangles.at(trianglesIndex++) = third;
                    triangles.push_back(first);
                    triangles.push_back(second);
                    triangles.push_back(third);

                    //remove the point in the middle, which is the predecessor to the current vertex iterator
                    if (!(currentVertex == loop.vertices_.begin()))
                        loop.vertices_.erase(--currentVertex);
                    else 
                        loop.vertices_.erase(--loop.vertices_.end());

                    //reset the robustness counter :-)
                    checkedVertices = 0;
                }
                else
                    checkedVertices++;
            }
            else
                checkedVertices++;

            if (checkedVertices == loop.vertices_.size()) {
                //TODO: what about the data in the triangles vector?!
                std::cout << "Triangulation failed. Possible reasons: not a two-manifold, intersecting cw loops (holes), numerical problems." << std::endl;
                std::cout << "Try to use robust mode that ignores holes to clip the mesh." << std::endl;
                return false;
            }
            
            //reflex vertex or vertex contained in triangle -> check next triangle
            if (!found) {
                currentVertex++;
                if (currentVertex == loop.vertices_.end())
                    currentVertex = loop.vertices_.begin();
                first = second;
                second = third;
            }
        }
    }

    return true;
}

//
// ------- Stuff for building a hierarchy of loops before triangulation
//  

struct LoopNode {

    LoopNode(Loop* loop) : loop_(loop) { }

    bool operator==(const LoopNode& a) const {
        return (loop_ == a.loop_);
    }

    Loop* loop_;
    std::list<LoopNode> children_;
};


//
// * Checks if a point is inside polygon, given by its vertices in correct order.
// * The algorithm uses the ray crossing method described in "Computational Geometry in C", ch. 7.4.2, p. 239
// * Implementation is based on the c code from http://www.ecse.rpi.edu/Homepages/wrf/Research/Short_Notes/pnpoly.html
// *
bool pointInPoly(const ClipVertex& v, const Loop& l) {
    bool c = false;
    tgt::vec2 vPos = v.pos2d_;
    
    tgt::vec2 pFirst = l.vertices_.back().pos2d_;
    tgt::vec2 pSecond;
    for (std::list<ClipVertex>::const_iterator iter = l.vertices_.begin(); iter != l.vertices_.end(); ++iter) {
        pSecond = iter->pos2d_;

        if (((pSecond.y > vPos.y) != (pFirst.y > vPos.y)) && (vPos.x < (pFirst.x - pSecond.x) * (vPos.y - pSecond.y) / (pFirst.y - pSecond.y) + pSecond.x))
            c = !c;

        pFirst = pSecond;
    }

    return c;
}

//
// Tests if loop a lies in loop b by checking all vertices using the winding number algorithm.
// 
bool allVerticesInLoop(const Loop& a, const Loop& b) {

    for (std::list<ClipVertex>::const_iterator iter = a.vertices_.begin(); iter != a.vertices_.end(); ++iter) {
        if (!pointInPoly(*iter, b))
            return false;
    }

    return true;
}

void placeInHierarchy(std::list<LoopNode>& level, Loop* loop) {
    
    // build new node for the loop
    LoopNode node(loop); 

    //check if loop is inside the current reference loop
    for (std::list<LoopNode>::iterator iter = level.begin(); iter != level.end(); ++iter) {
        if (allVerticesInLoop(*loop, *(iter->loop_))) {
            placeInHierarchy(iter->children_, loop);
            return;
        }
    }

    //check if loop encloses one of the current loops
    std::list<LoopNode>::iterator iter = level.begin();
    while (iter != level.end()) {
        if (allVerticesInLoop(*(iter->loop_), *loop)) {
            node.children_.push_back(*iter);
            level.remove(*iter);
            iter = level.begin();
        }
        else
            ++iter;
    }
    
    level.push_back(node);
}

// for sorting loops according to the x coordinate of their first (transformed) vertex in the vertex list
bool loopXSort(const Loop* l1, const Loop* l2) {
    return (l1->vertices_.front().pos2d_.x < l2->vertices_.front().pos2d_.x);
}

Loop* makePolygon(Loop* outerLoop, std::list<Loop*>& holes) {
    if (holes.empty())
        return outerLoop;

    //1. find vertex with minimal x coordinate for each hole and choose it as a starting point for the corresponding vertex list 
    for (std::list<Loop*>::iterator holeIter = holes.begin(); holeIter != holes.end(); ++holeIter) {
        std::list<ClipVertex>::iterator minIter = (*holeIter)->vertices_.begin();
        tgt::vec2 minPos = minIter->pos2d_;
        for (std::list<ClipVertex>::iterator vertexIter = (*holeIter)->vertices_.begin(); vertexIter != (*holeIter)->vertices_.end(); ++vertexIter) {
            tgt::vec2 currentPos = vertexIter->pos2d_;
            if (currentPos.x < minPos.x) {
                minIter = vertexIter;
                minPos = currentPos;
            }
            else if ((currentPos.x == minPos.x) && (currentPos.y < minPos.y)) {
                minIter = vertexIter;
                minPos = currentPos;
            }
        }

        if (minIter != (*holeIter)->vertices_.begin()) {
            std::list<ClipVertex> tmpList;
            tmpList.splice(tmpList.begin(), (*holeIter)->vertices_, (*holeIter)->vertices_.begin(), minIter);
            (*holeIter)->vertices_.splice((*holeIter)->vertices_.end(), tmpList);
        }
    } 

    //2. sort holes according to minimal x coordinate
    holes.sort(loopXSort);

    //3. for each hole: find contour bridge and insert the vertices into the outer loop
    //we know that we have to take the first vertex of the hole and only vertices to the left of it are relevant
    for (std::list<Loop*>::iterator holeIter = holes.begin(); holeIter != holes.end(); ++holeIter) {
        //find possible starting point for contour bridge
        tgt::vec2 endPos = (*holeIter)->vertices_.front().pos2d_;
        for (std::list<ClipVertex>::iterator candidateIter = outerLoop->vertices_.begin(); candidateIter != outerLoop->vertices_.end(); ++candidateIter) {
            tgt::vec2 candidatePos = candidateIter->pos2d_;
            if (candidatePos.x < endPos.x) {

                //possible starting point -> check new edge against all edges of the outer loop
                bool intersect = false;

                //iterate over edges and perform intersection test
                tgt::vec2 e1 = outerLoop->vertices_.back().pos2d_;
                tgt::vec2 e2;
                bool lastWasCandidate = false;
                for (std::list<ClipVertex>::iterator edgeIter = outerLoop->vertices_.begin(); edgeIter != outerLoop->vertices_.end(); ++edgeIter) {
                    //TODO: now we also check for intersections with the edges that have the candidate point as starting / ending point
                    //since this always leads to an intersection, the intersection check has been modified to only return the case
                    //of a "real" intersection (touching is ok), but instead the loop here should be changed to not check against the two problematic edges
                    e2 = edgeIter->pos2d_;
    
                    if (linesIntersect(e1, e2, candidatePos, endPos)) {
                        intersect = true;
                        break;
                    }

                    e1 = e2;
                }
                
                if (!intersect) {
                    // we found the contour bridge -> duplicate the bridge vertices
                    (*holeIter)->vertices_.push_back((*holeIter)->vertices_.front());
                    (*holeIter)->vertices_.push_back(*candidateIter);    
                    outerLoop->vertices_.splice(++candidateIter, (*holeIter)->vertices_);

                    break;
                }
            }
        }
    }

    return outerLoop;
}

void buildSimplePolygons(std::vector<Loop*>& polygons, std::list<LoopNode>& level) {
    for (std::list<LoopNode>::iterator iter = level.begin(); iter != level.end(); ++iter) {
        if (!(iter->loop_->ccw_)) {
            //current loop is a hole and should be handled above / not triangulated -> traverse child nodes
            buildSimplePolygons(polygons, iter->children_);
        }
        else {
            //current loop is ccw -> find holes
            std::list<Loop*> cwLoops;
            for (std::list<LoopNode>::iterator childIter = iter->children_.begin(); childIter != iter->children_.end(); ++childIter) 
                if (!(childIter->loop_->ccw_))
                    cwLoops.push_back(childIter->loop_);

            polygons.push_back(makePolygon(iter->loop_, cwLoops));
            
            //traverse children -> for cw loops, only children are traversed, for other loops triangulation is executed
            buildSimplePolygons(polygons, iter->children_);
        }
    }
}*/

} //namespace
