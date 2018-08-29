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

#ifndef VRN_GPUGEOMETRYCLIPPING_POLYGONCREATOR
#define VRN_GPUGEOMETRYCLIPPING_POLYGONCREATOR

#include "twomanifoldclippingutil.h"

namespace voreen {

// node in the data structure that is used to hierarchically order the loops
template<typename T>
struct LoopNode {

    LoopNode(Loop<T>* loop) : loop_(loop) { }

    bool operator==(const LoopNode<T>& a) const {
        return (loop_ == a.loop_);
    }

    Loop<T>* loop_;
    std::list<LoopNode<T> > children_;
};

// Checks if a point is inside a polygon given by its vertices in correct order.
// The algorithm uses the ray crossing method described in "Computational Geometry in C"
template<typename T>
bool pointInPoly(const ClipVertex<T>& v, const Loop<T>& l) {
    bool c = false;
    tgt::vec2 vPos = v.pos2d_;
    
    tgt::vec2 pFirst = l.vertices_.back().pos2d_;
    tgt::vec2 pSecond;
    for (typename std::list<ClipVertex<T> >::const_iterator iter = l.vertices_.begin(); iter != l.vertices_.end(); ++iter) {
        pSecond = iter->pos2d_;

        if (((pSecond.y > vPos.y) != (pFirst.y > vPos.y)) && (vPos.x < (pFirst.x - pSecond.x) * (vPos.y - pSecond.y) / (pFirst.y - pSecond.y) + pSecond.x))
            c = !c;

        pFirst = pSecond;
    }

    return c;
}

// Tests if loop a lies in loop b by checking all vertices using the crossings test. 
/*template<typename T>
bool allVerticesInLoop(const Loop<T>& a, const Loop<T>& b) {

    for (typename std::list<ClipVertex<T> >::const_iterator iter = a.vertices_.begin(); iter != a.vertices_.end(); ++iter) {
        if (!pointInPoly(*iter, b))
            return false;
    }

    return true;
}*/

//checks if bounding box A is enclosed by bounding box B
bool boundingBoxIsEnclosed(tgt::vec2 bbMinA, tgt::vec2 bbMaxA, tgt::vec2 bbMinB, tgt::vec2 bbMaxB);

// tests if loop a is enclosed by loop b by testing a vertex of a against b and testing for edge intersections
template<typename T>
bool loopIsEnclosed(const Loop<T>& a, const Loop<T>& b) {

    // use bounding box test first
    if (!boundingBoxIsEnclosed(a.bbMin_, a.bbMax_, b.bbMin_, b.bbMax_))
        return false;

    // 1. test a single vertex of a against b
    if (!pointInPoly(a.vertices_.front(), b))
        return false;

    // 2. peform edge intersection tests (TODO: this has quadratic complexity, could be realized much more efficient)
    tgt::vec2 posABegin = a.vertices_.back().pos2d_;
    tgt::vec2 posAEnd;
    tgt::vec2 posBBegin = b.vertices_.back().pos2d_;
    tgt::vec2 posBEnd;
    for (typename std::list<ClipVertex<T> >::const_iterator aEdgeEnd = a.vertices_.begin(); aEdgeEnd != a.vertices_.end(); ++aEdgeEnd) {
        posAEnd = aEdgeEnd->pos2d_;
        for (typename std::list<ClipVertex<T> >::const_iterator bEdgeEnd = b.vertices_.begin(); bEdgeEnd != b.vertices_.end(); ++bEdgeEnd) {
            posBEnd = bEdgeEnd->pos2d_;
            if (linesIntersect(posABegin, posAEnd, posBBegin, posBEnd))
                return false;
            posBBegin = posBEnd;
        }    
        posABegin = posAEnd;
    }

    return true;
}


/**
 * Abstract base class for classes that construct polygons from a set of loops
 */
template<typename T>
class PolygonCreator {
   
public:

    //PolygonCreator();

    /// constructs polygons from the set of loops and appends them to the list
    virtual void createPolygons(std::list<Loop<T> >& loops, std::list<Loop<T>* >& polygons) = 0;

};


// for sorting loops according to the x coordinate of their first (transformed) vertex in the vertex list
template<typename T>
bool loopXSort(const Loop<T>* l1, const Loop<T>* l2) {
    return (l1->vertices_.front().pos2d_.x < l2->vertices_.front().pos2d_.x);
}


/**
 * Implements a hierarchical odering of the loops and tries to construct simple polygons from them. 
 * Outer boundaries are assumed to be CCW, holes are CW. Contour bridges are inserted to connect the loops of a polygon, creating a (degenerated) simple polygon.
 * The polygon creation might fail if loops are intersecting. 
 */
template<typename T>
class SimplePolygonCreator : public PolygonCreator<T> {
   
public:

    /// If the parameter is set to true, only the top level CCW loops of the hierarchy are written out
    SimplePolygonCreator(bool onlyUseTopLevel = false)
            : PolygonCreator<T>()
            , onlyUseTopLevel_(onlyUseTopLevel)
    {  }



    /// constructs polygons from the set of loops and appends them to the list
    virtual void createPolygons(std::list<Loop<T> >& loops, std::list<Loop<T>* >& polygons) {

        // build hierarchy
        std::list<LoopNode<T> > topLevel;
        for (typename std::list<Loop<T> >::iterator i = loops.begin(); i != loops.end(); ++i) 
            placeInHierarchy(topLevel, &(*i));

        if (onlyUseTopLevel_) {
            //just put out the top level CCW loops
            for (typename std::list<LoopNode<T> >::iterator i = topLevel.begin(); i != topLevel.end(); ++i)
                if (i->loop_->ccw_)
                    polygons.push_back(i->loop_);
        
        }
        else {
            //traverse hierarchy and create simple polygons
            buildSimplePolygons(polygons, topLevel);

        }
    }

protected:

    /// places a loop in the hierarchy starting at the given level
    void placeInHierarchy(std::list<LoopNode<T> >& level, Loop<T>* loop) { 

        //check if loop is inside the current reference loop
        for (typename std::list<LoopNode<T> >::iterator iter = level.begin(); iter != level.end(); ++iter) {
            //top level mode -> heuristic -> check all points
            if (onlyUseTopLevel_) {
                if (loopIsEnclosed<T>(*loop, *(iter->loop_))) {
                    placeInHierarchy(iter->children_, loop);
                    return;
                }
            }
            else {      //only check one point
                if (boundingBoxIsEnclosed(loop->bbMin_, loop->bbMax_, iter->loop_->bbMin_, iter->loop_->bbMax_) && pointInPoly<T>(*(loop->vertices_.begin()), *(iter->loop_))) {
                    placeInHierarchy(iter->children_, loop);
                    return;
                }
            }
        }

        // build new node for the loop
        LoopNode<T> node(loop);

        //check if loop encloses one of the current loops
        typename std::list<LoopNode<T> >::iterator iter = level.begin();
        while (iter != level.end()) {
            if (onlyUseTopLevel_) {
                if (loopIsEnclosed<T>(*(iter->loop_), *loop)) {
                    node.children_.push_back(*iter);
                    level.remove(*iter);
                    iter = level.begin();
                }
                else
                    ++iter;
            }
            else {
                if (boundingBoxIsEnclosed(loop->bbMin_, loop->bbMax_, iter->loop_->bbMin_, iter->loop_->bbMax_) && pointInPoly<T>(*(iter->loop_->vertices_.begin()), *loop)) {
                    node.children_.push_back(*iter);
                    level.remove(*iter);
                    iter = level.begin();
                }
                else
                    ++iter;
            }
        }
    
        level.push_back(node);
    }

    void buildSimplePolygons(std::list<Loop<T>*>& polygons, std::list<LoopNode<T> >& level) {
        for (typename std::list<LoopNode<T> >::iterator iter = level.begin(); iter != level.end(); ++iter) {
            if (!(iter->loop_->ccw_)) {
                //current loop is a hole and should be handled above / not triangulated -> traverse child nodes
                buildSimplePolygons(polygons, iter->children_);
            }
            else {
                //current loop is ccw -> find holes
                std::list<Loop<T>*> cwLoops;
                for (typename std::list<LoopNode<T> >::iterator childIter = iter->children_.begin(); childIter != iter->children_.end(); ++childIter) 
                    if (!(childIter->loop_->ccw_))
                        cwLoops.push_back(childIter->loop_);

                polygons.push_back(makePolygon(iter->loop_, cwLoops));
            
                //traverse children -> for cw loops, only children are traversed, for other loops triangulation is executed
                buildSimplePolygons(polygons, iter->children_);
            }
        }
    }

    Loop<T>* makePolygon(Loop<T>* outerLoop, std::list<Loop<T>*>& holes) {
        if (holes.empty())
            return outerLoop;

        //1. find vertex with minimal x coordinate for each hole and choose it as a starting point for the corresponding vertex list 
        for (typename std::list<Loop<T>*>::iterator holeIter = holes.begin(); holeIter != holes.end(); ++holeIter) {
            typename std::list<ClipVertex<T> >::iterator minIter = (*holeIter)->vertices_.begin();
            tgt::vec2 minPos = minIter->pos2d_;
            for (typename std::list<ClipVertex<T> >::iterator vertexIter = (*holeIter)->vertices_.begin(); vertexIter != (*holeIter)->vertices_.end(); ++vertexIter) {
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
                std::list<ClipVertex<T> > tmpList;
                tmpList.splice(tmpList.begin(), (*holeIter)->vertices_, (*holeIter)->vertices_.begin(), minIter);
                (*holeIter)->vertices_.splice((*holeIter)->vertices_.end(), tmpList);
            }
        } 

        //2. sort holes according to minimal x coordinate
        holes.sort(loopXSort<T>);

        //3. for each hole: find contour bridge and insert the vertices into the outer loop
        //we know that we have to take the first vertex of the hole and only vertices to the left of it are relevant
        for (typename std::list<Loop<T> *>::iterator holeIter = holes.begin(); holeIter != holes.end(); ++holeIter) {
            //find possible starting point for contour bridge
            tgt::vec2 endPos = (*holeIter)->vertices_.front().pos2d_;
            for (typename std::list<ClipVertex<T> >::iterator candidateIter = outerLoop->vertices_.begin(); candidateIter != outerLoop->vertices_.end(); ++candidateIter) {
                tgt::vec2 candidatePos = candidateIter->pos2d_;
                if (candidatePos.x < endPos.x) {

                    //possible starting point -> check new edge against all edges of the outer loop
                    bool intersect = false;

                    //iterate over edges and perform intersection test
                    tgt::vec2 e1 = outerLoop->vertices_.back().pos2d_;
                    tgt::vec2 e2;
                    bool lastWasCandidate = false;
                    for (typename std::list<ClipVertex<T> >::iterator edgeIter = outerLoop->vertices_.begin(); edgeIter != outerLoop->vertices_.end(); ++edgeIter) {
                        //TODO: now we also check for intersections with the edges that have the candidate point as starting / ending point
                        //since this always leads to an intersection, the intersection check has been modified to only return the case
                        //of a "real" intersection (touching is ok), but instead the loop here should be changed to not check against the two problematic edges
                        e2 = edgeIter->pos2d_;
    
                        if (linesProperIntersect(e1, e2, candidatePos, endPos)) {
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

    bool onlyUseTopLevel_;

};



}

#endif
