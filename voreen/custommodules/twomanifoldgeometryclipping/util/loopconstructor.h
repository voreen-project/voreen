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

#ifndef VRN_GPUGEOMETRYCLIPPING_LOOPCONSTRUCTOR
#define VRN_GPUGEOMETRYCLIPPING_LOOPCONSTRUCTOR

#include "twomanifoldclippingutil.h"

#include "ext/tgt/logmanager.h"

namespace voreen {

 
// sorting function for sorting edges lexicographically with respect to the starting vertex positions using std::sort
template<typename T>
bool lexicographicalEdgeSort(const ClipEdge<T>& e1, const ClipEdge<T>& e2) {
    if (e1.start_.vertex_.pos_.x != e2.start_.vertex_.pos_.x)
        return (e1.start_.vertex_.pos_.x < e2.start_.vertex_.pos_.x);
    else if (e1.start_.vertex_.pos_.y != e2.start_.vertex_.pos_.y)
        return (e1.start_.vertex_.pos_.y < e2.start_.vertex_.pos_.y);
    else 
        return (e1.start_.vertex_.pos_.z < e2.start_.vertex_.pos_.z);
}  


// compare function for std::find_if to find the successor to an edge
template<typename T>
struct successor_compare : public std::unary_function<ClipEdge<T>, bool> {

    explicit successor_compare(const ClipEdge<T>& edge, float epsilon = 0.f) 
        : baseEdge_(edge) 
        , epsilon_(epsilon)  
    {}
        
    bool operator() (const ClipEdge<T>& arg) { 
        if (epsilon_ == 0.f)
            return (baseEdge_.end_.vertex_.pos_ == arg.start_.vertex_.pos_);
        else
            return (tgt::distance(baseEdge_.end_.vertex_.pos_, arg.start_.vertex_.pos_) <= epsilon_);
    }
          
    ClipEdge<T> baseEdge_;
    float epsilon_;
};


//computes the winding order of the (transformed) vertices (cw or ccw) and sets the coorresponding flag
//additionally computes the bounding box of the loop
template<typename T>
void computeLoopOrientationAndBoundingBox(Loop<T>& loop) {
    float area = 0;
    //begin with edge to the first vertex and then iterate through all edges
    tgt::vec2 begin = loop.vertices_.back().pos2d_;
    loop.bbMin_ = loop.bbMax_ = begin;
    tgt::vec2 end;
    for (typename std::list<ClipVertex<T> >::iterator i = loop.vertices_.begin(); i != loop.vertices_.end(); ++i) {
        end = i->pos2d_;
        loop.bbMin_ = tgt::min(loop.bbMin_, end);
        loop.bbMax_ = tgt::max(loop.bbMax_, end);
        area += begin.x * end.y - begin.y * end.x;
        begin = end;
    }

    loop.orientation_ = area;
    loop.ccw_ = (area > 0);
}

/**
 * Abstract base class for classes that construct loops from a set of edges
 */
template<typename T>
class LoopConstructor {
   
public:

    LoopConstructor(float epsilon = 0.f) 
        : epsilon_(epsilon)
    { }

    /// constructs loops from the set of edges and appends them to the list
    virtual void constructLoops(std::vector<ClipEdge<T> >& edges, std::list<Loop<T> >& loops) = 0;

protected:

    float epsilon_;     ///< epsilon for finding successor to an edge

};

// define for heuristic method
#define CLOSE_INCOMPLETE_LOOPS

template<typename T>
class SimpleLoopConstructor : LoopConstructor<T> {

public:
    SimpleLoopConstructor(float epsilon = 0.f)
        : LoopConstructor<T>(epsilon)
    { }

    virtual void constructLoops(std::vector<ClipEdge<T> >& edges, std::list<Loop<T> >& loops) {

        // sort edges lexicographically
        std::sort(edges.begin(), edges.end(), lexicographicalEdgeSort<T>);

        // find loops using binary search
#ifdef CLOSE_INCOMPLETE_LOOPS
        std::list<Loop<T> > incompleteLoops;
#endif

        while (!edges.empty()) {
            Loop<T> l;
            l.complete_ = true;
            typename std::vector<ClipEdge<T> >::iterator i = edges.begin();
            l.vertices_.push_back(i->start_);
            ClipEdge<T> currentEdge = *i;
            edges.erase(i);

            while (!(currentEdge.end_.vertex_.pos_ == l.vertices_.front().vertex_.pos_)) {
                typename std::vector<ClipEdge<T> >::iterator next = std::find_if(edges.begin(), edges.end(), successor_compare<T>(currentEdge, this->epsilon_)); 
                if (next == edges.end()) {
#ifdef CLOSE_INCOMPLETE_LOOPS

                    //no successor found -> check incomplete loops starting point
                    bool foundSuccessorLoop = false;
                    for (typename std::list<Loop<T> >::iterator incompleteIter = incompleteLoops.begin(); incompleteIter != incompleteLoops.end(); ++incompleteIter) {
                        if (tgt::distance(incompleteIter->vertices_.front().vertex_.pos_, currentEdge.end_.vertex_.pos_) <= this->epsilon_) {
                            foundSuccessorLoop = true;
                            //put l at the beginning of the found incomplete loop
                            incompleteIter->vertices_.splice(incompleteIter->vertices_.begin(), l.vertices_);
                            break;
                        }
                    } 
                    if (!foundSuccessorLoop) {
                        //no loop / edge with the end vertex of the current edge exists -> push back the end vertex as it is the end of the vertex list
                        l.vertices_.push_back(currentEdge.end_);
                        l.complete_ = false;
                        incompleteLoops.push_back(l);
                    }
                    break; 
#else
                    LWARNINGC("voreen.gpugeometryclipping.loopconstructor", "Found no successor to edge... loop is not complete!");
                    l.complete_ = false;
                    break;
#endif                
                }
                else {
                    currentEdge = *next;
                    edges.erase(next);
                    l.vertices_.push_back(currentEdge.start_);
                }
            } 

#ifdef CLOSE_INCOMPLETE_LOOPS
            if (!l.vertices_.empty() && l.complete_)
                loops.push_back(l);
#else      
            if (l.vertices_.size() > 2)
                loops.push_back(l);
            else
                LWARNINGC("voreen.gpugeometryclipping.loopconstructor", "Loop with 2 vertices");
#endif
        }

#ifdef CLOSE_INCOMPLETE_LOOPS
        //remove loops with only one edge, close the other ones and insert them into the complete loops
        typename std::list<Loop<T> >::iterator incompleteLoopIter = incompleteLoops.begin();
        size_t numIncompleteLoops = 0;
        while (!incompleteLoops.empty()) {
            //remove loops with only one edge 
            if (incompleteLoopIter->vertices_.size() == 2)
                incompleteLoops.erase(incompleteLoopIter);
            else {
                //other loops: set "complete" status (acts as if an edge exists from the last to the first vertex, ie. closing the loop
                incompleteLoopIter->complete_ = true;
                loops.splice(loops.begin(), incompleteLoops, incompleteLoopIter);
                numIncompleteLoops++;
            }   
            
            incompleteLoopIter = incompleteLoops.begin(); 
        }
       
        if (numIncompleteLoops > 0)
            LWARNINGC("voreen.gpugeometryclipping.loopconstructor", "Either due to non-manifold mesh or numerical issues " << numIncompleteLoops << " incomplete loops are closed.");
#endif

        //compute orientation and bounding box for every loop
        //TODO: this should be done during loop construction 
        for (typename std::list<Loop<T> >::iterator iter = loops.begin(); iter != loops.end(); ++iter) { 
            computeLoopOrientationAndBoundingBox<T>(*iter);
        }
    }

};    



} //namespace

#endif
