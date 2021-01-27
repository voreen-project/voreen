/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2021 University of Muenster, Germany,                        *
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

#ifndef VRN_GPUCLIPPINGSTRUCTS_H
#define VRN_GPUCLIPPINGSTRUCTS_H

#include <list>
#include <algorithm>
#include "voreen/core/datastructures/geometry/vertex.h"

namespace voreen {

//------------------------------
//data types
//------------------------------

/**
 * An extension of the Voreen vertex layouts.
 * Stores additional information needed for loop construction and triangulation.
 * T should be VertexBase or one of its derived classes.
 */
template<typename T>
struct ClipVertex {

    T vertex_;                  ///< the actual vertex 
    tgt::vec2 pos2d_;           ///< the transformed 2d position of the vertex
    
    bool earTipStatus_;         ///< needed for triangulation
    float minCos_;              ///< cos of the maximum angle in the ear of which the vertex is the ear tip

    /**
     * Compares two extended vertices by their original vertices.
     * Epsilon is set to zero.
     */
    bool operator==(const ClipVertex& a) const {
        //return (position_ == a.position_ && normal_ == a.normal_ && texCoords_ == a.texCoords_);
        return (vertex_.equals(a.vertex_, 0.0));  
    } 
};

/**
 * An edge formed by two clip vertices.
 */ 
template<typename T>
struct ClipEdge {
    ClipVertex<T> start_;
    ClipVertex<T> end_;
};

/**
 * A vertex loop.
 */ 
template<typename T>
struct Loop {
    std::list<ClipVertex<T> > vertices_;
    bool ccw_;
    bool complete_;
    float orientation_; ///< used as a temporary storage for computing the loop orientation during construction

    tgt::vec2 bbMin_;   ///< lower left of the 2d bounding box
    tgt::vec2 bbMax_;   ///< upper right of the 2d bounding box
};

/// orientation test for three points
float orientationTest(const tgt::vec2& pa, const tgt::vec2& pb, const tgt::vec2& pc);

/// checks if pc (strictly) lies to the left of the directed line through pa and pb
bool leftTurn(const tgt::vec2& pa, const tgt::vec2& pb, const tgt::vec2& pc);

/// checks if pc lies to the left of the directed line through pa and pb OR if the three points are collinear
bool leftTurnCollinear(const tgt::vec2& pa, const tgt::vec2& pb, const tgt::vec2& pc);

/** 
 * Performs the orientation test and compares to epsilon using >=
 * To increase the tolerance for regarding as a left turn, use a negative epsilon.
 */
bool leftTurnEpsilon(const tgt::vec2& pa, const tgt::vec2& pb, const tgt::vec2& pc, float epsilon);

/**
 * checks if the line segments (p1,p2) and (p3,p4) intersect
 * this only checks for "real" intersections, ie. the line segments may be touching
 */
bool linesProperIntersect(const tgt::vec2& p0, const tgt::vec2& p1, const tgt::vec2& p2, const tgt::vec2& p3);

/**
 * checks if the line segments (p1,p2) and (p3,p4) intersect
 * (also includes improper intersections)
 */ 
bool linesIntersect(const tgt::vec2& p0, const tgt::vec2& p1, const tgt::vec2& p2, const tgt::vec2& p3);

} //namespace


#endif
