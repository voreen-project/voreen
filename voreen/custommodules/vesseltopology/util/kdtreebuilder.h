/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2016 University of Muenster, Germany.                        *
 * Visualization and Computer Graphics Group <http://viscg.uni-muenster.de>        *
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
#pragma once

#include <vector>
#include "../ext/nanoflann/nanoflann.hpp"

template<typename E>
class KDTreeBuilder {
public:
    typedef typename E::VoxelType::ElemType coord_t;

    KDTreeBuilder()
        : points_()
    {
    }

    KDTreeBuilder(KDTreeBuilder&& other)
        : points_(std::move(other.points_))
    {
    }

    const std::vector<E>& points() const {
        return points_;
    }

    std::vector<E>& points() {
        return points_;
    }

    void push(E element) {
        points_.push_back(element);
    }

    // Must return the number of data points
    inline size_t kdtree_get_point_count() const {
        return points_.size();
    }

    // Returns the dim'th component of the idx'th point in the class:
    inline coord_t kdtree_get_pt(const size_t idx, int dim) const {
        return points_[idx].at(dim);
    }

    template <class BBOX>
    bool kdtree_get_bbox(BBOX&) const {
        return false;
    }
private:
    std::vector<E> points_;
};
