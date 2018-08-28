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

#ifndef VRN_PLOTZOOMSTATE_H
#define VRN_PLOTZOOMSTATE_H

#include "plotbase.h"
#include "interval.h"

#include "voreen/core/io/serialization/serializable.h"

namespace voreen {

/**
 * This class encapsulates three intervals (x,y and z axis). It is used
 * as a zoom state.
 */
class PlotZoomState : public Serializable {
public:
    PlotZoomState()
        : xZoom_(Interval<plot_t>(0, 0))
        , yZoom_(Interval<plot_t>(0, 0))
        , zZoom_(Interval<plot_t>(0, 0))
    {};

    PlotZoomState(Interval<plot_t> xZoom, Interval<plot_t> yZoom, Interval<plot_t> zZoom = Interval<plot_t>(0, 0))
        : xZoom_(xZoom), yZoom_(yZoom), zZoom_(zZoom)
    {};

    bool operator==(const PlotZoomState& rhs) const {
        if (xZoom_ != rhs.xZoom_)
            return false;
        if (yZoom_ != rhs.yZoom_)
            return false;
        if (zZoom_ != rhs.zZoom_)
            return false;
        return true;
    }

    void serialize(Serializer& s) const {
        s.serialize("xZoom", xZoom_);
        s.serialize("yZoom", yZoom_);
        s.serialize("zZoom", zZoom_);
    };

    void deserialize(Deserializer& s) {
        s.deserialize("xZoom", xZoom_);
        s.deserialize("yZoom", yZoom_);
        s.deserialize("zZoom", zZoom_);
    };

    Interval<plot_t> xZoom_;    ///< zoom interval for x axis
    Interval<plot_t> yZoom_;    ///< zoom interval for y axis
    Interval<plot_t> zZoom_;    ///< zoom interval for z axis
};



}

#endif // VRN_PLOTZOOMSTATE_H

