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

#include "multipassproxygeometry.h"
#include "voreen/core/datastructures/geometry/meshlistgeometry.h"
#include "voreen/core/datastructures/geometry/trianglemeshgeometry.h"
#include "voreen/core/datastructures/geometry/geometrysequence.h"
#include "tgt/glmath.h"
#include <vector>

namespace voreen {

MultiPassProxyGeometry::MultiPassProxyGeometry()
    : Processor()
    , camera_("camera", "Camera", tgt::Camera(), Processor::VALID)
    , test_("test", "test", 1, 0, 10/*, Processor::INVALID_PORTS*/)
    , loopInport_(Port::INPORT, "loop.inport", "Loop Inport")
    , inportProxyGeometry_(Port::INPORT, "proxygeometry_in", "Proxy Geometry"/*, Processor::INVALID_PORTS*/)
    , outportProxyGeometry_(Port::OUTPORT, "proxygeometry_out", "Proxy Geometry")
{
    loopInport_.setLoopPort(true);

    addPort(loopInport_);
    addPort(inportProxyGeometry_);
    addPort(outportProxyGeometry_);

    addProperty(test_);
    addProperty(camera_);
    camera_.onChange(MemberFunctionCallback<MultiPassProxyGeometry>(this, &MultiPassProxyGeometry::performDepthSort));
}

MultiPassProxyGeometry::~MultiPassProxyGeometry() {
}

Processor* MultiPassProxyGeometry::create() const {
    return new MultiPassProxyGeometry();
}

struct DepthSort {
    DepthSort(tgt::vec3 cameraPos) : cameraPos_(cameraPos) {}

    bool operator() (const Geometry* i, const Geometry* j) {
        float d1 = distance(cameraPos_, i->getBoundingBox().center());
        float d2 = distance(cameraPos_, j->getBoundingBox().center());
        return (d1>d2);
    }
protected:
    tgt::vec3 cameraPos_;
};

void MultiPassProxyGeometry::performDepthSort() {
    std::sort(sortedInput_.begin(), sortedInput_.end(), DepthSort(camera_.get().getPosition()));
}

void MultiPassProxyGeometry::beforeProcess() {
    Processor::beforeProcess();

    if((loopInport_.getLoopIteration() == 0) && inportProxyGeometry_.hasChanged()) {
        const Geometry* g = inportProxyGeometry_.getData();
        sortedInput_.clear();

        if(g) {
            const GeometrySequence* gs = dynamic_cast<const GeometrySequence*>(g);
            if(gs) {
                if(static_cast<size_t>(loopInport_.getNumLoopIterations()) != gs->getNumGeometries()) {
                    loopInport_.setNumLoopIterations(static_cast<int>(gs->getNumGeometries()));
                    invalidate(INVALID_PORTS);
                }

                for(size_t i=0; i<gs->getNumGeometries(); i++)
                    sortedInput_.push_back(gs->getGeometry(i));

                performDepthSort();
            }
            else if(loopInport_.getNumLoopIterations() != 1) {
                loopInport_.setNumLoopIterations(1);
                invalidate(INVALID_PORTS);
            }
        }
        else if(loopInport_.getNumLoopIterations() != 0) {
            loopInport_.setNumLoopIterations(0);
            invalidate(INVALID_PORTS);
        }

        inportProxyGeometry_.setValid();
    }
}

void MultiPassProxyGeometry::afterProcess() {
    Processor::afterProcess();
}

void MultiPassProxyGeometry::process() {
    int iteration = loopInport_.getLoopIteration();
    //LWARNING(iteration);

    const Geometry* g = inportProxyGeometry_.getData();
    if(g) {
        const GeometrySequence* gs = dynamic_cast<const GeometrySequence*>(g);
        if(gs) {
            //outportProxyGeometry_.setData(gs->getGeometry(iteration), false);
            if(iteration < static_cast<int>(sortedInput_.size()))
                outportProxyGeometry_.setData(sortedInput_.at(iteration), false);
        }
        else
            outportProxyGeometry_.setData(g, false);
    }
    else
        outportProxyGeometry_.setData(0);

    LGL_ERROR;
}

} // namespace
