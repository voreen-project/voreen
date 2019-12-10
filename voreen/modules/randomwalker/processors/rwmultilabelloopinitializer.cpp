/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2019 University of Muenster, Germany,                        *
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

#include "rwmultilabelloopinitializer.h"

#include "voreen/core/datastructures/volume/volumeram.h"
#include "voreen/core/datastructures/volume/volume.h"
#include "voreen/core/datastructures/volume/volumeatomic.h"
#include "voreen/core/datastructures/geometry/pointlistgeometry.h"
#include "voreen/core/datastructures/geometry/pointsegmentlistgeometry.h"
#include "voreen/core/voreenapplication.h"

#include "tgt/vector.h"

namespace voreen {

const std::string RWMultiLabelLoopInitializer::loggerCat_("voreen.RandomWalker.RWMultiLabelLoopInitializer");

RWMultiLabelLoopInitializer::RWMultiLabelLoopInitializer()
    : VolumeProcessor(),
      seedsInport_(Port::INPORT, "seeds.in", "seeds.in", true),
      foregroundSeedsOutport_(Port::OUTPORT, "foregroundSeeds.out"),
      backgroundSeedsOutport_(Port::OUTPORT, "backgroundSeeds.out"),
      loopInport_(Port::INPORT, "loop.in"),
      startButton_("startButton", "Start"),
      rwConnector_("rwConnector", "RW Connector"),
      running_(false)
{
    addPort(seedsInport_);
    addPort(foregroundSeedsOutport_);
    addPort(backgroundSeedsOutport_);

    loopInport_.setLoopPort(true);
    addPort(loopInport_);

    addProperty(startButton_);
    addProperty(rwConnector_);

    startButton_.onClick(MemberFunctionCallback<RWMultiLabelLoopInitializer>(this, &RWMultiLabelLoopInitializer::startButtonClicked));
}

RWMultiLabelLoopInitializer::~RWMultiLabelLoopInitializer() {
}

Processor* RWMultiLabelLoopInitializer::create() const {
    return new RWMultiLabelLoopInitializer();
}

void RWMultiLabelLoopInitializer::initialize() {
    VolumeProcessor::initialize();
}

void RWMultiLabelLoopInitializer::deinitialize() {
    VolumeProcessor::deinitialize();
}

void RWMultiLabelLoopInitializer::process() {

    std::vector<const Geometry*> inputGeom = seedsInport_.getAllData();

    std::vector<PointListGeometry<tgt::vec3>*> labelSeeds;

    for (size_t i=0; i<inputGeom.size(); i++) {
        const PointSegmentListGeometry<tgt::vec3>* seedList = dynamic_cast<const PointSegmentListGeometry<tgt::vec3>* >(
            inputGeom.at(i));
        if (!seedList)
            LWARNING("Invalid geometry. PointSegmentListGeometry<vec3> expected.");
        else {
            for (int j=0; j<seedList->getNumSegments(); j++) {
                if (seedList->getSegment(j).size() > 0) {
                    std::vector<tgt::vec3> segmentPoints = seedList->getSegment(j);
                    PointListGeometryVec3* segment = new PointListGeometryVec3();
                    segment->setData(segmentPoints);
                    labelSeeds.push_back(segment);
                }
            }
        }
    }

    if (labelSeeds.size() < 2) {
        LWARNING("at least two labels required");
        return;
    }
    //LINFO("Num labels: " << labelSeeds.size());
    if (labelSeeds.size() != (size_t)loopInport_.getNumLoopIterations())
        loopInport_.setNumLoopIterations(static_cast<int>(labelSeeds.size()));

    if (!running_)
        return;

    //currentLabel_.setMaxValue(labelSeeds.size()-1);
    size_t currentLabel = loopInport_.getLoopIteration();

    PointSegmentListGeometryVec3* foregroundSeeds = new PointSegmentListGeometryVec3();
    PointSegmentListGeometryVec3* backgroundSeeds = new PointSegmentListGeometryVec3();
    for (size_t i=0; i<labelSeeds.size(); i++) {
        if (i == currentLabel)
            foregroundSeeds->addSegment(labelSeeds.at(i)->getData());
        else
            backgroundSeeds->addSegment(labelSeeds.at(i)->getData());
    }

    foregroundSeedsOutport_.setData(foregroundSeeds);
    backgroundSeedsOutport_.setData(backgroundSeeds);

    rwConnector_.clicked();
    if (loopInport_.getLoopIteration() == loopInport_.getNumLoopIterations()-1)
        running_ = false;

}

void RWMultiLabelLoopInitializer::startButtonClicked() {
//    loopInport_.setLoopIteration(0);
    running_ = true;
    invalidate();
}

}   // namespace
