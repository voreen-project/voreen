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

#include "boundingboxratiolinker.h"

namespace voreen {

const std::string BoundingBoxRatioLinker::loggerCat_("voreen.OctreeRegionGrowing.BoundingBoxRatioLinker");

BoundingBoxRatioLinker::BoundingBoxRatioLinker() 
        : Processor()
        , inportA_(Port::INPORT, "volumeinportA", "Volume Inport A")
        , inportB_(Port::INPORT, "volumeinportB", "Volume Inport B")
        , clipRegionA_("clipRegionA", "Clip Region A", tgt::IntBounds(tgt::ivec3(0), tgt::ivec3(10000)), tgt::ivec3(0), tgt::ivec3(10000))
        , clipRegionB_("clipRegionB", "Clip Region B", tgt::IntBounds(tgt::ivec3(0), tgt::ivec3(10000)), tgt::ivec3(0), tgt::ivec3(10000))
{

    addPort(inportA_);
    addPort(inportB_);
    inportA_.onChange(MemberFunctionCallback<BoundingBoxRatioLinker>(this, &BoundingBoxRatioLinker::inputVolumesChanged));
    inportB_.onChange(MemberFunctionCallback<BoundingBoxRatioLinker>(this, &BoundingBoxRatioLinker::inputVolumesChanged));

    clipRegionA_.onChange(MemberFunctionCallback<BoundingBoxRatioLinker>(this, &BoundingBoxRatioLinker::clipRegionAChanged));
    clipRegionB_.onChange(MemberFunctionCallback<BoundingBoxRatioLinker>(this, &BoundingBoxRatioLinker::clipRegionBChanged));
    addProperty(clipRegionA_);
    addProperty(clipRegionB_);
}

Processor* BoundingBoxRatioLinker::create() const {
    return new BoundingBoxRatioLinker();
}

void BoundingBoxRatioLinker::process() {
    // do nothing :)
}

void BoundingBoxRatioLinker::inputVolumesChanged() {
    if (!inportA_.getData() || !inportB_.getData())
        return;

    // set max values
    clipRegionA_.setMaxValue(tgt::ivec3(inportA_.getData()->getDimensions()) - tgt::ivec3::one);
    clipRegionB_.setMaxValue(tgt::ivec3(inportB_.getData()->getDimensions()) - tgt::ivec3::one);
}

void BoundingBoxRatioLinker::clipRegionAChanged() {
    if (!inportA_.getData() || !inportB_.getData())
        return;
    
    tgt::vec3 llfRatio = tgt::vec3(clipRegionA_.get().getLLF()) / tgt::vec3(inportA_.getData()->getDimensions() - tgt::svec3::one);
    tgt::vec3 urbRatio = tgt::vec3(clipRegionA_.get().getURB()) / tgt::vec3(inportA_.getData()->getDimensions() - tgt::svec3::one);

    tgt::ivec3 resultLlf = tgt::iround(llfRatio * tgt::vec3(inportB_.getData()->getDimensions() - tgt::svec3::one));
    tgt::ivec3 resultUrb = tgt::iround(urbRatio * tgt::vec3(inportB_.getData()->getDimensions() - tgt::svec3::one));

    clipRegionB_.set(tgt::IntBounds(resultLlf, resultUrb));
}

void BoundingBoxRatioLinker::clipRegionBChanged() {
    if (!inportA_.getData() || !inportB_.getData())
        return;
    
    tgt::vec3 llfRatio = tgt::vec3(clipRegionB_.get().getLLF()) / tgt::vec3(inportB_.getData()->getDimensions() - tgt::svec3::one);
    tgt::vec3 urbRatio = tgt::vec3(clipRegionB_.get().getURB()) / tgt::vec3(inportB_.getData()->getDimensions() - tgt::svec3::one);

    tgt::ivec3 resultLlf = tgt::iround(llfRatio * tgt::vec3(inportA_.getData()->getDimensions() - tgt::svec3::one));
    tgt::ivec3 resultUrb = tgt::iround(urbRatio * tgt::vec3(inportA_.getData()->getDimensions() - tgt::svec3::one));

    clipRegionA_.set(tgt::IntBounds(resultLlf, resultUrb));
}

}   // namespace
