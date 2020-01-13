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

#include "slicenumberratiolinker.h"

namespace voreen {

const std::string SliceNumberRatioLinker::loggerCat_("voreen.OctreeRegionGrowing.SliceNumberRatioLinker");

SliceNumberRatioLinker::SliceNumberRatioLinker() 
        : Processor()
        , inportA_(Port::INPORT, "volumeinportA", "Volume Inport A")
        , inportB_(Port::INPORT, "volumeinportB", "Volume Inport B")
        , sliceNumberXSideA_("sliceNumberXSideA", "Slice Number X Side A", 0, 0, 10000)
        , sliceNumberYSideA_("sliceNumberYSideA", "Slice Number Y Side A", 0, 0, 10000)
        , sliceNumberZSideA_("sliceNumberZSideA", "Slice Number Z Side A", 0, 0, 10000)         
        , sliceNumberXSideB_("sliceNumberXSideB", "Slice Number X Side B", 0, 0, 10000)
        , sliceNumberYSideB_("sliceNumberYSideB", "Slice Number Y Side B", 0, 0, 10000)
        , sliceNumberZSideB_("sliceNumberZSideB", "Slice Number Z Side B", 0, 0, 10000)         
{

    addPort(inportA_);
    addPort(inportB_);
    inportA_.onChange(MemberFunctionCallback<SliceNumberRatioLinker>(this, &SliceNumberRatioLinker::inputVolumesChanged));
    inportB_.onChange(MemberFunctionCallback<SliceNumberRatioLinker>(this, &SliceNumberRatioLinker::inputVolumesChanged));

    sliceNumberXSideA_.onChange(MemberFunctionCallback<SliceNumberRatioLinker>(this, &SliceNumberRatioLinker::sliceNumbersAChanged));
    sliceNumberYSideA_.onChange(MemberFunctionCallback<SliceNumberRatioLinker>(this, &SliceNumberRatioLinker::sliceNumbersAChanged));
    sliceNumberZSideA_.onChange(MemberFunctionCallback<SliceNumberRatioLinker>(this, &SliceNumberRatioLinker::sliceNumbersAChanged));
    addProperty(sliceNumberXSideA_);
    addProperty(sliceNumberYSideA_);
    addProperty(sliceNumberZSideA_);

    sliceNumberXSideB_.onChange(MemberFunctionCallback<SliceNumberRatioLinker>(this, &SliceNumberRatioLinker::sliceNumbersBChanged));
    sliceNumberYSideB_.onChange(MemberFunctionCallback<SliceNumberRatioLinker>(this, &SliceNumberRatioLinker::sliceNumbersBChanged));
    sliceNumberZSideB_.onChange(MemberFunctionCallback<SliceNumberRatioLinker>(this, &SliceNumberRatioLinker::sliceNumbersBChanged));
    addProperty(sliceNumberXSideB_);
    addProperty(sliceNumberYSideB_);
    addProperty(sliceNumberZSideB_);
}

Processor* SliceNumberRatioLinker::create() const {
    return new SliceNumberRatioLinker();
}

void SliceNumberRatioLinker::process() {
    // do nothing :)
}

void SliceNumberRatioLinker::inputVolumesChanged() {
    if (!inportA_.getData() || !inportB_.getData())
        return;

    // set max values
    tgt::ivec3 dimA = tgt::ivec3(inportA_.getData()->getDimensions()) - tgt::ivec3::one;
    tgt::ivec3 dimB = tgt::ivec3(inportB_.getData()->getDimensions()) - tgt::ivec3::one;

    sliceNumberXSideA_.setMaxValue(dimA.x);
    sliceNumberYSideA_.setMaxValue(dimA.y);
    sliceNumberZSideA_.setMaxValue(dimA.z);

    sliceNumberXSideB_.setMaxValue(dimB.x);
    sliceNumberYSideB_.setMaxValue(dimB.y);
    sliceNumberZSideB_.setMaxValue(dimB.z);
}

void SliceNumberRatioLinker::sliceNumbersAChanged() {
    if (!inportA_.getData() || !inportB_.getData())
        return;

    tgt::ivec3 sliceNumbersA(sliceNumberXSideA_.get(), sliceNumberYSideA_.get(), sliceNumberZSideA_.get());
    tgt::vec3 sliceNumbersRatio = tgt::vec3(sliceNumbersA) / tgt::vec3(inportA_.getData()->getDimensions() - tgt::svec3::one);
    tgt::ivec3 result = tgt::iround(sliceNumbersRatio * tgt::vec3(inportB_.getData()->getDimensions() - tgt::svec3::one));

    sliceNumberXSideB_.set(std::min(result.x, sliceNumberXSideB_.getMaxValue()));
    sliceNumberYSideB_.set(std::min(result.y, sliceNumberYSideB_.getMaxValue()));    
    sliceNumberZSideB_.set(std::min(result.z, sliceNumberZSideB_.getMaxValue()));    
}

void SliceNumberRatioLinker::sliceNumbersBChanged() {
    if (!inportA_.getData() || !inportB_.getData())
        return;

    tgt::ivec3 sliceNumbersB(sliceNumberXSideB_.get(), sliceNumberYSideB_.get(), sliceNumberZSideB_.get());
    tgt::vec3 sliceNumbersRatio = tgt::vec3(sliceNumbersB) / tgt::vec3(inportB_.getData()->getDimensions() - tgt::svec3::one);
    tgt::ivec3 result = tgt::iround(sliceNumbersRatio * tgt::vec3(inportA_.getData()->getDimensions() - tgt::svec3::one));

    sliceNumberXSideA_.set(std::min(result.x, sliceNumberXSideA_.getMaxValue()));
    sliceNumberYSideA_.set(std::min(result.y, sliceNumberYSideA_.getMaxValue()));    
    sliceNumberZSideA_.set(std::min(result.z, sliceNumberZSideA_.getMaxValue()));    
}

}   // namespace
