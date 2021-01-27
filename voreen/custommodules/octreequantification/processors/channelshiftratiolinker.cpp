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

#include "channelshiftratiolinker.h"

namespace voreen {

const std::string ChannelShiftRatioLinker::loggerCat_("voreen.OctreeRegionGrowing.ChannelShiftRatioLinker");

ChannelShiftRatioLinker::ChannelShiftRatioLinker() 
        : Processor()
        , inportA_(Port::INPORT, "volumeinportA", "Volume Inport A")
        , inportB_(Port::INPORT, "volumeinportB", "Volume Inport B")
        , channelSelect_("channelselct", "Select Channel", 1, 1, 4)
        , applyChannelShift_("applyChannelShift", "Apply Channel Shift", false, Processor::INVALID_PROGRAM, Property::LOD_ADVANCED)
        , channelShiftA1_("channelShiftA0", "Channel Shift A1 (Input Channel 1)", tgt::vec3(0.f), tgt::vec3(-50.f), tgt::vec3(50.f), Processor::INVALID_RESULT, NumericProperty<tgt::vec3>::STATIC, Property::LOD_ADVANCED)
        , channelShiftA2_("channelShiftA1", "Channel Shift A2 (Input Channel 2)", tgt::vec3(0.f), tgt::vec3(-50.f), tgt::vec3(50.f), Processor::INVALID_RESULT, NumericProperty<tgt::vec3>::STATIC, Property::LOD_ADVANCED)
        , channelShiftA3_("channelShiftA2", "Channel Shift A3 (Input Channel 3)", tgt::vec3(0.f), tgt::vec3(-50.f), tgt::vec3(50.f), Processor::INVALID_RESULT, NumericProperty<tgt::vec3>::STATIC, Property::LOD_ADVANCED)
        , channelShiftA4_("channelShiftA3", "Channel Shift A4 (Input Channel 4)", tgt::vec3(0.f), tgt::vec3(-50.f), tgt::vec3(50.f), Processor::INVALID_RESULT, NumericProperty<tgt::vec3>::STATIC, Property::LOD_ADVANCED)
        , channelShiftB_("channelShiftB", "Channel Shift B (Output)", tgt::vec3(0.f), tgt::vec3(-50.f), tgt::vec3(50.f), Processor::INVALID_RESULT, NumericProperty<tgt::vec3>::STATIC, Property::LOD_ADVANCED)
{

    addPort(inportA_);
    addPort(inportB_);
    inportA_.onChange(MemberFunctionCallback<ChannelShiftRatioLinker>(this, &ChannelShiftRatioLinker::inputVolumesChanged));
    inportB_.onChange(MemberFunctionCallback<ChannelShiftRatioLinker>(this, &ChannelShiftRatioLinker::inputVolumesChanged));

    addProperty(channelSelect_);
    addProperty(applyChannelShift_);

    channelSelect_.onChange(MemberFunctionCallback<ChannelShiftRatioLinker>(this, &ChannelShiftRatioLinker::channelShiftAChanged));
    channelShiftA1_.onChange(MemberFunctionCallback<ChannelShiftRatioLinker>(this, &ChannelShiftRatioLinker::channelShiftAChanged));
    channelShiftA2_.onChange(MemberFunctionCallback<ChannelShiftRatioLinker>(this, &ChannelShiftRatioLinker::channelShiftAChanged));
    channelShiftA3_.onChange(MemberFunctionCallback<ChannelShiftRatioLinker>(this, &ChannelShiftRatioLinker::channelShiftAChanged));
    channelShiftA4_.onChange(MemberFunctionCallback<ChannelShiftRatioLinker>(this, &ChannelShiftRatioLinker::channelShiftAChanged));

    addProperty(channelShiftA1_);
    addProperty(channelShiftA2_);
    addProperty(channelShiftA3_);
    addProperty(channelShiftA4_);

    addProperty(channelShiftB_);
}

Processor* ChannelShiftRatioLinker::create() const {
    return new ChannelShiftRatioLinker();
}

void ChannelShiftRatioLinker::process() {
    // do nothing :)
}

void ChannelShiftRatioLinker::inputVolumesChanged() {
    if (!inportA_.getData() || !inportB_.getData())
        return;

    // -> compute
    channelShiftAChanged();
}

void ChannelShiftRatioLinker::channelShiftAChanged() {
    if (!inportA_.getData() || !inportB_.getData())
        return;

    // get current channel shift
    tgt::vec3 shift;
    switch (channelSelect_.get()) {
        case 1: { shift = channelShiftA1_.get(); break; }
        case 2: { shift = channelShiftA2_.get(); break; }
        case 3: { shift = channelShiftA3_.get(); break; }
        case 4: { shift = channelShiftA4_.get(); break; }
        default: shift = tgt::vec3(0.f);
    }

    // compute ratio and set value
    tgt::vec3 ratio = tgt::vec3(inportB_.getData()->getDimensions() - tgt::svec3::one) / tgt::vec3(inportA_.getData()->getDimensions() - tgt::svec3::one);
    channelShiftB_.set(shift * ratio);

}

}   // namespace
