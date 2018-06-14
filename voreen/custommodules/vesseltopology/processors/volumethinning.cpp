/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2015 University of Muenster, Germany.                        *
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

#include "volumethinning.h"

#include "voreen/core/datastructures/volume/volumeram.h"
#include "voreen/core/datastructures/volume/volumefactory.h"
#include "voreen/core/datastructures/volume/volumeminmax.h"
#include "voreen/core/io/progressreporter.h"

#include <chrono>

namespace voreen {

const std::string VolumeThinning::loggerCat_("voreen.vesseltopology.volumethinning");

VolumeThinning::VolumeThinning()
    : VolumeProcessor()
    , inport_(Port::INPORT, "volumethinning.inport", "Binary Volume Input")
    , sampleMaskInport_(Port::INPORT, "volumethinning.sampleMask.inport", "Sample Mask (optional)")
    , outport_(Port::OUTPORT, "volumethinning.outport", "Volume Output", false, Processor::VALID /* do not update if something connects to this port*/)
    , binarizationThreshold_("binarizationThreshold","Binarization threshold", 0.5, 0, std::numeric_limits<float>::max())
    , enabledProp_("enabledProp","Enabled",true)
    , thinningAlgorithm_("thinningAlgorithm", "Thinning Algorithm")
    , maxSteps_("maxSteps","Max thinning steps", 100, 0, 1000)
{
    addPort(inport_);
    addPort(sampleMaskInport_);
    addPort(outport_);
    addProperty(enabledProp_);
    addProperty(binarizationThreshold_);
    addProperty(thinningAlgorithm_);
        thinningAlgorithm_.addOption("improved", "Improved", VolumeMask::IMPROVED);
        thinningAlgorithm_.addOption("ma", "Ma", VolumeMask::MA);
        thinningAlgorithm_.addOption("chen", "Chen", VolumeMask::CHEN);
        thinningAlgorithm_.selectByValue(VolumeMask::IMPROVED);
    addProperty(maxSteps_);
}

VolumeThinning::~VolumeThinning() {
}

bool VolumeThinning::isReady() const {
    // sampleMaskInport_ is optional.
    return isInitialized() && inport_.isReady() && outport_.isReady();
}

void VolumeThinning::process() {
    const VolumeBase* maybeInvol = inport_.getData();
    const VolumeBase* sampleMask = sampleMaskInport_.getData(); //May be null
    if(!maybeInvol) {
        return;
    }
    const VolumeBase& invol = *maybeInvol;
    if(sampleMask && sampleMask->getDimensions() != invol.getDimensions()) {
        LERROR("Sample mask and input volume dimensions differ");
        tgtAssert(false, "Sample mask and input volume dimensions differ");
        return;
    }
    if(inport_.hasChanged()) {
        const VolumeMinMax* vmm = invol.getDerivedData<VolumeMinMax>();
        float currentVal = binarizationThreshold_.get();
        binarizationThreshold_.setMinValue(vmm->getMin());
        binarizationThreshold_.setMaxValue(vmm->getMax());
        binarizationThreshold_.adaptDecimalsToRange(3);
        if(vmm->getMin() > currentVal || currentVal > vmm->getMax()) {
            binarizationThreshold_.set((vmm->getMax() + vmm->getMin())/2);
        }
    }
    if(!enabledProp_.get()) {
        outport_.setData(&invol, false);
        return;
    }

    auto t_start = std::chrono::high_resolution_clock::now();

    float binarizationThreshold = invol.getRealWorldMapping().realWorldToNormalized(binarizationThreshold_.get());

    SubtaskProgressReporterCollection<3> subtaskReporters(*this, {{0.01,0.49,0.50}});

    VolumeMask m(
            std::move(binarizeVolume(invol, invol.getRealWorldMapping().realWorldToNormalized(binarizationThreshold_.get()),
                    sampleMask ? SubtaskProgressReporter(subtaskReporters.get<0>(), tgt::vec2(0,0.5)) : subtaskReporters.get<0>())),
            std::move(sampleMask
                ? boost::optional<LZ4SliceVolume<uint8_t>>(binarizeVolume(*sampleMask, 0.5, SubtaskProgressReporter(subtaskReporters.get<0>(), tgt::vec2(0.5,1.0))))
                : boost::none),
            NoFixedForeground(),
            subtaskReporters.get<1>());
    switch(thinningAlgorithm_.getValue()) {
        case VolumeMask::MA:
            m.skeletonize<VolumeMask::MA>(maxSteps_.get(), subtaskReporters.get<2>());
            break;
        case VolumeMask::CHEN:
            m.skeletonize<VolumeMask::CHEN>(maxSteps_.get(), subtaskReporters.get<2>());
            break;
        case VolumeMask::IMPROVED:
            m.skeletonize<VolumeMask::IMPROVED>(maxSteps_.get(), subtaskReporters.get<2>());
            break;
        case VolumeMask::IMPROVED_NO_LINE_PRESERVATION:
            m.skeletonize<VolumeMask::IMPROVED_NO_LINE_PRESERVATION>(maxSteps_.get(), subtaskReporters.get<2>());
            break;
    }
    VolumeBase* outvol = new Volume(m.toVolumeRAM("uint8"), &invol);

    auto t_end = std::chrono::high_resolution_clock::now();
    std::cout << "Time elapsed: " << std::chrono::duration<double>(t_end-t_start).count() << "s\n";

    outport_.setData(outvol);
}
VoreenSerializableObject* VolumeThinning::create() const {
    return new VolumeThinning();
}

} // namespace voreen
