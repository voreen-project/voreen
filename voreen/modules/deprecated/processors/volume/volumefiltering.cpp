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

#include "volumefiltering.h"

#include "voreen/core/datastructures/volume/volumeram.h"
#include "voreen/core/datastructures/volume/volume.h"
#include "voreen/core/datastructures/octree/volumeoctree.h"
#include "voreen/core/datastructures/callback/memberfunctioncallback.h"
#include "voreen/core/datastructures/volume/operators/volumeoperatormedian.h"
#include "voreen/core/datastructures/volume/operators/volumeoperatorgaussian.h"

#include "tgt/tgt_math.h"

namespace voreen {

const std::string VolumeFiltering::loggerCat_("voreen.VolumeFiltering");

VolumeFiltering::VolumeFiltering()
    : CachingVolumeProcessor()
    , inport_(Port::INPORT, "volumehandle.input", "Volume Input")
    , outport_(Port::OUTPORT, "volumehandle.output", "Volume Output", false)
    , enableProcessing_("enableProcessing", "Enable")
    , filteringOperator_("filteringOperator", "Operator")
    , kernelSize_("kernelSize", "Kernel Size")
    , deriveSigmaFromKernelSize_("sigmafromkernel", "Derive Sigma from Kernel Size", true)
    , sigma_("sigma", "Sigma", 1.f, 0.1f, 20.f)
    , forceUpdate_(true)
{
    addPort(inport_);
    addPort(outport_);

    filteringOperator_.addOption("median", "Median");
    filteringOperator_.addOption("gaussian", "Gaussian");
    ON_CHANGE(filteringOperator_, VolumeFiltering, adjustPropertyVisibility);

    kernelSize_.addOption("3",  "3x3x3",    3);
    kernelSize_.addOption("5",  "5x5x5",    5);
    kernelSize_.addOption("7",  "7x7x7",    7);
    kernelSize_.addOption("9",  "9x9x9",    9);
    kernelSize_.addOption("11", "11x11x11", 11);
    kernelSize_.addOption("13", "13x13x13", 13);
    kernelSize_.addOption("15", "15x15x15", 15);

    enableProcessing_.onChange(MemberFunctionCallback<VolumeFiltering>(this, &VolumeFiltering::forceUpdate));
    filteringOperator_.onChange(MemberFunctionCallback<VolumeFiltering>(this, &VolumeFiltering::forceUpdate));
    kernelSize_.onChange(MemberFunctionCallback<VolumeFiltering>(this, &VolumeFiltering::forceUpdate));
    inport_.onNewData(MemberFunctionCallback<VolumeFiltering>(this, &VolumeFiltering::forceUpdate));

    addProperty(enableProcessing_);
    addProperty(filteringOperator_);
    addProperty(kernelSize_);

    ON_CHANGE(kernelSize_, VolumeFiltering, adjustPropertyVisibility);
    ON_CHANGE(deriveSigmaFromKernelSize_, VolumeFiltering, adjustPropertyVisibility);
    ON_CHANGE(sigma_, VolumeFiltering, forceUpdate);
    addProperty(deriveSigmaFromKernelSize_);
    addProperty(sigma_);

    adjustPropertyVisibility();
}

VolumeFiltering::~VolumeFiltering() {
}

Processor* VolumeFiltering::create() const {
    return new VolumeFiltering();
}

void VolumeFiltering::process() {
    if (!enableProcessing_.get()) {
        outport_.setData(const_cast<VolumeBase*>(inport_.getData()), false);
    }
    else if (forceUpdate_) {

        // check for large data and an available VolumeRAM representation
        if (inport_.getData()->hasRepresentation<VolumeOctree>() && !inport_.getData()->hasRepresentation<VolumeRAM>()) {
            // volume octree that does not already provide a VolumeRAM representation -> do nothing
            LERROR("Volume has an octree representation, but not a VolumeRAM representation. Volume processing is not applicable.");
            outport_.clear();
            return;
        }

        applyOperator();
    }
}

// private methods
//

void VolumeFiltering::forceUpdate() {
    forceUpdate_ = true;
}

void VolumeFiltering::applyOperator() {
    if(!inport_.hasData())
        return;

    forceUpdate_ = false;

    if (inport_.getData()->getRepresentation<VolumeRAM>()) {
        const VolumeBase* input = inport_.getData();
        Volume* transformed = 0; // = input->clone();

        if (filteringOperator_.isSelected("median")) {
            //volOpMedian.setProgressBar(progressBar_);
            transformed = VolumeOperatorMedian::APPLY_OP(input, kernelSize_.getValue(), this);
        }
        else if (filteringOperator_.isSelected("gaussian")) {
            transformed = VolumeOperatorGaussian::APPLY_OP(input, static_cast<size_t>(kernelSize_.getValue()), sigma_.get(), this);
        }
        else {
            LERROR("Unknown operator: " << filteringOperator_.get());
        }

        outport_.setData(transformed);
    }
    else {
        outport_.setData(0);
    }
}

void VolumeFiltering::adjustPropertyVisibility() {
    deriveSigmaFromKernelSize_.setVisibleFlag(filteringOperator_.isSelected("gaussian"));
    sigma_.setVisibleFlag(filteringOperator_.isSelected("gaussian"));

    sigma_.setReadOnlyFlag(deriveSigmaFromKernelSize_.get());
    if (deriveSigmaFromKernelSize_.get()) {
        size_t kernelRadius = static_cast<size_t>(kernelSize_.getValue()) / 2;
        float sigma = static_cast<float>(kernelRadius) / 2.5f;
        sigma_.set(sigma);
    }
}

}   // namespace
