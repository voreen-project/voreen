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

#include "volumemorphology.h"

#include "voreen/core/datastructures/volume/volumeram.h"
#include "voreen/core/datastructures/volume/volume.h"
#include "voreen/core/datastructures/volume/operators/volumeoperatormorphology.h"

#include "tgt/tgt_math.h"

namespace voreen {

const std::string VolumeMorphology::loggerCat_("voreen.VolumeMorphology");

VolumeMorphology::VolumeMorphology()
    : CachingVolumeProcessor()
    , inport_(Port::INPORT, "volumehandle.input", "Volume Input")
    , outport_(Port::OUTPORT, "volumehandle.output", "Volume Output", false)
    , enableProcessing_("enableProcessing", "Enable")
    , morphologicOperator_("morphologicOperator", "Operator")
    , morphologicOperatorShape_("morphologicOperatorShape", "Shape")
    , kernelSize_("kernelSize", "Kernel Size")
    , forceUpdate_(true)
{
    addPort(inport_);
    addPort(outport_);

    morphologicOperator_.addOption("dilation", "Dilation");
    morphologicOperator_.addOption("erosion", "Erosion");

    morphologicOperatorShape_.addOption("cube", "Cube", Shape::CUBE);
    morphologicOperatorShape_.addOption("sphere", "Sphere", Shape::SPHERE);


    kernelSize_.addOption("3",  "3x3x3",    3);
    kernelSize_.addOption("5",  "5x5x5",    5);
    kernelSize_.addOption("7",  "7x7x7",    7);
    kernelSize_.addOption("9",  "9x9x9",    9);
    kernelSize_.addOption("15", "15x15x15", 15);
    kernelSize_.addOption("25", "25x25x25", 25);
    kernelSize_.addOption("35", "35x35x35", 35);
    kernelSize_.addOption("45", "45x45x45", 45);

    enableProcessing_.onChange(MemberFunctionCallback<VolumeMorphology>(this, &VolumeMorphology::forceUpdate));
    morphologicOperator_.onChange(MemberFunctionCallback<VolumeMorphology>(this, &VolumeMorphology::forceUpdate));
    morphologicOperatorShape_.onChange(MemberFunctionCallback<VolumeMorphology>(this, &VolumeMorphology::forceUpdate));
    kernelSize_.onChange(MemberFunctionCallback<VolumeMorphology>(this, &VolumeMorphology::forceUpdate));

    addProperty(enableProcessing_);
    addProperty(morphologicOperator_);
    addProperty(morphologicOperatorShape_);
    addProperty(kernelSize_);
    inport_.onNewData(MemberFunctionCallback<VolumeMorphology>(this, &VolumeMorphology::applyOperator));
}

VolumeMorphology::~VolumeMorphology() {}

Processor* VolumeMorphology::create() const {
    return new VolumeMorphology();
}

void VolumeMorphology::process() {
    if (!enableProcessing_.get()) {
        outport_.setData(const_cast<VolumeBase*>(inport_.getData()), false);
    }
    else if (forceUpdate_) {
        applyOperator();
    }
}

// private methods
//

void VolumeMorphology::forceUpdate() {
    forceUpdate_ = true;
}

void VolumeMorphology::applyOperator() {
    if(!inport_.hasData())
        return;

    forceUpdate_ = false;

    if (inport_.getData()->getRepresentation<VolumeRAM>()) {
        const VolumeBase* input = inport_.getData();
        Volume* transformed = 0;

        if (morphologicOperator_.isSelected("dilation")) {
            switch(morphologicOperatorShape_.getValue()) {
                case Shape::CUBE:
                    transformed = VolumeOperatorCubeDilation::APPLY_OP(input, kernelSize_.getValue(), (ProgressReporter*)this);
                    break;
                case Shape::SPHERE:
                    transformed = VolumeOperatorSphereDilation::APPLY_OP(input, kernelSize_.getValue(), (ProgressReporter*)this);
                    break;
                default:
                    tgtAssert(false, "Invalid shape");
            }
        }
        else if (morphologicOperator_.isSelected("erosion")) {
            switch(morphologicOperatorShape_.getValue()) {
                case Shape::CUBE:
                    transformed = VolumeOperatorCubeErosion::APPLY_OP(input, kernelSize_.getValue(), (ProgressReporter*)this);
                    break;
                case Shape::SPHERE:
                    transformed = VolumeOperatorSphereErosion::APPLY_OP(input, kernelSize_.getValue(), (ProgressReporter*)this);
                    break;
                default:
                    tgtAssert(false, "Invalid shape");
            }
        }
        else {
            LERROR("Unknown operator: " << morphologicOperator_.get());
        }

        outport_.setData(transformed);
    }
    else {
        outport_.setData(0);
    }
}
}   // namespace
