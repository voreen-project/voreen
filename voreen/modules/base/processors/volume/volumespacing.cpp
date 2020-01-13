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

#include "volumespacing.h"

#include "voreen/core/datastructures/volume/volumeram.h"
#include "voreen/core/datastructures/volume/volumedecorator.h"
#include "voreen/core/datastructures/meta/templatemetadata.h"

namespace voreen {

const std::string VolumeSpacing::loggerCat_("voreen.base.VolumeSpacing");

VolumeSpacing::VolumeSpacing()
    : VolumeProcessor(),
    inport_(Port::INPORT, "volumehandle.input", "Volume Input"),
    outport_(Port::OUTPORT, "volumehandle.output", "Volume Output", false),
    enableProcessing_("enabled", "Enable", false),
    mode_("mode", "Mode"),
    uniformSpacing_("uniformSpacing", "Uniform Spacing", false),
    spacingX_("spacingX", "Spacing X (mm)", 1.f, 0.000001f, 1000.f),
    spacingY_("spacingY", "Spacing Y (mm)", 1.f, 0.000001f, 1000.f),
    spacingZ_("spacingZ", "Spacing Z (mm)", 1.f, 0.000001f, 1000.f),
    reset_("reset", "Reset Spacing"),
    spacingDisplay_("spacingDisplay", "Resulting Spacing (mm)", tgt::vec3(1.0f), tgt::vec3(0.0f), tgt::vec3(1000.f))
{
    addPort(inport_);
    addPort(outport_);

    enableProcessing_.onChange(MemberFunctionCallback<VolumeSpacing>(this, &VolumeSpacing::adjustPropertyVisibility));
    addProperty(enableProcessing_);

    mode_.addOption("replace", "Replace");
    mode_.addOption("scale", "Scale");
    addProperty(mode_);

    addProperty(uniformSpacing_);

    spacingX_.setNumDecimals(6);
    spacingX_.setStepping(0.001f);
    spacingY_.setNumDecimals(6);
    spacingY_.setStepping(0.001f);
    spacingZ_.setNumDecimals(6);
    spacingZ_.setStepping(0.001f);
    addProperty(spacingX_);
    addProperty(spacingY_);
    addProperty(spacingZ_);
    reset_.onChange(MemberFunctionCallback<VolumeSpacing>(this, &VolumeSpacing::resetSpacing));
    addProperty(reset_);

    spacingDisplay_.setReadOnlyFlag(true);
    spacingDisplay_.setNumDecimals(6);
    addProperty(spacingDisplay_);

    spacingX_.onChange(
        MemberFunctionOneParameterCallback<VolumeSpacing, int>(this, &VolumeSpacing::spacingChanged, 0));
    spacingY_.onChange(
        MemberFunctionOneParameterCallback<VolumeSpacing, int>(this, &VolumeSpacing::spacingChanged, 1));
    spacingZ_.onChange(
        MemberFunctionOneParameterCallback<VolumeSpacing, int>(this, &VolumeSpacing::spacingChanged, 2));

    uniformSpacing_.onChange(
        MemberFunctionCallback<VolumeSpacing>(this, &VolumeSpacing::uniformScalingChanged));
}

Processor* VolumeSpacing::create() const {
    return new VolumeSpacing();
}

void VolumeSpacing::initialize() {
    VolumeProcessor::initialize();
    adjustPropertyVisibility();
}

void VolumeSpacing::process() {
    const VolumeBase* inputVolume = inport_.getData();

    if (!enableProcessing_.get()) {
        outport_.setData(inputVolume, false);
        spacingDisplay_.set(inputVolume->getSpacing());
        return;
    }

    tgt::vec3 spacing(spacingX_.get(), spacingY_.get(), spacingZ_.get());
    if (mode_.isSelected("scale"))
        spacing *= inputVolume->getSpacing();
    VolumeBase* outputVolume =
        new VolumeDecoratorReplaceSpacing(inputVolume, spacing);
    outport_.setData(outputVolume);
    spacingDisplay_.set(spacing);
}

void VolumeSpacing::spacingChanged(int dim) {

    if (!uniformSpacing_.get())
        return;

    if (dim == 0) {
        float xScale = spacingX_.get();
        spacingY_.set(xScale);
        spacingZ_.set(xScale);
    }
    else if (dim == 1) {
        float yScale = spacingY_.get();
        spacingX_.set(yScale);
        spacingZ_.set(yScale);
    }
    else if (dim == 2) {
        float zScale = spacingZ_.get();
        spacingX_.set(zScale);
        spacingY_.set(zScale);
    }
}

void VolumeSpacing::uniformScalingChanged() {
    if (uniformSpacing_.get())
        spacingChanged(0);
}

void VolumeSpacing::adjustPropertiesToInput() {
    const VolumeBase* inputVolume = inport_.getData();

    // only adjust spacing, if processor is disabled (= no effect on output volume)
    if (!enableProcessing_.get()) {
        if (mode_.isSelected("replace")) {
            spacingX_.set(inputVolume ? inputVolume->getSpacing().x : 1.f);
            spacingY_.set(inputVolume ? inputVolume->getSpacing().y : 1.f);
            spacingZ_.set(inputVolume ? inputVolume->getSpacing().z : 1.f);
        }
        else if (mode_.isSelected("scale")) {
            spacingX_.set(1.f);
            spacingY_.set(1.f);
            spacingZ_.set(1.f);
        }
    }
}

void VolumeSpacing::resetSpacing() {
    if (mode_.isSelected("replace")) {
        const VolumeBase* inputVolume = inport_.getData();
        if (inputVolume) {
            spacingX_.set(inputVolume->getSpacing().x);
            spacingY_.set(inputVolume->getSpacing().y);
            spacingZ_.set(inputVolume->getSpacing().z);
        }
        else {
            spacingX_.set(1.f);
            spacingY_.set(1.f);
            spacingZ_.set(1.f);
        }
    }
    else if (mode_.isSelected("scale")) {
        spacingX_.set(1.f);
        spacingY_.set(1.f);
        spacingZ_.set(1.f);
    }
}

void VolumeSpacing::adjustPropertyVisibility() {
    bool enabled = enableProcessing_.get();
    mode_.setReadOnlyFlag(!enabled);
    uniformSpacing_.setReadOnlyFlag(!enabled);
    spacingX_.setReadOnlyFlag(!enabled);
    spacingY_.setReadOnlyFlag(!enabled);
    spacingZ_.setReadOnlyFlag(!enabled);
    reset_.setReadOnlyFlag(!enabled);
}

}   // namespace
