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

#include "volumetransformation.h"

#include "voreen/core/datastructures/volume/volumeram.h"
#include "voreen/core/datastructures/volume/volume.h"
#include "voreen/core/datastructures/volume/volumedecorator.h"
#include "voreen/core/datastructures/meta/templatemetadata.h"
#include "voreen/core/datastructures/callback/lambdacallback.h"

#include "tgt/logmanager.h"

namespace voreen {

const std::string VolumeTransformation::loggerCat_("voreen.base.VolumeTransformation");

VolumeTransformation::VolumeTransformation()
    : VolumeProcessor()
    , inport_(Port::INPORT, "volumehandle.input", "Volume Input")
    , outport_(Port::OUTPORT, "volumehandle.output", "Volume Output", false)
    , enableProcessing_("enableProcessing", "Enable")
    , transformationMode_("transformationMode", "Transformation Mode")
    , sourceCoordinateSystem_("sourceCoordinateSystem", "Source Coordinate System")
    , transformationType_("transformationType", "Transformation Type")
    , transformMatrix_("transformMatrix", "Transformation Matrix", tgt::mat4::identity, tgt::mat4(-1e10), tgt::mat4(1e10))
    , changingTransformationInternally_(false)
    // Rotation properties
    , rotationAngle_("rotationAngle", "Rotation Angle", 0.0f, 0.0f, 2*tgt::PIf)
    , rotationAxis_("rotationAxis", "Rotation Axis", tgt::vec3(1.f,0.f,0.f), tgt::vec3(-1.f), tgt::vec3(1.f))
    , rotationReferencePoint_("rotationReferencePoint", "Rotation Reference Point", tgt::vec3(0.f), tgt::vec3(-1e10), tgt::vec3(1e10))
    , translationAmount_("translationAmount", "Translation Amount", tgt::vec3(0.f), tgt::vec3(-1e10), tgt::vec3(1e10))
{
    addPort(inport_);
    addPort(outport_);

        addProperty(enableProcessing_);
            enableProcessing_.setGroupID("general");

        addProperty(transformationMode_);
            transformationMode_.setGroupID("general");
            transformationMode_.addOption("concatenate", "Concatenate");
            transformationMode_.addOption("replace",     "Replace");
            ON_CHANGE(transformationMode_, VolumeTransformation, adaptToTransformationMode);

        addProperty(sourceCoordinateSystem_);
            sourceCoordinateSystem_.setGroupID("general");
            sourceCoordinateSystem_.addOption("voxel-coordinates",   "Voxel Coordinates");
            sourceCoordinateSystem_.addOption("texture-coordinates", "Texture Coordinates");
            sourceCoordinateSystem_.addOption("volume-coordinates",  "Volume/Physical Coordinates");
            sourceCoordinateSystem_.addOption("world-coordinates",   "World Coordinates");
            sourceCoordinateSystem_.select("volume-coordinates");

        addProperty(transformMatrix_);
            transformMatrix_.setGroupID("general");
            ON_CHANGE_LAMBDA(transformMatrix_, [this] {
                    if(!changingTransformationInternally_) {
                        transformationType_.selectByValue(ARBITRARY);
                    }
                    });

        addProperty(transformationType_);
            transformationType_.setGroupID("general");
            transformationType_.addOption("rotation", "Rotation", ROTATION);
            transformationType_.addOption("translation", "Translation", TRANSLATION);
            transformationType_.addOption("arbitrary", "Arbitrary", ARBITRARY);
            transformationType_.selectByValue(ARBITRARY); //Arbitrary Transformation by default
            ON_CHANGE(transformationType_, VolumeTransformation, adaptToTransformationType);
    setPropertyGroupGuiName("general", "General configuration");


        addProperty(rotationAngle_);
            rotationAngle_.setGroupID("rotation");
            ON_CHANGE(rotationAngle_, VolumeTransformation, setRotationMatrix);

        addProperty(rotationAxis_);
            rotationAxis_.setGroupID("rotation");
            ON_CHANGE(rotationAxis_, VolumeTransformation, setRotationMatrix);

        addProperty(rotationReferencePoint_);
            rotationReferencePoint_.setGroupID("rotation");
            ON_CHANGE(rotationReferencePoint_, VolumeTransformation, setRotationMatrix);

    setPropertyGroupGuiName("rotation", "Rotation configuration");


        addProperty(translationAmount_);
            translationAmount_.setGroupID("translation");
            ON_CHANGE(translationAmount_, VolumeTransformation, setTranslationMatrix);

    setPropertyGroupGuiName("translation", "Translation configuration");
}

VolumeTransformation::~VolumeTransformation() {}

Processor* VolumeTransformation::create() const {
    return new VolumeTransformation();
}

void VolumeTransformation::initialize() {
    VolumeProcessor::initialize();
    adaptToTransformationMode();
    adaptToTransformationType();
}

void VolumeTransformation::process() {
    const VolumeBase* inputVolume = inport_.getData();
    if (!enableProcessing_.get()) {
        outport_.setData(inputVolume, false);
        return;
    }

    tgt::mat4 outputTrafo = tgt::mat4::identity;
    if (transformationMode_.isSelected("concatenate")) {
        // transform volume into world-coords using current transformation, then apply specified transformation
        outputTrafo = transformMatrix_.get() * inputVolume->getPhysicalToWorldMatrix();
    }
    else if (transformationMode_.isSelected("replace")) {
        // transform from selected source coords into physical coords, then apply specified transformation
        tgt::mat4 toPhysical = tgt::mat4::identity;
        if (sourceCoordinateSystem_.isSelected("voxel-coordinates"))
            toPhysical = inputVolume->getPhysicalToVoxelMatrix();
        else if (sourceCoordinateSystem_.isSelected("volume-coordinates"))
            toPhysical = tgt::mat4::identity; // already in physical coordinates
        else if (sourceCoordinateSystem_.isSelected("texture-coordinates"))
            toPhysical = inputVolume->getPhysicalToTextureMatrix();
        else if (sourceCoordinateSystem_.isSelected("world-coordinates"))
            toPhysical = inputVolume->getPhysicalToWorldMatrix();
        else {
            LERROR("unknown source coordinates system selected: " << sourceCoordinateSystem_.get());
        }

        outputTrafo = transformMatrix_.get() * toPhysical;
    }
    else {
        LERROR("unknown transformationMode selected: " << transformationMode_.get());
    }

    VolumeBase* outputVolume =
        new VolumeDecoratorReplaceTransformation(inport_.getData(), outputTrafo);
    outport_.setData(outputVolume);
}

void VolumeTransformation::adaptToTransformationMode() {
    sourceCoordinateSystem_.setVisibleFlag(transformationMode_.isSelected("replace"));
}

void VolumeTransformation::adaptToTransformationType() {
    transformMatrix_.setReadOnlyFlag(transformationType_.getValue() != ARBITRARY);

    rotationReferencePoint_.setVisibleFlag(transformationType_.getValue() == ROTATION);
    rotationAxis_.setVisibleFlag(transformationType_.getValue() == ROTATION);
    rotationAngle_.setVisibleFlag(transformationType_.getValue() == ROTATION);
    setRotationMatrix(); //Update matrix if necessary

    translationAmount_.setVisibleFlag(transformationType_.getValue() == TRANSLATION);
    setTranslationMatrix(); //Update matrix if necessary
}

void VolumeTransformation::setTransformationMatrix(const tgt::mat4& transformation) {
    changingTransformationInternally_ = true;
    transformMatrix_.set(transformation);
    changingTransformationInternally_ = false;
}
void VolumeTransformation::setRotationMatrix() {
    if(transformationType_.getValue() != ROTATION) {
        // Do not change the value of transformMatrix if rotation is not selected
        return;
    }
    if(rotationAxis_.get() == tgt::vec3::zero) {
        LERROR("Rotation axis is zero");
        return;
    }
    // Axis will be normalized in createRotation
    tgt::mat4 rotation = tgt::mat4::createTranslation(-rotationReferencePoint_.get())
            * tgt::mat4::createRotation(rotationAngle_.get(), rotationAxis_.get())
            * tgt::mat4::createTranslation(rotationReferencePoint_.get());

    setTransformationMatrix(rotation);
}

void VolumeTransformation::setTranslationMatrix() {
    if(transformationType_.getValue() != TRANSLATION) {
        // Do not change the value of transformMatrix if rotation is not selected
        return;
    }
    // Axis will be normalized in createRotation
    tgt::mat4 translation = tgt::mat4::createTranslation(translationAmount_.get());

    setTransformationMatrix(translation);
}

}   // namespace
