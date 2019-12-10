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

#include "streamlinerotation.h"

#include "../../datastructures/streamlinelistdecorator.h"

namespace voreen {

const std::string StreamlineRotation::loggerCat_("voreen.flowreen.StreamlineRotation");

StreamlineRotation::StreamlineRotation()
    : Processor()
    // ports
    , streamlineInport_(Port::INPORT, "streamlineInport", "Streamline Input", false)
    , streamlineOutport_(Port::OUTPORT, "streamlineOutport", "Streamline Output", false)
    // properties
    , rotationAxisProp_("rotationAxisProp", "Axis of Rotation", tgt::vec3(1.f,0.f,0.f),tgt::vec3(-1.f),tgt::vec3(1.f))
    , rotationAngleProp_("rotationAngleProp", "Angle of Rotation", 0.f, 0.f, 359.99f)
    , linkAxisAndHandleProp_("linkAxisAndHandleProp","Link Handle with Axis?", true,Processor::VALID,Property::LOD_ADVANCED)
    // alignment
    , alignmentIsAppliedProp_("alignmentIsAppliedProp", "Is Alignment Applied?", false, Processor::VALID) //no invalidate here
    , currentAlignmentSettingsAreAppliedProp_("currentAlignmentSettingsAreAppliedProp", "Current Settings Applied?", false, Processor::VALID)  //no invalidate here
    , linkedHandleNormalProp_("linkedHandleNormalProp", "Link with PlaneManipulation",tgt::vec3(1.f,0.f,0.f),tgt::vec3(-1.f),tgt::vec3(1.f), Processor::VALID) //no invalidate here
    , alignToPlaneProp_("alignToPlaneProp", "Align to Plane",Processor::VALID) //no invalidate here
    , applyAlignmentProp_("applyAlignmentProp", "Apply Alignment")
    , clearAlignmentProp_("clearAlignmentProp", "Clear Alignment")
    , transMatOfLastStreamlineList_(tgt::mat4::identity)
    , currentRotationMatrix_(tgt::mat4::identity)
    , alignmentMatrix_(tgt::mat4::identity)
    , lastConfigurationVelocityRotation_(tgt::mat4::identity)
    , lastConfigurationStreamlineTransform_(tgt::mat4::identity)
    // to be linked
    , linkedMatrixProp_("linkedMatrixProp","Link with VolumeTransform",tgt::mat4::identity,tgt::mat4(-1000000),tgt::mat4(1000000))
{
    addPort(streamlineInport_);
        streamlineInport_.onNewData(MemberFunctionCallback<StreamlineRotation>(this, &StreamlineRotation::onNewDataStreamlinePort));
    addPort(streamlineOutport_);

    //generell rotation
    addProperty(rotationAxisProp_);
        rotationAxisProp_.setGroupID("rotate");
        rotationAxisProp_.onChange(MemberFunctionCallback<StreamlineRotation>(this, &StreamlineRotation::rotationAxisOnChange));
    addProperty(rotationAngleProp_);
        rotationAngleProp_.setGroupID("rotate");
        rotationAngleProp_.onChange(MemberFunctionCallback<StreamlineRotation>(this, &StreamlineRotation::rotationAngleOnChange));
    addProperty(linkAxisAndHandleProp_);
        linkAxisAndHandleProp_.setGroupID("rotate");
        linkAxisAndHandleProp_.onChange(MemberFunctionCallback<StreamlineRotation>(this, &StreamlineRotation::rotationAngleOnChange));
    setPropertyGroupGuiName("rotate", "Rotation");

    //align to plane
    addProperty(alignmentIsAppliedProp_);
        alignmentIsAppliedProp_.setGroupID("align");
        alignmentIsAppliedProp_.setReadOnlyFlag(true);
    addProperty(currentAlignmentSettingsAreAppliedProp_);
        currentAlignmentSettingsAreAppliedProp_.setGroupID("align");
        currentAlignmentSettingsAreAppliedProp_.setReadOnlyFlag(true);
    addProperty(linkedHandleNormalProp_);
        linkedHandleNormalProp_.setGroupID("align");
        linkedHandleNormalProp_.onChange(MemberFunctionCallback<StreamlineRotation>(this, &StreamlineRotation::alignmentSettingsOnChange));
        linkedHandleNormalProp_.onChange(MemberFunctionCallback<StreamlineRotation>(this, &StreamlineRotation::linkedHandleOnChange));
    addProperty(alignToPlaneProp_);
        alignToPlaneProp_.addOption("xNormal","YZ-Plane",YZ_PLANE,tgt::col4(255,0,0,255));
        alignToPlaneProp_.addOption("yNormal","XZ-Plane",XZ_PLANE,tgt::col4(0,255,0,255));
        alignToPlaneProp_.addOption("zNormal","XY-Plane",XY_PLANE,tgt::col4(0,0,255,255));
        alignToPlaneProp_.setGroupID("align");
        alignToPlaneProp_.onChange(MemberFunctionCallback<StreamlineRotation>(this, &StreamlineRotation::alignmentSettingsOnChange));
    addProperty(applyAlignmentProp_);
        applyAlignmentProp_.setGroupID("align");
        applyAlignmentProp_.onChange(MemberFunctionCallback<StreamlineRotation>(this, &StreamlineRotation::applyAlignmentOnChange));
    addProperty(clearAlignmentProp_);
        clearAlignmentProp_.setGroupID("align");
        clearAlignmentProp_.setReadOnlyFlag(true);
        clearAlignmentProp_.onChange(MemberFunctionCallback<StreamlineRotation>(this, &StreamlineRotation::clearAlignmentOnChange));
    setPropertyGroupGuiName("align", "Alignment");

    //result
    addProperty(linkedMatrixProp_);
        linkedMatrixProp_.setReadOnlyFlag(true);

    //get matrix right
    rotationAngleOnChange();
}

StreamlineRotation::~StreamlineRotation() {
}

void StreamlineRotation::serialize(Serializer& s) const {
    Processor::serialize(s);
    s.serialize("TransMatOfLastStreamlineList",transMatOfLastStreamlineList_);
    s.serialize("AlignmentMatrix",alignmentMatrix_);
    s.serialize("LastConfigurationVelocityRotation",lastConfigurationVelocityRotation_);
    s.serialize("LastConfigurationStreamlineTransform",lastConfigurationStreamlineTransform_);
}

void StreamlineRotation::deserialize(Deserializer& d) {
    Processor::deserialize(d);
    d.optionalDeserialize("TransMatOfLastStreamlineList",transMatOfLastStreamlineList_,tgt::mat4::identity);
    d.optionalDeserialize("AlignmentMatrix",alignmentMatrix_,tgt::mat4::identity);
    d.optionalDeserialize("LastConfigurationVelocityRotation",lastConfigurationVelocityRotation_,tgt::mat4::identity);
    d.optionalDeserialize("LastConfigurationStreamlineTransform",lastConfigurationStreamlineTransform_,tgt::mat4::identity);
}

void StreamlineRotation::process() {
    StreamlineListBase* inputStreamlines = const_cast<StreamlineListBase*>(streamlineInport_.getData());
    //calculate rotation matrix
    tgt::mat4 newRotationPart = currentRotationMatrix_ * alignmentMatrix_ * lastConfigurationVelocityRotation_;
    //set ouput
    streamlineOutport_.setData(new StreamlineListDecoratorReplaceTransformation(inputStreamlines, linkedMatrixProp_.get(),newRotationPart));
}

//---------------------------------------------------
//      Callbacks
//---------------------------------------------------
    //alignment
void StreamlineRotation::alignmentSettingsOnChange() {
    currentAlignmentSettingsAreAppliedProp_.set(false);
}

void StreamlineRotation::applyAlignmentOnChange() {
    //backup old alignment
    lastConfigurationStreamlineTransform_ = linkedMatrixProp_.get();
    lastConfigurationVelocityRotation_ = currentRotationMatrix_ * alignmentMatrix_ * lastConfigurationVelocityRotation_;

    //calc new alignment
    tgt::vec3 currentAxis;
    switch(alignToPlaneProp_.getValue()) {
    case YZ_PLANE:
        currentAxis = tgt::vec3(1.f,0.f,0.f);
        break;
    case XZ_PLANE:
        currentAxis = tgt::vec3(0.f,1.f,0.f);
        break;
    case XY_PLANE:
        currentAxis = tgt::vec3(0.f,0.f,1.f);
        break;
    default:
        tgtAssert(false,"unknown alignment");
        break;
    }

    tgt::vec3 crossProduct = tgt::cross(linkedHandleNormalProp_.get(),currentAxis);

    //check if product is valid
    if(crossProduct != tgt::vec3::zero) {
        alignmentMatrix_ = tgt::mat4::createRotation(std::acos(tgt::dot(linkedHandleNormalProp_.get(),currentAxis)),crossProduct);
    } else {
        tgt::vec3 rotationAxis(1.f,0.f,0.f);
        if(linkedHandleNormalProp_.get().x != 0.f) {
            rotationAxis = tgt::vec3(linkedHandleNormalProp_.get().y,-linkedHandleNormalProp_.get().x,0);
        }
        alignmentMatrix_ = tgt::mat4::createRotation(std::acos(tgt::dot(linkedHandleNormalProp_.get(),currentAxis)),rotationAxis);
    }

    //update handle axis
    linkedHandleNormalProp_.set((currentRotationMatrix_ * alignmentMatrix_ * lastConfigurationVelocityRotation_*tgt::vec4(linkedHandleNormalProp_.get(),1.f)).xyz());
    //rotationAxisProp_.blockCallbacks(true);
    //rotationAxisProp_.set(linkedHandleNormalProp_.get());
    //rotationAxisProp_.blockCallbacks(false);

    //force update
    rotationAngleOnChange();

    clearAlignmentProp_.setReadOnlyFlag(false);
    alignmentIsAppliedProp_.set(true);
    currentAlignmentSettingsAreAppliedProp_.set(true);
}

void StreamlineRotation::clearAlignmentOnChange() {
    clearAlignmentProp_.setReadOnlyFlag(true);
    alignmentIsAppliedProp_.set(false);
    currentAlignmentSettingsAreAppliedProp_.set(false);
    //reset all
    alignmentMatrix_ = tgt::mat4::identity;
    lastConfigurationVelocityRotation_ = streamlineInport_.getData()->getVelocityTransformMatrix();
    lastConfigurationStreamlineTransform_ = streamlineInport_.getData()->getListTransformMatrix();
    //update linked matrix
    updateLinkedMatrix();
}

void StreamlineRotation::rotationAngleOnChange() {
    currentRotationMatrix_ = tgt::mat4::createRotationDegree(rotationAngleProp_.get(), rotationAxisProp_.get());
    //update linked matrix
    updateLinkedMatrix();
}

void StreamlineRotation::rotationAxisOnChange() {
    rotationAngleProp_.set(0.f);
    if(linkAxisAndHandleProp_.get()) {
       linkedHandleNormalProp_.set(rotationAxisProp_.get());
    }
}

void StreamlineRotation::linkedHandleOnChange() {
    //if they are linked
    if(linkAxisAndHandleProp_.get()) {
        rotationAxisProp_.set(linkedHandleNormalProp_.get());
    }
}

void StreamlineRotation::updateLinkedMatrix() {
    //needed?
    if(!isInitialized()) return;

    const StreamlineListBase* inputList = streamlineInport_.getData();
    if(!inputList) {
        linkedMatrixProp_.set(tgt::mat4::identity);
        return;
    }

    //rotate around the bb center of the data
    tgt::vec3 dataCenter = inputList->getOriginalWorldBounds().transform(lastConfigurationStreamlineTransform_).center();
    tgt::mat4 dataCenterTransformRotationMatrix = tgt::mat4::createTranslation(dataCenter) * currentRotationMatrix_ * alignmentMatrix_ * tgt::mat4::createTranslation(-dataCenter) * lastConfigurationStreamlineTransform_;

    linkedMatrixProp_.set(dataCenterTransformRotationMatrix);
}

void StreamlineRotation::onNewDataStreamlinePort() {
    if(transMatOfLastStreamlineList_ == streamlineInport_.getData()->getListTransformMatrix()) {
        // simply process, which will be triggered trough the port invalidation
    } else {
        transMatOfLastStreamlineList_ = streamlineInport_.getData()->getListTransformMatrix();
        lastConfigurationVelocityRotation_ = streamlineInport_.getData()->getVelocityTransformMatrix();
        lastConfigurationStreamlineTransform_ = streamlineInport_.getData()->getListTransformMatrix();
        //reset all settings
        rotationAngleProp_.set(0.f);
        clearAlignmentOnChange();
    }
}

}   // namespace
