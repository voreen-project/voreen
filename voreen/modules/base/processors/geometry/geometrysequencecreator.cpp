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

#include "geometrysequencecreator.h"

#include "voreen/core/datastructures/geometry/geometrysequence.h"

namespace voreen {

const std::string GeometrySequenceCreator::loggerCat_("voreen.base.GeometrysequenceCreator");

GeometrySequenceCreator::GeometrySequenceCreator()
    : Processor()
    , inport_(Port::INPORT, "inport", "Input Geometries",true)
    , outport_(Port::OUTPORT, "outport", "Created Sequence Output")
    , invalidOutput_(false)
    , enableClippingProp_("enableClippingprop","Enable Clipping", false)
    , clippingBBoxProp_("clippingBBoxprop","Clipping Planes")
{
    addPort(inport_);
    addPort(outport_);

    addProperty(enableClippingProp_);
    addProperty(clippingBBoxProp_);

    ON_PROPERTY_CHANGE(enableClippingProp_,GeometrySequenceCreator,adjustPropertyVisibility);
    ON_PROPERTY_CHANGE(clippingBBoxProp_,GeometrySequenceCreator,invalidateOutput);

    adjustPropertyVisibility();
}

GeometrySequenceCreator::~GeometrySequenceCreator()
{}

Processor* GeometrySequenceCreator::create() const {
    return new GeometrySequenceCreator();
}

void GeometrySequenceCreator::adjustPropertiesToInput() {
    invalidOutput_ = true;
    //create tmp geometry
    GeometrySequence* tmpSequence = new GeometrySequence();
    std::vector<const Geometry*> vec = inport_.getAllData();
    for(std::vector<const Geometry*>::iterator it = vec.begin(); it != vec.end(); it++)
            tmpSequence->addGeometry((*it)->clone().release());
    //set range
    tgt::Bounds box = tmpSequence->getBoundingBox(false);
    clippingBBoxProp_.setMinValue(box.getLLF());
    clippingBBoxProp_.setMaxValue(box.getURB());
    //TODO: Handle float (epsilon) //Note: not working with optimized pg threads.
    //if(clippingBBoxProp_.get() != box)
    //    clippingBBoxProp_.set(tgt::Bounds(clippingBBoxProp_.getMinValue(),clippingBBoxProp_.getMaxValue()));
    delete tmpSequence;
}

void GeometrySequenceCreator::adjustPropertyVisibility() {
    clippingBBoxProp_.setVisibleFlag(enableClippingProp_.get());
    invalidateOutput();
}

void GeometrySequenceCreator::invalidateOutput() {
    invalidOutput_ = true;
}


void GeometrySequenceCreator::process() {
    if(invalidOutput_) {
        GeometrySequence* output = new GeometrySequence();
        std::vector<const Geometry*> vec = inport_.getAllData();

        //clone geometries into sequence
        for(std::vector<const Geometry*>::iterator it = vec.begin(); it != vec.end(); it++) {
            if(const GeometrySequence* sequence = dynamic_cast<const GeometrySequence*>(*it)) {
                for(size_t i = 0; i < sequence->getNumGeometries() ;i++) {
                    std::unique_ptr<Geometry> tmpGeo = sequence->getGeometry(i)->clone();
                    tmpGeo->setTransformationMatrix(sequence->getTransformationMatrix() * tmpGeo->getTransformationMatrix());
                    output->addGeometry(tmpGeo.release());
                }
            } else {
                output->addGeometry((*it)->clone().release());
            }
        }

        //clipping
        if(enableClippingProp_.get()) {
            tgt::vec3 llf = clippingBBoxProp_.get().getLLF();
            tgt::vec3 urb = clippingBBoxProp_.get().getURB();
            output->clip(tgt::plane(tgt::vec3(llf.x,0.f,0.f),tgt::vec3(llf.x,0.f,1.f),tgt::vec3(llf.x,1.f,0.f)));
            output->clip(tgt::plane(tgt::vec3(0.f,llf.y,0.f),tgt::vec3(1.f,llf.y,0.f),tgt::vec3(0.f,llf.y,1.f)));
            output->clip(tgt::plane(tgt::vec3(0.f,0.f,llf.z),tgt::vec3(0.f,1.f,llf.z),tgt::vec3(1.f,0.f,llf.z)));
            output->clip(tgt::plane(tgt::vec3(urb.x,0.f,0.f),tgt::vec3(urb.x,1.f,0.f),tgt::vec3(urb.x,0.f,1.f)));
            output->clip(tgt::plane(tgt::vec3(0.f,urb.y,0.f),tgt::vec3(0.f,urb.y,1.f),tgt::vec3(1.f,urb.y,0.f)));
            output->clip(tgt::plane(tgt::vec3(0.f,0.f,urb.z),tgt::vec3(1.f,0.f,urb.z),tgt::vec3(0.f,1.f,urb.z)));

        }

        //outport gets ownership
        outport_.setData(output);
        invalidOutput_ = false;
    }
}

}  //namespace
