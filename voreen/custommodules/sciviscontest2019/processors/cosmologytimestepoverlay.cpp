/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2024 University of Muenster, Germany,                        *
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

#include "cosmologytimestepoverlay.h"

#include "tgt/textureunit.h"
#include "tgt/event/mouseevent.h"
#include "tgt/immediatemode/immediatemode.h"
#include <sstream>

using tgt::TextureUnit;

namespace voreen {
    
CosmologyTimeStepOverlay::CosmologyTimeStepOverlay()
    : ImageProcessor("image/compositor")
    , imageInport_(Port::INPORT, "image.in", "Image Inpput")
    , inport_(Port::INPORT,   "particlehandle.output", "Particle Data Output")
    , outport_(Port::OUTPORT, "image.out", "Image Output")
	, enabled_("enabled", "Enable overlay", true)
    , timeStep_         ("timeStep",           "Time Step", 0.0f, 0.0f, 624.0f)
    , font_("font", "Font")
    , offset_("offset", "Offset", tgt::ivec2::zero, tgt::ivec2(-2000), tgt::ivec2(2000))
    , readOnly_("readOnly", "Read only", false)
    , multiView_("multiView", "Toggle Views independently", false)
    , overView_("overView", "Show in Overview", true)
    , view1_("view1", "Show in view 1", true)
    , view2_("view2", "Show in view 2", true)
    , view3_("view3", "Show in view 3", true)
    , view4_("view4", "Show in view 4", true)
	, currentSelectedView_("selectedView", "currently Selected View", 0, 0 , 4)

{
    addPort(imageInport_);
    addPort(inport_);
    addPort(outport_);

	addProperty(enabled_);
    addProperty(timeStep_);
    addProperty(font_);
    addProperty(offset_);
    addProperty(readOnly_);

    addProperty(multiView_);
	addProperty(currentSelectedView_);
    addProperty(overView_);
    addProperty(view1_);
    addProperty(view2_);
    addProperty(view3_);
    addProperty(view4_);

	currentSelectedView_.setVisibleFlag(false);

	overView_.setGroupName("Selected Views");
	overView_.setGroupID("selected");
	view1_.setGroupID("selected");
	view2_.setGroupID("selected");
	view3_.setGroupID("selected");
	view4_.setGroupID("selected");
	setPropertyGroupVisible("selected", false);

	ON_CHANGE_LAMBDA(multiView_, [this] {setPropertyGroupVisible("selected", multiView_.get()); });

	timeStep_.setStepping(1.0f);
	timeStep_.adaptDecimalsToRange(0);

    isDragging_ = false;
}

Processor* CosmologyTimeStepOverlay::create() const {
    return new CosmologyTimeStepOverlay();
}

void CosmologyTimeStepOverlay::initialize() {
    ImageProcessor::initialize();
}

void CosmologyTimeStepOverlay::deinitialize() {
    ImageProcessor::deinitialize();
}


static tgt::vec3 translate(tgt::mat4 m, tgt::vec2 v){
        tgt::vec4 v2 = m*tgt::vec4(v, 0.0f, 1.0);
        return tgt::vec3(v2.x*0.5f+0.5f, v2.y*0.5f+0.5f, 0.0f);
}

void CosmologyTimeStepOverlay::onEvent(tgt::Event *ev){
    if (readOnly_.get()){
        ImageProcessor::onEvent(ev);
        return;
    }
    tgt::MouseEvent *me = dynamic_cast<tgt::MouseEvent*>(ev);
    if (!me){
        ImageProcessor::onEvent(ev);
        return;
    }
    tgt::ivec2 coordi = me->coord();
    tgt::vec2 coord = tgt::vec2(coordi)/tgt::vec2(imageInport_.getSize())*tgt::vec2(2.0f, 2.0f)-tgt::vec2::one;
    coord.y*= -1.0f;

    MatStack.pushMatrix();
    
    
    MatStack.scale(0.8f, 0.5f, 1.0);
    MatStack.translate(0.0f, -0.6f, 0.0f);
    MatStack.translate(tgt::vec3(tgt::vec2(offset_.get())/tgt::vec2(outport_.getSize()), 0));

    tgt::mat4 M = MatStack.getModelViewMatrix();
    tgt::mat4 Mi;
    M.invert(Mi);
    if (me->button() == tgt::MouseEvent::MOUSE_BUTTON_LEFT && (me->action() == tgt::MouseEvent::RELEASED||me->action() == tgt::MouseEvent::EXIT) && isDragging_){
        isDragging_ = false;
        me->accept();
        MatStack.popMatrix();
        return;
    }


    float timestep = ((Mi*tgt::vec4(coord, 0.0f, 1.0f)).x + 1.0f) * 312.0f;
    float y = (Mi*tgt::vec4(coord, 0.0f, 1.0f)).y;
    MatStack.popMatrix();
    if (me->button() == tgt::MouseEvent::MOUSE_BUTTON_LEFT && me->action() == tgt::MouseEvent::PRESSED && y < -0.9 && y > -1.0f && timestep > -625.0f && timestep < 625.0f){
        isDragging_ = true;
    }
    if (timestep > -1 && timestep < 625 && isDragging_){
        timestep = tgt::clamp(timestep, 0.0f, 624.0f);
        timeStep_.set(timestep);
        me->accept();
        return;
    }
    ImageProcessor::onEvent(ev);

}

void CosmologyTimeStepOverlay::process() {

    tgtAssert(outport_.isReady(), "Outport not ready");
    tgtAssert(imageInport_.isReady(), "Inport not ready");

    

    //same code as in ImageOverlay
    //          |
    //          v
    outport_.activateTarget();
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    glActiveTexture(GL_TEXTURE0);
    glBindTexture(GL_TEXTURE_2D, imageInport_.getColorTexture()->getId());
    renderQuad();
    glBindTexture(GL_TEXTURE_2D, 0);

	if (!enabled_.get()){
		outport_.deactivateTarget();
		return;
	}

	bool enabledFlag = true;
	switch (currentSelectedView_.get()) {
	case 0:
		enabledFlag = overView_.get();
		break;
	case 1:
		enabledFlag = view1_.get();
		break;
	case 2:
		enabledFlag = view2_.get();
		break;
	case 3:
		enabledFlag = view3_.get();
		break;
	case 4:
		enabledFlag = view4_.get();
		break;
	}

	if (!enabledFlag) {
		outport_.deactivateTarget();
		return;
	}
		

    if (!inport_.getData())
        return;

    float lineWidth = 4.0f/imageInport_.getSize().x/0.8f;
    glDisable(GL_DEPTH_TEST);

    MatStack.pushMatrix();
    
    MatStack.scale(0.8f, 0.5f, 1.0);
    MatStack.translate(0.0f, -0.6f, 0.0f);
    MatStack.translate(tgt::vec3(tgt::vec2(offset_.get())/tgt::vec2(outport_.getSize()), 0));
    tgt::mat4 m = MatStack.getModelViewMatrix();
    
    glBegin(GL_QUADS);{
        glColor3f(47/255.0f,51/255.0f,47/255.0f);
        glVertex2f(-1, -0.9);
        glVertex2f(+1, -0.9);
        glVertex2f(+1, -1);
        glVertex2f(-1, -1);
    }glEnd();

    std::vector<CMParticleDataTimeSlice*> slices = inport_.getData()->particleDataTimeSlices();
    glBegin(GL_TRIANGLES);{
        glColor3f(94/255.0f,180/255.0f,214/255.0f);
        for(auto slice: slices){
            float ts = slice->getTimeStep();
            glVertex2f(ts*2-1-0.5f*lineWidth, -0.9); 
            glVertex2f(ts*2-1-0.5f*lineWidth, -1);
            glVertex2f(ts*2-1+0.5f*lineWidth, -1);
   
            glVertex2f(ts*2-1-0.5f*lineWidth, -0.9); 
            glVertex2f(ts*2-1+0.5f*lineWidth, -0.9); 
            glVertex2f(ts*2-1+0.5f*lineWidth, -1);
        }

    }glEnd();
    
    glLineWidth(4);
    glBegin(GL_LINE_LOOP);{
        
        glColor3f(0.4f, 0.4f, 0.4f);
        glVertex2f(-1, -0.9);
        glVertex2f(+1, -0.9);
        glVertex2f(+1, -1);
        glVertex2f(-1, -1);
    }glEnd();
    
	//draw slider
    glBegin(GL_TRIANGLES);{
        glColor3f(207/255.0f,255/255.0f,177/255.0f);
        glVertex2f(timeStep_.get() / 624.0f * 2-1 - 0.5f * 2.0f * lineWidth, -0.85);
        glVertex2f(timeStep_.get() / 624.0f * 2-1 - 0.5f * 2.0f * lineWidth, -1.05);
        glVertex2f(timeStep_.get() / 624.0f * 2-1 + 0.5f * 2.0f * lineWidth, -1.05);
   
        glVertex2f(timeStep_.get() / 624.0f * 2-1 - 0.5f * 2.0f * lineWidth, -0.85);
        glVertex2f(timeStep_.get() / 624.0f * 2-1 + 0.5f * 2.0f * lineWidth, -0.85);
        glVertex2f(timeStep_.get() / 624.0f * 2-1 + 0.5f * 2.0f * lineWidth, -1.05);
    }glEnd();
    
    
    MatStack.popMatrix();
    glLineWidth(1);

    tgt::Font* f= font_.get();
    tgt::vec2 screensize = imageInport_.getSize();
    f->setTextAlignment(tgt::Font::TextAlignment::BottomLeft);
    //f->setFontSize(15);
    f->setFontColor(tgt::vec4(94/255.0f,180/255.0f,214/255.0f, 1.0f));
	//labels below slider
    for(int i = 0; i <= 6; i++){
		//float effectiveFontSize = f->getEffektiveFontSize(15);
			

		//don't show particular label if slider is near
		float distance = (std::abs((i * 100.0f) - timeStep_.get()) / 625.0f) * screensize.x * 0.9f;

		if (distance < f->getFontSize() * 2.0f)
            continue;

		//set lavbel position
        float pos = (i * 100.f) / 625.0f * 2.0f - 1.0f;
		//set label text
        std::string text = std::to_string(i*100);

        tgt::vec2 s =  f->getSize(translate(m, tgt::vec2(pos, -1.0f))*tgt::vec3(screensize, 1.0f), text, screensize);
        f->render(translate(m, tgt::vec2(pos, -1.0f))*tgt::vec3(screensize, 1.0f)-tgt::vec3(0.5f*s.x, 0, 0), text, screensize);
    }

	//slider label
    {
        float pos = timeStep_.get() / 625.0f * 2.0f - 1.0f;
        std::stringstream ss;
        f->setFontColor(tgt::vec4(1.0f, 0.0f, 0.0f, 1.0));
        ss.precision(0);
        //ss << std::fixed;
        ss << std::round(timeStep_.get());
      
        std::string text = ss.str();
        tgt::vec2 s =  f->getSize(translate(m, tgt::vec2(pos, -1.0f))*tgt::vec3(screensize, 1.0f), text, screensize);
        f->render(translate(m, tgt::vec2(pos, -1.0f))*tgt::vec3(screensize, 1.0f)-tgt::vec3(0.5f*s.x, 0, 0), text, screensize);
    }

    glEnable(GL_DEPTH_TEST);
    IMode.color(tgt::vec4::one);

    outport_.deactivateTarget();
    TextureUnit::setZeroUnit();
    LGL_ERROR;
}

} // namespace voreen
