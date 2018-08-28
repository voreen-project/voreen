/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2018 University of Muenster, Germany,                        *
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

#include "lineprofile.h"

#include "voreen/core/voreenapplication.h"

#include "tgt/textureunit.h"

#include <sstream>

using tgt::TextureUnit;

namespace voreen {


//****************************************************************************************************************************
//            Invaidate
//****************************************************************************************************************************
void LineProfile::invalidate(int inv) {
    //test, if imports have changed
    if(calculateNew_ == CN_NON){
        if(volSingleInport_.hasChanged() || volMultiInport_.hasChanged())
                calculateNew_ = CN_PLOTDATA;
        else
            if(imgInport_.hasChanged())
                calculateNew_ = CN_IMAGE;
            else
                return;
    }
    //set to image, if no line is drawn
    if(!hasLine_) calculateNew_ = CN_IMAGE;
    //must be called
    PropertyOwner::invalidate(inv);
    //used, if plotdata change
    PlotData* outData = 0;
    std::vector<const VolumeBase* > vvh;
    //begin calculation
    switch(calculateNew_) {
        case CN_PLOTDATA: //set pVolData
            //get volume data and return, if it is empty
            if(!inportVolume(vvh)) return;
            //calculate plotdata
            setPlotData(vvh, *pVolData_);
            setMaxValue();
        case CN_FITFUNCTION: //setFitfunction
            delete pModData_;
            pModData_ = new PlotData(*pVolData_);
            setFitFunction();
        case CN_THRESHOLD:
            setThresholdValue();
            setThresholdLine();
            outData = new PlotData(*pModData_);
            plotOutport_.setData(outData);
            plotOutport_.invalidatePort();
        case CN_TEXTOUTPUT:
            setTextOutport();
        case CN_IMAGE:
            calculateNew_ = CN_NON;
            imgOutport_.invalidatePort();
        case CN_NON:
            break;
        default:
            break;
    }

}

//****************************************************************************************************************************
//        utils
//****************************************************************************************************************************
bool LineProfile::inportVolume(std::vector<const VolumeBase* >& vvh){
    //return, if inport is not ready
    if(!volSingleInport_.isReady()) return false;
    //push_back first volume (volume to work with)
    vvh.push_back(volSingleInport_.getData());
    //push_back other volumes, if connected
    if(volMultiInport_.isReady()){
        std::vector<const VolumeBase* > vvh2;
        vvh2 = volMultiInport_.getAllData();
        for(unsigned int i = 0; i < vvh2.size() ; i++){
            //test, if volume is already in single port
            if(vvh.at(0) != vvh2.at(i))
                vvh.push_back(vvh2.at(i));
        }
    }
    return true;
}

bool LineProfile::setPlotData(std::vector<const VolumeBase* >& vvh, PlotData& pData) {
    if(interactionMode())
        return LP_Measure::calculate(staSpherePos_,endSpherePos_,1.0f/interactionQuality_.get(),
                            vvh, fetchOption_.getValue(),&pData);
    else
        return LP_Measure::calculate(staSpherePos_,endSpherePos_,samplingRate_.get(),
                            vvh, fetchOption_.getValue(),&pData);
}

void LineProfile::setFitFunction() {
    if(!calculateLM_.get()) return;
    if(!interactionMode()){
        if(fitOption_.getValue() == LP_EigenFunctionFit::FFM_GAUSSIAN_SUM) {
            FFS_GAUSSIAN_SUM gSum = FFS_GAUSSIAN_SUM(fitOptionAddOn_.get(),singleRepresentation_.get());
            //set function fit
            LP_EigenFunctionFit::setFitFunction(fitOption_.getValue(),&gSum);
            LP_EigenFunctionFit::functionFit(pVolData_,1,pModData_);
        } else if (fitOption_.getValue() == LP_EigenFunctionFit::FFM_GAUSSIAN_SINGLE){
            FFS_GAUSSIAN_SINGLE gSingle = FFS_GAUSSIAN_SINGLE(fitOptionAddOn_.get(),singleRepresentation_.get());
            //set function fit
            LP_EigenFunctionFit::setFitFunction(fitOption_.getValue(),&gSingle);
            LP_EigenFunctionFit::functionFit(pVolData_,1,pModData_);
        }
    }
    else {
        std::string label;
        switch(fitOption_.getValue()) {
        case LP_EigenFunctionFit::FFM_GAUSSIAN_SUM:
            label = "Gaussian Sum Fit of " + pVolData_->getColumnLabel(1);
            break;
        case LP_EigenFunctionFit::FFM_GAUSSIAN_SINGLE:
            label = "Gaussian Single Fit of " + pVolData_->getColumnLabel(1);
            break;
        }
        LP_PlotModifier::addFloatColumn(pModData_,0.0f,label);
    }
}

void LineProfile::setMaxValue() {
    LP_MaxThreshold::setMaxValue(maxValue_,sphereMaxPosition_,pVolData_,1);
    setLineAbsMax();
}

void LineProfile::setThresholdValue(){
    if(useAbsLine_.get())
        threshold_ = thresholdLineAbs_.get();
    else
        threshold_ = maxValue_*thresholdLinePer_.get()/100.0f;
}

void LineProfile::setThresholdLine(){
    if(pModData_->getColumnLabel(pModData_->getColumnCount()-1).compare(0,9,"Threshold") == 0)
        LP_PlotModifier::removeColumn(pModData_,pModData_->getColumnCount()-1);
    if(calculateThreshold_.get()){
        LP_PlotModifier::addFloatColumn(pModData_,threshold_,"Threshold Line of " + pVolData_->getColumnLabel(1));
        if(!interactionMode()){
            if(useLM_.get())
                if(singleRepresentation_.get())
                    LP_MaxThreshold::setThresholdArea(threshold_,thresholdArea_,pModData_,pModData_->getColumnCount()-2-fitOptionAddOn_.get());
                else
                    LP_MaxThreshold::setThresholdArea(threshold_,thresholdArea_,pModData_,pModData_->getColumnCount()-2);
            else
                LP_MaxThreshold::setThresholdArea(threshold_,thresholdArea_,pVolData_,1);
        }
    }
}

void LineProfile::setTextOutport(){
    //make strstr ready
    strstr_.str("");
    strstr_.clear();
    if(!interactionMode()){
        //set strstr
        if (showMax_.get())
            LP_TextEditor::setTextMaxValue(strstr_,maxValue_);
        if(calculateThreshold_.get()){
            LP_TextEditor::setTextThresholdValue(strstr_, threshold_);
            LP_TextEditor::setTextThresholdArea(strstr_,thresholdArea_);
        }
    }
    //make outport ready
    textOutport_.setData(strstr_.str());
}

//****************************************************************************************************************************
//            Process and Renering
//****************************************************************************************************************************
void LineProfile::process() {
    // set modelview and projection matrices
    MatStack.matrixMode(tgt::MatrixStack::PROJECTION);
    MatStack.loadMatrix(camera_.get().getProjectionMatrix(imgOutport_.getSize()));
    MatStack.matrixMode(tgt::MatrixStack::MODELVIEW);
    MatStack.loadMatrix(camera_.get().getViewMatrix());
    LGL_ERROR;

    //Test, if line is in image plane
    tgt::ivec2 size = imgInport_.getSize();
    tgt::vec3 A = fhpInport_.getRenderTarget()->getColorAtPos(tgt::ivec2(size.x/2,size.y/2)).xyz();
    tgt::vec3 B = fhpInport_.getRenderTarget()->getColorAtPos(tgt::ivec2(size.x/4,size.y/2)).xyz();
    tgt::vec3 C = fhpInport_.getRenderTarget()->getColorAtPos(tgt::ivec2(size.x/2,size.y/4)).xyz();
    tgt::vec3 N = tgt::cross(B-A,C-A);
    float res = tgt::dot(N,A);

//RENDER IN TEMP BUFFER THE ARROW
    //clear Buffer
    imgTempport_.activateTarget();
    glClearDepth(1.0);
    glClear(GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT);
    //Render
    if (arrowVisible_.get() && (tgt::dot(N,staSpherePos_) == res) && (tgt::dot(N,endSpherePos_) == res)) {
        drawLine();
    }
    imgTempport_.deactivateTarget();
    LGL_ERROR;

    // restore matrices
    MatStack.matrixMode(tgt::MatrixStack::PROJECTION);
    MatStack.loadIdentity();
    MatStack.matrixMode(tgt::MatrixStack::MODELVIEW);
    MatStack.loadIdentity();
    LGL_ERROR;

//COMBINE BOTH RESULTS
    // if inport is connected, combine both results:
   if (imgInport_.isReady()) {
        imgOutport_.activateTarget();
        imgOutport_.clearTarget();
        glClearDepth(1.0);
        glClearColor(0.0, 0.0, 0.0, 0.0);
        glClear(GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT);
        LGL_ERROR;

        imgInport_.bindTextures(GL_TEXTURE0, GL_TEXTURE1);
        imgTempport_.bindTextures(GL_TEXTURE2, GL_TEXTURE3);

        program_->activate();

        tgt::Camera cam = camera_.get();
        setGlobalShaderParameters(program_, &cam);
        program_->setUniform("colorTex0_", 0);
        program_->setUniform("depthTex0_", 1);
        imgInport_.setTextureParameters(program_, "textureParameters0_");

        program_->setUniform("colorTex1_", 2);
        program_->setUniform("depthTex1_", 3);
        imgTempport_.setTextureParameters(program_, "textureParameters1_");

        renderQuad();
        program_->deactivate();
        imgOutport_.deactivateTarget();
        LGL_ERROR;
    }
    TextureUnit::setZeroUnit();
    LGL_ERROR;

}

void LineProfile::drawLine() {
    float M[16];
    //get help vectors
    tgt::vec3 orien = endSpherePos_-staSpherePos_;
    tgt::vec3 Norien = tgt::normalize(orien);
    tgt::vec3 z(0,0,1);
    tgt::vec3 v = tgt::cross(z,Norien);
    //calculate rotation matrix
    if(tgt::dot(v,v) == 0) {
        if (orien.z > 0) {
            M[0] = 1; M[1] = 0; M[2] = 0; M[3] = 0;
            M[4] = 0; M[5] = 1; M[6] = 0; M[7] = 0;
            M[8] = 0; M[9] = 0; M[10] = 1; M[11] = 0;
            M[12] = 0; M[13] =  0; M[14] = 0; M[15]= 1;
        }
        else {
            M[0] = -1; M[1] = 0; M[2] = 0; M[3] = 0;
            M[4] = 0; M[5] = 1; M[6] = 0; M[7] = 0;
            M[8] = 0; M[9] = 0; M[10] = -1; M[11] = 0;
            M[12] = 0; M[13] =  0; M[14] = 0; M[15]= 1;
        }
    }
    else {
        float c = tgt::dot(z,Norien);
        float h = (1-c)/tgt::dot(v,v);
        M[0] = c+h*v.x*v.x; M[1] = h*v.x*-v.y+v.z; M[2] = h*v.x*v.z+v.y; M[3] = 0;
        M[4] = h*v.x*-v.y-v.z; M[5] = c+h*v.y*v.y; M[6] = h*-v.y*v.z+v.x; M[7] = 0;
        M[8] = h*v.x*v.z+-v.y; M[9] = h*-v.y*v.z-v.x; M[10] = c+h*v.z*v.z; M[11] = 0;
        M[12] = 0; M[13] =  0; M[14] = 0; M[15]= 1;
    }
    //draw differnt parts of the arrow
    MatStack.pushMatrix();
        MatStack.translate(staSpherePos_.x, staSpherePos_.y, staSpherePos_.z);
        MatStack.multMatrix(tgt::mat4(M));

        //draw marked arrow
        for(int i = 0; i < selectionProp_.size(); i++)
            LP_Graphic::drawLineInterval((float)selectionProp_.getSelectionAt(i).selection_.getSelection().at(0).second->getIntervalRepresentation().getLeft(),
            (float)selectionProp_.getSelectionAt(i).selection_.getSelection().at(0).second->getIntervalRepresentation().getRight(),tgt::length(orien),arrowSize_.get(),arrowMarkColor_.get());
        //draw threshold arrow
        if(showIntervals_.get() && !interactionMode())
            for(std::vector<std::pair<float,float> >::iterator it = thresholdArea_.begin(); it != thresholdArea_.end() ; it++)
                LP_Graphic::drawLineInterval(it->first,it->second,tgt::length(orien),arrowSize_.get(),arrowHotColor_.get());
        //draw normal arrow
        LP_Graphic::drawLineInterval(0,tgt::length(orien),tgt::length(orien),arrowSize_.get(),arrowColor_.get());
    MatStack.popMatrix();

    //draw spheres
    if(sphereMoveVisible_.get() || drawStaSphere_) {
        LP_Graphic::drawSphere(staSpherePos_,sphereMoveSize_.get(),sphereMoveColor_.get());
    }
    if(sphereMoveVisible_.get() || drawMidSphere_) {
        LP_Graphic::drawSphere(midSpherePos_,sphereMoveSize_.get(),sphereMoveColor_.get());
    }
    if(sphereMoveVisible_.get() || drawEndSphere_) {
        LP_Graphic::drawSphere(endSpherePos_,sphereMoveSize_.get(),sphereMoveColor_.get());
    }
    if(showMax_.get() && sphereMaxPosition_ >= 0 && !interactionMode()){
        LP_Graphic::drawSphere((staSpherePos_+tgt::normalize(endSpherePos_-staSpherePos_)*sphereMaxPosition_),sphereMaxSize_.get(),sphereMaxColor_.get());
    }
}


//****************************************************************************************************************************
//            General Functions
//****************************************************************************************************************************
LineProfile::LineProfile()
    : ImageProcessor("image/compositor")
    , hasLine_(false)
    //the ports
    , volSingleInport_(Port::INPORT,"volume.single.inport","volume.single.inport",false)
    , volMultiInport_(Port::INPORT, "volume.multi.inport", "volume.multi.inport", true)
    , imgInport_(Port::INPORT, "image.inport")
    , fhpInport_(Port::INPORT, "fhp.inport", "fhp.inport", false, Processor::INVALID_PROGRAM, RenderPort::RENDERSIZE_DEFAULT, GL_RGBA16F_ARB)
    , imgOutport_(Port::OUTPORT, "image.output")
    , imgTempport_(Port::OUTPORT, "image.tempport")
    , plotOutport_(Port::OUTPORT,"plotdata.outport")
    , textOutport_(Port::OUTPORT,"text.outport")
    //mouse events
    , mouseDownDrawLine_(false)
    , mouseDownMoveSphere_(false)
    //the camera
    , camera_("camera", "Camera", tgt::Camera(tgt::vec3(0.f, 0.f, 3.5f), tgt::vec3(0.f, 0.f, 0.f), tgt::vec3(0.f, 1.f, 0.f)))
    //properties: plotting
    , fetchOption_("fetchOption","Voxel Fetched By")
    , samplingRate_("samplingRate","Sampling Rate",2,1,1000)
    , interactionQuality_("interactionQuality","Interaction Quality",15,1,1000)
    //lm-algorithm
    , calculateLM_("calculateLM","Calculate Fit?",true)
    , useLM_("useLM","Use Fit For Threshold",true)
    , singleRepresentation_("singleRepresentation","Single Representation",true)
    , fitRepresentation_("fitRepresentation","Function Term","")
    , fitOption_("fitOption","Fit Function")
    , fitOptionAddOn_("foao","Number of Gaussian",2,1,5)
    //maxthreshold
    , calculateThreshold_("calculateThreshold","Use Threshold?",true)
    , showIntervals_("showIntervals","Show Intervals",true)
    , useAbsLine_("absLine","Absolut Value?",false)
    , thresholdLineAbs_("tla","Absolut Value",1.0f,0.0f,100000.0f)
    , thresholdLinePer_("tlp","Percent Value",50.0f,0.0f,100.0f)
    //properties: arrow
    , arrowVisible_("arrowVisible","Visible",true)
    , arrowSize_("arrowSize", "Size", 0.005f, 0.001f, 1.0f)
    , arrowColor_("arrowColor", "Arrow Color", tgt::vec4(1.0f, 1.0f, 1.0f, 1.0f))
    , arrowHotColor_("arrowHotColor", "Threshold Color", tgt::vec4(1.0f, 1.0f, 1.0f, 1.0f))
    , arrowMarkColor_("arrowMarkColor", "Marked Color", tgt::vec4(1.0f, 1.0f, 1.0f, 1.0f))
        //private poperties
        , selectionProp_("selectionProp","Selection",new PlotData(0,0),*(new PlotEntitiesProperty()),false)
    //properties: move sphere
    , sphereMoveVisible_("sphereMoveVisible","Allways Visible", false)
    , sphereMoveSize_("sphereMoveSize", "Size", 0.04f, 0.001f, 1.0f)
    , sphereMoveColor_("sphereMoveColor", "Color", tgt::vec4(1.0f, 1.0f, 1.0f, 1.0f))
    //properties: max sphere
    , showMax_("showMax","Show Maximum",true)
    , sphereMaxSize_("sphereMaxSize", "Size", 0.04f, 0.001f, 1.0f)
    , sphereMaxColor_("sphereMaxColor", "Color", tgt::vec4(1.0f, 1.0f, 1.0f, 1.0f))
    //other attributs
    , pVolData_(new PlotData(1,1)), pModData_(0), sphereMaxPosition_(-1.0f)
    , graped_(NON)
    , drawStaSphere_(false), drawMidSphere_(false), drawEndSphere_(false)
    , staSpherePos_(tgt::vec3(0.0f)), midSpherePos_(tgt::vec3(0.0f)), endSpherePos_(tgt::vec3(0.0f))
    , calculateNew_(CN_NON)

{
//adding ports to processor
    addPort(volSingleInport_);
    addPort(volMultiInport_);
    addPort(imgInport_);
    addPort(fhpInport_);
    addPort(imgOutport_);
    addPort(plotOutport_);
    addPort(textOutport_);
    addPrivateRenderPort(imgTempport_);

//adding the properties
    addProperty(camera_);
        camera_.setVisible(false);
    //plotting
    addProperty(fetchOption_);
        //set options
        fetchOption_.addOption("fo_normal","Nearest",LP_Measure::FO_NORMAL);
        fetchOption_.addOption("fo_linear","Linear",LP_Measure::FO_LINEAR);
        fetchOption_.addOption("fo_cubic","Cubic",LP_Measure::FO_CUBIC);
        fetchOption_.onChange(CallMemberAction<LineProfile>(this, &LineProfile::plotdataMakeNew));
    addProperty(samplingRate_);
        samplingRate_.onChange(CallMemberAction<LineProfile>(this, &LineProfile::plotdataMakeNew));
    addProperty(interactionQuality_);
        interactionQuality_.onChange(CallMemberAction<LineProfile>(this, &LineProfile::plotdataMakeNew));
    //fitting
    addProperty(calculateLM_);
        calculateLM_.onChange(CallMemberAction<LineProfile>(this, &LineProfile::calculateLMOnChange));
    addProperty(useLM_);
        useLM_.onChange(CallMemberAction<LineProfile>(this, &LineProfile::thresholdMakeNew));
    addProperty(singleRepresentation_);
        singleRepresentation_.onChange(CallMemberAction<LineProfile>(this, &LineProfile::fitFunctionMakeNew));
    addProperty(fitOption_);
        fitOption_.addOption("ff_gaussian_sum","Gauss Sum Function",LP_EigenFunctionFit::FFM_GAUSSIAN_SUM);
        fitOption_.addOption("ff_gaussian_single","Gauss Single Function",LP_EigenFunctionFit::FFM_GAUSSIAN_SINGLE);
        fitOption_.onChange(CallMemberAction<LineProfile>(this, &LineProfile::fitOptionOnChange));
    addProperty(fitOptionAddOn_);
        fitOptionAddOn_.onChange(CallMemberAction<LineProfile>(this, &LineProfile::fitOptionOnChange));
    addProperty(fitRepresentation_);
        fitRepresentation_.setWidgetsEnabled(false);
    //maxthreshold
    addProperty(calculateThreshold_);
        calculateThreshold_.onChange( CallMemberAction<LineProfile>(this, &LineProfile::calculateThresholdOnChange) );
    addProperty(showIntervals_);
        showIntervals_.onChange( CallMemberAction<LineProfile>(this, &LineProfile::imageMakeNew) );
    addProperty(useAbsLine_);
        useAbsLine_.onChange( CallMemberAction<LineProfile>(this, &LineProfile::useAbsLineOnChange) );
    addProperty(thresholdLineAbs_);
        thresholdLineAbs_.onChange( CallMemberAction<LineProfile>(this, &LineProfile::thresholdMakeNew) );
    addProperty(thresholdLinePer_);
        thresholdLinePer_.onChange( CallMemberAction<LineProfile>(this, &LineProfile::thresholdMakeNew) );
    //arrow
    addProperty(arrowVisible_);
        arrowVisible_.onChange( CallMemberAction<LineProfile>(this, &LineProfile::imageMakeNew) );
    addProperty(arrowSize_);
        arrowSize_.setNumDecimals(3);
        arrowSize_.onChange( CallMemberAction<LineProfile>(this, &LineProfile::imageMakeNew) );
    addProperty(arrowColor_);
        arrowColor_.setViews(Property::COLOR);
        arrowColor_.onChange( CallMemberAction<LineProfile>(this, &LineProfile::imageMakeNew) );
    addProperty(arrowHotColor_);
        arrowHotColor_.setViews(Property::COLOR);
        arrowHotColor_.onChange( CallMemberAction<LineProfile>(this, &LineProfile::imageMakeNew) );
    addProperty(arrowMarkColor_);
        arrowMarkColor_.setViews(Property::COLOR);
        arrowMarkColor_.onChange( CallMemberAction<LineProfile>(this, &LineProfile::imageMakeNew) );
    addProperty(selectionProp_);
        selectionProp_.setVisible(false);
        selectionProp_.onChange( CallMemberAction<LineProfile>(this, &LineProfile::imageMakeNew) );
    //move spheres
    addProperty(sphereMoveVisible_);
    sphereMoveVisible_.onChange( CallMemberAction<LineProfile>(this, &LineProfile::imageMakeNew) );
    addProperty(sphereMoveSize_);
        sphereMoveSize_.setNumDecimals(3);
        sphereMoveSize_.onChange( CallMemberAction<LineProfile>(this, &LineProfile::imageMakeNew) );
    addProperty(sphereMoveColor_);
        sphereMoveColor_.setViews(Property::COLOR);
        sphereMoveColor_.onChange( CallMemberAction<LineProfile>(this, &LineProfile::imageMakeNew) );
    //max sphere
    addProperty(showMax_);
        showMax_.onChange( CallMemberAction<LineProfile>(this, &LineProfile::textOutputMakeNew) );
    addProperty(sphereMaxSize_);
        sphereMaxSize_.setNumDecimals(3);
        sphereMaxSize_.onChange( CallMemberAction<LineProfile>(this, &LineProfile::imageMakeNew) );
    addProperty(sphereMaxColor_);
        sphereMaxColor_.setViews(Property::COLOR);
        sphereMaxColor_.onChange( CallMemberAction<LineProfile>(this, &LineProfile::imageMakeNew) );

//set GUI
    fetchOption_.setGroupID("Plotting");
    samplingRate_.setGroupID("Plotting");
    interactionQuality_.setGroupID("Plotting");
    calculateLM_.setGroupID("Function Fit");
    useLM_.setGroupID("Function Fit");
    singleRepresentation_.setGroupID("Function Fit");
    fitOption_.setGroupID("Function Fit");
    fitOptionAddOn_.setGroupID("Function Fit");
    fitRepresentation_.setGroupID("Function Fit");
    calculateThreshold_.setGroupID("Threshold Line");
    showIntervals_.setGroupID("Threshold Line");
    useAbsLine_.setGroupID("Threshold Line");
    thresholdLineAbs_.setGroupID("Threshold Line");
    thresholdLinePer_.setGroupID("Threshold Line");
    arrowVisible_.setGroupID("Arrow");
    arrowSize_.setGroupID("Arrow");
    arrowColor_.setGroupID("Arrow");
    arrowHotColor_.setGroupID("Arrow");
    arrowMarkColor_.setGroupID("Arrow");
    sphereMoveVisible_.setGroupID("Spheres");
    sphereMoveSize_.setGroupID("Spheres");
    sphereMoveColor_.setGroupID("Spheres");
    showMax_.setGroupID("Maximum Sphere");
    sphereMaxSize_.setGroupID("Maximum Sphere");
    sphereMaxColor_.setGroupID("Maximum Sphere");

    setPropertyGroupGuiName("Plotting","Plotting");
    setPropertyGroupGuiName("Function Fit","Function Fit");
    setPropertyGroupGuiName("Threshold Line","Threshold Line");
    setPropertyGroupGuiName("Arrow","Arrow");
    setPropertyGroupGuiName("Spheres","Spheres");
    setPropertyGroupGuiName("Maximum Sphere","Maximum Sphere");

//setting the mouse events
    mouseEventDrawLine_ = new EventProperty<LineProfile>(
        "mouseEvent.drawLine", "Draw Line", this,
        &LineProfile::drawLine,
        tgt::MouseEvent::MOUSE_BUTTON_LEFT,
        tgt::MouseEvent::MOTION | tgt::MouseEvent::PRESSED | tgt::MouseEvent::RELEASED,
        tgt::Event::CTRL, false);

    mouseEventMoveSphere_ = new EventProperty<LineProfile>(
        "mouseEvent.movesShere", "Move Sphere", this,
        &LineProfile::moveSphere,
         tgt::MouseEvent::MOUSE_BUTTON_NONE,
        tgt::MouseEvent::MOTION, tgt::Event::MODIFIER_NONE, false);

    mouseEventPressSphere_ = new EventProperty<LineProfile>(
        "mouseEvent.pressSphere", "Press Sphere",
        this, &LineProfile::moveSphere,
        tgt::MouseEvent::MOUSE_BUTTON_LEFT,
        tgt::MouseEvent::PRESSED | tgt::MouseEvent::RELEASED | tgt::MouseEvent::MOTION, tgt::Event::MODIFIER_NONE,
        false);

    addEventProperty(mouseEventDrawLine_);
    addEventProperty(mouseEventMoveSphere_);
    addEventProperty(mouseEventPressSphere_);

    //preset all properties
    preSettings();
}

LineProfile::~LineProfile() {
    //delete selectionProp_ attributes
    delete selectionProp_.getPlotData();
    delete &(selectionProp_.getEntitiesProperty());
    delete pVolData_;
    delete pModData_;
    //delete mouse events
    delete mouseEventDrawLine_;
    delete mouseEventMoveSphere_;
    delete mouseEventPressSphere_;
}

Processor* LineProfile::create() const {
    return new LineProfile();
}

bool LineProfile::isReady() const {
    if (!isInitialized() || !imgInport_.isReady() || !fhpInport_.isReady() || !imgOutport_.isReady() || !volSingleInport_.isReady() || !plotOutport_.isReady())
        return false;
    return true;
}

std::string LineProfile::generateHeader(const tgt::GpuCapabilities::GlVersion* /*version*/) {
    std::string header = RenderProcessor::generateHeader();
    header += "#define MODE_ALPHA_COMPOSITING\n";
    return header;
}

void LineProfile::interactionModeToggled(){
    if(!interactionMode()){
        //FIXME: update everthing?
        calculateNew_ = CN_PLOTDATA;
        invalidate();
    }
}

//****************************************************************************************************************************
//        On Change Functions
//****************************************************************************************************************************
void LineProfile::plotdataMakeNew() {
    if(isInitialized()){
        calculateNew_ = CN_PLOTDATA;
        invalidate();
    }
}

void LineProfile::fitFunctionMakeNew() {
    if(isInitialized()){
        calculateNew_ = CN_FITFUNCTION;
        invalidate();
    }
}

void LineProfile::thresholdMakeNew() {
    if(isInitialized()){
        calculateNew_ = CN_THRESHOLD;
        invalidate();
    }
}

void LineProfile::textOutputMakeNew() {
    if(isInitialized()){
        calculateNew_ = CN_TEXTOUTPUT;
        invalidate();
    }
}

void LineProfile::imageMakeNew() {
    if(isInitialized()){
        calculateNew_ = CN_IMAGE;
        invalidate();
    }
}

void LineProfile::preSettings(){
    fitOptionOnChange();
    calculateLMOnChange();
    calculateThresholdOnChange();
    useAbsLineOnChange();
}

void LineProfile::calculateLMOnChange() {
    if(calculateLM_.get()){
        useLM_.setWidgetsEnabled(true);
        singleRepresentation_.setWidgetsEnabled(true);
        fitOption_.setWidgetsEnabled(true);
        fitOptionAddOn_.setWidgetsEnabled(true);
    } else {
        useLM_.setWidgetsEnabled(false);
        useLM_.set(false);
        singleRepresentation_.setWidgetsEnabled(false);
        fitOption_.setWidgetsEnabled(false);
        fitOptionAddOn_.setWidgetsEnabled(false);
    }
    fitFunctionMakeNew();
}

void LineProfile::fitOptionOnChange() {
    switch(fitOption_.getValue()){
    case LP_EigenFunctionFit::FFM_GAUSSIAN_SUM:
        fitRepresentation_.set("Sum[i=0;n][a_i*e^-((x-m_i)/s_i)^2]+b_1");
        break;
    case LP_EigenFunctionFit::FFM_GAUSSIAN_SINGLE:
        fitRepresentation_.set("a_i*e^-((x-m_i)/s_i)^2+b_i");
        break;
    }
    fitFunctionMakeNew();
}

void LineProfile::calculateThresholdOnChange(){
    if(calculateThreshold_.get()){
        showIntervals_.setWidgetsEnabled(true);
        useAbsLine_.setWidgetsEnabled(true);
        thresholdLineAbs_.setWidgetsEnabled(true);
        thresholdLinePer_.setWidgetsEnabled(true);
    } else{
        showIntervals_.setWidgetsEnabled(false);
        showIntervals_.set(false);
        useAbsLine_.setWidgetsEnabled(false);
        thresholdLineAbs_.setWidgetsEnabled(false);
        thresholdLinePer_.setWidgetsEnabled(false);
    }
    thresholdMakeNew();
}

void LineProfile::useAbsLineOnChange(){
    if(useAbsLine_.get()){
        thresholdLineAbs_.setVisible(true);
        thresholdLinePer_.setVisible(false);
    } else{
        thresholdLineAbs_.setVisible(false);
        thresholdLinePer_.setVisible(true);
    }
    thresholdMakeNew();
}

void LineProfile::setLineAbsMax(){
    if(interactionMode()) return;
    if(thresholdLineAbs_.get() > ((int)maxValue_+1))
        thresholdLineAbs_.set(maxValue_);
    thresholdLineAbs_.setMaxValue((const float)((int)maxValue_+1));
}

//****************************************************************************************************************************
//        Mouse Events
//****************************************************************************************************************************
tgt::ivec2 LineProfile::cropToViewport(tgt::ivec2 mousePos) {
    tgt::ivec2 result = mousePos;
    tgt::ivec2 size = imgInport_.getSize();
    if (result.x < 0) result.x = 0;
    else if (result.x > size.x-1) result.x = size.x-1;
    if (result.y < 0) result.y = 0;
    else if (result.y > size.y-1) result.y = size.y-1;
    return result;
}

void LineProfile::drawLine(tgt::MouseEvent* e) {
    //mouse is pressed
    if (e->action() & tgt::MouseEvent::PRESSED && !mouseDownDrawLine_) {
        //crop mouse position to canvas size
        tgt:: ivec2 mousePos = cropToViewport(tgt::ivec2(e->coord().x, e->viewport().y-e->coord().y));
        //get world position
        tgt::vec3 Color = fhpInport_.getRenderTarget()->getColorAtPos(mousePos).xyz();
        //set position of spheres
        endSpherePos_ = Color; staSpherePos_ = Color; midSpherePos_ = Color;
        //mouse is pressed and arrow position has changed
        mouseDownDrawLine_ = true; calculateNew_ = CN_PLOTDATA;
        //calculation can start
        hasLine_ = true;
        //activate interaction mode
        toggleInteractionMode(true,this);
        e->accept();
    }
    //mouse is pressed and moves
    if (e->action() & tgt::MouseEvent::MOTION && mouseDownDrawLine_) {
        //crop mouse position to canvas size
        tgt::ivec2 mousePos = cropToViewport(tgt::ivec2(e->coord().x, e->viewport().y-e->coord().y));
        //get world position
        tgt::vec3 Color = fhpInport_.getRenderTarget()->getColorAtPos(mousePos).xyz();
        //set position of spheres
        endSpherePos_ = Color;
        midSpherePos_ = (endSpherePos_-staSpherePos_)*(tgt::vec3(0.5))+staSpherePos_;
        //arrow position has changed
        calculateNew_ = CN_PLOTDATA;
        e->accept();
    }
    //mouse is released
    if ((e->action() & tgt::MouseEvent::RELEASED) && (mouseDownDrawLine_ || mouseDownMoveSphere_)) {
        //mouse is no longer pressed
        mouseDownDrawLine_ = false;    mouseDownMoveSphere_ = false; graped_ = NON;
        //deactivate interaction mode and set hasChanged
        toggleInteractionMode(false,this);
        calculateNew_ = CN_PLOTDATA;
        e->accept();
    }
    //if arrow position has changed, invalidate
    if(calculateNew_ == CN_PLOTDATA) invalidate();
}

void LineProfile::moveSphere(tgt::MouseEvent* e) {
    //mouse moved
    if (e->action() & tgt::MouseEvent::MOTION ) {
        //crop mouse position to canvas size
        tgt::ivec2 mouseTempPos = cropToViewport(tgt::ivec2(e->coord().x, e->viewport().y-e->coord().y));
        //get world position
        tgt::vec3 Color = fhpInport_.getRenderTarget()->getColorAtPos(mouseTempPos).xyz();
        tgt::vec3 moveVec; //used in cas MID

        switch(graped_) {
        //if start sphere is graped
        case STA:
            e->accept();
            staSpherePos_ = Color;
            midSpherePos_ = (endSpherePos_-staSpherePos_)*(tgt::vec3(0.5))+staSpherePos_;
            calculateNew_ = CN_PLOTDATA;
            break;
        //if end sphere is graped
        case END:
            e->accept();
            endSpherePos_ = Color;
            midSpherePos_ = (endSpherePos_-staSpherePos_)*(tgt::vec3(0.5))+staSpherePos_;
            calculateNew_ = CN_PLOTDATA;
            break;
        //if middle sphere is graped
        case MID:
            e->accept();
            moveVec = Color-midSpherePos_;
            staSpherePos_ = staSpherePos_+moveVec;
            endSpherePos_ = endSpherePos_+moveVec;
            midSpherePos_ = (endSpherePos_-staSpherePos_)*(tgt::vec3(0.5))+staSpherePos_;
            calculateNew_ = CN_PLOTDATA;
            break;
        //if no sphere is graped, test, if the mouse is over one sphere
        case NON:
            //start sphere
            if(tgt::distance(Color,staSpherePos_) <= sphereMoveSize_.get()) {
                if(!drawStaSphere_) {
                    drawStaSphere_ = true;
                    calculateNew_ = CN_IMAGE;
                }
            }
            else {
                if(drawStaSphere_) {
                    drawStaSphere_ = false;
                    calculateNew_ = CN_IMAGE;
                }
            }
            //end sphere
            if(tgt::distance(Color,endSpherePos_) <= sphereMoveSize_.get()) {
                if(!drawEndSphere_) {
                    drawEndSphere_ = true;
                    calculateNew_ = CN_IMAGE;
                }
            }
            else {
                if(drawEndSphere_) {
                    drawEndSphere_ = false;
                    calculateNew_ = CN_IMAGE;
                }
            }
            //mid sphere
            if(tgt::distance(Color,midSpherePos_) <= sphereMoveSize_.get()) {
                if(!drawMidSphere_) {
                    drawMidSphere_ = true;
                    calculateNew_ = CN_IMAGE;
                }
            }
            else {
                if(drawMidSphere_) {
                    drawMidSphere_ = false;
                    calculateNew_ = CN_IMAGE;
                }
            }
            break;
        }
        //FIXME: problem should be solved in 2D camera
        if(mouseDownDrawLine_){
            e->accept();
        }
    }
    //mouse pressed
    if (e->action() & tgt::MouseEvent::PRESSED && !mouseDownMoveSphere_) {
        //crop mouse position to canvas size
        tgt::ivec2 mouseTempPos_ = cropToViewport(tgt::ivec2(e->coord().x, e->viewport().y-e->coord().y));
        //get world position
        tgt::vec3 Color = fhpInport_.getRenderTarget()->getColorAtPos(mouseTempPos_).xyz();
        //test, if a sphere was clicked
        if (tgt::distance(Color,staSpherePos_) <= sphereMoveSize_.get()) {
            graped_ = STA;
            mouseDownMoveSphere_ = true;
            e->accept();
        }
        else if (tgt::distance(Color,midSpherePos_) <= sphereMoveSize_.get()) {
            graped_ = MID;
            mouseDownMoveSphere_ = true;
            e->accept();
        }
        else if (tgt::distance(Color,endSpherePos_) <= sphereMoveSize_.get()) {
            graped_ = END;
            mouseDownMoveSphere_ = true;
            e->accept();
        }
        //if a sphere is graped, switch to interaction mode
        if(graped_ != NON) {
            toggleInteractionMode(true,this);
            calculateNew_ = CN_PLOTDATA;
        }
    }

    //mouse is released
    if (e->action() & tgt::MouseEvent::RELEASED && (mouseDownDrawLine_ || mouseDownMoveSphere_)) {
        //mouse is no longer pressed
        mouseDownDrawLine_ = false; mouseDownMoveSphere_ = false; graped_ = NON;
        //deactivate interaction mode
        toggleInteractionMode(false,this);
        calculateNew_ = CN_PLOTDATA;
        e->accept();
    }
    // if image/plotdata have to be updated
    if (calculateNew_ != CN_NON) invalidate();
}

} // namespace voreen
