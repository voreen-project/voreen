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

#include "vtkVoreenActor.h"
#include "vtkdummycanvas.h"

#include "vtkObjectFactory.h"
#include "vtkRenderWindow.h"
#include "vtkRendererCollection.h"

#include "tgt/shadermanager.h"
#include "tgt/camera.h"

#include "voreen/core/utils/voreenpainter.h"
#include "voreen/core/network/networkevaluator.h"
#include "voreen/core/network/processornetwork.h"
#include "voreen/core/network/workspace.h"
#include "voreen/core/processors/canvasrenderer.h"
#include "voreen/core/io/serialization/serialization.h"
#include "voreen/modules/base/basemodule.h"
#include "voreen/modules/base/processors/proxygeometry/cubemeshproxygeometry.h"
#include "voreen/modules/base/processors/geometry/geometryprocessor.h"
#include "voreen/core/datastructures/volume/volumecontainer.h"
#include "voreen/core/datastructures/volume/volumeatomic.h"
#include "voreen/core/datastructures/volume/volume.h"
#include "voreen/modules/base/processors/render/singlevolumeraycaster.h"
#include "voreen/modules/base/processors/render/simpleraycaster.h"
#include "voreen/modules/base/processors/render/glslraycaster.h"
#include "voreen/modules/base/processors/geometry/boundingboxrenderer.h"
#include "voreen/modules/base/processors/datasource/volumesource.h"
#include "voreen/modules/base/processors/image/grayscale.h"

#include "voreen/core/properties/cameraproperty.h"
#include "voreen/core/properties/stringproperty.h"
#include "voreen/core/datastructures/transfunc/transfuncmappingkey.h"

#include "voreen/core/voreenapplication.h"

// Convert define into a string, compare
// http://gcc.gnu.org/onlinedocs/cpp/Stringification.html
#define VRN_XSTRINGIFY(s) VRN_STRINGIFY(s)
#define VRN_STRINGIFY(s) #s

using namespace voreen;

vtkCxxRevisionMacro(vtkVoreenActor, "$Revision: 1.2 $");
vtkStandardNewMacro(vtkVoreenActor);
NetworkEvaluator* networkEvaluator_ = 0;
bool networkChanged_ = false;
bool singleVolumeRaycaster_ = true;
bool simpleRaycaster_ = false;
bool glslRaycaster_ = false;
std::string raycaster_ = "SingleVolumeRaycaster";
std::string geometry_ = "CubeMeshProxyGeometry";
std::string compositing1_ = "dvr";
std::string selectedPlane_ = "left";
int selectedKey_ = 0;
int shading_ = 5;
voreen::VtkDummyCanvas* canvas_ = 0;
tgt::Camera* cam_ = 0;
vtkRenderer *renderer_ = 0;
Workspace* workspace_ = 0;
std::stringstream outstream_;
vtkTextActor* statusText_ = 0;
vtkRenderer* textRenderer_ = 0;
vtkImageReader2* inputReader_ = 0;

vtkVoreenActor::vtkVoreenActor() {
    voreen_initialized_ = false;
}

vtkVoreenActor::~vtkVoreenActor() {
}

void vtkVoreenActor::printSelf(ostream& os, vtkIndent indent) {
    cout << "Status of vtkVoreenActor:\n";
    this->Superclass::PrintSelf(os,indent);
    cout << "==================\n";
}

void vtkVoreenActor::voreen_init() {
    voreen_initialized_ = true;

    VoreenApplication* app = new VoreenApplication("vtkVoreen", "vtkVoreen", 0, 0);
    
    // add base module
    app->addModule(new BaseModule());

    //window control is done by vtk, but the canvasRenderer needs a canvas with its painter to work 
    canvas_ = new VtkDummyCanvas("VtkDummyCanvas", tgt::ivec2(1,1), tgt::GLCanvas::RGBADD);
    canvas_->init();
    
    //init voreen application
    app->init();
    app->initGL();

    //set shader paths
    ShdrMgr.addPath(app->getShaderPath());
    ShdrMgr.addPath(app->getShaderPath() + "/util");

    std::string voreenDir = VRN_XSTRINGIFY(VRN_DIR);
    ShdrMgr.addPath(voreenDir + "/src/core/glsl");
    ShdrMgr.addPath(voreenDir + "/src/modules/base/glsl");
    ShdrMgr.addPath("../../src/core/glsl");
    ShdrMgr.addPath("../../src/modules/base/glsl");
    
    //create new workspace and load standard workspace
    workspace_ = new Workspace();
    networkEvaluator_ = new NetworkEvaluator();
    PrintMessage("[vtkVoreenActor] : Successfully initialized");
}


void vtkVoreenActor::loadWorkspace(std::string workspace) {

    if (!isInitialized()) {
        PrintMessage("[vtkVoreenActor]: Not initialized.");
        return;
    }
    tgtAssert(workspace_, "No workspace");

    workspace_->load(workspace);
    int useCanvas = 0;  //currently we only use the first canvas
    //get the canvas renderer from workspace processors
    ProcessorNetwork* network = workspace_->getProcessorNetwork();
    if (!network) {
        PrintMessage("[vtkVoreenActor]: Workspace loading failed: " + workspace);
        networkEvaluator_->setProcessorNetwork(0);
        return;
    }

    std::vector<CanvasRenderer*> canvasRenderer = network->getProcessorsByType<CanvasRenderer>();
    if(canvasRenderer.size() > useCanvas) {
        //create painter and connect painter, network evaluator and canvas renderer
        VoreenPainter* painter = new VoreenPainter(canvas_, networkEvaluator_, canvasRenderer[useCanvas]);
        canvas_->setPainter(painter);
        canvasRenderer[useCanvas]->setCanvas(canvas_);
        //add network to evaluator
        networkEvaluator_->lock();
        networkChanged_ = true;
        networkEvaluator_->setProcessorNetwork(network);
    }
    PrintMessage("[vtkVoreenActor]: Workspace loaded: \n" + workspace);
}

void vtkVoreenActor::connectToRenderer(vtkRenderer* ren) {
    renderer_ = ren;
    createStatusOutput();
    voreen_init();
    renderer_->AddVolume(this);
}

void vtkVoreenActor::createStatusOutput() {

    vtkRenderer *textRen = vtkRenderer::New();
    textRen->SetViewport(0.7,0,1,1);
    textRen->SetBackground(0,0,0);
    renderer_->GetRenderWindow()->AddRenderer(textRen);
    vtkTextActor* status = vtkTextActor::New();
    status->SetInput("Status");
    vtkTextProperty* statusprop = status->GetTextProperty();
    statusprop->SetLineSpacing(1.0);
    statusprop->SetShadow(0);
    statusprop->SetFontSize(14);
    statusprop->SetBold(0);
    statusprop->SetShadow(0);
    statusprop->SetShadowOffset(0.1,0.1);
    statusprop->SetFontFamilyToCourier();
    statusprop->SetColor(1,1,1);
    status->SetDisplayPosition(textRen->GetViewport()[0]*renderer_->GetRenderWindow()->GetSize()[0]+10,10);
    renderer_->SetViewport(0.0,0.0,0.7,1.0);
    textRen->AddActor2D(status);
    textRen->Render();
    setStatusOutput(textRen,status);
    PrintMessage("Press h for command list  \nPress s to toggle status window \n========================");
}

void vtkVoreenActor::setStatusOutput(vtkRenderer* textRen, vtkTextActor* status) {
    statusText_ = status;
    textRenderer_ = textRen;
}

void vtkVoreenActor::syncCamera(vtkCamera* cam){

    if (!isInitialized()) {
        PrintMessage("[vtkVoreenActor]: Not initialized.");
        return;
    }
    tgtAssert(workspace_, "No workspace");

    //Get values from tgt camera
    double pos[] = {0,0,0};
    cam->GetPosition(pos);
    double up[] = {0,0,0};
    cam->GetViewUp(up);
    double foc[] = {0,0,0};
    cam->GetFocalPoint(foc);
    double angle = cam->GetViewAngle();
    double clip[]= {0,0};
    cam->GetClippingRange(clip);
    
    //Get all cameras from workspace and set them to the vtk values
    std::vector<CameraProperty*> camprops = workspace_->getProcessorNetwork()->getPropertiesByType<CameraProperty>();

    for (int i = 0; i < camprops.size();i++){
        tgt::Camera* tgt_cam = tgt::Camera(tgt::vec3(pos[0],pos[1],pos[2]),tgt::vec3(foc[0],foc[1],foc[2]), tgt::vec3(up[0],up[1],up[2]),angle,1,clip[0],clip[1]);
        camprops.at(i)->set(*tgt_cam);
    }
}

VtkDummyCanvas* vtkVoreenActor::getCanvas() {
    return canvas_;
}

std::string vtkVoreenActor::getWorkspaceName(){
    if (workspace_)
        return workspace_->getFilename();
    else
        return "";
}

bool vtkVoreenActor::isInitialized() {
    return voreen_initialized_;
}

void vtkVoreenActor::listProcessors() {

    tgtAssert(workspace_, "No workspace");
    if (!workspace_->getProcessorNetwork()) {
        PrintMessage("[vtkVoreenActor]: No processor network");
        return;
    }

    std::vector<Processor*> processorlist = workspace_->getProcessorNetwork()->getProcessors();
    outstream_ << "======================== \n";
    outstream_ << "Processors in Workspace: \n";
    for (size_t i=0; i < processorlist.size(); i++) {
        outstream_ << "Prozessor ";
        outstream_ << i;
        outstream_ << ": ";
        outstream_ << processorlist[i]->getName();
        outstream_ << " type: ";
        outstream_ << processorlist[i]->getCategory();
        outstream_ << "\n";
    }
    
    StatusOutput();
}

void vtkVoreenActor::listProperties(std::string procName) {

    tgtAssert(workspace_, "No network evaluator");
    if (!workspace_->getProcessorNetwork()) {
        PrintMessage("[vtkVoreenActor]: No processor network");
        return;
    }

    Processor* processor = 0;
    processor = workspace_->getProcessorNetwork()->getProcessor(procName);
    
    if (!processor) {
        outstream_ << "[vtkVoreenActor]: No processor named " << procName << " found!" << endl;
        return;
    }
    outstream_ << "======================== \n";
    outstream_ << "Selected Processor:"  << processor->getName() << endl;
    std::vector<Property*> props = processor->getProperties();
    for (size_t i=0; i < props.size(); i++) {
        FloatProperty* floatProp = 0;
        IntProperty* intProp = 0;
        BoolProperty* boolProp = 0;
        StringProperty* stringProp = 0;
        StringOptionProperty* stringOptionProp = 0;
        floatProp = dynamic_cast<FloatProperty*>(props.at(i));      
        intProp = dynamic_cast<IntProperty*>(props.at(i));
        boolProp = dynamic_cast<BoolProperty*>(props.at(i));
        stringProp = dynamic_cast<StringProperty*>(props.at(i));
        stringOptionProp = dynamic_cast<StringOptionProperty*>(props.at(i));
        outstream_ << "Property ID ";
        outstream_ << i;
        outstream_ << ": ";
        outstream_ << props.at(i)->getID();
        outstream_ << " Type&Value: ";
        if (floatProp != 0) {
            outstream_ << "FloatProperty (" << floatProp->get() << ") " ;
            outstream_ << "(" << floatProp->getMinValue() << " " << floatProp->getStepping() << " " << floatProp->getMaxValue() << ")";
        }
        else if (intProp != 0) {
            outstream_ << "IntProperty (" << intProp->get() << ") " ;
            outstream_ << "(" << intProp->getMinValue() << " " << intProp->getStepping() << " " << intProp->getMaxValue() << ")";
        }
        else if (boolProp != 0) {
            outstream_ << "BoolProperty (" << boolProp->get() << ") " ;
        }
        else if (stringProp != 0) {
            outstream_ << "StringProperty (" << stringProp->get() << ") " ;
        }
       else if (stringOptionProp != 0) {
            outstream_ << "StringOptionProperty (" << stringOptionProp->get() << ") " ;
        }
        else 
            outstream_ << "undefined Property-type."; 
        outstream_ << "\n";
    }

    StatusOutput();
}

void vtkVoreenActor::listProperties(int procID) {

    tgtAssert(workspace_, "No workspace");
    if (!workspace_->getProcessorNetwork()) {
        PrintMessage("[vtkVoreenActor]: No processor network");
        return;
    }

    Processor* processor = 0;
    std::vector<Processor*> processors = workspace_->getProcessorNetwork()->getProcessors();
   
    if (procID < processors.size()) 
        listProperties(processors[procID]->getName());
    else {
        outstream_ << "[vtkVoreenActor]: No processor with ID " << procID << " in network." << endl;
        StatusOutput();
    }  
}

void vtkVoreenActor::StatusOutput() {

    if (statusText_ && textRenderer_){
        statusText_->SetInput(outstream_.str().c_str());
        textRenderer_->GetRenderWindow()->Render();
    } 
    else {
        cout << "[vtkVoreenActor] : Not connected to a status window!";
    }
}

void vtkVoreenActor::PrintMessage(std::string s) {

    if (isInitialized()) {
        LINFOC("vtkVoreenActor", s);
    }

    if (outstream_.str().size() > 4000) {
        outstream_.str("");
        outstream_.clear();
        outstream_ << "\n=========Outstream resettet=========\n";
    }
    outstream_ << "======================== \n";
    outstream_ << s << "\n";

    StatusOutput();
}


void vtkVoreenActor::changeIntPropByName(std::string procName, std::string propName, bool increase) {
    
    tgtAssert(workspace_, "No workspace");
    if (!workspace_->getProcessorNetwork()) {
        PrintMessage("[vtkVoreenActor]: No processor network");
        return;
    }

    Processor* proxy = 0;
    proxy = workspace_->getProcessorNetwork()->getProcessor(procName);

    if (!proxy) {
        PrintMessage("[vtkVoreenActor]: No processor named " + procName + " found!");
        return;
    }

    std::vector<IntProperty*> intVec = proxy->getPropertiesByType<IntProperty>();

    if (!increase) {
        for (size_t i=0; i<intVec.size(); i++ ) {
            if (intVec.at(i)->getID() == propName) {
                intVec.at(i)->decrease();
                renderer_->GetVTKWindow()->Render();
                return;
            }    
        }
        PrintMessage("[vtkVoreenActor]: No property named " + propName + " found!");
    }
    else {
        for (size_t i=0; i<intVec.size(); i++) {
            if (intVec.at(i)->getID() == propName) {
                intVec.at(i)->increase();
                renderer_->GetVTKWindow()->Render();
                return;
            }
        }
        PrintMessage("[vtkVoreenActor]: No property named " + propName + " found!");
    }
}


void vtkVoreenActor::changeIntPropByName(std::string procName, std::string propName, int value, int index) {
    
    tgtAssert(workspace_, "No workspace");
    if (!workspace_->getProcessorNetwork()) {
        PrintMessage("[vtkVoreenActor]: No processor network");
        return;
    }
    
    Processor* proxy = 0;
    proxy = workspace_->getProcessorNetwork()->getProcessor(procName);

    if (!proxy) {
        PrintMessage("[vtkVoreenActor]: No processor named " + procName + " found!");
        return;
    }

    std::vector<IntProperty*> intVec = proxy->getPropertiesByType<IntProperty>();

    for (size_t i=0; i<intVec.size(); i++ ) {
        if (intVec.at(i)->getID() == propName) {
            if (index == 0) 
                intVec.at(i)->set(value);
            if (index == 1)
                intVec.at(i)->setStepping(value);
            if (index == 2)
                intVec.at(i)->setMinValue(value);
            if (index == 3)
                intVec.at(i)->setMaxValue(value);
            renderer_->GetVTKWindow()->Render();
            return;
        }    
   }
   PrintMessage("[vtkVoreenActor]: No property named " + propName + " found!");  
}


void vtkVoreenActor::changeFloatPropByName(std::string procName, std::string propName, bool increase) {
    
    tgtAssert(workspace_, "No workspace");
    if (!workspace_->getProcessorNetwork()) {
        PrintMessage("[vtkVoreenActor]: No processor network");
        return;
    }
    
    Processor* proxy = 0;
    proxy = workspace_->getProcessorNetwork()->getProcessor(procName);

    if (!proxy) {
        PrintMessage("[vtkVoreenActor]: No processor named " + procName + " found!");
        return;
    }

    std::vector<FloatProperty*> intVec = proxy->getPropertiesByType<FloatProperty>();

    if (!increase) {
        for (size_t i=0; i<intVec.size(); i++ ) {
            if (intVec.at(i)->getID() == propName) {
                intVec.at(i)->decrease();
                renderer_->GetVTKWindow()->Render();
                return;
           }    
        }
        PrintMessage("[vtkVoreenActor]: No property named " + propName + " found!");
    }
    else {
        for (size_t i=0; i<intVec.size(); i++ ) {
            if (intVec.at(i)->getID() == propName) {
                intVec.at(i)->increase();
                renderer_->GetVTKWindow()->Render();
                return;
            }
        }
        PrintMessage("[vtkVoreenActor]: No property named " + propName + " found!");
    }
}


void vtkVoreenActor::changeFloatPropByName(std::string procName, std::string propName, float value, int index) {
    
    tgtAssert(workspace_, "No workspace");
    if (!workspace_->getProcessorNetwork()) {
        PrintMessage("[vtkVoreenActor]: No processor network");
        return;
    }

    Processor* proxy = 0;
    proxy = workspace_->getProcessorNetwork()->getProcessor(procName);

    if (!proxy) {
        PrintMessage("[vtkVoreenActor]: No processor named " + procName + " found!");
        return;
    }

    std::vector<FloatProperty*> floatVec = proxy->getPropertiesByType<FloatProperty>();

    for (size_t i=0; i<floatVec.size(); i++ ) {
        if (floatVec.at(i)->getID() == propName) {
            if (index == 0) 
                floatVec.at(i)->set(value);
            if (index == 1)
                floatVec.at(i)->setStepping(value);
            if (index == 2)
                floatVec.at(i)->setMinValue(value);
            if (index == 3)
                floatVec.at(i)->setMaxValue(value);
            renderer_->GetVTKWindow()->Render();
            return;
        }    
   }
   PrintMessage("[vtkVoreenActor]: No property named " + propName + " found!");
}


void vtkVoreenActor::changeStringOptionPropByName(std::string procName, std::string propName, std::string key) {

    tgtAssert(workspace_, "No workspace");
    if (!workspace_->getProcessorNetwork()) {
        PrintMessage("[vtkVoreenActor]: No processor network");
        return;
    }

    Processor* proxy = 0;
    proxy = workspace_->getProcessorNetwork()->getProcessor(procName);

    if (!proxy) {
        PrintMessage("[vtkVoreenActor]: No processor named "  + procName + " found!");
        return;
    }

    std::vector<StringOptionProperty*> propVec = proxy->getPropertiesByType<StringOptionProperty>();
    
    for (size_t i=0; i<propVec.size(); i++ ) {
        if (propName.compare(propVec.at(i)->getID())== 0){
            std::vector<std::string> allowedKeys = propVec.at(i)->getKeys(); 
            bool allowed = false;
            for (size_t keyIndex=0; keyIndex<allowedKeys.size(); keyIndex++){
                if (key.compare(allowedKeys.at(keyIndex)) == 0) allowed = true;
            }
            if (allowed)
                propVec.at(i)->selectByKey(key);
            else {
               PrintMessage("[vtkVoreenActor]: Value not allowed!"); 
            }
            renderer_->GetVTKWindow()->Render();
            return;
        }
        
    }    
    PrintMessage("[vtkVoreenActor]: No property named " + propName + " found!"); 
} 


void vtkVoreenActor::changeStringOptionPropByName(std::string procName, std::string propName, int key_id) {
    
    tgtAssert(workspace_, "No workspace");
    if (!workspace_->getProcessorNetwork()) {
        PrintMessage("[vtkVoreenActor]: No processor network");
        return;
    }
    
    Processor* proxy = 0;
    proxy = workspace_->getProcessorNetwork()->getProcessor(procName);

    if (!proxy) {
        PrintMessage("[vtkVoreenActor]: No processor named " + procName + " found!");
        return;
    }

    std::vector<StringOptionProperty*> propVec = proxy->getPropertiesByType<StringOptionProperty>();
    
    for (size_t i=0; i<propVec.size(); i++ ) {
        if (propName.compare(propVec.at(i)->getID())== 0){
            std::vector<std::string> allowedKeys = propVec.at(i)->getKeys(); 
            if (key_id < allowedKeys.size())
                propVec.at(i)->selectByKey(allowedKeys[key_id]);
            else{
                PrintMessage("[vtkVoreenActor]: Value not allowed!");
            }
            //PrintMessage("[StringProp]: Changed " + procName + ":" +propName+"( " + allowedKeys[key_id] + " )");
            renderer_->GetVTKWindow()->Render();
            return;
        }
        
    }    
    PrintMessage("[vtkVoreenActor]: No property named " + propName + " found!");

}


void vtkVoreenActor::resetTransferFunction(std::string procName) {
    
    tgtAssert(workspace_, "No workspace");
    if (!workspace_->getProcessorNetwork()) {
        PrintMessage("[vtkVoreenActor]: No processor network");
        return;
    }

    if (procName == "all") {
        std::vector<TransFuncProperty*> tfVec = workspace_->getProcessorNetwork()->getPropertiesByType<TransFuncProperty>();
        for (int i = 0; i < tfVec.size(); i++)         {
            TransFuncIntensity* tfIntensity = dynamic_cast<TransFuncIntensity*>(tfVec[i]->get());
            if (tfIntensity) {
                tfIntensity->createStdFunc();
                tfVec[i]->invalidate();
            }
        }
        PrintMessage("[vtkVoreenActor]: All TFs resettet");
        selectTransferfunctionKey(0);
        renderer_->GetVTKWindow()->Render();
        return;
    } 
    else {
        Processor* proxy = 0;
        proxy = workspace_->getProcessorNetwork()->getProcessor(procName);

        if (!proxy) {
            PrintMessage("[vtkVoreenActor]: No processor named " + procName + " found!");
            return;
        }
        std::vector<TransFuncProperty*> tfVec = proxy->getPropertiesByType<TransFuncProperty>();
        TransFuncIntensity* tfIntensity = dynamic_cast<TransFuncIntensity*>(tfVec[0]->get());
        if (tfIntensity) {
            tfIntensity->createStdFunc();
            tfVec[0]->invalidate();
            PrintMessage("[vtkVoreenActor]: Reset to Standard TF");
            selectTransferfunctionKey(0);
            renderer_->GetVTKWindow()->Render();
            
        }
    }
}

void vtkVoreenActor::clearTransferFunction(std::string procName) {

    tgtAssert(workspace_, "No workspace");
    if (!workspace_->getProcessorNetwork()) {
        PrintMessage("[vtkVoreenActor]: No processor network");
        return;
    }
    
    if (procName == "all")
        procName = raycaster_;
    Processor* proxy = 0;
    proxy = workspace_->getProcessorNetwork()->getProcessor(procName);

    if (!proxy) {
        PrintMessage("[vtkVoreenActor]: No processor named " + procName + " found!");
        return;
    }
    std::vector<TransFuncProperty*> tfVec = proxy->getPropertiesByType<TransFuncProperty>();
    TransFuncIntensity* tfIntensity = dynamic_cast<TransFuncIntensity*>(tfVec[0]->get());
    if (tfIntensity) {
        tfIntensity->clearKeys();
        tfVec[0]->invalidate();
        PrintMessage("[vtkVoreenActor]: Removed all Keys from TF");
        selectTransferfunctionKey(-1);
        renderer_->GetVTKWindow()->Render();
    }
}

void vtkVoreenActor::loadTransferFunction(std::string procName, std::string file){

    tgtAssert(workspace_, "No workspace");
    if (!workspace_->getProcessorNetwork()) {
        PrintMessage("[vtkVoreenActor]: No processor network");
        return;
    }
    
    if (procName == "all") {
        std::vector<TransFuncProperty*> tfVec = workspace_->getProcessorNetwork()->getPropertiesByType<TransFuncProperty>();  
        for (int i=0; i<tfVec.size(); i++) { 
            TransFuncIntensity* tfIntensity = dynamic_cast<TransFuncIntensity*>(tfVec[i]->get());
            if (tfIntensity) {
                tfIntensity->clearKeys();
                tfIntensity->load(file);
                tfVec[i]->set(tfIntensity);
                tfVec[i]->invalidate();
               
            }
        }
        PrintMessage("[vtkVoreenActor]: Loaded transfer function to all renderers");
        renderer_->GetVTKWindow()->Render();
    }
    else {
        Processor* proxy = 0;
        proxy = workspace_->getProcessorNetwork()->getProcessor(procName);

        if (!proxy) {
            PrintMessage("[vtkVoreenActor]: No processor named " + procName + " found!");
            return;
        }
        std::vector<TransFuncProperty*> tfVec = proxy->getPropertiesByType<TransFuncProperty>();
        TransFuncIntensity* tfIntensity = dynamic_cast<TransFuncIntensity*>(tfVec[0]->get());
        if (tfIntensity) {
            tfIntensity->clearKeys();
            tfIntensity->load(file);
            tfVec[0]->set(tfIntensity);
            tfVec[0]->invalidate();
            PrintMessage("[vtkVoreenActor]: Loaded transfer function");
            selectTransferfunctionKey(-1);
            renderer_->GetVTKWindow()->Render();
        }
    }
}

void vtkVoreenActor::selectTransferfunctionKey(int id) {

    tgtAssert(workspace_, "No workspace");
    if (!workspace_->getProcessorNetwork()) {
        PrintMessage("[vtkVoreenActor]: No processor network");
        return;
    }
    
    if (id == -1) {
        PrintMessage("[vtkVoreenActor]: Selected no Key");
        selectedKey_ = -1;
    }
    else {
        std::vector<TransFuncProperty*> tfVec = workspace_->getProcessorNetwork()->getPropertiesByType<TransFuncProperty>();
        TransFuncIntensity* tfIntensity = dynamic_cast<TransFuncIntensity*>(tfVec[0]->get());
        selectedKey_ = id;
        if (id > tfIntensity->getNumKeys()-1)
            selectedKey_ = 0;
        outstream_ << "[vtkVoreenActor]: Selected tf key " << selectedKey_ << ".\n";
        StatusOutput();
    }
}

void vtkVoreenActor::moveTransferFunctionKey(int direction) {

    tgtAssert(workspace_, "No workspace");
    if (!workspace_->getProcessorNetwork()) {
        PrintMessage("[vtkVoreenActor]: No processor network");
        return;
    }
    
    if (selectedKey_ == -1) { 
        PrintMessage("[vtkVoreenActor]: No Key Selected");
    } 
    else {
        double stepping = 0.05;
        std::vector<TransFuncProperty*> tfVec = workspace_->getProcessorNetwork()->getPropertiesByType<TransFuncProperty>();
        TransFuncIntensity* tfIntensity = dynamic_cast<TransFuncIntensity*>(tfVec[0]->get());

        TransFuncMappingKey* selectedKey = tfIntensity->getKey(selectedKey_);

        double max_left = 0;
        double max_right = 1;

        if (selectedKey_ > 0)
            max_left = tfIntensity->getKey(selectedKey_-1)->getIntensity();
        if (selectedKey_ < tfVec.size())
            max_right = tfIntensity->getKey(selectedKey_+1)->getIntensity();
        
        float i = selectedKey->getIntensity();
        float a_l = selectedKey->getAlphaL();
        float a_r = selectedKey->getAlphaR();
        //0 = left  1 = up  2 = right  3 = down
        if (direction == 0){
            i = i - stepping;
            if (i <= max_left){
                i = max_left+0.001;
            }
        } 
        else if (direction == 2) {
            i = i + stepping;
            if (i >= max_right){
                i = max_right-0.001;
            }
        } 
        else if (direction == 1){
            a_l = a_l + stepping;
            a_r = a_r + stepping;
            if (a_l > 1)
                a_l = 1;
            if (a_r > 1)
                a_r = 1;
        } 
        else if (direction == 3) {
            a_l = a_l - stepping;
            a_r = a_r - stepping;
            if (a_l < 0)
                a_l = 0;
            if (a_r < 0)
                a_r = 0;
        }

        selectedKey->setIntensity(i);
        selectedKey->setAlphaL(a_l);
        selectedKey->setAlphaR(a_r);

        tfIntensity->updateKey(selectedKey);
        
        tfVec[0]->invalidate();
        renderer_->GetVTKWindow()->Render();

        outstream_ << "[vtkVoreenActor]: Moved Key " << selectedKey_ << "(" << i << " " << a_l << ")\n";
        StatusOutput();
    }
}

int vtkVoreenActor::getSelectedKey() {
    return selectedKey_;
}

void vtkVoreenActor::addTransferFunctionKey(std::string procName, float intensity, 
                                            float red, float green, float blue, float opacity) {
    
    tgtAssert(workspace_, "No workspace");
    if (!workspace_->getProcessorNetwork()) {
        PrintMessage("[vtkVoreenActor]: No processor network");
        return;
    }
    
    if (procName == "all")
        procName = raycaster_;
    Processor* proxy = 0;
    proxy = workspace_->getProcessorNetwork()->getProcessor(procName);

    if (!proxy) {
       PrintMessage("[vtkVoreenActor]: No processor named " + procName + " found!");
        return;
    }
    std::vector<TransFuncProperty*> tfVec = proxy->getPropertiesByType<TransFuncProperty>();
    TransFuncIntensity* tfIntensity = dynamic_cast<TransFuncIntensity*> (tfVec[0]->get());
    if (tfIntensity) {
        std::vector<VolumeSource*> sources = workspace_->getProcessorNetwork()->getProcessorsByType<VolumeSource>();
        double max = sources.front()->getVolumeHandle()->getRepresentation<VolumeRAM>()->elementRange()[1];
        if (inputReader_ != 0)
             max = inputReader_->GetOutput()->GetScalarTypeMax();
        opacity = opacity*max;
        red = 255*red;
        blue = 255*blue;
        green = 255*green;
        TransFuncMappingKey* myKey = new TransFuncMappingKey(intensity, tgt::col4(red,green,blue,opacity));
        tfIntensity->addKey(myKey);
        std::vector<TransFuncMappingKey *> keys = tfIntensity->getKeys();
        tfIntensity->invalidateTexture();
        tfVec.at(0)->invalidate();
        std::stringstream tmp;
        tmp << "[vtkVoreenActor]: TF-Key added at " << intensity << ":" << opacity << "( " << red << " " << green << " " << blue << " )";
        PrintMessage(tmp.str());
        renderer_->GetVTKWindow()->Render();
    } 
    else {
        PrintMessage("[vtkVoreenActor]: Could not get transfer function.");
    }  
}

void vtkVoreenActor::setInput(vtkImageReader2* rawReader){

    if (!isInitialized()) {
        PrintMessage("[vtkVoreenActor]: Voreen not initialized");
        return;
    }
    if (inputReader_ != rawReader)
        inputReader_ = rawReader;
    else {
        vtkImageData* data = rawReader->GetOutput();        //ImageData from Reader
        std::string type = data->GetScalarTypeAsString();   //Data Type (short/char/unsigned char...)
        //Read extents and spacing
        int range[] = {0,0,0,0,0,0};
        rawReader->GetDataExtent(range);
        int x = range[1]+1;
        int y = range[3]+1;
        int z = range[5]+1;
        double spacing[] = {0.0,0.0,0.0};
        rawReader->GetDataSpacing(spacing);
        float sp_x = spacing[0];
        float sp_y = spacing[1];
        float sp_z = spacing[2]; 

        void *ptr = data->GetScalarPointer();               //Pointer to raw data

        VolumeRAM* volume = 0;

        //Cast into the nesseccary volume depending on the data type
        if (type == "unsigned char") {
            volume = new VolumeRAM_UInt8(reinterpret_cast<uint8_t*>(ptr), tgt::ivec3(x,y,z), tgt::vec3(sp_x, sp_y, sp_z));
        }
        else if (type == "unsigned short") {
            volume = new VolumeRAM_UInt16(reinterpret_cast<uint16_t*>(ptr), tgt::ivec3(x,y,z), tgt::vec3(sp_x, sp_y, sp_z));
        }
        else if (type == "char") {
            volume = new VolumeRAM_Int8(reinterpret_cast<int8_t*>(ptr), tgt::ivec3(x,y,z), tgt::vec3(sp_x, sp_y, sp_z));
        }
        else if (type == "short") {
            volume = new VolumeRAM_Int16(reinterpret_cast<int16_t*>(ptr), tgt::ivec3(x,y,z), tgt::vec3(sp_x, sp_y, sp_z));
        }
        else if (type == "float") {
            volume = new VolumeRAM_Float(reinterpret_cast<float*>(ptr), tgt::ivec3(x,y,z), tgt::vec3(sp_x, sp_y, sp_z));
        }
        else if (type == "double") {
            volume = new VolumeRAM_Double(reinterpret_cast<double*>(ptr), tgt::ivec3(x,y,z), tgt::vec3(sp_x, sp_y, sp_z));
        }
        else {
            PrintMessage("[vtkVoreenActor]: Unknown data type: " + type);
            return;
        }
        
        if (volume) {
            //Create new volume with new volume and add it to the volume source processor
            Volume* handle = new Volume(volume);
            ProcessorNetwork* network = workspace_->getProcessorNetwork();
            std::vector<VolumeSource*> sources = network->getProcessorsByType<VolumeSource>();

            if (sources.empty()) {
                PrintMessage("[vtkVoreenActor]: No volume source processor");
                return;
            }
            else {
                sources.front()->setVolumeHandle(handle);
                //Set clipping planes to new extent
                std::vector<CubeMeshProxyGeometry*> tmp = workspace_->getProcessorNetwork()->getProcessorsByType<CubeMeshProxyGeometry>();
                std::string cubeproxyname = tmp.at(0)->getName();
                
                changeFloatPropByName(cubeproxyname,"leftClippingPlane",range[0],0);
                changeFloatPropByName(cubeproxyname,"rightClippingPlane",range[1],0);
                changeFloatPropByName(cubeproxyname,"bottomClippingPlane",range[2],0);
                changeFloatPropByName(cubeproxyname,"topClippingPlane",range[3],0);
                changeFloatPropByName(cubeproxyname,"frontClippingPlane",range[4],0);
                changeFloatPropByName(cubeproxyname,"backClippingPlane",range[5],0);

                PrintMessage("[vtkVoreenActor]: Successfully imported Data");
            }
        }

        resetTransferFunction("all");
    }
}


bool vtkVoreenActor::isSingleVolumeRaycaster(){
    return singleVolumeRaycaster_;
}

bool vtkVoreenActor::isSimpleRaycaster(){
    return simpleRaycaster_;
}

bool vtkVoreenActor::isGlslRaycaster(){
    return glslRaycaster_;
}

void vtkVoreenActor::setRendererToSingleVolumeRayCaster() {
    
    tgtAssert(workspace_, "No workspace");
    if (!workspace_->getProcessorNetwork()) {
        PrintMessage("[vtkVoreenActor]: No processor network");
        return;
    }
    if (workspace_->getProcessorNetwork()->empty()) {
        PrintMessage("[vtkVoreenActor]: Empty processor network");
        return;
    }
    if (!workspace_->getProcessorNetwork()->getProcessor(raycaster_)) {
        PrintMessage("[vtkVoreenActor]: No processor named " + raycaster_ + " found.");
        return;
    }
    
    SingleVolumeRaycaster* raycaster = new SingleVolumeRaycaster();
    workspace_->getProcessorNetwork()->replaceProcessor(workspace_->getProcessorNetwork()->getProcessor(raycaster_),raycaster);
    if (inputReader_)
        setInput(inputReader_);
    singleVolumeRaycaster_ = true;
    simpleRaycaster_ = false;
    glslRaycaster_ = false;
    raycaster_ = "SingleVolumeRaycaster";

    PrintMessage("[vtkVoreenActor]: SingleVolumeRaycaster activated");
}

void vtkVoreenActor::setRendererToSimpleRayCaster() {
    
    tgtAssert(workspace_, "No workspace");
    if (!workspace_->getProcessorNetwork()) {
        PrintMessage("[vtkVoreenActor]: No processor network");
        return;
    }
    if (workspace_->getProcessorNetwork()->empty()) {
        PrintMessage("[vtkVoreenActor]: Empty processor network");
        return;
    }
    if (!workspace_->getProcessorNetwork()->getProcessor(raycaster_)) {
        PrintMessage("[vtkVoreenActor]: No processor named " + raycaster_ + " found.");
        return;
    }

    SimpleRaycaster* raycaster = new SimpleRaycaster();
    workspace_->getProcessorNetwork()->replaceProcessor(workspace_->getProcessorNetwork()->getProcessor(raycaster_),raycaster);
    if (inputReader_)
        setInput(inputReader_);
    singleVolumeRaycaster_ = false;
    simpleRaycaster_ = true;
    glslRaycaster_ = false;
    raycaster_ = "SimpleRaycaster";
    PrintMessage("[vtkVoreenActor]: Simple Raycaster activated");
}

void vtkVoreenActor::setRendererToGlslRayCaster() {
    
    tgtAssert(workspace_, "No workspace");
    if (!workspace_->getProcessorNetwork()) {
        PrintMessage("[vtkVoreenActor]: No processor network");
        return;
    }
    if (workspace_->getProcessorNetwork()->empty()) {
        PrintMessage("[vtkVoreenActor]: Empty processor network");
        return;
    }
    if (!workspace_->getProcessorNetwork()->getProcessor(raycaster_)) {
        PrintMessage("[vtkVoreenActor]: No processor named " + raycaster_ + " found.");
        return;
    }

    GLSLRaycaster* raycaster = new GLSLRaycaster();
    workspace_->getProcessorNetwork()->replaceProcessor(workspace_->getProcessorNetwork()->getProcessor(raycaster_),raycaster);
    if (inputReader_)
        setInput(inputReader_);
    singleVolumeRaycaster_ = false;
    simpleRaycaster_ = false;
    glslRaycaster_ = true;
    raycaster_ = "GLSLRaycaster";
    PrintMessage("[vtkVoreenActor]: GLSL-Raycaster activated");
}

void vtkVoreenActor::setCompositingDVR() {
    raycaster_ = "SingleVolumeRaycaster";
    PrintMessage("[vtkVoreenActor]: using DVR compositing");
    changeStringOptionPropByName(raycaster_,"compositing","dvr");
    compositing1_ = "dvr";
}

void vtkVoreenActor::setCompositingMIP() {
    raycaster_ = "SingleVolumeRaycaster";
    PrintMessage("[vtkVoreenActor]: using MIP compositing");
    changeStringOptionPropByName(raycaster_,"compositing","mip");
    compositing1_ = "mip";
}

void vtkVoreenActor::setCompositingISO() {
    raycaster_ = "SingleVolumeRaycaster";
    PrintMessage("[vtkVoreenActor]: using ISO compositing");
    changeStringOptionPropByName(raycaster_,"compositing","iso");
    compositing1_ = "iso";
}

void vtkVoreenActor::setCompositingFHP() {
    raycaster_ = "SingleVolumeRaycaster";
    PrintMessage("[vtkVoreenActor]: using FHP compositing");
    changeStringOptionPropByName(raycaster_,"compositing","fhp");
    compositing1_ = "fhp";
}

void vtkVoreenActor::setCompositingFHN() {
    raycaster_ = "SingleVolumeRaycaster";
    PrintMessage("[vtkVoreenActor]: using FHN compositing");
    changeStringOptionPropByName(raycaster_,"compositing","fhn");
    compositing1_ = "fhn";
}

std::string vtkVoreenActor::getCompositeStyle() {
    return compositing1_;
}

void vtkVoreenActor::setShading(int i) {
    shading_ = i;
    switch (i){
        case 0: PrintMessage("[vtkVoreenActor]: using no shading"); break;
        case 1: PrintMessage("[vtkVoreenActor]: using phong shading (diffuse)"); break;
        case 2: PrintMessage("[vtkVoreenActor]: using phong shading (specular)"); break;
        case 3: PrintMessage("[vtkVoreenActor]: using phong shading (diff+ambient)"); break;
        case 4: PrintMessage("[vtkVoreenActor]: using phong shading (diff+specular)"); break;
        case 5: PrintMessage("[vtkVoreenActor]: using phong shading (full)"); break;
        case 6: PrintMessage("[vtkVoreenActor]: using toon shading"); break;
    }
    changeStringOptionPropByName(raycaster_,"shading",i);
}

int vtkVoreenActor::getShading(){
    return shading_;
}

void vtkVoreenActor::selectClippingPlane(std::string plane) {
    if (plane=="top" || plane=="bottom" || plane=="left" || plane=="right" || plane == "front" || plane == "back") {
        selectedPlane_ = plane;
        PrintMessage("[vtkVorenActor]: selected " + plane + "  ClippingPlane");
    }
    else
        PrintMessage("[vtkVoreenActor]: Wrong Clipping plane");
        
}

void vtkVoreenActor::moveClippingPlane(bool increase){
    changeFloatPropByName(geometry_, selectedPlane_ + "ClippingPlane", increase);
}

void vtkVoreenActor::setClippingStepping(int value){
    changeIntPropByName(geometry_, selectedPlane_ + "ClippingPlane",value,1);
}

void vtkVoreenActor::setClippingStepping(bool increase){
    if (increase){
        changeFloatPropByName(geometry_, selectedPlane_ + "ClippingPlane",getClippingStepping()+1,1);
    } 
    else {
        changeFloatPropByName(geometry_, selectedPlane_ + "ClippingPlane",getClippingStepping()-1,1);
    }
    outstream_ << "[vtkVoreenActor]: " << selectedPlane_ << " clipping plane stepping set to " << getClippingStepping() << "\n";
    StatusOutput();
}

int vtkVoreenActor::getClippingStepping(){

    tgtAssert(workspace_, "No workspace");
    if (!workspace_->getProcessorNetwork()) {
        PrintMessage("[vtkVoreenActor]: No processor network");
        return -1;
    }
    
    Processor* processor = 0;
    processor = workspace_->getProcessorNetwork()->getProcessor(geometry_);
    
    if (!processor) {
        outstream_ << "[vtkVoreenActor]: No processor named " << geometry_ << " found!" << endl;
        return -1;
    }
    std::vector<Property*> props = processor->getProperties();
    for (size_t i=0; i < props.size(); i++) {
        IntProperty* intProp = 0; 
        if (props.at(i)->getID() == (selectedPlane_ + "ClippingPlane")){    
            intProp = dynamic_cast<IntProperty*>(props.at(i));
            return intProp->getStepping();  
        }
    }   

    outstream_ << "[vtkVoreenActor]: No property named " << selectedPlane_ + "ClippingPlane found!" << endl;
    return -1;
}

std::string vtkVoreenActor::getSelectedClippingPlane(){
    return selectedPlane_;
}

void vtkVoreenActor::increaseSamplingRate(){
    changeFloatPropByName(raycaster_, "samplingRate",0.5,1);
    changeFloatPropByName(raycaster_,"samplingRate",true);
    //PrintMessage("[vtkVoreenActor]: Increased Sampling Rate");
}

void vtkVoreenActor::decreaseSamplingRate(){
    changeFloatPropByName(raycaster_, "samplingRate",0.5,1);
    changeFloatPropByName(raycaster_,"samplingRate",false);
    //PrintMessage("[vtkVoreenActor]: Decreased Sampling Rate");
}


void vtkVoreenActor::setSamplingRate(float val){
    changeFloatPropByName(raycaster_,"samplingRate",val,0);
    outstream_ << "[vtkVoreenActor]: Set Sampling Rate to " << val << "\n";
    StatusOutput();

}

std::string vtkVoreenActor::getRaycasterName() {
    return raycaster_;
}

void vtkVoreenActor::enableBoundingBox() {

    tgtAssert(workspace_, "No workspace");
    if (!workspace_->getProcessorNetwork()) {
        PrintMessage("[vtkVoreenActor]: No processor network");
        return;
    }
    
    ProcessorNetwork* network = workspace_->getProcessorNetwork();
    Processor* processor = 0;
    processor = workspace_->getProcessorNetwork()->getProcessor("BoundingBox");
    if (processor == 0) { //create new objects
        BoundingBoxRenderer* bounding_box = new BoundingBoxRenderer();
        GeometryProcessor* geom_proc = new GeometryProcessor();
        network->addProcessorInConnection(network->getProcessor("SingleVolumeRaycaster")->getOutports()[0],network->getProcessor("Background")->getInports()[0],geom_proc);
        network->addProcessorInConnection(network->getProcessor("VolumeSource")->getOutports()[0],geom_proc->getCoProcessorInports()[0],bounding_box);
        network->connectPorts(bounding_box->getCoProcessorOutports()[0], geom_proc->getCoProcessorInports()[0]);
        network->linkProperties<CameraProperty>(std::vector<Processor*>(), std::vector<std::string>());
        renderer_->GetRenderWindow()->Render();
    
    }
        
}

void vtkVoreenActor::statusToggle() {
    double main_viewport[4];
    renderer_->GetViewport(main_viewport);
    if (main_viewport[2] != 1){
        renderer_->SetViewport(0,0,1,1);
        textRenderer_->SetViewport(0,0,0,0);
    } 
    else {
        renderer_->SetViewport(0,0,0.7,1);
        textRenderer_->SetViewport(0.7,0,1,1);
        statusText_->SetDisplayPosition(textRenderer_->GetViewport()[0]*textRenderer_->GetRenderWindow()->GetSize()[0]+10,10);
        
        textRenderer_->Render();
        renderer_->Render();
    }
    renderer_->GetRenderWindow()->Render();    
}


int vtkVoreenActor::RenderVolumetricGeometry(vtkViewport* view) {
    
    if (!isInitialized()){
        voreen_init();
    }

    if (!networkEvaluator_) {
        PrintMessage("[vtkVoreenActor]: No network evaluator");
        return 0;
    }

    //Check current window size
    int canvas_width = canvas_->getWidth();
    int canvas_height = canvas_->getHeight();
    int win_width = renderer_->GetSize()[0];
    int win_height = renderer_->GetSize()[1];
    if (canvas_width != win_width || canvas_height != win_height){
        canvas_->sizeChanged(tgt::ivec2(renderer_->GetSize()[0],renderer_->GetSize()[1]));
        if (statusText_)
            statusText_->SetDisplayPosition( (1 - renderer_->GetViewport()[0])*win_width+10,10);
    }
    //Check if Network was changed
    if (networkChanged_){
        networkEvaluator_->unlock();
        networkChanged_ = false;
        networkEvaluator_->onNetworkChange();
    }
    //save and reset OpenGL states
    glPushAttrib(GL_ALL_ATTRIB_BITS);
    glMatrixMode(GL_MODELVIEW);
    glPushMatrix();
    glLoadIdentity();
    glMatrixMode(GL_PROJECTION);
    glPushMatrix();
    glLoadIdentity();
    glClearColor(0,0,0,0);
    glMatrixMode(GL_MODELVIEW);  
    glDisable(GL_BLEND);
    //perform voreen rendering
    networkEvaluator_->process();
    //restore OpenGL states
    glMatrixMode(GL_PROJECTION);
    glPopMatrix();
    glMatrixMode(GL_MODELVIEW);
    glPopMatrix();
    glPopAttrib();
    return 1;
 }
