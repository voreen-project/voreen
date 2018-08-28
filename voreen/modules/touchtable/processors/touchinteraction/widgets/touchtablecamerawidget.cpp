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

#include "touchtablecamerawidget.h"

#include "voreen/core/datastructures/datetime.h"
#include "voreen/core/utils/voreenqualitymode.h"
#include <fstream>

namespace voreen {

const std::string TouchTableCameraWidget::loggerCat_("voreen.touchtable.TouchTableCameracWidget");


TouchTableCameraWidget::TouchTableCameraWidget()
    : TouchTableMenuWidgetScrollable("videocamera.png")
    , cam_("cam", "Camera")
    , enableOrientationProp_("enableOrientationProp","Enable Overlay (to be linked)",true)
    , orientationTypeProp_("orientationTypeProp","Overlay Type (to be linked)")
    , presetDirectory_("presetdirectory", "Preset Directory", "Set preset directory...", VoreenApplication::app()->getCoreResourcePath(), "", FileDialogProperty::DIRECTORY)
    , coronal_(25, tgt::ivec2(180, 175), 0, 1.f, tgt::vec4(0.2f), tgt::vec4(0.f), 0, false, false)
    , sagittal_(25, tgt::ivec2(225, 175), 0, 1.f, tgt::vec4(0.2f), tgt::vec4(0.f), 0, false, false)
    , axial_(25, tgt::ivec2(270, 175), 0, 1.f, tgt::vec4(0.2f), tgt::vec4(0.f), 0, false, false)
    , coronalReverse_(25, tgt::ivec2(180, 175), 0, 1.f, tgt::vec4(0.2f), tgt::vec4(0.f), 0, false, false)
    , sagittalReverse_(25, tgt::ivec2(225, 175), 0, 1.f, tgt::vec4(0.2f), tgt::vec4(0.f), 0, false, false)
    , axialReverse_(25, tgt::ivec2(270, 175), 0, 1.f, tgt::vec4(0.2f), tgt::vec4(0.f), 0, false, false)
    , orientationBox_(25, tgt::ivec2(300, 175))
    , orientationAxes_(25, tgt::ivec2(300, 200))
    , loadPreset_(25, tgt::ivec2(615, 175))
    , savePreset_(25, tgt::ivec2(180, 175), 0, 1.f, tgt::vec4(0.2f), tgt::vec4(0.f), 0, false, false)
    , rotateXElem_(25, tgt::ivec2(615, 175))
    , rotateYElem_(25, tgt::ivec2(615, 175))
    , rotateZElem_(25, tgt::ivec2(615, 175))
    , coronalTex_(0)
    , coronalInvTex_(0)
    , sagittalTex_(0)
    , sagittalInvTex_(0)
    , axialTex_(0)
    , axialInvTex_(0)
    , orientationBoxTex_(0)
    , orientationAxesTex_(0)
    , rotationXTex_(0)
    , rotationYTex_(0)
    , rotationZTex_(0)
    , presetTex_(0)
    , saveTex_(0)
    , orientationTimer_(0)
    , rotationTimer_(0)
    , messageTimer_(0)
    , eventHandler_()
    , rotateX_(false)
    , rotateY_(false)
    , rotateZ_(false)
{

    //set timer
    eventHandler_.addListenerToBack(this);
    orientationTimer_ = VoreenApplication::app()->createTimer(&eventHandler_);
    rotationTimer_ = VoreenApplication::app()->createTimer(&eventHandler_);
    messageTimer_ = VoreenApplication::app()->createTimer(&eventHandler_);

    scrollableMenuIsOpen_.set(false);

    addProperty(cam_);

    addProperty(enableOrientationProp_);
    addProperty(orientationTypeProp_);
        orientationTypeProp_.addOption("axis3d","Axis 3D",OrientationOverlay::OT_AXES_3D);
        orientationTypeProp_.addOption("col_cube","Colored Cube",OrientationOverlay::OT_COLOR_CUBE);
    enableOrientationProp_.onChange(MemberFunctionCallback<TouchTableCameraWidget>(this, &TouchTableCameraWidget::updateComponents));
    orientationTypeProp_.onChange(MemberFunctionCallback<TouchTableCameraWidget>(this, &TouchTableCameraWidget::updateComponents));

    //addProperty(renderLightSource_);
    //renderLightSource_.onChange(CallMemberAction<TouchTableCameraWidget>(this, &TouchTableCameraWidget::updateComponents));

    addProperty(presetDirectory_);

    //update content of menu if preset directory changes or menu is opened
    presetDirectory_.onChange(MemberFunctionCallback<TouchTableCameraWidget>(this, &TouchTableCameraWidget::updateScrollableMenuContent));
    scrollableMenuIsOpen_.onChange(MemberFunctionCallback<TouchTableCameraWidget>(this, &TouchTableCameraWidget::updateScrollableMenuContent));

    menuDimensions_.set(tgt::ivec2(520, 200));
}

TouchTableCameraWidget::~TouchTableCameraWidget() {
    delete orientationTimer_;
    delete rotationTimer_;
    delete messageTimer_;
}

void TouchTableCameraWidget::initialize() {
    TouchTableMenuWidget::initialize();

    coronalTex_ = TexMgr.load("y_up.png");
    coronalInvTex_ = TexMgr.load("y_down.png");
    sagittalTex_ = TexMgr.load("x_up.png");
    sagittalInvTex_ = TexMgr.load("x_down.png");
    axialTex_ = TexMgr.load("z_down.png");
    axialInvTex_ = TexMgr.load("z_up.png");

    orientationBoxTex_ = TexMgr.load("orientation_box.png");
    orientationAxesTex_ = TexMgr.load("orientation_axes.png");

    rotationXTex_ = TexMgr.load("rotate_x.png");
    rotationYTex_ = TexMgr.load("rotate_y.png");
    rotationZTex_ = TexMgr.load("rotate_z.png");

    presetTex_ = TexMgr.load("document_open.png");
    saveTex_ = TexMgr.load("document_save.png");

    sagittal_.setSymbolTexture(sagittalTex_);
    coronal_.setSymbolTexture(coronalTex_);
    axial_.setSymbolTexture(axialTex_);
    sagittalReverse_.setSymbolTexture(sagittalInvTex_);
    coronalReverse_.setSymbolTexture(coronalInvTex_);
    axialReverse_.setSymbolTexture(axialInvTex_);

    orientationBox_.setSymbolTexture(orientationBoxTex_);
    orientationAxes_.setSymbolTexture(orientationAxesTex_);

    rotateXElem_.setSymbolTexture(rotationXTex_);
    rotateYElem_.setSymbolTexture(rotationYTex_);
    rotateZElem_.setSymbolTexture(rotationZTex_);

    loadPreset_.setSymbolTexture(presetTex_);
    savePreset_.setSymbolTexture(saveTex_);
}

void TouchTableCameraWidget::deinitialize() {

    orientationTimer_->stop();
    rotationTimer_->stop();
    messageTimer_->stop();

    //clean up
    sagittal_.setSymbolTexture(0);
    coronal_.setSymbolTexture(0);
    axial_.setSymbolTexture(0);
    sagittalReverse_.setSymbolTexture(0);
    coronalReverse_.setSymbolTexture(0);
    axialReverse_.setSymbolTexture(0);

    orientationBox_.setSymbolTexture(0);
    orientationAxes_.setSymbolTexture(0);

    rotateXElem_.setSymbolTexture(0);
    rotateYElem_.setSymbolTexture(0);
    rotateZElem_.setSymbolTexture(0);

    loadPreset_.setSymbolTexture(0);
    savePreset_.setSymbolTexture(0);


    TexMgr.dispose(symbolTex_);     symbolTex_ = 0;
    TexMgr.dispose(sagittalTex_);       sagittalTex_ = 0;
    TexMgr.dispose(sagittalInvTex_);    sagittalInvTex_ = 0;
    TexMgr.dispose(coronalTex_);       coronalTex_ = 0;
    TexMgr.dispose(coronalInvTex_);     coronalInvTex_ = 0;
    TexMgr.dispose(axialTex_);    axialTex_ = 0;
    TexMgr.dispose(axialInvTex_); axialInvTex_ = 0;

    TexMgr.dispose(orientationBoxTex_);    orientationBoxTex_ = 0;
    TexMgr.dispose(orientationAxesTex_);    orientationAxesTex_ = 0;

    TexMgr.dispose(rotationXTex_);  rotationXTex_ = 0;
    TexMgr.dispose(rotationYTex_);  rotationYTex_ = 0;
    TexMgr.dispose(rotationZTex_);  rotationZTex_ = 0;

    TexMgr.dispose(presetTex_); presetTex_ = 0;
    TexMgr.dispose(saveTex_); saveTex_ = 0;

    //delete textures of old preset files
    for (int i = 0; i < cameraPresetFiles_.size(); ++i) {
        TexMgr.dispose(cameraPresetFiles_.at(i).second);
        cameraPresetFiles_.at(i).second = 0;
    }
    cameraPresetFiles_.clear();

    TouchTableMenuWidget::deinitialize();
}

void TouchTableCameraWidget::handleTouchPoints(const std::deque<tgt::TouchPoint>& tp) {

    for (std::deque<tgt::TouchPoint>::const_iterator currentTp = tp.begin(); currentTp != tp.end(); ++currentTp) {
        //TODO: currently only released points are handled
        if (currentTp->state() != tgt::TouchPoint::TouchPointReleased)
            continue;

        //convert to menu coordinates
        tgt::vec2 pos = tgt::vec2(convertToMenuCoordinates(currentTp->pos()));

        //check if one of the control elements has been hit
        if (hitsControlElement(pos, sagittal_))
            orientationChanged(SAGITTAL);

        if (hitsControlElement(pos, coronal_))
            orientationChanged(CORONAL);

        if (hitsControlElement(pos, axial_))
            orientationChanged(AXIAL);

        if (hitsControlElement(pos, sagittalReverse_))
            orientationChanged(SAGITTAL_INVERSE);

        if (hitsControlElement(pos, coronalReverse_))
            orientationChanged(CORONAL_INVERSE);

        if (hitsControlElement(pos, axialReverse_))
            orientationChanged(AXIAL_INVERSE);

        if (hitsControlElement(pos, orientationBox_)) {
            if(orientationTypeProp_.getValue() == OrientationOverlay::OT_COLOR_CUBE) {
                enableOrientationProp_.set(!enableOrientationProp_.get());
            } else {
                orientationTypeProp_.selectByValue(OrientationOverlay::OT_COLOR_CUBE);
                enableOrientationProp_.set(true);
            }
        }

        if (hitsControlElement(pos, orientationAxes_)) {
            if(orientationTypeProp_.getValue() == OrientationOverlay::OT_AXES_3D) {
                enableOrientationProp_.set(!enableOrientationProp_.get());
            } else {
                orientationTypeProp_.selectByValue(OrientationOverlay::OT_AXES_3D);
                enableOrientationProp_.set(true);
            }
        }

        if (hitsControlElement(pos, rotateXElem_)) {
            rotateX_ = !rotateX_;
            updateComponents();
            invalidate();

            setRotationTimerState();
        }

        if (hitsControlElement(pos, rotateYElem_)) {
            rotateY_ = !rotateY_;
            updateComponents();
            invalidate();

            setRotationTimerState();
        }

        if (hitsControlElement(pos, rotateZElem_)) {
            rotateZ_ = !rotateZ_;
            updateComponents();
            invalidate();

            setRotationTimerState();
        }

        if (hitsControlElement(pos, loadPreset_))
            scrollableMenuIsOpen_.set(!scrollableMenuIsOpen_.get());

        if (hitsControlElement(pos, savePreset_))
            saveCamera();
    }
}

void TouchTableCameraWidget::renderComponents() {

    //render control elements
    overlay_->renderControlElement(&sagittal_, menu_.getLL());
    overlay_->renderControlElement(&coronal_, menu_.getLL());
    overlay_->renderControlElement(&axial_, menu_.getLL());

    overlay_->renderControlElement(&sagittalReverse_, menu_.getLL());
    overlay_->renderControlElement(&coronalReverse_, menu_.getLL());
    overlay_->renderControlElement(&axialReverse_, menu_.getLL());

    overlay_->renderControlElement(&orientationBox_, menu_.getLL());
    overlay_->renderControlElement(&orientationAxes_, menu_.getLL());

    overlay_->renderControlElement(&rotateXElem_, menu_.getLL());
    overlay_->renderControlElement(&rotateYElem_, menu_.getLL());
    overlay_->renderControlElement(&rotateZElem_, menu_.getLL());

    overlay_->renderControlElement(&loadPreset_, menu_.getLL());
    overlay_->renderControlElement(&savePreset_, menu_.getLL());

    //render message box
    if (fileMessage_.display_ == true) {

        //render frame and box
        overlay_->renderColoredQuad(tgt::ivec2(8) + menu_.getLL(), tgt::ivec2(menuDimensions_.get().x - 8, scrollableMenu_.getRowHeight() * 2 + 12) + menu_.getLL(), tgt::vec4(tgt::vec3(0.85f), 1.f));
        overlay_->renderColoredQuad(tgt::ivec2(10) + menu_.getLL(), tgt::ivec2(menuDimensions_.get().x - 10, scrollableMenu_.getRowHeight() * 2 + 10) + menu_.getLL(), tgt::vec4(tgt::vec3(0.2f), 1.f));

        //save settings
        int width = static_cast<int>(fontProp_.get()->getLineWidth());
        int size = fontProp_.get()->getFontSize();
        fontProp_.get()->setLineWidth(static_cast<float>(menuDimensions_.get().x - 60));
        fontProp_.get()->setFontSize(size * 2 / 3);

        //render text
        tgt::ivec2 messagePos = tgt::ivec2(20, scrollableMenu_.getRowHeight() * 2 - fontProp_.get()->getFontSize()- 4);
        overlay_->renderMenuWidgetMessage(messagePos, menu_.getLL(), fileMessage_.message_, tgt::vec4(1.f), fontProp_.get());

        //restore settings
        fontProp_.get()->setLineWidth(static_cast<float>(width));
        fontProp_.get()->setFontSize(size);
    }
}

void TouchTableCameraWidget::updateComponents() {

    //update radius
    sagittal_.setRadius(controlElementRadius_.get());
    coronal_.setRadius(controlElementRadius_.get());
    axial_.setRadius(controlElementRadius_.get());

    sagittalReverse_.setRadius(controlElementRadius_.get());
    coronalReverse_.setRadius(controlElementRadius_.get());
    axialReverse_.setRadius(controlElementRadius_.get());

    orientationBox_.setRadius(controlElementRadius_.get());
    orientationAxes_.setRadius(controlElementRadius_.get());

    rotateXElem_.setRadius(controlElementRadius_.get());
    rotateYElem_.setRadius(controlElementRadius_.get());
    rotateZElem_.setRadius(controlElementRadius_.get());

    loadPreset_.setRadius(controlElementRadius_.get());
    savePreset_.setRadius(controlElementRadius_.get());

    orientationBox_.setActive(enableOrientationProp_.get() && (orientationTypeProp_.getValue() == OrientationOverlay::OT_COLOR_CUBE));
    orientationAxes_.setActive(enableOrientationProp_.get() && (orientationTypeProp_.getValue() == OrientationOverlay::OT_AXES_3D));

    rotateXElem_.setActive(rotateX_);
    rotateYElem_.setActive(rotateY_);
    rotateZElem_.setActive(rotateZ_);

    loadPreset_.setActive(scrollableMenuIsOpen_.get());

    //update control element coordinates within the menu
    int yVal = menuDimensions_.get().y - 10 - controlElementRadius_.get();
    sagittal_.setPosition(tgt::ivec2(10 + controlElementRadius_.get(), yVal));
    coronal_.setPosition(tgt::ivec2(15 + 3 * controlElementRadius_.get(), yVal));
    axialReverse_.setPosition(tgt::ivec2(20 + 5 * controlElementRadius_.get(), yVal));

    orientationBox_.setPosition(tgt::ivec2(menuDimensions_.get().x / 2, yVal));
    orientationAxes_.setPosition(tgt::ivec2(menuDimensions_.get().x / 2 + 5 + 2 * controlElementRadius_.get(), yVal));

    loadPreset_.setPosition(menuDimensions_.get() - tgt::ivec2(controlElementRadius_.get() + 10));
    savePreset_.setPosition(menuDimensions_.get() - tgt::ivec2(controlElementRadius_.get() + 10) - tgt::ivec2(2 * controlElementRadius_.get() + 5, 0));

    yVal = yVal - 5 - 2* controlElementRadius_.get();
    sagittalReverse_.setPosition(tgt::ivec2(10 + controlElementRadius_.get(), yVal));
    coronalReverse_.setPosition(tgt::ivec2(15 + 3 * controlElementRadius_.get(), yVal));
    axial_.setPosition(tgt::ivec2(20 + 5 * controlElementRadius_.get(), yVal));


    rotateXElem_.setPosition(tgt::ivec2(menuDimensions_.get().x / 2, yVal));
    rotateYElem_.setPosition(tgt::ivec2(menuDimensions_.get().x / 2 + 5 + 2 * controlElementRadius_.get(), yVal));
    rotateZElem_.setPosition(tgt::ivec2(menuDimensions_.get().x / 2 + 10 + 4 * controlElementRadius_.get(), yVal));

}

void TouchTableCameraWidget::orientationChanged(CameraOrientation o) {

    //set current rotation
    currentRotation_ = o;
    //store old quat
    oldQuat_ = cam_.get().getQuat();

    // start timer for animation
    orientationTimer_->stop();
    orientationTimer_->start(50, 10);

}

void TouchTableCameraWidget::setRotationTimerState() {

    if (/*rotationTimer_->isActive()*/!rotationTimer_->isStopped()) {
        if (!(rotateX_ || rotateY_ || rotateZ_)) {
            QualityMode.requestQualityMode(VoreenQualityMode::RQ_DEFAULT, this);
            rotationTimer_->stop();
        }
    }
    else {
        if (rotateX_ || rotateY_ || rotateZ_) {
            QualityMode.requestQualityMode(VoreenQualityMode::RQ_INTERACTIVE, this);
            rotationTimer_->start(25, 0);
        }
    }

}

void TouchTableCameraWidget::timerEvent(tgt::TimeEvent* te) {

    //check orientation timer
    if (te->getTimer() == orientationTimer_) {
        int count = orientationTimer_->getCount();

        if (count < 10)
            QualityMode.requestQualityMode(VoreenQualityMode::RQ_INTERACTIVE, this);

        VoreenTrackball& track = cam_.getTrackball();

        tgt::quat newQuat;
        const float c = 0.5f * sqrtf(2.0f);

        switch(currentRotation_) {
            case AXIAL:
                newQuat = tgt::quat(0.f, 0.f, 0.f, 1.f);
                break;
            case SAGITTAL:
                newQuat = tgt::quat(-0.5f, 0.5f, 0.5f, -0.5f);
                break;
            case CORONAL:
                newQuat = tgt::quat(c, 0.0f, 0.0f, c);
                break;
            case AXIAL_INVERSE:
                newQuat = tgt::quat(1.0f, 0.f, 0.f, 0.f);
                break;
            case SAGITTAL_INVERSE:
                newQuat = tgt::quat(-0.5f, -0.5f, -0.5f, -0.5f);
                break;
            case CORONAL_INVERSE:
                newQuat = tgt::normalize(tgt::quat(0.0f, c, c, 0.0f));
                break;
        }

        float t = 0.1f * static_cast<float>(count);

        tgt::quat initQuat = cam_.get().getQuat();
        initQuat.invert();

        tgt::quat tmp = tgt::slerpQuat(oldQuat_, newQuat, std::min(t, 1.f));

        track.rotate(initQuat);
        track.rotate(tmp);

        if ((count == 10) && !(rotateX_ || rotateY_ || rotateZ_))
            QualityMode.requestQualityMode(VoreenQualityMode::RQ_DEFAULT, this);

        cam_.invalidate();
    }

    if (te->getTimer() == rotationTimer_) {
        //axis rotation

        VoreenTrackball& track = cam_.getTrackball();

        if (rotateX_)
            track.rotate(tgt::vec3(1.f, 0.f, 0.f), 0.05f);

        if (rotateY_)
            track.rotate(tgt::vec3(0.f, 1.f, 0.f), 0.05f);

        if (rotateZ_)
            track.rotate(tgt::vec3(0.f, 0.f, 1.f), 0.05f);

        cam_.invalidate();
    }

    if (te->getTimer() == messageTimer_) {
        //clear message
        fileMessage_.message_ = "";
        fileMessage_.display_ = false;
        invalidate();
    }

}

void TouchTableCameraWidget::handleScrollableMenuSelection(std::string file) {

    //stop orientation change
    orientationTimer_->stop();

    if (!(rotateX_ || rotateY_ || rotateZ_))
            QualityMode.requestQualityMode(VoreenQualityMode::RQ_DEFAULT, this);

    std::string url = presetDirectory_.get() + "/" + file + ".vcm";

    XmlDeserializer d(url);
    std::fstream stream(url.c_str(), std::ios_base::in);
    d.read(stream);
    d.deserialize("Property", cam_);
    stream.close();

    scrollableMenu_.clearSelection();
    scrollableMenuIsOpen_.set(false);

    cam_.invalidate();
    //checkCameraState();

    fileMessage_.message_ = "Loaded camera preset from file:\n" + url;
    fileMessage_.display_ = true;
    messageTimer_->start(3000, 1);
}

void TouchTableCameraWidget::updateScrollableMenuPosition() {

        scrollableMenu_.setAnchor(getPosition());

        scrollableMenu_.setLL(tgt::ivec2(10, 10));
        scrollableMenu_.setUR(menuDimensions_.get() - tgt::ivec2(2* controlElementRadius_.get() + 15, 10));

        scrollableMenu_.setTextColumnWidth((scrollableMenu_.getUR().x - scrollableMenu_.getLL().x) * 3 / 4);
        scrollableMenu_.setRowHeight(fontProp_.get()->getFontSize() + 10);

}

void TouchTableCameraWidget::updateScrollableMenuContent() {
    loadPreset_.setActive(scrollableMenuIsOpen_.get());

    if (scrollableMenuIsOpen_.get()) {

        if(!tgt::FileSystem::dirExists(presetDirectory_.get())) {
            LERROR("Preset Directory not found. Please select valid directory.");
            scrollableMenuIsOpen_.set(false);
            return;
        }

        //the menu has been opened -> get the list of camera presets

        //first delete textures of old preset files
        for (int i = 0; i < cameraPresetFiles_.size(); ++i) {
            TexMgr.dispose(cameraPresetFiles_.at(i).second);
            cameraPresetFiles_.at(i).second = 0;
        }

        cameraPresetFiles_.clear();
        std::string path = presetDirectory_.get();

        //content for preset menu
        tgt::FileSystem sys;
        std::vector<std::string> filenames = sys.listFiles(path, true);

        for (int i = 0; i < filenames.size(); ++i) {

            if (endsWith(filenames.at(i), ".vcm")) {
                tgt::Texture* screenshot = 0;
                std::string url = presetDirectory_.get() + "/" + filenames.at(i);
                std::string imageUrl = url.substr(0, url.length() - 4) + ".png";
                if (tgt::FileSystem::fileExists(imageUrl))
                    screenshot = TexMgr.load(imageUrl);

                std::string displayname = filenames.at(i).substr(0, filenames.at(i).length() - 4);

                cameraPresetFiles_.push_back(std::make_pair(displayname, screenshot));

            }
        }

        scrollableMenu_.setContent(cameraPresetFiles_);
    }
}

void TouchTableCameraWidget::saveCamera() {

    DateTime dt = DateTime::now();
    std::stringstream s;
    s << dt.getMonth() << "-" << dt.getDay() << "-" << dt.getYear() << "_" << dt.getHour() << "." << dt.getMinute() << "." << dt.getSecond();

    std::string url = presetDirectory_.get() + "/" + s.str() + ".vcm";

    std::fstream stream(url.c_str(), std::ios_base::out);
    if(stream.fail()) {
        LERROR("Preset Directory does not exist. Please select valid directory.");
        return;
    }

    XmlSerializer d(url);
    d.serialize("Property", cam_);
    d.write(stream);

    stream.close();

    //save screenshot
    std::vector<Port*> inports = overlay_->getInports();
    if(!inports.empty())
        ((RenderPort*) inports.back())->saveToImage(presetDirectory_.get() + "/" + s.str() + ".png");
    else
        LERROR("Inport conntection problem with TouchTableOverlay");

    fileMessage_.message_ = "Saved camera preset to file:\n" + url;
    fileMessage_.display_ = true;
    messageTimer_->start(3000, 1);

    LINFO("Saved camera to " << url);

    invalidate();

}


} // namespace
