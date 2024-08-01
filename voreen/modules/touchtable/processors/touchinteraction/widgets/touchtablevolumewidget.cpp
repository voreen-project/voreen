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

#include "touchtablevolumewidget.h"

namespace voreen {

    const std::string TouchTableVolumeWidget::loggerCat_("voreen.touchtable.TouchTableVolumeWidget");

    TouchTableVolumeWidget::TouchTableVolumeWidget()
        : TouchTableMenuWidgetScrollable("volumes.png")
        , volumesDirectory_("volumesdirectory", "Volumes Directory", "Set volumes directory...", VoreenApplication::app()->getCoreResourcePath(), "", FileDialogProperty::DIRECTORY)
        , fontProp_("fontprop", "Font")
        //, volumeMenu_(tgt::ivec2(100,100), tgt::vec4(0.4f, 0.4f, 0.4f, 0.5f))
        , volumeURL_("volumeURL", "Volume","")
        , volumeInfo_("volumeinfo", "")
    {
        /*fontProp_.get()->setSize(25);
        fontProp_.get()->setLineWidth(400);
        fontProp_.get()->setVerticalTextAlignment(tgt::Font::VerticalTextAlignment::Middle);*/
        addProperty(volumeURL_);
        addProperty(volumeInfo_);
        addProperty(fontProp_);
        addProperty(volumesDirectory_);
        //volumeURL_.setVisible(false);
        volumeURL_.addInfoProperty(&volumeInfo_);
        //scrollableMenu_.setFont(fontProp_.get());
        //scrollableMenuIsOpen_.set(true);
        //scrollableMenu_.setSelectionHandler(std::bind1st(std::mem_fun(&TouchTableVolumeWidget::handleScrollableMenuSelection), this));
        fontProp_.get()->setFontColor(tgt::vec4(0.65f,0.65f,0.65f,1.f));
        fontProp_.onChange(MemberFunctionCallback<TouchTableVolumeWidget>(this, &TouchTableVolumeWidget::updateFontProp));
        volumesDirectory_.onChange(MemberFunctionCallback<TouchTableVolumeWidget>(this, &TouchTableVolumeWidget::updateScrollableMenuContent));
        updateMenuCoordinates();
        validFileExtensions_.push_back("nrrd");
        validFileExtensions_.push_back("vvd");

    }

    void TouchTableVolumeWidget::renderComponents(){

        //overlay_->renderScrollableMenu(&scrollableMenu_, getMenuLL());

    }

    void TouchTableVolumeWidget::updateComponents() {



    }

    void TouchTableVolumeWidget::initialize() {
        TouchTableMenuWidgetScrollable::initialize();
    }

    void TouchTableVolumeWidget::deinitialize() {
        TexMgr.dispose(symbolTex_);     symbolTex_ = 0;
        TouchTableMenuWidgetScrollable::deinitialize();
    }

    void TouchTableVolumeWidget::pressed() {
        std::string path = volumesDirectory_.get();
        tgt::FileSystem sys;

        if(!sys.dirExists(path)) {
            LERROR("Volume Directory not found. Please select valid directory.");
            if (overlay_)
                overlay_->setMenuWidget(0);
            setActive(false);
            return;
        }

        TouchTableMenuWidgetScrollable::pressed();
    }

    void TouchTableVolumeWidget::handleTouchPoints(const std::deque<tgt::TouchPoint>& tp) {

        //scrollableMenu_.handleTouchPoints(convertTouchPointsToMenuCoordinates(tp));
        //invalidate();
    }

    void TouchTableVolumeWidget::updateScrollableMenuPosition() {

        scrollableMenu_.setLL(tgt::ivec2(0,0));
        scrollableMenu_.setUR(menuDimensions_.get());
        scrollableMenu_.setTextColumnWidth(scrollableMenu_.getUR().x - scrollableMenu_.getLL().x);
        fontProp_.get()->setLineWidth(scrollableMenu_.getTextColumnWidth());
        //also update control element position and preset menu
        scrollableMenu_.placeContentMenu();
        scrollableMenu_.placeSlider();
    }

    void TouchTableVolumeWidget::updateScrollableMenuContent() {

        std::string path = volumesDirectory_.get();
        tgt::FileSystem sys;

        if(!sys.dirExists(path)) {
            return;
        }

        std::vector<std::string> tmpVectors = sys.listFiles(path, true);
        std::vector<std::pair<std::string, tgt::Texture*> > tmpPairs;
        std::vector<std::string>::iterator iter;

        for(iter = tmpVectors.begin(); iter != tmpVectors.end(); ++iter){
            if(isExtensionValid(tgt::FileSystem::fileExtension(*iter)))
                tmpPairs.push_back(std::make_pair(*iter, static_cast<tgt::Texture*>(0)));
        }
        scrollableMenu_.clearSelection();
        scrollableMenu_.setContent(tmpPairs);
    }

    void TouchTableVolumeWidget::updateFontProp(){
        scrollableMenu_.setFont(fontProp_.get());
    }

    void TouchTableVolumeWidget::handleScrollableMenuSelection(std::string file){

        std::string path = volumesDirectory_.get();
        path.append("\\" + file);
        volumeURL_.set(path);
        volumeURL_.invalidate();
        LINFO("Load Volume: " << path);

        //check links
        const std::vector<PropertyLink*>& links = volumeURL_.getLinks();
        for (std::vector<PropertyLink*>::const_iterator i = links.begin(); i != links.end(); ++i) {
            //get the linked property, cast it to VolumeURLProperty and load volume
            Property* src = (*i)->getSourceProperty();
            Property* dest = (*i)->getDestinationProperty();

            VolumeURLProperty* prop = (src == &volumeURL_) ? dynamic_cast<VolumeURLProperty*>(dest) : dynamic_cast<VolumeURLProperty*>(src);

            if (prop) {
                //LINFO("could cast");
                prop->set(path);
                prop->loadVolume();
            }

        }
    }
    void TouchTableVolumeWidget::setValidFileExtensions(std::vector<std::string> extensions){
        validFileExtensions_ = extensions;
    }
    std::vector<std::string> TouchTableVolumeWidget::getValidFileExtensions(){
        return validFileExtensions_;
    }

    bool TouchTableVolumeWidget::isExtensionValid(std::string extension){
        std::vector<std::string>::iterator iter;
        bool isValid = false;
        for(iter = validFileExtensions_.begin(); iter != validFileExtensions_.end(); ++iter){
            if(extension == *iter)
                isValid = true;
        }
        return isValid;
    }
}

