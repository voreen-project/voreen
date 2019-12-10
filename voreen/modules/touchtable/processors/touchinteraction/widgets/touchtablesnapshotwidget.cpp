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

#include "touchtablesnapshotwidget.h"
#include <time.h>
namespace voreen {


    const std::string TouchTableSnapShotWidget::loggerCat_("voreen.touchtable.TouchTableSnapShotWidget");


    TouchTableSnapShotWidget::TouchTableSnapShotWidget()
        : TouchTableMenuWidgetScrollable("camera.png")
        , snapshotElem_(25, tgt::ivec2(10 + 20, 175),0, 1.0f, tgt::vec4(0.0f), tgt::vec4(0.0f), 0, false, false)
        , deleteElem_(25, tgt::ivec2(10 + 20, 175),0, 1.0f, tgt::vec4(0.0f), tgt::vec4(0.0f), 1, false, false)
        , snapshotDirectory_("snapshotDirectory", "Snapshot Directory", "Set snapshot directory...", VoreenApplication::app()->getCoreResourcePath(), "", FileDialogProperty::DIRECTORY, Processor::INVALID_RESULT, Property::LOD_DEFAULT, VoreenFileWatchListener::ALWAYS_OFF)
        , snapshotIconSize_("snapshotIconSize", "Snapshot IconSize", 60, 15, 200)
        , snapshotMenuHeight_("snapshotMenuHeight", "SnapshotMenu Height", 200, 100, 1000)
        , snapshotMenuWidth_("previewLength", "SnapshotMenu Length", 200, 100, 1000)
        , snapshotPreview_(tgt::ivec2(0),0)
        , fontProperty_("fontPropertry" , "Font Propertry")
        , messageTimer_(0)
        , saveSnapshotTimer_(0)
        , eventHandler_()
        , currentFilename_("")

    {
        //set timer
        eventHandler_.addListenerToBack(this);
        messageTimer_ = VoreenApplication::app()->createTimer(&eventHandler_);
        saveSnapshotTimer_ = VoreenApplication::app()->createTimer(&eventHandler_);

        addProperty(snapshotDirectory_);
        addProperty(snapshotIconSize_);
        addProperty(snapshotMenuHeight_);
        addProperty(snapshotMenuWidth_);
        addProperty(fontProperty_);

    }

    TouchTableSnapShotWidget::~TouchTableSnapShotWidget(){
        delete messageTimer_;
        delete saveSnapshotTimer_;
    }

    void TouchTableSnapShotWidget::pressed() {
        std::string path = snapshotDirectory_.get();
        tgt::FileSystem sys;

        if(!sys.dirExists(path)) {
            LERROR("SnapShot Directory not found. Please select valid directory.");
            if (overlay_)
                overlay_->setMenuWidget(0);
            setActive(false);
            return;
        }

        TouchTableMenuWidgetScrollable::pressed();
    }

    void TouchTableSnapShotWidget::initialize() {
        TouchTableMenuWidget::initialize();

        //load textures for control elements
        snapshotTex_ = TexMgr.load("camera.png");
        deleteTex_ = TexMgr.load("delete.png");
        //set the textures to the control elements accordingly
        snapshotElem_.setSymbolTexture(snapshotTex_);
        deleteElem_.setSymbolTexture(deleteTex_);
        //load existing screenshots from harddrive into snapshotmenu
        snapshots_.clear();
        tgt::FileSystem filesystem;
        std::vector<std::string> filenames = filesystem.listFiles(snapshotDirectory_.get(), true);
        for(std::vector<std::string>::iterator filenameIter= filenames.begin(); filenameIter != filenames.end(); ++filenameIter){
            if(filesystem.fileExtension(*filenameIter) ==  "png"){
                std::string filename = *filenameIter;
                snapshots_.push_back(TouchTablePreviewIcon(tgt::ivec2(0),0,0,TexMgr.load(snapshotDirectory_.get()+ "/" +filename), filename, filesystem.baseName(filename)));
            }
        }

        //default settings
        fontProperty_.get()->setFontSize(23);
        fontProperty_.get()->setFontColor(tgt::vec4(0.65f,0.65f,0.65f,1.f));
        menuDimensions_.set(tgt::ivec2(719,288));
        snapshotMenuHeight_.set(194);
        snapshotMenuWidth_.set(355);
        updateSnapshotMenu();
        scrollableMenu_.setFont(fontProperty_.get());
        scrollableMenu_.setRowHeight(40);
        scrollableMenuIsOpen_.set(true);

    }

    void TouchTableSnapShotWidget::deinitialize() {

        messageTimer_->stop();
        saveSnapshotTimer_->stop();

        //clean up

        //dispose Textures
        TexMgr.dispose(snapshotTex_); snapshotTex_ = 0;
        TexMgr.dispose(deleteTex_); deleteTex_ = 0;
        for(std::vector<TouchTablePreviewIcon>::iterator snapshotIter = snapshots_.begin(); snapshotIter != snapshots_.end(); ++ snapshotIter){
            TexMgr.dispose(snapshotIter->getTexture());
            snapshotIter->setTexture(0);
        }

        TouchTableMenuWidget::deinitialize();
    }



    void TouchTableSnapShotWidget::handleTouchPoints(const std::deque<tgt::TouchPoint>& tp) {

        for (std::deque<tgt::TouchPoint>::const_iterator currentTp = tp.begin(); currentTp != tp.end(); ++currentTp) {

            tgt::vec2 tpPos = convertToMenuCoordinates(currentTp->pos());

            if (currentTp->state() == tgt::TouchPoint::TouchPointPressed) {


            }

            if (currentTp->state() == tgt::TouchPoint::TouchPointMoved) {

            }

            if(currentTp->state() == tgt::TouchPoint::TouchPointReleased){
                //snapshot save hit
                if(hitsControlElement(tpPos, snapshotElem_)){
                    //get current time and date
                    time_t t = std::time(0);
                    struct tm * now = std::localtime(& t);
                    std::string min = now->tm_min > 9 ? itos(now->tm_min) : "0" + itos(now->tm_min);
                    std::string sec = now->tm_sec > 9 ? itos(now->tm_sec) : "0" + itos(now->tm_sec);
                    std::string caption = itos(now->tm_hour) + "." + min + "."+  sec +  "_" + itos(now->tm_mday) + "." + itos(now->tm_mon + 1) + "." + itos(now->tm_year -100);
                    std::string fileEnding= "/" + caption + ".png";
                    currentFilename_ = caption;
                    //get texture of geometry from inport TouchTableOverlay and save to file
                    setLoadingMessage( "Save snapshot to file: \n" + snapshotDirectory_.get() + fileEnding, true);

                    saveSnapshotTimer_->start(100,1);

                    continue;
                }

                //delete element hit
                if(hitsControlElement(tpPos, deleteElem_)){
                    for(std::vector<TouchTablePreviewIcon>::iterator snapshotIter = snapshots_.begin(); snapshotIter != snapshots_.end();){
                        if(snapshotIter->isSelected()){
                            std::string filename = snapshotIter->getFilename();
                            //delete selected snapshot from snapshots-vector
                            TexMgr.dispose(snapshotIter->getTexture()); snapshotIter->setTexture(0);
                            snapshotIter = snapshots_.erase(snapshotIter);
                            //delete selected snapshot from harddrive
                            tgt::FileSystem filesystem;
                            filesystem.deleteFile(snapshotDirectory_.get() + "/" + filename);
                            //delete snapshot from preview
                            if(snapshotIter != snapshots_.end()){
                                snapshotPreview_.previewTexture_ = snapshotIter->getTexture();
                                snapshotIter->setIsSelected(true);
                            }
                            else
                                snapshotPreview_.previewTexture_ = 0;

                            updateSnapshotMenu();
                            invalidate(); break;
                        }else
                            snapshotIter++;
                    }
                }

            }
        }//end for
    }

    void TouchTableSnapShotWidget::setLoadingMessage(std::string dialog, bool show){
        fileMessage_.message_ =dialog;
        fileMessage_.display_ = show;
        messageTimer_->start(3000,1);
        invalidate();
    }

    void TouchTableSnapShotWidget::timerEvent(tgt::TimeEvent* te){
        if(te->getTimer() == messageTimer_){
            //clearMessage that as been set
            fileMessage_.message_ ="";
            fileMessage_.display_ = false;
            messageTimer_->stop();
            invalidate();
        }
        if(te->getTimer() == saveSnapshotTimer_){
            saveSnapshot(currentFilename_);
            saveSnapshotTimer_->stop();
        }
    }

    void TouchTableSnapShotWidget::saveSnapshot(std::string filename){

        std::vector<Port*> inports = overlay_->getInports();
        if(!inports.empty())
            ((RenderPort*) inports.back())->saveToImage(snapshotDirectory_.get() + "/" + filename + ".png");
        else
            LERROR("Inport conntection problem with TouchTableOverlay");

        //save snapshot as a PreviewIcon
        snapshots_.push_back(TouchTablePreviewIcon(tgt::ivec2(0,0), 0 ,0 , TexMgr.load(snapshotDirectory_.get() + "/" + filename + ".png"), "/" + filename + ".png", currentFilename_));
        //show last taken snapshot in bigger preview
        snapshotPreview_.previewTexture_ = snapshots_.back().getTexture();
        snapshots_.back().setIsSelected(true);
        invalidate();

        updateSnapshotMenu();
    }

    void TouchTableSnapShotWidget::renderComponents() {

        //render control elements
        overlay_->renderControlElement(&snapshotElem_, menu_.getLL());
        if(snapshotPreview_.previewTexture_)
            overlay_->renderControlElement(&deleteElem_, menu_.getLL());
        else
            overlay_->renderControlElement(&deleteElem_, menu_.getLL(), 0.7f); //gray out delete, when no snapshot selected

        //render snapshot preview
        overlay_->renderColoredQuad(snapshotPreview_.posLL_, snapshotPreview_.posUR_, tgt::vec4(0.f, 0.f, 0.f, 0.8f));
        if(snapshotPreview_.previewTexture_) {
            tgt::ivec2 tmpLL, tmpUR; int wT = 2, h = 1;
            //resize texture to fit preview nicely
            tgt::vec2 preSize = tgt::vec2(snapshotPreview_.posUR_ - snapshotPreview_.posLL_);
            tgt::vec2 texSize = tgt::vec2(snapshotPreview_.previewTexture_->getDimensions().xy());
            float ratioPre = preSize.x / preSize.y;
            float ratioTex = texSize.x / texSize.y;
            if(ratioPre < ratioTex) {
                int border = (preSize.y - (preSize.x/ratioTex))/2;
                tmpLL = tgt::ivec2(snapshotPreview_.posLL_.x,snapshotPreview_.posLL_.y + border);
                tmpUR = tgt::ivec2(snapshotPreview_.posUR_.x,snapshotPreview_.posUR_.y - border);
            } else {
                int border = (preSize.x - (preSize.y*ratioTex))/2;
                tmpLL = tgt::ivec2(snapshotPreview_.posLL_.x + border,snapshotPreview_.posLL_.y);
                tmpUR = tgt::ivec2(snapshotPreview_.posUR_.x - border,snapshotPreview_.posUR_.y);
            }

            overlay_->renderTexturedQuad(tmpLL, tmpUR, snapshotPreview_.previewTexture_);
        }

        //render message box
        if (fileMessage_.display_ == true) {
            //render frame and box
            overlay_->renderColoredQuad(tgt::ivec2(menu_.getLL().x, menu_.getUR().y) + tgt::ivec2(16+ 2*controlElementRadius_.get(),-15 - 2*controlElementRadius_.get() ), tgt::ivec2(menu_.getLL().x + deleteElem_.getPosition().x- controlElementRadius_.get() - 10, menu_.getUR().y-5), tgt::vec4(tgt::vec3(0.85f), 1.f));
            overlay_->renderColoredQuad(tgt::ivec2(menu_.getLL().x, menu_.getUR().y) + tgt::ivec2(18+ 2*controlElementRadius_.get(), -13 - 2*controlElementRadius_.get()), tgt::ivec2(menu_.getLL().x + deleteElem_.getPosition().x- controlElementRadius_.get() - 12, menu_.getUR().y-7), tgt::vec4(tgt::vec3(0.2f), 1.f));

            //save settings
            int width = static_cast<int>(fontProperty_.get()->getLineWidth());
            int size = fontProperty_.get()->getFontSize();
            //set settings
            fontProperty_.get()->setLineWidth(static_cast<float>(deleteElem_.getPosition().x) - 3*controlElementRadius_.get() - 20);
            fontProperty_.get()->setFontSize(size * 2 / 3);
            glColor4f(1.f,1.f,1.f,1.f);

            //render text
            tgt::ivec2 messagePos = tgt::ivec2(controlElementRadius_.get()*2+20, (menu_.getUR().y-menu_.getLL().y) - 15 - fontProperty_.get()->getFontSize());
            overlay_->renderMenuWidgetMessage(messagePos, menu_.getLL(), fileMessage_.message_, tgt::vec4(1.f), fontProperty_.get());

            //restore settings
            fontProperty_.get()->setLineWidth(static_cast<float>(width));
            fontProperty_.get()->setFontSize(size);
        }


    }

    void TouchTableSnapShotWidget::updateComponents(){
        //set positions of control elements
        int yPos = menu_.getUR().y - menu_.getLL().y - 10  - controlElementRadius_.get();
        snapshotElem_.setPosition( tgt::ivec2(10 + controlElementRadius_.get(), yPos));
        deleteElem_.setPosition(tgt::ivec2(menu_.getUR().x - menu_.getLL().x - 10 - controlElementRadius_.get(), yPos) );

        //set radius of control elements
        snapshotElem_.setRadius(controlElementRadius_.get());
        deleteElem_.setRadius(controlElementRadius_.get());

        //update taken screenshotes in snapshotmenu
        updateScrollableMenuPosition();

        //update position and size of snapshot preview
        snapshotPreview_.posLL_ = menu_.getLL() + scrollableMenu_.getLL() + tgt::ivec2(snapshotMenuWidth_.get() + 10, 0);
        snapshotPreview_.posUR_ = tgt::ivec2(menu_.getUR().x - 5, menu_.getLL().y + scrollableMenu_.getUR().y);
    }

    void TouchTableSnapShotWidget::updateScrollableMenuPosition(){
        //set position of menu
        scrollableMenu_.setAnchor(getPosition());
        int yPos = menu_.getUR().y - menu_.getLL().y - 2*controlElementRadius_.get() - 25;
        scrollableMenu_.setLL(tgt::ivec2(10, yPos - snapshotMenuHeight_.get()));
        scrollableMenu_.setUR(tgt::ivec2(10+snapshotMenuWidth_.get(), yPos));


        scrollableMenu_.placeContentMenu();
        scrollableMenu_.placeSlider();

        scrollableMenu_.setTextColumnWidth((scrollableMenu_.getUR().x - scrollableMenu_.getLL().x) * 3 / 4);
        scrollableMenu_.setRowHeight(static_cast<int>(fontProp_.get()->getFontSize() + 10));

        fontProperty_.get()->setLineWidth((float) scrollableMenu_.getTextColumnWidth());
        fontProperty_.invalidate();
    }

    void TouchTableSnapShotWidget::updateSnapshotMenu(){
        //set content of menu
        std::vector<std::pair<std::string, tgt::Texture*> > menuPairs;
        for(std::vector<TouchTablePreviewIcon>::iterator snapshotIter = snapshots_.begin(); snapshotIter != snapshots_.end(); ++snapshotIter){
            menuPairs.push_back(std::make_pair(snapshotIter->getCaption(), snapshotIter->getTexture()));
        }
        scrollableMenu_.setContent(menuPairs);
    }

    void TouchTableSnapShotWidget::handleScrollableMenuSelection(std::string filename){

        for(std::vector<TouchTablePreviewIcon>::iterator snapshotIter = snapshots_.begin() ; snapshotIter != snapshots_.end(); ++snapshotIter){
            if(snapshotIter->getCaption() == filename){
                //select new hit snapshot and show it in bigger preview
                snapshotIter->setIsSelected(true);
                snapshotPreview_.previewTexture_= snapshotIter->getTexture();
                invalidate();
            }else{
                //deselect current snapshot in Preview if necessary
                snapshotIter->setIsSelected(false);
            }
        }
    }


} // namespace

