/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2018 University of Muenster, Germany.                        *
 * Visualization and Computer Graphics Group <http://viscg.uni-muenster.de>        *
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

#include "touchtabletransfuncwidget.h"

#include "voreen/core/datastructures/transfunc/1d/preintegrationtable.h"

namespace voreen {


const std::string TouchTableTransFuncWidget::loggerCat_("voreen.touchtable.TouchTableTransFuncWidget");


TouchTableTransFuncWidget::TouchTableTransFuncWidget()
    : TouchTableMenuWidgetScrollable("colorMap.png")
    , transfunc_("transfunc", "Transfer Function")
    , internalTransFunc_("internaltf", "Internal TF (do not link!)")
    , volumePort_(Port::INPORT, "volumeport", "Volume Port")
    , presetDirectory_("presetdirectory", "Preset Directory", "Set preset directory...", VoreenApplication::app()->getCoreResourcePath(), "", FileDialogProperty::DIRECTORY)
    , renderBackground_("renderbackground", "Render Background", true)
    , renderTF_("rendertf", "Render TF", true)
    , renderHistogram_("renderhistogram", "Render Histogram", true)
    , viewLeftRight_("viewlr", "View L/R (Histogram)", tgt::vec2(0.f, 1.f), tgt::vec2(-32000.f)/*tgt::vec2(FLT_MIN)*/, tgt::vec2(32000.f)/*tgt::vec2(FLT_MAX)*/)
    , backgroundElem_(25, tgt::ivec2(10 + 20, 175))
    , tfElem_(25, tgt::ivec2(15 + 60, 175))
    , histogramElem_(25, tgt::ivec2(20 + 100, 175))
    , zoomIn_(25, tgt::ivec2(180, 175), 0, 1.f, tgt::vec4(0.f), tgt::vec4(0.f), 0, false, false)
    , zoomOut_(25, tgt::ivec2(225, 175), 0, 1.f, tgt::vec4(0.f), tgt::vec4(0.f), 0, false, false)
    , zoomReset_(25, tgt::ivec2(270, 175), 0, 1.f, tgt::vec4(0.f), tgt::vec4(0.f), 0, false, false)
    , fitToData_(25, tgt::ivec2(400, 175), 0, 1.f, tgt::vec4(0.f), tgt::vec4(0.f), 0, false, false)
    , loadPreset_(25, tgt::ivec2(615, 175))
    , domainSlider_(tgt::ivec2(30, 20), 630, 15, HORIZONTAL, 25, 20)
    , histogram_(0)
{
    addPort(volumePort_);

    addProperty(presetDirectory_);
    presetDirectory_.onChange(MemberFunctionCallback<TouchTableTransFuncWidget>(this, &TouchTableTransFuncWidget::updateScrollableMenuContent));

    addProperty(transfunc_);
    transfunc_.onChange(MemberFunctionCallback<TouchTableTransFuncWidget>(this, &TouchTableTransFuncWidget::adaptSliderPositionsToDomain));
    addProperty(renderBackground_);
    addProperty(renderTF_);
    addProperty(renderHistogram_);

    //properties set control element settings if changed
    renderBackground_.onChange(MemberFunctionCallback<TouchTableTransFuncWidget>(this, &TouchTableTransFuncWidget::updateComponents));
    renderTF_.onChange(MemberFunctionCallback<TouchTableTransFuncWidget>(this, &TouchTableTransFuncWidget::updateComponents));
    renderHistogram_.onChange(MemberFunctionCallback<TouchTableTransFuncWidget>(this, &TouchTableTransFuncWidget::updateComponents));

    //scrollableMenuIsOpen_.onChange(CallMemberAction<TouchTableTransFuncWidget>(this, &TouchTableTransFuncWidget::updateScrollableMenuContent));

    addProperty(viewLeftRight_);
    viewLeftRight_.setVisibleFlag(false);

    //symbol should not be changed: set property invisible
    symbolFile_.setVisibleFlag(false);

    addProperty(internalTransFunc_);
    internalTransFunc_.setVisibleFlag(false);

    //set default value for menu dimensions
    menuDimensions_.set(tgt::ivec2(650, 225));
}

void TouchTableTransFuncWidget::initialize() {
    TouchTableMenuWidget::initialize();

    //load textures for control elements
    gridTex_ = TexMgr.load("gridtex.png");
    mappingTex_ = TexMgr.load("mappingtex.png");
    histTex_ = TexMgr.load("histogramtex.png");

    zoomInTex_ = TexMgr.load("zoomIn.png");
    zoomOutTex_ = TexMgr.load("zoomOut.png");
    zoomResetTex_ = TexMgr.load("zoomReset.png");

    fitToDataTex_ = TexMgr.load("fitToData.png");

    presetTex_ = TexMgr.load("document_open.png");

    backgroundElem_.setSymbolTexture(gridTex_);
    tfElem_.setSymbolTexture(mappingTex_);
    histogramElem_.setSymbolTexture(histTex_);

    zoomIn_.setSymbolTexture(zoomInTex_);
    zoomOut_.setSymbolTexture(zoomOutTex_);
    zoomReset_.setSymbolTexture(zoomResetTex_);

    fitToData_.setSymbolTexture(fitToDataTex_);

    loadPreset_.setSymbolTexture(presetTex_);
}

void TouchTableTransFuncWidget::deinitialize() {

    //clean up

    backgroundElem_.setSymbolTexture(0);
    tfElem_.setSymbolTexture(0);
    histogramElem_.setSymbolTexture(0);

    TexMgr.dispose(symbolTex_);     symbolTex_ = 0;
    TexMgr.dispose(gridTex_);       gridTex_ = 0;
    TexMgr.dispose(mappingTex_);    mappingTex_ = 0;
    TexMgr.dispose(histTex_);       histTex_ = 0;
    TexMgr.dispose(zoomInTex_);     zoomInTex_ = 0;
    TexMgr.dispose(zoomOutTex_);    zoomOutTex_ = 0;
    TexMgr.dispose(zoomResetTex_);  zoomResetTex_ = 0;
    TexMgr.dispose(fitToDataTex_);  fitToDataTex_ = 0;
    TexMgr.dispose(presetTex_);     presetTex_ = 0;

    //delete transfuncs of preset files
    for (int i = 0; i < tfPresetFiles_.size(); ++i) {
        //TexMgr.dispose(tfPresetFiles_.at(i).second());
        delete tfPresetFiles_.at(i).second;
    }

    TouchTableMenuWidget::deinitialize();
}

void TouchTableTransFuncWidget::handleTouchPoints(const std::deque<tgt::TouchPoint>& tp) {

    for (std::deque<tgt::TouchPoint>::const_iterator currentTp = tp.begin(); currentTp != tp.end(); ++currentTp) {

        //convert to menu coordinate space
        tgt::ivec2 pos = convertToMenuCoordinates(tgt::ivec2(currentTp->pos()));

        if (currentTp->state() == tgt::TouchPoint::TouchPointPressed) {

            if (domainSlider_.checkIndicatorHit(pos, currentTp->id()))
                domainSliderIDs_.push_back(currentTp->id());

            continue;
        }

        if (currentTp->state() == tgt::TouchPoint::TouchPointMoved) {
            std::vector<int>::iterator i = std::find(domainSliderIDs_.begin(), domainSliderIDs_.end(), currentTp->id());

            if (i != domainSliderIDs_.end()) {

                domainSlider_.updateIndicatorPosition(pos, currentTp->id());

                float minRW = viewLeftRight_.get().x + (viewLeftRight_.get().y - viewLeftRight_.get().x) * domainSlider_.getMinMax().x;
                float maxRW = viewLeftRight_.get().x + (viewLeftRight_.get().y - viewLeftRight_.get().x) * domainSlider_.getMinMax().y;

                if (minRW >= maxRW)
                    maxRW = minRW + 0.0001f;

                transfunc_.get()->setDomain(tgt::vec2(minRW, maxRW));
                transfunc_.invalidate();

                //invalidate();
            }

            continue;
        }

        //delete slider touch points on release
        std::vector<int>::iterator i = std::find(domainSliderIDs_.begin(), domainSliderIDs_.end(), currentTp->id());
        if (i != domainSliderIDs_.end()) {
            domainSliderIDs_.erase(i);
            domainSlider_.releaseIndicator(currentTp->id());
            continue;
        }

        //check if one of the control elements has been hit
        if (hitsControlElement(pos, backgroundElem_))
            renderBackground_.set(!renderBackground_.get());

        if (hitsControlElement(pos, tfElem_))
            renderTF_.set(!renderTF_.get());

        if (hitsControlElement(pos, histogramElem_))
            renderHistogram_.set(!renderHistogram_.get());

        if (hitsControlElement(pos, zoomIn_))
            zoomIn();

        if (hitsControlElement(pos, zoomOut_))
            zoomOut();

        if (hitsControlElement(pos, zoomReset_))
            resetZoom();

        if (hitsControlElement(pos, fitToData_)) {
            transfunc_.applyDomainFromData();
            resetZoom();
        }

        if (hitsControlElement(pos, loadPreset_))
            scrollableMenuIsOpen_.set(!scrollableMenuIsOpen_.get());
    }
}

void TouchTableTransFuncWidget::process() {

    TouchTableMenuWidget::process();

    if (volumePort_.hasChanged()) {
        histogram_ = 0; transfunc_.setVolume(0);

        if (volumePort_.getData()) {
            //set the new volume to the transfer function
            transfunc_.setVolume(volumePort_.getData());
            //get histogram
            VolumeHistogramIntensity* tmp = volumePort_.getData()->getDerivedData<VolumeHistogramIntensity>();
            histogram_ = &(tmp->getHistogram());
            resetZoom();
        }
    }
}

void TouchTableTransFuncWidget::renderComponents() {

    domainSlider_.setBarOrigin(tgt::ivec2(controlElementRadius_.get() + controlElementRadius_.get()/2, controlElementRadius_.get() * 3 / 4 + 5));

    //render histogram etc.
    if (!scrollableMenuIsOpen_.get())
        renderHistogram(menu_.getLL() + tgt::ivec2(domainSlider_.getBarLL().x, 10 + domainSlider_.getIndicatorHeight()), tgt::ivec2(menu_.getLL().x + domainSlider_.getBarUR().x, menu_.getUR().y - 2 * controlElementRadius_.get() - 15), overlay_->getOutportSize(), renderBackground_.get(), renderTF_.get(), renderHistogram_.get());

    //render control elements
    overlay_->renderControlElement(&backgroundElem_, menu_.getLL());
    overlay_->renderControlElement(&tfElem_, menu_.getLL());
    overlay_->renderControlElement(&histogramElem_, menu_.getLL());

    overlay_->renderControlElement(&zoomIn_, menu_.getLL());
    overlay_->renderControlElement(&zoomOut_, menu_.getLL());
    overlay_->renderControlElement(&zoomReset_, menu_.getLL());

    overlay_->renderControlElement(&fitToData_, menu_.getLL());

    overlay_->renderControlElement(&loadPreset_, menu_.getLL());

    //render slider
    overlay_->renderSlider(&domainSlider_, menu_.getLL());
}

void TouchTableTransFuncWidget::updateComponents() {

    //update radius
    backgroundElem_.setRadius(controlElementRadius_.get());
    tfElem_.setRadius(controlElementRadius_.get());
    histogramElem_.setRadius(controlElementRadius_.get());

    zoomIn_.setRadius(controlElementRadius_.get());
    zoomOut_.setRadius(controlElementRadius_.get());
    zoomReset_.setRadius(controlElementRadius_.get());

    fitToData_.setRadius(controlElementRadius_.get());

    loadPreset_.setRadius(controlElementRadius_.get());

    //update active / inactive settings
    backgroundElem_.setActive(renderBackground_.get());
    tfElem_.setActive(renderTF_.get());
    histogramElem_.setActive(renderHistogram_.get());

    loadPreset_.setActive(scrollableMenuIsOpen_.get());

    //update control element coordinates within the menu
    int yVal = menuDimensions_.get().y - 10 - controlElementRadius_.get();
    backgroundElem_.setPosition(tgt::ivec2(10 + controlElementRadius_.get(), yVal));
    tfElem_.setPosition(tgt::ivec2(15 + 3 * controlElementRadius_.get(), yVal));
    histogramElem_.setPosition(tgt::ivec2(20 + 5 * controlElementRadius_.get(), yVal));

    int xMid = menuDimensions_.get().x / 2;
    zoomIn_.setPosition(tgt::ivec2(xMid - (2 * controlElementRadius_.get() + 5), yVal));
    zoomOut_.setPosition(tgt::ivec2(xMid, yVal));
    zoomReset_.setPosition(tgt::ivec2(xMid + (2 * controlElementRadius_.get() + 5), yVal));

    fitToData_.setPosition(tgt::ivec2((menuDimensions_.get().x - 10 - controlElementRadius_.get() + xMid + (2 * controlElementRadius_.get() + 5)) / 2, yVal));

    loadPreset_.setPosition(tgt::ivec2(menuDimensions_.get().x - 10 - controlElementRadius_.get(), yVal));

    //update domain slider size
    domainSlider_.setIndicatorHeight(controlElementRadius_.get() + controlElementRadius_.get()/2);
    domainSlider_.setIndicatorWidth(controlElementRadius_.get());
    domainSlider_.setBarOrigin(tgt::ivec2(controlElementRadius_.get() + controlElementRadius_.get()/2, controlElementRadius_.get() * 3 / 4 + 5));
    domainSlider_.setBarWidth(controlElementRadius_.get() * 3 / 4);
    domainSlider_.setBarLength(menuDimensions_.get().x - 3 * controlElementRadius_.get());
}

void TouchTableTransFuncWidget::updateScrollableMenuContent() {

    loadPreset_.setActive(scrollableMenuIsOpen_.get());

    //if (!scrollableMenuIsOpen_.get())
    //    return;

    scrollableMenu_.clearSelection();

    if (scrollableMenuIsOpen_.get()) {

        //the menu has been opened -> get the list of TF presets

        //first delete transfuncs of preset files
        for (int i = 0; i < tfPresetFiles_.size(); ++i) {
            //TexMgr.dispose(tfPresetFiles_.at(i).second);
            //tfPresetFiles_.at(i).second = 0;
            delete tfPresetFiles_.at(i).second;
            tfPresetFiles_.at(i).second = 0;
        }

        tfPresetFiles_.clear();
        std::string path = presetDirectory_.get();

        //content for preset menu
        std::vector<std::pair<std::string, tgt::Texture*> > contentVector;

        tgt::FileSystem sys;
        std::vector<std::string> filenames = sys.listFiles(path, true);

        for (int i = 0; i < filenames.size(); ++i) {
#ifdef WIN32
            internalTransFunc_.get()->load(presetDirectory_.get() + "\\" + filenames.at(i));
#else
            internalTransFunc_.get()->load(presetDirectory_.get() + "/" + filenames.at(i));
#endif
            TransFunc1DKeys* tf = 0;
            if (internalTransFunc_.get())
                tf = internalTransFunc_.get()->clone();
            tfPresetFiles_.push_back(std::make_pair(filenames.at(i), tf));

            if (dynamic_cast<TransFunc1DKeys*>(tf))
                contentVector.push_back(std::make_pair(filenames.at(i), tf->getTexture()));
            else
                contentVector.push_back(std::make_pair(filenames.at(i), static_cast<tgt::Texture*>(0)));
        }

        scrollableMenu_.setContent(contentVector);
    }
}

void TouchTableTransFuncWidget::handleScrollableMenuSelection(std::string file) {
    transfunc_.get()->load(presetDirectory_.get() + "/" + file);
    transfunc_.invalidate();

    scrollableMenuIsOpen_.set(false);

    transfunc_.applyDomainFromData();
    resetZoom();
}

void TouchTableTransFuncWidget::updateScrollableMenuPosition() {
        scrollableMenu_.setAnchor(getPosition());

        scrollableMenu_.setLL(tgt::ivec2(10, 5));
        scrollableMenu_.setUR(menuDimensions_.get() - tgt::ivec2(10, 2* controlElementRadius_.get() + 20));


        scrollableMenu_.setTextColumnWidth((scrollableMenu_.getUR().x - scrollableMenu_.getLL().x) * 3 / 4);
        scrollableMenu_.setRowHeight(fontProp_.get()->getFontSize() + 10);
}

void TouchTableTransFuncWidget::zoomIn() {
    if(transfunc_.get()) {
        TransFunc1DKeys* tf = transfunc_.get();
        float domCenter = (tf->getDomain().x + tf->getDomain().y) * 0.5f;
        float domSize = (tf->getDomain().y - tf->getDomain().x);
        float viewSize = viewLeftRight_.get().y - viewLeftRight_.get().x;
        viewSize *= 0.5f;
        if(viewSize < domSize)
            viewSize = domSize;

        viewLeftRight_.set(tgt::vec2(domCenter - (viewSize * 0.5f), domCenter + (viewSize * 0.5f)));

        adaptSlidersToZoom();
    }
}

void TouchTableTransFuncWidget::zoomOut() {
    if(transfunc_.get()) {
        TransFunc1DKeys* tf = transfunc_.get();
        float domCenter = (tf->getDomain().x + tf->getDomain().y) * 0.5f;
        float domSize = (tf->getDomain().y - tf->getDomain().x);
        float viewSize = viewLeftRight_.get().y - viewLeftRight_.get().x;
        float viewCenter = domCenter;

        viewSize *= 2.0f;
        if(viewSize > domSize)
            viewCenter = domCenter;

        viewLeftRight_.set(tgt::vec2(viewCenter - (viewSize * 0.5f), viewCenter + (viewSize * 0.5f)));

        adaptSlidersToZoom();
    }
}

void TouchTableTransFuncWidget::resetZoom() {
    if(transfunc_.get()) {
        TransFunc1DKeys* tf = transfunc_.get();
        viewLeftRight_.set(tgt::vec2(tf->getDomain().x, tf->getDomain().y));

        adaptSlidersToZoom();
    }
}

void TouchTableTransFuncWidget::adaptSlidersToZoom() {

    float vl = viewLeftRight_.get().x;
    float vs = viewLeftRight_.get().y - viewLeftRight_.get().x;
    float dl = transfunc_.get()->getDomain().x;
    float dr = transfunc_.get()->getDomain().y;

    float minNorm = (dl - vl) / vs;
    float maxNorm = (dr - vl) / vs;

    domainSlider_.setIndicatorPositions(tgt::vec2(tgt::clamp(minNorm, 0.f, 1.f), tgt::clamp(maxNorm, 0.f, 1.f)));

    invalidate();
}

void TouchTableTransFuncWidget::adaptSliderPositionsToDomain() {
    float scale = viewLeftRight_.get().y - viewLeftRight_.get().x;
    float offset = viewLeftRight_.get().x;

    float min = (transfunc_.get()->getDomain().x - offset) / scale;
    float max = (transfunc_.get()->getDomain().y - offset) / scale;

    domainSlider_.setIndicatorPositions(tgt::vec2(min, max));
}

void TouchTableTransFuncWidget::renderHistogram(tgt::ivec2 ll, tgt::ivec2 ur, tgt::ivec2 viewportSize, bool renderBackground, bool renderTF, bool renderHistogram) {

    if (!volumePort_.getData() || !transfunc_.get())
        return;

    tgt::ivec2 tmp = ll;
    ll = tgt::min(ll, ur);
    ur = tgt::max(tmp, ur);

    bool logarithmic = true;
    int width = (ur.x - ll.x);
    int height = (ur.y - ll.y);
    float viewLeft = viewLeftRight_.get().x;
    float viewRight = viewLeftRight_.get().y;

    //get the tf
    TransFunc1DKeys* tf1d = transfunc_.get();

    if (!tf1d) {
        LWARNING("Not a 1D TF!");
        return;
    }

    //set transformation to use pixel coordinates
    MatStack.matrixMode(tgt::MatrixStack::PROJECTION);
    MatStack.loadMatrix(tgt::mat4::createOrtho(0, viewportSize.x, 0, viewportSize.y, -1, 1));
    MatStack.matrixMode(tgt::MatrixStack::MODELVIEW);
    MatStack.loadIdentity();

    if (isMenuInverted()) {

        //set OpenGL transformation
        float xt = - (static_cast<float>(getMenuLL().x) + static_cast<float>(getMenuDimensions().x) / 2.f);
        float yt = - (static_cast<float>(getMenuLL().y) + static_cast<float>(getMenuDimensions().y) / 2.f);

        MatStack.translate(-xt,-yt, 0.f);
        MatStack.rotate(180.f, 0.f, 0.f, 1.f);
        MatStack.translate(xt, yt, 0.f);

    }
    MatStack.translate(static_cast<float>(ll.x), static_cast<float>(ll.y), 0.f);


    if (renderBackground) {
        glDisable(GL_DEPTH_TEST);

        float inc = static_cast<float>(width) / 10.f;
        // paint checkerboard
        for (int i = 0 ; i < 10 ; ++i) {
            glBegin(GL_QUADS);
                // Front Face
                if (i % 2)
                    glColor3f(0.6f, 0.6f, 0.6f);
                else
                    glColor3f(1.f, 1.f, 1.f);
                glVertex3f( i      * inc, 0.0f * height,  -0.5f);  // Bottom Left
                glVertex3f((i + 1) * inc, 0.0f * height,  -0.5f);  // Bottom Right
                glVertex3f((i + 1) * inc, 0.5f * height,  -0.5f);  // Top Right
                glVertex3f( i      * inc, 0.5f * height,  -0.5f);  // Top Left
            glEnd();
            glBegin(GL_QUADS);
                // Front Face
                if (i % 2)
                    glColor3f(1.f, 1.f, 1.f);
                else
                    glColor3f(0.6f, 0.6f, 0.6f);
                glVertex3f( i      * inc, 0.5f * height,  -0.5f);  // Bottom Left
                glVertex3f((i + 1) * inc, 0.5f * height,  -0.5f);  // Bottom Right
                glVertex3f((i + 1) * inc, 1.0f * height,  -0.5f);  // Top Right
                glVertex3f( i      * inc, 1.0f * height,  -0.5f);  // Top Left
            glEnd();
        }

        glEnable(GL_DEPTH_TEST);

    }

    MatStack.scale(1.0f / (viewRight - viewLeft), 1.0f, 1.0f);
    MatStack.translate(-(viewLeft * width), 0.0f, 0.0f);

    if (renderTF) {

        tgt::vec2 dom = tf1d->getDomain();

        // paint transfer function
        glEnable(GL_BLEND);
        glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
        glDisable(GL_DEPTH_TEST);

        glActiveTexture(GL_TEXTURE0);
        glEnable(GL_TEXTURE_1D);
        glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_COMBINE_EXT);
        glTexEnvf(GL_TEXTURE_ENV, GL_COMBINE_RGB_EXT, GL_REPLACE);
        tf1d->getTexture()->bind();

        MatStack.matrixMode(tgt::MatrixStack::TEXTURE);
        MatStack.loadIdentity();
        MatStack.matrixMode(tgt::MatrixStack::MODELVIEW);

        glBegin(GL_QUADS);
            glColor4f(1.f, 1.f, 1.f, 1.f);

            // Quad left of domain:
            glTexCoord1f(0.f);
            glVertex3f(static_cast<float>(ll.x), 0.f, -0.5f);
            glVertex3f(dom.x * width, 0.f, -0.5f);
            glVertex3f(dom.x * width, 1.f * height, -0.5f);
            glVertex3f(static_cast<float>(ll.x), 1.f * height, -0.5f);

            // Inside domain:
            glTexCoord1f(0.f);
            glVertex3f(dom.x * width, 0.f, -0.5f);
            glTexCoord1f(1.f);
            glVertex3f(dom.y * width , 0.f, -0.5f);
            glVertex3f(dom.y * width , 1.f * height, -0.5f);
            glTexCoord1f(0.f);
            glVertex3f(dom.x * width , 1.f * height, -0.5f);

            // Quad right of domain:
            glTexCoord1f(1.f);
            glVertex3f(dom.y * width, 0.f, -0.5f);
            glVertex3f(viewRight * width, 0.f, -0.5f);
            glVertex3f(viewRight * width, 1.f * height, -0.5f);
            glVertex3f(dom.y * width, 1.f * height, -0.5f);
        glEnd();

        glBlendFunc(GL_ONE, GL_ZERO);
        glDisable(GL_BLEND);
        glDisable(GL_TEXTURE_1D);
        glEnable(GL_DEPTH_TEST);
    }

    if (renderHistogram && histogram_) {
        // paint histogram
        glDisable(GL_DEPTH_TEST);

        int nBuckets = (int) histogram_->getNumBuckets();

        glColor3f(1.f, 0.f, 0.f);
        float bucketSize = (histogram_->getMaxValue() - histogram_->getMinValue()) / nBuckets;
        float offset = histogram_->getMinValue();
        float max = static_cast<float>(histogram_->getMaxBucket());
        if(max > 0) {
            if(logarithmic)
                max = logf(max);
            glBegin(GL_QUADS);
            for (int i = 0 ; i < nBuckets ; ++i) {
                float y = 0.0f;
                int bucket = static_cast<int>(histogram_->getBucket(i));
                if(bucket > 0) {
                    if(logarithmic) {
                        y = logf(static_cast<float>(bucket)) / max;
                    }
                    else {
                        y = bucket / max;
                    }
                }

                float leftVal = (offset + (i * bucketSize)) * width;
                float rightVal = (offset + ((i+1) * bucketSize)) * width;

                //clamp to window size
                float leftBorder = viewLeft * width;
                if (rightVal < leftBorder)
                    continue;
                else if (leftVal < leftBorder)
                    leftVal = leftBorder;

                float rightBorder = width * (viewRight-viewLeft) + viewLeft * width;
                if (leftVal > rightBorder)
                    break;
                else if (rightVal > rightBorder)
                    rightVal = rightBorder;

                glVertex3f(leftVal, 0.0f,  -0.5f);  // Bottom Left
                glVertex3f(rightVal, 0.0f,  -0.5f);  // Bottom Right
                glVertex3f(rightVal, y * height,  -0.5f);  // Top Right
                glVertex3f(leftVal, y * height,  -0.5f);  // Top Left
            }
            glEnd();
        }

        glEnable(GL_DEPTH_TEST);
    }

    MatStack.matrixMode(tgt::MatrixStack::PROJECTION);
    MatStack.loadIdentity();
    MatStack.matrixMode(tgt::MatrixStack::MODELVIEW);
    MatStack.loadIdentity();
    glColor3f(1.f,1.f,1.f);
    glEnable(GL_DEPTH_TEST);

    LGL_ERROR;
}

} // namespace
