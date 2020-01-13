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

#include "smartlabel.h"

#include "tgt/immediatemode/immediatemode.h"

#include <algorithm>
#include <set>
#include <cmath>
#include <limits>

#include "tgt/font.h"
#include "tgt/vector.h"
#include "tgt/tgt_gl.h"

namespace voreen {

const std::string SmartLabelGroupBase::loggerCat_("voreen.plotting.SmartLabel");

namespace {
struct sorterY {
    bool operator() (const SmartLabel& lhs, const SmartLabel& rhs) {
        return (lhs.boundingBox_.center().y < rhs.boundingBox_.center().y);
    };
};

struct sorterX {
    bool operator() (const SmartLabel& lhs, const SmartLabel& rhs) {
        return (lhs.boundingBox_.center().x < rhs.boundingBox_.center().x);
    };
};
};

SmartLabel::SmartLabel(std::string text, tgt::vec3 pos, Alignment align, tgt::Color color, int size)
    : text_(text)
    , position_(pos)
    , align_(align)
    , color_(color)
    , size_(size)
{}

SmartLabelGroupBase::SmartLabelGroupBase(tgt::Font* font, float padding, tgt::Bounds bounds)
    : font_(font)
    , padding_(padding)
    , bounds_(bounds)
{}

void SmartLabelGroupBase::addLabel(std::string text, tgt::vec3 position, tgt::Color color, int size, SmartLabel::Alignment align) {
    SmartLabel sm(text, position, align, color, size);
    calculateBoundingBox(sm);
    labels_.push_back(sm);
}

void SmartLabelGroupBase::addLabel(std::string text, tgt::dvec2 position, tgt::Color color, int size, SmartLabel::Alignment align) {
    addLabel(text, tgt::vec3((float)position.x, (float)position.y, 0), color, size, align);
}

void SmartLabelGroupBase::reset() {
    labels_.clear();
}

void SmartLabelGroupBase::setBounds(tgt::Bounds bounds) {
    bounds_ = bounds;
}

void SmartLabelGroupBase::render() {
    for (std::vector<SmartLabel>::const_iterator it = labels_.begin(); it < labels_.end(); ++it)
        renderLabel(*it);
}



///////////////////////////////////////////////////////////////////////////////////////////
//
//                                 OPENGL - PART
//
//
///////////////////////////////////////////////////////////////////////////////////////////

SmartLabelGroupBaseOpenGl::SmartLabelGroupBaseOpenGl(tgt::Font* font, float padding, tgt::Bounds bounds, tgt::ivec2 screensize)
    : SmartLabelGroupBase(font,padding,bounds)
    , screensize_(screensize)
{
}

void SmartLabelGroupBaseOpenGl::setScreenSize(tgt::ivec2 screensize){
    screensize_ = screensize;
}

void SmartLabelGroupBaseOpenGl::renderSingleLabel(tgt::Font* font, const std::string& text, const tgt::vec3& position,
                                            const tgt::Color& color, int size, SmartLabel::Alignment align,
                                            float padding, tgt::ivec2 screensize) {
    font->setFontSize(size);
    font->setTextAlignment(tgt::Font::TextAlignment::TopLeft);
    font->setFontColor(color);
    tgt::vec3 drawingPos = position;
    IMode.color(color.r, color.g, color.b, color.a);
    font->setLineWidth(100000);
    tgt::vec2 labelSize = font->getSize(position, text, screensize);
    float textHeight = labelSize.y;
    float textWidth  = labelSize.x;

    switch (align) {
        case SmartLabel::TOPLEFT:
            drawingPos.x += padding;
            drawingPos.y += padding;
            break;
        case SmartLabel::TOPCENTERED:
            drawingPos.x -= textWidth/2;
            drawingPos.y += padding;
            break;
        case SmartLabel::TOPRIGHT:
            drawingPos.x -= textWidth + padding;
            drawingPos.y += padding;
            break;
        case SmartLabel::MIDDLELEFT:
            drawingPos.y -= textHeight/2;
            drawingPos.x += padding;
            break;
        case SmartLabel::CENTERED:
            drawingPos.x -= textWidth/2;
            drawingPos.y -= textHeight/2;
            break;
        case SmartLabel::MIDDLERIGHT:
            drawingPos.x -= textWidth + padding;
            drawingPos.y -= textHeight/2;
            break;
        case SmartLabel::BOTTOMLEFT:
            drawingPos.x += padding;
            drawingPos.y -= textHeight + padding;
            break;
        case SmartLabel::BOTTOMCENTERED:
            drawingPos.x -= textWidth/2;
            drawingPos.y -= textHeight + padding;
            break;
        case SmartLabel::BOTTOMRIGHT:
            drawingPos.x -= textWidth + padding;
            drawingPos.y -= textHeight + padding;
            break;
    }
    font->render(drawingPos, text, screensize);
}

void SmartLabelGroupBaseOpenGl::calculateBoundingBox(SmartLabel& sm) {

    font_->setFontSize(sm.size_);

    tgt::vec3 drawingPos = sm.position_;
    tgt::vec3 enlargement(padding_/2, padding_/2, 0);
    font_->setLineWidth(10000);
    tgt::vec2 labelSize = font_->getSize(sm.position_, sm.text_, screensize_);


    float textHeight = labelSize.y;
    float textWidth  = labelSize.x;

    switch (sm.align_) {
        case SmartLabel::TOPLEFT:
            drawingPos.x += padding_;
            drawingPos.y += padding_;
            break;
        case SmartLabel::TOPCENTERED:
            drawingPos.x -= textWidth/2;
            drawingPos.y += padding_;
            break;
        case SmartLabel::TOPRIGHT:
            drawingPos.x -= textWidth + padding_;
            drawingPos.y += padding_;
            break;
        case SmartLabel::MIDDLELEFT:
            drawingPos.x += padding_;
            drawingPos.y -= textHeight/2;
            break;
        case SmartLabel::CENTERED:
            drawingPos.x -= textWidth/2;
            drawingPos.y -= textHeight/2;
            break;
        case SmartLabel::MIDDLERIGHT:
            drawingPos.x -= textWidth + padding_;
            drawingPos.y -= textHeight/2;
            break;
        case SmartLabel::BOTTOMLEFT:
            drawingPos.x += padding_;
            drawingPos.y -= textHeight + padding_;
            break;
        case SmartLabel::BOTTOMCENTERED:
            drawingPos.x -= textWidth/2;
            drawingPos.y -= textHeight + padding_;
            break;
        case SmartLabel::BOTTOMRIGHT:
            drawingPos.x -= textWidth + padding_;
            drawingPos.y -= textHeight + padding_;
            break;
    }

    // check label position within bounds
    if (drawingPos.x + textWidth > bounds_.getURB().x - padding_)
        drawingPos.x = bounds_.getURB().x - textWidth - padding_;
    if (drawingPos.y + textHeight > bounds_.getURB().y - padding_)
        drawingPos.y = bounds_.getURB().y - textHeight - padding_;
    if (drawingPos.x < bounds_.getLLF().x + padding_)
        drawingPos.x = bounds_.getLLF().x + padding_;
    if (drawingPos.y < bounds_.getLLF().y + padding_)
        drawingPos.y = bounds_.getLLF().y + padding_;

    sm.boundingBox_ = tgt::Bounds(drawingPos - enlargement, drawingPos + enlargement + tgt::vec3(textWidth, textHeight, 0));
    sm.position_ = drawingPos;
}


void SmartLabelGroupBaseOpenGl::renderLabel(const SmartLabel& sm) {
    // note: the font size may not be smaller than 8
    font_->setFontSize(sm.size_);

    //IMode.color(sm.color_.r, sm.color_.g, sm.color_.b, sm.color_.a);
    font_->setFontColor(sm.color_);
    font_->setLineWidth(100000); // don't EVER break the line
    font_->render(sm.position_, sm.text_, screensize_);
}


/* ------------------------------------------------------------------------- */

SmartLabelGroupOpenGlVerticalMoving::SmartLabelGroupOpenGlVerticalMoving(tgt::Font* font, float padding, tgt::Bounds bounds, tgt::ivec2 screensize)
    : SmartLabelGroupBaseOpenGl(font, padding, bounds, screensize)
{}

void SmartLabelGroupOpenGlVerticalMoving::performLayout() {
    if (labels_.empty())
        return;

    // first sort array by y scale
    sorterY sorter;
    std::sort(labels_.begin(), labels_.end(), sorter);

    std::vector<SmartLabel>::iterator it;
    for (it = labels_.begin(); it < labels_.end() && (it+1) < labels_.end(); ++it) {
        tgt::Bounds bottom = it->boundingBox_;
        tgt::Bounds top = (it+1)->boundingBox_;
        float diff = bottom.getURB().y - top.getLLF().y + padding_;
        if (diff > 0) {
            (it+1)->position_.y += diff;
            calculateBoundingBox(*(it+1));
        }
    }

    // some labels may have been pushed out of bounds_, we try to fix this:
    tgt::Bounds top = labels_.rbegin()->boundingBox_;
    float diff = top.getURB().y - bounds_.getURB().y + padding_;
    if (diff > 0) {
        labels_.rbegin()->position_.y -= diff;
        calculateBoundingBox(*(labels_.rbegin()));
        for (std::vector<SmartLabel>::reverse_iterator rit = labels_.rbegin(); rit < labels_.rend() && (rit+1) < labels_.rend(); ++rit) {
            tgt::Bounds top = rit->boundingBox_;
            tgt::Bounds bottom = (rit+1)->boundingBox_;
            diff = bottom.getURB().y - top.getLLF().y + padding_;
            if (diff <= 0) {
                break;
            }
            else {
                (rit+1)->position_.y -= diff;
                calculateBoundingBox(*(rit+1));
            }
        }
    }
}


/* ------------------------------------------------------------------------- */

SmartLabelGroupOpenGlHorizontalMoving::SmartLabelGroupOpenGlHorizontalMoving(tgt::Font* font, float padding, tgt::Bounds bounds, tgt::ivec2 screensize)
: SmartLabelGroupBaseOpenGl(font, padding, bounds, screensize)
{}

void SmartLabelGroupOpenGlHorizontalMoving::performLayout() {
    if (labels_.empty())
        return;

    // first sort array by x scale
    sorterX sorter;
    std::sort(labels_.begin(), labels_.end(), sorter);

    std::vector<SmartLabel>::iterator it;

    // first get maximum count of intersecting labels
    int intersectionLevel = 1;
    for (it = labels_.begin(); it < labels_.end(); ++it) {
        tgt::Bounds current = it->boundingBox_;

        // check with how many of the following labels current label intersects
        int i = 1;
        while ((it+1) < labels_.end() && current.intersects((it+1)->boundingBox_)) {
            ++i;
            ++it;
        }
        intersectionLevel = std::max(intersectionLevel, i);
    }

    if (intersectionLevel > 1) {
        float step = (bounds_.getURB().y - bounds_.getLLF().y - (2*padding_))/intersectionLevel;
        float start = bounds_.getLLF().y + padding_;

        float* positions = new float[intersectionLevel];
        for (int i=0; i<intersectionLevel; ++i)
            positions[i] = start + step*(intersectionLevel-i-1);

        int count = 0;
        // now positon the labels in intersectionLevel different steps
        for (it = labels_.begin(); it < labels_.end(); ++it) {
            it->position_.y = positions[count];
            //calculateBoundingBox(*it);
            count = (count+1) % intersectionLevel;
        }

        delete[] positions;
    }
}



/* ------------------------------------------------------------------------- */

SmartLabelGroupOpenGlNoLayout::SmartLabelGroupOpenGlNoLayout(tgt::Font* font, float padding, tgt::ivec2 screensize)
    : SmartLabelGroupBaseOpenGl(font,
                          padding,
                          tgt::Bounds(tgt::Vector(-std::numeric_limits<float>::max(), -std::numeric_limits<float>::max(), -std::numeric_limits<float>::max()),
                                      tgt::Vector( std::numeric_limits<float>::max(),  std::numeric_limits<float>::max(),  std::numeric_limits<float>::max())),
                          screensize)
{}

void SmartLabelGroupOpenGlNoLayout::performLayout() {
    /* its a feature, not a bug: nothing happens here - really */
}


/* ------------------------------------------------------------------------- */

SmartLabelGroupOpenGlNoLayoutWithBackground::SmartLabelGroupOpenGlNoLayoutWithBackground(tgt::Font* font, float padding, tgt::Bounds bounds, tgt::Color backgroundColor, tgt::ivec2 screensize)
    : SmartLabelGroupBaseOpenGl(font, padding, bounds, screensize)
    , backgroundColor_(backgroundColor)
{}

void SmartLabelGroupOpenGlNoLayoutWithBackground::performLayout() {
    /* its a feature, not a bug: nothing happens here - really */
}

void SmartLabelGroupOpenGlNoLayoutWithBackground::render() {
    for (std::vector<SmartLabel>::const_iterator it = labels_.begin(); it < labels_.end(); ++it) {
        // first render background bounding box
        IMode.color(backgroundColor_);
        float bbWidth = it->boundingBox_.getURB().x - it->boundingBox_.getLLF().x;
        IMode.begin(tgt::ImmediateMode::TRIANGLE_STRIP);
            IMode.vertex(it->boundingBox_.getLLF().xy());
            IMode.vertex(it->boundingBox_.getLLF().x + bbWidth, it->boundingBox_.getLLF().y);
            IMode.vertex(it->boundingBox_.getURB().xy());
            IMode.vertex(it->boundingBox_.getURB().x - bbWidth, it->boundingBox_.getURB().y);
        IMode.end();

        // then render label
        renderLabel(*it);
    }
}


/* ------------------------------------------------------------------------- */
///////////////////////////////////////////////////////////////////////////////////////////
//
//                                 SVG - PART
//
//
///////////////////////////////////////////////////////////////////////////////////////////

SmartLabelGroupBaseSvg::SmartLabelGroupBaseSvg(tgt::Font* font, float padding, tgt::Bounds bounds)
    : SmartLabelGroupBase(font,padding,bounds)
{
}

const std::vector<SmartLabel>& SmartLabelGroupBaseSvg::getLabels() const {
    return labels_;
}

tgt::Font* SmartLabelGroupBaseSvg::getFont() const {
    return font_;
}

float SmartLabelGroupBaseSvg::getpadding() const {
    return padding_;
}

void SmartLabelGroupBaseSvg::renderSingleLabel(tgt::Font* font, const std::string& /*text*/, int size,
                                               SmartLabel::Alignment align, double padding, tgt::dvec3& position) {

    font->setFontSize(size);

   // tgt::vec3 drawingPos = position;
    //tgt::Bounds labelBounds = font->getBounds(position, text);
    //float textHeight = labelBounds.getURB().y - labelBounds.getLLF().y;
    //float textWidth  = labelBounds.getURB().x - labelBounds.getLLF().x;

    switch (align) {
        case SmartLabel::TOPLEFT:
            position.x += padding;
            position.y -= padding;
            break;
        case SmartLabel::TOPCENTERED:
            position.x += padding;
            position.y -= padding;
            break;
        case SmartLabel::TOPRIGHT:
            position.x -= padding;
            position.y += padding;
            break;
        case SmartLabel::MIDDLELEFT:
            position.y += - size/2.0;
            position.x += padding;
            break;
        case SmartLabel::CENTERED:
            //position.x += padding;
            position.y += - size/2.0;
            break;
        case SmartLabel::MIDDLERIGHT:
            position.x -= padding;
            position.y += - size/2.0;
            break;
        case SmartLabel::BOTTOMLEFT:
            position.x += padding;
            position.y -= padding + size;
            break;
        case SmartLabel::BOTTOMCENTERED:
            position.x += padding;
            position.y -= padding + size;
            break;
        case SmartLabel::BOTTOMRIGHT:
            position.x += -padding;
            position.y -= padding + size;
            break;
    }
}

void SmartLabelGroupBaseSvg::calculateBoundingBox(SmartLabel& sm) {
    font_->setFontSize(sm.size_);

    tgt::vec3 drawingPos = sm.position_;
    tgt::vec3 enlargement(padding_/2, padding_/2, 0);
    //tgt::Bounds labelBounds = font_->getBounds(sm.position_, sm.text_);

    //tgt::vec3 urb = labelBounds.getURB();
    //tgt::vec3 llf = labelBounds.getLLF();

    //float textHeight = urb.y - llf.y;
    //float textWidth  = urb.x - llf.x;
    //HACK(apv): calculate font size in some meaningful way
    float textHeight = sm.size_;
    float textWidth  = sm.size_*sm.text_.size()*0.6f;

    switch (sm.align_) {
        case SmartLabel::TOPLEFT:
            drawingPos.x += padding_;
            drawingPos.y += padding_;
            break;
        case SmartLabel::TOPCENTERED:
            //drawingPos.x -= textWidth/2;
            drawingPos.y += padding_;
            break;
        case SmartLabel::TOPRIGHT:
            drawingPos.x -= /*textWidth*/ + padding_;
            drawingPos.y += padding_;
            break;
        case SmartLabel::MIDDLELEFT:
            drawingPos.x += padding_;
            drawingPos.y -= sm.size_/2.f;//textHeight/2;
            break;
        case SmartLabel::CENTERED:
            drawingPos.x -= textWidth/3.f;
            drawingPos.y -= sm.size_/2.f;//textHeight/2;
            break;
        case SmartLabel::MIDDLERIGHT:
            drawingPos.x -= /*textWidth*/ + padding_;
            drawingPos.y -= sm.size_/2.f;//textHeight/2;
            break;
        case SmartLabel::BOTTOMLEFT:
            drawingPos.x += padding_;
            drawingPos.y -= /*textHeight*/ sm.size_ + padding_;
            break;
        case SmartLabel::BOTTOMCENTERED:
            //drawingPos.x -= textWidth/2;
            drawingPos.y -= /*textHeight*/ sm.size_ + padding_;
            break;
        case SmartLabel::BOTTOMRIGHT:
            drawingPos.x -= /*textWidth*/ + padding_;
            drawingPos.y -= /*textHeight*/ sm.size_ + padding_;
            break;
    }

    // check label position within bounds
    if (drawingPos.x + textWidth > bounds_.getURB().x - padding_)
        drawingPos.x = bounds_.getURB().x - textWidth - padding_;
    if (drawingPos.y + textHeight > bounds_.getURB().y - padding_)
        drawingPos.y = bounds_.getURB().y - textHeight - padding_;
    if (drawingPos.x < bounds_.getLLF().x + padding_)
        drawingPos.x = bounds_.getLLF().x + padding_;
    if (drawingPos.y < bounds_.getLLF().y + padding_)
        drawingPos.y = bounds_.getLLF().y + padding_;

    sm.boundingBox_ = tgt::Bounds(drawingPos - enlargement, drawingPos + enlargement + tgt::vec3(textWidth, textHeight, 0));
    sm.position_ = drawingPos;
}

void SmartLabelGroupBaseSvg::renderLabel(const SmartLabel& sm) {
    // note: the font size may not be smaller than 8

    font_->setFontSize(sm.size_);
}
/* ------------------------------------------------------------------------- */

SmartLabelGroupSvgHorizontalMoving::SmartLabelGroupSvgHorizontalMoving(tgt::Font* font, float padding, tgt::Bounds bounds)
    : SmartLabelGroupBaseSvg(font, padding, bounds)
{}

void SmartLabelGroupSvgHorizontalMoving::performLayout() {
    if (labels_.empty())
        return;

    // first sort array by x scale
    sorterX sorter;
    std::sort(labels_.begin(), labels_.end(), sorter);

    std::vector<SmartLabel>::iterator it;

    // first get maximum count of intersecting labels
    int intersectionLevel = 1;
    for (it = labels_.begin(); it < labels_.end(); ++it) {
        tgt::Bounds current = it->boundingBox_;

        // check with how many of the following labels current label intersects
        int i = 1;
        while ((it+1) < labels_.end() && current.intersects((it+1)->boundingBox_)) {
            ++i;
            ++it;
        }
        intersectionLevel = std::max(intersectionLevel, i);
    }

    if (intersectionLevel > 1) {
        float step = (bounds_.getURB().y - bounds_.getLLF().y - (2*padding_))/intersectionLevel;
        float start = bounds_.getLLF().y + padding_;

        float* positions = new float[intersectionLevel];
        for (int i=0; i<intersectionLevel; ++i)
            positions[i] = start + step*(intersectionLevel-i-1);

        int count = 0;
        // now positon the labels in intersectionLevel different steps
        for (it = labels_.begin(); it < labels_.end(); ++it) {
            it->position_.y = positions[count];
            //calculateBoundingBox(*it);
            count = (count+1) % intersectionLevel;
        }

        delete[] positions;
    }
}


/* ------------------------------------------------------------------------- */

SmartLabelGroupSvgVerticalMoving::SmartLabelGroupSvgVerticalMoving(tgt::Font* font, float padding, tgt::Bounds bounds)
    : SmartLabelGroupBaseSvg(font, padding, bounds)
{}

void SmartLabelGroupSvgVerticalMoving::performLayout() {
    if (labels_.empty())
        return;

    // first sort array by y scale
    sorterY sorter;
    std::sort(labels_.begin(), labels_.end(), sorter);

    std::vector<SmartLabel>::iterator it;
    for (it = labels_.begin(); it < labels_.end() && (it+1) < labels_.end(); ++it) {
        tgt::Bounds bottom = it->boundingBox_;
        tgt::Bounds top = (it+1)->boundingBox_;
        float diff = bottom.getURB().y - top.getLLF().y + padding_;
        if (diff > 0) {
            (it+1)->position_.y += diff;
            calculateBoundingBox(*(it+1));
        }
    }

    // some labels may have been pushed out of bounds_, we try to fix this:
    tgt::Bounds top = labels_.rbegin()->boundingBox_;
    float diff = top.getURB().y - bounds_.getURB().y + padding_;
    if (diff > 0) {
        labels_.rbegin()->position_.y -= diff;
        calculateBoundingBox(*(labels_.rbegin()));
        for (std::vector<SmartLabel>::reverse_iterator rit = labels_.rbegin(); rit < labels_.rend() && (rit+1) < labels_.rend(); ++rit) {
            tgt::Bounds top = rit->boundingBox_;
            tgt::Bounds bottom = (rit+1)->boundingBox_;
            diff = bottom.getURB().y - top.getLLF().y + padding_;
            if (diff <= 0) {
                break;
            }
            else {
                (rit+1)->position_.y -= diff;
                calculateBoundingBox(*(rit+1));
            }
        }
    }
}

/* ------------------------------------------------------------------------- */

SmartLabelGroupSvgNoLayout::SmartLabelGroupSvgNoLayout(tgt::Font* font, float padding)
    : SmartLabelGroupBaseSvg(font,padding,
                          tgt::Bounds(tgt::Vector(-std::numeric_limits<float>::max(), -std::numeric_limits<float>::max(), -std::numeric_limits<float>::max()),
                                      tgt::Vector( std::numeric_limits<float>::max(),  std::numeric_limits<float>::max(),  std::numeric_limits<float>::max())))
{}

void SmartLabelGroupSvgNoLayout::performLayout() {
    /* its a feature, not a bug: nothing happens here - really */
}

/* ------------------------------------------------------------------------- */

SmartLabelGroupSvgNoLayoutWithBackground::SmartLabelGroupSvgNoLayoutWithBackground(tgt::Font* font, float padding, tgt::Bounds bounds, tgt::Color backgroundColor)
    : SmartLabelGroupBaseSvg(font, padding, bounds)
    , backgroundColor_(backgroundColor)
    , polygonPoints_()
{}

void SmartLabelGroupSvgNoLayoutWithBackground::render() {
    polygonPoints_.clear();
    std::vector< tgt::Vector3<double> > points;
    points.resize(4);
    for (std::vector<SmartLabel>::const_iterator it = labels_.begin(); it < labels_.end(); ++it) {
        // first render background bounding box
        double bbWidth = it->boundingBox_.getURB().x - it->boundingBox_.getLLF().x;
        //IMode.begin(tgt::ImmediateMode::TRIANGLE_STRIP);
        //    IMode.vertex(*(tgt::vec2*)it->boundingBox_.getLLF().elem);
        //    IMode.vertex(it->boundingBox_.getLLF().x + bbWidth, it->boundingBox_.getLLF().y);
        //    IMode.vertex(*(tgt::vec2*)it->boundingBox_.getURB().elem);
        //    IMode.vertex(it->boundingBox_.getURB().x - bbWidth, it->boundingBox_.getURB().y);
        //IMode.end();
        points[0] = tgt::Vector3<double>(it->boundingBox_.getLLF());
        points[1] = tgt::Vector3<double>(it->boundingBox_.getLLF().x + bbWidth, it->boundingBox_.getLLF().y,0);
        points[2] = tgt::Vector3<double>(it->boundingBox_.getURB());
        points[3] = tgt::Vector3<double>(it->boundingBox_.getURB().x - bbWidth, it->boundingBox_.getURB().y,0);
        polygonPoints_.push_back(points);
        // then render label
        renderLabel(*it);
    }
}

void SmartLabelGroupSvgNoLayoutWithBackground::performLayout() {
    /* its a feature, not a bug: nothing happens here - really */
}

const std::vector< std::vector< tgt::Vector3<double> > >& SmartLabelGroupSvgNoLayoutWithBackground::getPolygonPoints() const {
    return polygonPoints_;
}

const tgt::Color& SmartLabelGroupSvgNoLayoutWithBackground::getBackgroundColor() const {
    return backgroundColor_;
}



/* ------------------------------------------------------------------------- */
///////////////////////////////////////////////////////////////////////////////////////////
//
//                                 Latex - PART
//
//
///////////////////////////////////////////////////////////////////////////////////////////

SmartLabelGroupBaseLatex::SmartLabelGroupBaseLatex(tgt::Font* font, float padding, tgt::Bounds bounds)
    : SmartLabelGroupBase(font,padding,bounds)
{
}

const std::vector<SmartLabel>& SmartLabelGroupBaseLatex::getLabels() const {
    return labels_;
}

tgt::Font* SmartLabelGroupBaseLatex::getFont() const {
    return font_;
}

float SmartLabelGroupBaseLatex::getpadding() const {
    return padding_;
}

void SmartLabelGroupBaseLatex::renderSingleLabel(tgt::Font* font, const std::string& /*text*/, int size,
                                               SmartLabel::Alignment align, double padding, tgt::dvec3& position) {

    font->setFontSize(size);

   // tgt::vec3 drawingPos = position;
    //tgt::Bounds labelBounds = font->getBounds(position, text);
    //float textHeight = labelBounds.getURB().y - labelBounds.getLLF().y;
    //float textWidth  = labelBounds.getURB().x - labelBounds.getLLF().x;

    switch (align) {
        case SmartLabel::TOPLEFT:
            position.x += padding;
            position.y -= padding;
            break;
        case SmartLabel::TOPCENTERED:
            position.x += padding;
            position.y -= padding;
            break;
        case SmartLabel::TOPRIGHT:
            position.x -= padding;
            position.y += padding;
            break;
        case SmartLabel::MIDDLELEFT:
            position.y += - size/2.0;
            position.x += padding;
            break;
        case SmartLabel::CENTERED:
            //position.x += padding;
            position.y += - size/2.0;
            break;
        case SmartLabel::MIDDLERIGHT:
            position.x -= padding;
            position.y += - size/2.0;
            break;
        case SmartLabel::BOTTOMLEFT:
            position.x += padding;
            position.y -= padding + size;
            break;
        case SmartLabel::BOTTOMCENTERED:
            position.x += padding;
            position.y -= padding + size;
            break;
        case SmartLabel::BOTTOMRIGHT:
            position.x += -padding;
            position.y -= padding + size;
            break;
    }
}

void SmartLabelGroupBaseLatex::calculateBoundingBox(SmartLabel& sm) {
    font_->setFontSize(sm.size_);

    tgt::vec3 drawingPos = sm.position_;
    tgt::vec3 enlargement(padding_/2, padding_/2, 0);

    //HACK(apv): calculate font size in some meaningful way
    float textHeight = sm.size_;
    float textWidth  = sm.size_*sm.text_.size()*0.6f;

    //float textHeight = urb.y - llf.y;
    //float textWidth  = urb.x - llf.x;

    switch (sm.align_) {
        case SmartLabel::TOPLEFT:
            drawingPos.x += padding_;
            drawingPos.y += padding_;
            break;
        case SmartLabel::TOPCENTERED:
            //drawingPos.x -= textWidth/2;
            drawingPos.y += padding_;
            break;
        case SmartLabel::TOPRIGHT:
            drawingPos.x -= /*textWidth*/ + padding_;
            drawingPos.y += padding_;
            break;
        case SmartLabel::MIDDLELEFT:
            drawingPos.x += padding_;
            drawingPos.y -= sm.size_/2.f;//textHeight/2;
            break;
        case SmartLabel::CENTERED:
            drawingPos.x -= 2*textWidth/5.f;
            drawingPos.y -= sm.size_/2.f;//textHeight/2;
            break;
        case SmartLabel::MIDDLERIGHT:
            drawingPos.x -= /*textWidth*/ + padding_;
            drawingPos.y -= sm.size_/2.f;//textHeight/2;
            break;
        case SmartLabel::BOTTOMLEFT:
            drawingPos.x += padding_;
            drawingPos.y -= /*textHeight*/ sm.size_ + padding_;
            break;
        case SmartLabel::BOTTOMCENTERED:
            //drawingPos.x -= textWidth/2;
            drawingPos.y -= /*textHeight*/ sm.size_ + padding_;
            break;
        case SmartLabel::BOTTOMRIGHT:
            drawingPos.x -= /*textWidth*/ + padding_;
            drawingPos.y -= /*textHeight*/ sm.size_ + padding_;
            break;
    }

    // check label position within bounds
    if (drawingPos.x + textWidth > bounds_.getURB().x - padding_)
        drawingPos.x = bounds_.getURB().x - textWidth - padding_;
    if (drawingPos.y + textHeight > bounds_.getURB().y - padding_)
        drawingPos.y = bounds_.getURB().y - textHeight - padding_;
    if (drawingPos.x < bounds_.getLLF().x + padding_)
        drawingPos.x = bounds_.getLLF().x + padding_;
    if (drawingPos.y < bounds_.getLLF().y + padding_)
        drawingPos.y = bounds_.getLLF().y + padding_;

    sm.boundingBox_ = tgt::Bounds(drawingPos - enlargement, drawingPos + enlargement + tgt::vec3(textWidth, textHeight, 0));
    sm.position_ = drawingPos;
}

void SmartLabelGroupBaseLatex::renderLabel(const SmartLabel& sm) {
    // note: the font size may not be smaller than 8
    font_->setFontSize(sm.size_);
}
/* ------------------------------------------------------------------------- */

SmartLabelGroupLatexHorizontalMoving::SmartLabelGroupLatexHorizontalMoving(tgt::Font* font, float padding, tgt::Bounds bounds)
    : SmartLabelGroupBaseLatex(font, padding, bounds)
{}

void SmartLabelGroupLatexHorizontalMoving::performLayout() {
    if (labels_.empty())
        return;

    // first sort array by x scale
    sorterX sorter;
    std::sort(labels_.begin(), labels_.end(), sorter);

    std::vector<SmartLabel>::iterator it;

    // first get maximum count of intersecting labels
    int intersectionLevel = 1;
    for (it = labels_.begin(); it < labels_.end(); ++it) {
        tgt::Bounds current = it->boundingBox_;

        // check with how many of the following labels current label intersects
        int i = 1;
        while ((it+1) < labels_.end() && current.intersects((it+1)->boundingBox_)) {
            ++i;
            ++it;
        }
        intersectionLevel = std::max(intersectionLevel, i);
    }

    if (intersectionLevel > 1) {
        float step = (bounds_.getURB().y - bounds_.getLLF().y - (2*padding_))/intersectionLevel;
        float start = bounds_.getLLF().y + padding_;

        float* positions = new float[intersectionLevel];
        for (int i=0; i<intersectionLevel; ++i)
            positions[i] = start + step*(intersectionLevel-i-1);

        int count = 0;
        // now positon the labels in intersectionLevel different steps
        for (it = labels_.begin(); it < labels_.end(); ++it) {
            it->position_.y = positions[count];
            //calculateBoundingBox(*it);
            count = (count+1) % intersectionLevel;
        }

        delete[] positions;
    }
}


/* ------------------------------------------------------------------------- */

SmartLabelGroupLatexVerticalMoving::SmartLabelGroupLatexVerticalMoving(tgt::Font* font, float padding, tgt::Bounds bounds)
    : SmartLabelGroupBaseLatex(font, padding, bounds)
{}

void SmartLabelGroupLatexVerticalMoving::performLayout() {
    if (labels_.empty())
        return;

    // first sort array by y scale
    sorterY sorter;
    std::sort(labels_.begin(), labels_.end(), sorter);

    std::vector<SmartLabel>::iterator it;
    for (it = labels_.begin(); it < labels_.end() && (it+1) < labels_.end(); ++it) {
        tgt::Bounds bottom = it->boundingBox_;
        tgt::Bounds top = (it+1)->boundingBox_;
        float diff = bottom.getURB().y - top.getLLF().y + padding_;
        if (diff > 0) {
            (it+1)->position_.y += diff;
            calculateBoundingBox(*(it+1));
        }
    }

    // some labels may have been pushed out of bounds_, we try to fix this:
    tgt::Bounds top = labels_.rbegin()->boundingBox_;
    float diff = top.getURB().y - bounds_.getURB().y + padding_;
    if (diff > 0) {
        labels_.rbegin()->position_.y -= diff;
        calculateBoundingBox(*(labels_.rbegin()));
        for (std::vector<SmartLabel>::reverse_iterator rit = labels_.rbegin(); rit < labels_.rend() && (rit+1) < labels_.rend(); ++rit) {
            tgt::Bounds top = rit->boundingBox_;
            tgt::Bounds bottom = (rit+1)->boundingBox_;
            diff = bottom.getURB().y - top.getLLF().y + padding_;
            if (diff <= 0) {
                break;
            }
            else {
                (rit+1)->position_.y -= diff;
                calculateBoundingBox(*(rit+1));
            }
        }
    }
}

/* ------------------------------------------------------------------------- */

SmartLabelGroupLatexNoLayout::SmartLabelGroupLatexNoLayout(tgt::Font* font, float padding)
    : SmartLabelGroupBaseLatex(font,padding,
                          tgt::Bounds(tgt::Vector(-std::numeric_limits<float>::max(), -std::numeric_limits<float>::max(), -std::numeric_limits<float>::max()),
                                      tgt::Vector( std::numeric_limits<float>::max(),  std::numeric_limits<float>::max(),  std::numeric_limits<float>::max())))
{}

void SmartLabelGroupLatexNoLayout::performLayout() {
    /* its a feature, not a bug: nothing happens here - really */
}

/* ------------------------------------------------------------------------- */

SmartLabelGroupLatexNoLayoutWithBackground::SmartLabelGroupLatexNoLayoutWithBackground(tgt::Font* font, float padding, tgt::Bounds bounds, tgt::Color backgroundColor)
    : SmartLabelGroupBaseLatex(font, padding, bounds)
    , backgroundColor_(backgroundColor)
    , polygonPoints_()
{}

void SmartLabelGroupLatexNoLayoutWithBackground::render() {
    polygonPoints_.clear();
    std::vector< tgt::Vector3<double> > points;
    points.resize(4);
    for (std::vector<SmartLabel>::const_iterator it = labels_.begin(); it < labels_.end(); ++it) {
        // first render background bounding box
        double bbWidth = it->boundingBox_.getURB().x - it->boundingBox_.getLLF().x;
        //IMode.begin(tgt::ImmediateMode::TRIANGLE_STRIP);
        //    IMode.vertex(*(tgt::vec2*)it->boundingBox_.getLLF().elem);
        //    IMode.vertex(it->boundingBox_.getLLF().x + bbWidth, it->boundingBox_.getLLF().y);
        //    IMode.vertex(*(tgt::vec2*)it->boundingBox_.getURB().elem);
        //    IMode.vertex(it->boundingBox_.getURB().x - bbWidth, it->boundingBox_.getURB().y);
        //IMode.end();
        points[0] = tgt::Vector3<double>(it->boundingBox_.getLLF());
        points[1] = tgt::Vector3<double>(it->boundingBox_.getLLF().x + bbWidth, it->boundingBox_.getLLF().y,0);
        points[2] = tgt::Vector3<double>(it->boundingBox_.getURB());
        points[3] = tgt::Vector3<double>(it->boundingBox_.getURB().x - bbWidth, it->boundingBox_.getURB().y,0);
        polygonPoints_.push_back(points);
        // then render label
        renderLabel(*it);
    }
}

void SmartLabelGroupLatexNoLayoutWithBackground::performLayout() {
    /* its a feature, not a bug: nothing happens here - really */
}

const std::vector< std::vector< tgt::Vector3<double> > >& SmartLabelGroupLatexNoLayoutWithBackground::getPolygonPoints() const {
    return polygonPoints_;
}

const tgt::Color& SmartLabelGroupLatexNoLayoutWithBackground::getBackgroundColor() const {
    return backgroundColor_;
}


} // namespace

