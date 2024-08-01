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

#ifndef VRN_TOUCHTABLESCROLLABLETIMELINEMENU_H
#define VRN_TOUCHTABLESCROLLABLETIMELINEMENU_H

#include "touchtablemenuframe.h"
#include "touchtableslider.h"
#include "voreen/core/properties/vectorproperty.h"
#include "tgt/font.h"
#include "tgt/texture.h"
#include "touchtabletimeline.h"
#include "boost/tuple/tuple.hpp"

namespace voreen{

class TouchTableScrollableTimelineMenu : public  TouchTableMenuFrame {

public:
    TouchTableScrollableTimelineMenu(tgt::ivec2 anchor, tgt::vec4 color, int duration, int stepSize, tgt::Font* font = 0);

    virtual void setUR(tgt::ivec2 ur);
    virtual void setLL(tgt::ivec2 ll);
    virtual void handleTouchPoints(const std::deque<tgt::TouchPoint>& tp);
    virtual void setSelectionHandler(boost::function<void (boost::tuple<float,bool, float>)>);    ///< handle gets called if content menu got hit with corresponding string as argument
    virtual void placeSlider();
    virtual void placeContentMenu();
    virtual  TouchTableMenuFrame& getContentMenu();
    virtual TouchTableSlider& getSlider();
    virtual tgt::Font* getFont() const;
    virtual void setFont(tgt::Font* font);
    virtual int getContentOffset() const;
    virtual TouchTableTimeline& getTimeline();

    /**
    * sets the slider indicator to the given position. Offset of content menu is set accordingly
    *
    *@param pos position to set the slider indicator at (x-coordinate btw 0 and 1)
    */
    virtual void setIndicatorPosition(float pos);

    /**
    * adds keyValue to the timeline GUI at current position of the time cursor
    *
    *@return bool returns true if keyValue could be added, false otherwise
    */
    virtual bool addKeyValue();

protected:

    void handleContentMenu(tgt::TouchPoint& tp);
    void computeContentMenuProperties();

    int contentOffset_;
    TouchTableSlider menuSlider_;            ///< slider to scroll through menu
    tgt::ivec2 sliderMargin_;                ///< distance of slider to border of  main menu
    tgt::ivec2 menuMargin_;                    ///< distance of menu with content to border of  main menu
     TouchTableMenuFrame contentMenu_;            ///< menu with content to scroll through
    int sliderID_;                         ///< if slider got hit by Touchpoint, this is the TPs ID. Otherwise -1
    boost::function<void (boost::tuple<float,bool,float>)> handle_; ///< function getting flaot time, bool keyValuehit, int oldKeyValuePos
    tgt::Font* font_;
    TouchTableTimeline timeline_;

    std::map<int,std::pair<float, KeyValue*> > keyValueTouchPointMappping_;    ///< maps hit keyvalue with according touchpoint and saves time of keyvalue when hit

};

}

#endif
