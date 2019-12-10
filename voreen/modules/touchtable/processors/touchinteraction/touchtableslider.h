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

#ifndef VRN_TOUCHTABLESLIDER_H
#define VRN_TOUCHTABLESLIDER_H

#include "tgt/vector.h"
#include "tgt/bounds.h"
#include "tgt/event/touchpoint.h"

#include <deque>

namespace voreen{

    enum Orientation{
        HORIZONTAL, VERTICAL
    };

    struct Indicator{

        int height_;                    ///<height of indicator
        int width_;                        ///<width of Indicator
        float position_;                ///<position of the indicator
        tgt::Bounds bounds_;            ///<bounds of indcator to check if indicator got hit
        int associatedTouchPointID_;    ///<ID of the touchpoint currently moving the indicator
        bool isLeft_;                    ///<for sliders with two indicators this flag is used to identify the left indicator

    };

    struct Bar{

        tgt::ivec2 origin_;            ///<origin of the bar. If bar is horizontal it is most left point, if is vertical is most upper point
        int length_;                ///<length of bar
        int width_;                    ///<width of bar
        Orientation orientation_;    ///<bar is either horizontal or vertical orientated

    };

class TouchTableSlider{

public:

    TouchTableSlider(tgt::ivec2 origin, int length, int width, Orientation orientation, int indHeight = 0, int indWidth = 0, float indPos = 0.5f);

    virtual void setIndicatorPos(float pos);
    virtual float getIndicatorPos() const;
    virtual void setBarOrigin(tgt::ivec2 origin);
    virtual tgt::ivec2 getBarOrigin() const;
    virtual void setBarLength(int length);
    virtual int getBarLength() const;
    virtual Orientation getOrientation() const;
    virtual void setOrientation(Orientation orientation);
    virtual void setIndicatorHeight(int height);
    virtual void setIndicatorWidth(int width);
    virtual tgt::ivec2 getIndicatorLL() const;
    virtual tgt::ivec2 getIndicatorUR() const;
    virtual tgt::ivec2 getBarLL() const;
    virtual tgt::ivec2 getBarUR() const;
    virtual int getIndicatorHeight() const;
    virtual int getIndicatorWidth() const;
    virtual int getBarWidth() const;
    virtual void setBarWidth(int width);

    /**
    * Updates indicator position on bar corresponing with moving touchpoint
    *
    * @param newPos new position of touchpoint associated with indicator
    */
    void updateIndicatorPosition(tgt::ivec2 newPos);

    /**
    * check if touchpoint hits indicator
    *
    * @return returns true if indicator was hit
    */
    virtual bool checkIndicatorHit(tgt::ivec2 tp);

protected:

    virtual tgt::ivec2 interpolate();


    virtual void updateIndicatorBounds();

    Indicator indicator_;
    Bar bar_;

    static const std::string loggerCat_;

};

class TouchTableTwoIndicatorSlider : public TouchTableSlider{

public:

    TouchTableTwoIndicatorSlider(tgt::ivec2 origin, int length, int width, Orientation orientation, int indHeight, int indWidth);

    virtual void setIndicatorPositions(tgt::vec2 pos);
    virtual bool checkIndicatorHit(tgt::vec2 pos, int id);
    virtual void setBarOrigin(tgt::ivec2 origin);
    virtual void setBarLength(int length);
    virtual void setIndicatorHeight(int height);
    virtual void setIndicatorWidth(int width);
    virtual int getIndicatorHeight() const;
    virtual int getIndicatorWidth() const;
    tgt::vec2 getIndicatorDimensions() const;
    std::deque<tgt::ivec2> getIndicatorPositions();

    /**
    * Updates indicator position on bar corresponing with moving touchpoint
    *
    * @param pos new position of touchpoint associated with indicator
    * @param id the id of the indicator to update
    */
    void updateIndicatorPosition(tgt::vec2 pos, int id);

    /**
    * called when touchpoint associated with indicator is released
    *
    * @param id the id of the indicator to release
    */
    virtual void releaseIndicator(int id);


    //get positions of left and right indicator
    virtual tgt::vec2 getMinMax() const;

protected:

    Indicator leftIndicator_;
    Indicator rightIndicator_;

    virtual tgt::ivec2 getIndicatorLL(Indicator* indicator) const;
    virtual tgt::ivec2    getIndicatorUR(Indicator* indicator) const;
    virtual tgt::ivec2 interpolate(Indicator* indicator);

    /**
    * updates bounds of indicator depending on its position
    *
    * @param indicator inditcator to update bounds from
    */
    virtual void updateIndicatorBounds(Indicator* indicator);
    static const std::string loggerCat_;

};

}

#endif
