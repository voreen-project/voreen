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

#ifndef VRN_TOUCHTABLETIMELINE_H
#define VRN_TOUCHTABLETIMELINE_H

#include "tgt/bounds.h"
#include "voreen/core/animation/propertykeyvalue.h"

namespace voreen {

    /**
    * KeyValue represents a keyvalue in the timeline GUI. It has a position on the timeline as well as
    * boundaries to check it has been hit
    */
    struct KeyValue{
        KeyValue(int pos, int width, int height, float time)
            : pos_(pos)
            , height_(height)
            , width_(width)
            , lastValidPos_(pos)
            , isSelected_(false)
            , time_(time)
        {
        }

        int pos_;
        int width_;
        int height_;
        int lastValidPos_;
        bool isSelected_;
        float time_;
    };
    // widget menu
class TouchTableTimeline {

public:
    TouchTableTimeline(int duration, int stepSize);
    virtual int getStepSize();
    virtual void setStepSize(int stepSize);
    virtual int getDuration();
    virtual void setDuration(int duration);
    virtual int getTimeOffset() const;
    virtual void setTimeOffset(int offset);
    virtual int getCurrentPos();
    virtual void setCurrentPos(int pos);
    virtual int getMainLine();
    virtual void setMainLine(int pos);
    virtual float getCurrentTime();
    virtual void setCurrentTime(float time);
    virtual int getTimePerStep();
    virtual std::vector<KeyValue>& getKeyValues();

    /**
    * Time at the current position of cursor in the timeline
    *
    *@return std::string in format XX:XX (min:sec)
    */
    virtual std::string getCurrentTimeAsString();

    /**
    * updates the current time accoring to the position of the cursor
    */
    virtual void updateCurrentTime();

    /**
    * calculates the time for a given position in the timeline.
    *
    *@param pos x-position of the given position in the timeline (without offset)
    *@return std::string time in format XX:XX (min:sec)
    */
    virtual std::string getTimeForPosAsString(int pos);

    /**
    * calculates the time for a given position in the timeline.
    *
    *@param pos x-position of the given position in the timeline (without offset)
    *@return float time in seconds
    */
    virtual float getTimeForPos(int pos);

    /**
    * adds a keyvalue to the timeline GUI, if no keyvalue already exists at that position
    *
    *@param pos keyValue Position as x-coordinate on the timeline
    *@param height of the keyvalue, default 20
    *@param width of the keyvalue, default 10
    *@return bool true if keyvalue could be added, false otherwise
    */
    virtual bool addKeyValueAt(int pos, int height = 30, int width = 20);

    /**
    * checks whether the given position hits a keyvalue and returns its position
    *
    *@ param xPos x-position to check for keyValue
    *@ param yPos y-position ot check for keyValue
    *@ return ´KeyValue* pointer to the hit keyValue, if no keyValue hit: 0
    */
    virtual KeyValue* keyValueHit(int xPos, int yPos);

    /**
    * returns the selected value if one is selected, 0 otherwise
    */
    virtual KeyValue* getSelectedKeyValue();

    /**
    * erases selected keyValue, if no keyValue is selected nothing will be deleted
    *
    */
    virtual void eraseSelectedKeyValue();

    /**
    * checks whether or not one of the keyValues is selected
    *
    *@return bool true if a keyValue is selected, false otherwise
    */
    virtual bool anyKeyValueSelected();

    /**
    * checks if the given position colides with the position of another keyValue
    *
    *@param pos position of the keyValue, whose position needs to be validated
    */
    virtual bool validateKeyValuePos(int pos, int lastValidPos);

protected:

    //duration of animation
    int duration_;

    //current time in seconds
    float currentTime_;

    //current position, represented by a line at the according position on the timeline
    int currentPosition_;

    //starting time visible in scrollable menu
    int timeOffset_;

    //step size to render lines representing the time steps
    int stepSize_;

    //duration of time for one stepSize in seconds
    int timePerStep_;

    //positions of keyvalues in timeline
    std::vector<KeyValue> keyValues_;

    //y-position of the main line in the timeline
    int mainLinePos_;

    std::map<int,PropertyKeyValueBase*> propertyKeyValueMapping_;            ///< mapping of gui keyvalues and intern keyvalues

};

}



#endif

