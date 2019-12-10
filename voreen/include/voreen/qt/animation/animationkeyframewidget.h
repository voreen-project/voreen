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

#ifndef ANIMATIONKEYFRAMEWIDGET_H
#define ANIMATIONKEYFRAMEWIDGET_H

#include "voreen/qt/animation/animationinfowidgetbase.h"

namespace voreen {

class AnimationTimelineWidget;
class AnimationKeyframeWidget;

//----------------------------------------------------
// AnimationKeyframeLineWidget
//----------------------------------------------------

class AnimationKeyframeLineWidget : public AnimationLineWidgetBase {
    Q_OBJECT

public:
    AnimationKeyframeLineWidget(AnimationKeyframeWidget* parent, PropertyTimeline* tl);

    // deletes the property copy and clears all qt objects
    ~AnimationKeyframeLineWidget();

    void updateProperty(PropertyKeyValueBase* key);

    void updateKeyframe(PropertyKeyValueBase* key);

    /// enables / disables the gui components according to the activation status of the property timeline in the core
    virtual void updateActivationStatus();

    /// checks if construction of the widget worked, ie. if the property and property widget are defined
    virtual bool isValid() const;

protected:

    //property and widgets for the key
    Property* property_;
    PropertyWidget* propertyWidget_;

    //remember the current key everytime the property is updated by the key
    PropertyKeyValueBase* currentKey_;

protected slots:

    void updateKeyframe();

    //void tfWidgetZoomLevelChanged();

    //void adjustTfWidgetZoomLevel();

};

//----------------------------------------------------
// AnimationKeyframeWidget
//----------------------------------------------------

class AnimationKeyframeWidget : public AnimationInfoWidgetBase {
    Q_OBJECT
public:
    AnimationKeyframeWidget(AnimationTimelineWidget* = 0);

    /// creates an internal copy of the property associated with the PropertyTimeline and adds a corresponding property widget
    virtual void addPropertyTimeline(PropertyTimeline* prop);

    /// Remove all widgets and delete internal copies of properties
    virtual void clear();

    /// Updates the values of the properties (if a key has been selected)
    virtual void updateProperties(const std::vector<std::pair<PropertyKeyValueBase*, PropertyTimeline*> >& propertyValues);

    /// this enables / disables the gui components according to the activation status of the property timelines in the core
    virtual void updateActivationStatus();

protected:

friend class AnimationKeyframeLineWidget;

    virtual void createHeader();

    virtual void resizeHeader(int);

    /// associated property timelines to their corresponding line widgets
    std::map<PropertyTimeline*, AnimationKeyframeLineWidget*> timelineWidgets_;
};

} // namespace voreen

#endif

