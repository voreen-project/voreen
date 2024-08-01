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

#ifndef ANIMATIONPROPERTYLISTWIDGET_H
#define ANIMATIONPROPERTYLISTWIDGET_H

#include "voreen/qt/animation/animationinfowidgetbase.h"

namespace voreen {

class AnimationTimelineWidget;
class AnimationPropertyListWidget;

//----------------------------------------------------
// AnimationPropertyListLineWidget
//----------------------------------------------------

class AnimationPropertyListLineWidget : public AnimationLineWidgetBase {
    Q_OBJECT

public:
    AnimationPropertyListLineWidget(AnimationPropertyListWidget* parent, PropertyTimeline* tl);

    virtual ~AnimationPropertyListLineWidget();

    virtual bool isValid() const {
        return true;
    }

public slots:

    virtual void setFixedWidthSlot(int);

};

//----------------------------------------------------
// AnimationPropertyListWidget
//----------------------------------------------------

class AnimationPropertyListWidget : public AnimationInfoWidgetBase {
    Q_OBJECT
public:
    AnimationPropertyListWidget(AnimationTimelineWidget* = 0);

    /// creates an internal copy of the property associated with the PropertyTimeline and adds a corresponding property widget
    virtual void addPropertyTimeline(PropertyTimeline* prop);

    /// Remove all widgets and delete internal copies of properties
    virtual void clear();

    /// this enables / disables the gui components according to the activation status of the property timelines in the core
    virtual void updateActivationStatus();

protected:

    friend class AnimationPropertyListLineWidget;

    virtual void createHeader();

    virtual void resizeHeader(int);

    /// associated property timelines to their corresponding line widgets
    std::map<PropertyTimeline*, AnimationPropertyListLineWidget*> timelineWidgets_;

};

} // namespace voreen

#endif

