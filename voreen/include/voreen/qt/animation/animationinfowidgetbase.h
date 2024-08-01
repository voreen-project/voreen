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

#ifndef ANIMATIONINFOWIDGETBASE_H
#define ANIMATIONINFOWIDGETBASE_H

#include <QWidget>
#include <QVBoxLayout>
#include <map>
#include <QGroupBox>
#include <QLabel>

#include "voreen/core/network/workspace.h"
#include "voreen/core/animation/animation.h"
#include "voreen/core/animation/propertykeyvalue.h"

class QVBoxLayout;
class QPushButton;

namespace voreen {

//----------------------------------------------------
// AnimationLineWidgetBase
//----------------------------------------------------

class AnimationInfoWidgetBase;

class AnimationLineWidgetBase : public QWidget {
    Q_OBJECT

public:
    AnimationLineWidgetBase(AnimationInfoWidgetBase* parent, PropertyTimeline* tl);

    // deletes the property copy and clears all qt objects
    virtual ~AnimationLineWidgetBase();

    /** Enables / disables the gui components (ie. the activation button) according to the activation status of the property timeline in the core.
     * Subclasses have to overwrite this method to activate / deactivate additional gui components, but still should call this method.
     */
    virtual void updateActivationStatus();

    /// checks if construction of the widget worked, ie. if the properties and property widgets are defined
    virtual bool isValid() const = 0;

public slots:

    /** For setting the fixed width. necessary because of the use of QScrollArea in SimpleTimelineWidget.
     *  Subclasses have to overwrite this method the set the width of additional gui components, but should still call this method.
     */
    virtual void setFixedWidthSlot(int);

    // if the property timeline is activated / deactivated by using the checkbox, this updates the core
    virtual void timelineActivated();
    virtual void deleteButtonClicked();

signals:

    void deleteTimelineSignal(PropertyTimeline*);

protected:

    QGridLayout* lineLayout_;

    QPushButton* deleteButton_;
    QPushButton* activateTimelineButton_;
    QLabel* nameLabel_;

    PropertyTimeline* propertyTimeline_;

    AnimationInfoWidgetBase* informationWidget_;    ///< the "parent" widget for updating keys etc.

//protected slots:

};

//----------------------------------------------------
// AnimationInfoWidgetBase
//----------------------------------------------------

class AnimationInfoWidgetBase : public QWidget {
    Q_OBJECT
public:
    AnimationInfoWidgetBase(QWidget* = 0);

    /// adds a property timeline (and the corresponding line widget) to the widget
    virtual void addPropertyTimeline(PropertyTimeline* prop) = 0;

    /// Remove all widgets and delete internal copies of properties
    virtual void clear() = 0;

    /// this enables / disables the gui components according to the activation status of the property timelines in the core
    virtual void updateActivationStatus() = 0;

    /// sets the current application mode configuration
    virtual void setApplicationModeConfig(ApplicationModeConfiguration* appConfig);

    /// helper method that sets the value of the given key to the property and updates the PropertyWidget from the property afterwards (if all types correspond)
    virtual void updatePropertyFromKey(Property*, PropertyWidget*, const PropertyKeyValueBase*, PropertyTimeline* tl);

    virtual void updateKeyFromWidget(PropertyTimeline* timeline, const PropertyKeyValueBase* key, Property* prop);

public slots:
    /// for setting the fixed width. Sadly there seems to be no easy way to do the correct resizing within qscrollareas
    virtual void setFixedWidthSlot(int);

protected:

    //helper method for sorting timelines by their priority within an application mode group
    static bool sortTimelinesByPriority(const std::pair<int, PropertyTimeline*>& a, const std::pair<int, PropertyTimeline*>& b);

    virtual QSize sizeHint() const;

    virtual void removeWidgetsFromLayout(QLayout*);

    virtual void createHeader() = 0;

    virtual void resizeHeader(int) = 0;

    QVBoxLayout* mainLayout_;

    QWidget* headerWidget_;
    QGridLayout* headerLayout_;

    /// associated application mode property groups with their corresponding group boxes in the gui
    std::map<std::string, QGroupBox*> applicationModeGroups_;

    /// associates application mode property groups with an ordered list of their property timelines
    std::map<std::string, std::vector<std::pair<int, PropertyTimeline*> > > groupTimelines_;

    ApplicationModeConfiguration* appConfig_;

signals:

    void deleteTimelineSignal(PropertyTimeline*);
    void setFixedWidthSignal(int);
};

} // namespace voreen

#endif

