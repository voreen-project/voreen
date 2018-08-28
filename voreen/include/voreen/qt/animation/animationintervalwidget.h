/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2018 University of Muenster, Germany,                        *
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

#ifndef ANIMATIONINTERVALWIDGET_H
#define ANIMATIONINTERVALWIDGET_H

#include "voreen/qt/animation/animationinfowidgetbase.h"

#include "voreen/core/animation/propertykeyvalue.h"
#include "voreen/core/animation/interpolationfunctionfactory.h"
#include "voreen/core/animation/templatepropertytimeline.h"

#include <QMenu>
#include <QPushButton>
#include <QVBoxLayout>
#include <QHBoxLayout>
#include <QWidget>

namespace voreen {

class AnimationIntervalWidget;

//----------------------------------------------------
// AnimationIntervalLineWidget
//----------------------------------------------------

class AnimationIntervalLineWidget : public AnimationLineWidgetBase {
    Q_OBJECT

public:
    AnimationIntervalLineWidget(AnimationIntervalWidget* parent, PropertyTimeline* tl);

    // deletes the property copy and clears all qt objects
    virtual ~AnimationIntervalLineWidget();

    void updateProperties(PropertyKeyValueBase* lKey, PropertyKeyValueBase* rKey);

    void updateKeyframes(PropertyKeyValueBase* lKey, PropertyKeyValueBase* rKey);

    /// this enables / disables the gui components according to the activation status of the property timelines in the core
    virtual void updateActivationStatus();

    /// checks if construction of the widget worked, ie. if all components are defined
    virtual bool isValid() const;

public slots:

    /// for setting the fixed width. necessary because of the use of QScrollArea in SimpleTimelineWidget
    virtual void setFixedWidthSlot(int);

protected:

    void populateInterpolationMenu();

    template<class T>
    void populateTemplateInterpolationMenu() {
        InterpolationFunctionFactory* iff = new InterpolationFunctionFactory();
        std::vector<InterpolationFunction<T>*> listOfFunctions = iff->getListOfFunctions<T>();
        typename std::vector<InterpolationFunction<T>*>::iterator it;
        std::map<std::string, QMenu*> categories;
        it = listOfFunctions.begin();
        while(it != listOfFunctions.end()) {
            std::string id = (*it)->getCategory();
            if (categories.find(id) == categories.end()) {
                QMenu* categoryMenu = new QMenu(QString::fromStdString(id), this);
                categories[id] = categoryMenu;
                interpolationMenu_->addMenu(categoryMenu);
            }
            it++;
        }

        it = listOfFunctions.begin();
        while(it != listOfFunctions.end()) {
            QMenu* categoryMenu = (*categories.find((*it)->getCategory())).second;
            QAction* action  = new QAction(QString::fromStdString((*it)->getGuiName()), this);
            categoryMenu->addAction(action);
            QActionInterpolationFunctionMap_[action] = (*it)->create();
            it++;
        }
        delete iff;

    }

    template<class T>
    void changeInterpolationTemplate(QAction* action) {
        TemplatePropertyTimeline<T>* tl = dynamic_cast<TemplatePropertyTimeline<T>*>(propertyTimeline_);
        PropertyKeyValue<T>* lK = dynamic_cast<PropertyKeyValue<T>*>(leftKey_);
        PropertyKeyValue<T>* rK = dynamic_cast<PropertyKeyValue<T>*>(rightKey_);
        InterpolationFunction<T>* func = dynamic_cast<InterpolationFunction<T>*>(QActionInterpolationFunctionMap_[action]);

        tl->setInterpolationFunctionAfter(func->create(), lK);
        tl->setInterpolationFunctionBefore(func->create(), rK);

        //TODO: what about the smooth stuff?! see templatepropertytimelinewidget, l. 482 - 509

        // update the interpolation menu button
        interpolationFunction_->setText(QString::fromStdString(func->getGuiName()));
    }

    //properties and widgets for left and right key
    Property* lProperty_;
    Property* rProperty_;
    PropertyWidget* lPropertyWidget_;
    PropertyWidget* rPropertyWidget_;

    QPushButton* interpolationFunction_;

    QMenu* interpolationMenu_;

    //remember the left and right key everytime the properties are updated by the keys
    PropertyKeyValueBase* leftKey_;
    PropertyKeyValueBase* rightKey_;

    /// links the string with the corresponding InterpolationFunction
    std::map<QAction*, InterpolationFunctionBase*> QActionInterpolationFunctionMap_;

protected slots:

    void updateLeftKeyframe();
    void updateRightKeyframe();

    void changeInterpolation(QAction*);

};


//----------------------------------------------------
// AnimationIntervalWidget
//----------------------------------------------------

class AnimationIntervalWidget : public AnimationInfoWidgetBase {
    Q_OBJECT
public:
    AnimationIntervalWidget(QWidget* = 0);

    /// creates an internal copy of the property associated with the PropertyTimeline and adds a corresponding property widget
    virtual void addPropertyTimeline(PropertyTimeline* prop);

    /// Remove all widgets and delete internal copies of properties
    virtual void clear();

    /// Updates the values of the properties (if a key has been selected)
    virtual void updateProperties(const std::vector<std::pair<PropertyKeyValueBase*, PropertyTimeline*> >& leftPropertyValues,
                            const std::vector<std::pair<PropertyKeyValueBase*, PropertyTimeline*> >& rightPropertyValues);

    /// this enables / disables the gui components according to the activation status of the property timelines in the core
    virtual void updateActivationStatus();

protected:

friend class AnimationIntervalLineWidget;

    virtual void createHeader();

    virtual void resizeHeader(int);

    std::map<PropertyTimeline*, AnimationIntervalLineWidget*> timelineWidgets_;
};

} // namespace voreen

#endif

