/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2021 University of Muenster, Germany,                        *
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

#ifndef VRN_GROUPPROPERTYWIDGET_H
#define VRN_GROUPPROPERTYWIDGET_H

#include "voreen/qt/widgets/property/qpropertywidget.h"
#include <QGroupBox>

class QTabWidget;
class QPushButton;

namespace voreen {

class GroupPropertyWidget : public QPropertyWidget {
    Q_OBJECT
public:
    GroupPropertyWidget(Property*, bool, std::string, QWidget* = 0);
    void addWidget(QPropertyWidget*, const QString);
    void addWidget(QPropertyWidget*, QWidget*, const QString);

    Property* getProperty();

    std::string getGroupName() const;

    /// Check if there is a visible property (used to hide grouping widget in user mode if necessary)
    bool isAnyPropertyVisible(Property::LevelOfDetail lod);

    virtual CustomLabel* getOrCreateNameLabel() const;
protected slots:
    virtual void updateFromPropertySlot();
    void hideGroup(bool);

private:
    QPushButton* hideWidgetButton_;
    QGridLayout* gridLayout_;
    std::vector<QPropertyWidget*> propertyWidgets_;
    QTabWidget* tabWidget_;
    QGroupBox* groupBox_;
    Property* prop_;
    bool tabbed_;
    std::string groupName_;

};

} // namespace voreen

#endif // VRN_GROUPPROPERTYWIDGET_H
