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

#ifndef VRN_COMPACTPROPERTYVECTORWIDGET_H
#define VRN_COMPACTPROPERTYVECTORWIDGET_H

#include "voreen/qt/widgets/property/qpropertywidget.h"
#include "voreen/core/properties/propertyvector.h"

class QGridLayout;

namespace voreen {

class PropertyVector;
class QPropertyWidgetFactory;

class PropertyVectorWidget : public QPropertyWidget {
Q_OBJECT
public:

    enum Properties {
        BOOL,
        BUTTON,
        FLOATMAT2,
        FLOATMAT3,
        FLOATMAT4,
        FLOATVEC2,
        FLOATVEC3,
        FLOATVEC4,
        INTVEC2,
        INTVEC3,
        INTVEC4,
        /*CAMERA,*/
        COLOR,
        /*FILEDIALOG,*/
        FLOAT,
        INT,
        LIGHT,
        OPTION,
        /*SHADER,*/
        STRING
    };

    PropertyVectorWidget(PropertyVector* prop, QWidget* parent = 0);

public slots:
    void setProperty(PropertyVector* value);

protected:

    /**
     * Creates the property widget for the passed property using the passed factory and
     * adds it to the property layout
     */
    void createAndAddPropertyWidget(Property* prop);
    void createAndAddPropertyWidgetByAction(QAction*);

    void initializePropertyMenu();

    /// The PropertyVector assigned to this widget
    PropertyVector* property_;

    /// Contains the property widgets.
    QGridLayout* propertiesLayout_;

    virtual void contextMenuEvent (QContextMenuEvent*);

    void propertyAdded();

    std::map<QAction*, int> propertyMap_;

    QMenu* propertyMenu_;

protected slots:
    virtual void updateFromPropertySlot();

};

} // namespace

#endif // VRN_COMPACTBOOLPROPERTYWIDGET_H
