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

#ifndef VRN_QPROPERTYWIDGET_H
#define VRN_QPROPERTYWIDGET_H

#include "voreen/core/properties/property.h"
#include "voreen/core/properties/propertywidget.h"

#include "voreen/qt/voreenqtapi.h"

#include <QWidget>
#include <QBoxLayout>
#include <QToolButton>

class QLabel;

namespace voreen {

class CustomLabel;

class VRN_QT_API QPropertyWidget : public QWidget, public PropertyWidget {
    Q_OBJECT;

    friend class CustomLabel;

public:
    QPropertyWidget(Property* prop, QWidget* parent = 0, bool showNameLabel = true);
    virtual ~QPropertyWidget();
    virtual QSize sizeHint() const;

    /**
     * Public method called by the owning property. Delegates the call to updateFromPropertySlot()
     * via the signal/slot mechanism. This redirection is necessary to allow background threads
     * to call updateFromProperty().

     * @note: Do not overwrite this method, but updateFromPropertySlot() instead!
     */
    virtual void updateFromProperty();

    virtual void setEnabled(bool enabled);
    virtual void setVisible(bool state);
    virtual void updateViewFlags(Property::ViewFlags flags);

    virtual void disconnect();

    std::string getPropertyGuiName();
    virtual void setPropertyGuiName(std::string);
    virtual CustomLabel* getOrCreateNameLabel() const;

    // this is a static variable for the font size used in all propertywidgets
    static const int fontSize_;

public slots:
    /** Used to toggle the interaction mode. @see TransFuncPropertyWidget(Painter) */
    virtual void toggleInteractionMode(bool im);
    /** Used to invalidate the associated property. */
    virtual void invalidateProperty();
    //virtual void showNameLabel(bool); //remove?

signals:
    void valueModifiedByUser(); //TODO: deprecated. Used by timelines and should be removed
    void mouseClicked();
    void widgetChanged();

    void checkGroupVisibility();
    void checkVisibility();

protected:
    void addWidget(QWidget* widget);
    void addLayout(QLayout* layout);
    void mouseMoveEvent(QMouseEvent*);

    bool disconnected_;
    QBoxLayout* layout_;

    mutable CustomLabel* nameLabel_;
    bool showNameLabel_;

protected slots:
    /**
     * Override this method to perform the update operations
     * that would usually be placed in updateFromProperty().
     *
     * updateFromProperty() calls are redirected to this method via a signal.
     */
    virtual void updateFromPropertySlot() = 0;

signals:
    void updateFromPropertySignal();
};

} // namespace

#endif // VRN_QPROPERTYWIDGET_H
