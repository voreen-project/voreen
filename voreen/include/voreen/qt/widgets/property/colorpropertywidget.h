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

#ifndef VRN_COLORPROPERTYWIDGET_H
#define VRN_COLORPROPERTYWIDGET_H

#include "voreen/core/properties/colorproperty.h"
#include "voreen/qt/widgets/colorselectorwidget.h"

#include "voreen/qt/widgets/property/qpropertywidget.h"

#include "tgt/vector.h"

namespace voreen {

/**
 * The property widget associated with the color property.
 */
class ColorPropertyWidget : public QPropertyWidget {
    Q_OBJECT
public:
    /** Constructor */
    ColorPropertyWidget(ColorProperty* prop, QWidget* parent = 0);

protected slots:
    /** @override QPropertyWidget::updateFromPropertySlot() */
    virtual void updateFromPropertySlot() override;

public slots:
    /** Triggered by colorSelectorWidget_ on a color change. */
    void colorChangedSlot(QColor color);

    //------------------
    //  Helper
    //------------------
protected:
    /** Transform a QColor into a tgt::Color. */
    tgt::Color toTgtColor(const QColor& color);
    /** Transform a tgt::Color into a QColor. */
    QColor toQColor(const tgt::Color& color);

    //------------------
    //  Members
    //------------------
private:
    ColorProperty* property_;                   ///< pointer to the associated property
    ColorSelectorWidget* colorSelectorWidget_;  ///< used for the color selection
};

} // namespace

#endif // VRN_COLORPROPERTYWIDGET_H
