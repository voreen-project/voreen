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

#ifndef VRN_FONTPROPERTYWIDGET_H
#define VRN_FONTPROPERTYWIDGET_H

#include <QWidget>
#include <QComboBox>
#include <QGroupBox>
#include <QSlider>
#include <QLCDNumber>
#include <math.h>
#include <tgt/vector.h>
#include "voreen/qt/widgets/property/qpropertywidget.h"
#include "voreen/core/properties/fontproperty.h"
#include "voreen/core/properties/colorproperty.h"
#include "voreen/qt/widgets/property/colorpropertywidget.h"
#include "voreen/qt/widgets/sliderspinboxwidget.h"

namespace voreen {
class FontPropertyWidget : public QPropertyWidget {
Q_OBJECT
public:
    FontPropertyWidget(FontProperty* prop, QWidget* parent = 0);
public slots:
    void updateProperty();

protected:
    FontProperty* property_;

protected slots:
    virtual void updateFromPropertySlot();

private:
    QGroupBox* groupBox_;
    QComboBox* tgtFontName_;
    ColorProperty colorProperty_;
    ColorPropertyWidget* colorPropertyWidget_;
    SliderSpinBoxWidget* fontSizeSlider_;
};

} // namespace voreen

#endif // VRN_LIGHTPROPTERYWIDGET_H
