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

#ifndef VRN_COLORMAPPROPERTYWIDGET_H
#define VRN_COLORMAPPROPERTYWIDGET_H

#include "../../datastructures/colormap.h"
#include "voreen/qt/widgets/property/qpropertywidgetwithtoolwindow.h"
#include "tgt/vector.h"

#include <QLabel>
#include <QComboBox>

namespace voreen {

class ColorMapProperty;

class ClickableColorMapLabel : public QLabel {
    Q_OBJECT
public:
    ClickableColorMapLabel(const QString& text, QWidget* parent=0, Qt::WindowFlags f=0);
    void setColorMap(const ColorMap& cm);

signals:
    void clicked();

protected:
    void mousePressEvent(QMouseEvent* e);
    void paintEvent(QPaintEvent* event);
    tgt::Color toTgtColor(QColor color);
    QColor toQColor(tgt::Color color);

    ColorMap cm_;
};



// --------------------------------------------------------------------------------------



class ColorMapPropertyWidget : public QPropertyWidget {
    Q_OBJECT
public:
    ColorMapPropertyWidget(ColorMapProperty* prop, QWidget* parent = 0);

protected slots:
    virtual void updateFromPropertySlot();

private slots:
    void selectColormap(int index);

protected:

    void updateColorLabel();
    tgt::Color toTgtColor(QColor color);
    QColor toQColor(tgt::Color color);

private:
    ColorMapProperty* property_;
    QComboBox* cbMaps_;
    ClickableColorMapLabel* colorLbl_;
};

} // namespace

#endif // VRN_COLORMAPPROPERTYWIDGET_H
