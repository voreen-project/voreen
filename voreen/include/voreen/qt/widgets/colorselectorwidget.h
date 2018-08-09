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

#ifndef VRN_COLORSELECTORWIDGET_H
#define VRN_COLORSELECTORWIDGET_H

#include "voreen/qt/widgets/property/qpropertywidget.h"

#include "tgt/vector.h"

#include <QLabel>
#include <QColor>

namespace voreen {

/**
 * Widget used for color selection, e.g., in the ColorPropertyWidget.
 */
class ColorSelectorWidget : public QLabel {
    Q_OBJECT
public:
    /** Constructor */
    ColorSelectorWidget(const QString& text, QWidget* parent=0, bool useAlphaChannel=true);

    /** Setter */
    void setColor(const QColor& color);
    /** Getter */
    QColor getColor() const;

signals:
    /** Signal emitted on color changed*/
    void colorChangedSignal(QColor);

protected:
    /**
     * Opens a color dialog.
     * @see mousePressEvent */
    void openColorDialog();

    //------------------
    //  Events
    //------------------
    /**
     * @override QLabel::mousePressEvent
     * Opens a color dialog on mouse click.
     */
    virtual void mousePressEvent(QMouseEvent* e) override;
    /**
     * @override QLabel::mousePressEvent
     * Opens a color dialog on mouse click.
     */
    virtual void paintEvent(QPaintEvent* event) override;

    //------------------
    //  Members
    //------------------
    QColor currentColor_;   ///< currently selected color
    bool useAlphaChannel_;  ///< determines, if the alpha channel should be selectable
};
} // namespace

#endif // VRN_COLORSELECTORWIDGET_H
