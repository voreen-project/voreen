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

#ifndef CUSTOMLABEL_H
#define CUSTOMLABEL_H

#include "voreen/qt/widgets/property/qpropertywidget.h"

#include <QLabel>
#include <QWidget>

class QLineEdit;

namespace voreen {

class CustomLabel : public QLabel {
Q_OBJECT
public:
    CustomLabel(const char* text, QPropertyWidget* pw = 0, QWidget* parent = 0,
                Qt::WindowFlags f = 0, bool disallowWordwrap = false, bool editable = false);
    CustomLabel(QPropertyWidget* pw = 0, QWidget* parent = 0, Qt::WindowFlags f = 0, bool editable = false);
    ~CustomLabel();

    void init();

protected slots:
    void editingFinished();

protected:
    void initFont();
    void processText();
    virtual void contextMenuEvent (QContextMenuEvent*);

    QPropertyWidget* propertyWidget_;
    QLineEdit* edit_;
    bool disallowWordwrap_;
    bool editable_;
};

} // namespace

#endif // CUSTOMLABEL_H
