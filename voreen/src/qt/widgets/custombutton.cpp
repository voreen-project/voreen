/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2020 University of Muenster, Germany,                        *
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

#include "voreen/qt/widgets/custombutton.h"
#include "voreen/qt/widgets/property/qpropertywidget.h"

//#include "voreen/core/properties/transfuncproperty.h"

#include <QContextMenuEvent>
#include <QFont>
#include <QLineEdit>
#include <QMenu>
#include <QWidget>

namespace voreen {

CustomButton::CustomButton(const QString& text, QPropertyWidget* pw, bool editable, QWidget* parent)
    : QPushButton(text, parent)
    , propertyWidget_(pw)
    , editable_(editable)
{
    init();
}

CustomButton::~CustomButton() {
}

void CustomButton::init() {
    initFont();
    edit_ = new QLineEdit(this);
    edit_->hide();
    connect(edit_, SIGNAL(editingFinished()), this, SLOT(editingFinished()));
}

void CustomButton::editingFinished() {
    setText(edit_->text());

    if(propertyWidget_ != 0) {
        propertyWidget_->setPropertyGuiName(edit_->text().toStdString());
    }
    edit_->hide();
}

void CustomButton::contextMenuEvent(QContextMenuEvent* e) {
    if (editable_) {
        QMenu* men = new QMenu(this);
        men->addAction("Rename");
        QAction* ac = men->exec(e->globalPos());
        if(ac != 0){
            if(ac->iconText().compare("Rename") == 0){
                edit_->setText(text());
                edit_->setFocus();
                edit_->setCursorPosition(edit_->text().length());
                edit_->resize(size());
                edit_->show();
            }
        }
    }
}

void CustomButton::initFont() {
    QFontInfo fontInfo(font());
    setFont(QFont(fontInfo.family(), QPropertyWidget::fontSize_));
    setToolTip(text());
}

} // namespace
