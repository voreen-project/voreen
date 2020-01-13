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

#include "voreen/qt/widgets/property/stringpropertywidget.h"
#include "voreen/qt/widgets/customlabel.h"

#include <QLineEdit>

namespace voreen {

StringPropertyWidget::StringPropertyWidget(StringProperty* prop, QWidget* parent)
    : QPropertyWidget(prop, parent)
    , property_(prop)
{
    QString st = QString::fromUtf8(property_->get().c_str());
    lineEdit_ = new QLineEdit(st);
    addWidget(lineEdit_);
    lineEdit_->setReadOnly(!property_->isEditable() || property_->isReadOnlyFlagSet());

    connect(lineEdit_, SIGNAL(textChanged(QString)), this, SLOT(textHasChanged(QString)));
    connect(lineEdit_, SIGNAL(editingFinished()), this, SLOT(hasFinishedEditing()));
    connect(lineEdit_, SIGNAL(textChanged(QString)), this, SIGNAL(widgetChanged()));
}

void StringPropertyWidget::updateViewFlags(Property::ViewFlags flags) {
    // StringPropertyWidget is a special case and doesn't gray out the text field, since
    // this reduces visibility and readability of contained text.
    if (nameLabel_)
        nameLabel_->setEnabled(!(flags & Property::ViewFlags::VF_READ_ONLY));
    emit checkVisibility();
}

CustomLabel* StringPropertyWidget::getOrCreateNameLabel() const {
    QPropertyWidget::getOrCreateNameLabel();
    nameLabel_->setEnabled(!prop_->isReadOnlyFlagSet());
    return nameLabel_;
}

void StringPropertyWidget::updateFromPropertySlot() {
    lineEdit_->blockSignals(true);
    //QString st = QString::fromStdString(property_->get());
    QString st = QString::fromUtf8(property_->get().c_str());
    if(st != lineEdit_->text())
        lineEdit_->setText(st);
    lineEdit_->blockSignals(false);

    lineEdit_->setReadOnly(!property_->isEditable() || property_->isReadOnlyFlagSet());
}

void StringPropertyWidget::setProperty(const QString& text) {
    if (!disconnected_) {
        QByteArray bytes = text.toUtf8();
        std::string str(bytes.data(), bytes.length());

        property_->set(str);
        //property_->set(text.toStdString());
        emit valueModifiedByUser();
    }
    else
        updateFromProperty();
}

void StringPropertyWidget::textHasChanged(const QString& text)
{
    if (property_->getInstantUpdate()){
        setProperty(text);
    }
}

void StringPropertyWidget::hasFinishedEditing()
{
    setProperty(lineEdit_->text());
}

} // namespace
