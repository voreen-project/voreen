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

#include "voreen/qt/widgets/property/temppathpropertywidget.h"

#include "voreen/core/properties/temppathproperty.h"

#include <QCheckBox>
#include <QPushButton>

#include "tgt/filesystem.h"

namespace voreen {

TempPathPropertyWidget::TempPathPropertyWidget(TempPathProperty* prop, QWidget* parent)
    : FileDialogPropertyWidget(prop, parent)
    , property_(prop)
    , useGeneratedPathCheckBox_(nullptr)
{
    openFileDialogBtn_->setDisabled(prop->getUseGeneratedPath());
    useGeneratedPathCheckBox_ = new QCheckBox(tr("Generate and Use Temporary Path"), this);
    useGeneratedPathCheckBox_->setChecked(prop->getUseGeneratedPath());
    connect(useGeneratedPathCheckBox_, SIGNAL(clicked(bool)), this, SLOT(changeUseGeneratedPathOption(bool)));

    mainLayout_->addWidget(useGeneratedPathCheckBox_);
}

void TempPathPropertyWidget::changeUseGeneratedPathOption(bool checked) {
    property_->setUseGeneratedPath(checked);
    openFileDialogBtn_->setDisabled(checked);
}

void TempPathPropertyWidget::updateFromPropertySlot() {
    updateButtonText(property_->get());
    useGeneratedPathCheckBox_->setChecked(property_->getUseGeneratedPath());
}

} // namespace
