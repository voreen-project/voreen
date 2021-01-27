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

#include "voreen/qt/widgets/property/volumeinfopropertywidget.h"

#include "voreen/core/datastructures/volume/volume.h"

#include "voreen/qt/widgets/volumeviewhelper.h"

#include <QLabel>
#include <QVBoxLayout>

namespace voreen {


const std::string VolumeInfoPropertyWidget::loggerCat_("voreen.qt.VolumeInfoPropertyWidget");

VolumeInfoPropertyWidget::VolumeInfoPropertyWidget(VolumeInfoProperty* volumeInfoProp, QWidget* parent)
    : QPropertyWidget(volumeInfoProp, parent, false)
{
    if (!volumeInfoProp) {
        tgtAssert(false, "No volume property");
        LERROR("No volume property");
        return;
    }

    QVBoxLayout* mainLayout = new QVBoxLayout();
    mainLayout->setContentsMargins(0, 2, 0, 0);

    layout_->addLayout(mainLayout);

    QHBoxLayout* previewLayout = new QHBoxLayout();
    previewLayout->setContentsMargins(2, 0, 2, 2);
    QGridLayout* infoLayout = new QGridLayout();
    infoLayout->setContentsMargins(2, 0, 2, 2);
    infoLayout->setColumnStretch(0, 1);
    infoLayout->setColumnStretch(1, 2);

    previewLabel_ = new QLabel(this);
    previewLabel_->setMaximumWidth(75);

    dimensionLabelCaption_ = new CustomLabel(this);
    spacingLabelCaption_ = new CustomLabel(this);
    memSizeLabelCaption_ = new CustomLabel(this);

    dimensionLabelCaption_->setText(" Dimensions:");
    spacingLabelCaption_->setText(" Spacing:");
    memSizeLabelCaption_->setText(" MemSize:");

    volumeNameLabel_ = new CustomLabel(this);
    pathLabel_ = new CustomLabel(this);
    dimensionLabel_ = new CustomLabel(this);
    spacingLabel_ = new CustomLabel(this);
    memSizeLabel_ = new CustomLabel(this);

    volumeNameLabel_->setTextInteractionFlags(Qt::TextSelectableByMouse);
    volumeNameLabel_->setAlignment(static_cast<Qt::Alignment>(Qt::AlignHCenter | Qt::AlignVCenter));
    pathLabel_->setTextInteractionFlags(Qt::TextSelectableByMouse);
    pathLabel_->setWordWrap(true);
    dimensionLabel_->setTextInteractionFlags(Qt::TextSelectableByMouse);
    dimensionLabel_->setWordWrap(true);
    spacingLabel_->setTextInteractionFlags(Qt::TextSelectableByMouse);
    spacingLabel_->setWordWrap(true);
    memSizeLabel_->setTextInteractionFlags(Qt::TextSelectableByMouse);
    memSizeLabel_->setWordWrap(true);

    QHBoxLayout* volumeLayout = new QHBoxLayout();
    volumeLayout->setContentsMargins(0, 0, 4, 0);
    volumeLayout->setSpacing(4);
    volumeLayout->setMargin(0);

    mainLayout->addLayout(volumeLayout);

    previewLayout->addWidget(previewLabel_);
    previewLayout->addLayout(infoLayout);
    infoLayout->addWidget(pathLabel_, 0, 0, 1, 2, 0);
    infoLayout->addWidget(dimensionLabelCaption_, 1, 0);
    infoLayout->addWidget(spacingLabelCaption_, 2, 0);
    infoLayout->addWidget(memSizeLabelCaption_, 3, 0);

    infoLayout->addWidget(dimensionLabel_, 1, 1);
    infoLayout->addWidget(spacingLabel_, 2, 1);
    infoLayout->addWidget(memSizeLabel_, 3, 1);

    QHBoxLayout* separatorLayout = new QHBoxLayout();
    QFrame* frame = new QFrame();
    frame->setFrameShape(QFrame::HLine);
    separatorLayout->addWidget(frame);
    separatorLayout->addWidget(volumeNameLabel_);
    frame = new QFrame();
    frame->setFrameShape(QFrame::HLine);
    separatorLayout->addWidget(frame);

    mainLayout->addLayout(separatorLayout);
    mainLayout->addLayout(previewLayout);

    updateFromProperty();
}

const VolumeBase* VolumeInfoPropertyWidget::getVolumeFromProperty() const {

    VolumeInfoProperty* infoProp = dynamic_cast<VolumeInfoProperty*>(prop_);
    if (!infoProp) {
        LWARNING("No volume property");
        return 0;
    }

    return infoProp->getVolume();
}

void VolumeInfoPropertyWidget::resizeEvent(QResizeEvent* e) {
    QPropertyWidget::resizeEvent(e);
    pathLabel_->setText(QString::fromStdString(clipText(pathLabel_,pathLabel_->toolTip())));
}

void VolumeInfoPropertyWidget::updateFromPropertySlot() {
    const VolumeBase* handle = getVolumeFromProperty();
    if (handle) {
        dimensionLabel_->show();
        pathLabel_->show();
        spacingLabel_->show();
        memSizeLabel_->show();
        previewLabel_->show();
        dimensionLabelCaption_->show();
        spacingLabelCaption_->show();
        memSizeLabelCaption_->show();

        std::string name = VolumeViewHelper::getStrippedVolumeName(handle);
        std::string path = VolumeViewHelper::getVolumePath(handle);
        if(name.size() > 30) {
            volumeNameLabel_->setToolTip(QString::fromStdString(name));
            int end = static_cast<int>(name.size());
            std::string startString;
            std::string endString;
            for(size_t i = 0; i < 14; i++){
                 startString += name.at(i);
                 endString += name.at(end-14+i);
            }
            name = startString+"..."+endString;
        }
        pathLabel_->setToolTip(QString::fromStdString(path));
        path = clipText(pathLabel_,QString::fromStdString(path));

        volumeNameLabel_->setText(QString::fromStdString(" " + name + " (" + handle->getFormat() + ") "));

        pathLabel_->setText(QString::fromStdString(" "+path));
        dimensionLabel_->setText(QString::fromStdString(VolumeViewHelper::getVolumeDimension(handle)));
        spacingLabel_->setText(QString::fromStdString(VolumeViewHelper::getVolumeSpacing(handle)));
        memSizeLabel_->setText(QString::fromStdString(VolumeViewHelper::getVolumeMemorySize(handle)));
        previewLabel_->setPixmap(VolumeViewHelper::generateBorderedPreview(handle, 70, 0));

    }
    else {
        volumeNameLabel_->setText(tr(" no volume"));
        volumeNameLabel_->adjustSize();

        pathLabel_->hide();
        previewLabel_->setPixmap(QPixmap());
        dimensionLabel_->hide();
        spacingLabel_->hide();
        memSizeLabel_->hide();
        previewLabel_->hide();
        dimensionLabelCaption_->hide();
        spacingLabelCaption_->hide();
        memSizeLabelCaption_->hide();
    }
}

std::string VolumeInfoPropertyWidget::clipText(CustomLabel* l, QString fullText) {
    if(fullText.isEmpty()) return fullText.toStdString();

    int textWidth = l->fontMetrics().boundingRect(fullText).width();
    int labelWidth = l->width();

    if (textWidth+10 > labelWidth) {
        int estimatedPixelPerLetter = textWidth/fullText.size() + 1;
        int lettersToRemove =(textWidth+10-labelWidth)/ estimatedPixelPerLetter + 3;
        int count = (fullText.size() - lettersToRemove) / 2;

        int end = static_cast<int>(fullText.size());
        QString startString;
        QString endString;

        for(int i = 0; i < count; i++){
             startString += fullText.at(i);
             endString += fullText.at(end-count+i);
        }
        fullText = startString+"..."+endString;
    }

    return fullText.toStdString();
}

void VolumeInfoPropertyWidget::showNameLabel(bool) {
    if (nameLabel_)
        nameLabel_->hide();
}

CustomLabel* VolumeInfoPropertyWidget::getOrCreateNameLabel() const {
    return 0;
}

} //namespace voreen
