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

#include "voreen/qt/widgets/property/transfunc/1d/1dgaussian/transfunc1dgaussianpropertywidget.h"

#include "voreen/qt/widgets/property/transfunc/1d/1dgaussian/transfunc1dgaussianpropertyeditor.h"

#include "voreen/core/datastructures/transfunc/1d/1dgaussian/transfunc1dgaussian.h"

#include <QPixmap>
#include <QPainter>
#include <QAction>

namespace voreen {

TransFunc1DGaussianPropertyWidget::TransFunc1DGaussianPropertyWidget(TransFunc1DGaussianProperty* prop, QWidget* parent)
    : TransFunc1DPropertyWidget(prop, parent)
{}

TransFunc1DGaussianPropertyWidget::~TransFunc1DGaussianPropertyWidget() {
}

//-----------------------------------------------------------------------------------
//    Override Functions
//-----------------------------------------------------------------------------------
QWidget* TransFunc1DGaussianPropertyWidget::createToolWindowWidget()  {
    editor_ = new TransFunc1DGaussianPropertyEditor(dynamic_cast<TransFunc1DGaussianProperty*>(property_));
    editor_->initialize();
    return editor_;
}

const QPixmap TransFunc1DGaussianPropertyWidget::createPresetPreviewIcon(std::string& presetPath) const {
    TransFunc1DGaussian preview;
    preview.load(presetPath);
    QPixmap icon(32,24);
    QPainter painter(&icon);
    for(int i = 0; i < 32; i++) {
        tgt::Color color = preview.getTexture()->texelAsFloat(static_cast<size_t>((static_cast<float>(i)/31.f)*(preview.getDimensions().x-1)));
        painter.setPen(QColor(color.r*255,color.g*255,color.b*255));
        painter.drawLine(i,0,i,23);
    }
    return icon;
}

void TransFunc1DGaussianPropertyWidget::loadPresetSlot(QAction* action) {
    QString data = action->data().toString();
    if(property_ && dynamic_cast<TransFunc1DGaussian*>(property_->get())) {
        TransFunc1DGaussian* tf1DGaussian = dynamic_cast<TransFunc1DGaussian*>(property_->get());
        tgt::vec2  oldDomain = tf1DGaussian->getDomain();
        TransFuncBase::AlphaMode alphaMode = property_->get()->getAlphaMode();
        if (property_->get()->load(data.toStdString())) {
            if(tf1DGaussian->getDomain() == tgt::vec2(0.0f, 1.0f)) {
                tf1DGaussian->setDomain(oldDomain);
            }
            property_->get()->setAlphaMode(alphaMode);
            property_->invalidate();
        }
        else {
            //TODO: error
        }
    }
}

} // namespace voreen
