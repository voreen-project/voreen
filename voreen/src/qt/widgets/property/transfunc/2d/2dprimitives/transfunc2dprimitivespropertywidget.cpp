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

#include "voreen/qt/widgets/property/transfunc/2d/2dprimitives/transfunc2dprimitivespropertywidget.h"

#include "voreen/qt/widgets/property/transfunc/2d/2dprimitives/transfunc2dprimitivespropertyeditor.h"
#include "voreen/core/datastructures/transfunc/2d/2dprimitives/transfunc2dprimitives.h"

#include <QPixmap>
#include <QPainter>
#include <QAction>

namespace voreen {

TransFunc2DPrimitivesPropertyWidget::TransFunc2DPrimitivesPropertyWidget(TransFunc2DPrimitivesProperty* prop, QWidget* parent)
    : TransFunc2DPropertyWidget(prop, parent)
{}

TransFunc2DPrimitivesPropertyWidget::~TransFunc2DPrimitivesPropertyWidget() {
}

//-----------------------------------------------------------------------------------
//    Override Functions
//-----------------------------------------------------------------------------------
QWidget* TransFunc2DPrimitivesPropertyWidget::createToolWindowWidget()  {
    editor_ = new TransFunc2DPrimitivesPropertyEditor(dynamic_cast<TransFunc2DPrimitivesProperty*>(property_));
    editor_->initialize();
    return editor_;
}

const QPixmap TransFunc2DPrimitivesPropertyWidget::createPresetPreviewIcon(std::string& presetPath) const {
    TransFunc2DPrimitives preview;
    preview.load(presetPath);
    QPixmap icon(32,32);
    QPainter painter(&icon);
    for(int x = 0; x < 32; x++) {
        for(int y = 0; y < 32; y++) {
            tgt::Color color = preview.getTexture()->texelAsFloat(tgt::svec2((tgt::vec2(x,y)/tgt::vec2(31.f))*tgt::vec2(preview.getDimensions().xy() - tgt::ivec2::one)));
            painter.setPen(QColor(color.r*255,color.g*255,color.b*255));
            painter.drawLine(x,y,x,y);
        }
    }
    return icon;
}

void TransFunc2DPrimitivesPropertyWidget::loadPresetSlot(QAction* action) {
    QString data = action->data().toString();
    if(property_ && dynamic_cast<TransFunc2DPrimitives*>(property_->get())) {
        TransFunc2DPrimitives* tf2DPrimitives = dynamic_cast<TransFunc2DPrimitives*>(property_->get());
        tgt::vec2  oldDomainX = tf2DPrimitives->getDomain(0);
        tgt::vec2  oldDomainY = tf2DPrimitives->getDomain(1);
        TransFuncBase::AlphaMode alphaMode = property_->get()->getAlphaMode();
        if (property_->get()->load(data.toStdString())) {
            if(tf2DPrimitives->getDomain(0) == tgt::vec2(0.0f, 1.0f)) {
                tf2DPrimitives->setDomain(oldDomainX,0);
                tf2DPrimitives->setDomain(oldDomainY,1);
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
