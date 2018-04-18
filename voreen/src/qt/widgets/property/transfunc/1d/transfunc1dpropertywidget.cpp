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

#include "voreen/qt/widgets/property/transfunc/1d/transfunc1dpropertywidget.h"
#include "voreen/core/datastructures/transfunc/1d/transfunc1d.h"

#include "voreen/core/voreenapplication.h"
#include "voreen/core/utils/stringutils.h"
#include "voreen/core/datastructures/meta/templatemetadata.h"
#include "voreen/core/datastructures/volume/histogram.h"

#include "tgt/vector.h"
#include "tgt/qt/qtcanvas.h"

#include <QPixmap>
#include <QMenu>

namespace voreen {

TransFunc1DPropertyWidget::TransFunc1DPropertyWidget(TransFuncPropertyBase* prop, QWidget* parent)
    : TransFuncPropertyWidgetBase(prop, parent)
{}

TransFunc1DPropertyWidget::~TransFunc1DPropertyWidget() {
}

//-----------------------------------------------------------------------------------
//    Override Functions
//-----------------------------------------------------------------------------------
TransFuncPropertyWidgetBase::WidgetBaseLayout TransFunc1DPropertyWidget::getBaseLayout() const {
    return WBL_1D;
}

TransFuncPropertyWidgetPainterBase* TransFunc1DPropertyWidget::createPreviewCanvasPainter(tgt::QtCanvas* canvas) const {
    return new TransFunc1DPropertyWidgetPainter(canvas);
}

void TransFunc1DPropertyWidget::storeZoomMetaData() const {
    property_->getMetaDataContainer().addMetaData("TransFunc1DPropertyWidgetPainterZoom",
                                  new Vec2MetaData(tgt::vec2(dynamic_cast<TransFunc1DPropertyWidgetPainter*>(previewPainter_)->getMinVisibleDomainValue(),
                                                             dynamic_cast<TransFunc1DPropertyWidgetPainter*>(previewPainter_)->getMaxVisibleDomainValue())));
}

void TransFunc1DPropertyWidget::restoreZoomMetaData() const {
    if(property_->getMetaDataContainer().hasMetaData("TransFunc1DPropertyWidgetPainterZoom")) {
        tgt::vec2 range = static_cast<Vec2MetaData*>(property_->getMetaDataContainer().getMetaData("TransFunc1DPropertyWidgetPainterZoom"))->getValue();
        dynamic_cast<TransFunc1DPropertyWidgetPainter*>(previewPainter_)->setVisibleDomainValues(range);
    }
}


Histogram* TransFunc1DPropertyWidget::getHistogramFromVolume(const VolumeBase* vb, const size_t channel) const {
    if(vb && vb->hasDerivedData<VolumeHistogramIntensity>())
        return  &vb->getDerivedData<VolumeHistogramIntensity>()->getHistogram(channel);
    else {
        // only compute the histogram if this functionality is activated
        if (property_->getComputeHistogram())
            vb->getDerivedDataThreaded<VolumeHistogramIntensity>();
    }
    return 0;
}

void TransFunc1DPropertyWidget::createDomainMenuSlot() {
    TransFuncPropertyWidgetBase::createDomainMenuSlot();

    domainMenu_->addSeparator();

    std::string applicationPresetPath = VoreenApplication::app()->getCoreResourcePath("transferfuncs/windows/" + getPathSubFolder());

    std::vector<std::string> subdirs = FileSys.listSubDirectories(applicationPresetPath, true);
    for(size_t i=0; i<subdirs.size(); i++) {
        QMenu* subMenu = domainMenu_->addMenu(QString::fromStdString(subdirs[i]));

        std::vector<std::string> presets = FileSys.listFiles(applicationPresetPath + "/" + subdirs[i], true);
        for(size_t j=0; j<presets.size(); i++) {
            QAction* propAction = new QAction(QString::fromStdString(FileSys.baseName(presets[j])), subMenu);
            propAction->setData(QVariant(QString::fromStdString(applicationPresetPath + "/" + subdirs[i] + "/" + presets[j])));
            subMenu->addAction(propAction);
        }
    }
}

void TransFunc1DPropertyWidget::loadDomainSlot(QAction* action) {
    TransFuncPropertyWidgetBase::loadDomainSlot(action);
    QString data = action->data().toString();
    tgt::File* file = FileSys.open(data.toStdString());
    if(file) {
        std::string content = file->getAsString();
        std::vector<std::string> expl = strSplit(content, '\n');
        if(expl.size() == 2) {
             std::cout << expl[0] << expl[1];
             float width = stof(expl[0]);
             float center = stof(expl[1]);

             tgt::vec2 domain = tgt::vec2(center - (width * 0.5f), center + (width * 0.5f));
             if(TransFunc1D* tf1D = dynamic_cast<TransFunc1D*>(property_->get())) {
                 tf1D->setDomain(domain);
                 property_->invalidate();
              }
              else {
                 LERRORC("voreen.TransFuncPropertyWidget", "Failed to parse window preset file");
             }
        }
    }

}

} // namespace voreen
