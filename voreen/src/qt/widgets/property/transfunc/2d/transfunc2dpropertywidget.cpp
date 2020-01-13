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

#include "voreen/qt/widgets/property/transfunc/2d/transfunc2dpropertywidget.h"
#include "voreen/core/datastructures/transfunc/2d/transfunc2d.h"

#include "voreen/core/voreenapplication.h"
#include "voreen/core/utils/stringutils.h"
#include "voreen/core/datastructures/volume/histogram.h"

#include "tgt/filesystem.h"

#include <QPixmap>
#include <QMenu>

namespace voreen {

TransFunc2DPropertyWidget::TransFunc2DPropertyWidget(TransFuncPropertyBase* prop, QWidget* parent)
    : TransFuncPropertyWidgetBase(prop, parent)
{}

TransFunc2DPropertyWidget::~TransFunc2DPropertyWidget() {
}

//-----------------------------------------------------------------------------------
//    Override Functions
//-----------------------------------------------------------------------------------
TransFuncPropertyWidgetBase::WidgetBaseLayout TransFunc2DPropertyWidget::getBaseLayout() const {
    return WBL_2D;
}

TransFuncPropertyWidgetPainterBase* TransFunc2DPropertyWidget::createPreviewCanvasPainter(tgt::QtCanvas* canvas) const {
    return new TransFunc2DPropertyWidgetPainter(canvas, QWidget::palette().color(QWidget::backgroundRole()));
}

void TransFunc2DPropertyWidget::storeZoomMetaData() const {
    property_->getMetaDataContainer().addMetaData("TransFunc2DPropertyWidgetPainterZoom",
                                  new Vec4MetaData(tgt::vec4(dynamic_cast<TransFunc2DPropertyWidgetPainter*>(previewPainter_)->getMinVisibleDomainValue(0),
                                                             dynamic_cast<TransFunc2DPropertyWidgetPainter*>(previewPainter_)->getMaxVisibleDomainValue(0),
                                                             dynamic_cast<TransFunc2DPropertyWidgetPainter*>(previewPainter_)->getMinVisibleDomainValue(1),
                                                             dynamic_cast<TransFunc2DPropertyWidgetPainter*>(previewPainter_)->getMaxVisibleDomainValue(1))));
}

void TransFunc2DPropertyWidget::restoreZoomMetaData() const {
    if(property_->getMetaDataContainer().hasMetaData("TransFunc2DPropertyWidgetPainterZoom")) {
        tgt::vec4 range = static_cast<Vec4MetaData*>(property_->getMetaDataContainer().getMetaData("TransFunc2DPropertyWidgetPainterZoom"))->getValue();
        dynamic_cast<TransFunc2DPropertyWidgetPainter*>(previewPainter_)->setVisibleDomainValues(range.xy(),range.zw());
    }
}

Histogram* TransFunc2DPropertyWidget::getHistogramFromVolume(const VolumeBase* vb, const size_t channel) const {
    if(vb && vb->hasDerivedData<VolumeHistogramIntensityGradient >())
        return  &vb->getDerivedData<VolumeHistogramIntensityGradient>()->getHistogram(channel);
    else {
        // only auto-compute histogram if the option is activated
        if (property_->getComputeHistogram())
            vb->getDerivedDataThreaded<VolumeHistogramIntensityGradient>();
    }
    return 0;
}

void TransFunc2DPropertyWidget::createDomainMenuSlot() {
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

void TransFunc2DPropertyWidget::loadDomainSlot(QAction* action) {
    TransFuncPropertyWidgetBase::loadDomainSlot(action);
    QString data = action->data().toString();
    tgt::File* file = FileSys.open(data.toStdString());
    if(file) {
        std::string content = file->getAsString();
        std::vector<std::string> expl = strSplit(content, '\n');
        if(expl.size() == 4) {
             std::cout << expl[0] << expl[1];
             float width1 = stof(expl[0]);
             float center1 = stof(expl[1]);
             float width2 = stof(expl[2]);
             float center2 = stof(expl[3]);

             tgt::vec2 domain1 = tgt::vec2(center1 - (width1 * 0.5f), center1 + (width1 * 0.5f));
             tgt::vec2 domain2 = tgt::vec2(center2 - (width2 * 0.5f), center2 + (width2 * 0.5f));
             if(TransFunc2D* tf2D = dynamic_cast<TransFunc2D*>(property_->get())) {
                 tf2D->setDomain(domain1,0); tf2D->setDomain(domain2,1);
                 property_->invalidate();
              }
              else {
                 LERRORC("voreen.TransFuncPropertyWidget", "Failed to parse window preset file");
             }
        }
    }

}








//
//
//
//
//
//
//
//
//
//
//
//
//
//
//
//
//
//
//
//void TransFuncPropertyWidget::setVisible(bool visible) {
//    QPropertyWidgetWithEditorWindow::setVisible(visible);
//    textureCanvas_->update();
//}
//
//void TransFuncPropertyWidget::updateFromPropertySlot() {
//    if (plugin_)
//        plugin_->updateFromProperty();
//
//    const VolumeBase* vb = property_->getVolumeHandle();
//    if(property_->getMetaDataContainer().hasMetaData("TransfuncPropertyWidgetPainterZoom")) {
//        tgt::vec2 range = static_cast<Vec2MetaData*>(property_->getMetaDataContainer().getMetaData("TransfuncPropertyWidgetPainterZoom"))->getValue();
//        texturePainter_->setTransFunc(dynamic_cast<TransFunc1DKeys*>(property_->get()), range.x, range.y);
//    } else {
//        texturePainter_->setTransFunc(dynamic_cast<TransFunc1DKeys*>(property_->get()));
//    }
//    // We want the view to reset when the TF is deserialized.
//    // To achieve this we set the widget as not initialized as long as only the standard function has been set
//    const TransFuncBase* tf = property_->get();
//    /*if(!viewInitialized_) {
//        if(tf) {
//            const TransFunc1DKeys* tf1d = dynamic_cast<const TransFunc1DKeys*>(tf);
//            if(tf1d) {
//                if(!tf1d->isStandardFunc())
//                    viewInitialized_ = true;
//            }
//            else
//                viewInitialized_ = true;
//        }
//
//        if(viewInitialized_)
//            resetZoom();
//    }*/
//
//    clearObserveds();
//    if(vb) {
//        vb->addObserver(this);
//
//        if(vb->hasDerivedData<VolumeHistogramIntensity>())
//            texturePainter_->setHistogram(&(vb->getDerivedData<VolumeHistogramIntensity>()->getHistogram(property_->getVolumeChannel())));
//        else
//            vb->getDerivedDataThreaded<VolumeHistogramIntensity>();
//    }
//    else {
//        texturePainter_->setHistogram(0);
//    }
//
//    if(tf) {
//        alphaBt_->blockSignals(true);
//        updateAlpha(tf->getAlphaMode());
//        alphaBt_->blockSignals(false);
//    }
//
//    textureCanvas_->update();
//}
//
//void TransFuncPropertyWidget::showToolTipSlot(QPoint pos, QString tip) {
//    QToolTip::showText(mapToGlobal(pos),tip);
//}
//
//void TransFuncPropertyWidget::hideToolTipSlot() {
//    QToolTip::hideText();
//}
//
//
//void TransFuncPropertyWidget::setProperty() {
//    if (!disconnected_) {
//        // lazy instantiation of transfunc editor window
//        if (!editorWindow_) {
//            createEditorWindow(Qt::NoDockWidgetArea);
//            tgtAssert(editorWindow_, "Transfunc editor not instantiated");
//        }
//
//        if (editorWindow_->isVisible()) {
//            //close widget
//            editorWindow_->close();
//        }
//        else {
//            //open Widget
//            editorWindow_->showNormal();
//            plugin_->updateFromProperty();
//        }
//    }
//}
//
//void TransFuncPropertyWidget::disconnect() {
//    disconnected_ = true;
//    if (plugin_)
//        plugin_->disconnect();
//}
//
//QWidget* TransFuncPropertyWidget::createEditorWindowWidget() {
//    plugin_ = new TransFuncPluginXXX(property_, parentWidget(), Qt::Horizontal);
//    plugin_->createWidgets();
//    plugin_->createConnections();
//    //connect(plugin_, SIGNAL(transferFunctionChanged()), this, SIGNAL(modified()));
//
//    return plugin_;
//}
//
//void TransFuncPropertyWidget::customizeEditorWindow() {
//    editorWindow_->setAllowedAreas(Qt::NoDockWidgetArea);
//    editorWindow_->setFloating(true);
//}
//
//Property* TransFuncPropertyWidget::getProperty() {
//    return property_;
//}
//
//void TransFuncPropertyWidget::showNameLabel(bool visible) {
//    if (nameLabel_) {
//        if (visible)
//            nameLabel_->show();
//        else
//            nameLabel_->hide();
//    }
//}
//
//CustomLabel* TransFuncPropertyWidget::getOrCreateNameLabel() const {
//    return 0;
//}
//
//void TransFuncPropertyWidget::toggleInteractionMode(bool on) {
//    QualityMode.requestQualityMode((on? VoreenQualityMode::RQ_INTERACTIVE : VoreenQualityMode::RQ_DEFAULT), this);
//}
//
//void TransFuncPropertyWidget::invalidateProperty() {
//    property_->invalidate();
//}
//
//void TransFuncPropertyWidget::updateZoomMeta() {
//    property_->getMetaDataContainer().addMetaData("TransfuncPropertyWidgetPainterZoom",
//                                      new Vec2MetaData(tgt::vec2(texturePainter_->getMinDomainValue(),texturePainter_->getMaxDomainValue())));
//    property_->invalidate();
//}
//
//void TransFuncPropertyWidget::volumeDelete(const VolumeBase* source) {
//    texturePainter_->setHistogram(0);
//    textureCanvas_->update();
//}
//
//void TransFuncPropertyWidget::volumeChange(const VolumeBase* source) {
//    texturePainter_->setHistogram(0);
//    textureCanvas_->update();
//}
//
//void TransFuncPropertyWidget::derivedDataThreadFinished(const VolumeBase* source) {
//    if(source->hasDerivedData<VolumeHistogramIntensity>()) {
//        texturePainter_->setHistogram(&(static_cast<const VolumeHistogramIntensity*>(derivedData)->getHistogram(property_->getVolumeChannel())));
//        textureCanvas_->update();
//    }
//}


} // namespace voreen
