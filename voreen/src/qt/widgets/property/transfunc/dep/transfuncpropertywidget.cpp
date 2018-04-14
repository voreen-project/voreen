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

#include "voreen/qt/widgets/property/transfuncpropertywidget.h"

#include "voreen/qt/widgets/voreentoolwindow.h"
#include "voreen/qt/widgets/transfunc/transfuncpluginxxx.h"
#include "voreen/qt/widgets/transfunc/utils/transfunciohelperqt.h"

#include "voreen/core/processors/processor.h"
#include "voreen/core/voreenapplication.h"
#include "voreen/core/datastructures/volume/volume.h"
#include "voreen/core/datastructures/volume/histogram.h"
#include "voreen/core/properties/transfuncproperty.h"
#include "voreen/core/datastructures/transfunc/1d/transfunc1d.h"
#include "voreen/core/datastructures/transfunc/1d/1dkeys/transfunc1dkeys.h"
#include "voreen/core/utils/voreenqualitymode.h"
#include "voreen/qt/widgets/transfunc/1d/transfuncpropertywidgetpainter.h"
#include "voreen/qt/widgets/customlabel.h"

#include <QToolButton>
#include <QMenu>
#include <QToolTip>
#include "tgt/qt/qtcanvas.h"
#include "tgt/filesystem.h"

namespace voreen {

TransFuncPropertyWidget::TransFuncPropertyWidget(TransFuncProperty* prop, QWidget* parent)
    : QPropertyWidgetWithEditorWindow(prop, parent, false)
    , plugin_(0)
    , property_(prop)
    , texturePainter_(0)
    , ignoreSlideUpdates_(false)
    , viewInitialized_(false)
    , advancedBt_(new QToolButton())
    , alphaBt_(new QToolButton())
    , tfBt_(new QToolButton())
    , windowBt_(new QToolButton())
    , zoomInBt_(new QToolButton())
    , zoomOutBt_(new QToolButton())
{
    //main layout
    QVBoxLayout* mainLayout = new QVBoxLayout();
    mainLayout->setAlignment(Qt::AlignLeft);
    mainLayout->setContentsMargins(0, 2, 0, 0);
    layout_->addLayout(mainLayout);
    //button layout
    QGridLayout* buttonLayout = new QGridLayout();
    buttonLayout->setSpacing(4);
    buttonLayout->setMargin(0);
    buttonLayout->setColumnStretch(0,4);
    buttonLayout->setColumnStretch(1,1);
    buttonLayout->setColumnStretch(2,1);
    buttonLayout->setColumnStretch(3,1);
    buttonLayout->setColumnStretch(4,1);
    buttonLayout->setColumnStretch(5,1);
    buttonLayout->setColumnStretch(6,1);
    buttonLayout->setColumnStretch(7,1);
    buttonLayout->setColumnStretch(8,1);
    mainLayout->addLayout(buttonLayout);
    //texture canvas
    textureCanvas_ = new tgt::QtCanvas("", tgt::ivec2(10, 10), tgt::GLCanvas::RGBADD, 0, true);
    texturePainter_ = new TransFuncPropertyWidgetPainter(textureCanvas_); //sets the painter to the canvas
    textureCanvas_->init();
    if(property_->getMetaDataContainer().hasMetaData("TransfuncPropertyWidgetPainterZoom")) {
        tgt::vec2 range = static_cast<Vec2MetaData*>(property_->getMetaDataContainer().getMetaData("TransfuncPropertyWidgetPainterZoom"))->getValue();
        texturePainter_->setTransFunc(dynamic_cast<TransFunc1D*>(property_->get()), range.x, range.y);
    } else {
        texturePainter_->setTransFunc(dynamic_cast<TransFunc1D*>(property_->get()));
    }
    textureCanvas_->setFixedHeight(30);
    textureCanvas_->setSizePolicy(QSizePolicy::MinimumExpanding, QSizePolicy::Fixed);
    textureCanvas_->getEventHandler()->addListenerToBack(texturePainter_);
    mainLayout->addWidget(textureCanvas_);
    connect(texturePainter_,SIGNAL(changedGamma()),this,SLOT(invalidateProperty()));
    connect(texturePainter_,SIGNAL(changedDomain()),this,SLOT(updateZoomMeta()));
    connect(texturePainter_,SIGNAL(interaction(bool)),this,SLOT(toggleInteractionMode(bool)));
    connect(texturePainter_,SIGNAL(showInfoToolTip(QPoint, QString)),this,SLOT(showToolTipSlot(QPoint, QString)));
    connect(texturePainter_,SIGNAL(hideInfoToolTip()),this,SLOT(hideToolTipSlot()));

     //tf button
    loadTFMenu_ = new QMenu();
    tfBt_->setIcon(QPixmap(":/qt/icons/colorize.png"));
    tfBt_->setToolTip("Color Map Presets");
    tfBt_->setMenu(loadTFMenu_);
    tfBt_->setPopupMode(QToolButton::InstantPopup);
    connect(loadTFMenu_, SIGNAL(aboutToShow()), this, SLOT(populateLoadTFMenu()));
    connect(loadTFMenu_, SIGNAL(triggered(QAction*)), this, SLOT(loadTF(QAction*)));
    //fitting button
    loadWindowMenu_ = new QMenu();
    windowBt_->setIcon(QPixmap(":/qt/icons/histogram_fit.png"));
    windowBt_->setToolTip("Window Fitting/Presets");
    windowBt_->setMenu(loadWindowMenu_);
    windowBt_->setPopupMode(QToolButton::InstantPopup);
    connect(loadWindowMenu_, SIGNAL(aboutToShow()), this, SLOT(populateLoadWindowMenu()));
    connect(loadWindowMenu_, SIGNAL(triggered(QAction*)), this, SLOT(loadWindow(QAction*)));
    //alpha button
    alphaMenu_ = new QMenu();
    alphaBt_->setIcon(QPixmap(":/qt/icons/alpha.png"));
    alphaBt_->setToolTip("Malipulate Alpha Channel");
    alphaBt_->setMenu(alphaMenu_);
    alphaBt_->setPopupMode(QToolButton::InstantPopup);
    connect(alphaMenu_, SIGNAL(aboutToShow()), this, SLOT(populateAlphaMenu()));
    connect(alphaMenu_, SIGNAL(triggered(QAction*)), this, SLOT(setAlpha(QAction*)));
    tgtAssert(prop, "No property passed");
    if(prop && prop->get())
        updateAlpha(prop->get()->getAlphaMode());

    //zoom buttons
    zoomInBt_->setIcon(QPixmap(":/qt/icons/viewmag+.png"));
    zoomInBt_->setToolTip("Zoom In");
    zoomOutBt_->setIcon(QPixmap(":/qt/icons/viewmag_.png"));
    zoomOutBt_->setToolTip("Zoom Out");
    connect(zoomInBt_, SIGNAL(clicked()), this, SLOT(zoomIn()));
    connect(zoomOutBt_, SIGNAL(clicked()), this, SLOT(zoomOut()));

    //advanced button
    advancedMenu_ = new QMenu();
    advancedBt_->setIcon(QPixmap(":/qt/icons/mapping-function.png"));
    advancedBt_->setToolTip("Advanced Settings");
    advancedBt_->setMenu(advancedMenu_);
    advancedBt_->setPopupMode(QToolButton::InstantPopup);
    connect(advancedMenu_, SIGNAL(aboutToShow()), this, SLOT(populateAdvancedMenu()));
    connect(advancedMenu_, SIGNAL(triggered(QAction*)), this, SLOT(doAdvancedAction(QAction*)));
    if (!prop->getLazyEditorInstantiation() || editorVisibleOnStartup())
        createEditorWindow(Qt::NoDockWidgetArea);

    //set button layout
    buttonLayout->addWidget(QPropertyWidgetWithEditorWindow::getOrCreateNameLabel(),0,0,1,1);

    buttonLayout->addWidget(tfBt_,0,1,1,1);
    buttonLayout->addWidget(windowBt_,0,2,1,1);
    buttonLayout->addWidget(alphaBt_,0,3,1,1);

    buttonLayout->addWidget(zoomInBt_,0,5,1,1);
    buttonLayout->addWidget(zoomOutBt_,0,6,1,1);

    buttonLayout->addWidget(advancedBt_,0,8,1,1);

    //stuff
    QFontInfo fontInfo(font());
    advancedBt_->setFont(QFont(fontInfo.family(), QPropertyWidget::fontSize_));

    updateFromProperty();
}

//-----------------------------------------------------------------------------------
//    tf menu
//-----------------------------------------------------------------------------------
void TransFuncPropertyWidget::populateLoadTFMenu() {
    loadTFMenu_->clear();

    std::string applicationPresetPath = VoreenApplication::app()->getCoreResourcePath("transferfuncs/presets");
    populateTFMenuFromDirectory(loadTFMenu_, applicationPresetPath);

    loadTFMenu_->addSeparator();

    std::string userPresetPath = VoreenApplication::app()->getUserDataPath("transferfuncs/presets");
    populateTFMenuFromDirectory(loadTFMenu_, userPresetPath);
}

void TransFuncPropertyWidget::populateTFMenuFromDirectory(QMenu* menu, std::string directory) {
    int currentSubSection = 0;
    std::vector<std::string> presets = FileSys.listFiles(directory, true);
    for(size_t i=0; i<presets.size(); i++) {
        //add seperator after sub section
        if((std::isdigit(FileSys.baseName(presets[i])[0])) && (FileSys.baseName(presets[i])[0] != currentSubSection)) {
            menu->addSeparator();
            currentSubSection = FileSys.baseName(presets[i])[0];
        }
        //create action and cut off prefix
        QAction* propAction;
        if(std::isdigit(FileSys.baseName(presets[i])[0]))
            propAction = new QAction(QString::fromStdString(FileSys.baseName(presets[i]).substr(4)), menu);
        else
            propAction = new QAction(QString::fromStdString(FileSys.baseName(presets[i])), menu);
        propAction->setData(QVariant(QString::fromStdString(directory + "/" + presets[i])));
        //generate icon
        TransFunc1DKeys preview;
        preview.load(directory + "/" + presets[i]);
        QPixmap icon(32,24);
        QPainter painter(&icon);
        for(int i = 0; i < 32; i++) {
            tgt::Color color = preview.getTexture()->texelAsFloat(static_cast<size_t>((static_cast<float>(i)/31.f)*(preview.getDimensions().x-1)));
            painter.setPen(QColor(color.r*255,color.g*255,color.b*255));
            painter.drawLine(i,0,i,23);
        }
        propAction->setIcon(icon);
        //add tf to menu
        menu->addAction(propAction);
    }
}

void TransFuncPropertyWidget::loadTF(QAction* action) {
    QString data = action->data().toString();
    if(property_ && dynamic_cast<TransFunc1DKeys*>(property_->get())) {
        TransFunc1DKeys* tf1DKeys = dynamic_cast<TransFunc1DKeys*>(property_->get());
        tgt::vec2  oldDomain = tf1DKeys->getDomain();
        TransFuncBase::AlphaMode alphaMode = property_->get()->getAlphaMode();
        if (property_->get()->load(data.toStdString())) {
            if(tf1DKeys->getDomain() == tgt::vec2(0.0f, 1.0f)) {
                tf1DKeys->setDomain(oldDomain);
            }

            property_->get()->setAlphaMode(alphaMode);
            property_->invalidate();
        }
        else {
            //TODO: error
        }
    }
}

//-----------------------------------------------------------------------------------
//    window menu
//-----------------------------------------------------------------------------------
void TransFuncPropertyWidget::populateLoadWindowMenu() {
    loadWindowMenu_->clear();

    QAction* propAction = new QAction(QString::fromStdString("Fit to Data"), loadWindowMenu_);
    propAction->setData(QVariant(QString::fromStdString("fitToData")));
    loadWindowMenu_->addAction(propAction);

    loadWindowMenu_->addSeparator();

    std::string applicationPresetPath = VoreenApplication::app()->getCoreResourcePath("transferfuncs/windows");
    populateWindowMenuFromDirectory(loadWindowMenu_, applicationPresetPath);

    loadWindowMenu_->addSeparator();

    std::string userPresetPath = VoreenApplication::app()->getUserDataPath("transferfuncs/windows");
    populateWindowMenuFromDirectory(loadWindowMenu_, userPresetPath);
}

void TransFuncPropertyWidget::populateWindowMenuFromDirectory(QMenu* menu, std::string directory) {
    std::vector<std::string> subdirs = FileSys.listSubDirectories(directory, true);
    for(size_t i=0; i<subdirs.size(); i++) {
        QMenu* subMenu = menu->addMenu(QString::fromStdString(subdirs[i]));
        populateWindowMenuFromDirectory(subMenu, directory + "/" + subdirs[i]);
    }

    std::vector<std::string> presets = FileSys.listFiles(directory, true);
    for(size_t i=0; i<presets.size(); i++) {
        QAction* propAction = new QAction(QString::fromStdString(FileSys.baseName(presets[i])), menu);
        propAction->setData(QVariant(QString::fromStdString(directory + "/" + presets[i])));
        menu->addAction(propAction);
    }
}

void TransFuncPropertyWidget::loadWindow(QAction* action) {
    QString data = action->data().toString();
    if(property_ && property_->get()) {
        if(data == "fitToData") {
            fitToData();
        }
        else {
            tgt::File* file = FileSys.open(data.toStdString());
            if(file) {
                std::string content = file->getAsString();
                std::vector<std::string> expl = strSplit(content, '\n');
                if(expl.size() == 2) {
                    std::cout << expl[0] << expl[1];
                    float width = stof(expl[0]);
                    float level = stof(expl[1]);

                    tgt::vec2 domain = tgt::vec2(level - (width * 0.5f), level + (width * 0.5f));
                    if(TransFunc1DKeys* tf1DKeys = dynamic_cast<TransFunc1DKeys*>(property_->get()))
                        tf1DKeys->setDomain(domain);
                    property_->invalidate();
                }
                else {
                    LERRORC("voreen.TransFuncPropertyWidget", "Failed to parse window preset file");
                }
            }
        }
    }
    else {
        LERRORC("voreen.TransFuncPropertyWidget", "Failed to open window preset file");
    }
}

//-----------------------------------------------------------------------------------
//    alpha menu
//-----------------------------------------------------------------------------------
void TransFuncPropertyWidget::populateAlphaMenu() {
    alphaMenu_->clear();

    QActionGroup* group = new QActionGroup(alphaMenu_);

    QAction* useAction = new QAction(QString::fromStdString("Use Alpha"), alphaMenu_);
    useAction->setData(QVariant(TransFuncBase::TF_USE_ALPHA));
    useAction->setCheckable(true);
    alphaMenu_->addAction(useAction);
    group->addAction(useAction);

    alphaMenu_->addSeparator();

    QAction* zeroAction = new QAction(QString::fromStdString("Transparent"), alphaMenu_);
    zeroAction->setData(QVariant(TransFuncBase::TF_ZERO_ALPHA));
    zeroAction->setCheckable(true);
    alphaMenu_->addAction(zeroAction);
    group->addAction(zeroAction);

    QAction* oneAction = new QAction(QString::fromStdString("Opaque"), alphaMenu_);
    oneAction->setData(QVariant(TransFuncBase::TF_ONE_ALPHA));
    oneAction->setCheckable(true);
    alphaMenu_->addAction(oneAction);
    group->addAction(oneAction);

    if(property_ && property_->get()) {
        switch(property_->get()->getAlphaMode()) {
        case TransFuncBase::TF_ZERO_ALPHA:
            zeroAction->setChecked(true);
            break;
        case TransFuncBase::TF_USE_ALPHA:
            useAction->setChecked(true);
            break;
        case TransFuncBase::TF_ONE_ALPHA:
            oneAction->setChecked(true);
            break;
        }
    }
}

void TransFuncPropertyWidget::setAlpha(QAction* action) {
    TransFuncBase::AlphaMode mode = static_cast<TransFuncBase::AlphaMode>(action->data().toInt());
    updateAlpha(mode);
    property_->invalidate();
    textureCanvas_->update();
}

void TransFuncPropertyWidget::updateAlpha(TransFuncBase::AlphaMode mode) {
    if(property_ && property_->get()) {
        property_->get()->setAlphaMode(mode);
    }
    switch(mode) {
    case TransFuncBase::TF_ZERO_ALPHA:
        alphaBt_->setIcon(QPixmap(":/qt/icons/alpha_trans.png"));
        break;
    case TransFuncBase::TF_USE_ALPHA:
        alphaBt_->setIcon(QPixmap(":/qt/icons/alpha_use.png"));
        break;
    case TransFuncBase::TF_ONE_ALPHA:
        alphaBt_->setIcon(QPixmap(":/qt/icons/alpha_opaque.png"));
        break;
    }
}

//-----------------------------------------------------------------------------------
//    advanced menu
//-----------------------------------------------------------------------------------
void TransFuncPropertyWidget::populateAdvancedMenu() {
    advancedMenu_->clear();

    QAction* editorAction = new QAction(QString::fromStdString("Show/Hide Editor"), advancedMenu_);
    editorAction->setData(QVariant(0));
    editorAction->setIcon(QPixmap(":/qt/icons/mapping-function.png"));
    advancedMenu_->addAction(editorAction);

    advancedMenu_->addSeparator();

    QAction* invertAction = new QAction(QString::fromStdString("Invert Color Map"), advancedMenu_);
    invertAction->setData(QVariant(1));
    invertAction->setIcon(QPixmap(":/qt/icons/arrow-leftright.png"));
    advancedMenu_->addAction(invertAction);

    advancedMenu_->addSeparator();

    QAction* loadAction = new QAction(QString::fromStdString("Load Transfer Function"), advancedMenu_);
    loadAction->setData(QVariant(2));
    loadAction->setIcon(QPixmap(":/qt/icons/open.png"));
    advancedMenu_->addAction(loadAction);

    QAction* saveAction = new QAction(QString::fromStdString("Save Transfer Function"), advancedMenu_);
    saveAction->setData(QVariant(3));
    saveAction->setIcon(QPixmap(":/qt/icons/save.png"));
    advancedMenu_->addAction(saveAction);

}

void TransFuncPropertyWidget::doAdvancedAction(QAction* action) {
    int data = action->data().toInt();
    switch (data) {
    case 0: //open editor
        setProperty();
        widgetChanged();
        break;
    case 1: //invert color map
        if(property_ && property_->get()) {
            if(TransFunc1DKeys* func = dynamic_cast<TransFunc1DKeys*>(property_->get())) {
                func->invertKeys();
                property_->invalidate();
                textureCanvas_->update();
            } else {
                LINFOC("TransFuncPropertyWidget","Invert color map is not supported for 2d color maps");
            }
        }
        break;
    case 2: // load tf
        if(property_ && property_->get() && TransFuncIOHelperQt::loadTransferFunction(property_->get())) {
            property_->invalidate();
            textureCanvas_->update();
        }
        break;
    case 3: // save tf
        if(property_ && property_->get()) {
            TransFuncIOHelperQt::saveTransferFunction(property_->get());
        }
        break;
    default:
        tgtAssert(false,"Unknown Action!");
        break;
    }
}

void TransFuncPropertyWidget::setVisible(bool visible) {
    QPropertyWidgetWithEditorWindow::setVisible(visible);
    textureCanvas_->update();
}

void TransFuncPropertyWidget::updateFromPropertySlot() {
    if (plugin_)
        plugin_->updateFromProperty();

    const VolumeBase* vb = property_->getVolumeHandle();
    if(property_->getMetaDataContainer().hasMetaData("TransfuncPropertyWidgetPainterZoom")) {
        tgt::vec2 range = static_cast<Vec2MetaData*>(property_->getMetaDataContainer().getMetaData("TransfuncPropertyWidgetPainterZoom"))->getValue();
        texturePainter_->setTransFunc(dynamic_cast<TransFunc1DKeys*>(property_->get()), range.x, range.y);
    } else {
        texturePainter_->setTransFunc(dynamic_cast<TransFunc1DKeys*>(property_->get()));
    }
    // We want the view to reset when the TF is deserialized.
    // To achieve this we set the widget as not initialized as long as only the standard function has been set
    const TransFuncBase* tf = property_->get();
    /*if(!viewInitialized_) {
        if(tf) {
            const TransFunc1DKeys* tf1d = dynamic_cast<const TransFunc1DKeys*>(tf);
            if(tf1d) {
                if(!tf1d->isStandardFunc())
                    viewInitialized_ = true;
            }
            else
                viewInitialized_ = true;
        }

        if(viewInitialized_)
            resetZoom();
    }*/

    clearObserveds();
    if(vb) {
        vb->addObserver(this);

        if(vb->hasDerivedData<VolumeHistogramIntensity>())
            texturePainter_->setHistogram(&(vb->getDerivedData<VolumeHistogramIntensity>()->getHistogram(property_->getVolumeChannel())));
        else
            vb->getDerivedDataThreaded<VolumeHistogramIntensity>();
    }
    else {
        texturePainter_->setHistogram(0);
    }

    if(tf) {
        alphaBt_->blockSignals(true);
        updateAlpha(tf->getAlphaMode());
        alphaBt_->blockSignals(false);
    }

    textureCanvas_->update();
}

void TransFuncPropertyWidget::showToolTipSlot(QPoint pos, QString tip) {
    QToolTip::showText(mapToGlobal(pos),tip);
}

void TransFuncPropertyWidget::hideToolTipSlot() {
    QToolTip::hideText();
}


void TransFuncPropertyWidget::zoomIn() {
    texturePainter_->zoomIn();
    textureCanvas_->update();
    property_->getMetaDataContainer().addMetaData("TransfuncPropertyWidgetPainterZoom",
                                      new Vec2MetaData(tgt::vec2(texturePainter_->getMinDomainValue(),texturePainter_->getMaxDomainValue())));
    emit zoomLevelChanged();
}

void TransFuncPropertyWidget::zoomOut() {
    texturePainter_->zoomOut();
    textureCanvas_->update();
    property_->getMetaDataContainer().addMetaData("TransfuncPropertyWidgetPainterZoom",
                                      new Vec2MetaData(tgt::vec2(texturePainter_->getMinDomainValue(),texturePainter_->getMaxDomainValue())));
    emit zoomLevelChanged();
}

void TransFuncPropertyWidget::resetZoom() {
    texturePainter_->resetZoom();
    textureCanvas_->update();
    property_->getMetaDataContainer().addMetaData("TransfuncPropertyWidgetPainterZoom",
                                      new Vec2MetaData(tgt::vec2(texturePainter_->getMinDomainValue(),texturePainter_->getMaxDomainValue())));
    emit zoomLevelChanged();
}

void TransFuncPropertyWidget::fitToData() {
    property_->fitDomainToData();
    resetZoom();
    property_->getMetaDataContainer().addMetaData("TransfuncPropertyWidgetPainterZoom",
                                      new Vec2MetaData(tgt::vec2(texturePainter_->getMinDomainValue(),texturePainter_->getMaxDomainValue())));
}

void TransFuncPropertyWidget::setProperty() {
    if (!disconnected_) {
        // lazy instantiation of transfunc editor window
        if (!editorWindow_) {
            createEditorWindow(Qt::NoDockWidgetArea);
            tgtAssert(editorWindow_, "Transfunc editor not instantiated");
        }

        if (editorWindow_->isVisible()) {
            //close widget
            editorWindow_->close();
        }
        else {
            //open Widget
            editorWindow_->showNormal();
            plugin_->updateFromProperty();
        }
    }
}

void TransFuncPropertyWidget::disconnect() {
    disconnected_ = true;
    if (plugin_)
        plugin_->disconnect();
}

QWidget* TransFuncPropertyWidget::createEditorWindowWidget() {
    plugin_ = new TransFuncPluginXXX(property_, parentWidget(), Qt::Horizontal);
    plugin_->createWidgets();
    plugin_->createConnections();
    //connect(plugin_, SIGNAL(transferFunctionChanged()), this, SIGNAL(modified()));

    return plugin_;
}

void TransFuncPropertyWidget::customizeEditorWindow() {
    editorWindow_->setAllowedAreas(Qt::NoDockWidgetArea);
    editorWindow_->setFloating(true);
}

Property* TransFuncPropertyWidget::getProperty() {
    return property_;
}

void TransFuncPropertyWidget::showNameLabel(bool visible) {
    if (nameLabel_) {
        if (visible)
            nameLabel_->show();
        else
            nameLabel_->hide();
    }
}

CustomLabel* TransFuncPropertyWidget::getOrCreateNameLabel() const {
    return 0;
}

void TransFuncPropertyWidget::toggleInteractionMode(bool on) {
    QualityMode.requestQualityMode((on? VoreenQualityMode::RQ_INTERACTIVE : VoreenQualityMode::RQ_DEFAULT), this);
}

void TransFuncPropertyWidget::invalidateProperty() {
    property_->invalidate();
}

void TransFuncPropertyWidget::updateZoomMeta() {
    property_->getMetaDataContainer().addMetaData("TransfuncPropertyWidgetPainterZoom",
                                      new Vec2MetaData(tgt::vec2(texturePainter_->getMinDomainValue(),texturePainter_->getMaxDomainValue())));
    property_->invalidate();
}

void TransFuncPropertyWidget::volumeDelete(const VolumeBase* source) {
    texturePainter_->setHistogram(0);
    textureCanvas_->update();
}

void TransFuncPropertyWidget::volumeChange(const VolumeBase* source) {
    texturePainter_->setHistogram(0);
    textureCanvas_->update();
}

void TransFuncPropertyWidget::derivedDataThreadFinished(const VolumeBase* source, const VolumeDerivedData* derivedData) {
    if(dynamic_cast<const VolumeHistogramIntensity*>(derivedData)) {
        texturePainter_->setHistogram(&(static_cast<const VolumeHistogramIntensity*>(derivedData)->getHistogram(property_->getVolumeChannel())));
        textureCanvas_->update();
    }
}


} // namespace voreen
