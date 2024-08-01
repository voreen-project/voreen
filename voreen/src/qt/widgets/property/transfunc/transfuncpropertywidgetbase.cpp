/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2024 University of Muenster, Germany,                        *
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

#include "voreen/qt/widgets/property/transfunc/transfuncpropertywidgetbase.h"

#include "voreen/core/voreenapplication.h"
#include "voreen/core/properties/transfunc/transfuncpropertybase.h"

#include "voreen/qt/widgets/property/transfunc/utils/transfunciohelperqt.h"
#include "voreen/qt/widgets/property/transfunc/transfuncpropertyeditorbase.h"
#include "voreen/qt/widgets/property/transfunc/transfuncpropertywidgetpainterbase.h"
#include "voreen/qt/widgets/customlabel.h"
#include "voreen/qt/widgets/voreentoolwindow.h"

#include "tgt/qt/qtcanvas.h"
#include "tgt/filesystem.h"
#include "tgt/glcontextmanager.h"

#include <QMenu>
#include <QToolButton>
#include <QToolTip>


namespace voreen {

TransFuncPropertyWidgetBase::TransFuncPropertyWidgetBase(TransFuncPropertyBase* prop, QWidget* parent)
    : QPropertyWidgetWithToolWindow(prop, parent, true, true)
    , presetButton_(0), presetMenu_(0), domainButton_(0), domainMenu_(0)
    , alphaButton_(0), alphaMenu_(0), zoomInButton_(0), zoomOutButton_(0)
    , advancedButton_(0), advancedMenu_(0)
    , previewCanvas_(0), previewPainter_(0)
    , property_(prop), editor_(0)
{}

TransFuncPropertyWidgetBase::~TransFuncPropertyWidgetBase() {
}

//----------------------------------------------------------
//      Gui Functions
//----------------------------------------------------------
CustomLabel* TransFuncPropertyWidgetBase::getOrCreateNameLabel() const {
    return 0;
}

void TransFuncPropertyWidgetBase::initialize() {
    //initialize gui components
    initializeBaseLayout();
    //initialize gui values
    updateFromProperty();
}

void TransFuncPropertyWidgetBase::initializeBaseLayout() {
    //preset button
    presetMenu_ = new QMenu();
    presetButton_ = new QToolButton();
    presetButton_->setIcon(QPixmap(":/qt/icons/colorize.png"));
    presetButton_->setToolTip("Color Map Presets");
    presetButton_->setMenu(presetMenu_);
    presetButton_->setPopupMode(QToolButton::InstantPopup);
    connect(presetMenu_, SIGNAL(aboutToShow()), this, SLOT(createPresetMenuSlot()));
    connect(presetMenu_, SIGNAL(triggered(QAction*)), this, SLOT(loadPresetSlot(QAction*)));

    //domain button
    domainMenu_ = new QMenu();
    domainButton_ = new QToolButton();
    domainButton_->setIcon(QPixmap(":/qt/icons/histogram_fit.png"));
    domainButton_->setToolTip("Window Fitting/Presets");
    domainButton_->setMenu(domainMenu_);
    domainButton_->setPopupMode(QToolButton::InstantPopup);
    connect(domainMenu_, SIGNAL(aboutToShow()), this, SLOT(createDomainMenuSlot()));
    connect(domainMenu_, SIGNAL(triggered(QAction*)), this, SLOT(loadDomainSlot(QAction*)));

    //alpha button
    alphaMenu_ = new QMenu();
    alphaButton_ = new QToolButton();
    alphaButton_->setIcon(QPixmap(":/qt/icons/alpha.png"));
    alphaButton_->setToolTip("Malipulate Alpha Channel");
    alphaButton_->setMenu(alphaMenu_);
    alphaButton_->setPopupMode(QToolButton::InstantPopup);
    connect(alphaMenu_, SIGNAL(aboutToShow()), this, SLOT(createAlphaMenuSlot()));
    connect(alphaMenu_, SIGNAL(triggered(QAction*)), this, SLOT(loadAlphaSlot(QAction*)));
    if(property_ && property_->get())
        updateAlpha(property_->get()->getAlphaMode());

    //zoom buttons
    zoomInButton_ = new QToolButton();
    zoomInButton_->setIcon(QPixmap(":/qt/icons/viewmag+.png"));
    zoomInButton_->setToolTip("Zoom In");
    zoomOutButton_ = new QToolButton();
    zoomOutButton_->setIcon(QPixmap(":/qt/icons/viewmag_.png"));
    zoomOutButton_->setToolTip("Zoom Out");
    connect(zoomInButton_, SIGNAL(clicked()), this, SLOT(zoomInSlot()));
    connect(zoomOutButton_, SIGNAL(clicked()), this, SLOT(zoomOutSlot()));

    //advanced button
    advancedMenu_ = new QMenu();
    advancedButton_ = new QToolButton();
    advancedButton_->setIcon(QPixmap(":/qt/icons/mapping-function.png"));
    advancedButton_->setToolTip("Advanced Settings");
    advancedButton_->setMenu(advancedMenu_);
    advancedButton_->setPopupMode(QToolButton::InstantPopup);
    connect(advancedMenu_, SIGNAL(aboutToShow()), this, SLOT(createAdvancedMenuSlot()));
    connect(advancedMenu_, SIGNAL(triggered(QAction*)), this, SLOT(doAdvancedActionSlot(QAction*)));
    if (isToolWindowVisibleOnStartup())
        createToolWindow(Qt::NoDockWidgetArea);

    //preview canvas
    previewCanvas_ = new tgt::QtCanvas("TransFunc", tgt::ivec2(10, 10), tgt::GLCanvas::RGBADD, 0);
    previewPainter_ = createPreviewCanvasPainter(previewCanvas_);
    previewPainter_->setTransFunc(property_->get());
    previewCanvas_->getEventHandler()->addListenerToBack(previewPainter_);
    connect(previewPainter_,SIGNAL(changedGamma()),this,SLOT(invalidateProperty()));
    connect(previewPainter_,SIGNAL(changedDomain()),this,SLOT(invalidateProperty()));
    connect(previewPainter_,SIGNAL(changedDomain()),this,SLOT(storeZoomMetaData()));
    connect(previewPainter_,SIGNAL(storeZoomMetaDataSignal()),this,SLOT(storeZoomMetaData()));
    connect(previewPainter_,SIGNAL(interaction(bool)),this,SLOT(toggleInteractionMode(bool)));
    connect(previewPainter_,SIGNAL(showInfoToolTip(QPoint, QString)),this,SLOT(showToolTipSlot(QPoint, QString)));
    connect(previewPainter_,SIGNAL(hideInfoToolTip()),this,SLOT(hideToolTipSlot()));

    //layout components
    //main layout
    QVBoxLayout* mainLayout = new QVBoxLayout();
    mainLayout->setAlignment(Qt::AlignLeft);
    mainLayout->setContentsMargins(0, 0, 0, 0);
    mainLayout->setSpacing(2);

    //button layout
    QGridLayout* buttonLayout = new QGridLayout();
    buttonLayout->setSpacing(2);

    buttonLayout->setColumnStretch(0,4); //name label is much bigger
    buttonLayout->setColumnStretch(1,1);
    buttonLayout->setColumnStretch(2,1);
    buttonLayout->setColumnStretch(3,1);
    buttonLayout->setColumnStretch(4,1);
    buttonLayout->setColumnStretch(5,1);
    buttonLayout->setColumnStretch(6,1);
    buttonLayout->setColumnStretch(7,1);
    buttonLayout->setColumnStretch(8,1);

    //set button layout
    buttonLayout->addWidget(QPropertyWidget::getOrCreateNameLabel(),0,0,1,1);
    nameLabel_->setParent(this);

    buttonLayout->addWidget(presetButton_,0,1,1,1);
    buttonLayout->addWidget(domainButton_,0,2,1,1);
    buttonLayout->addWidget(alphaButton_,0,3,1,1);

    buttonLayout->addWidget(zoomInButton_,0,5,1,1);
    buttonLayout->addWidget(zoomOutButton_,0,6,1,1);

    buttonLayout->addWidget(advancedButton_,0,8,1,1);

    mainLayout->addLayout(buttonLayout);

    switch(getBaseLayout()) {
    case WBL_1D:
        initialize1DLayout(mainLayout);
        break;
    case WBL_2D:
        initialize2DLayout(mainLayout);
        break;
    default:
        tgtAssert(false,"undefined layout");
        break;
    }
}

void TransFuncPropertyWidgetBase::initialize1DLayout(QLayout* mainLayout) {
    //canvas
    previewCanvas_->setFixedHeight(30);
    previewCanvas_->setSizePolicy(QSizePolicy::MinimumExpanding, QSizePolicy::Fixed);

    mainLayout->addWidget(previewCanvas_);
    layout_->addLayout(mainLayout);
    previewCanvas_->init();
}

void TransFuncPropertyWidgetBase::initialize2DLayout(QLayout* mainLayout) {
    //canvas
    previewCanvas_->setFixedHeight(200);
    previewCanvas_->setFixedWidth(200);
    previewCanvas_->setSizePolicy(QSizePolicy::Fixed, QSizePolicy::Fixed);
    QHBoxLayout* tmpLayout = new QHBoxLayout();
    tmpLayout->setAlignment(Qt::AlignCenter);
    tmpLayout->addWidget(previewCanvas_);
    mainLayout->addItem(tmpLayout);
    layout_->addLayout(mainLayout);
    previewCanvas_->init();
}

//----------------------------------------------------------
//      Menu Handling
//----------------------------------------------------------
    //--------------------
    //    preset menu
    //--------------------
void TransFuncPropertyWidgetBase::createPresetMenuSlot() {
    presetMenu_->clear();

    //hard coded preset path
    std::string applicationPresetPath = VoreenApplication::app()->getCoreResourcePath("transferfuncs/presets/" + getPathSubFolder());

    //the first number of the preset name encodes the subsection
    int currentSubSection = 0;
    std::vector<std::string> presets = FileSys.listFiles(applicationPresetPath, true);
    for(size_t i=0; i<presets.size(); i++) {
        //add seperator after sub section
        if(FileSys.baseName(presets[i])[0] != currentSubSection) {
            presetMenu_->addSeparator();
            currentSubSection = FileSys.baseName(presets[i])[0];
        }
        //create action and cut off prefix
        QAction* propAction;
        propAction = new QAction(QString::fromStdString(FileSys.baseName(presets[i]).substr(4)), presetMenu_);
        propAction->setData(QVariant(QString::fromStdString(applicationPresetPath + "/" + presets[i])));
        //generate icon
        std::string iconPath = applicationPresetPath + "/" + presets[i];
        propAction->setIcon(createPresetPreviewIcon(iconPath));
        //add tf to menu
        presetMenu_->addAction(propAction);
    }
}

    //--------------------
    //    domain menu
    //--------------------
void TransFuncPropertyWidgetBase::createDomainMenuSlot() {
    domainMenu_->clear();

    QAction* propAction = new QAction(QString::fromStdString("Fit to Data"), domainMenu_);
    propAction->setData(QVariant(QString::fromStdString("fitToData")));
    domainMenu_->addAction(propAction);
}

void TransFuncPropertyWidgetBase::loadDomainSlot(QAction* action) {
    if(property_ && property_->get()) {
        QString data = action->data().toString();
        if(data == "fitToData") {
            fitToDataSlot();
        }
    }
}

void TransFuncPropertyWidgetBase::fitToDataSlot() {
    property_->applyDomainFromData();
    resetZoomSlot();
    storeZoomMetaData();
}

    //--------------------
    //    alpha menu
    //--------------------
void TransFuncPropertyWidgetBase::createAlphaMenuSlot() {
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

void TransFuncPropertyWidgetBase::loadAlphaSlot(QAction* action) {
    TransFuncBase::AlphaMode mode = static_cast<TransFuncBase::AlphaMode>(action->data().toInt());
    updateAlpha(mode);
    property_->invalidate();
}

    //--------------------
    //    zoom menu
    //--------------------
void TransFuncPropertyWidgetBase::zoomInSlot() {
    previewPainter_->zoomIn();
    storeZoomMetaData();
}

void TransFuncPropertyWidgetBase::zoomOutSlot() {
    previewPainter_->zoomOut();
    storeZoomMetaData();
}

void TransFuncPropertyWidgetBase::resetZoomSlot() {
    previewPainter_->resetZoom();
    storeZoomMetaData();
}

    //--------------------
    //    advanced menu
    //--------------------
enum AdvancedActions {
    AA_OPEN_EDITOR = 0,
    AA_CHANGE_PLUGIN = 5,
    AA_LOAD_TF = 20,
    AA_SAVE_TF = 30
};

void TransFuncPropertyWidgetBase::createAdvancedMenuSlot() {
    advancedMenu_->clear();

    QAction* editorAction = new QAction(QString::fromStdString("Show/Hide Editor"), advancedMenu_);
    editorAction->setData(QVariant(AA_OPEN_EDITOR));
    editorAction->setIcon(QPixmap(":/qt/icons/mapping-function.png"));
    advancedMenu_->addAction(editorAction);

    advancedMenu_->addSeparator();

    QAction* loadAction = new QAction(QString::fromStdString("Load Transfer Function"), advancedMenu_);
    loadAction->setData(QVariant(AA_LOAD_TF));
    loadAction->setIcon(QPixmap(":/qt/icons/open.png"));
    advancedMenu_->addAction(loadAction);

    QAction* saveAction = new QAction(QString::fromStdString("Save Transfer Function"), advancedMenu_);
    saveAction->setData(QVariant(AA_SAVE_TF));
    saveAction->setIcon(QPixmap(":/qt/icons/save.png"));
    advancedMenu_->addAction(saveAction);
}

void TransFuncPropertyWidgetBase::doAdvancedActionSlot(QAction* action) {
    int data = action->data().toInt();
    switch (data) {
    case AA_OPEN_EDITOR: //open editor
        //lazy instanziation
        if(!toolWindow_) {
            createToolWindow(Qt::NoDockWidgetArea); //< @see QPropertyWidgetWithToolWindow
            toolWindow_->setVisible(true);
        }else
            toggleToolWindow(); //< @see QPropertyWidgetWithToolWindow
        //widgetChanged(); // y?
        break;
    case AA_LOAD_TF: // load tf
        if(property_ && property_->get() && TransFuncIOHelperQt::loadTransferFunction(property_->get())) {
            property_->invalidate();
            previewCanvas_->update();
        }
        break;
    case AA_SAVE_TF: // save tf
        if(property_ && property_->get()) {
            TransFuncIOHelperQt::saveTransferFunction(property_->get());
        }
        break;
    default:
        tgtAssert(false,"Unknown Action!");
        break;
    }
}

//----------------------------------------------------------
//      Update Functions
//----------------------------------------------------------
void TransFuncPropertyWidgetBase::updateAlpha(TransFuncBase::AlphaMode mode) {
    if(property_ && property_->get()) {
        property_->get()->setAlphaMode(mode);
    }
    switch(mode) {
    case TransFuncBase::TF_ZERO_ALPHA:
        alphaButton_->setIcon(QPixmap(":/qt/icons/alpha_trans.png"));
        break;
    case TransFuncBase::TF_USE_ALPHA:
        alphaButton_->setIcon(QPixmap(":/qt/icons/alpha_use.png"));
        break;
    case TransFuncBase::TF_ONE_ALPHA:
        alphaButton_->setIcon(QPixmap(":/qt/icons/alpha_opaque.png"));
        break;
    }

}

void TransFuncPropertyWidgetBase::customizeToolWindow() {
    if(toolWindow_) {
        toolWindow_->setAllowedAreas(Qt::NoDockWidgetArea);
        toolWindow_->setFloating(true);
    }
}

void TransFuncPropertyWidgetBase::showToolTipSlot(QPoint pos, QString tip) {
    QToolTip::showText(mapToGlobal(pos),tip);
}

void TransFuncPropertyWidgetBase::hideToolTipSlot() {
    QToolTip::hideText();
}

void TransFuncPropertyWidgetBase::updateFromPropertySlot() {
    if(!property_) return;

    //Update editor Window
    if (editor_)
        editor_->updateFromProperty();

    //Updates histogram from volume changes
    const VolumeBase* vb = property_->getVolume();
    if(vb) {
        Histogram* histo = getHistogramFromVolume(vb, property_->getVolumeChannel());
        previewPainter_->setHistogram(histo);
    } else {
        previewPainter_->setHistogram(0);
    }

    //Updates from tf
    TransFuncBase* tf = property_->get();
    if(tf) {
        //use meta data correctly
        restoreZoomMetaData();
        //set tf
        previewPainter_->setTransFunc(tf);
        //udate alpha button
        alphaButton_->blockSignals(true);
        updateAlpha(tf->getAlphaMode());
        alphaButton_->blockSignals(false);
    }

    //force repaint
    previewCanvas_->update();
}

} // namespace voreen
