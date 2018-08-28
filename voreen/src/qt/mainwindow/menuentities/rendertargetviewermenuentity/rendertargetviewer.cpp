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

#include "voreen/qt/mainwindow/menuentities/rendertargetviewermenuentity/rendertargetviewer.h"
#include "voreen/qt/mainwindow/menuentities/rendertargetviewermenuentity/rendertargetviewerpainter.h"

#include "voreen/core/voreenapplication.h"
#include "voreen/core/network/processornetwork.h"
#include "voreen/core/network/networkevaluator.h"
#include "voreen/core/ports/renderport.h"
#include "voreen/core/datastructures/rendertarget/rendertarget.h"
#include "voreen/qt/utils/voreenqtworkspacehandler.h"

#include "voreen/qt/widgets/voreentoolwindow.h"

#include "tgt/qt/qtcanvas.h"
#include "tgt/font.h"
#include "tgt/vector.h"
#include "tgt/gpucapabilities.h"
#include "tgt/logmanager.h"

#include <math.h>
#include <QMouseEvent>
#include <QString>
#include <QMenu>
#include <QFileDialog>
#include <QApplication>
#include <QUrl>
#include <QMessageBox>
#include <QMainWindow>
#include <QLayout>
#include <QDesktopServices>
#include <QWidget>

#include <iostream>
#include <sstream>
#include <iomanip>

using namespace std;
using tgt::Texture;


namespace voreen {

const std::string RenderTargetViewer::loggerCat_ = "voreen.RenderTargetViewer";

RenderTargetViewer::RenderTargetViewer()
    : QWidget()
    , evaluator_(0)
    , canvas_(0)
    , selectedRenderPortIndex_(-1)
    , maximizeOnePort_(false)
    , mouseX_(0)
    , mouseY_(0)
    , mouseIsInside_(false)
    , zoomScale_(1)
    , zoomTranslateX_(0)
    , zoomTranslateY_(0)
    , zoomMouseX_(0)
    , zoomMouseY_(0)
    , zoomOffsetX_(0)
    , zoomOffsetY_(0)

{
    setMouseTracking(true);
    setFocusPolicy(Qt::StrongFocus);

    contextMenuMEN_ = new QMenu(this);
    contextMenuMEN_->setParent(this);
    backgroundBlackACT_ = contextMenuMEN_->addAction("Background: Black");
    backgroundBlackACT_->setCheckable(true);
    backgroundBlackACT_->setChecked(true);
    backgroundWhiteACT_ = contextMenuMEN_->addAction("Background: White");
    backgroundWhiteACT_->setCheckable(true);
    backgroundCheckerboardPatternACT_ = contextMenuMEN_->addAction("Background: Checkerboard Pattern");
    backgroundCheckerboardPatternACT_->setCheckable(true);

    contextMenuMEN_->addSeparator();

    colorRACT_ = contextMenuMEN_->addAction("Activate R Channel");
    colorRACT_->setCheckable(true);
    colorRACT_->setChecked(true);
    colorGACT_ = contextMenuMEN_->addAction("Activate G Channel");
    colorGACT_->setCheckable(true);
    colorGACT_->setChecked(true);
    colorBACT_ = contextMenuMEN_->addAction("Activate B Channel");
    colorBACT_->setCheckable(true);
    colorBACT_->setChecked(true);
    colorAACT_ = contextMenuMEN_->addAction("Activate A Channel");
    colorAACT_->setCheckable(true);
    colorAACT_->setChecked(true);

    contextMenuMEN_->addSeparator();

    colorRGBAACT_ = contextMenuMEN_->addAction("Color Channel");
    colorRGBAACT_->setCheckable(true);
    alphaOnlyACT_ = contextMenuMEN_->addAction("Alpha Channel");
    alphaOnlyACT_->setCheckable(true);
    depthOnlyACT_ = contextMenuMEN_->addAction("Depth Buffer");
    depthOnlyACT_->setCheckable(true);
    hOnlyACT_ = contextMenuMEN_->addAction("Hue Channel");
    hOnlyACT_->setCheckable(true);
    sOnlyACT_ = contextMenuMEN_->addAction("Saturation Channel");
    sOnlyACT_->setCheckable(true);
    vOnlyACT_ = contextMenuMEN_->addAction("Value Channel");
    vOnlyACT_->setCheckable(true);

    contextMenuMEN_->addSeparator();

    filterPortyBySelectedProcessorsACT_ = contextMenuMEN_->addAction("Display Selected Processor Ports Only");
    filterPortyBySelectedProcessorsACT_->setCheckable(true);
    filterPortyBySelectedProcessorsACT_->setChecked(false);

    keepAspectRatioACT_ = contextMenuMEN_->addAction("Keep Aspect Ratio");
    keepAspectRatioACT_->setCheckable(true);
    keepAspectRatioACT_->setChecked(true);

    contextMenuMEN_->addSeparator();

    showInfosACT_ = contextMenuMEN_->addAction("Show Infos");
    showInfosACT_->setCheckable(true);
    showInfosACT_->setChecked(true);
    showInfosDetailsACT_ = contextMenuMEN_->addAction("Show Info Details");
    showInfosDetailsACT_->setCheckable(true);
    showInfosDetailsACT_->setChecked(true);

    contextMenuMEN_->addSeparator();

    saveScreenshotACT_ = contextMenuMEN_->addAction("Save As Image...");
    saveScreenshotWithOverlayACT_ = contextMenuMEN_->addAction("Save As Image With Overlay...");

    contextMenuMEN_->addSeparator();

    backgroundToShowACG_ = new QActionGroup(this);
    backgroundToShowACG_->addAction(backgroundBlackACT_);
    backgroundToShowACG_->addAction(backgroundWhiteACT_);
    backgroundToShowACG_->addAction(backgroundCheckerboardPatternACT_);

    typeToShowACG_ = new QActionGroup(this);
    typeToShowACG_->addAction(colorRGBAACT_);
    typeToShowACG_->addAction(alphaOnlyACT_);
    typeToShowACG_->addAction(depthOnlyACT_);
    typeToShowACG_->addAction(hOnlyACT_);
    typeToShowACG_->addAction(sOnlyACT_);
    typeToShowACG_->addAction(vOnlyACT_);
    setFocusPolicy(Qt::StrongFocus);

    //create rendering canvas
    canvas_ = new tgt::QtCanvas("RenderTargetViewer", tgt::ivec2(10, 10), tgt::GLCanvas::RGBADD, 0);
    QHBoxLayout* mainLayout = new QHBoxLayout();
    mainLayout->addWidget(canvas_);
    setLayout(mainLayout);

    RenderTargetViewerPainter* painter = new RenderTargetViewerPainter(canvas_, this);
    canvas_->setMouseTracking(true);
    canvas_->init();

    mouseIsInside_ = canvas_->underMouse();
}

RenderTargetViewer::~RenderTargetViewer() {
    deinitialize();
}

void RenderTargetViewer::afterNetworkProcess() {
    if (isVisible())
        canvas_->update();
}

void RenderTargetViewer::setWorkspace(Workspace* workspace) {
    if (isVisible())
        canvas_->update();
}

void RenderTargetViewer::setEvaluator(NetworkEvaluator* evaluator) {
    if (evaluator_)
        evaluator_->removeObserver(this);

    evaluator_ = evaluator;

    if (evaluator_)
        evaluator_->addObserver(this);

    if (isVisible())
        canvas_->update();
}

void RenderTargetViewer::initialize() {
}

void RenderTargetViewer::deinitialize() {
    if (evaluator_)
        evaluator_->removeObserver(this);
}

//-----------------------------------------------
//  Helper Functions
//-----------------------------------------------
struct CompProcID {
    inline bool operator() (const Port* p1, const Port* p2) {
        return (p1->getProcessor()->getID() < p2->getProcessor()->getID());
    }
};

std::vector<RenderPort*> RenderTargetViewer::collectRenderPorts() {
    tgtAssert(evaluator_, "no evaluator");

    // In case the evaluator is locked, we break here and request a new update.
    if (evaluator_->isLocked()) {
        canvas_->update();
        return std::vector<RenderPort*>();
    }

    std::vector<RenderPort*> collectedRenderPorts = evaluator_->collectRenderPorts();

    //return all ports if no ports are selected
    if(collectedRenderPorts.empty() || !filterPortyBySelectedProcessorsACT_->isChecked()) {
        std::sort(collectedRenderPorts.begin(), collectedRenderPorts.end(), CompProcID());
        return collectedRenderPorts;
    }

    // return only selected ports
    std::vector<RenderPort*> selectedPorts;
    for(unsigned int i=0; i<collectedRenderPorts.size(); i++) {
        if(std::find(selectedProcessors_.begin(), selectedProcessors_.end(), collectedRenderPorts[i]->getProcessor()) != selectedProcessors_.end())
            selectedPorts.push_back(collectedRenderPorts[i]);
    }
    std::sort(selectedPorts.begin(), selectedPorts.end(), CompProcID());

    return selectedPorts;
}

void RenderTargetViewer::processorsSelected(const QList<Processor*>& processors) {
    selectedProcessors_ = processors;
    maximizeOnePort_ = false;
    canvas_->update();
}

void RenderTargetViewer::takeScreenshot(bool withOverlay) {
    QString path;

    //tgtAssert(initialized_, "Not initialized");
    tgtAssert(evaluator_, "No evaluator");

    if (!maximizeOnePort_) {
        // should not get here
        QMessageBox::critical(this, tr("Error saving screenshot"), tr("Only supported if one port is visible."));
        return;
    }

    std::vector<RenderPort*> renderPorts = collectRenderPorts();

    tgtAssert(selectedRenderPortIndex_ >= 0 && selectedRenderPortIndex_ < (int)renderPorts.size(), "Invalid render port index");

    QFileDialog filedialog(this);
    filedialog.setWindowTitle(tr("Save Screenshot"));
    filedialog.setDirectory(VoreenApplication::app()->getUserDataPath("screenshots").c_str());
    filedialog.setDefaultSuffix(tr("png"));

    QStringList filter;
    filter << tr("PNG image (*.png)");
    filter << tr("JPEG image (*.jpg)");
    filedialog.setNameFilters(filter);
    filedialog.setAcceptMode(QFileDialog::AcceptSave);
    filedialog.setOption(QFileDialog::DontUseNativeDialog);

    QList<QUrl> urls;
    urls << QUrl::fromLocalFile(VoreenApplication::app()->getUserDataPath("screenshots").c_str());
    for (auto f : QStandardPaths::standardLocations(QStandardPaths::DesktopLocation))
        urls << QUrl::fromLocalFile(f);
    for (auto f : QStandardPaths::standardLocations(QStandardPaths::HomeLocation))
        urls << QUrl::fromLocalFile(f);
    filedialog.setSidebarUrls(urls);

    struct tm* Tm;
    time_t currentTime = time(NULL);
    Tm = localtime(&currentTime);
    std::stringstream timestamp;
    timestamp << "screenshot " << (Tm->tm_year+1900) << "-" << (Tm->tm_mon+1) << "-" << Tm->tm_mday << "-" << Tm->tm_hour << "-" << Tm->tm_min << "-" << Tm->tm_sec;
    timestamp << ".png";
    filedialog.selectFile(tr(timestamp.str().c_str()));

    QStringList fileList;
    if (filedialog.exec())
        fileList = filedialog.selectedFiles();
    if (fileList.empty())
        return;

    path = filedialog.directory().absolutePath();

    if (!fileList.at(0).endsWith(".jpg") && !fileList.at(0).endsWith(".png")) {
        std::string text = "Screenshot could not be saved.\n";
        int index = fileList[0].lastIndexOf(".");
        if ((index == -1) || (index+1 == fileList[0].size()))
            text += "No file extension specified.";
        else
            text += "Invalid file extension: " + fileList[0].right(fileList[0].size()-index-1).toStdString();

        QMessageBox::critical(this, tr("Error saving screenshot"), tr(text.c_str()));
        return;
    }

    QApplication::setOverrideCursor(QCursor(Qt::WaitCursor));
    qApp->processEvents();
    try {
        try {
            if(withOverlay) {
                if(!canvas_->grabFramebuffer().save(tr(fileList.at(0).toStdString().c_str())))
                    throw VoreenException("Could not grab framebuffer.");
            } else {
                // repaint the widget without overlay if necessary, then grab screenshot
                bool showInfos = false;
                if(showInfosACT_->isChecked()) {
                    showInfos = true;
                    showInfosACT_->setChecked(false);
                    canvas_->update();
                }

                if(!canvas_->grabFramebuffer().save(tr(fileList.at(0).toStdString().c_str())))
                    throw VoreenException("Could not grab framebuffer.");

                if(showInfos)
                    showInfosACT_->setChecked(true);

            }
        }
        catch (const VoreenException& e) {
            QString text = tr("Screenshot could not be saved:\n%1").arg(e.what());
            QMessageBox::warning(this, tr("Error saving screenshot"), text);
        }
    }
    catch (const std::exception& e) {
        QString text = tr("Screenshot could not be saved:\n%1").arg(e.what());
        QMessageBox::warning(this, tr("Error saving screenshot"), text);
    }
    QApplication::restoreOverrideCursor();
}

//-----------------------------------------------
//  Event Handling
//-----------------------------------------------

void RenderTargetViewer::mousePressEvent(QMouseEvent* e) {
    updateMousePosition(e);
    if (!evaluator_) {
        selectedRenderPortIndex_ = -1;
        return;
    }

    if (!maximizeOnePort_)
        updateSelected();

    if (selectedRenderPortIndex_ < 0)
        return;

    if (e->button() == Qt::RightButton) {
        std::vector<RenderPort*> renderPorts = collectRenderPorts();
        if (selectedRenderPortIndex_ >= 0 && selectedRenderPortIndex_ < (int)renderPorts.size()) {
            RenderTarget* rt = renderPorts[selectedRenderPortIndex_]->getRenderTarget();
            colorRGBAACT_->setEnabled(false);
            alphaOnlyACT_->setEnabled(false);
            depthOnlyACT_->setEnabled(false);
            hOnlyACT_->setEnabled(false);
            sOnlyACT_->setEnabled(false);
            vOnlyACT_->setEnabled(false);

            colorRGBAACT_->setVisible(true);
            alphaOnlyACT_->setVisible(true);
            depthOnlyACT_->setVisible(true);
            hOnlyACT_->setVisible(true);
            sOnlyACT_->setVisible(true);
            vOnlyACT_->setVisible(true);

            backgroundBlackACT_->setChecked((showType_[selectedRenderPortIndex_] & BackgroundBlack) != 0);
            backgroundWhiteACT_->setChecked((showType_[selectedRenderPortIndex_] & BackgroundWhite) != 0);
            backgroundCheckerboardPatternACT_->setChecked((showType_[selectedRenderPortIndex_] & CheckerboardPattern) != 0);

            if (rt->getColorTexture()) {
                colorRGBAACT_->setEnabled(true);
                alphaOnlyACT_->setEnabled(true);
                hOnlyACT_->setEnabled(true);
                sOnlyACT_->setEnabled(true);
                vOnlyACT_->setEnabled(true);
            }

            if (rt->getDepthTexture())
                depthOnlyACT_->setEnabled(true);

            if ((showType_[selectedRenderPortIndex_] & Color) != 0)
                colorRGBAACT_->setChecked(true);
            if (showType_[selectedRenderPortIndex_] & Alpha)
                alphaOnlyACT_->setChecked(true);
            if (showType_[selectedRenderPortIndex_] & Depth)
                depthOnlyACT_->setChecked(true);
            if (showType_[selectedRenderPortIndex_] & H)
                hOnlyACT_->setChecked(true);
            if (showType_[selectedRenderPortIndex_] & S)
                sOnlyACT_->setChecked(true);
            if (showType_[selectedRenderPortIndex_] & V)
                vOnlyACT_->setChecked(true);
        }
        else {
            colorRGBAACT_->setVisible(false);
            alphaOnlyACT_->setVisible(false);
            depthOnlyACT_->setVisible(false);
            hOnlyACT_->setVisible(false);
            sOnlyACT_->setVisible(false);
            vOnlyACT_->setVisible(false);
        }

        showInfosDetailsACT_->setVisible(maximizeOnePort_);
        saveScreenshotACT_->setVisible(maximizeOnePort_);
        saveScreenshotWithOverlayACT_->setVisible(maximizeOnePort_);

        QAction* currentAction = contextMenuMEN_->exec(e->globalPos());
        if (currentAction == saveScreenshotACT_)
            takeScreenshot(false);
        else if (currentAction == saveScreenshotWithOverlayACT_)
            takeScreenshot(true);

        if (colorRGBAACT_->isChecked())
            showType_[selectedRenderPortIndex_] =
            Color
            | (colorRACT_->isChecked() ? R : 0)
            | (colorGACT_->isChecked() ? G : 0)
            | (colorBACT_->isChecked() ? B : 0)
            | (colorAACT_->isChecked() ? A : 0);
        else if (alphaOnlyACT_->isChecked())
            showType_[selectedRenderPortIndex_] = Alpha;
        else if (depthOnlyACT_->isChecked())
            showType_[selectedRenderPortIndex_] = Depth;
        else if (hOnlyACT_->isChecked())
            showType_[selectedRenderPortIndex_] = H;
        else if (sOnlyACT_->isChecked())
            showType_[selectedRenderPortIndex_] = S;
        else if (vOnlyACT_->isChecked())
            showType_[selectedRenderPortIndex_] = V;


        if (backgroundBlackACT_->isChecked())
            showType_[selectedRenderPortIndex_] |= BackgroundBlack;
        if (backgroundWhiteACT_->isChecked())
            showType_[selectedRenderPortIndex_] |= BackgroundWhite;
        if (backgroundCheckerboardPatternACT_->isChecked())
            showType_[selectedRenderPortIndex_] |= CheckerboardPattern;
    }

    if (maximizeOnePort_ && e->button() == Qt::MidButton) {
        zoomScale_ = 1.0f;
    }

    if (e->button() == Qt::LeftButton) {
        maximizeOnePort_ = !maximizeOnePort_;
        updateSelected();
    }

    canvas_->update();
}

void RenderTargetViewer::mouseReleaseEvent(QMouseEvent* /*e*/) {
}

void RenderTargetViewer::mouseMoveEvent(QMouseEvent* e) {
    updateMousePosition(e);
    mouseIsInside_ = canvas_->underMouse();
    if (!maximizeOnePort_)
        updateSelected();
    canvas_->update();
}

void RenderTargetViewer::wheelEvent(QWheelEvent* e) {

    if(zoomScale_ == 1.0) {
        zoomOffsetX_ = mouseX_;
        zoomOffsetY_ = mouseY_;
        zoomMouseX_ = -1.0;
        zoomMouseY_ = -1.0;
    }

    if(zoomMouseX_ != mouseX_ || zoomMouseY_ != mouseY_) {
        zoomOffsetX_ += (mouseX_ - zoomOffsetX_) / zoomScale_;
        zoomOffsetY_ += (mouseY_ - zoomOffsetY_) / zoomScale_;
        zoomMouseX_ = mouseX_;
        zoomMouseY_ = mouseY_;
    }
    zoomScale_ = std::max(1.0f, zoomScale_ + (e->delta()/1200.0f));

    zoomTranslateX_ = -((zoomScale_*zoomOffsetX_)-mouseX_)/zoomScale_;
    zoomTranslateY_ = -((zoomScale_*zoomOffsetY_)-mouseY_)/zoomScale_;

    canvas_->update();
}

void RenderTargetViewer::closeEvent(QCloseEvent *event) {
    event->accept();
}

void RenderTargetViewer::keyPressEvent(QKeyEvent* e) {
    if (maximizeOnePort_) {
        if (e->key() == Qt::Key_Plus)
            zoomScale_ += 0.5;
        if (e->key() == Qt::Key_Minus)
            zoomScale_ = max(1.0, zoomScale_ - 0.5);
        if (e->key() == Qt::Key_Escape) {
            zoomScale_ = 1;
            maximizeOnePort_ = false;
        }
    }
    canvas_->update();
}

//-----------------------------------------------
//  Event Helper
//-----------------------------------------------
void RenderTargetViewer::updateMousePosition(QMouseEvent* e) {
    QPoint pos = canvas_->mapFromParent(e->pos());
    mouseX_ = pos.x();
    mouseY_ = canvas_->height() - pos.y();
}

void RenderTargetViewer::updateSelected() {
    //set default value
    selectedRenderPortIndex_ = -1;

    // return if no evaluator is present
    if (!evaluator_) return;

    //get all render ports
    std::vector<RenderPort*> renderPorts = collectRenderPorts();
    // return, if no render ports are collected
    if (renderPorts.empty()) return;
    if (mouseX_ < 0 || mouseY_ < 0) return;

    // update port grid
    int countX = static_cast<int>(std::ceil(std::sqrt(static_cast<float>(renderPorts.size()))));
    int countY = static_cast<int>(std::ceil(static_cast<float>(renderPorts.size()) / countX));

    int pixelPerPortWidth = static_cast<int>(std::ceil(static_cast<float>(canvas_->getSize().x) / countX));
    int pixelPerPortHeight = static_cast<int>(std::ceil(static_cast<float>(canvas_->getSize().y) / countY));

    size_t index = mouseX_/pixelPerPortWidth + (countY - 1 - mouseY_/pixelPerPortHeight) * countX;

    // set right value
    if (index < renderPorts.size())
        selectedRenderPortIndex_ = static_cast<int>(index);
}


} // namespace voreen
