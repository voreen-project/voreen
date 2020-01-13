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

#include "voreen/qt/mainwindow/menuentities/screenshotmenuentity/canvasrendererscreenshotplugin.h"

#include "modules/core/processors/output/canvasrenderer.h" //< core module is always available
#include "modules/core/qt/processor/canvasrendererwidget.h"

#include <QSpinBox>
#include <QComboBox>
#include <QApplication>
#include <QMessageBox>

namespace voreen {

CanvasRendererScreenshotPlugin::CanvasRendererScreenshotPlugin(QWidget* parent, Workspace* workspace, CanvasRenderer* canvasRenderer)
    : ScreenshotPluginBase(parent, workspace)
    , canvasRenderer_(canvasRenderer)
{
    //override title and prefix
    if (canvasRenderer_) {
        setWindowTitle(QString::fromStdString("Screenshot: " + canvasRenderer_->getGuiName()));
        metaDataPrefix_ = canvasRenderer_->getGuiName();
    }
    //load meta data
    loadMetaData();
    //connect native
    connect(static_cast<CanvasRendererWidget*>(canvasRenderer->getProcessorWidget()),SIGNAL(canvasRendererWidget_sizeChanged_Signal()),this,SLOT(nativeSizeHasChangedSlot()));
}

void CanvasRendererScreenshotPlugin::nativeSizeHasChangedSlot() {
    userDefinedWidthSpinBox_->blockSignals(true);
    userDefinedHeightSpinBox_->blockSignals(true);

    // if native
    if (canvasRenderer_ && resolutionComboBox_->currentIndex() == 0) { //"native"
        tgt::ivec2 size = canvasRenderer_->getCanvas()->getSize();
        userDefinedWidthSpinBox_->setValue(size.x);
        userDefinedHeightSpinBox_->setValue(size.y);
    }

    userDefinedWidthSpinBox_->blockSignals(false);
    userDefinedHeightSpinBox_->blockSignals(false);
}

void CanvasRendererScreenshotPlugin::saveScreenshot(const QString& filename, int width, int height) {
    if (!canvasRenderer_)
        return;

    tgt::ivec2 size(width, height);

    QApplication::setOverrideCursor(QCursor(Qt::WaitCursor));
    qApp->processEvents();
    try {
        bool success = canvasRenderer_->renderToImage(filename.toStdString(), size);
        QApplication::restoreOverrideCursor();
        if (!success) {
            QMessageBox::warning(this, tr("Error saving snapshot"),
                                 tr("Screenshot could not be saved:\n%1").arg(
                                     QString::fromStdString(canvasRenderer_->getRenderToImageError())));
        }
    }
    catch (const std::exception& e) {
        QApplication::restoreOverrideCursor();
        QMessageBox::warning(this, tr("Error saving snapshot"),
                             tr("Screenshot could not be saved:\n%1").arg(e.what()));
    }
}

void CanvasRendererScreenshotPlugin::canvasRendererHasBeenRenamed() {
    if (!canvasRenderer_)
        return; //should not happen
    removeMetaData();
    setWindowTitle(QString::fromStdString("Screenshot: " + canvasRenderer_->getGuiName()));
    metaDataPrefix_ = canvasRenderer_->getGuiName();
    saveMetaData();
}

} // namespace voreen
