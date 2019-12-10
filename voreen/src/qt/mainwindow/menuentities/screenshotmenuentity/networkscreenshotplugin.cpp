/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2019 University of Muenster, Germany,                        *
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

#include "voreen/qt/mainwindow/menuentities/screenshotmenuentity/networkscreenshotplugin.h"

#include "voreen/qt/networkeditor/networkeditor.h"

#include <QSpinBox>
#include <QComboBox>

namespace voreen {

NetworkScreenshotPlugin::NetworkScreenshotPlugin(QWidget* parent, Workspace* workspace, NetworkEditor* networkEditorWidget)
    : ScreenshotPluginBase(parent, workspace)
    , networkEditorWidget_(networkEditorWidget)
{
    metaDataPrefix_ = "NetworkScreenshotPlugin";
    //override title
    setWindowTitle(tr("Network Screenshot"));
    //load meta data
    loadMetaData();
    //connect
    connect(networkEditorWidget_,SIGNAL(networkEditor_visibleSceneSizeChanged_Signal()),this,SLOT(nativeSizeHasChangedSlot()));
}

void NetworkScreenshotPlugin::nativeSizeHasChangedSlot() {
    userDefinedWidthSpinBox_->blockSignals(true);
    userDefinedHeightSpinBox_->blockSignals(true);

    //if native is selected, the user defined spinboxes are set
    if(networkEditorWidget_ && resolutionComboBox_->currentIndex() == 0) { // "native"
        QRectF visibleRect;
        foreach (QGraphicsItem* item, networkEditorWidget_->scene()->items()) {
            if (item->isVisible()) {
                QRectF iRect = item->mapRectToScene(item->boundingRect());
                visibleRect = visibleRect.united(iRect);
            }
        }
        userDefinedWidthSpinBox_->setValue(visibleRect.width());
        userDefinedHeightSpinBox_->setValue(visibleRect.height());
    }

    userDefinedWidthSpinBox_->blockSignals(false);
    userDefinedHeightSpinBox_->blockSignals(false);
}

void NetworkScreenshotPlugin::saveScreenshot(const QString& filename, int width, int height) {
    if (!networkEditorWidget_)
        return;

    //determine current size
    QRectF visibleRect;
    foreach (QGraphicsItem* item, networkEditorWidget_->scene()->items()) {
        if (item->isVisible()) {
            QRectF iRect = item->mapRectToScene(item->boundingRect());
            visibleRect = visibleRect.united(iRect);
        }
    }

    //prepare pixmap (fill used to activate transparency)
    QPixmap pixmap(width, height);
    pixmap.fill(QColor(255, 255, 255, 0));

    //prepare painter
    QPainter painter(&pixmap);
    painter.setRenderHint(QPainter::Antialiasing, true);
    painter.setRenderHint(QPainter::TextAntialiasing, true);

    //render and save scene
    networkEditorWidget_->scene()->render(&painter,QRectF(),visibleRect);
    pixmap.save(filename);
}

void NetworkScreenshotPlugin::saveScreenshot(const QString& filename) {
    QRectF visibleRect;
    foreach (QGraphicsItem* item, networkEditorWidget_->scene()->items()) {
        if (item->isVisible()) {
            QRectF iRect = item->mapRectToScene(item->boundingRect());
            visibleRect = visibleRect.united(iRect);
        }
    }
    //save screenshot with native settings
    saveScreenshot(filename, visibleRect.width(), visibleRect.height());
}

} // namespace voreen
