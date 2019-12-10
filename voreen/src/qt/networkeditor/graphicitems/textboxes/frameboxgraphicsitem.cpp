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

#include "voreen/qt/networkeditor/graphicitems/textboxes/frameboxgraphicsitem.h"
#include "voreen/qt/networkeditor/styles/nwestyle_base.h"

#include <QAction>

namespace voreen {

FrameBoxGraphicsItem::FrameBoxGraphicsItem(NetworkEditor* nwe)
    : TextBoxBaseGraphicsItem(nwe)
{
    currentSize_ = QSizeF(150.f,110.f);
    minimalSize_ = QSizeF(100.f,100.f);
    baseColor_ = Qt::gray;

    //no content
    delete contentItem_;
    contentItem_ = 0;
    contentEditor_ = 0; // is deletet as child of contentItem

    //add actions after child creation
    setContextMenuActions();

    //adjust children again
    captionItem_->setPlainText("Right Click -> Rename");
    captionItem_->setDefaultTextColor(fontColor_);
    captionChangedSlot();
}

FrameBoxGraphicsItem::~FrameBoxGraphicsItem() {
}

void FrameBoxGraphicsItem::initializePaintSettings() {
    return currentStyle()->FrameBoxGI_initializePaintSettings(this);
}

void FrameBoxGraphicsItem::mainPaint(QPainter* painter, const QStyleOptionGraphicsItem* option, QWidget* widget, NWEItemSettings& setting) {
    return currentStyle()->FrameBoxGI_paint(this,painter,option,widget,setting);
}

void FrameBoxGraphicsItem::setContextMenuActions() {
    TextBoxBaseGraphicsItem::setContextMenuActions();
    //visible
    toggleCaptionAction_ = new QAction(tr("Show/Hide Header"), this);
    toggleCaptionAction_->setCheckable(true);
    toggleCaptionAction_->setChecked(showCaption_);
    connect(toggleCaptionAction_,SIGNAL(triggered()),this,SLOT(toggleCaptionSlot()));
    customContextMenu_->addAction(toggleCaptionAction_);
    customContextMenu_->addSeparator();

    //rename
    renameCaptionAction_ = new QAction(QIcon(":/qt/icons/rename.png"), tr("Rename"), this);
    connect(renameCaptionAction_,SIGNAL(triggered()),this,SLOT(renameCaptionSlot()));
    renameCaptionAction_->setVisible(showCaption_);
    customContextMenu_->addAction(renameCaptionAction_);

    //change color of font
    QPixmap pixmap(32,32);

    pixmap.fill(fontColor_);
    fontColorAction_ = new QAction(QIcon(pixmap), tr("Set Font Color ..."), this);
    connect(fontColorAction_,SIGNAL(triggered()),this,SLOT(changeFontColorSlot()));
    fontColorAction_->setVisible(showCaption_);
    customContextMenu_->addAction(fontColorAction_);
        //font size
    QMenu* fontMenu = customContextMenu_->addMenu(tr("Font Size"));

    delete decreaseCaptionFontSizeAction_;
    decreaseCaptionFontSizeAction_ = new QAction(/*QIcon(":/qt/icons/document-icon.png"),*/ tr("Smaller"), this);
    connect(decreaseCaptionFontSizeAction_,SIGNAL(triggered()),this,SLOT(decreaseCaptionFontSizeSlot()));
    fontMenu->addAction(decreaseCaptionFontSizeAction_);

    delete increaseCaptionFontSizeAction_;
    increaseCaptionFontSizeAction_ = new QAction(/*QIcon(":/qt/icons/document-icon.png"),*/ tr("Larger"), this);
    connect(increaseCaptionFontSizeAction_,SIGNAL(triggered()),this,SLOT(increaseCaptionFontSizeSlot()));
    fontMenu->addAction(increaseCaptionFontSizeAction_);

    customContextMenu_->addSeparator();

    pixmap.fill(baseColor_);
    baseColorAction_ = new QAction(QIcon(pixmap), tr("Set Frame Color ..."), this);
    connect(baseColorAction_,SIGNAL(triggered()),this,SLOT(changeBaseColorSlot()));
    customContextMenu_->addAction(baseColorAction_);
}

} // namespace
