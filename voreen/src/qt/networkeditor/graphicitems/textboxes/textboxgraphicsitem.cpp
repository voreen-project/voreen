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

#include "voreen/qt/networkeditor/graphicitems/textboxes/textboxgraphicsitem.h"
#include "voreen/qt/networkeditor/styles/nwestyle_base.h"

#include <QAction>

namespace voreen {

TextBoxGraphicsItem::TextBoxGraphicsItem(NetworkEditor* nwe)
    : TextBoxBaseGraphicsItem(nwe)
{
    currentSize_ = QSizeF(121.f,35.f);
    minimalSize_ = QSizeF(100.f,30.f);

    //no caption
    delete captionItem_;
    captionItem_ = 0;
    //content
    contentEditor_->setPlainText("Press \"right click -> edit\" to modify the text note.");
    //add actions after child creation
    setContextMenuActions();
    //update items
    layoutChildItems();  

    // Setup shadows.
    resetPaintInitialization();
}

TextBoxGraphicsItem::~TextBoxGraphicsItem() {
}

void TextBoxGraphicsItem::resetPaintInitialization() {
    enableShadows(currentStyle()->getShadowsEnabled());
    TextBoxBaseGraphicsItem::resetPaintInitialization();
}
void TextBoxGraphicsItem::initializePaintSettings() {
    currentStyle()->TextBoxGI_initializePaintSettings(this);
}
void TextBoxGraphicsItem::mainPaint(QPainter* painter, const QStyleOptionGraphicsItem* option, QWidget* widget, NWEItemSettings& setting) {
    currentStyle()->TextBoxGI_paint(this,painter,option,widget,setting);
}

void TextBoxGraphicsItem::setContextMenuActions() {
    TextBoxBaseGraphicsItem::setContextMenuActions();
    QAction* editAction = new QAction(QIcon(":/qt/icons/document-icon.png"), tr("Edit"), this);
    connect(editAction,SIGNAL(triggered()),this,SLOT(switchContentEditModeSlot()));
    customContextMenu_->addAction(editAction);

    //font size
    QMenu* fontMenu = customContextMenu_->addMenu(tr("Font Size"));

    delete decreaseContentFontSizeAction_;
    decreaseContentFontSizeAction_ = new QAction(/*QIcon(":/qt/icons/document-icon.png"),*/ tr("Smaller"), this);
    connect(decreaseContentFontSizeAction_,SIGNAL(triggered()),this,SLOT(decreaseContentFontSizeSlot()));
    fontMenu->addAction(decreaseContentFontSizeAction_);

    delete increaseContentFontSizeAction_;
    increaseContentFontSizeAction_ = new QAction(/*QIcon(":/qt/icons/document-icon.png"),*/ tr("Larger"), this);
    connect(increaseContentFontSizeAction_,SIGNAL(triggered()),this,SLOT(increaseContentFontSizeSlot()));
    fontMenu->addAction(increaseContentFontSizeAction_);
}

} // namespace
