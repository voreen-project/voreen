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

#include "voreen/qt/networkeditor/graphicitems/textboxes/textboxbasegraphicsitem.h"
#include "voreen/qt/networkeditor/meta/textboxmetadata.h"
#include "voreen/qt/networkeditor/styles/nwestyle_base.h"

#include "voreen/qt/networkeditor/graphicitems/utils/renamabletextgraphicsitem.h"

#include <QApplication>
#include <QGraphicsSceneMouseEvent>
#include <QAction>
#include <QColorDialog>

namespace voreen {

TextBoxBaseGraphicsItem::TextBoxBaseGraphicsItem(NetworkEditor* nwe)
    : NWEBaseGraphicsItem(nwe)
    //resize
    , resizeDirection_(RD_NONE)
    //caption
    , toggleCaptionAction_(0), renameCaptionAction_(0), fontColorAction_(0), baseColorAction_(0)
    , decreaseCaptionFontSizeAction_(0), increaseCaptionFontSizeAction_(0)
    //content
    , contentEditDisabled_(true)
    , lastEditStyleSheet_("")
    //persist members
    , showCaption_(true)
    , currentSize_(150.f,150.f)
    , minimalSize_(100.f,100.f)
    , resizeBorder_(3.f)
    , fontColor_(Qt::black)
    , baseColor_(nwe->getCurrentStyle()->NWEStyle_TextBoxBaseMainColor)
    , captionItem_(0)
    , contentItem_(0), contentEditor_(0), decreaseContentFontSizeAction_(0), increaseContentFontSizeAction_(0)
    {
    tgtAssert(nwe, "passed null pointer");

    setFlag(ItemIsSelectable);
    setFlag(ItemIsMovable,true);
    //setFlag(ItemIsPanel);
    setAcceptHoverEvents(true);
    setZValue(ZValuesTextBoxBaseGraphicsItem);

    createChildItems();
}

TextBoxBaseGraphicsItem::~TextBoxBaseGraphicsItem() {
    deleteChildItems();
}

void TextBoxBaseGraphicsItem::copyFromMeta(TextBoxMetaData* meta) {
    tgtAssert(meta,"no meta data passed");
    prepareGeometryChange();
    currentSize_.setWidth(meta->getWidth());
    currentSize_.setHeight(meta->getHeight());
    fontColor_ = meta->getFontColor();
    if(fontColorAction_) {
        QPixmap pmap(32,32);
        pmap.fill(fontColor_);
        fontColorAction_->setIcon(QIcon(pmap));
    }
    baseColor_ = meta->getBaseColor();
    if(baseColorAction_) {
        QPixmap pmap(32,32);
        pmap.fill(baseColor_);
        baseColorAction_->setIcon(QIcon(pmap));
    }
    if(captionItem_) {
        captionItem_->setDefaultTextColor(fontColor_);
        if(meta->getCaptionFontSize() != -1) {
            QFont font = captionItem_->font();
            font.setPointSize(meta->getCaptionFontSize());
            captionItem_->setFont(font);
        }
        captionItem_->setPlainText(meta->getCaption().c_str());
    }
    if(contentEditor_) {
        if(meta->getContentFontSize() != -1) {
            QFont font = contentEditor_->font();
            font.setPointSize(meta->getContentFontSize());
            contentEditor_->setFont(font);
        }
        contentEditor_->setPlainText(meta->getContent().c_str());
    }
    showCaption_ = meta->getShowCaption();
    if(toggleCaptionAction_) {
        toggleCaptionAction_->setChecked(showCaption_);
        toggleCaptionSlot();
    }
    setPos(QPointF(meta->getX(),meta->getY()));
}

//---------------------------------------------------------------------------------------------------------------
//                  getter and setter functions
//---------------------------------------------------------------------------------------------------------------
void TextBoxBaseGraphicsItem::setPos(const QPointF & pos) {
    NWEBaseGraphicsItem::setPos(pos);
    layoutChildItems();
}

QSizeF TextBoxBaseGraphicsItem::getCurrentSize() const {
    return currentSize_;
}

qreal TextBoxBaseGraphicsItem::getResizeBorder() const {
    return resizeBorder_;
}

QColor TextBoxBaseGraphicsItem::getFontColor() const {
    return fontColor_;
}

QColor TextBoxBaseGraphicsItem::getBaseColor() const {
    return baseColor_;
}

bool TextBoxBaseGraphicsItem::isEditDisabled() const {
    return contentEditDisabled_;
}
bool TextBoxBaseGraphicsItem::isCaptionEnabled() const {
    return showCaption_;
}


const QGraphicsTextItem* TextBoxBaseGraphicsItem::getCaptionItem() const {
    return captionItem_;
}

QTextEdit* TextBoxBaseGraphicsItem::getContentEditor() const {
    return contentEditor_;
}

//---------------------------------------------------------------------------------------------------------------
//                  nwebasegraphicsitem functions
//---------------------------------------------------------------------------------------------------------------
QRectF TextBoxBaseGraphicsItem::boundingRect() const {
    return currentStyle()->TextBoxBaseGI_boundingRect(this);
}
QPainterPath TextBoxBaseGraphicsItem::shape() const {
    return currentStyle()->TextBoxBaseGI_shape(this);
}

void TextBoxBaseGraphicsItem::layoutChildItems() {
    prepareGeometryChange();
    if(captionItem_) {
        captionItem_->setPos(resizeBorder_,resizeBorder_);
        if(contentItem_)
            contentItem_->setPos(resizeBorder_,3*resizeBorder_+captionItem_->boundingRect().height());
    } else {
        if(contentItem_)
            contentItem_->setPos(resizeBorder_,resizeBorder_);
    }
    if(contentEditor_) {
        contentEditor_->setMaximumSize(currentSize_.width()-2*resizeBorder_,currentSize_.height()-2*resizeBorder_);
        contentEditor_->setMinimumSize(currentSize_.width()-2*resizeBorder_,currentSize_.height()-2*resizeBorder_);
        contentEditor_->update();
    }
}

void TextBoxBaseGraphicsItem::createChildItems() {
    //caption items
    captionItem_ = new RenamableTextGraphicsItem("Caption",this);
    captionItem_->setParentItem(this);
    captionItem_->setAcceptHoverEvents(false);
    QFont font = captionItem_->font();
    font.setPointSize(9);
    captionItem_->setFont(font);
    connect(captionItem_, SIGNAL(renameFinished()), this, SLOT(renameCaptionFinishedSlot()));
    connect(captionItem_, SIGNAL(textChanged()), this, SLOT(captionChangedSlot()));
    //content items
    contentItem_ = new QGraphicsProxyWidget(this);
    contentItem_->setAcceptHoverEvents(true);
    contentEditor_ = new QTextEdit(nullptr);
    contentEditor_->setMouseTracking(true);
    contentEditor_->setFrameStyle(QFrame::NoFrame);
    contentEditor_->installEventFilter(this);
    font = contentEditor_->font();
    font.setPointSize(6);
    contentEditor_->setFont(font);
    contentItem_->setWidget(contentEditor_);
    contentItem_->setOpacity(0.99999); // HACK: Workaround for Qt5 bug: QTBUG-55070

    //update editor
    contentEditDisabled_ = false; //is switched back in next call
    switchContentEditModeSlot();
    layoutChildItems();
}

void TextBoxBaseGraphicsItem::deleteChildItems() {
    delete captionItem_;
    captionItem_ = 0;
    delete contentItem_;
    contentItem_ = 0;
    contentEditor_ = 0;
}

//---------------------------------------------------------------------------------------------------------------
//                  resize event functions
//---------------------------------------------------------------------------------------------------------------
void TextBoxBaseGraphicsItem::mousePressEvent(QGraphicsSceneMouseEvent * event) {
    if(contentEditDisabled_) {
        resizeDirection_ = determineResizeDirection(event->pos());
        if(resizeDirection_ != RD_NONE) {
            setFlag(ItemIsMovable,false);
        }
    }
    NWEBaseGraphicsItem::mousePressEvent(event);
}

void TextBoxBaseGraphicsItem::mouseMoveEvent (QGraphicsSceneMouseEvent* event){
    if(resizeDirection_ != RD_NONE) {
        prepareGeometryChange();
        qreal minWidth = qMax(minimalSize_.width(),(captionItem_ ? captionItem_->boundingRect().width()+2*resizeBorder_ : 0.f));
        qreal minHeight = qMax(minimalSize_.height(),(captionItem_ ? captionItem_->boundingRect().height()+6*resizeBorder_ : 0.f));
        switch(resizeDirection_){
        case RD_LEFT: {
            qreal shift = (currentSize_.width()-event->pos().x() < minWidth ? 0.f : event->pos().x());
            currentSize_.setWidth(currentSize_.width()-shift);
            setPos(pos()+(mapFromItem(this,QPointF(shift,0.f))));
                      } break;
        case RD_RIGHT:
            currentSize_.setWidth(qMax(minWidth,event->pos().x()));
            break;
        case RD_UP: {
            qreal shift = (currentSize_.height()-event->pos().y() < minHeight ? 0.f : event->pos().y());
            currentSize_.setHeight(currentSize_.height()-shift);
            setPos(pos()+(mapFromItem(this,QPointF(0.f,shift))));
                      } break;
            break;
        case RD_DOWN:
            currentSize_.setHeight(qMax(minHeight,event->pos().y()));
            break;
        case RD_LEFT_UP: {
            qreal shiftX = (currentSize_.width()-event->pos().x() < minWidth ? 0.f : event->pos().x());
            qreal shiftY = (currentSize_.height()-event->pos().y() < minHeight ? 0.f : event->pos().y());
            currentSize_ = currentSize_ - QSizeF(shiftX,shiftY);
            setPos(pos()+(mapFromItem(this,QPointF(shiftX,shiftY))));
                         } break;
        case RD_LEFT_DOWN: {
            qreal shiftX = (currentSize_.width()-event->pos().x() < minWidth ? 0.f : event->pos().x());
            currentSize_.setWidth(currentSize_.width()-shiftX);
            currentSize_.setHeight(qMax(minHeight,event->pos().y()));
            setPos(pos()+(mapFromItem(this,QPointF(shiftX,0.f))));
                           } break;
        case RD_RIGHT_UP: {
            qreal shiftY = (currentSize_.height()-event->pos().y() < minHeight ? 0.f : event->pos().y());
            currentSize_.setHeight(currentSize_.height()-shiftY);
            currentSize_.setWidth(qMax(minWidth,event->pos().x()));
            setPos(pos()+(mapFromItem(this,QPointF(0.f,shiftY))));
                          } break;
        case RD_RIGHT_DOWN:
            currentSize_.setWidth(qMax(minWidth,event->pos().x()));
            currentSize_.setHeight(qMax(minHeight,event->pos().y()));
            break;
        default:
            tgtAssert(false,"should not get here!");
        }
        layoutChildItems();
    } else {
        //mouse is pressed
        if(!QApplication::overrideCursor() || QApplication::overrideCursor()->shape() != Qt::SizeAllCursor) {
            QApplication::restoreOverrideCursor(); //propably not needed
            QApplication::setOverrideCursor(Qt::SizeAllCursor);
        }
    }
    NWEBaseGraphicsItem::mouseMoveEvent(event);
}

void TextBoxBaseGraphicsItem::mouseReleaseEvent(QGraphicsSceneMouseEvent * event) {
    if(resizeDirection_ != RD_NONE) {
        setFlag(ItemIsMovable,true);
        resizeDirection_ = RD_NONE;
    }
    QApplication::restoreOverrideCursor();
    NWEBaseGraphicsItem::mouseReleaseEvent(event);
}

void TextBoxBaseGraphicsItem::hoverEnterEvent (QGraphicsSceneHoverEvent* event){
    NWEBaseGraphicsItem::hoverEnterEvent(event);
}

void TextBoxBaseGraphicsItem::hoverMoveEvent (QGraphicsSceneHoverEvent* event){
    if(contentEditDisabled_) {
        switch(determineResizeDirection(event->pos())) {
        case RD_LEFT:
        case RD_RIGHT:
            if(!QApplication::overrideCursor() || QApplication::overrideCursor()->shape() != Qt::SizeHorCursor) {
                QApplication::restoreOverrideCursor();
                QApplication::setOverrideCursor(Qt::SizeHorCursor);
            }
            break;
        case RD_UP:
        case RD_DOWN:
            if(!QApplication::overrideCursor() || QApplication::overrideCursor()->shape() != Qt::SizeVerCursor) {
                QApplication::restoreOverrideCursor();
                QApplication::setOverrideCursor(Qt::SizeVerCursor);
            }
            break;
        case RD_LEFT_UP:
        case RD_RIGHT_DOWN:
            if(!QApplication::overrideCursor() || QApplication::overrideCursor()->shape() != Qt::SizeFDiagCursor) {
                QApplication::restoreOverrideCursor();
                QApplication::setOverrideCursor(Qt::SizeFDiagCursor);
            }
            break;
        case RD_RIGHT_UP:
        case RD_LEFT_DOWN:
            if(!QApplication::overrideCursor() || QApplication::overrideCursor()->shape() != Qt::SizeBDiagCursor) {
                QApplication::restoreOverrideCursor();
                QApplication::setOverrideCursor(Qt::SizeBDiagCursor);
            }
            break;
        default:
            QApplication::restoreOverrideCursor();
            break;
        }
    }
    NWEBaseGraphicsItem::hoverMoveEvent(event);
}

void TextBoxBaseGraphicsItem::hoverLeaveEvent (QGraphicsSceneHoverEvent* event){
    //restore cursor
    QApplication::restoreOverrideCursor();
    NWEBaseGraphicsItem::hoverLeaveEvent(event);
}

TextBoxBaseGraphicsItem::ResizeDirection TextBoxBaseGraphicsItem::determineResizeDirection(QPointF pos) {
    if(pos.x() <= resizeBorder_) {
        if(pos.y() <= resizeBorder_)
            return RD_LEFT_UP;
        else if(pos.y() >= currentSize_.height()-resizeBorder_)
            return RD_LEFT_DOWN;
        else
            return RD_LEFT;
    }
    else if(pos.x() >= currentSize_.width()-resizeBorder_) {
        if(pos.y() <= resizeBorder_)
            return RD_RIGHT_UP;
        else if(pos.y() >= currentSize_.height()-resizeBorder_)
            return RD_RIGHT_DOWN;
        else
            return RD_RIGHT;
    }
    else if(pos.y() <= resizeBorder_)
        return RD_UP;
    else if(pos.y() >= currentSize_.height()-resizeBorder_)
        return RD_DOWN;
    else
        return RD_NONE;
}

bool TextBoxBaseGraphicsItem::eventFilter(QObject* watched, QEvent* event) {
    if(watched == contentEditor_) {
        //LERRORC("Debug Event",event->type());
        if(contentEditDisabled_) {
            event->ignore();
            if(event->type() == QEvent::MouseMove) {
                QGraphicsSceneHoverEvent* sceneEvent = new QGraphicsSceneHoverEvent();
                sceneEvent->setPos(static_cast<QMouseEvent*>(event)->pos()+QPoint(resizeBorder_,resizeBorder_));
                hoverMoveEvent(sceneEvent);
                delete sceneEvent;
            }
            return NWEBaseGraphicsItem::eventFilter(watched,event);
        } else {
            if(event->type() == QEvent::MouseMove) {
                QApplication::restoreOverrideCursor();
            }
            if(!this->isSelected()) { //if selection is lost, stop editting
                switchContentEditModeSlot();
                return true;
            }
        }
    }
    return NWEBaseGraphicsItem::eventFilter(watched,event);
}

//---------------------------------------------------------------------------------------------------------------
//                  caption edit / show functions
//---------------------------------------------------------------------------------------------------------------
void TextBoxBaseGraphicsItem::renameCaptionSlot() {
    if(captionItem_) {
        captionItem_->setTextInteractionFlags(Qt::TextEditorInteraction);
        captionItem_->setFlag(QGraphicsItem::ItemIsFocusable, true);
        captionItem_->setFocus();
        QTextCursor cursor = captionItem_->textCursor();
        cursor.select(QTextCursor::Document);
        captionItem_->setTextCursor(cursor);
    }
}

void TextBoxBaseGraphicsItem::captionChangedSlot() {
    if(captionItem_) {
        currentSize_ = currentSize_.expandedTo(captionItem_->boundingRect().size() + QSizeF(2*resizeBorder_,6*resizeBorder_));
        layoutChildItems();
        update();
    }
}

void TextBoxBaseGraphicsItem::renameCaptionFinishedSlot() {
    if(captionItem_) {
        captionItem_->setTextInteractionFlags(Qt::NoTextInteraction);
        captionItem_->setFlag(QGraphicsItem::ItemIsFocusable, false);
        captionItem_->setFlag(QGraphicsItem::ItemIsSelectable, false);
        //update scene
        captionChangedSlot();
    }
}

void TextBoxBaseGraphicsItem::toggleCaptionSlot() {
    if(toggleCaptionAction_) {
        showCaption_ = toggleCaptionAction_->isChecked();
        if(fontColorAction_) fontColorAction_->setVisible(showCaption_);
        if(renameCaptionAction_) renameCaptionAction_->setVisible(showCaption_);
        if(captionItem_) captionItem_->setVisible(showCaption_);
        update();
    }
}

void TextBoxBaseGraphicsItem::decreaseCaptionFontSizeSlot() {
    if(captionItem_) {
        QFont font = captionItem_->font();
        int oldSize = font.pointSize();
        int newSize = tgt::iround((float)oldSize / 1.5f);
        if(newSize >= 4) {
            font.setPointSize(newSize);
            captionItem_->setFont(font);
            captionChangedSlot(); //adjust size of graphics item
            //enable/disable don't make the action grey. it is only unclickable
            if(increaseCaptionFontSizeAction_ && !increaseCaptionFontSizeAction_->isVisible())
                increaseCaptionFontSizeAction_->setVisible(true);
            if(decreaseCaptionFontSizeAction_ && (tgt::iround((float)newSize / 1.5f) < 4))
                decreaseCaptionFontSizeAction_->setVisible(false);
        } else {
            tgtAssert(false,"should not get here!");
        }
        //LERRORC("Debug",oldSize << " / " << newSize);
    }
}

void TextBoxBaseGraphicsItem::increaseCaptionFontSizeSlot() {
    if(captionItem_) {
        QFont font = captionItem_->font();
        int oldSize = font.pointSize();
        int newSize = tgt::iround((float)oldSize * 1.5f);
        if(newSize <= 50) {
            font.setPointSize(newSize);
            captionItem_->setFont(font);
            captionChangedSlot(); //adjust size of graphics item
            //enable/disable don't make the action grey. it is only unclickable
            if(decreaseCaptionFontSizeAction_ && !decreaseCaptionFontSizeAction_->isVisible())
                decreaseCaptionFontSizeAction_->setVisible(true);
            if(increaseCaptionFontSizeAction_ && (tgt::iround((float)newSize * 1.5f) > 50))
                increaseCaptionFontSizeAction_->setVisible(false);
        } else {
            tgtAssert(false,"should not get here!");
        }
        //LERRORC("Debug",oldSize << " / " << newSize);
    }
}


//---------------------------------------------------------------------------------------------------------------
//                  content edit functions
//---------------------------------------------------------------------------------------------------------------
void TextBoxBaseGraphicsItem::switchContentEditModeSlot() {
    //LERRORC("Debug","toggle edit");
    contentEditDisabled_ = !contentEditDisabled_;
    if(contentEditor_) {
        if(contentEditDisabled_) {
            contentEditor_->setReadOnly(true);
            contentEditor_->setTextInteractionFlags(Qt::NoTextInteraction);
            if(!lastEditStyleSheet_.isEmpty())
                contentEditor_->setStyleSheet(lastEditStyleSheet_);
            contentEditor_->viewport()->setCursor(QCursor(Qt::ArrowCursor));
            QTextCursor cur = contentEditor_->textCursor();
            cur.clearSelection();
            contentEditor_->setTextCursor(cur);
        } else {
            lastEditStyleSheet_ = contentEditor_->styleSheet();
            contentEditor_->setStyleSheet("background-color: white");
            contentEditor_->moveCursor(QTextCursor::End);
            contentEditor_->setReadOnly(false);
            contentEditor_->viewport()->setCursor(QCursor(Qt::IBeamCursor));
            contentEditor_->selectAll();
            //force new focus
            contentEditor_->setFocus();
        }
    }
}

void TextBoxBaseGraphicsItem::decreaseContentFontSizeSlot() {
    if(contentEditor_) {
        QFont font = contentEditor_->font();
        int oldSize = font.pointSize();
        int newSize = tgt::iround((float)oldSize / 1.5f);
        if(newSize >= 4) {
            font.setPointSize(newSize);
            contentEditor_->setFont(font);
            //enable/disable don't make the action grey. it is only unclickable
            if(increaseContentFontSizeAction_ && !increaseContentFontSizeAction_->isVisible())
                increaseContentFontSizeAction_->setVisible(true);
            if(decreaseContentFontSizeAction_ && (tgt::iround((float)newSize / 1.5f) < 4))
                decreaseContentFontSizeAction_->setVisible(false);
        } else {
            tgtAssert(false,"should not get here!");
        }
        //LERRORC("Debug",oldSize << " / " << newSize);
    }
}

void TextBoxBaseGraphicsItem::increaseContentFontSizeSlot() {
    if(contentEditor_) {
        QFont font = contentEditor_->font();
        int oldSize = font.pointSize();
        int newSize = tgt::iround((float)oldSize * 1.5f);
        if(newSize <= 50) {
            font.setPointSize(newSize);
            contentEditor_->setFont(font);
            //enable/disable don't make the action grey. it is only unclickable
            if(decreaseContentFontSizeAction_ && !decreaseContentFontSizeAction_->isVisible())
                decreaseContentFontSizeAction_->setVisible(true);
            if(increaseContentFontSizeAction_ && (tgt::iround((float)newSize * 1.5f) > 50))
                increaseContentFontSizeAction_->setVisible(false);
        } else {
            tgtAssert(false,"should not get here!");
        }
        //LERRORC("Debug",oldSize << " / " << newSize);
    }
}

//---------------------------------------------------------------------------------------------------------------
//                  color edit functions
//---------------------------------------------------------------------------------------------------------------
void TextBoxBaseGraphicsItem::changeFontColorSlot() {
    QObject* obj = sender();
    QAction* act = dynamic_cast<QAction*>(obj);
    changeColorHelper(act,true);
}

void TextBoxBaseGraphicsItem::changeBaseColorSlot() {
    QObject* obj = sender();
    QAction* act = dynamic_cast<QAction*>(obj);
    changeColorHelper(act,false);
}

void TextBoxBaseGraphicsItem::changeColorHelper(QAction* act, bool changeFontColor) {
    tgtAssert(act,"no action passed");
    QColor tmpColor; QString tmpString;
    if(changeFontColor) {
        tmpColor = fontColor_;
        tmpString = "Select Font Color";
    } else {
        tmpColor = baseColor_;
        tmpString = "Select Frame Color";
    }
    //get color
    tmpColor = QColorDialog::getColor(tmpColor,0,tmpString);
    //check, if cancel has been called
    if(tmpColor.isValid()) {
        QPixmap pixmap(32,32);
        pixmap.fill(tmpColor);
        act->setIcon(QIcon(pixmap));
        if(changeFontColor) {
            fontColor_ = tmpColor;
            captionItem_->setDefaultTextColor(fontColor_);
        } else {
            baseColor_ = tmpColor;
        }
        update();
    }
}

} // namespace
