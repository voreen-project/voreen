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

#ifndef VRN_TEXTBOXBASEGRAPHICSITEM_H
#define VRN_TEXTBOXBASEGRAPHICSITEM_H

#include "../nwebasegraphicsitem.h"

#include <QTextEdit>
#include <QGraphicsProxyWidget>

namespace voreen {

    class TextBoxMetaData;
    class RenamableTextGraphicsItem;

    /**
     * Base item of all text box items.
     * This is the super class of TextBoxGraphicsItem and FrameBoxGraphicsItem.
     * The class controls the header visibility and color controls.
     */
class TextBoxBaseGraphicsItem : public NWEBaseGraphicsItem {
Q_OBJECT
public:
    //constructor + destructor
    TextBoxBaseGraphicsItem(NetworkEditor* nwe);
    virtual ~TextBoxBaseGraphicsItem();

    /**
     * Restores all informations from the passed meta data.
     */
    void copyFromMeta(TextBoxMetaData* meta);

    //--------------------------------------------------------------------------------
    //      getter and setter functions
    //--------------------------------------------------------------------------------
    /// sets the position of the text item
    virtual void setPos(const QPointF & pos);
    /// returns the current size
    QSizeF getCurrentSize() const;
    /// returns the resize/highlight border
    qreal getResizeBorder() const;
    /// returns the font color.
    QColor getFontColor() const;
    /// returns the base color.
    QColor getBaseColor() const;
    /// checks, if the edit function is enabled
    bool isEditDisabled() const;
    /// checks, if the caption should be shown
    bool isCaptionEnabled() const;
    /// returns the header label
    const QGraphicsTextItem* getCaptionItem() const;
    /// returns the text editor
    QTextEdit* getContentEditor() const;

    //--------------------------------------------------------------------------------
    //      nwebasegraphicsitem functions
    //--------------------------------------------------------------------------------
    virtual int type() const {return UserTypesTextBoxBaseGraphicsItem;}
    virtual void updateNWELayerAndCursor() {}
    virtual QRectF boundingRect() const;
    virtual QPainterPath shape() const;
    //handling child items (no children)
public:
    void layoutChildItems();
protected:
    void createChildItems();
    void deleteChildItems();

    //--------------------------------------------------------------------------------
    //   resize events
    //--------------------------------------------------------------------------------
protected:
    virtual void mousePressEvent(QGraphicsSceneMouseEvent * event);
    virtual void mouseMoveEvent(QGraphicsSceneMouseEvent* event);
    virtual void mouseReleaseEvent(QGraphicsSceneMouseEvent * event);
    virtual void hoverEnterEvent(QGraphicsSceneHoverEvent* event);
    virtual void hoverMoveEvent(QGraphicsSceneHoverEvent* event);
    virtual void hoverLeaveEvent(QGraphicsSceneHoverEvent* event);
    /** Filter to convert events from QTextEdit to GraphicsItem */
    virtual bool eventFilter(QObject* watched, QEvent* event );
    /** Directions of the current resize event */
    enum ResizeDirection {
        RD_NONE, RD_LEFT, RD_RIGHT, RD_UP, RD_DOWN,
        RD_LEFT_UP, RD_LEFT_DOWN, RD_RIGHT_UP, RD_RIGHT_DOWN
    };
    /** Returns the resize direction depending on the mouse position */
    ResizeDirection determineResizeDirection(QPointF pos);

    ResizeDirection resizeDirection_;  //< containing the direction determined during press event

    //--------------------------------------------------------------------------------
    //   edit/show caption slots
    //--------------------------------------------------------------------------------
protected slots:
    /** Enters the caption renaming */
    void renameCaptionSlot();
    /** Adjust label during caption change */
    void captionChangedSlot();
    /** Finishes the caption renaming */
    void renameCaptionFinishedSlot();

    /** enables and disables the caption */
    void toggleCaptionSlot();

    /** Decreases the font size of the caption item. */
    void decreaseCaptionFontSizeSlot();
    /** increases the font size of the caption item. */
    void increaseCaptionFontSizeSlot();

protected:
    QAction* toggleCaptionAction_; //< member action to modify visibility in toggleCaptionSlot
    QAction* renameCaptionAction_; //< member action to modify visibility in toggleCaptionSlot
    QAction* fontColorAction_;     //< member action to modify visibility in toggleCaptionSlot
    QAction* baseColorAction_;     //< member action to modify visibility in toggleCaptionSlot

    QAction* decreaseCaptionFontSizeAction_;     ///< member action to modify the caption font size
    QAction* increaseCaptionFontSizeAction_;     ///< member action to modify the caption font size

    //--------------------------------------------------------------------------------
    //   edit content slots
    //--------------------------------------------------------------------------------
protected slots:
    /** enables / disables the text editor */
    void switchContentEditModeSlot();
    /** Decreases the font size of the content editor. */
    void decreaseContentFontSizeSlot();
    /** increases the font size of the content editor. */
    void increaseContentFontSizeSlot();

protected:
    bool contentEditDisabled_;       //< stores, if we are in edit mode
    QString lastEditStyleSheet_;     //< contains the last used editor style sheet

    QAction* decreaseContentFontSizeAction_;     ///< member action to modify the content font size
    QAction* increaseContentFontSizeAction_;     ///< member action to modify the content font size

    //--------------------------------------------------------------------------------
    //   edit color slots
    //--------------------------------------------------------------------------------
protected slots:
    /** changes the color of the font */
    void changeFontColorSlot();
    /** changes the color of a base item */
    void changeBaseColorSlot();
private:
    /** helper used in both slots */
    void changeColorHelper(QAction* act, bool changeFontColor);

    //--------------------------------------------------------------------------------
    //   member
    //--------------------------------------------------------------------------------
protected:
     bool showCaption_;     //< used to disable the caption
    //sizes
    QSizeF currentSize_;    //< current size of the text item
    QSizeF minimalSize_;    //< minimal size of the text item
    qreal resizeBorder_;    //< border size to resize the item and to draw it
    //colors
    QColor fontColor_;      //< color of the font
    QColor baseColor_;      //< base color (can be modified)
    //base items
    RenamableTextGraphicsItem* captionItem_; //< caption (header) of the text box
    QGraphicsProxyWidget* contentItem_;      //< content of the text box
    QTextEdit* contentEditor_;               //< editor contained in the contentItem
};

} // namespace

#endif // VRN_TEXTBOXBASEGRAPHICSITEM_H




