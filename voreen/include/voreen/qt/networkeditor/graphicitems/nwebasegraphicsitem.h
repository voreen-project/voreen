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

#ifndef VRN_NWEBASEGRAPHICSITEM_H
#define VRN_NWEBASEGRAPHICSITEM_H

#include "../networkeditor.h"

#include <QGraphicsItem>

#include <QMenu>

class QAction;
class QGraphicsDropShadowEffect;

namespace voreen {

class NWEStyle_Base;
struct NWEItemSettings;
class ToolTipBaseGraphicsItem;

/**
 * Base GraphicsItem used for all NetworkEditor elements. This class is abstract.
 * It provides basic functions for context menu and tooltip handling.
 * Functions used for painting and arranging child items have to be reimplemented in the sub-classes.
 */
class NWEBaseGraphicsItem : public QGraphicsObject {
Q_OBJECT
public:
    /**
     * Contructor
     * @param nwe the networkeditor used to visualize the item. MUST NOT BE ZERO.
     */
    NWEBaseGraphicsItem(NetworkEditor* nwe);
    virtual ~NWEBaseGraphicsItem();

    void enableShadows(bool enable);

//-- member --
protected:
    NetworkEditor* networkEditor_;           //< the NetworkEditor assigned in the constructor
    bool paintHasBeenInitialized_;           //< stores, if the paint settings have been initialized
    QMenu* customContextMenu_;               //< custom context menu of this item
    ToolTipBaseGraphicsItem* toolTipItem_;   //< the tooltip graphicsitem if defined
    QGraphicsDropShadowEffect* shadowEffect_;//< optional shadow effect

//-------------------------------------------------------------------------------------
//                  getter and setter
//-------------------------------------------------------------------------------------
public:
    /// returns the NetworkEditor
    NetworkEditor* getNetworkEditor() const;
    ///returns the current layer of the NetworkEditor
    NetworkEditorLayer currentLayer() const;
    ///returns the current cursor mode of the NetworkEditor
    NetworkEditorCursorMode currentCursorMode() const;
    ///returns the current style used in the NetworkEditor
    NWEStyle_Base* currentStyle() const;
    ///returns the current tooltip mode of the NetworkEditor
    bool currentToolTipMode() const;
//---paint getter and setter --
    ///returns, if the paint settings have been initialized
    bool isPaintInitialized() const;
    ///sets the paint settings to not initialized
    void resetPaintInitialization();
//-- tooltips --
    ///set a new tooltip graphicsitem
    void setToolTipGraphicsItem(ToolTipBaseGraphicsItem* tooltip);
    ///returns the tooltip graphicsitem
    ToolTipBaseGraphicsItem* getToolTipGraphicsItem() const;
//-- help function --
    ///just makes this function public
    virtual void prepareGeometryChange();

//-------------------------------------------------------------------------------------
//                  paint
//-------------------------------------------------------------------------------------
public:
    /**
     * This funtion is called by qt to paint the graphicsitem. If 'paintHasBeenInitialized_' is false,
     * initializePaintSettings() is called.
     */
    void paint(QPainter* painter, const QStyleOptionGraphicsItem* option, QWidget* widget);
    /**
     * Gets the NWEItemSettings from this.
     */
    NWEItemSettings beforePaint(QPainter* painter, const QStyleOptionGraphicsItem* option);
//-------------------------------------------------------------------------------------
//                  context menu
//-------------------------------------------------------------------------------------
public:
    /**
     * This function adds all actions in 'contextMenuActionList_' to the parameter menu.
     * The list can be updated by 'setContextMenuActions'.
     * @param menu the context menu the actions will be added to.
     */
    void addActionsToContextMenu(QMenu* menu);
protected:
    /**
     * This function must be reimplemented by each sub class. It should add all GraphicsItem specific actions to
     * 'contextMenuActionList_'. The default implementation only deletes the old list.
     */
    virtual void setContextMenuActions();

//-------------------------------------------------------------------------------------
//                  events (tooltips)
//-------------------------------------------------------------------------------------
protected:
    /**
     * This function starts the tooltip timer if the 'currentToolTipMode' is true.
     */
    virtual void hoverEnterEvent (QGraphicsSceneHoverEvent* event);
    /**
     * This function hides the tooltip if the 'currentToolTipMode' is false.
     */
    virtual void mouseMoveEvent (QGraphicsSceneMouseEvent* event);
    /**
     * This function hides the tooltip.
     */
    virtual void hoverLeaveEvent (QGraphicsSceneHoverEvent* event);

//-------------------------------------------------------------------------------------
//                  NEEDS TO BE IMPLEMENTED
//-------------------------------------------------------------------------------------
public:
    ///works in Qt like 'getClassName'. The different types are defined in editor_settings.h.
    virtual int type() const = 0;
    /**
     * Should update the state of the GraphicsItem, i.e. selectable etc., based on the current
     * layer/cursor mode of the NetworkEditor.
     */
    virtual void updateNWELayerAndCursor() = 0;
//-- paint and style --
    /**
     * Is used by Qt. Should call currentStyle()->THIS_boundingRect();
     */
    virtual QRectF boundingRect() const = 0;
    /**
     * Is used by Qt. Should call currentStyle()->THIS_shape();
     */
    virtual QPainterPath shape() const = 0;
protected:
    /**
     * Is called in the first paint and only ones. Should call currentStyle()->THIS_initializePaintSettings();
     */
    virtual void initializePaintSettings() = 0;
    /**
     * Is used in NWEBaseGraphicsItem::paint(). Should call currentStyle()->THIS_paint();
     */
    virtual void mainPaint(QPainter* painter, const QStyleOptionGraphicsItem* option, QWidget* widget, NWEItemSettings& setting) = 0;
//-- handling  child elements --
public:
    /**
     * Arrange the child items.
     */
    virtual void layoutChildItems() = 0;
protected:
    /**
     * Creates the child items.
     */
    virtual void createChildItems() = 0;
    /**
     * Deletes the child items.
     */
    virtual void deleteChildItems() = 0;
};

} //namespace voreen

#endif // VRN_NWEBASEGRAPHICSITEM_H

