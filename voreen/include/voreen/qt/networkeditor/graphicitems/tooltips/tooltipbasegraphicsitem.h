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

#ifndef VRN_TOOLTIPBASEGRAPHICSITEM_H
#define VRN_TOOLTIPBASEGRAPHICSITEM_H

#include "../nwebasegraphicsitem.h"
#include "../../editor_settings.h"

class QGraphicsTextItem;
class QTimer;

namespace voreen {

/**
 * The basic tooltip class. Each NWEBaseGraphicsItem can have its own toolip inherited from this class.
 */
class ToolTipBaseGraphicsItem : public NWEBaseGraphicsItem {
Q_OBJECT
public:
    /**
     * Constructor
     * @param item the graphicsitem belonging to this tooltip.
     */
    ToolTipBaseGraphicsItem(NWEBaseGraphicsItem* item);

protected:
//-- member --
    QGraphicsTextItem* toolTipText_;            //< the textitem used in the tooltip
    QTimer* timer_;                             //< the timer to show the tooltip with delay
    QPointF toolTipTimerTriggertMousePosition_; //< the mouse position when triggering the tooltip event
//---------------------------------------------------------------------------------
//              nwebasegraphicsitem functions
//---------------------------------------------------------------------------------
public:
    /**
     * At the moment tooltips dont change during layer changes.
     */
    virtual void updateNWELayerAndCursor(){};
//handling child elements
    virtual void layoutChildItems(){};
protected:
    virtual void createChildItems(){};
    virtual void deleteChildItems(){};

//---------------------------------------------------------------------------------
//              tooltipgraphicsitem functions
//---------------------------------------------------------------------------------
public:
    /**
     * Upadates the position of the tooltip.
     * @param x position of the mouse cursor.
     * @param y position of the mouse cursor.
     * @return should only set visible, if returns true.
     */
    virtual bool updateToolTip(qreal x, qreal y);
//-- getter and setter --
    /// retruns the textitem.
    const QGraphicsTextItem* getToolTipTextItem() const;
    /// sets the mouse position used to set the tooltip position.
    void setToolTipTimerTriggertMousePosition(QPointF& point);
//-- timer --
    /// starts the timer
    void startTimer();
    /// stops the timer
    void stopTimer();
protected slots:
    /**
     * On timer timeout the position of the tooltip will be set and the tooltip will become visible.
     */
    void timerTimeout();
};

} //namespace voreen

#endif // VRN_TOOLTIPBASEGRAPHICSITEM_H

