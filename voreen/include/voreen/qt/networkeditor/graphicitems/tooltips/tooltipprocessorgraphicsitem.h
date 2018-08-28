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

#ifndef VRN_TOOLTIPPROCESSORGRAPHICSITEM_H
#define VRN_TOOLTIPPROCESSORGRAPHICSITEM_H

#include "tooltipbasegraphicsitem.h"
#include "voreen/qt/networkeditor/editor_settings.h"

namespace voreen {

class ProcessorGraphicsItem;

/**
 * The tooltip class of all processor items. At the moment it shows only the errorState of the processor.
 * If the processor has no error, no tooltip will be shown.
 */
class ToolTipProcessorGraphicsItem : public ToolTipBaseGraphicsItem {
public:
    /**
     * Constructor
     */
    ToolTipProcessorGraphicsItem(ProcessorGraphicsItem* item);

//---------------------------------------------------------------------------------
//      nwebasegraphicsitem functions
//---------------------------------------------------------------------------------
    ///the type
    int type() const {return UserTypesToolTipProcessorGraphicsItem;}
    //style
    QRectF boundingRect() const;
    QPainterPath shape() const;
protected:
    void initializePaintSettings();
    void mainPaint(QPainter* painter, const QStyleOptionGraphicsItem* option, QWidget* widget, NWEItemSettings& setting);

//---------------------------------------------------------------------------------
//      tooltipbase functions
//---------------------------------------------------------------------------------
public:
    /**
     * Sets the tooltip text.
     * @see updateToolTip of base class.
     */
    bool updateToolTip(qreal x, qreal y);
};

} //namespace voreen

#endif // VRN_TOOLTIPPROCESSORGRAPHICSITEM_H

