/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2021 University of Muenster, Germany,                        *
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

#ifndef VRN_TOOLTIPPORTGRAPHICSITEM_H
#define VRN_TOOLTIPPORTGRAPHICSITEM_H

#include "tooltipbasegraphicsitem.h"

class QImage;

namespace voreen {

class PortGraphicsItem;

/**
 * The tooltip class of all port items. It supports an image used for volume or render ports.
 * The shown text is defined by Port::getContentDescriptionHTML().
 */
class ToolTipPortGraphicsItem : public ToolTipBaseGraphicsItem {
public:
    /**
     * Constructor
     */
    ToolTipPortGraphicsItem(PortGraphicsItem* item);

protected:
//-- member --
    QImage* image_;         //< the image used by specific port classes
    float imageAspect_;    //< image aspect

//---------------------------------------------------------------------------------
//      nwebasegraphicsitem functions
//---------------------------------------------------------------------------------
    ///the type
    int type() const {return UserTypesToolTipPortGraphicsItem;}
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
     * Sets the tooltip text and the image.
     * @see updateToolTip of base class.
     */
    bool updateToolTip(qreal x, qreal y);
//-- getter and setter --
    /// returns the image
    const QImage* getImage() const;
    /// returns the aspect
    const float getImageAspect() const;
};

} //namespace voreen

#endif // VRN_TOOLTIPPORTGRAPHICSITEM_H

