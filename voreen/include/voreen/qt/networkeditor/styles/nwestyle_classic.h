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

#ifndef VRN_NWESTYLE_CLASSIC_H
#define VRN_NWESTYLE_CLASSIC_H

#include "nwestyle_base.h"

namespace voreen {

/**
 * The classic look of Voreen.
 */
class NWEStyle_Classic : public NWEStyle_Base{
public:

    /*********************************************************************
     *                       General Color Defines
     ********************************************************************/
    virtual QColor getButtonBackgroundColor() const;
    virtual QBrush getBackgroundBrush() const;
    virtual QColor getSelectionColor() const;
    virtual QColor getHoverColor() const;
    virtual QColor getConnectionYes() const;
    virtual QColor getConnectionMaybe() const;
    virtual QColor getConnectionNo() const;
    virtual QColor getProcessorColor1() const;
    virtual QColor getPortArrowColor() const;
    virtual QColor getPortOwnerLinkArrowColor() const;
    virtual QColor getPropertyLinkArrowColor() const;
    virtual QColor getPortSizeLinkArrowColor() const;
    virtual QColor getToolTipBackgroundColor() const;
    virtual QColor getTextBoxBaseMainColor() const;
    virtual bool getShadowsEnabled() const;


    NWEStyle_Classic(NetworkEditor* networkeditor);
    virtual ~NWEStyle_Classic();

    /*********************************************************************
     *                       Core Elements
     ********************************************************************/

    //port
    virtual QRectF PortGI_boundingRect(const PortGraphicsItem* item) const;
    virtual QPainterPath PortGI_shape(const PortGraphicsItem* item) const;
    virtual void PortGI_initializePaintSettings(PortGraphicsItem* item);
    virtual void PortGI_paint(PortGraphicsItem* item, QPainter* painter, const QStyleOptionGraphicsItem* option, QWidget* widget, NWEItemSettings& setting);

    //processor
    virtual QRectF ProcessorGI_boundingRect(const ProcessorGraphicsItem* item) const;
    virtual QPainterPath ProcessorGI_shape(const ProcessorGraphicsItem* item) const;
    virtual void ProcessorGI_initializePaintSettings(ProcessorGraphicsItem* item);
    virtual void ProcessorGI_paint(ProcessorGraphicsItem* item, QPainter* painter, const QStyleOptionGraphicsItem* option, QWidget* widget, NWEItemSettings& setting);

    //property
    virtual QRectF PropertyGI_boundingRect(const PropertyGraphicsItem* item) const;
    virtual QPainterPath PropertyGI_shape(const PropertyGraphicsItem* item) const;
    virtual void PropertyGI_initializePaintSettings(PropertyGraphicsItem* item);
    virtual void PropertyGI_paint(PropertyGraphicsItem* item, QPainter* painter, const QStyleOptionGraphicsItem* option, QWidget* widget, NWEItemSettings& setting);

    //property list
    virtual QRectF PropertyListGI_boundingRect(const PropertyListGraphicsItem* item) const;
    virtual QPainterPath PropertyListGI_shape(const PropertyListGraphicsItem* item) const;
    virtual void PropertyListGI_initializePaintSettings(PropertyListGraphicsItem* item);
    virtual void PropertyListGI_paint(PropertyListGraphicsItem* item, QPainter* painter, const QStyleOptionGraphicsItem* option, QWidget* widget, NWEItemSettings& setting);

    /*********************************************************************
     *                       Util Elements
     ********************************************************************/
    //ProgressBarGraphicsItem
    virtual QRectF ProgressBarGI_boundingRect(const ProgressBarGraphicsItem* item) const;
    virtual QPainterPath ProgressBarGI_shape(const ProgressBarGraphicsItem* item) const;
    virtual void ProgressBarGI_initializePaintSettings(ProgressBarGraphicsItem* item);
    virtual void ProgressBarGI_paint(ProgressBarGraphicsItem* item, QPainter* painter, const QStyleOptionGraphicsItem* option, QWidget* widget, NWEItemSettings& setting);

    //PropertyListButtonGraphicsItem
    virtual QRectF PropertyListButtonGI_boundingRect(const PropertyListButtonGraphicsItem* item) const;
    virtual QPainterPath PropertyListButtonGI_shape(const PropertyListButtonGraphicsItem* item) const;
    virtual void PropertyListButtonGI_initializePaintSettings(PropertyListButtonGraphicsItem* item);
    virtual void PropertyListButtonGI_paint(PropertyListButtonGraphicsItem* item, QPainter* painter, const QStyleOptionGraphicsItem* option, QWidget* widget, NWEItemSettings& setting);

    //WidgetToggleButtonGraphicsItem
    virtual QRectF WidgetToggleButtonGI_boundingRect(const WidgetToggleButtonGraphicsItem* item) const;
    virtual QPainterPath WidgetToggleButtonGI_shape(const WidgetToggleButtonGraphicsItem* item) const;
    virtual void WidgetToggleButtonGI_initializePaintSettings(WidgetToggleButtonGraphicsItem* item);
    virtual void WidgetToggleButtonGI_paint(WidgetToggleButtonGraphicsItem* item, QPainter* painter, const QStyleOptionGraphicsItem* option, QWidget* widget, NWEItemSettings& setting);

    /*********************************************************************
     *                       Connection Elements
     ********************************************************************/
    //port arrow
    virtual QRectF PortArrowGI_boundingRect(const PortArrowGraphicsItem* item) const;
    virtual QPainterPath PortArrowGI_shape(const PortArrowGraphicsItem* item) const;
    virtual void PortArrowGI_initializePaintSettings(const PortArrowGraphicsItem* item);
    virtual void PortArrowGI_paint(PortArrowGraphicsItem* item, QPainter* painter, const QStyleOptionGraphicsItem* option, QWidget* widget, NWEItemSettings& setting);

    //port size link arrow
    virtual QRectF PortSizeLinkArrowGI_boundingRect(const PortSizeLinkArrowGraphicsItem* item) const;
    virtual QPainterPath PortSizeLinkArrowGI_shape(const PortSizeLinkArrowGraphicsItem* item) const;
    virtual void PortSizeLinkArrowGI_initializePaintSettings(PortSizeLinkArrowGraphicsItem* item);
    virtual void PortSizeLinkArrowGI_paint(PortSizeLinkArrowGraphicsItem* item, QPainter* painter, const QStyleOptionGraphicsItem* option, QWidget* widget, NWEItemSettings& setting);


    //processor link arrow
    virtual QRectF PortOwnerLinkArrowGI_boundingRect(const PortOwnerLinkArrowGraphicsItem* item) const;
    virtual QPainterPath PortOwnerLinkArrowGI_shape(const PortOwnerLinkArrowGraphicsItem* item) const;
    virtual void PortOwnerLinkArrowGI_initializePaintSettings(PortOwnerLinkArrowGraphicsItem* item);
    virtual void PortOwnerLinkArrowGI_paint(PortOwnerLinkArrowGraphicsItem* item, QPainter* painter, const QStyleOptionGraphicsItem* option, QWidget* widget, NWEItemSettings& setting);

    //property link arrow
    virtual QRectF PropertyLinkArrowGI_boundingRect(const PropertyLinkArrowGraphicsItem* item) const;
    virtual QPainterPath PropertyLinkArrowGI_shape(const PropertyLinkArrowGraphicsItem* item) const;
    virtual void PropertyLinkArrowGI_initializePaintSettings(const PropertyLinkArrowGraphicsItem* item);
    virtual void PropertyLinkArrowGI_paint(PropertyLinkArrowGraphicsItem* item, QPainter* painter, const QStyleOptionGraphicsItem* option, QWidget* widget, NWEItemSettings& setting);

    /*********************************************************************
     *                       ToolTip Elements
     ********************************************************************/
    //port tooltip
    virtual QRectF ToolTipPortGI_boundingRect(const ToolTipPortGraphicsItem* item) const;
    virtual QPainterPath ToolTipPortGI_shape(const ToolTipPortGraphicsItem* item) const;
    virtual void ToolTipPortGI_initializePaintSettings(ToolTipPortGraphicsItem* item);
    virtual void ToolTipPortGI_paint(ToolTipPortGraphicsItem* item, QPainter* painter, const QStyleOptionGraphicsItem* option, QWidget* widget, NWEItemSettings& setting);
    //processor tooltip
    virtual QRectF ToolTipProcessorGI_boundingRect(const ToolTipProcessorGraphicsItem* item) const;
    virtual QPainterPath ToolTipProcessorGI_shape(const ToolTipProcessorGraphicsItem* item) const;
    virtual void ToolTipProcessorGI_initializePaintSettings(ToolTipProcessorGraphicsItem* item);
    virtual void ToolTipProcessorGI_paint(ToolTipProcessorGraphicsItem* item, QPainter* painter, const QStyleOptionGraphicsItem* option, QWidget* widget, NWEItemSettings& setting);
    /*********************************************************************
     *                       TextBoxes Elements
     ********************************************************************/
    //text box base
    virtual QRectF TextBoxBaseGI_boundingRect(const TextBoxBaseGraphicsItem* item) const;
    virtual QPainterPath TextBoxBaseGI_shape(const TextBoxBaseGraphicsItem* item) const;
    virtual void TextBoxGI_initializePaintSettings(TextBoxGraphicsItem* item);
    virtual void TextBoxGI_paint(TextBoxGraphicsItem* item, QPainter* painter, const QStyleOptionGraphicsItem* option, QWidget* widget, NWEItemSettings& setting);
    virtual void FrameBoxGI_initializePaintSettings(FrameBoxGraphicsItem* item);
    virtual void FrameBoxGI_paint(FrameBoxGraphicsItem* item, QPainter* painter, const QStyleOptionGraphicsItem* option, QWidget* widget, NWEItemSettings& setting);
};

} //namespace voreen

#endif // VRN_NWESTYLE_CLASSIC_H

