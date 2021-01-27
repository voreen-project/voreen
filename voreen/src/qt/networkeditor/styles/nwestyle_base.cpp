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

#include "voreen/qt/networkeditor/styles/nwestyle_base.h"

namespace voreen {

NWEItemSettings NWEStyle_Base::getItemSettings(const NWEBaseGraphicsItem* item, const QStyleOptionGraphicsItem* option) {
    QColor color1 = Qt::black;
    qreal opacity = 1.0;
    qreal zValue = 0.0;
    ErrorState error = ES_NO_ERROR;
    //switch over all graphicsitem types
    switch (item->type()) {
    /********************************************
     *      Port                                *
     ********************************************/
    case UserTypesPortGraphicsItem: {
        const PortGraphicsItem* portItem = dynamic_cast<const PortGraphicsItem*>(item);
        color1 = portItem->getColor();
        if (!portItem->getPort()->areConditionsMet())
            error = ES_ERROR_T1;
        switch (item->currentLayer()) {
            case NetworkEditorLayerDataFlow:
                opacity = 1.0;
                break;
            case NetworkEditorLayerGeneralLinking:
            case NetworkEditorLayerCameraLinking:
                opacity = 0.2;
                break;
            case NetworkEditorLayerPortSizeLinking:
                if (!portItem->getPort()->getPropertiesByType<RenderSizeOriginProperty>().empty() ||
                    !portItem->getPort()->getPropertiesByType<RenderSizeReceiveProperty>().empty())
                    opacity = 1.0;
                else
                    opacity = 0.2;
                break;
            default:
                tgtAssert(false, "should not get here");
                break;
        }

        return NWEItemSettings(color1, opacity, zValue, error);
        break;
    }
        /********************************************
         *      Processor                           *
         ********************************************/
    case UserTypesProcessorGraphicsItem: {
        const ProcessorGraphicsItem* procItem = dynamic_cast<const ProcessorGraphicsItem*>(item);
        //get aggregation color and depth
        color1 = getProcessorColor1();
        // frame indicates selected process
        if (option->state & QStyle::State_Selected) {
            color1 = getSelectionColor();
            zValue = ZValuesSelectedPortOwnerGraphicsItem;
        }
            // hover effect
        else if (option->state & QStyle::State_MouseOver) {
            color1 = getHoverColor();
            zValue = ZValuesSelectedPortOwnerGraphicsItem;
        }
            // no effect
        else {
            zValue = ZValuesPortOwnerGraphicsItem;
        }
        //set error
        switch (procItem->getProcessor()->getProcessorState()) {
            case Processor::PROCESSOR_STATE_NOT_INITIALIZED:
                error = ES_ERROR_T2;
                break;
            case Processor::PROCESSOR_STATE_NOT_READY:
                error = ES_ERROR_T1;
                break;
            case Processor::PROCESSOR_STATE_READY:
                error = ES_NO_ERROR;
                break;
            default:
                tgtAssert(false, "Unknown processor state!!!");
                break;
        }
        //get opacity by layer
        bool hasSpecialProp = false;
        switch (item->currentLayer()) {
            case NetworkEditorLayerCameraLinking:
                        foreach (Processor* processor, procItem->getProcessors()) {
                        if (!processor->getPropertiesByType<CameraProperty>().empty()) {
                            hasSpecialProp = true;
                            break;
                        }
                    }
                if (hasSpecialProp)
                    opacity = 1.0;
                else
                    opacity = 0.20;
                break;
            case NetworkEditorLayerPortSizeLinking:
                        foreach (Port* port, procItem->getPorts()) {
                        if (!port->getPropertiesByType<RenderSizeOriginProperty>().empty() ||
                            !port->getPropertiesByType<RenderSizeReceiveProperty>().empty()) {
                            hasSpecialProp = true;
                            break;
                        }
                    }
                if (hasSpecialProp)
                    opacity = 1.0;
                else
                    opacity = 0.20;
                break;
            default:
                opacity = 1.0;
                break;
        }
        return NWEItemSettings(color1, opacity, zValue, error);
        break;
    }
        /********************************************
         *      Property                            *
         ********************************************/
    case UserTypesPropertyGraphicsItem: {
        //nothing yet
        return NWEItemSettings();
        break;
    }
        /********************************************
         *      PropertyList                        *
         ********************************************/
    case UserTypesPropertyListGraphicsItem: {
        //nothing yet
        return NWEItemSettings();
        break;
    }
        /********************************************
         *      ProgressBar                         *
         ********************************************/
    case UserTypesProgressBarGraphicsItem: {
        //nothing yet
        return NWEItemSettings();
        break;
    }
        /********************************************
         *      PropertyListButton                  *
         ********************************************/
    case UserTypesPropertyListButtonGraphicsItem: {
        //nothing yet
        return NWEItemSettings();
        break;
    }
        /********************************************
         *      WidgetToggleButton                  *
         ********************************************/
    case UserTypesWidgetToggleButtonGraphicsItem: {
        //nothing yet
        return NWEItemSettings();
        break;
    }
        /********************************************
         *      PortArrow                           *
         ********************************************/
    case UserTypesPortArrowGraphicsItem: {
        const PortArrowGraphicsItem* arrowItem = dynamic_cast<const PortArrowGraphicsItem*>(item);
        switch (arrowItem->getColorConnectableMode()) {
            case ConnectionBaseGraphicsItem::CCM_DEFAULT:
                if (item->isSelected()) {
                    color1 = getSelectionColor();
                    zValue = ZValuesSelectedPortArrowGraphicsItem;
                } else {
                    color1 = getPortArrowColor();
                    zValue = ZValuesPortArrowGraphicsItem;
                }

                // In case the data condition of the dest port has failed, we mark this connection as invalid.
                if (PortGraphicsItem* portItem = dynamic_cast<PortGraphicsItem*>(arrowItem->getDestinationItem()))
                    if (!portItem->getPort()->areConditionsMet())
                        color1 = getConnectionNo();

                if ((option->state & QStyle::State_MouseOver) || arrowItem->getIsHovered()) {
                    if (color1 == Qt::black) {    // Qt is unable to brighten up Qt::black
                        color1 = getHoverColor();
                    } else {
                        color1 = color1.lighter(130);
                    }
                }
                break;
            case ConnectionBaseGraphicsItem::CCM_NO:
                color1 = getConnectionNo();
                break;
            case ConnectionBaseGraphicsItem::CCM_MAYBE:
                color1 = getConnectionMaybe();
                break;
            case ConnectionBaseGraphicsItem::CCM_YES:
                color1 = getConnectionYes();
                break;
        }
        //push to front
        if (!arrowItem->getDestinationItem())
            zValue = ZValuesGrabbedPortArrowGraphicsItem;
        //get opacity
        if (item->currentLayer() == NetworkEditorLayerDataFlow)
            opacity = 1.0;
        else
            opacity = 0.2;
        //get
        return NWEItemSettings(color1, opacity, zValue, error);
        break;
    }
        /********************************************
         *      PortOwnerLinkArrow                  *
         ********************************************/
    case UserTypesPortOwnerLinkArrowGraphicsItem: {
        const PortOwnerLinkArrowGraphicsItem* polaItem = dynamic_cast<const PortOwnerLinkArrowGraphicsItem*>(item);
        switch (polaItem->getColorConnectableMode()) {
            case ConnectionBaseGraphicsItem::CCM_DEFAULT:
                if (item->isSelected())
                    color1 = getSelectionColor();
                else
                    color1 = getPortOwnerLinkArrowColor();

                if (option->state & QStyle::State_MouseOver) {
                    if (color1 == Qt::lightGray)
                        color1 = getHoverColor();
                    else
                        color1 = color1.lighter();
                }
                break;
            case ConnectionBaseGraphicsItem::CCM_NO:
                color1 = getConnectionNo();
                break;
            case ConnectionBaseGraphicsItem::CCM_MAYBE:
                color1 = getConnectionMaybe();
                break;
            case ConnectionBaseGraphicsItem::CCM_YES:
                color1 = getConnectionYes();
                break;
        }
        return NWEItemSettings(color1, opacity, zValue, error);
        break;
    }
        /********************************************
         *      PropertyLinkArrow                   *
         ********************************************/
    case UserTypesPropertyLinkArrowGraphicsItem: {
        const PropertyLinkArrowGraphicsItem* plaItem = dynamic_cast<const PropertyLinkArrowGraphicsItem*>(item);
        switch (plaItem->getColorConnectableMode()) {
            case ConnectionBaseGraphicsItem::CCM_DEFAULT:
                if (item->isSelected()) {
                    color1 = getSelectionColor();
                    zValue = ZValuesSelectedPropertyLinkArrowGraphicsItem;
                } else {
                    color1 = getPropertyLinkArrowColor();
                    zValue = ZValuesPropertyLinkArrowGraphicsItem;
                }

                if (option->state & QStyle::State_MouseOver || plaItem->getIsHovered()) {
                    if (color1 == Qt::lightGray)
                        color1 = getHoverColor();
                    else
                        color1 = color1.lighter();
                }
                break;
            case ConnectionBaseGraphicsItem::CCM_NO:
                color1 = getConnectionNo();
                break;
            case ConnectionBaseGraphicsItem::CCM_MAYBE:
                color1 = getConnectionMaybe();
                break;
            case ConnectionBaseGraphicsItem::CCM_YES:
                color1 = getConnectionYes();
                break;
        }
        //push to front
        if (option->state & QStyle::State_MouseOver || plaItem->getIsHovered()) {
            zValue = ZValuesSelectedPropertyLinkArrowGraphicsItem;
        }

        return NWEItemSettings(color1, opacity, zValue, error);
        break;
    }
        /********************************************
         *      PortSizeLinkArrow                   *
         ********************************************/
    case UserTypesPortSizeLinkArrowGraphicsItem: {
        const PortSizeLinkArrowGraphicsItem* pslaItem = dynamic_cast<const PortSizeLinkArrowGraphicsItem*>(item);
        switch (pslaItem->getColorConnectableMode()) {
            case ConnectionBaseGraphicsItem::CCM_DEFAULT:
                if (item->isSelected())
                    color1 = getSelectionColor();
                else
                    color1 = getPortSizeLinkArrowColor();

                if ((option->state & QStyle::State_MouseOver) || pslaItem->getIsHovered()) {
                    if (color1 == Qt::black)
                        color1 = getHoverColor();
                    else
                        color1 = color1.lighter();
                }
                break;
            case ConnectionBaseGraphicsItem::CCM_NO:
                color1 = getConnectionNo();
                break;
            case ConnectionBaseGraphicsItem::CCM_MAYBE:
                color1 = getConnectionMaybe();
                break;
            case ConnectionBaseGraphicsItem::CCM_YES:
                color1 = getConnectionYes();
                break;
        }
        return NWEItemSettings(color1, opacity, zValue, error);
        break;
    }
        /********************************************
         *      ToolTipPort                         *
         ********************************************/
    case UserTypesToolTipPortGraphicsItem: {
        //nothing yet
        return NWEItemSettings(getToolTipBackgroundColor(), opacity, zValue, error);
        break;
    }
        /********************************************
         *      ToolTipProcessor                    *
         ********************************************/
    case UserTypesToolTipProcessorGraphicsItem: {
        //nothing yet
        return NWEItemSettings(getToolTipBackgroundColor(), opacity, zValue, error);
        break;
    }
        /********************************************
         *      Text Boxes                          *
         ********************************************/
        //case UserTypesTextBoxBaseGraphicsItem:
    case UserTypesTextBoxGraphicsItem: {
        //const TextBoxGraphicsItem* tbItem = dynamic_cast<const TextBoxGraphicsItem*>(item);
        //get aggregation color and depth
        color1 = getTextBoxBaseMainColor();
        zValue = ZValuesTextBoxGraphicsItem;
        // frame indicates selected process
        if (option->state & QStyle::State_Selected) {
            color1 = getSelectionColor();
            zValue = ZValuesSelectedPortOwnerGraphicsItem;
        }//ignore hover
        return NWEItemSettings(color1, opacity, zValue, error);
        break;
    }
    case UserTypesFrameBoxGraphicsItem: {
        const FrameBoxGraphicsItem* frameItem = dynamic_cast<const FrameBoxGraphicsItem*>(item);
        //get aggregation color and depth
        color1 = frameItem->getBaseColor();
        zValue = ZValuesFrameBoxGraphicsItem;
        // frame indicates selected process
        if (option->state & QStyle::State_Selected) {
            color1 = getSelectionColor();
        }
        //ignore hover
        return NWEItemSettings(color1, opacity, zValue, error);
        break;
    }
        /********************************************
         *      Default                             *
         ********************************************/
    default:
        tgtAssert(false, "NWEBaseGraphicsItem::Type is not known!!!");
        LERRORC("NWEStyle_Base::getActualColorAndDepthSetting",
                "NWEBaseGraphicsItem::Type " << item->type() << "is not known!!!");
        break;
    }
    //should not get here
    return NWEItemSettings();
}

}