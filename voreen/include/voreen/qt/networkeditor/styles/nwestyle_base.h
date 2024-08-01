/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2024 University of Muenster, Germany,                        *
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

#ifndef VRN_NWESTYLE_BASE_H
#define VRN_NWESTYLE_BASE_H

#include "voreen/qt/networkeditor/editor_settings.h"

//gi
#include "voreen/qt/networkeditor/graphicitems/core/portgraphicsitem.h"
#include "voreen/qt/networkeditor/graphicitems/core/processorgraphicsitem.h"
#include "voreen/qt/networkeditor/graphicitems/core/propertygraphicsitem.h"
#include "voreen/qt/networkeditor/graphicitems/core/propertylistgraphicsitem.h"
#include "voreen/qt/networkeditor/graphicitems/connections/portarrowgraphicsitem.h"
#include "voreen/qt/networkeditor/graphicitems/connections/portsizelinkarrowgraphicsitem.h"
#include "voreen/qt/networkeditor/graphicitems/connections/portownerlinkarrowgraphicsitem.h"
#include "voreen/qt/networkeditor/graphicitems/connections/propertylinkarrowgraphicsitem.h"
#include "voreen/qt/networkeditor/graphicitems/utils/progressbargraphicsitem.h"
#include "voreen/qt/networkeditor/graphicitems/utils/propertylistbuttongraphicsitem.h"
#include "voreen/qt/networkeditor/graphicitems/utils/widgettogglebuttongraphicsitem.h"
#include "voreen/qt/networkeditor/graphicitems/tooltips/tooltipportgraphicsitem.h"
#include "voreen/qt/networkeditor/graphicitems/tooltips/tooltipprocessorgraphicsitem.h"
#include "voreen/qt/networkeditor/graphicitems/textboxes/textboxbasegraphicsitem.h"
#include "voreen/qt/networkeditor/graphicitems/textboxes/textboxgraphicsitem.h"
#include "voreen/qt/networkeditor/graphicitems/textboxes/frameboxgraphicsitem.h"

//core
#include "voreen/core/ports/port.h"
#include "voreen/core/ports/coprocessorport.h"
#include "voreen/core/ports/renderport.h"
#include "voreen/core/processors/processor.h"
#include "voreen/core/properties/cameraproperty.h"
#include "voreen/core/processors/processorwidget.h"

//qt
#include <QPainter>
#include <QStyle>
#include <QStyleOption>
#include <QColor>
#include <QtSvg/QSvgRenderer>

namespace voreen {
    /// Struct used to return all needed informations to draw a graphicsitem in the right state
    struct NWEItemSettings;

/**
 * Abstract base class of NWEStyles.
 * All functions have to be implemented in the concrete used style. An example is NWESytle_Classis.
 * Each graphicsitem has four functions, which are required for painting.
 * @note It is useful to call getActualColorSetting in each function.
 *
 *  \**
 *   * Returns the bounding rect of the given graphics item. This functions is used by QGraphicsScene etc...
 *   * @param item the actual item.
 *   * @return the bounding box of 'item'.
 *   *\
 *  QRectF xxxGI_boundingRect(const xxxGraphicsItem* item) const
 *
 *  \**
 *   * This function returns the exact shape of the graphicsitem. It is for instance used in the item selection.
 *   * The shape has to fit into the bounding box.
 *   * @param item the actual item.
 *   * @return the shape of 'item'.
 *   *\
 *  QPainterPath xxxGI_shape(const xGraphicsItem* item) const
 *
 *  \**
 *   * This function is called once before the first paint call. Settings like the font color and size should be set here.
 *   * @param item the actual item.
 *  void xxxGI_initializePaintSettings(xxxGraphicsItem* item)
 *
 *  \**
 *   * The main paint function used by Qt.
 *   * @param item the actual item.
 *   * @param painter the painter used in this function.
 *   * @param option the style option of 'item'.
 *   * @param widget a special widget, whose default value is null.
 *   * @param setting struct containing informations needed to draw the item appropriate.
 *  void xxxGI_paint(xxxGraphicsItem* item, QPainter* painter, const QStyleOptionGraphicsItem* option, QWidget* widget, NWEItemSettings& setting)
 */
class NWEStyle_Base {
public:
    /*********************************************************************
     *                       General Color Defines
     ********************************************************************/
    //editor
    virtual QColor getButtonBackgroundColor() const = 0;     //< Color of NetworkEditor buttons
    virtual QBrush getBackgroundBrush() const = 0;           //< Color/Brush of the NetworkEditor background

    //general graphicsitems
    virtual QColor getSelectionColor() const = 0;            //< Selection color of each graphicsitem
    virtual QColor getHoverColor() const = 0;                //< Hover color of each graphicsitem

    virtual QColor getConnectionYes() const = 0;             //< Color of potential connections if they are valid
    virtual QColor getConnectionMaybe() const = 0;           //< Color of potential connections if they are partly valid
    virtual QColor getConnectionNo() const = 0;              //< Color of potential connections if they are not valid

    //port
    //processor
    virtual QColor getProcessorColor1() const = 0;           //< Default color of the processor graphicsitem
    //property
    //propertylist
    //progressbar
    //propertylistbutton
    //widgettogglebutton
    //portarrow
    virtual QColor getPortArrowColor() const = 0;            //< Default color of the portarrow graphicsitem
    //portownerlinkarrow
    virtual QColor getPortOwnerLinkArrowColor() const = 0;   //< Default color of the portownrlinkarrow graphicsitem
    //propertylinkarrow
    virtual QColor getPropertyLinkArrowColor() const = 0;    //< Default color of the propertylinkarrow graphicsitem
    //portsizelinkarrow
    virtual QColor getPortSizeLinkArrowColor() const = 0;    //< Default color of the portsizelinkarrow graphicsitem
    //tool tips
    virtual QColor getToolTipBackgroundColor() const = 0;    //< Background color of all tooltips
    //text boxes
    virtual QColor getTextBoxBaseMainColor() const = 0;      //< Main color of text boxes
    //shadows
    virtual bool getShadowsEnabled() const = 0;              //< Enable shadows to items that support it

    static QSvgRenderer NWEStyle_Error1Renderer;     //< Renderer for error type 1 sign
    static QSvgRenderer NWEStyle_Error2Renderer;     //< Renderer for error type 2 sign
    static const QString NWEStyle_Error1SVGPath;     //< Path to the error 1 svg element
    static const QString NWEStyle_Error2SVGPath;     //< Path to the error 2 svg element

    /**
     * Enum used in 'ItemSetting'. At the moment only supported in processor graphicsitem.
     */
    enum ErrorState {
        ES_NO_ERROR = 0,
        ES_ERROR_T1 = 1,
        ES_ERROR_T2 = 2
    };
protected:
    NetworkEditor* nwe_;    ///< pointer to the editor. Used for font size.
public:
    NWEStyle_Base(NetworkEditor* networkeditor){
        NWEStyle_Error1Renderer.load(NWEStyle_Error1SVGPath);
        NWEStyle_Error2Renderer.load(NWEStyle_Error2SVGPath);
        nwe_ = networkeditor;
    };

    virtual ~NWEStyle_Base(){}

    /*********************************************************************
     *                       General Color Function
     ********************************************************************/

    /**
     * Retruns the color and depth settings of the given graphicsitem. Should be called at the beginning of the paint function.
     * @param item the nwebasegraphicsitem of which the color and depth settings are requested.
     * @param option QStyleOption used for settings calculation.
     * @return the color and depth settings of the 'item'.
     */
    NWEItemSettings getItemSettings(const NWEBaseGraphicsItem* item, const QStyleOptionGraphicsItem* option);

    /*********************************************************************
     *                       Core Elements
     ********************************************************************/
    //port
    virtual QRectF PortGI_boundingRect(const PortGraphicsItem* item) const = 0;
    virtual QPainterPath PortGI_shape(const PortGraphicsItem* item) const = 0;
    virtual void PortGI_initializePaintSettings(PortGraphicsItem* item) = 0;
    virtual void PortGI_paint(PortGraphicsItem* item, QPainter* painter, const QStyleOptionGraphicsItem* option, QWidget* widget, NWEItemSettings& setting) = 0;

    //processor
    virtual QRectF ProcessorGI_boundingRect(const ProcessorGraphicsItem* item) const = 0;
    virtual QPainterPath ProcessorGI_shape(const ProcessorGraphicsItem* item) const = 0;
    virtual void ProcessorGI_initializePaintSettings(ProcessorGraphicsItem* item) = 0;
    virtual void ProcessorGI_paint(ProcessorGraphicsItem* item, QPainter* painter, const QStyleOptionGraphicsItem* option, QWidget* widget, NWEItemSettings& setting) = 0;

    //property
    virtual QRectF PropertyGI_boundingRect(const PropertyGraphicsItem* item) const = 0;
    virtual QPainterPath PropertyGI_shape(const PropertyGraphicsItem* item) const = 0;
    virtual void PropertyGI_initializePaintSettings(PropertyGraphicsItem* item) = 0;
    virtual void PropertyGI_paint(PropertyGraphicsItem* item, QPainter* painter, const QStyleOptionGraphicsItem* option, QWidget* widget, NWEItemSettings& setting) = 0;

    //property list
    virtual QRectF PropertyListGI_boundingRect(const PropertyListGraphicsItem* item) const = 0;
    virtual QPainterPath PropertyListGI_shape(const PropertyListGraphicsItem* item) const = 0;
    virtual void PropertyListGI_initializePaintSettings(PropertyListGraphicsItem* item) = 0;
    virtual void PropertyListGI_paint(PropertyListGraphicsItem* item, QPainter* painter, const QStyleOptionGraphicsItem* option, QWidget* widget, NWEItemSettings& setting) = 0;

    /*********************************************************************
     *                       Util Elements
     ********************************************************************/
    //ProgressBarGraphicsItem
    virtual QRectF ProgressBarGI_boundingRect(const ProgressBarGraphicsItem* item) const = 0;
    virtual QPainterPath ProgressBarGI_shape(const ProgressBarGraphicsItem* item) const = 0;
    virtual void ProgressBarGI_initializePaintSettings(ProgressBarGraphicsItem* item) = 0;
    virtual void ProgressBarGI_paint(ProgressBarGraphicsItem* item, QPainter* painter, const QStyleOptionGraphicsItem* option, QWidget* widget, NWEItemSettings& setting) = 0;

    //PropertyListButtonGraphicsItem
    virtual QRectF PropertyListButtonGI_boundingRect(const PropertyListButtonGraphicsItem* item) const = 0;
    virtual QPainterPath PropertyListButtonGI_shape(const PropertyListButtonGraphicsItem* item) const = 0;
    virtual void PropertyListButtonGI_initializePaintSettings(PropertyListButtonGraphicsItem* item) = 0;
    virtual void PropertyListButtonGI_paint(PropertyListButtonGraphicsItem* item, QPainter* painter, const QStyleOptionGraphicsItem* option, QWidget* widget, NWEItemSettings& setting) = 0;

    //WidgetToggleButtonGraphicsItem
    virtual QRectF WidgetToggleButtonGI_boundingRect(const WidgetToggleButtonGraphicsItem* item) const = 0;
    virtual QPainterPath WidgetToggleButtonGI_shape(const WidgetToggleButtonGraphicsItem* item) const = 0;
    virtual void WidgetToggleButtonGI_initializePaintSettings(WidgetToggleButtonGraphicsItem* item) = 0;
    virtual void WidgetToggleButtonGI_paint(WidgetToggleButtonGraphicsItem* item, QPainter* painter, const QStyleOptionGraphicsItem* option, QWidget* widget, NWEItemSettings& setting) = 0;

    /*********************************************************************
     *                       Connection Elements
     ********************************************************************/
    //port arrow
    virtual QRectF PortArrowGI_boundingRect(const PortArrowGraphicsItem* item) const = 0;
    virtual QPainterPath PortArrowGI_shape(const PortArrowGraphicsItem* item) const = 0;
    virtual void PortArrowGI_initializePaintSettings(const PortArrowGraphicsItem* item) = 0;
    virtual void PortArrowGI_paint(PortArrowGraphicsItem* item, QPainter* painter, const QStyleOptionGraphicsItem* option, QWidget* widget, NWEItemSettings& setting) = 0;

    //port size link arrow
    virtual QRectF PortSizeLinkArrowGI_boundingRect(const PortSizeLinkArrowGraphicsItem* item) const = 0;
    virtual QPainterPath PortSizeLinkArrowGI_shape(const PortSizeLinkArrowGraphicsItem* item) const = 0;
    virtual void PortSizeLinkArrowGI_initializePaintSettings(PortSizeLinkArrowGraphicsItem* item) = 0;
    virtual void PortSizeLinkArrowGI_paint(PortSizeLinkArrowGraphicsItem* item, QPainter* painter, const QStyleOptionGraphicsItem* option, QWidget* widget, NWEItemSettings& setting) = 0;

    //processor link arrow
    virtual QRectF PortOwnerLinkArrowGI_boundingRect(const PortOwnerLinkArrowGraphicsItem* item) const = 0;
    virtual QPainterPath PortOwnerLinkArrowGI_shape(const PortOwnerLinkArrowGraphicsItem* item) const = 0;
    virtual void PortOwnerLinkArrowGI_initializePaintSettings(PortOwnerLinkArrowGraphicsItem* item) = 0;
    virtual void PortOwnerLinkArrowGI_paint(PortOwnerLinkArrowGraphicsItem* item, QPainter* painter, const QStyleOptionGraphicsItem* option, QWidget* widget, NWEItemSettings& setting) = 0;

    //property link arrow
    virtual QRectF PropertyLinkArrowGI_boundingRect(const PropertyLinkArrowGraphicsItem* item) const = 0;
    virtual QPainterPath PropertyLinkArrowGI_shape(const PropertyLinkArrowGraphicsItem* item) const = 0;
    virtual void PropertyLinkArrowGI_initializePaintSettings(const PropertyLinkArrowGraphicsItem* item) = 0;
    virtual void PropertyLinkArrowGI_paint(PropertyLinkArrowGraphicsItem* item, QPainter* painter, const QStyleOptionGraphicsItem* option, QWidget* widget, NWEItemSettings& setting) = 0;

    /*********************************************************************
     *                       ToolTip Elements
     ********************************************************************/
    //port tooltip
    virtual QRectF ToolTipPortGI_boundingRect(const ToolTipPortGraphicsItem* item) const = 0;
    virtual QPainterPath ToolTipPortGI_shape(const ToolTipPortGraphicsItem* item) const = 0;
    virtual void ToolTipPortGI_initializePaintSettings(ToolTipPortGraphicsItem* item) = 0;
    virtual void ToolTipPortGI_paint(ToolTipPortGraphicsItem* item, QPainter* painter, const QStyleOptionGraphicsItem* option, QWidget* widget, NWEItemSettings& setting) = 0;
    //port tooltip
    virtual QRectF ToolTipProcessorGI_boundingRect(const ToolTipProcessorGraphicsItem* item) const = 0;
    virtual QPainterPath ToolTipProcessorGI_shape(const ToolTipProcessorGraphicsItem* item) const = 0;
    virtual void ToolTipProcessorGI_initializePaintSettings(ToolTipProcessorGraphicsItem* item) = 0;
    virtual void ToolTipProcessorGI_paint(ToolTipProcessorGraphicsItem* item, QPainter* painter, const QStyleOptionGraphicsItem* option, QWidget* widget, NWEItemSettings& setting) = 0;
    /*********************************************************************
     *                       TextBoxes Elements
     ********************************************************************/
    //text boxes
    virtual QRectF TextBoxBaseGI_boundingRect(const TextBoxBaseGraphicsItem* item) const = 0;
    virtual QPainterPath TextBoxBaseGI_shape(const TextBoxBaseGraphicsItem* item) const = 0;
    virtual void TextBoxGI_initializePaintSettings(TextBoxGraphicsItem* item) = 0;
    virtual void TextBoxGI_paint(TextBoxGraphicsItem* item, QPainter* painter, const QStyleOptionGraphicsItem* option, QWidget* widget, NWEItemSettings& setting) = 0;
    virtual void FrameBoxGI_initializePaintSettings(FrameBoxGraphicsItem* item) = 0;
    virtual void FrameBoxGI_paint(FrameBoxGraphicsItem* item, QPainter* painter, const QStyleOptionGraphicsItem* option, QWidget* widget, NWEItemSettings& setting) = 0;

};

/**
 * Struct containing all information needed to paint the GraphicsItem.
 * Will be automaticly calculated assigned to a parameter in the 'mainPaint' function of the item.
 */
struct NWEItemSettings{
    QColor color1_; //< primary color
    qreal opacity_;
    qreal zValue_;
    NWEStyle_Base::ErrorState error_; //only used in processor graphicsitem
    NWEItemSettings(QColor color1 = Qt::black, qreal opacity = 1.0, qreal zValue = 0.0, NWEStyle_Base::ErrorState error = NWEStyle_Base::ES_NO_ERROR)
        : color1_(color1)
        , opacity_(opacity)
        , zValue_(zValue)
        , error_(error)
    {}
};

} //namespace voreen

#endif // VRN_NWESTYLE_BASE_H

