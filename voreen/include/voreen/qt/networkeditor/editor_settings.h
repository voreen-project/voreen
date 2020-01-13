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

#ifndef VRN_NETWORKEDITOR_COMMON_H
#define VRN_NETWORKEDITOR_COMMON_H

#include "tgt/matrix.h"

#include <vector>
#include <string>

#include <QtGlobal>
#include <QList>
#include <QGraphicsItem>
#include <QColor>

namespace voreen {

enum NWEBaseGraphicsItemUserTypes {
    UserTypesProcessorGraphicsItem           = QGraphicsItem::UserType +  1,
    UserTypesPortGraphicsItem                = QGraphicsItem::UserType +  2,
    UserTypesPropertyGraphicsItem            = QGraphicsItem::UserType +  3,
    UserTypesPropertyListGraphicsItem        = QGraphicsItem::UserType +  4,
    UserTypesPortArrowGraphicsItem           = QGraphicsItem::UserType +  5,
    UserTypesWidgetToggleButtonGraphicsItem  = QGraphicsItem::UserType +  6,
    UserTypesPropertyListButtonGraphicsItem  = QGraphicsItem::UserType +  7,
    UserTypesPortOwnerLinkArrowGraphicsItem  = QGraphicsItem::UserType +  8,
    UserTypesPropertyLinkArrowGraphicsItem   = QGraphicsItem::UserType +  9,
    UserTypesPortSizeLinkArrowGraphicsItem   = QGraphicsItem::UserType + 10,
    //UserTypesAggregationGraphicsItem         = QGraphicsItem::UserType + 11, // No longer supported
    UserTypesRenamableTextGraphicsItem       = QGraphicsItem::UserType + 12,
    UserTypesProgressBarGraphicsItem         = QGraphicsItem::UserType + 13,
    UserTypesToolTipProcessorGraphicsItem    = QGraphicsItem::UserType + 20,
    UserTypesToolTipPortGraphicsItem         = QGraphicsItem::UserType + 21,
    UserTypesTextBoxBaseGraphicsItem         = QGraphicsItem::UserType + 30,
    UserTypesTextBoxGraphicsItem             = QGraphicsItem::UserType + 31,
    UserTypesFrameBoxGraphicsItem            = QGraphicsItem::UserType + 32
};

enum NWEBaseGraphicsItemZValues {
    ZValuesTextBoxBaseGraphicsItem = 1,
    ZValuesFrameBoxGraphicsItem = 2,
    ZValuesTextBoxGraphicsItem = 3,

    ZValuesPortArrowGraphicsItem = 5,
    ZValuesSelectedPortArrowGraphicsItem = 6,

    ZValuesPortOwnerGraphicsItem = 10,
    ZValuesPortGraphicsItem = 11,
    ZValuesPropertyListButtonGraphicsItem = 11,

    ZValuesPortOwnerLinkArrowGraphicsItem = 33,
    ZValuesPropertyLinkArrowGraphicsItem = 34,
    ZValuesPropertyListGraphicsItem = 35,
    ZValuesPropertyGraphicsItem = 35,

    ZValuesPortSizeLinkArrowGraphicsItem = 40,

    ZValuesSelectedPortOwnerGraphicsItem = 50,
    ZValuesGrabbedPortArrowGraphicsItem = 51,
    ZValuesSelectedPropertyLinkArrowGraphicsItem = 52,
    ZValuesSelectedPortOwnerLinkArrowGraphicsItem = 53,

    ZValuesToolTipGraphicsItem = 100
};

enum NetworkEditorStyles {
    NWESTYLE_CLASSIC = 0,
    NWESTYLE_CLASSIC_PRINT = 1
};

enum NetworkEditorGraphLayouts {
    NWEGL_SUGIYAMA = 0
};

    /**
     * NetworkEditor Settings
     */
    /// editor settings
    const QSize minimumSize = QSize(200, 200);
    //const qreal tooltiptimerRadius = 3.f;
    const QRectF sceneRectSpacing = QRect(-100000, -100000, 200000, 200000);
    const qreal scaleFactorFactor = 0.9f;
    const float minimalScaleFactor = 1.25f;
    const qreal keyPressScaleFactor = pow(2.f, 0.35f);
    // editor button settings
    const QSize NWEMainButtonSize = QSize(42, 42);
    const QSize NWESub1ButtonSize = QSize(36, 36);
    const QSize NWESub2ButtonSize = QSize(30, 30);
    const int NWEButtonBackgroundMargin = 4;
    const int NWEMarginLayerToLinking = 20;
    const int NWEButtonBorderSpacingX = 10;
    const int NWEButtonBorderSpacingY = 10;
    //ToolTip settings
    const int NWEToolTipTime = 500; //msec

    /**
     * PortOwner Settings
     */
    const qreal buttonsOffsetX = 6.0;
    const qreal buttonsOffsetY = 6.0;
    const qreal drawingRectMinimumWidth = 80.0;
    const qreal drawingRectMinimumHeight = 60.0;
    const qreal minimumDistanceToStartDrawingArrow = 10.0;
    /**
     * Processor Settings
     */
    const qreal utilityOpacityValue = 0.5;


    /**
     * ConnectionBase Settings
     */
    const qreal cbArrowHeadSize = 10.0;

    /**
     * PortArrow Settings
     */
    const qreal portArrowHoverSize = 2.5f;
    const QColor paTransparentColor(80, 80, 80, 100);
    const qreal paEventRadius = 15.f;







static tgt::Matrix3d QTransformToMat3(const QTransform& transform) {
    return tgt::mat3(transform.m11(), transform.m12(), transform.m13(),
                     transform.m21(), transform.m22(), transform.m23(),
                     transform.m31(), transform.m32(), transform.m33());
}

static QTransform Mat3ToQTransform(const tgt::Matrix3d& matrix) {
    return QTransform(matrix.t00, matrix.t01, matrix.t02,
                      matrix.t10, matrix.t11, matrix.t12,
                      matrix.t20, matrix.t21, matrix.t22);
}

template<class T>
QList<T> stdVectorToQList(const std::vector<T>& vec) {
    QList<T> result;

    foreach (T i, vec)
        result.append(i);

    return result;
}

template<class T>
std::vector<T> qListToStdVector(const QList<T>& vec) {
    std::vector<T> result;

    foreach (T i, vec)
        result.push_back(i);

    return result;
}

template<class T, class U>
QList<U> convertQList(const QList<T>& list, bool ignoreFailedConversion = true) {
    QList<U> result;

    foreach (T item, list) {
        U castItem = dynamic_cast<U>(item);

        if (!ignoreFailedConversion) {
            tgtAssert(castItem, "dynamic_cast failed");
        }

        if (castItem)
            result.append(castItem);
    }

    return result;
}

} // namespace

#endif // VRN_NETWORKEDITOR_COMMON_H
