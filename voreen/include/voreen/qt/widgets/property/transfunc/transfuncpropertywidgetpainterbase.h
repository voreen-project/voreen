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

#ifndef VRN_TRANSFUNCPROPERTYWIDGETPAINTERBASE_H
#define VRN_TRANSFUNCPROPERTYWIDGETPAINTERBASE_H

#include "voreen/qt/widgets/property/transfunc/utils/transfuncpropertywidgetpainterslider.h"

#include "voreen/qt/voreenqtapi.h"

#include "tgt/painter.h"
#include "tgt/event/eventlistener.h"
#include "tgt/vector.h"

#include <QObject>
#include <QPoint>
#include <QString>

namespace voreen {

class TransFuncBase;
class Histogram;

/**
 * The base painter used in transfer function property widgets.
 * It displays the tf texture and can be used to adjust domains and gamma value.
 */
class VRN_QT_API TransFuncPropertyWidgetPainterBase : public QObject, public tgt::Painter, public tgt::EventListener {
    Q_OBJECT
public:
    /**
     * @param canvas canvas that belongs to this painter
     */
    TransFuncPropertyWidgetPainterBase(tgt::GLCanvas* canvas);
    ~TransFuncPropertyWidgetPainterBase();

    //--------------------
    //  getter and setter
    //--------------------
    /** Sets the transfer function, which will be displayed. */
    virtual void setTransFunc(TransFuncBase* tf) = 0;
    /** Sets the histogram. */
    virtual void setHistogram(const Histogram* histogram) = 0;

signals:
    void changedGamma();            ///< signal which is emited, if the gamma value has been changed
    void changedDomain();           ///< signal which is emited, if the domain sliders have changed
    void storeZoomMetaDataSignal(); ///< forces the widget to store the current zoom level
    void interaction(bool inter); ///< toggles interaction mode
    void showInfoToolTip(QPoint pos, QString tip); ///< signal to show the info tooltip on slider changes
    void hideInfoToolTip();       ///< signal to hide the info tooltip

    //--------------------
    //  zoom functions
    //--------------------
public:
    /** Zooms in by changing the minimal and maximal domain values. */
    virtual void zoomIn() = 0;
    /** Zooms in by changing the minimal and maximal domain values. */
    virtual void zoomOut() = 0;
    /** Fits min/maxDomainValue to tf_ domain. */
    virtual void resetZoom() = 0;
};

} // namespace voreen

#endif // VRN_TRANSFUNCPROPERTYWIDGETPAINTERBASE_H
