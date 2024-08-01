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

#ifndef VRN_TRANSFUNC2DPROPERTYWIDGET_H
#define VRN_TRANSFUNC2DPROPERTYWIDGET_H

#include "voreen/qt/widgets/property/transfunc/transfuncpropertywidgetbase.h"
#include "voreen/qt/widgets/property/transfunc/2d/transfunc2dpropertywidgetpainter.h"

#include "voreen/core/properties/transfunc/transfuncpropertybase.h"

#include "tgt/qt/qtcanvas.h"
namespace voreen {

/**
 * Property widget used to represent 1D transfer functions.
 */
class TransFunc2DPropertyWidget : public TransFuncPropertyWidgetBase {
    Q_OBJECT
public:
    /** Constructor */
    TransFunc2DPropertyWidget(TransFuncPropertyBase* prop, QWidget* parent = 0);
    /** Destructor */
    ~TransFunc2DPropertyWidget();

    //--------------------
    // override functions
    //--------------------
protected:
    virtual WidgetBaseLayout getBaseLayout() const;
    virtual TransFuncPropertyWidgetPainterBase* createPreviewCanvasPainter(tgt::QtCanvas* canvas) const;
    virtual Histogram* getHistogramFromVolume(const VolumeBase* vb, const size_t channel) const;
    virtual void restoreZoomMetaData() const;
protected slots:
    virtual void storeZoomMetaData() const;
    virtual void createDomainMenuSlot();
    virtual void loadDomainSlot(QAction* action);

};


} // namespace voreen

#endif // VRN_TRANSFUNC2DPROPERTYWIDGET_H
