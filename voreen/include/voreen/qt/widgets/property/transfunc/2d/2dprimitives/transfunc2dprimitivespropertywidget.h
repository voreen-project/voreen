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

#ifndef VRN_TRANSFUNC2DPRIMITIVESPROPERTYWIDGET_H
#define VRN_TRANSFUNC2DPRIMITIVESPROPERTYWIDGET_H

#include "voreen/qt/widgets/property/transfunc/2d/transfunc2dpropertywidget.h"

#include "voreen/core/properties/transfunc/2d/2dprimitives/transfunc2dprimitivesproperty.h"

namespace voreen {

/**
 * Property widget used to represent 2DPrimitives transfer functions.
 */
class TransFunc2DPrimitivesPropertyWidget : public TransFunc2DPropertyWidget {
    Q_OBJECT
public:
    /** Constructor */
    TransFunc2DPrimitivesPropertyWidget(TransFunc2DPrimitivesProperty* prop, QWidget* parent = 0);
    /** Destructor */
    ~TransFunc2DPrimitivesPropertyWidget();

    //--------------------
    // override functions
    //--------------------
protected:
    virtual QWidget* createToolWindowWidget();

    virtual std::string getPathSubFolder() const {return "2dprimitives";}
    virtual const QPixmap createPresetPreviewIcon(std::string& presetPath) const;

protected slots:
    virtual void loadPresetSlot(QAction* action);



};

} // namespace voreen

#endif // VRN_TRANSFUNC2DPRIMITIVESEDITOR_H
