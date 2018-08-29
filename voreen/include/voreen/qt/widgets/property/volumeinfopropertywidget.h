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

#ifndef VRN_VOLUMEINFOPROPERTYWIDGET_H
#define VRN_VOLUMEINFOPROPERTYWIDGET_H

#include "voreen/qt/widgets/property/qpropertywidget.h"
#include "voreen/core/properties/volumeinfoproperty.h"

#include "voreen/qt/widgets/customlabel.h"

class QLabel;

namespace voreen {

class VolumeInfoPropertyWidget : public QPropertyWidget {
    Q_OBJECT

public:
    VolumeInfoPropertyWidget(VolumeInfoProperty* volumeInfoProp, QWidget* parent);

    /// Returns the null pointer, since this widget does not need a separate label.
    virtual CustomLabel* getOrCreateNameLabel() const;

protected:
    /// @see QPropertyWidget
    void showNameLabel(bool);

    virtual void resizeEvent(QResizeEvent* e);
protected slots:
    virtual void updateFromPropertySlot();

private:
    const VolumeBase* getVolumeFromProperty() const;
    std::string clipText(CustomLabel* l, QString fullText);

    QLabel* previewLabel_;

    CustomLabel* volumeNameLabel_;
    CustomLabel* pathLabel_;
    CustomLabel* dimensionLabel_;
    CustomLabel* spacingLabel_;
    CustomLabel* memSizeLabel_;

    CustomLabel* dimensionLabelCaption_;
    CustomLabel* spacingLabelCaption_;
    CustomLabel* memSizeLabelCaption_;

    static const std::string loggerCat_;

};

} // namespace voreen

#endif // VRN_VOLUMEINFOPROPERTYWIDGET_H
